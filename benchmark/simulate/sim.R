# Simulated dataset with full ground truth (PLAN.md 2.1).
# Biology: MIDASim (parametric mode) fitted to a real template; DA taxa planted by shifting their mean relative abundance in cases.
# Batch: injected on top of the true relative abundances, so the batch-free truth (oracle.tsv) is known exactly.
# Usage: Rscript simulate/sim.R <template_counts.tsv> <outdir> [key=value ...]   (keys: see `p` below)
#        Rscript simulate/sim.R <template_counts.tsv> --fit     (fit + cache the template only; run once before a grid)
suppressMessages(library(MIDASim))
a <- commandArgs(TRUE)
template <- a[1]; outdir <- a[2]
p <- list(n_batch = 3,      # number of batches
          n_per = 50,       # mean samples per batch
          imbalance = 1,    # largest / smallest batch size
          conf = 0,         # batch-phenotype confounding: 0 balanced .. 1 fully confounded
          da_prop = 0.1,    # fraction of taxa that are DA (0 = null scenario)
          da_log2fc = 2,    # DA effect size, |log2 fold change| of mean relative abundance
          bias_sd = 1,      # multiplicative per-batch per-taxon bias, sd on log scale (McLaren model)
          affected = 1,     # fraction of taxa hit by the bias in each batch (<1 = non-systematic)
          scale_sd = 0,     # per-batch per-taxon variance scaling, sd on log scale (ComBat-style)
          dropout = 0,      # max extra per-batch dropout probability of present taxa
          depth_fold = 1,   # max ratio of sequencing depth between batches
          seed = 1)
for (kv in a[-(1:2)]) { kv <- strsplit(kv, "=")[[1]]; if (!kv[1] %in% names(p)) stop("unknown key ", kv[1]); p[[kv[1]]] <- as.numeric(kv[2]) }
set.seed(p$seed)
# --- Fit the template once, cache next to it ----------------------------------------------------
cache <- paste0(template, ".midasim.rds")
fit <- if (file.exists(cache)) readRDS(cache) else {
  tab <- as.matrix(read.delim(template, row.names = 1, check.names = FALSE))
  f <- MIDASim.setup(tab[, colSums(tab > 0) > 0], mode = "parametric")  # parametric: n and mean abundances can change freely
  saveRDS(f, cache); f
}
if (outdir == "--fit") quit(save = "no")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
taxa <- fit$taxa.names
P <- length(taxa)

# --- Design: batch sizes, phenotype with controlled confounding --------------------------------
w <- exp(seq(0, log(p$imbalance), length.out = p$n_batch))  # geometric spread of batch sizes
sizes <- round(p$n_per * p$n_batch * w / sum(w))
batch <- rep(seq_len(p$n_batch), sizes); N <- length(batch)
sgn <- ifelse(seq_len(p$n_batch) %% 2 == 1, 1, -1)       # odd batches lean case, even lean control
pheno <- rbinom(N, 1, 0.5 + p$conf / 2 * sgn[batch])
if (p$conf >= 1 && p$n_batch == 1) stop("full confounding needs >= 2 batches")

# DA taxa among reasonably prevalent ones, random direction
prev <- fit$taxa.1.prop
cand <- which(prev >= 0.2)
# da_prop is relative to taxa that survive the >= 10% prevalence filter below, not to all template taxa
da <- sample(cand, min(length(cand), round(p$da_prop * sum(prev >= 0.1))))
fc <- rep(1, P); fc[da] <- 2^(p$da_log2fc * sample(c(-1, 1), length(da), replace = TRUE))

# --- True relative abundances (MIDASim), per phenotype group -------------------------------------
sim_rel <- function(n, mult) {
  mra <- fit$mean.rel.abund * mult
  m <- MIDASim.modify(fit, lib.size = rep(1e5, n), mean.rel.abund = mra / sum(mra))  # n samples = length(lib.size)
  MIDASim(m, only.rel = TRUE)$sim_rel
}
rel <- matrix(0, N, P)
if (any(pheno == 0)) rel[pheno == 0, ] <- sim_rel(sum(pheno == 0), rep(1, P))
if (any(pheno == 1)) rel[pheno == 1, ] <- sim_rel(sum(pheno == 1), fc)
lib <- sample(fit$lib.size, N, replace = TRUE)

# --- Batch effects on the true relative abundances ------------------------------------------------
obs <- rel
depth <- exp(runif(p$n_batch, -1, 1) * log(p$depth_fold) / 2)
for (b in seq_len(p$n_batch)) {
  i <- batch == b
  hit <- runif(P) < p$affected
  w <- ifelse(hit, exp(rnorm(P, 0, p$bias_sd)), 1)
  x <- obs[i, , drop = FALSE]
  if (p$scale_sd > 0) {                                     # stretch log-abundance around each taxon's batch mean
    s <- ifelse(hit, exp(rnorm(P, 0, p$scale_sd)), 1)
    l <- log(x); l[!is.finite(l)] <- NA
    mu <- colMeans(l, na.rm = TRUE)
    l <- sweep(sweep(sweep(l, 2, mu), 2, s, "*"), 2, mu, "+")
    x <- ifelse(is.na(l), 0, exp(l))
  }
  x <- sweep(x, 2, w, "*")
  if (p$dropout > 0) x[x > 0 & matrix(runif(length(x)) < runif(1, 0, p$dropout), nrow(x))] <- 0
  obs[i, ] <- x
}
draw <- function(r, n) t(sapply(seq_len(N), function(k) if (sum(r[k, ]) > 0) rmultinom(1, n[k], r[k, ]) else rep(0, P)))
counts <- draw(obs, round(lib * depth[batch]))
oracle <- draw(rel, lib)

# --- Same filter as real data: prevalence >= 10%, drop empty samples ------------------------------
keep_t <- colMeans(counts > 0) >= 0.10
keep_s <- rowSums(counts[, keep_t]) > 0 & rowSums(oracle[, keep_t]) > 0
ids <- sprintf("S%04d", seq_len(N))[keep_s]
save_tab <- function(m, f) {
  m <- m[keep_s, keep_t, drop = FALSE]; dimnames(m) <- list(ids, taxa[keep_t])
  write.table(data.frame(sample_id = ids, m, check.names = FALSE), file.path(outdir, f),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
save_tab(counts, "counts.tsv")
save_tab(oracle, "oracle.tsv")
meta <- data.frame(sample_id = ids, batch = paste0("B", batch[keep_s]), phenotype = pheno[keep_s])
write.table(meta, file.path(outdir, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(intersect(taxa[da], taxa[keep_t]), file.path(outdir, "truth_da.txt"))

tb <- table(meta$batch, meta$phenotype)
p$cramers_v <- if (min(dim(tb)) > 1) sqrt(suppressWarnings(chisq.test(tb, correct = FALSE))$statistic / sum(tb)) else 1
p$n_samples <- nrow(meta); p$n_taxa <- sum(keep_t); p$n_da <- length(intersect(taxa[da], taxa[keep_t]))
p$template <- basename(template)
jsonlite::write_json(lapply(p, unname), file.path(outdir, "params.json"), auto_unbox = TRUE, pretty = TRUE)
