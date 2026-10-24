# Score one wrapper output. Usage: Rscript metrics/compute.R <datadir> <out.tsv> <result.tsv>
# <datadir> has meta.tsv, and for simulations oracle.tsv (batch-free truth) and truth_da.txt (true DA taxa, may be empty).
# Writes a one-row TSV. Metrics that do not apply to the output kind are NA, never dropped.
suppressMessages({ library(vegan); library(cluster) })
a <- commandArgs(TRUE)
if (length(a) != 3) stop("usage: Rscript compute.R <datadir> <out.tsv> <result.tsv>")
meta <- read.delim(file.path(a[1], "meta.tsv"), row.names = 1)
x <- as.matrix(read.delim(a[2], row.names = 1, check.names = FALSE))
kind <- jsonlite::fromJSON(paste0(a[2], ".json"))$kind
stopifnot(identical(rownames(x), rownames(meta)), all(is.finite(x)))
if (kind != "embedding")
  stopifnot(identical(colnames(x), colnames(read.delim(file.path(a[1], "counts.tsv"), row.names = 1,
                                                       nrows = 1, check.names = FALSE))))
batch <- factor(meta$batch); pheno <- factor(meta$phenotype)

# --- Put every output kind into the spaces the metrics need -------------------------------------
prop <- function(m) { m <- pmax(m, 0); m / rowSums(m) }
clr_of_prop <- function(p) {  # zeros -> half the smallest positive proportion (same rule for every method)
  p[p == 0] <- min(p[p > 0]) / 2
  l <- log(p); l - rowMeans(l)
}
softmax <- function(l) { e <- exp(l - apply(l, 1, max)); e / rowSums(e) }

eucl <- switch(kind, counts = , relabund = clr_of_prop(prop(x)), x)          # Aitchison for compositions
comp <- switch(kind, counts = , relabund = prop(x), clr = , log = softmax(x),  # for Bray-Curtis
               percentile = x, embedding = NULL)
pa <- if (kind %in% c("counts", "relabund")) (x > 0) * 1 else NULL             # for Jaccard

D <- list(ait = dist(eucl),
          bc = if (!is.null(comp)) vegdist(comp, "bray"),
          jac = if (!is.null(pa)) vegdist(pa, "jaccard", binary = TRUE))

r2 <- function(d, g) if (is.null(d)) NA else adonis2(d ~ g, permutations = 0)$R2[1]
asw <- function(d, g) mean(silhouette(as.integer(g), d)[, 3])

# ponytail: unweighted kNN LISI (inverse Simpson of neighbour labels), not the perplexity-weighted original
lisi <- function(d, g, k = 15) {
  m <- as.matrix(d); k <- min(k, nrow(m) - 1)
  mean(apply(m, 1, function(r) {
    nb <- g[order(r)[2:(k + 1)]]
    1 / sum((table(nb) / k)^2)
  }))
}

res <- data.frame(method = sub("\\.tsv$", "", basename(a[2])), kind = kind)
for (nm in names(D)) {
  res[[paste0("R2_batch_", nm)]] <- r2(D[[nm]], batch)
  res[[paste0("R2_pheno_", nm)]] <- r2(D[[nm]], pheno)
}
res$ASW_batch <- asw(D$ait, batch)
res$ASW_pheno <- asw(D$ait, pheno)
res$iLISI <- (lisi(D$ait, batch) - 1) / (nlevels(batch) - 1)  # 1 = perfectly mixed batches
res$cLISI <- (nlevels(pheno) - lisi(D$ait, pheno)) / (nlevels(pheno) - 1)  # 1 = phenotypes fully separated

# --- Ground truth (simulations only) -----------------------------------------------------------
# One Wilcoxon pass over a table, in whatever space the caller hands it.
da_test <- function(space, group) {
  p <- apply(space, 2, function(v) suppressWarnings(wilcox.test(v ~ group, exact = FALSE)$p.value))
  p[is.na(p)] <- 1                                   # constant taxa
  list(p = p, sig = names(p)[p.adjust(p, "BH") < 0.05])
}
# Error rate, power and average precision of one call set against one definition of the truth.
against <- function(called, p, truth)
  list(FDR = if (length(called)) mean(!called %in% truth) else 0,
       power = if (length(truth)) mean(truth %in% called) else NA,
       AP = if (length(truth)) { hit <- names(sort(p)) %in% truth
                                 mean((cumsum(hit) / seq_along(hit))[hit]) } else NA)

oracle_f <- file.path(a[1], "oracle.tsv")
res$oracle_mantel <- res$oracle_taxon_rho <- NA
oracle_sig <- NULL
if (file.exists(oracle_f)) {
  o <- clr_of_prop(prop(as.matrix(read.delim(oracle_f, row.names = 1, check.names = FALSE))))
  stopifnot(identical(rownames(o), rownames(x)))
  res$oracle_mantel <- cor(c(D$ait), c(dist(o)))  # Mantel r with no permutations = correlation of distances
  if (kind != "embedding")
    res$oracle_taxon_rho <- mean(sapply(colnames(o), function(j) cor(eucl[, j], o[, j], method = "spearman")),
                                 na.rm = TRUE)
  if (nlevels(pheno) == 2) oracle_sig <- da_test(o, pheno)$sig
}

# Two truths, because they answer different questions.
#
# `truth_da.txt` is what the simulation planted: did the method recover the biology? On its own it
# overstates the error rate, because planting a fold change renormalises the community -- the taxa
# that were *not* planted really do shift, and counting them as errors blames the method for the
# simulation (measured: the planted effects change total mass by 1.5-2x, which moves every other
# taxon's prevalence by ~0.04).
#
# The oracle's own calls are what the same test finds on the batch-free table: did the method recover
# the answer the data would have given with no batch effect at all? That is the question a
# batch-correction benchmark actually asks, and it is well defined whatever the simulation did to the
# composition. By construction the oracle scores perfectly against it, so it bounds the scale rather
# than acting as a comparator.
truth_f <- file.path(a[1], "truth_da.txt")
res$DA_FDR <- res$DA_power <- res$DA_AP <- res$DA_nsig <- NA
res$DA_oracle_FDR <- res$DA_oracle_power <- res$DA_oracle_AP <- NA
if (kind != "embedding" && nlevels(pheno) == 2) {
  called <- da_test(eucl, pheno)
  res$DA_nsig <- length(called$sig)
  if (file.exists(truth_f)) {
    planted <- against(called$sig, called$p, readLines(truth_f))
    res$DA_FDR <- planted$FDR; res$DA_power <- planted$power; res$DA_AP <- planted$AP
  }
  if (!is.null(oracle_sig)) {
    recovered <- against(called$sig, called$p, oracle_sig)
    res$DA_oracle_FDR <- recovered$FDR
    res$DA_oracle_power <- recovered$power
    res$DA_oracle_AP <- recovered$AP
  }
}

write.table(res, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
