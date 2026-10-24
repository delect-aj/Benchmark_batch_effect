# Known-composition scoring: how far is a corrected table from the composition that was actually pipetted?
#
# Applies to the designed mock datasets, which ship a truth.tsv of expected relative abundances
# (samples x taxa). Only the samples listed there are scored, so a dataset that mixes mock and
# biological specimens (jmbc_extraction) contributes its mock samples alone -- and on that subset every
# difference between batches is spurious by construction, which is what makes R2_batch and false_da
# below interpretable as false structure rather than as biology.
#
# Embeddings have no taxon axis and cannot be compared with a composition; the workflow does not ask
# for them here (see the truth rule in the Snakefile).
#
# Usage: Rscript metrics/truth.R <datadir> <corrected.tsv> <summary.tsv> <per_taxon.tsv>
suppressMessages(library(vegan))

as_matrix <- function(path) as.matrix(read.delim(path, row.names = 1, check.names = FALSE))

# Put each output kind back on the simplex. Percentile tables are not compositions; rescaling them to
# sum to one compares their shape only, which is the same compromise metrics/compute.R makes.
as_proportion <- function(values, kind) switch(kind,
  counts = , relabund = , percentile = { v <- pmax(values, 0); v / rowSums(v) },
  clr = , log = { e <- exp(values - apply(values, 1, max)); e / rowSums(e) })

# One rule for every method: a zero where the truth says the taxon is present is a miss, not a
# structural absence, so it is replaced by half the smallest positive value in that sample rather than
# dropped. Taxa the truth puts at zero are genuinely absent and are excluded from the sample's support.
closed_on_support <- function(proportions, support) {
  p <- proportions[support]
  p[p == 0] <- min(p[p > 0]) / 2
  p / sum(p)
}

clr_rows <- function(p) { l <- log(p); l - rowMeans(l) }

aitchison <- function(observed, expected) {
  lo <- log(observed) - mean(log(observed))
  le <- log(expected) - mean(log(expected))
  sqrt(sum((lo - le)^2))
}

bray <- function(observed, expected) sum(abs(observed - expected)) / 2

# Per sample: distance to the known composition, the mass that landed outside it, and the per-taxon
# log2 ratio that the multiplicative bias model predicts should be constant within a batch.
score_sample <- function(sample_proportions, truth_row) {
  support <- truth_row > 0
  observed <- closed_on_support(sample_proportions[names(truth_row)], support)
  expected <- truth_row[support] / sum(truth_row[support])
  list(aitchison = aitchison(observed, expected),
       bray = bray(observed, expected),
       off_truth_mass = 1 - sum(sample_proportions[names(truth_row)[support]]),
       log2_ratio = log2(observed / expected))
}

main <- function(args) {
  datadir <- args[1]; corrected <- args[2]; summary_out <- args[3]; taxon_out <- args[4]
  dataset <- sub("^results/data/", "", datadir)
  method <- sub("\\.tsv$", "", basename(corrected))

  truth <- as_matrix(file.path(datadir, "truth.tsv"))
  meta <- read.delim(file.path(datadir, "meta.tsv"), row.names = 1, check.names = FALSE)
  kind <- jsonlite::fromJSON(paste0(corrected, ".json"))$kind
  proportions <- as_proportion(as_matrix(corrected), kind)

  samples <- rownames(truth)
  taxa <- colnames(truth)
  stopifnot(all(samples %in% rownames(proportions)), all(taxa %in% colnames(proportions)))
  scored <- lapply(samples, function(s) score_sample(proportions[s, ], truth[s, taxa]))
  names(scored) <- samples

  # Per-taxon bias: the mean log2 ratio over the samples where the taxon was pipetted in, and how much
  # that ratio moves between batches -- the residual the correction was supposed to remove.
  batch <- factor(meta[samples, "batch"])
  ratios <- lapply(taxa, function(t) {
    present <- vapply(scored, function(s) t %in% names(s$log2_ratio), TRUE)
    data.frame(taxon = t,
               n_samples = sum(present),
               log2_ratio = mean(vapply(scored[present], function(s) s$log2_ratio[[t]], 0)),
               batch_sd = sd(tapply(vapply(scored[present], function(s) s$log2_ratio[[t]], 0),
                                    droplevels(batch[present]), mean)))
  })
  per_taxon <- cbind(dataset = dataset, method = method, do.call(rbind, ratios))

  # False structure: with one material and one pipetted composition, any difference between batches is
  # an artefact. Kruskal-Wallis across all batches avoids inventing a two-group split. Taken over the
  # whole truth taxon set, including taxa the truth puts at zero, where anything observed is noise.
  observed_full <- proportions[samples, taxa, drop = FALSE]
  observed_full[observed_full == 0] <- min(observed_full[observed_full > 0]) / 2
  observed_clr <- clr_rows(observed_full / rowSums(observed_full))
  kw <- apply(observed_clr, 2, function(v) kruskal.test(v ~ batch)$p.value)
  summary_row <- data.frame(
    dataset = dataset, method = method, kind = kind,
    n_samples = length(samples), n_taxa = length(taxa), n_batch = nlevels(batch),
    aitchison_to_truth = mean(vapply(scored, function(s) s$aitchison, 0)),
    bray_to_truth = mean(vapply(scored, function(s) s$bray, 0)),
    bias_abs_mean = mean(abs(per_taxon$log2_ratio), na.rm = TRUE),   # NA where the truth never contains the taxon
    bias_batch_sd = mean(per_taxon$batch_sd, na.rm = TRUE),
    off_truth_mass = mean(vapply(scored, function(s) s$off_truth_mass, 0)),
    R2_batch = adonis2(dist(observed_clr) ~ batch, permutations = 0)$R2[1],
    false_da_n = sum(p.adjust(kw, "BH") < 0.05))

  write.table(summary_row, summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(per_taxon, taxon_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat(sprintf("%s / %s: Aitchison to truth %.2f, |bias| %.2f, false DA %d/%d\n", dataset, method,
              summary_row$aitchison_to_truth, summary_row$bias_abs_mean,
              summary_row$false_da_n, summary_row$n_taxa))
}

main(commandArgs(TRUE))
