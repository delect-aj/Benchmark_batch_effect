# Validate every wrapper output against the contract, then print batch / phenotype PERMANOVA R2.
# Usage: Rscript check.R <datadir> <outdir>
# Fails (non-zero) if a present output breaks the contract, or if limma does not reduce batch R2 vs raw.
a <- commandArgs(TRUE)
meta <- read.delim(file.path(a[1], "meta.tsv"), row.names = 1)
counts <- as.matrix(read.delim(file.path(a[1], "counts.tsv"), row.names = 1))

clr <- function(x) { l <- log(pmax(x, 0) + 0.5); l - rowMeans(l) }
# Euclidean on CLR / embedding (Aitchison-like). ponytail: one distance only; the full metric set lives in metrics/
to_space <- function(x, kind) switch(kind,
  counts = , relabund = clr(if (kind == "relabund") x * 1e4 else x),
  x)

r2 <- function(x, g) vegan::adonis2(dist(x) ~ g, permutations = 0)$R2[1]

res <- list()
for (f in list.files(a[2], "\\.tsv$", full.names = TRUE)) {
  m <- sub("\\.tsv$", "", basename(f))
  x <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
  kind <- jsonlite::fromJSON(paste0(f, ".json"))$kind
  stopifnot(identical(rownames(x), rownames(meta)), all(is.finite(x)))
  if (kind != "embedding") stopifnot(identical(colnames(x), colnames(counts)))
  z <- to_space(x, kind)
  res[[m]] <- data.frame(method = m, kind = kind,
                         R2_batch = r2(z, meta$batch), R2_pheno = r2(z, factor(meta$phenotype)))
}
res <- do.call(rbind, res)
print(res[order(res$R2_batch), ], row.names = FALSE, digits = 3)

stopifnot("raw and limma outputs required" = all(c("raw", "limma") %in% res$method))
stopifnot("limma should lower batch R2 vs raw" = res["limma", "R2_batch"] < res["raw", "R2_batch"])
cat("check passed\n")
