# Shared I/O for every R wrapper. Contract: see benchmark/README.md
# Usage inside a wrapper:  d <- read_input();  ...;  write_output(x, "clr", d$out)

KINDS <- c("counts", "relabund", "clr", "log", "percentile", "embedding")

read_input <- function() {
  a <- commandArgs(TRUE)
  if (length(a) != 3) stop("usage: Rscript <wrapper>.R counts.tsv meta.tsv out.tsv")
  counts <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
  meta <- read.delim(a[2], row.names = 1, check.names = FALSE)
  if (!identical(rownames(counts), rownames(meta))) stop("counts and meta sample_id differ / out of order")
  if (!all(c("batch", "phenotype") %in% names(meta))) stop("meta needs columns: batch, phenotype")
  meta$batch <- factor(meta$batch)
  meta$phenotype <- factor(meta$phenotype)
  list(counts = counts, meta = meta, out = a[3])
}

write_output <- function(x, kind, out) {
  if (!kind %in% KINDS) stop("unknown kind: ", kind)
  x <- as.matrix(x)
  if (kind == "embedding") colnames(x) <- paste0("dim", seq_len(ncol(x)))
  write.table(data.frame(sample_id = rownames(x), x, check.names = FALSE),
              out, sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(sprintf('{"kind": "%s"}', kind), paste0(out, ".json"))
}

clr <- function(counts, pseudo = 0.5) {
  l <- log(counts + pseudo)
  l - rowMeans(l)
}

# Protected covariate design, used identically by every method that accepts one
bio_design <- function(meta) model.matrix(~ phenotype, meta)
