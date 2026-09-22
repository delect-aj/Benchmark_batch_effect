# Seed microbiome, 16S (Foxx & Rivers 2025, Mendeley Data 10.17632/5xrfg5dym6.1): 9 studies, 5 amplicon regions.
# The confounded case: phenotype = host plant (multi-class), and most plants appear in only some studies.
# Features are the published ASVs (sequence hashes, no taxonomy), kept as published. Batch = study.
# Usage: Rscript data/seed_16s.R <outdir> <raw download dir>
a <- commandArgs(TRUE); outdir <- a[1]; raw <- a[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE); dir.create(raw, recursive = TRUE, showWarnings = FALSE)
base <- "https://data.mendeley.com/public-files/datasets/5xrfg5dym6/files/"
get <- function(id, f) { p <- file.path(raw, f); if (!file.exists(p)) download.file(paste0(base, id, "/file_downloaded"), p, mode = "wb"); p }
x <- read.csv(get("7b71e396-90aa-45a6-9e38-c768778e1ab6", "count.csv"), check.names = FALSE)[, -1]
m <- read.csv(get("91155828-5f29-45ad-8c71-aa531f169f08", "metadata.csv"), check.names = FALSE)[, -1]
m <- m[!duplicated(m$SampleID), ]
rownames(x) <- x$SampleID; x <- as.matrix(x[, -1])
ids <- intersect(rownames(x), m$SampleID)
counts <- x[ids, ]; m <- m[match(ids, m$SampleID), ]
meta <- data.frame(batch = m$study_id, phenotype = m$host_plant, region = m$gene_region, row.names = ids)

# Same filter for every dataset (PLAN.md 2.3): prevalence >= 10%, drop empty samples
counts <- counts[, colMeans(counts > 0) >= 0.10, drop = FALSE]
ok <- rowSums(counts) > 0
counts <- counts[ok, ]; meta <- meta[ok, ]
w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(as.data.frame(counts), "counts.tsv"); w(meta, "meta.tsv")
print(table(meta$batch, meta$phenotype))
cat("seed_16s:", nrow(counts), "samples x", ncol(counts), "ASVs\n")
