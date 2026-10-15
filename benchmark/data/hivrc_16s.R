# HIVRC: 17 HIV gut 16S studies (HIV Re-analysis Consortium; the DEBIAS-M HIV benchmark).
# Source files are downloaded manually from Synapse (needs an account) into results/raw/hivrc/ — see
# data/MANUAL_DATASETS.md. Upstream processing is the consortium's: one common Resphera Insight pipeline.
# Usage: Rscript data/hivrc_16s.R <outdir> <raw dir with the Synapse files>
a <- commandArgs(TRUE); outdir <- a[1]; raw <- a[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Genus table: rows = full lineage, columns = samples, already rarefied to 20,000 reads per sample
g <- read.delim(file.path(raw, "insight.merged_genus.txt"), row.names = 1, check.names = FALSE)
m <- as.data.frame(readxl::read_excel(file.path(raw, "SupplementaryMaterial.xlsx"), sheet = 1, skip = 1))
m <- m[, !grepl("^\\.\\.\\.", names(m)), drop = FALSE]
m$SeqID <- as.character(m$SeqID)

ids <- intersect(colnames(g), m$SeqID)
stopifnot("no samples shared between the abundance table and the manifest" = length(ids) > 0)
counts <- round(t(as.matrix(g[, ids, drop = FALSE])))
m <- m[match(ids, m$SeqID), ]
keep_cols <- intersect(c("Age", "gender", "msm", "ARTuse", "vl", "cd4"), names(m))
meta <- data.frame(batch = m$Study, phenotype = as.integer(m$hivstatus), m[, keep_cols, drop = FALSE],
                   row.names = ids)

# Same filter for every dataset (PLAN.md 2.3): prevalence >= 10%, drop empty samples
counts <- counts[, colMeans(counts > 0) >= 0.10, drop = FALSE]
ok <- rowSums(counts) > 0
counts <- counts[ok, ]; meta <- meta[ok, ]

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
tax <- do.call(rbind, lapply(strsplit(colnames(counts), ";"), function(l) {
  v <- sub("^[a-z]__", "", l)[seq_along(ranks)]
  v[is.na(v) | v == ""] <- "unclassified"
  v
}))
tax <- data.frame(tax, row.names = colnames(counts)); names(tax) <- ranks

w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(as.data.frame(counts), "counts.tsv"); w(meta, "meta.tsv"); w(tax, "taxonomy.tsv")
print(table(meta$batch, meta$phenotype))
cat("hivrc_16s:", nrow(counts), "samples x", ncol(counts), "genera,", nlevels(factor(meta$batch)), "studies\n")
