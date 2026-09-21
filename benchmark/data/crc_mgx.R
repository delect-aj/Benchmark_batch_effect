# Real pilot dataset: CRC vs control shotgun metagenomes, one batch per study (Wirbel 2019 cohorts in curatedMetagenomicData).
# Writes counts.tsv, meta.tsv, taxonomy.tsv in the wrapper contract. No oracle / truth: real data is scored on proxies only.
# Usage: Rscript data/crc_mgx.R <outdir>
suppressMessages(library(curatedMetagenomicData))
outdir <- commandArgs(TRUE)[1]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

studies <- c("FengQ_2015", "ZellerG_2014", "YuJ_2015", "VogtmannE_2016", "WirbelJ_2018",
             "ThomasAM_2018a", "ThomasAM_2018b", "YachidaS_2019", "HanniganGD_2017")
tse <- lapply(studies, function(s) {
  x <- curatedMetagenomicData(paste0(s, ".relative_abundance"), dryrun = FALSE, counts = TRUE, rownames = "short")[[1]]
  cd <- SummarizedExperiment::colData(x)
  keep <- cd$study_condition %in% c("CRC", "control")
  if ("days_from_first_collection" %in% names(cd))  # one sample per subject
    keep <- keep & (is.na(cd$days_from_first_collection) | cd$days_from_first_collection == 0)
  x[, keep & !duplicated(cd$subject_id)]
})
names(tse) <- studies

taxa <- Reduce(union, lapply(tse, rownames))
counts <- do.call(rbind, lapply(tse, function(x) {
  m <- matrix(0, ncol(x), length(taxa), dimnames = list(colnames(x), taxa))
  m[, rownames(x)] <- t(as.matrix(SummarizedExperiment::assay(x)))
  m
}))
meta <- data.frame(batch = rep(studies, sapply(tse, ncol)),
                   phenotype = as.integer(unlist(lapply(tse, function(x) x$study_condition == "CRC"))),
                   row.names = rownames(counts))

# Same filter for every dataset (PLAN.md 2.3): prevalence >= 10%, drop empty samples
counts <- counts[, colMeans(counts > 0) >= 0.10]
ok <- rowSums(counts) > 0
counts <- round(counts[ok, ]); meta <- meta[ok, ]

# Key on an explicit taxon column: rbind() silently renames duplicate row names, which broke a row-name lookup
tax <- do.call(rbind, lapply(tse, function(x)
  data.frame(taxon = rownames(x), as.data.frame(SummarizedExperiment::rowData(x)), row.names = NULL)))
tax <- tax[!duplicated(tax$taxon), ]
rownames(tax) <- tax$taxon
tax <- tax[colnames(counts), -1, drop = FALSE]
stopifnot("taxonomy lookup failed" = !anyNA(tax[[ncol(tax)]]))

w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(as.data.frame(counts), "counts.tsv"); w(meta, "meta.tsv"); w(tax, "taxonomy.tsv")
print(table(meta$batch, meta$phenotype))
cat(nrow(counts), "samples x", ncol(counts), "taxa\n")
