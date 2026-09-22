# Real shotgun datasets from curatedMetagenomicData, one batch per study. Writes counts.tsv, meta.tsv, taxonomy.tsv
# in the wrapper contract. Real data has no oracle/truth: it is scored on proxies only. Runs in the bench-data env.
# Usage: Rscript data/cmd.R <crc_mgx|ibd_mgx|ici_mgx> <outdir>
suppressMessages(library(curatedMetagenomicData))
a <- commandArgs(TRUE); name <- a[1]; outdir <- a[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Per dataset: studies, which samples to keep, and the binary phenotype (1 = case)
sets <- list(
  crc_mgx = list(studies = c("FengQ_2015", "ZellerG_2014", "YuJ_2015", "VogtmannE_2016", "WirbelJ_2018",   # Wirbel 2019 cohorts
                             "ThomasAM_2018a", "ThomasAM_2018b", "YachidaS_2019", "HanniganGD_2017"),
                 keep = function(cd) cd$study_condition %in% c("CRC", "control"),
                 pheno = function(cd) cd$study_condition == "CRC"),
  ibd_mgx = list(studies = c("HMP_2019_ibdmdb", "HallAB_2017", "IjazUZ_2017", "LiJ_2014", "NielsenHB_2014"),  # cohorts with cases and controls
                 keep = function(cd) cd$study_condition %in% c("IBD", "control"),
                 pheno = function(cd) cd$study_condition == "IBD"),
  ici_mgx = list(studies = c("FrankelAE_2017", "LeeKA_2022", "WindTT_2020"),  # immune-checkpoint response (ORR)
                 keep = function(cd) cd$ORR %in% c("yes", "no"),
                 pheno = function(cd) cd$ORR == "yes"))
s <- sets[[name]]; if (is.null(s)) stop("unknown dataset ", name)

tse <- lapply(s$studies, function(st) {
  x <- curatedMetagenomicData(paste0(st, ".relative_abundance"), dryrun = FALSE, counts = TRUE, rownames = "short")[[1]]
  cd <- as.data.frame(SummarizedExperiment::colData(x))
  keep <- s$keep(cd)
  # one sample per subject: the earliest collection (longitudinal cohorts such as HMP2)
  t <- if ("days_from_first_collection" %in% names(cd)) cd$days_from_first_collection else rep(NA, nrow(cd))
  o <- order(!keep, is.na(t), t)
  first <- rep(FALSE, nrow(cd)); first[o[!duplicated(cd$subject_id[o])]] <- TRUE
  x[, keep & first]
})
names(tse) <- s$studies

taxa <- Reduce(union, lapply(tse, rownames))
counts <- do.call(rbind, lapply(tse, function(x) {
  m <- matrix(0, ncol(x), length(taxa), dimnames = list(colnames(x), taxa))
  m[, rownames(x)] <- t(as.matrix(SummarizedExperiment::assay(x)))
  m
}))
meta <- data.frame(batch = rep(s$studies, sapply(tse, ncol)),
                   phenotype = as.integer(unlist(lapply(tse, function(x) s$pheno(as.data.frame(SummarizedExperiment::colData(x)))))),
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
tax[is.na(tax)] <- "unclassified"  # MetaDICT compares all higher ranks; NA ranks make its neighbour count NA

w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(as.data.frame(counts), "counts.tsv"); w(meta, "meta.tsv"); w(tax, "taxonomy.tsv")
print(table(meta$batch, meta$phenotype))
cat(name, ":", nrow(counts), "samples x", ncol(counts), "taxa\n")
