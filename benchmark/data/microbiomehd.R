# Real 16S datasets from MicrobiomeHD (Duvallet et al. 2017, Zenodo record 569601), one batch per study.
# OTUs are called de novo per study, so they are not comparable across studies: aggregate to genus via the RDP
# lineage (PLAN.md 2.3). Study-level sample filters follow the authors' db/dataset_info.yaml.
# Usage: Rscript data/microbiomehd.R <crc_16s|ibd_16s|hiv_16s> <outdir> <raw download dir>
a <- commandArgs(TRUE); name <- a[1]; outdir <- a[2]; raw <- a[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE); dir.create(raw, recursive = TRUE, showWarnings = FALSE)

# Case/control labels follow the MicrobiomeHD conventions; anything else (e.g. adenoma "nonCRC") is dropped
sets <- list(crc_16s = list(studies = c("crc_baxter", "crc_xiang", "crc_zackular", "crc_zeller", "crc_zhao"),
                            case = "CRC", control = "H"),
             ibd_16s = list(studies = c("ibd_alm", "ibd_engstrand_maxee", "ibd_gevers_2014", "ibd_huttenhower"),
                            case = c("CD", "UC"), control = c("H", "nonIBD")),
             hiv_16s = list(studies = c("hiv_dinh", "hiv_lozupone", "hiv_noguerajulian"),
                            case = "HIV", control = "H"))
s <- sets[[name]]; if (is.null(s)) stop("unknown dataset ", name)
# Per-study sample filters, copied from the authors' db/dataset_info.yaml (microbiomeHD repo); for these 12
# studies the disease label column is always DiseaseState. (That yaml has duplicate keys R's parser rejects.)
filters <- list(hiv_lozupone = list(time_point = "1"), hiv_noguerajulian = list(cohort = c("BCN0", "STK")))

read_study <- function(st) {
  tgz <- file.path(raw, paste0(st, "_results.tar.gz"))
  if (!file.exists(tgz))
    download.file(sprintf("https://zenodo.org/records/569601/files/%s_results.tar.gz?download=1", st), tgz, mode = "wb")
  untar(tgz, exdir = raw)
  d <- file.path(raw, paste0(st, "_results"))
  otu <- read.delim(list.files(file.path(d, "RDP"), "rdp_assigned$", full.names = TRUE), row.names = 1,
                    check.names = FALSE)                       # OTUs x samples, row names = RDP lineage
  m <- read.delim(list.files(d, "metadata.txt$", full.names = TRUE), row.names = 1, colClasses = "character",
                  check.names = FALSE)  # no re-encoding: latin1 conversion truncates some files in a C locale
  lab <- "DiseaseState"
  for (col in names(filters[[st]])) m <- m[m[[col]] %in% filters[[st]][[col]], , drop = FALSE]
  st_lab <- trimws(m[[lab]])
  m <- m[st_lab %in% c(s$case, s$control), , drop = FALSE]
  ids <- intersect(colnames(otu), rownames(m))
  x <- t(as.matrix(otu[, ids, drop = FALSE]))
  x <- x[rowSums(x) >= 100, , drop = FALSE]                    # authors' own floor: >= 100 reads per sample
  # genus key = lineage up to g__ (unclassified genera keep their deepest known rank)
  genus <- sub(";s__.*$", "", colnames(x))
  g <- t(rowsum(t(x), genus))
  list(counts = g, phenotype = as.integer(trimws(m[rownames(g), lab]) %in% s$case))
}
parts <- lapply(s$studies, read_study); names(parts) <- s$studies

taxa <- Reduce(union, lapply(parts, function(p) colnames(p$counts)))
counts <- do.call(rbind, lapply(names(parts), function(st) {
  p <- parts[[st]]
  m <- matrix(0, nrow(p$counts), length(taxa), dimnames = list(paste0(st, "_", rownames(p$counts)), taxa))
  m[, colnames(p$counts)] <- p$counts
  m
}))
meta <- data.frame(batch = rep(names(parts), sapply(parts, function(p) nrow(p$counts))),
                   phenotype = unlist(lapply(parts, `[[`, "phenotype")), row.names = rownames(counts))

# Same filter for every dataset (PLAN.md 2.3): prevalence >= 10%, drop empty samples
counts <- counts[, colMeans(counts > 0) >= 0.10]
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
cat(name, ":", nrow(counts), "samples x", ncol(counts), "genera\n")
