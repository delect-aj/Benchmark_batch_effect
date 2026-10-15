# Storage/handling experiment on pig faeces and sewage (ENA PRJEB31650): 4 biological matrices held under 21
# storage conditions with technical replicates, across 3 library preparations and 12 flow-cell lanes, half of
# them spiked with an 8-organism mock community. Published as Poulsen et al. 2021 Microbiol Spectr 9:e01387-21;
# re-used by Fachrul et al. 2022 Sci Rep 12:22236 (the "Salim 2022" entry in PLAN.md).
# No raw-read processing needed: the authors published the MGmapper genus count table on figshare (2.6 MB) and
# the per-sample conditions as Table S1 of the paper. Table S1 is taken from Europe PMC, which serves the
# supplementary files as one zip; PMC itself gates downloads behind a proof-of-work and ASM returns 403.
# For calibration only (PLAN.md 2.2) -- not a case-control dataset, so not in the method run list.
# Usage: Rscript data/salim_pig.R <outdir> <raw download dir>
a <- commandArgs(TRUE); outdir <- a[1]; raw <- a[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE); dir.create(raw, recursive = TRUE, showWarnings = FALSE)
get <- function(url, f) { p <- file.path(raw, f); if (!file.exists(p)) download.file(url, p, mode = "wb"); p }

# Mock community as mixed into each matrix: Table S3, "added cells to samples" (cells/mL), keyed by the genus
# names MGmapper reports. Pig_feces_1 got half as much P. freudenreichii and S. aureus as the other three.
mock <- c(Propionibacterium = 1e9, Bacteroides = 1e8, Staphylococcus = 1e9, Fusobacterium = 1e8,
          Escherichia = 1e8, Salmonella = 1e8, Saccharomyces = 2e7, Cryptosporidium = 2e6)

cnt <- read.delim(get("https://ndownloader.figshare.com/files/22065585", "Item_A_TaxonomyRaw.txt"),
                  row.names = 1, check.names = FALSE)   # genera x samples, raw counts (halved: mapped as proper pairs)
tax <- read.delim(get("https://ndownloader.figshare.com/files/22065591", "Item_C_FeatureShort.txt"),
                  row.names = 1, check.names = FALSE)   # same genera in a different order, so join on the genus name
names(tax) <- c("kingdom", "phylum", "class", "order", "family", "genus")
tax <- tax[match(rownames(cnt), tax$genus), ]; rownames(tax) <- tax$genus
stopifnot("genus names do not line up between the count and feature tables" = !anyNA(tax$genus))

s1 <- "spectrum.01387-21-s0001.xlsx"
unzip(get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8510183/supplementaryFiles", "pmc_suppl.zip"),
      files = s1, exdir = raw, overwrite = FALSE)
m <- as.data.frame(readxl::read_excel(file.path(raw, s1), sheet = 2, col_types = "text"))  # sheet 1 is the caption
m <- m[!is.na(m$Sample) & m$SpikedUnspiked != "EC", ]   # drop the 30 DNA-extraction blanks and 2 footnote rows

# Keep the microbial domains, plus the two eukaryotes of the mock community (yeast and Cryptosporidium);
# MGmapper also maps against human, animal and plant references, which are not part of the community.
counts <- t(as.matrix(cnt[tax$kingdom %in% c("Bacteria", "Archaea") | tax$genus %in% names(mock), ]))
ids <- intersect(rownames(counts), m$Sample)            # one count column comes from a later run absent from Table S1
counts <- counts[ids, , drop = FALSE]; m <- m[match(ids, m$Sample), ]
temp <- m[["Temperature (degree C)"]]
meta <- data.frame(batch = paste0(m$Time, "h_", ifelse(temp == "NA", "fresh", paste0(temp, "C"))),
                   phenotype = m$SampleMatrix,          # Pig_feces_1/2, Sewage_1/2: the real biological signal
                   replicate = m$Replicate_Groups,      # same biological sample, repeated: what RUV-III-NB needs
                   spiked = m$SpikedUnspiked, library_prep = m$Library_preparation, lane = m$Flow_cell_lane,
                   time_h = m$Time, temp_c = temp, freeze_thaw = m$Freeze_thaw_cycles, row.names = ids)

# Same filter for every dataset (PLAN.md 2.3): prevalence >= 10%, drop empty samples
counts <- counts[, colMeans(counts > 0) >= 0.10, drop = FALSE]
ok <- rowSums(counts) > 0
counts <- counts[ok, ]; meta <- meta[ok, ]

# Truth for the pure mock samples: cell proportions, as with brooks_mock. Genome-size differences stay inside
# the measured bias, which is what the McLaren per-taxon multiplicative model is meant to capture.
pure <- rownames(meta)[meta$spiked == "pureMC"]
truth <- t(sapply(pure, function(id) {
  v <- mock
  if (meta[id, "phenotype"] == "Pig_feces_1") v[c("Propionibacterium", "Staphylococcus")] <- 5e8
  v / sum(v)
}))

w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(as.data.frame(counts), "counts.tsv"); w(meta, "meta.tsv"); w(as.data.frame(truth), "truth.tsv")
w(tax[colnames(counts), ], "taxonomy.tsv")
print(table(meta$phenotype, meta$spiked))
cat("salim_pig:", nrow(counts), "samples x", ncol(counts), "genera,", nlevels(factor(meta$batch)),
    "storage conditions,", nlevels(factor(meta$replicate)), "replicate groups,", length(pure), "pure mock\n")
