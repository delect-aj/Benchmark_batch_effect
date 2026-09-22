# Mock communities (Brooks et al. 2015, bundled with metacal): known mixtures of 7 vaginal species on 6 plates.
# For calibration only (PLAN.md 2.2): checks whether per-batch multiplicative bias (the simulation's batch model)
# describes a real protocol effect. Not a case-control dataset, so it is not in the method run list.
# Writes counts.tsv (observed reads), truth.tsv (actual proportions), meta.tsv (batch = plate).
# Usage: Rscript data/brooks_mock.R <outdir>
outdir <- commandArgs(TRUE)[1]; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ext <- system.file("extdata", package = "metacal")
obs <- read.csv(file.path(ext, "brooks2015-observed.csv"), row.names = 1, check.names = FALSE)
act <- read.csv(file.path(ext, "brooks2015-actual.csv"), row.names = 1, check.names = FALSE)
sd <- read.csv(file.path(ext, "brooks2015-sample-data.csv"), row.names = 1)
sp <- setdiff(intersect(colnames(obs), colnames(act)), "Other")
ids <- Reduce(intersect, list(rownames(obs), rownames(act), rownames(sd)))
meta <- data.frame(batch = paste0("plate", sd[ids, "Plate"]), mixture_type = sd[ids, "Mixture_type"],
                   num_species = sd[ids, "Num_species"], row.names = ids)
w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(obs[ids, sp], "counts.tsv"); w(act[ids, sp], "truth.tsv"); w(meta, "meta.tsv")
print(table(meta$batch, meta$mixture_type))
cat("brooks_mock:", length(ids), "samples x", length(sp), "species\n")
