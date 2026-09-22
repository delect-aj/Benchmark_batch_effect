source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
suppressMessages(library(doParallel))  # ConQuR uses foreach %do% without importing it; its vignette attaches doParallel
# Sensitivity variant: ConQuR untuned, reference batch = first level (the main entry is methods/conqur.R = Tune_ConQuR)
ref <- levels(d$meta$batch)[1]
# ConQuR bug: its fallback for sparse taxa (simple_QQ) reads a global `batchid` it never defines; reference level first
batchid <- relevel(d$meta$batch, ref = ref)
x <- ConQuR::ConQuR(tax_tab = d$counts, batchid = d$meta$batch,
                    covariates = d$meta[, "phenotype", drop = FALSE],
                    batch_ref = ref)
write_output(x, "counts", d$out)
