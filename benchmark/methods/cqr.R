source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
suppressMessages(library(doParallel))  # ConQuR uses foreach %do% without importing it; its vignette attaches doParallel
# Park & Park 2025 (Front Microbiol 16:1484183). The released code does not run as-is: the NB step reads a
# non-existent "Intercept" coefficient and the quantile step relies on undefined globals (batchid, standard_name).
# Reimplemented from the paper: reference batch = lowest robust CV (MAD/median) of library size, then composite
# quantile regression over k = 19 quantiles (5%..95%) via ConQuR, whose "composite" mode is the same CQR step.
# The NB mean-adjustment step is omitted. Report this as a partial reimplementation.
lib <- rowSums(d$counts)
rcv <- tapply(lib, d$meta$batch, function(x) mad(x, constant = 1) / median(x))
ref <- names(which.min(rcv))
# ConQuR bug: its fallback for sparse taxa (simple_QQ) reads a global `batchid` it never defines; reference level first
batchid <- relevel(d$meta$batch, ref = ref)
x <- ConQuR::ConQuR(tax_tab = d$counts, batchid = d$meta$batch,
                    covariates = d$meta[, "phenotype", drop = FALSE],
                    batch_ref = ref,
                    quantile_type = "composite", taus = seq(0.05, 0.95, by = 0.05))
write_output(x, "counts", d$out)
