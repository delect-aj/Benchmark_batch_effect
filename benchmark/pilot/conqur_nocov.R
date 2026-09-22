# Diagnostic A1(b): ConQuR exactly as methods/conqur.R but WITHOUT the protected phenotype covariate.
# If false DA calls in the null scenario disappear, the inflation comes from conditioning on phenotype.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "../methods/_common.R"))
d <- read_input()
suppressMessages(library(doParallel))
ref <- levels(d$meta$batch)[1]
batchid <- relevel(d$meta$batch, ref = ref)
x <- ConQuR::ConQuR(tax_tab = d$counts, batchid = d$meta$batch,
                    covariates = NULL, batch_ref = ref)  # NULL: X = data.frame(batchid) inside ConQuR; a constant dummy would make every fit singular and silently return raw counts
write_output(x, "counts", d$out)
