# Sensitivity variant: ConQuR default (methods/conqur_default.R), but the covariate is a PERMUTED phenotype (ConQuR requires a
# covariate). Same distribution, no link to the real phenotype: isolates the effect of conditioning on phenotype.
# If false DA calls in the null scenario disappear, the inflation comes from that conditioning.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
suppressMessages(library(doParallel))
set.seed(1)
ref <- levels(d$meta$batch)[1]
batchid <- relevel(d$meta$batch, ref = ref)
x <- ConQuR::ConQuR(tax_tab = d$counts, batchid = d$meta$batch,
                    covariates = data.frame(phenotype_perm = sample(d$meta$phenotype)), batch_ref = ref)
write_output(x, "counts", d$out)
