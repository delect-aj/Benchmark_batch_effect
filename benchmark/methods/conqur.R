source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
# ponytail: reference batch = first level; authors suggest picking one deliberately, revisit per dataset
x <- ConQuR::ConQuR(tax_tab = d$counts, batchid = d$meta$batch,
                    covariates = d$meta[, "phenotype", drop = FALSE],
                    batch_ref = levels(d$meta$batch)[1])
write_output(x, "counts", d$out)
