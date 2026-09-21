source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
x <- MMUPHin::adjust_batch(feature_abd = t(d$counts / rowSums(d$counts)), batch = "batch",
                              covariates = "phenotype", data = d$meta)$feature_abd_adj
write_output(t(x), "relabund", d$out)
