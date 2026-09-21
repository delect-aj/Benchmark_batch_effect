source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
x <- PLSDAbatch::PLSDA_batch(X = clr(d$counts), Y.trt = d$meta$phenotype, Y.bat = d$meta$batch,
                              ncomp.trt = 1, ncomp.bat = nlevels(d$meta$batch) - 1, balance = TRUE)$X.nobatch
write_output(x, "clr", d$out)
