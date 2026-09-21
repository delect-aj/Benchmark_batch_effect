source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
x <- batchelor::fastMNN(t(clr(d$counts)), batch = d$meta$batch, d = min(20, nrow(d$counts) - 1))
write_output(SingleCellExperiment::reducedDim(x, "corrected"), "embedding", d$out)
