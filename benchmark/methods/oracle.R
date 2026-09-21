source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
# Upper bound, simulations only: the batch-free truth. Fails on real data (no oracle.tsv), by design.
f <- side_file(d, "oracle.tsv")
if (is.null(f)) stop("no oracle.tsv: oracle exists only for simulated data")
write_output(as.matrix(read.delim(f, row.names = 1, check.names = FALSE)), "counts", d$out)
