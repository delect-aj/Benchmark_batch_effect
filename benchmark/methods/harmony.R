source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
pcs <- prcomp(clr(d$counts), rank. = min(20, nrow(d$counts) - 1))$x
x <- harmony::RunHarmony(pcs, d$meta, vars_use = "batch")
write_output(x, "embedding", d$out)
