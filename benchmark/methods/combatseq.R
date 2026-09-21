source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
x <- sva::ComBat_seq(t(d$counts), batch = d$meta$batch, group = d$meta$phenotype)
write_output(t(x), "counts", d$out)
