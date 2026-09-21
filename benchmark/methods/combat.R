source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
x <- sva::ComBat(t(clr(d$counts)), batch = d$meta$batch, mod = bio_design(d$meta))
write_output(t(x), "clr", d$out)
