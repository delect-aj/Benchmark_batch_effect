source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
# Gibbons 2018: within each batch, each taxon -> percentile of that batch control distribution
ra <- d$counts / rowSums(d$counts)
ra[ra == 0] <- runif(sum(ra == 0), 0, 1e-9)  # jitter zeros as in the original implementation
ctrl <- d$meta$phenotype == "0"
out <- ra
for (b in levels(d$meta$batch)) {
  i <- d$meta$batch == b
  if (!any(i & ctrl)) stop("batch ", b, " has no controls")
  for (j in seq_len(ncol(ra))) {
    cd <- sort(ra[i & ctrl, j])
    x <- ra[i, j]
    out[i, j] <- 50 * (findInterval(x, cd, left.open = TRUE) + findInterval(x, cd)) / length(cd)  # scipy kind="mean"
  }
}
write_output(out, "percentile", d$out)
