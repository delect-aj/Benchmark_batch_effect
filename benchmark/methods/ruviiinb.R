source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
# Replicate sets: meta$replicate if the dataset has technical replicates, else phenotype groups
# (the authors allow known sub-populations to define M). Negative controls: the 20% of taxa least associated with phenotype.
grp <- if ("replicate" %in% names(d$meta)) d$meta$replicate else d$meta$phenotype
M <- model.matrix(~ 0 + factor(grp))
pv <- apply(clr(d$counts), 2, function(x) wilcox.test(x ~ d$meta$phenotype, exact = FALSE)$p.value)
b <- as.numeric(d$meta$batch)
fit <- ruvIIInb::ruvIII.nb(Y = t(d$counts), M = M, ctl = pv >= quantile(pv, 0.8), k = 2, batch = b,
                           zeroinf = rep(TRUE, nrow(d$counts)))  # authors recommend ZINB for non-UMI data
# ruvIII.nb returns Mb per replicate group, but get.res slices it per sample; expand it the way ruvIII.nb does
# internally (Mb[, apply(M, 1, which)]) before asking for percentile-adjusted counts
fit$Mb <- fit$Mb[, apply(fit$M, 1, which), drop = FALSE]
x <- t(as.matrix(ruvIIInb::get.res(fit, type = "quantile", batch = b)))
dimnames(x) <- dimnames(d$counts)
write_output(x, "counts", d$out)
