# Toy data for the pilot only: tests the wrapper plumbing, not method quality.
# Real simulations (MIDASim / SparseDOSSA2 + batch injection) replace this. See PLAN.md section 2.1.
# Model: log-normal community + 10 DA taxa (4x in phenotype 1) + McLaren multiplicative per-batch per-taxon bias.
# Usage: Rscript toy_sim.R <outdir> [confounding 0..1]
a <- commandArgs(TRUE)
outdir <- a[1]
conf <- if (length(a) > 1) as.numeric(a[2]) else 0
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

n_batch <- 3; n_per <- 40; p <- 100; N <- n_batch * n_per
batch <- rep(seq_len(n_batch), each = n_per)
# conf = P(phenotype follows batch parity) beyond chance; 0 = balanced, 1 = fully confounded
pheno <- ifelse(runif(N) < conf, batch %% 2, rbinom(N, 1, 0.5))

base <- log(rgamma(p, 0.3) + 1e-6)
lp <- matrix(base, N, p, byrow = TRUE) + matrix(rnorm(N * p), N)
da <- 1:10
lp[pheno == 1, da] <- lp[pheno == 1, da] + log(4)
bias <- matrix(rnorm(n_batch * p), n_batch)

draw <- function(logit) {
  pr <- exp(logit); pr <- pr / rowSums(pr)
  depth <- rpois(N, 2e4)
  t(sapply(seq_len(N), function(i) rmultinom(1, depth[i], pr[i, ])))
}
ids <- sprintf("S%03d", seq_len(N)); taxa <- sprintf("t%03d", seq_len(p))
save_tab <- function(m, f) {
  dimnames(m) <- list(ids, taxa)
  write.table(data.frame(sample_id = ids, m), file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
}
save_tab(draw(lp + bias[batch, ]), "counts.tsv")
save_tab(draw(lp), "oracle.tsv")  # same biology, no batch bias
write.table(data.frame(sample_id = ids, batch = paste0("B", batch), phenotype = pheno),
            file.path(outdir, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(taxa[da], file.path(outdir, "truth_da.txt"))
