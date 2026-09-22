# A1 summary: false DA calls in the null scenario (10 reps) and CRC batch R2 for ConQuR variants.
# Usage: Rscript pilot/conqur_diag.R <results dir>
r <- commandArgs(TRUE)[1]
f <- list.files(r, "\\.score$", recursive = TRUE, full.names = TRUE)
f <- f[grepl("null_16s_conf0/|crc_mgx/", f) & grepl("/(raw|oracle|conqur|conqur_permcov|conqur_tuned)\\.score$", f)]
s <- do.call(rbind, lapply(f, function(x) { d <- read.delim(x); if (ncol(d) < 3) return(NULL)
  d$ds <- basename(dirname(x)); d$scen <- basename(dirname(dirname(x))); d[, c("method", "ds", "scen", "R2_batch_ait", "R2_pheno_ait", "iLISI", "DA_nsig")] }))
nl <- s[s$scen == "null_16s_conf0", ]
cat("== null_16s_conf0: false DA calls per rep (truth = none)\n")
print(tapply(nl$DA_nsig, list(nl$ds, nl$method), identity))
cat("\nmean false calls:\n"); print(round(tapply(nl$DA_nsig, nl$method, mean), 1))
cat("\nmean batch R2:\n"); print(round(tapply(nl$R2_batch_ait, nl$method, mean), 4))
cat("\n== CRC\n"); print(s[s$ds == "crc_mgx", c("method", "R2_batch_ait", "R2_pheno_ait", "iLISI")], row.names = FALSE, digits = 3)
