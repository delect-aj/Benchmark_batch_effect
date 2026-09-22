# Diagnostic A1(c): ConQuR with the authors' tuning (Tune_ConQuR, pools as in the package vignette).
# Reference-batch pool: all batches, or the 3 largest when there are more than 5 (CRC: each fit ~13 min).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "../methods/_common.R"))
d <- read_input()
suppressMessages(library(doParallel))
tb <- sort(table(d$meta$batch), decreasing = TRUE)
pool <- names(tb)[seq_len(if (length(tb) > 5) 3 else length(tb))]
batchid <- d$meta$batch
res <- ConQuR::Tune_ConQuR(tax_tab = d$counts, batchid = d$meta$batch,
                           covariates = d$meta[, "phenotype", drop = FALSE],
                           batch_ref_pool = pool, logistic_lasso_pool = FALSE,
                           quantile_type_pool = c("standard", "lasso"), simple_match_pool = FALSE,
                           lambda_quantile_pool = c(NA, "2p/n"), interplt_pool = FALSE,
                           frequencyL = 0, frequencyU = 1, num_core = 8)
print(res$method_final)
write_output(res$tax_final, "counts", d$out)
