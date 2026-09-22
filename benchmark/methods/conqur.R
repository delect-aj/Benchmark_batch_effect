source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
suppressMessages(library(doParallel))  # ConQuR uses foreach %do% without importing it; its vignette attaches doParallel
# Main ConQuR entry: the authors' tuning (Tune_ConQuR), search pools as in the package vignette.
# Reference-batch pool: the 3 largest batches (all if fewer), to bound runtime. Pre-run check (PLAN.md): the untuned
# default with an arbitrary reference fills in zeros of zero-rich batches and inflates zero-sensitive batch distances.
tb <- sort(table(d$meta$batch), decreasing = TRUE)
pool <- names(tb)[seq_len(min(3, length(tb)))]
batchid <- d$meta$batch  # ConQuR bug: its sparse-taxon fallback (simple_QQ) reads a global `batchid`
# Tune_ConQuR calls vegan::adonis (removed in vegan 2.7) and only reads the batch-term R2 ($aov.tab[1, 5]).
# Shim: same R2 from adonis2 (Bray-Curtis, the old default), in the old shape, where ConQuR finds its vegan imports.
shim <- function(formula, ...) {
  e <- environment(formula)  # adonis2 does not look up the LHS here, so evaluate both sides ourselves
  Y <- eval(formula[[2]], e); g <- eval(formula[[3]], e)
  r <- vegan::adonis2(vegan::vegdist(Y, "bray") ~ g, permutations = 0)
  list(aov.tab = data.frame(Df = r$Df, SumsOfSqs = r$SumOfSqs, MeanSqs = NA, F.Model = NA, R2 = r$R2))
}
imp <- parent.env(asNamespace("ConQuR")); unlockBinding("adonis", imp); assign("adonis", shim, envir = imp)
res <- ConQuR::Tune_ConQuR(tax_tab = d$counts, batchid = d$meta$batch,
                           covariates = d$meta[, "phenotype", drop = FALSE],
                           batch_ref_pool = pool, logistic_lasso_pool = FALSE,
                           quantile_type_pool = c("standard", "lasso"), simple_match_pool = FALSE,
                           lambda_quantile_pool = c(NA, "2p/n"), interplt_pool = FALSE,
                           frequencyL = 0, frequencyU = 1, num_core = 2)
print(res$method_final)
write_output(res$tax_final, "counts", d$out)
