# Scenario table for the simulation grid (PLAN.md 2.1). Usage: Rscript simulate/design.R <n_random> <out.tsv>
# 1. sweep_*: confounding sweep, everything else at defaults, both templates (the main figure)
# 2. null_*:  no biological signal (over-correction / leakage check)
# 3. lhs_*:   balanced random design: every factor's levels appear equally often, levels permuted
#             independently per factor (discrete Latin hypercube). Covers interactions without the full factorial.
a <- commandArgs(TRUE); n <- as.integer(a[1]); set.seed(2026)

levels <- list(template = c("16s", "mgx"), n_batch = c(2, 5, 10), n_per = c(20, 50, 200), imbalance = c(1, 3),
               conf = c(0, 0.3, 0.6, 0.9, 1), da_prop = c(0.05, 0.1), da_log2fc = c(1, 2, 3),
               bias_sd = c(0.5, 1, 2), affected = c(0.1, 0.5, 1), scale_sd = c(0, 0.5),
               dropout = c(0, 0.2), depth_fold = c(1, 4))
default <- list(template = "16s", n_batch = 5, n_per = 50, imbalance = 1, conf = 0, da_prop = 0.1,
                da_log2fc = 2, bias_sd = 1, affected = 1, scale_sd = 0, dropout = 0, depth_fold = 1)
row <- function(id, ...) { r <- modifyList(default, list(...)); data.frame(scenario = id, r) }

sweep <- do.call(rbind, unlist(lapply(levels$template, function(t) lapply(levels$conf, function(c)
  row(sprintf("sweep_%s_conf%s", t, c), template = t, conf = c))), recursive = FALSE))
null <- do.call(rbind, lapply(c(0, 0.6, 1), function(c) row(sprintf("null_16s_conf%s", c), conf = c, da_prop = 0)))
lhs <- as.data.frame(lapply(levels, function(l) sample(rep_len(l, n))), stringsAsFactors = FALSE)
lhs <- cbind(scenario = sprintf("lhs_%03d", seq_len(n)), lhs)

d <- rbind(sweep, null, lhs)
write.table(d, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(d), "scenarios ->", a[2], "\n")
