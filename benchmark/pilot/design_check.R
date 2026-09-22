# A2/A3: is the simulation design sane? Uses raw + oracle scores of every scenario (rep1).
# Usage: Rscript pilot/design_check.R results/scores.tsv
s <- read.delim(commandArgs(TRUE)[1])
s <- s[grepl("^sim/", s$dataset) & grepl("/rep1$", s$dataset) & s$status == "ok", ]
raw <- s[s$method == "raw", ]; ora <- s[s$method == "oracle", ]
cat("scenarios with raw scores:", nrow(raw), " oracle:", nrow(ora), "\n")
q <- function(v) round(quantile(v, c(0, .1, .25, .5, .75, .9, 1), na.rm = TRUE), 3)
cat("\n== A2: raw batch R2 (Aitchison), all scenarios. Real: CRC 0.077; HIVRC 0.12-0.30 (literature)\n"); print(q(raw$R2_batch_ait))
by <- function(d, v, f) round(tapply(d[[v]], d[[f]], median, na.rm = TRUE), 3)
for (f in c("bias_sd", "affected", "scale_sd", "dropout", "depth_fold", "n_batch", "template"))
  { cat("\nraw batch R2 median by", f, "\n"); print(by(raw, "R2_batch_ait", f)) }
cat("\nshare of scenarios with raw batch R2 > 0.30:", round(mean(raw$R2_batch_ait > 0.30), 3), "\n")
cat("\n== A3: oracle DA power (no batch at all) by log2FC x n_per\n")
print(round(tapply(ora$DA_power, list(ora$da_log2fc, ora$n_per), median, na.rm = TRUE), 3))
cat("\n== B1 preview: oracle DA FDR by confounding (should be ~0.05 if the truth definition is right)\n")
print(round(tapply(ora$DA_FDR, ora$conf, median, na.rm = TRUE), 3))
cat("\noracle DA FDR by log2FC\n"); print(round(tapply(ora$DA_FDR, ora$da_log2fc, median, na.rm = TRUE), 3))
