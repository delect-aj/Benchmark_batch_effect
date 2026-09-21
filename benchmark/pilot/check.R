# Boundary checks on pilot scores (PLAN.md section 7). Usage: Rscript check.R <scores.tsv>
s <- read.delim(commandArgs(TRUE)[1]); rownames(s) <- s$method
print(s[order(s$R2_batch_ait), c("method", "kind", "R2_batch_ait", "R2_pheno_ait", "R2_batch_bc",
                                 "iLISI", "oracle_mantel", "DA_FDR", "DA_power")], row.names = FALSE, digits = 3)
stopifnot("raw, limma and oracle scores required" = all(c("raw", "limma", "oracle") %in% s$method),
          "limma should lower batch R2 vs raw" = s["limma", "R2_batch_ait"] < s["raw", "R2_batch_ait"],
          "oracle should have less batch R2 than raw" = s["oracle", "R2_batch_ait"] < s["raw", "R2_batch_ait"],
          "oracle should match itself" = s["oracle", "oracle_mantel"] > 0.99,
          "raw should be further from oracle than oracle" = s["raw", "oracle_mantel"] < s["oracle", "oracle_mantel"])
cat("check passed\n")
