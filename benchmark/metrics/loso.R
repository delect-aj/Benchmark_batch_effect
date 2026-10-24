# Leave-one-study-out cross-study prediction, scored on the table the method already produced.
#
# Protocol. The correction is run once, on the whole dataset, by the ordinary `run` rule; this script only
# re-splits the corrected table by study, trains a classifier on the remaining studies and predicts the one
# held out. Correction therefore saw every label, including the held-out study's, so these auROCs are an
# optimistic **upper bound** on what cross-study prediction can deliver, not an estimate of it. They are
# labelled as such throughout.
#
# Two things make the bound useful rather than merely cheap. First, `raw` involves no correction and so no
# leakage, which makes it an unbiased reference: if a corrected table does not beat raw even with the
# leakage working in its favour, correction did not help, and that conclusion is safe. Second, the auROC is
# computed on exactly the table that the R2, ASW and LISI metrics are computed on, so the rank comparison
# between "batch removed" and "predicts across studies" is like for like.
#
# Restricted to datasets whose batch is the study and whose phenotype is a disease contrast (config `loso:`).
# Leaving one study out is undefined when batches are laboratories or protocols, and a phenotype whose
# classes are absent from the training studies cannot be predicted at all -- in seed_16s no study contains
# all five host plants.
#
# Usage: Rscript metrics/loso.R <data dir> <corrected table> <out.tsv>
suppressMessages({library(glmnet); library(randomForest); library(pROC)})
a <- commandArgs(TRUE)
datadir <- a[1]; corrected <- a[2]; out <- a[3]
ds <- basename(datadir); mth <- sub("\\.tsv$", "", basename(corrected))

meta <- read.delim(file.path(datadir, "meta.tsv"), row.names = 1, check.names = FALSE)
y <- factor(meta$phenotype)
stopifnot("LOSO needs a binary phenotype" = nlevels(y) == 2)

x <- as.matrix(read.delim(corrected, row.names = 1, check.names = FALSE))
stopifnot("corrected table does not cover the sample set" = all(rownames(meta) %in% rownames(x)))
x <- x[rownames(meta), , drop = FALSE]
kind <- sub('.*"kind"[^"]*"([^"]+)".*', "\\1", readLines(paste0(corrected, ".json"), warn = FALSE)[1])
if (kind %in% c("counts", "relabund")) { l <- log(x + 0.5); x <- l - rowMeans(l) }   # onto a scale a classifier can use

auc_of <- function(truth, score)
  tryCatch(as.numeric(pROC::auc(pROC::roc(truth, score, quiet = TRUE, direction = "<"))),
           error = function(e) NA_real_)
row <- function(s, cl, v, test, note) data.frame(
  dataset = ds, method = mth, held_out = s, classifier = cl, auroc = v,
  n_train = sum(!test), n_test = sum(test), n_pos_test = sum(y[test] == levels(y)[2]), status = note)

rows <- list()
for (s in sort(unique(meta$batch))) {
  test <- meta$batch == s
  note <- if (nlevels(droplevels(y[test])) < 2) "held-out study has one class" else
          if (nlevels(droplevels(y[!test])) < 2) "training studies have one class" else NA
  if (!is.na(note)) { rows[[length(rows) + 1]] <- row(s, NA, NA_real_, test, note); next }

  xi <- x[, apply(x[!test, , drop = FALSE], 2, sd) > 0, drop = FALSE]   # constant in training -> uninformative
  yt <- droplevels(y[!test]); ys <- droplevels(y[test])
  fits <- list(
    # lambda.min, not lambda.1se: the one-standard-error rule is tuned for parsimony and on these data it
    # selects no features at all, leaving an intercept-only model whose predictions are constant. That scores
    # exactly 0.5 for every method and would erase the comparison. A fold where even lambda.min keeps nothing
    # is reported as degenerate rather than as chance performance -- no model was fitted, so there is no
    # performance to report.
    L1_logistic = tryCatch({
      cv <- cv.glmnet(xi[!test, ], yt, family = "binomial", alpha = 1, nfolds = 5)
      if (sum(coef(cv, s = "lambda.min")[-1] != 0) == 0) NA
      else as.numeric(predict(cv, xi[test, , drop = FALSE], s = "lambda.min", type = "response"))
    }, error = function(e) NULL),
    random_forest = tryCatch({
      rf <- randomForest(xi[!test, ], yt, ntree = 500)
      predict(rf, xi[test, , drop = FALSE], type = "prob")[, levels(yt)[2]]
    }, error = function(e) NULL))

  for (cl in names(fits)) {
    p <- fits[[cl]]
    st <- if (is.null(p)) "classifier failed" else if (length(p) == 1 && is.na(p)) "degenerate: no features selected" else "ok"
    rows[[length(rows) + 1]] <- row(s, cl, if (st == "ok") auc_of(ys, p) else NA_real_, test, st)
  }
}

res <- do.call(rbind, rows)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
ok <- res[res$status == "ok" & !is.na(res$auroc), ]
cat(sprintf("%s / %s: %d folds, %d scored, mean auROC %s\n", ds, mth,
            length(unique(res$held_out)), nrow(ok),
            if (nrow(ok)) sprintf("%.3f", mean(ok$auroc)) else "-"))
