# Leave-one-study-out cross-study prediction.
#
# What leaks and what does not. A batch-correction method cannot place an unseen study on a common scale
# without seeing data from it -- the batch effect of a study is only estimable from that study's samples. A
# strictly inductive protocol ("fit on the training studies, then transform a study the method has never
# seen") is therefore impossible for essentially every method here, and demanding it would mean benchmarking
# a setting nobody can use. What must not leak is the *labels*: the correction may use the held-out study's
# features, but never its phenotype. This driver enforces exactly that by masking the held-out phenotype to
# "unknown" before the wrapper is called, so methods that take a protected biological covariate receive one
# that carries no information about the held-out samples.
#
# Restricted to datasets whose batch variable is the study and whose phenotype is a disease contrast
# (config `loso:`). Leaving one study out is undefined when batches are laboratories or protocols, and a
# phenotype whose classes are not all present in the training studies cannot be predicted at all -- in
# seed_16s, for instance, no study contains all five host plants.
#
# Usage: Rscript metrics/loso.R <data dir> <exe> <wrapper script> <out.tsv> [timeout seconds]
suppressMessages({library(glmnet); library(randomForest); library(pROC)})
a <- commandArgs(TRUE)
datadir <- a[1]; exe <- a[2]; script <- a[3]; out <- a[4]
tmo <- if (length(a) > 4) a[5] else "3600"
ds <- basename(datadir); mth <- sub("\\.(R|py)$", "", basename(script))

counts <- as.matrix(read.delim(file.path(datadir, "counts.tsv"), row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(datadir, "meta.tsv"), row.names = 1, check.names = FALSE)
stopifnot(identical(rownames(counts), rownames(meta)))
y <- factor(meta$phenotype)
stopifnot("LOSO needs a binary phenotype" = nlevels(y) == 2)
studies <- sort(unique(meta$batch))

work <- file.path(tempdir(), paste0("loso_", ds, "_", mth))
dir.create(work, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(datadir, "counts.tsv"), file.path(work, "counts.tsv"), overwrite = TRUE)
for (s in c("taxonomy.tsv", "tree.nwk"))                       # side inputs wrappers may look for
  if (file.exists(file.path(datadir, s))) file.copy(file.path(datadir, s), file.path(work, s), overwrite = TRUE)

clr <- function(x, pseudo = 0.5) { l <- log(x + pseudo); l - rowMeans(l) }

# Corrected tables come back on different scales; put them all on one the classifiers can use
as_features <- function(x, kind) if (kind %in% c("counts", "relabund")) clr(x) else x

auc_of <- function(truth, score)
  tryCatch(as.numeric(pROC::auc(pROC::roc(truth, score, quiet = TRUE, direction = "<"))),
           error = function(e) NA_real_)

rows <- list()
for (s in studies) {
  test <- meta$batch == s
  note <- if (nlevels(droplevels(y[test])) < 2) "held-out study has one class" else
          if (nlevels(droplevels(y[!test])) < 2) "training studies have one class" else NA
  if (!is.na(note)) {
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, method = mth, held_out = s, classifier = NA,
                                           auroc = NA, n_train = sum(!test), n_test = sum(test),
                                           status = note)
    next
  }

  # The correction sees the held-out study's counts but not its labels
  m2 <- meta
  m2$phenotype <- as.character(m2$phenotype)
  m2$phenotype[test] <- "unknown"
  write.table(data.frame(sample_id = rownames(m2), m2, check.names = FALSE),
              file.path(work, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  o <- file.path(work, "out.tsv"); unlink(c(o, paste0(o, ".json")))
  st <- system2("timeout", c(tmo, exe, script, file.path(work, "counts.tsv"), file.path(work, "meta.tsv"), o),
                stdout = FALSE, stderr = FALSE)

  if (st != 0 || !file.exists(o)) {
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, method = mth, held_out = s, classifier = NA,
                                           auroc = NA, n_train = sum(!test), n_test = sum(test),
                                           status = paste0("correction failed (exit ", st, ")"))
    next
  }
  x <- as.matrix(read.delim(o, row.names = 1, check.names = FALSE))
  kind <- sub('.*"kind"[^"]*"([^"]+)".*', "\\1", readLines(paste0(o, ".json"), warn = FALSE)[1])
  x <- as_features(x[rownames(meta), , drop = FALSE], kind)
  x <- x[, apply(x[!test, , drop = FALSE], 2, sd) > 0, drop = FALSE]   # constant in training -> uninformative

  yt <- droplevels(y[!test]); ys <- droplevels(y[test])
  fits <- list(
    L1_logistic = tryCatch({
      cv <- cv.glmnet(x[!test, ], yt, family = "binomial", alpha = 1, nfolds = 5)
      as.numeric(predict(cv, x[test, , drop = FALSE], s = "lambda.1se", type = "response"))
    }, error = function(e) NULL),
    random_forest = tryCatch({
      rf <- randomForest(x[!test, ], yt, ntree = 500)
      predict(rf, x[test, , drop = FALSE], type = "prob")[, levels(yt)[2]]
    }, error = function(e) NULL))

  for (cl in names(fits)) {
    p <- fits[[cl]]
    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds, method = mth, held_out = s, classifier = cl,
      auroc = if (is.null(p)) NA_real_ else auc_of(ys, p),
      n_train = sum(!test), n_test = sum(test),
      status = if (is.null(p)) "classifier failed" else "ok")
  }
}

res <- do.call(rbind, rows)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
ok <- res[res$status == "ok" & !is.na(res$auroc), ]
cat(sprintf("%s / %s: %d folds, %d scored, mean auROC %s\n", ds, mth, length(studies), nrow(ok),
            if (nrow(ok)) sprintf("%.3f", mean(ok$auroc)) else "-"))
