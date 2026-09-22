# Collect every score into one long table: dataset, status, runtime, scenario params, metrics.
# Usage: Rscript report/aggregate.R <results dir> <file listing the .score paths> <out.tsv>
a <- commandArgs(TRUE); res <- a[1]
files <- readLines(a[2])
rows <- lapply(files, function(f) {
  ds <- dirname(sub(paste0("^", file.path(res, "out"), "/"), "", f))
  s <- read.delim(f)
  s$dataset <- ds
  s$status <- readLines(sub("\\.score$", ".status", f))[1]
  s$status[s$status == "ok" && ncol(s) <= 3] <- "invalid output"  # ran, but compute.R rejected the output
  b <- sub("^out/", "bench/", sub("\\.score$", ".txt", sub(paste0("^", res, "/"), "", f)))
  b <- file.path(res, b)
  if (file.exists(b)) { bt <- read.delim(b); s$seconds <- bt$s[1]; s$max_rss_mb <- bt$max_rss[1] }
  pj <- file.path(res, "data", ds, "params.json")
  if (file.exists(pj)) s <- cbind(s, as.data.frame(jsonlite::fromJSON(pj)))
  s
})
cols <- unique(unlist(lapply(rows, names)))
out <- do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(out), "rows,", sum(out$status == "ok"), "ok ->", a[3], "\n")
