# MBQC (Microbiome Quality Control, Sinha et al. 2017): 22 specimens sent to 13 handling labs, each doing its
# own DNA extraction and 16S V4 sequencing. The consortium's integrated OTU table is gone (the server that held
# it is offline), so we re-call ASVs ourselves from the raw forward reads with one pipeline -- which is the
# point: every batch here is pure wet-lab difference, with the same specimens repeated across all 13 labs.
#
# Forward reads only (the labs ran different configurations -- 2x126 up to 2x294 -- so reverse reads would not
# overlap consistently). Each lab is denoised on its own, which is also DADA2's rule for separate runs; ASVs
# therefore differ in length between labs and are NOT comparable as sequences, so the shared currency is genus,
# exactly as for the MicrobiomeHD datasets.
#
# Rscript is broken on the compute nodes (old prefix compiled in), so run it through the R wrapper:
#   . ~/renv.sh && $RR -f data/mbqc.R --args <outdir> <raw dir> <silva trainset .fa.gz>
a <- commandArgs(TRUE)
outdir <- a[1]; raw <- a[2]
silva <- if (length(a) > 2) a[3] else "/beegfs/dongbiao/silva_138.2/silva_138.2_dada2_genus.fa.gz"
suppressMessages({library(dada2); library(ShortRead)})
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ck <- file.path(raw, "checkpoints"); dir.create(ck, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# The checkpoints are on shared storage and this is launched by hand on whichever node is free, so guard
# against two builds racing on the same RDS files. The lock is refreshed as each lab finishes; one older than
# three hours is assumed stale and overridden.
lock <- file.path(ck, ".lock")
if (file.exists(lock) && difftime(Sys.time(), file.mtime(lock), units = "hours") < 3)
  stop("another build holds the lock (", readLines(lock, warn = FALSE)[1], "); delete ", lock, " to override")
hold <- function() writeLines(paste(Sys.info()[["nodename"]], Sys.getpid(), format(Sys.time())), lock)
hold(); on.exit(unlink(lock), add = TRUE)

# Keep the whole manifest: the sample metadata at the end is built from it, and a lab that already has a
# checkpoint needs no reads on disk at all -- which is what lets the 13 GB of FastQ be deleted once the
# per-lab sequence tables exist.
man <- read.delim(file.path(raw, "mbqc_manifest.tsv"), colClasses = "character")
man$fq <- file.path(raw, "fastq", paste0(man$run, ".fastq.gz"))
man$have <- file.exists(man$fq) & file.size(man$fq) > 5000         # empty/failed downloads count as absent
say(sum(man$have), "of", nrow(man), "runs have reads on disk, across",
    length(unique(man$lab)), "labs")

# Truncate where the per-cycle median quality drops below 25, but never past the length the reads actually
# have: filterAndTrim discards every read shorter than truncLen, and the labs ran anything from 126 to 294 bp.
# The quality matrix pads short reads with NA, so the trailing cycles' medians come from the few long reads
# and would otherwise pull truncLen past the bulk of the data -- cap it at the 10th percentile of read length.
trunc_for <- function(files) {
  qs <- list(); lens <- integer()
  for (f in head(files, 4)) {
    fq <- tryCatch(yield(FastqSampler(f, 2000)), error = function(e) NULL)
    if (is.null(fq) || length(fq) == 0) next
    qs[[length(qs) + 1L]] <- apply(as(quality(fq), "matrix"), 2, median, na.rm = TRUE)
    lens <- c(lens, width(fq))
  }
  if (!length(qs)) return(c(trunc = 100L, p10 = NA_integer_, med = NA_integer_))
  n <- min(lengths(qs)); med <- rowMeans(vapply(qs, function(x) x[seq_len(n)], numeric(n)))
  ok <- which(med >= 25)
  p10 <- as.integer(quantile(lens, 0.10, names = FALSE))
  c(trunc = max(100L, min(if (length(ok)) max(ok) else n, p10)), p10 = p10, med = as.integer(median(lens)))
}

labs <- sort(unique(man$lab))
for (lab in labs) {
  f_rds <- file.path(ck, paste0("seqtab_", lab, ".rds"))
  if (file.exists(f_rds)) { say("lab", lab, "already done, skipping"); next }
  sub <- man[man$lab == lab & man$have, ]
  if (!nrow(sub)) { say("lab", lab, ": no checkpoint and no reads on disk, skipping"); next }
  tl <- trunc_for(sub$fq)
  say("lab", lab, ":", nrow(sub), "samples, read length median", tl[["med"]], "p10", tl[["p10"]],
      "-> truncLen", tl[["trunc"]])

  filt <- file.path(raw, "filtered", lab, paste0(sub$run, ".fastq.gz"))
  dir.create(dirname(filt[1]), recursive = TRUE, showWarnings = FALSE)
  # One malformed FastQ would otherwise abort the whole build, so filter file by file and skip the bad ones
  st <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i)
    tryCatch(filterAndTrim(sub$fq[i], filt[i], truncLen = tl[["trunc"]], maxN = 0, maxEE = 2, truncQ = 2,
                           rm.phix = TRUE, compress = TRUE, multithread = FALSE, verbose = FALSE),
             error = function(e) { say("   unreadable, skipped:", sub$run[i], "-", conditionMessage(e))
                                   matrix(0L, 1, 2, dimnames = list(sub$run[i], c("reads.in", "reads.out"))) })))
  keep <- file.exists(filt) & file.size(filt) > 100
  rate <- sum(st[, "reads.out"]) / max(1, sum(st[, "reads.in"]))
  say("  kept", sum(keep), "of", nrow(sub), "samples;", round(100 * rate), "% of reads passed the filter")
  # Median quality can sit just above 25 for the whole read and still blow the maxEE budget over 250 cycles,
  # which costs a lab most of its data. Keep maxEE the same everywhere and shorten that lab's reads instead.
  for (try in 1:2) {
    if (rate >= 0.50 || tl[["trunc"]] <= 100L) break
    tl[["trunc"]] <- max(100L, as.integer(tl[["trunc"]] * 0.8))
    say("  retrying", lab, "with truncLen", tl[["trunc"]])
    st <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i)
      tryCatch(filterAndTrim(sub$fq[i], filt[i], truncLen = tl[["trunc"]], maxN = 0, maxEE = 2, truncQ = 2,
                             rm.phix = TRUE, compress = TRUE, multithread = FALSE, verbose = FALSE),
               error = function(e) matrix(0L, 1, 2, dimnames = list(sub$run[i], c("reads.in", "reads.out"))))))
    keep <- file.exists(filt) & file.size(filt) > 100
    rate <- sum(st[, "reads.out"]) / max(1, sum(st[, "reads.in"]))
    say("  kept", sum(keep), "of", nrow(sub), "samples;", round(100 * rate), "% of reads passed the filter")
  }
  if (sum(keep) < 2) { say("  too few samples, skipping lab"); next }
  # A lab that still loses nearly everything is broken; fail loudly instead of checkpointing junk
  if (rate < 0.20) { say("  pass rate too low -- not checkpointing", lab, "; check truncLen"); next }
  filt <- filt[keep]; names(filt) <- sub$run[keep]

  err <- learnErrors(filt, nbases = 5e7, multithread = FALSE, randomize = TRUE, verbose = FALSE)
  dd <- dada(derepFastq(filt, verbose = FALSE), err = err, multithread = FALSE, verbose = FALSE)
  if (!is.list(dd) || inherits(dd, "dada")) dd <- list(dd)      # dada() unwraps a single sample
  names(dd) <- names(filt)
  tab <- removeBimeraDenovo(makeSequenceTable(dd), method = "consensus", multithread = FALSE, verbose = FALSE)
  say("  ", ncol(tab), "ASVs after chimera removal,", sum(tab), "reads")
  saveRDS(tab, f_rds); hold()
}

tabs <- lapply(labs, function(l) { f <- file.path(ck, paste0("seqtab_", l, ".rds")); if (file.exists(f)) readRDS(f) })
names(tabs) <- labs; tabs <- tabs[!vapply(tabs, is.null, TRUE)]
stopifnot("no lab produced a sequence table" = length(tabs) > 0)

# One taxonomy pass over the union of ASVs: loading the 451k-sequence reference is the expensive part.
tax_rds <- file.path(ck, "taxonomy.rds")
seqs <- unique(unlist(lapply(tabs, colnames), use.names = FALSE))
say("assigning taxonomy to", length(seqs), "unique ASVs (this is the slow step)")
cached <- if (file.exists(tax_rds)) readRDS(tax_rds) else NULL
tx <- if (!is.null(cached) && all(seqs %in% rownames(cached))) cached[seqs, , drop = FALSE] else {
  t <- assignTaxonomy(seqs, silva, minBoot = 50, multithread = FALSE, verbose = FALSE)
  saveRDS(t, tax_rds); t          # a lab added later brings new ASVs, so the cache is only good if it covers them
}
ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
colnames(tx) <- ranks
lin <- apply(tx, 1, function(r) { r[is.na(r)] <- "unclassified"; paste(r, collapse = ";") })

# Collapse each lab's ASVs to genus, then stack the labs on the union of genera
per_lab <- lapply(names(tabs), function(l) {
  tb <- tabs[[l]]
  g <- t(rowsum(t(tb), lin[colnames(tb)]))
  rownames(g) <- paste0(l, "_", rownames(tb))
  g
})
taxa <- Reduce(union, lapply(per_lab, colnames))
counts <- do.call(rbind, lapply(per_lab, function(g) {
  m <- matrix(0L, nrow(g), length(taxa), dimnames = list(rownames(g), taxa))
  m[, colnames(g)] <- g; m
}))

run_of <- sub("^[^_]+_", "", rownames(counts))
mi <- man[match(run_of, man$run), ]
sm <- read.delim(file.path(raw, "mbqc_samples.tsv"), colClasses = "character")
si <- sm[match(mi$specimen, sm$specimen), ]
meta <- data.frame(batch = mi$lab,              # handling lab: the only thing that differs between replicates
                   phenotype = mi$specimen,     # 22 specimens, each sent to every lab: the biological signal
                   replicate = mi$specimen,     # same material measured repeatedly, for RUV-style methods
                   run = mi$run, specimen_type = si$specimen_type, health_status = si$health_status,
                   sex = si$sex, age = si$age, artificial = si$artificial,
                   row.names = rownames(counts))

# Same filter for every dataset (PLAN.md 2.3): prevalence >= 10%, drop empty samples
counts <- counts[, colMeans(counts > 0) >= 0.10, drop = FALSE]
ok <- rowSums(counts) > 0
counts <- counts[ok, ]; meta <- meta[ok, ]

tax <- do.call(rbind, strsplit(colnames(counts), ";"))
tax <- data.frame(tax, row.names = colnames(counts)); names(tax) <- ranks

w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(as.data.frame(counts), "counts.tsv"); w(meta, "meta.tsv"); w(tax, "taxonomy.tsv")
print(table(meta$batch))
say("mbqc:", nrow(counts), "samples x", ncol(counts), "genera,",
    nlevels(factor(meta$batch)), "handling labs,", nlevels(factor(meta$phenotype)), "specimens")
