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

man <- read.delim(file.path(raw, "mbqc_manifest.tsv"), colClasses = "character")
man$fq <- file.path(raw, "fastq", paste0(man$run, ".fastq.gz"))
man <- man[file.exists(man$fq) & file.size(man$fq) > 5000, ]       # skip empty/failed downloads
say(nrow(man), "runs across", length(unique(man$lab)), "labs")

# Truncate where the per-cycle median quality drops below 25, judged from a few files per lab: the labs used
# different read lengths, so one global truncLen would throw away whole labs (filterAndTrim drops any read
# shorter than truncLen).
trunc_for <- function(files) {
  q <- lapply(head(files, 3), function(f) {
    fq <- yield(FastqSampler(f, 2000)); if (length(fq) == 0) return(NULL)
    m <- as(quality(fq), "matrix"); apply(m, 2, median, na.rm = TRUE)
  })
  q <- q[!vapply(q, is.null, TRUE)]
  if (!length(q)) return(100L)
  n <- min(lengths(q)); med <- rowMeans(vapply(q, function(x) x[seq_len(n)], numeric(n)))
  ok <- which(med >= 25)
  max(100L, if (length(ok)) max(ok) else n)
}

labs <- sort(unique(man$lab))
for (lab in labs) {
  f_rds <- file.path(ck, paste0("seqtab_", lab, ".rds"))
  if (file.exists(f_rds)) { say("lab", lab, "already done, skipping"); next }
  sub <- man[man$lab == lab, ]
  tl <- trunc_for(sub$fq)
  say("lab", lab, ":", nrow(sub), "samples, truncLen =", tl)

  filt <- file.path(raw, "filtered", lab, paste0(sub$run, ".fastq.gz"))
  dir.create(dirname(filt[1]), recursive = TRUE, showWarnings = FALSE)
  st <- filterAndTrim(sub$fq, filt, truncLen = tl, maxN = 0, maxEE = 2, truncQ = 2,
                      rm.phix = TRUE, compress = TRUE, multithread = FALSE, verbose = FALSE)
  keep <- file.exists(filt) & file.size(filt) > 100
  say("  kept", sum(keep), "of", nrow(sub), "samples;",
      round(100 * sum(st[, "reads.out"]) / max(1, sum(st[, "reads.in"]))), "% of reads passed the filter")
  if (sum(keep) < 2) { say("  too few samples, skipping lab"); next }
  filt <- filt[keep]; names(filt) <- sub$run[keep]

  err <- learnErrors(filt, nbases = 5e7, multithread = FALSE, randomize = TRUE, verbose = FALSE)
  dd <- dada(derepFastq(filt, verbose = FALSE), err = err, multithread = FALSE, verbose = FALSE)
  if (!is.list(dd) || inherits(dd, "dada")) dd <- list(dd)      # dada() unwraps a single sample
  names(dd) <- names(filt)
  tab <- removeBimeraDenovo(makeSequenceTable(dd), method = "consensus", multithread = FALSE, verbose = FALSE)
  say("  ", ncol(tab), "ASVs after chimera removal,", sum(tab), "reads")
  saveRDS(tab, f_rds)
}

tabs <- lapply(labs, function(l) { f <- file.path(ck, paste0("seqtab_", l, ".rds")); if (file.exists(f)) readRDS(f) })
names(tabs) <- labs; tabs <- tabs[!vapply(tabs, is.null, TRUE)]
stopifnot("no lab produced a sequence table" = length(tabs) > 0)

# One taxonomy pass over the union of ASVs: loading the 451k-sequence reference is the expensive part.
tax_rds <- file.path(ck, "taxonomy.rds")
seqs <- unique(unlist(lapply(tabs, colnames), use.names = FALSE))
say("assigning taxonomy to", length(seqs), "unique ASVs (this is the slow step)")
tx <- if (file.exists(tax_rds)) readRDS(tax_rds) else {
  t <- assignTaxonomy(seqs, silva, minBoot = 50, multithread = FALSE, verbose = FALSE)
  saveRDS(t, tax_rds); t
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
