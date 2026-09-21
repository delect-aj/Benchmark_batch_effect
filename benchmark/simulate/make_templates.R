# Simulation templates (samples x taxa counts). Usage: Rscript simulate/make_templates.R <16s|mgx> <out.tsv> [crc_dir]
#   16s: MIDASim's bundled IBD 16S table (146 samples x 614 taxa)
#   mgx: control samples of the CRC metagenome dataset (data/crc_mgx.R output)
a <- commandArgs(TRUE)
tab <- if (a[1] == "16s") {
  data(count.ibd, package = "MIDASim", envir = environment()); count.ibd
} else {
  x <- as.matrix(read.delim(file.path(a[3], "counts.tsv"), row.names = 1, check.names = FALSE))
  m <- read.delim(file.path(a[3], "meta.tsv"), row.names = 1)
  x[m$phenotype == 0, ]
}
tab <- tab[rowSums(tab) > 0, colSums(tab > 0) > 0]
dir.create(dirname(a[2]), recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(sample_id = rownames(tab), tab, check.names = FALSE), a[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(a[1], "template:", nrow(tab), "samples x", ncol(tab), "taxa\n")
