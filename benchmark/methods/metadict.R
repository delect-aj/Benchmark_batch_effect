source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "_common.R"))
d <- read_input()
# Taxon similarity: tree.nwk > taxonomy.tsv (rows = taxa, columns = ranks) > none
args <- list(count = t(d$counts), meta = d$meta[, c("batch", "phenotype")], covariates = "phenotype", verbose = FALSE)
if (!is.null(f <- side_file(d, "tree.nwk"))) {
  args$tree <- ape::read.tree(f)
} else if (!is.null(f <- side_file(d, "taxonomy.tsv"))) {
  args$taxonomy <- read.delim(f, row.names = 1)
  args$tax_level <- tail(colnames(args$taxonomy), 1)
} else {
  # ponytail: no taxonomy -> flat distances, MetaDICT loses its smoothness prior; give real data a taxonomy.tsv
  p <- ncol(d$counts)
  args$distance_matrix <- matrix(1, p, p, dimnames = list(colnames(d$counts), colnames(d$counts))) - diag(p)
}
x <- t(as.matrix(do.call(MetaDICT::MetaDICT, args)$count))
dimnames(x) <- dimnames(d$counts)
write_output(x, "counts", d$out)
