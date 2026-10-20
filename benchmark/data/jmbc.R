# JMBC (Tourlousse et al. 2021, Microbiome 9:95): the Japanese consortium's protocol comparison, rebuilt from
# raw reads because the paper ships nothing but a PDF. Two phase-I experiments, kept apart because their
# batch axes are different things:
#
#   extraction  9 DNA extraction protocols (L-T) applied to a 20-strain cell mock and to feces.
#               batch = protocol, phenotype = material. This is the only dataset in the benchmark where the
#               batch variable is the DNA extraction method, which the literature calls the dominant source
#               of protocol bias (Costea 2017, DEBIAS-M).
#   library     11 library kits x 3 input amounts, all on one DNA mock. Nothing biological varies at all, so
#               it is a null set: whatever structure a method leaves behind is batch, and whatever it invents
#               is a false signal. Calibration only, not a case-control dataset.
#
# Usage: Rscript data/jmbc.R <extraction|library> <outdir> <raw dir>
a <- commandArgs(TRUE); which <- a[1]; outdir <- a[2]; raw <- a[3]
stopifnot("first argument must be extraction or library" = which %in% c("extraction", "library"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

mock <- read.delim("data/jmbc/jmbc_mock.tsv", colClasses = "character")
runs <- read.delim("data/jmbc/jmbc_runs.tsv", colClasses = "character")
long <- read.delim(file.path(raw, "counts_long.tsv"), header = FALSE,
                   col.names = c("run", "assembly", "n"))

runs <- runs[runs$experiment == which & runs$run %in% long$run, ]
long <- long[long$run %in% runs$run, ]
stopifnot("no runs quantified for this experiment" = nrow(runs) > 0)

# reads x strain -> samples x strains, naming columns by species rather than accession
m <- matrix(0L, nrow(runs), nrow(mock), dimnames = list(runs$run, mock$species))
idx <- cbind(match(long$run, runs$run), match(long$assembly, mock$assembly))
ok <- !is.na(idx[, 1]) & !is.na(idx[, 2])
m[idx[ok, ]] <- as.integer(long$n[ok])

meta <- data.frame(batch = runs$protocol,        # extraction protocol L-T, or library kit + input amount
                   phenotype = runs$material,    # cell mock vs feces; constant in the library experiment
                   experiment = runs$experiment, detail = runs$detail,
                   read_pairs_total = runs$read_pairs, row.names = runs$run)

# Known composition, from Table S1. The cell mock omits Megamonas funiformis (0.0%), the DNA mock is even.
pct <- as.numeric(mock[[if (which == "library") "pct_dna_mock" else "pct_cell_mock"]])
names(pct) <- mock$species
mock_rows <- rownames(meta)[meta$phenotype %in% c("cell_mock", "dna_mock")]
truth <- matrix(pct / sum(pct), length(mock_rows), length(pct), byrow = TRUE,
                dimnames = list(mock_rows, names(pct)))

keep <- colSums(m) > 0
m <- m[, keep, drop = FALSE]; truth <- truth[, keep, drop = FALSE]
ok <- rowSums(m) > 0
m <- m[ok, , drop = FALSE]; meta <- meta[ok, , drop = FALSE]
truth <- truth[rownames(truth) %in% rownames(m), , drop = FALSE]

tax <- data.frame(genus = sub(" .*", "", colnames(m)), species = colnames(m),
                  gram = mock$gram[match(colnames(m), mock$species)],
                  genome_size = mock$genome_size[match(colnames(m), mock$species)],
                  row.names = colnames(m))

w <- function(df, f) write.table(data.frame(sample_id = rownames(df), df, check.names = FALSE),
                                 file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
w(as.data.frame(m), "counts.tsv"); w(meta, "meta.tsv"); w(tax, "taxonomy.tsv")
if (nrow(truth)) w(as.data.frame(truth), "truth.tsv")
print(table(meta$batch, meta$phenotype))
cat("jmbc_", which, ": ", nrow(m), " samples x ", ncol(m), " strains, ",
    nlevels(factor(meta$batch)), " protocols, median ", round(median(rowSums(m))), " mapped reads\n", sep = "")
