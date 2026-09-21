# Install R packages not provided by envs/r.yaml. Idempotent: skips what is already installed.
# Called by envs/build.sh. Usage (inside bench-r): Rscript envs/install_r_extras.R <dir with GitHub clones>
# GitHub is not reachable from compute nodes, so GitHub packages are cloned elsewhere and installed from disk.
src <- commandArgs(TRUE)[1]
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
need <- function(p) !requireNamespace(p, quietly = TRUE)

cran <- c("cqrReg", "GUniFrac", "MIDASim")
bioc <- c("MMUPHin", "PLSDAbatch", "MetaDICT", "Maaslin2", "phyloseq")
local <- c(ConQuR = "ConQuR", ruvIIInb = "ruvIIInb", metacal = "metacal",
           MetaDICT = "MetaDICT")  # MetaDICT: Bioconductor first, GitHub clone if this R is too old

for (p in cran[sapply(cran, need)]) install.packages(p)
todo <- bioc[sapply(bioc, need)]
if (length(todo)) BiocManager::install(todo, update = FALSE, ask = FALSE)
# BiocManager only warns on failure, so anything still missing falls through to the local clone
for (p in names(local)[sapply(names(local), need)])
  remotes::install_local(file.path(src, local[[p]]), dependencies = NA, upgrade = "never")  # NA: no Suggests

# ruvIII.nb calls DescTools::Winsorize(probs = ...), an argument removed in DescTools 0.99.50
if (packageVersion("DescTools") >= "0.99.50") remotes::install_version("DescTools", "0.99.49", upgrade = "never")

missing <- c(cran, bioc, names(local))[sapply(c(cran, bioc, names(local)), need)]
if (length(missing)) stop("still missing: ", paste(missing, collapse = ", "))
cat("all R packages installed\n")
