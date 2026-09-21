#!/bin/bash
# Build both benchmark envs. Source of truth for dependencies (replaces env yamls: `conda env create -f`
# cannot override the user's configured channels, and mixing r/defaults with bioconda stalls the solver).
# Built for CentOS 7 / glibc 2.17 nodes: everything that needs compiling comes as a conda-forge binary,
# because the system cmake (2.8) and system libs are too old for source builds.
# Usage: bash envs/build.sh <dir with GitHub clones: ConQuR ruvIIInb metacal MetaDICT DEBIAS-M>
set -euo pipefail
src=$1
source ~/software/miniconda3/etc/profile.d/conda.sh
CH="--override-channels -c conda-forge -c bioconda"

conda env list | grep -q '^bench-r ' || conda create -y -n bench-r $CH r-base=4.4
# r-cvxr<1: ANCOMBC 2.8 calls CVXR::solve, which CVXR 1.x no longer exports
conda install -y -n bench-r $CH \
  compilers make cmake pkg-config libuv libxml2 libcurl nlopt cairo fontconfig freetype harfbuzz fribidi \
  r-biocmanager r-remotes r-jsonlite r-vegan r-harmony r-quantreg r-glmnet r-dplyr r-doparallel r-gplots \
  r-ade4 r-compositions r-randomforest r-rocr r-ape r-fastdummies r-mass r-pscl r-ggpubr r-rmarkdown \
  r-fs r-xml2 r-systemfonts r-nloptr r-haven r-httpuv r-rcurl r-cairo r-rcppparallel r-osqp r-igraph \
  bioconductor-limma bioconductor-sva bioconductor-batchelor bioconductor-singlecellexperiment \
  bioconductor-edger bioconductor-biomformat \
  r-matrix r-lme4 r-tidyverse r-ragg gsl r-gsl r-energy \
  bioconductor-ancombc bioconductor-treesummarizedexperiment "r-cvxr<1" \
  bioconductor-scater bioconductor-singler   # ruvIIInb imports
# Separate step: bioconda's curatedMetagenomicData 3.14 fails its post-link against current rbiom (unifrac no
# longer exported), and a failed post-link rolls back the whole install. Only real-data prep needs it.
# rbiom 1.x binaries link against pre-2021 TBB (tbb::task), hence the tbb pin.
conda install -y -n bench-r $CH bioconductor-curatedmetagenomicdata "r-rbiom<2" "tbb<2021" || echo "WARN: curatedMetagenomicData not installed"
conda run -n bench-r Rscript envs/install_r_extras.R "$src"
echo R_ENV_DONE

conda env list | grep -q '^bench-py ' && conda env remove -y -n bench-py
conda create -y -n bench-py $CH python=3.11 numpy=1.26.4 pandas scikit-learn h5py anndata scvi-tools \
  'pytorch=*=cpu*' lightning pip snakemake-minimal
conda run -n bench-py pip install --no-deps lightning-bolts "$src/DEBIAS-M"
conda run -n bench-py python -c "import debiasm, scvi, anndata; print('python imports ok')"
echo PY_ENV_DONE
