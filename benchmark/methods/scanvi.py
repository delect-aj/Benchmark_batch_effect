import os
import sys

import anndata as ad
import pandas as pd
import scvi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_input, write_output

counts, meta, out = read_input()
scvi.settings.seed = 0
adata = ad.AnnData(counts.values.astype("float32"), obs=meta[["batch", "phenotype"]].astype(str))
scvi.model.SCVI.setup_anndata(adata, batch_key="batch", labels_key="phenotype")
vae = scvi.model.SCVI(adata)
vae.train()
lvae = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="Unknown")  # all samples labelled, as in the MetaDICT benchmark
lvae.train()
write_output(pd.DataFrame(lvae.get_latent_representation(), index=counts.index), "embedding", out)
