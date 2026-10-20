import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_input, write_output

# debiasm 0.0.2's multi-class path (multiclass.py:93) still passes the pytorch-lightning <1.7 argument names,
# while its own log-additive twin twenty lines below already uses the current ones -- the author renamed one
# function and not the other. Translate them here rather than pinning lightning back, which would disturb the
# binary path that works as is. Nothing about the method changes: all three switches mean the same thing
# ("no checkpointing", "no model summary", "no progress bar") under either spelling.
import pytorch_lightning as pl

_Trainer = pl.Trainer
_RENAMED = (("checkpoint_callback", "enable_checkpointing"),
            ("weights_summary", "enable_model_summary"),
            ("progress_bar_refresh_rate", "enable_progress_bar"))


def _trainer(*args, **kw):
    for old, new in _RENAMED:
        if old in kw:
            kw[new] = bool(kw.pop(old))   # False / None / 0 all mean "off"
    return _Trainer(*args, **kw)


pl.Trainer = _trainer

from debiasm import DebiasMClassifier

counts, meta, out = read_input()
# DEBIAS-M input: first column = integer batch index (from 0), then read counts
batch_idx = pd.factorize(meta["batch"])[0]
X = np.hstack([batch_idx[:, None], counts.values])
# factorize, not astype(int): phenotypes are labels, and some datasets use names rather than 0/1
y = pd.factorize(meta["phenotype"])[0]

model = DebiasMClassifier(x_val=X)  # transductive: all samples used to align batches
model.fit(X, y)
x = pd.DataFrame(model.transform(X), index=counts.index, columns=counts.columns)
write_output(x, "relabund", out)
