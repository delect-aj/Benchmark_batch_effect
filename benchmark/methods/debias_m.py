import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_input, write_output
from debiasm import DebiasMClassifier

counts, meta, out = read_input()
# DEBIAS-M input: first column = integer batch index (from 0), then read counts
batch_idx = pd.factorize(meta["batch"])[0]
X = np.hstack([batch_idx[:, None], counts.values])
y = meta["phenotype"].astype(int).values

model = DebiasMClassifier(x_val=X)  # transductive: all samples used to align batches
model.fit(X, y)
x = pd.DataFrame(model.transform(X), index=counts.index, columns=counts.columns)
write_output(x, "relabund", out)
