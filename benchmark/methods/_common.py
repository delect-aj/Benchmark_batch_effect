"""Shared I/O for every Python wrapper. Contract: see benchmark/README.md"""
import json
import sys

import pandas as pd

KINDS = {"counts", "relabund", "clr", "log", "percentile", "embedding"}


def read_input():
    if len(sys.argv) != 4:
        sys.exit("usage: python <wrapper>.py counts.tsv meta.tsv out.tsv")
    counts = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
    meta = pd.read_csv(sys.argv[2], sep="\t", index_col=0)
    if not counts.index.equals(meta.index):
        sys.exit("counts and meta sample_id differ / out of order")
    if not {"batch", "phenotype"} <= set(meta.columns):
        sys.exit("meta needs columns: batch, phenotype")
    return counts, meta, sys.argv[3]


def write_output(x: pd.DataFrame, kind: str, out: str):
    if kind not in KINDS:
        raise ValueError(f"unknown kind: {kind}")
    if kind == "embedding":
        x.columns = [f"dim{i + 1}" for i in range(x.shape[1])]
    x.index.name = "sample_id"
    x.to_csv(out, sep="\t")
    with open(out + ".json", "w") as f:
        json.dump({"kind": kind}, f)
