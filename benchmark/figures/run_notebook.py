#!/usr/bin/env python
"""Execute one figure notebook in place.

The notebooks are the figure code now: each panel is defined and drawn in its own cell, and the
last cell assembles them and writes results/figures/figN.pdf. Running a notebook here keeps the
rendered panels inside the .ipynb, so the file a reader opens already shows every panel without
needing a kernel or the 13 GB of results behind it.

Needs the `ir` kernel (IRkernel in the bench-r environment); see the repository README.

Usage: python figures/run_notebook.py figures/fig1.ipynb
"""
import sys

import nbformat
from nbclient import NotebookClient

path = sys.argv[1]
notebook = nbformat.read(path, as_version=4)
# cwd must be benchmark/ -- the panels read config/ and results/ by relative path
NotebookClient(notebook, timeout=3600, kernel_name="ir",
               resources={"metadata": {"path": "."}}).execute()
nbformat.write(notebook, path)

drawn = sum("image/png" in output.get("data", {})
            for cell in notebook.cells for output in cell.get("outputs", []))
print(f"executed {path}: {len(notebook.cells)} cells, {drawn} figures drawn")
