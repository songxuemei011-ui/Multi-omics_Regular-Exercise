#!/usr/bin/env python
# ============================================
# Prepare .loom file for PySCENIC
# ============================================

import anndata
import loompy as lp
import numpy as np

# Input/Output paths
INPUT_H5AD = "/data/work/pyscenic/CD8CTL.h5ad"
OUTPUT_LOOM = "/data/work/pyscenic/CD8CTL_scenic.loom"

# Load data
print("Loading:", INPUT_H5AD)
adata = anndata.read_h5ad(INPUT_H5AD)
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

# Prepare row attributes (genes)
row_attrs = {
    "Gene": np.array(adata.var.index, dtype=str).astype('S')
}

# Prepare column attributes (cells)
col_attrs = {
    "CellID": np.array(adata.obs.index, dtype=str).astype('S'),
    "nGene": np.array(np.sum(adata.X > 0, axis=1)).flatten(),
    "nUMI": np.array(np.sum(adata.X, axis=1)).flatten(),
}

# Create loom file
print("Creating loom:", OUTPUT_LOOM)
lp.create(OUTPUT_LOOM, adata.X.transpose(), row_attrs, col_attrs)
print("Done!")

# Verify
with lp.connect(OUTPUT_LOOM) as ds:
    print(f"Loom file created: {ds.shape[0]} genes, {ds.shape[1]} cells")
    print(f"Example genes: {ds.ra.Gene[:5]}")
