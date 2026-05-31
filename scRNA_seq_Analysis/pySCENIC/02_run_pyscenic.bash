#!/bin/bash
# ============================================
# PySCENIC Pipeline: GRN + Regulon + AUCell
# ============================================

# Set working directory
cd /data/work/pyscenic/

# ============================================
# Step 1: GRN inference (GRNBoost2)
# ============================================
echo "Step 1: Running GRN inference..."
pyscenic grn \
  --num_workers 8 \
  --output CD8CTL_scenic_grn.tsv \
  --method grnboost2 \
  CD8CTL_scenic.loom \
  hs_hgnc_tfs.txt

# ============================================
# Step 2: Regulon identification (CisTarget)
# ============================================
echo "Step 2: Identifying regulons..."
pyscenic ctx \
  CD8CTL_scenic_grn.tsv \
  hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
  --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname CD8CTL_scenic.loom \
  --mode dask_multiprocessing \
  --output CD8CTL_scenic_ctx.csv \
  --num_workers 8 \
  --mask_dropouts

# ============================================
# Step 3: AUCell calculation
# ============================================
echo "Step 3: Calculating AUCell scores..."
pyscenic aucell \
  CD8CTL_scenic.loom \
  CD8CTL_scenic_ctx.csv \
  --output CD8CTL_scenic_auc.loom \
  --num_workers 8

echo "Done!"
