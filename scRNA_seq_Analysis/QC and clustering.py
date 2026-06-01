"""
Quality Control and Gene Filtering for Single-cell RNA-seq Data
Data source:
    - Raw expression matrix downloaded from CIMA database (https://cima.cngb.org/)
    - The matrix was already pre-processed and cell-type annotated by CIMA
    - Sample list and metadata used in this study are provided in Supplementary Table 1
This script performs:
    - Cell filtering (genes, UMIs, mitochondrial percentage)
    - Doublet removal (Scrublet)
    - Gene filtering (HB, ncRNA, AC/AL pseudogenes)
    - Normalization and log transformation
    - HVG selection and batch correction (Harmony)
    - QC plots (violin + boxplot by group)

Prerequisites:
    - adata object must be loaded
    - adata.obs must contain: exercise_group

Output:
    - qc_plots.pdf (violin + boxplot for QC metrics)
    - adata_processed.h5ad
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scrublet as scr
import re
import warnings
warnings.filterwarnings('ignore')

# ============================================
# 0. Set up plotting style
# ============================================
sc.settings.set_figure_params(dpi=100, figsize=(8, 6))
COLOR_SED = "#7EBFC9"
COLOR_EX = "#BCDF7A"

# ============================================
# 1. Quality control (cell filtering)
# ============================================
print("=" * 50)
print("Step 1: Quality Control")
print("=" * 50)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

print(f"Before filtering: {adata.n_obs} cells, {adata.n_vars} genes")

# Filter cells
sc.pp.filter_cells(adata, min_genes=500, max_genes=6000)
adata = adata[adata.obs['total_counts'] >= 1000, :]      # MIN_UMI
adata = adata[adata.obs['total_counts'] <= 25000, :]     # MAX_UMI
adata = adata[adata.obs['pct_counts_mt'] < 10, :]
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")

# ============================================
# 2. QC plots (before gene filtering)
# ============================================
print("\n" + "=" * 50)
print("Step 2: QC Plots")
print("=" * 50)

# Prepare data for plotting
qc_df = adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'exercise_group']].copy()
qc_df = qc_df.dropna()

# Create figure with subplots
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# Violin plot for n_genes
sns.violinplot(x='exercise_group', y='n_genes_by_counts', data=qc_df, 
               palette=[COLOR_SED, COLOR_EX], ax=axes[0], cut=0, inner=None)
sns.boxplot(x='exercise_group', y='n_genes_by_counts', data=qc_df, 
            width=0.2, ax=axes[0], color='white', boxprops=dict(alpha=0.5))
axes[0].set_title('Genes per cell')
axes[0].set_xlabel('')
axes[0].set_ylabel('Number of genes')

# Violin plot for UMIs
sns.violinplot(x='exercise_group', y='total_counts', data=qc_df, 
               palette=[COLOR_SED, COLOR_EX], ax=axes[1], cut=0, inner=None)
sns.boxplot(x='exercise_group', y='total_counts', data=qc_df, 
            width=0.2, ax=axes[1], color='white', boxprops=dict(alpha=0.5))
axes[1].set_title('UMIs per cell')
axes[1].set_xlabel('')
axes[1].set_ylabel('Number of UMIs')

# Violin plot for mitochondrial percentage
sns.violinplot(x='exercise_group', y='pct_counts_mt', data=qc_df, 
               palette=[COLOR_SED, COLOR_EX], ax=axes[2], cut=0, inner=None)
sns.boxplot(x='exercise_group', y='pct_counts_mt', data=qc_df, 
            width=0.2, ax=axes[2], color='white', boxprops=dict(alpha=0.5))
axes[2].set_title('Mitochondrial %')
axes[2].set_xlabel('')
axes[2].set_ylabel('Percentage')

plt.tight_layout()
plt.savefig('qc_metrics_by_group.pdf', dpi=300, bbox_inches='tight')
plt.savefig('qc_metrics_by_group.png', dpi=300, bbox_inches='tight')
plt.show()
print("QC plots saved")

# ============================================
# 3. Doublet removal (Scrublet)
# ============================================
print("\n" + "=" * 50)
print("Step 3: Doublet Removal")
print("=" * 50)

raw_counts = adata.X.copy()
scrub = scr.Scrublet(raw_counts, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets
adata = adata[~predicted_doublets, :].copy()
print(f"After doublet removal: {adata.n_obs} cells")

# ============================================
# 4. Remove unwanted genes (HB, ncRNA, AC/AL)
# ============================================
print("\n" + "=" * 50)
print("Step 4: Gene Filtering")
print("=" * 50)

# === FIX: Mark genes before removing ===
# Mark hemoglobin genes
adata.var['hb'] = adata.var_names.str.contains('^HB[AB]', regex=True, case=False)

# Mark ncRNA genes
adata.var['ncRNA'] = adata.var_names.str.contains('^MIR|^SNOR|^LINC|^RN', regex=True, case=False)

# Remove hemoglobin genes
hb_genes = adata.var.index[adata.var['hb'] == True].tolist()
if hb_genes:
    adata = adata[:, ~adata.var.index.isin(hb_genes)].copy()
    print(f"Removed {len(hb_genes)} hemoglobin genes")
else:
    print("No hemoglobin genes found")

# Remove ncRNA genes
ncrna_genes = adata.var.index[adata.var['ncRNA'] == True].tolist()
if ncrna_genes:
    adata = adata[:, ~adata.var.index.isin(ncrna_genes)].copy()
    print(f"Removed {len(ncrna_genes)} ncRNA genes")
else:
    print("No ncRNA genes found")

# Remove AC/AL pseudogenes (improved regex)
ac_al_genes = [gene for gene in adata.var_names 
               if re.match(r"^AC\d{6}\.\d+$", gene) or 
                  re.match(r"^AL\d{6}\.\d+$", gene)]
if ac_al_genes:
    adata = adata[:, ~adata.var_names.isin(ac_al_genes)].copy()
    print(f"Removed {len(ac_al_genes)} AC/AL pseudogenes")
else:
    print("No AC/AL pseudogenes found")

print(f"Final gene count: {adata.n_vars}")

# ============================================
# 5. Normalization and log transformation
# ============================================
print("\n" + "=" * 50)
print("Step 5: Normalization")
print("=" * 50)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("Data normalized and log-transformed")

# ============================================
# 6. Highly variable genes
# ============================================
print("\n" + "=" * 50)
print("Step 6: HVG Selection")
print("=" * 50)

sc.pp.highly_variable_genes(adata, n_top_genes=3000)
print(f"Number of HVGs: {sum(adata.var.highly_variable)}")

# Plot HVG
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig('hvg_plot.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("HVG plot saved")

# ============================================
# 7. Scaling and PCA
# ============================================
print("\n" + "=" * 50)
print("Step 7: Scaling and PCA")
print("=" * 50)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
print("PCA completed")

# Plot PCA variance
sc.pl.pca_variance_ratio(adata, n_pcs=30, show=False)
plt.savefig('pca_variance.pdf', dpi=300, bbox_inches='tight')
plt.close()
print("PCA variance plot saved")

# ============================================
# 8. Batch correction (Harmony)
# ============================================
print("\n" + "=" * 50)
print("Step 8: Batch Correction")
print("=" * 50)

if 'batch' in adata.obs.columns:
    sc.external.pp.harmony_integrate(adata, key='batch')
    print("Batch correction completed with Harmony")
else:
    print("No batch column found, skipping Harmony")

# ============================================
# 9. Save processed data
# ============================================
print("\n" + "=" * 50)
print("Step 9: Saving Data")
print("=" * 50)

adata.write("adata_processed.h5ad", compression='gzip')
print("Processed data saved to: adata_processed.h5ad")
print("\nDone!")
# ============================================
# 10.  Cell Proportion Analysis
# ============================================
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')
# Prepare data
meta = adata.obs[['sample', 'celltype', 'exercise_group', 'sex', 'BMI']].copy()
meta = meta.dropna()

# Calculate cell type proportions per sample
prop_list = []
for sample in meta['sample'].unique():
    sample_data = meta[meta['sample'] == sample]
    total = len(sample_data)
    group = sample_data['exercise_group'].iloc[0]
    sex = sample_data['sex'].iloc[0]
    bmi = sample_data['BMI'].iloc[0]
    for ct in meta['celltype'].unique():
        count = len(sample_data[sample_data['celltype'] == ct])
        prop = count / total * 100
        prop_list.append([sample, ct, group, sex, bmi, prop])

df = pd.DataFrame(prop_list, columns=['sample', 'cell_type', 'group', 'sex', 'BMI', 'proportion'])

# Statistical test
cell_types = df['cell_type'].unique()
pvals = {}
for ct in cell_types:
    df_ct = df[df['cell_type'] == ct].copy()
    if len(df_ct) > 10:
        try:
            df_ct['group_num'] = (df_ct['group'] == 'EX').astype(int)
            model = smf.ols('proportion ~ group_num + C(sex) + BMI', data=df_ct).fit()
            pvals[ct] = model.pvalues['group_num']
        except:
            pvals[ct] = np.nan
    else:
        pvals[ct] = np.nan

# BH correction
valid_ct = [ct for ct in cell_types if not np.isnan(pvals[ct])]
valid_p = [pvals[ct] for ct in valid_ct]
if len(valid_p) > 0:
    _, padj, _, _ = multipletests(valid_p, method='fdr_bh')
    padj_dict = dict(zip(valid_ct, padj))
else:
    padj_dict = {}

def sig_label(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'

# Plot
COLOR_SED = "#7EBFC9"
COLOR_EX = "#BCDF7A"

fig, ax = plt.subplots(figsize=(12, 6))
n_ct = len(cell_types)
x_pos = np.arange(n_ct)
width = 0.35

for i, ct in enumerate(cell_types):
    df_ct = df[df['cell_type'] == ct]
    sed_vals = df_ct[df_ct['group'] == 'SED']['proportion'].values
    ex_vals = df_ct[df_ct['group'] == 'EX']['proportion'].values
    
    # SED
    if len(sed_vals) > 0:
        q1, med, q3 = np.percentile(sed_vals, [25, 50, 75])
        iqr = q3 - q1
        lower = max(np.min(sed_vals), q1 - 1.5 * iqr)
        upper = min(np.max(sed_vals), q3 + 1.5 * iqr)
        rect = plt.Rectangle((x_pos[i] - width - 0.05, q1), width*0.8, q3 - q1,
                              facecolor=COLOR_SED, edgecolor='black', linewidth=1)
        ax.add_patch(rect)
        ax.hlines(med, x_pos[i] - width - 0.05, x_pos[i] - width - 0.05 + width*0.8, 
                  color='black', linewidth=1.5)
        ax.vlines(x_pos[i] - width - 0.05 + width*0.4, lower, q1, color='black', linewidth=1)
        ax.vlines(x_pos[i] - width - 0.05 + width*0.4, q3, upper, color='black', linewidth=1)
        np.random.seed(42)
        jitter = np.random.uniform(-0.08, 0.08, len(sed_vals))
        ax.scatter(x_pos[i] - width - 0.05 + width*0.4 + jitter, sed_vals,
                   color=COLOR_SED, edgecolors='black', linewidths=0.5, s=20, zorder=5, alpha=0.7)
    
    # EX
    if len(ex_vals) > 0:
        q1, med, q3 = np.percentile(ex_vals, [25, 50, 75])
        iqr = q3 - q1
        lower = max(np.min(ex_vals), q1 - 1.5 * iqr)
        upper = min(np.max(ex_vals), q3 + 1.5 * iqr)
        rect = plt.Rectangle((x_pos[i] + 0.05, q1), width*0.8, q3 - q1,
                              facecolor=COLOR_EX, edgecolor='black', linewidth=1)
        ax.add_patch(rect)
        ax.hlines(med, x_pos[i] + 0.05, x_pos[i] + 0.05 + width*0.8, 
                  color='black', linewidth=1.5)
        ax.vlines(x_pos[i] + 0.05 + width*0.4, lower, q1, color='black', linewidth=1)
        ax.vlines(x_pos[i] + 0.05 + width*0.4, q3, upper, color='black', linewidth=1)
        np.random.seed(42)
        jitter = np.random.uniform(-0.08, 0.08, len(ex_vals))
        ax.scatter(x_pos[i] + 0.05 + width*0.4 + jitter, ex_vals,
                   color=COLOR_EX, edgecolors='black', linewidths=0.5, s=20, zorder=5, alpha=0.7)
    
    # Significance
    p_adj = padj_dict.get(ct, np.nan)
    if not np.isnan(p_adj):
        label = sig_label(p_adj)
        y_max = max(np.max(sed_vals) if len(sed_vals) > 0 else 0,
                    np.max(ex_vals) if len(ex_vals) > 0 else 0)
        ax.text(x_pos[i], y_max + 2, label, ha='center', va='bottom', fontsize=11,
                fontweight='bold' if label != 'ns' else 'normal')

ax.set_xticks(x_pos)
ax.set_xticklabels(cell_types, rotation=45, ha='right', fontsize=10)
ax.set_ylabel('Percentage of total cells (%)', fontsize=12)
ax.set_xlim(x_pos[0] - 0.6, x_pos[-1] + 0.6)
ax.spines[['top', 'right']].set_visible(False)

legend_handles = [
    Patch(facecolor=COLOR_SED, edgecolor='black', label='SED'),
    Patch(facecolor=COLOR_EX, edgecolor='black', label='EX'),
]
ax.legend(handles=legend_handles, title='Group', loc='upper right', frameon=True)

plt.tight_layout()
plt.savefig('cell_proportion_boxplot.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Print results
print("\n=== Statistical Results (adjusted for sex and BMI, BH corrected) ===")
for ct in cell_types:
    p_adj = padj_dict.get(ct, np.nan)
    if not np.isnan(p_adj):
        print(f"{ct}: p_adj = {p_adj:.4f} {sig_label(p_adj)}")



