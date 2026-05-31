import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

# ============================================
# User settings - MODIFY THESE
# ============================================

CELL_TYPE = "cMono"
TARGET_GENES = ['HLA-DQA1', 'HLA-DQB1', 'C3AR1', 'EGR1', 'CX3CR1', 'RGS1', 'NR4A2', 'JUN']
OUTPUT_DIR = Path("./violin_plots/")
GROUP_COL = "exercise_group"
CONTROL_GROUP = "SED"
TREATMENT_GROUP = "EX"

# Colors
COLOR_CONTROL = '#7EBFC9'
COLOR_TREATMENT = '#BCDF7A'

# ============================================
# DO NOT MODIFY BELOW
# ============================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    'font.family': 'Arial',
    'axes.linewidth': 1.2,
})
sc.settings.set_figure_params(dpi=300, frameon=True, figsize=(2, 1.5))
sns.set_style("whitegrid")

# Extract cell type
adata_subset = adata[adata.obs['cell_type_l3'] == CELL_TYPE, :].copy()
print(f"Cells in {CELL_TYPE}: {adata_subset.n_obs}")

# Prepare data
plot_data = pd.DataFrame({
    'group': adata_subset.obs[GROUP_COL]
})

for gene in TARGET_GENES:
    if gene in adata_subset.var_names:
        plot_data[gene] = adata_subset[:, gene].X.toarray().flatten()
    else:
        print(f"Warning: {gene} not found")
        plot_data[gene] = np.nan

# Create PDF
output_file = OUTPUT_DIR / f'{CELL_TYPE}_violin_plots.pdf'
print(f"Saving to: {output_file}")

with PdfPages(output_file) as pdf:
    for gene in TARGET_GENES:
        if gene not in plot_data.columns or plot_data[gene].isna().all():
            continue
        
        fig, ax = plt.subplots(figsize=(4, 4.5))
        
        # Violin plot
        sns.violinplot(
            data=plot_data,
            x='group',
            y=gene,
            order=[CONTROL_GROUP, TREATMENT_GROUP],
            palette={CONTROL_GROUP: COLOR_CONTROL, TREATMENT_GROUP: COLOR_TREATMENT},
            cut=0,
            scale="width",
            inner=None,
            linewidth=0.7,
            alpha=1,
            ax=ax
        )
        
        # Boxplot overlay
        sns.boxplot(
            data=plot_data,
            x='group',
            y=gene,
            order=[CONTROL_GROUP, TREATMENT_GROUP],
            width=0.08,
            boxprops={'facecolor':'white', 'edgecolor':'black', 'linewidth':0.6, 'alpha':1},
            medianprops={'color':'black', 'linewidth':0.6},
            whiskerprops={'color':'black', 'linewidth':0.6},
            capprops={'color':'black', 'linewidth':0},
            showfliers=False,
            ax=ax
        )
        
        ax.set_title(f"{CELL_TYPE}: {gene}", pad=15, fontsize=12,
                    bbox=dict(facecolor='#E0E0E0', edgecolor='none', pad=3, boxstyle='square,pad=0.5'))
        ax.set_xlabel('')
        ax.set_ylabel('Expression Level', labelpad=10)
        ax.grid(axis='y', visible=False)
        
        for spine in ax.spines.values():
            spine.set_color('black')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Statistical test and log2FC
        ctrl_vals = plot_data.loc[plot_data['group'] == CONTROL_GROUP, gene].dropna()
        trt_vals = plot_data.loc[plot_data['group'] == TREATMENT_GROUP, gene].dropna()
        
        if len(ctrl_vals) > 0 and len(trt_vals) > 0:
            _, p_val = mannwhitneyu(ctrl_vals, trt_vals)
            
            # Calculate log2FC
            mean_ctrl = np.mean(ctrl_vals)
            mean_trt = np.mean(trt_vals)
            log2fc = np.log2(mean_trt / mean_ctrl) if mean_ctrl > 0 else np.nan
            
            y_max = plot_data[gene].max()
            y_min = plot_data[gene].min()
            y_range = y_max - y_min
            
            line_height = y_max + y_range * 0.05
            text_height = line_height + y_range * 0.02
            label_height = line_height - y_range * 0.03
            
            if p_val < 0.001:
                sig_symbol = '***'
            elif p_val < 0.01:
                sig_symbol = '**'
            elif p_val < 0.05:
                sig_symbol = '*'
            else:
                sig_symbol = 'ns'
            
            # Draw significance line
            ax.plot([0, 1], [line_height, line_height], color='black', lw=1)
            # Add significance symbol above the line
            ax.text(0.5, text_height, sig_symbol, ha='center', va='bottom', fontsize=12)
            # Add log2FC below the line (between line and violin top)
            if not np.isnan(log2fc):
                ax.text(0.5, label_height, f'log2FC = {log2fc:.3f}', 
                       ha='center', va='top', fontsize=8, color='black')
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()
        print(f"Plotted: {gene}")

print(f"\nDone! Output: {output_file}")
