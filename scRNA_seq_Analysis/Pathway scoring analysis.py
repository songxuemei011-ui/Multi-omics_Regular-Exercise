# ============================================
# AddModule Score Pathway Analysis for Myeloid Cells
# ============================================

import glob
import os
import scanpy as sc
import pandas as pd

# ============================================
# USER SETTINGS - MODIFY THESE
# ============================================

# Input data
# adata = sc.read_h5ad("data/processed/scRNA-seq.h5ad")

# Path to gene set files (TSV format)
GENE_SET_DIR = "data/gene_sets/Myeloid_Pathway/"

# Output directory
OUTPUT_DIR = "results/PathwayScore/Myeloid/"

# Cell types to include
MYELOID_SUBTYPES = ['cDC', 'cMono', 'intMono', 'ncMono', 'pDC']

# Metadata column names
SAMPLE_COL = "sample"
GROUP_COL = "exercise_group"
CELL_TYPE_L1_COL = "cell_type_l1"
CELL_TYPE_L3_COL = "cell_type_l3"

# ============================================
# DO NOT MODIFY BELOW
# ============================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load gene sets
gene_sets = {}
tsv_files = glob.glob(os.path.join(GENE_SET_DIR, '*.tsv'))
print(f"Found {len(tsv_files)} gene set files")

for file_path in tsv_files:
    with open(file_path, 'r') as file:
        gene_set_name = os.path.splitext(os.path.basename(file_path))[0]
        gene_list = []
        for line in file:
            line = line.strip()
            if line.startswith("GENE_SYMBOLS\t"):
                genes = line.split("\t")[1].strip()
                if genes:
                    gene_list = [gene.strip() for gene in genes.split(",") if gene.strip()]
                break
        if gene_list:
            gene_sets[gene_set_name] = gene_list
            print(f"Loaded: {gene_set_name} ({len(gene_list)} genes)")

# Subset myeloid cells
myeloid_subset = adata[adata.obs[CELL_TYPE_L3_COL].isin(MYELOID_SUBTYPES)].copy()

# Normalize if needed
if myeloid_subset.X.max() > 50:
    sc.pp.normalize_total(myeloid_subset, target_sum=1e4)
    sc.pp.log1p(myeloid_subset)

# Calculate pathway scores
valid_pathways = []
for gene_set_name, gene_list in gene_sets.items():
    valid_genes = [g for g in gene_list if g in myeloid_subset.var_names]
    if len(valid_genes) < 3:
        continue
    sc.tl.score_genes(myeloid_subset, gene_list=valid_genes, 
                      score_name=gene_set_name, use_raw=False)
    valid_pathways.append(gene_set_name)

# Save results
result_df = pd.DataFrame({
    SAMPLE_COL: myeloid_subset.obs[SAMPLE_COL],
    GROUP_COL: myeloid_subset.obs[GROUP_COL],
    CELL_TYPE_L1_COL: myeloid_subset.obs[CELL_TYPE_L1_COL],
    CELL_TYPE_L3_COL: myeloid_subset.obs[CELL_TYPE_L3_COL]
})

for pathway in valid_pathways:
    result_df[pathway] = myeloid_subset.obs[pathway]

output_file = os.path.join(OUTPUT_DIR, 'Myeloid_pathway_scores.csv')
result_df.to_csv(output_file, index=False)
print(f"Results saved to: {output_file}")
# ============================================
# Pathway Score Boxplot Visualization
# BH-corrected Mann-Whitney U test between groups
# ============================================

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from matplotlib.backends.backend_pdf import PdfPages
import os

# ============================================
# USER SETTINGS - MODIFY THESE
# ============================================

# Input and output paths
INPUT_FILE = "results/PathwayScore/Myeloid/Myeloid_pathway_scores.csv"
OUTPUT_DIR = "results/PathwayScore/Myeloid/plots/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Pathways to plot (select from your CSV columns)
SELECTED_PATHWAYS = [
    'GOBP_POSITIVE_REGULATION_OF_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY.v2025.1.Hs',
    'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION.v2025.1.Hs',
    'GOBP_INNATE_IMMUNE_RESPONSE.v2025.1.Hs',
    'GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX.v2025.1.Hs',
    'GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_I_PROTEIN_COMPLEX.v2025.1.Hs',
    'HALLMARK_HYPOXIA.v2025.1.Hs'
]

# Cell type order (left to right on x-axis)
CELL_TYPE_ORDER = ['cMono', 'intMono', 'ncMono', 'cDC', 'pDC']

# Group names
CONTROL_GROUP = "Inactive"
TREATMENT_GROUP = "Regularly Active"
DISPLAY_CONTROL = "SED"
DISPLAY_TREATMENT = "EX"

# Colors
COLOR_CONTROL = '#7EBFC9'   # Light blue for SED/Inactive
COLOR_TREATMENT = '#BCDF7A' # Light green for EX/Active

# Plot settings
FIGURE_SIZE = (8, 4)
DPI = 300
FONT_FAMILY = 'Arial'

# Statistical test
ALPHA = 0.05  # Significance level

# ============================================
# DO NOT MODIFY BELOW
# ============================================

# Load data
print("=" * 60)
print("Pathway Score Visualization")
print("=" * 60)

score_df = pd.read_csv(INPUT_FILE)
print(f"Data loaded: {score_df.shape[0]} cells")

# Preprocess data
score_df['cell_type_l3'] = pd.Categorical(
    score_df['cell_type_l3'], 
    categories=CELL_TYPE_ORDER, 
    ordered=True
)
score_df = score_df.sort_values('cell_type_l3')

# Set plot style
sns.set_style("whitegrid")
plt.rcParams.update({
    'font.family': FONT_FAMILY,
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.linewidth': 1.0,
    'grid.color': '0.9'
})

# Color palette
palette = {
    CONTROL_GROUP: COLOR_CONTROL, 
    TREATMENT_GROUP: COLOR_TREATMENT
}

# Function: Perform all Mann-Whitney U tests
def perform_all_tests(data, pathways):
    results = []
    for pathway in pathways:
        if data[pathway].isnull().all():
            continue
            
        for cell_type in CELL_TYPE_ORDER:
            group_data = data[data['cell_type_l3'] == cell_type]
            if len(group_data) < 3:
                continue
                
            control = group_data[group_data['exercise_group'] == CONTROL_GROUP][pathway]
            treatment = group_data[group_data['exercise_group'] == TREATMENT_GROUP][pathway]
            
            try:
                _, p_val = mannwhitneyu(control, treatment, alternative='two-sided')
                y_max = group_data[pathway].max()
                y_min = group_data[pathway].min()
                y_range = y_max - y_min
                line_height = y_max + 0.001 * y_range
                text_height = line_height + 0.0001 * y_range
                
                results.append({
                    'pathway': pathway,
                    'cell_type': cell_type,
                    'p_val': p_val,
                    'line_height': line_height,
                    'text_height': text_height
                })
            except:
                pass
    return results

# Function: BH correction for multiple testing
def adjust_pvalues(results):
    p_values = [r['p_val'] for r in results]
    _, corrected_p, _, _ = multipletests(p_values, method='fdr_bh')
    
    for i, result in enumerate(results):
        result['corrected_p'] = corrected_p[i]
    return results

# Function: Plot boxplot with significance
def plot_celltype_pathway_boxplot(data, pathway, ax, corrected_results):
    # Draw boxplot
    sns.boxplot(
        data=data,
        x='cell_type_l3',
        y=pathway,
        hue='exercise_group',
        hue_order=[CONTROL_GROUP, TREATMENT_GROUP],
        palette=palette,
        width=0.6,
        gap=0.2,
        linewidth=1.0,
        showfliers=False,
        ax=ax
    )
    
    # Style boxes
    for i, artist in enumerate(ax.artists):
        artist.set_edgecolor('black')
        artist.set_facecolor(palette[CONTROL_GROUP if i % 2 == 0 else TREATMENT_GROUP])
        artist.set_alpha(0.8)
        artist.set_linewidth(0.8)
        ax.lines[i*6+4].set_color('black')
        ax.lines[i*6+4].set_linewidth(1.5)
    
    # Add significance annotations
    for i, cell_type in enumerate(CELL_TYPE_ORDER):
        result = next((r for r in corrected_results 
                      if r['pathway'] == pathway and r['cell_type'] == cell_type), None)
        
        if result and result['corrected_p'] < ALPHA:
            if result['corrected_p'] < 0.001:
                sig_symbol = '***'
            elif result['corrected_p'] < 0.01:
                sig_symbol = '**'
            else:
                sig_symbol = '*'
            
            # Draw significance line and text
            ax.plot([i-0.2, i+0.2], [result['line_height']]*2, color='black', lw=1.2)
            ax.text(i, result['text_height'], sig_symbol, 
                   ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    # Labels and title
    # Shorten pathway name for display
    short_name = pathway.replace('.v2025.1.Hs', '').replace('GOBP_', '').replace('HALLMARK_', '')
    ax.set_title(short_name, pad=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Pathway Score', labelpad=8)
    
    # Legend with SED/EX labels
    handles = ax.get_legend_handles_labels()[0]
    ax.legend(
        handles=handles, 
        labels=[DISPLAY_CONTROL, DISPLAY_TREATMENT],
        title='Group', 
        bbox_to_anchor=(1.02, 1), 
        loc='upper left'
    )
    
    ax.grid(True, axis='y', alpha=0.3)
    ax.grid(False, axis='x')
    plt.setp(ax.get_xticklabels(), rotation=0, ha='center')

# Main execution
print("\n" + "-" * 40)
print("Step 1: Filtering pathways")
available_pathways = [p for p in SELECTED_PATHWAYS if p in score_df.columns]
filtered_pathways = [p for p in available_pathways if not score_df[p].isnull().all()]
print(f"Pathways to plot: {len(filtered_pathways)}")

print("\nStep 2: Running statistical tests")
all_results = perform_all_tests(score_df, filtered_pathways)
print(f"Tests performed: {len(all_results)}")

print("\nStep 3: BH correction")
corrected_results = adjust_pvalues(all_results)
significant = sum(1 for r in corrected_results if r['corrected_p'] < ALPHA)
print(f"Significant after correction: {significant}")

print("\nStep 4: Generating plots")
output_pdf = os.path.join(OUTPUT_DIR, 'Myeloid_pathway_boxplots_BH_corrected.pdf')

with PdfPages(output_pdf) as pdf:
    for pathway in filtered_pathways:
        fig, ax = plt.subplots(figsize=FIGURE_SIZE)
        plot_celltype_pathway_boxplot(score_df, pathway, ax, corrected_results)
        plt.tight_layout()
        pdf.savefig(fig, dpi=DPI, bbox_inches='tight')
        plt.close()
        print(f"  Plotted: {pathway}")

print("\n" + "=" * 60)
print(f"SUCCESS! Plot saved to: {output_pdf}")
print("=" * 60)




