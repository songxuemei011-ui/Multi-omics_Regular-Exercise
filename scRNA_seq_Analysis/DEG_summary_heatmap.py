#!/usr/bin/env python
# ============================================
# DEG Summary Heatmap
# Count up/down regulated genes per cell type
# ============================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ============================================
# USER SETTINGS - MODIFY THESE
# ============================================

# Path to DEG results directory
RESULTS_DIR = "./DEG_results/"

# Cell types to include (change to your list)
CELL_TYPES = [
    'CD8_CTL', 'CD8_Tcm', 'CD8_Tem', 'CD8_Tn'
]

# DEG thresholds
P_ADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 0.25

# Output
OUTPUT_FILE = "DEG_summary_heatmap.pdf"

# ============================================
# DO NOT MODIFY BELOW
# ============================================

# Collect results
results = {'Cell Type': [], 'Upregulated': [], 'Downregulated': [], 'Total DEGs': []}

for file in os.listdir(RESULTS_DIR):
    if file.endswith('_significant.csv') or file.endswith('_DEGs.csv'):
        # Extract cell type name
        celltype = file.split('_')[0]
        
        if celltype in CELL_TYPES:
            file_path = os.path.join(RESULTS_DIR, file)
            
            try:
                df = pd.read_csv(file_path)
                
                # Count up/down regulated
                up = df[(df['p_adj'] < P_ADJ_THRESHOLD) & (df['log2FC'] > LOG2FC_THRESHOLD)]
                down = df[(df['p_adj'] < P_ADJ_THRESHOLD) & (df['log2FC'] < -LOG2FC_THRESHOLD)]
                
                results['Cell Type'].append(celltype)
                results['Upregulated'].append(len(up))
                results['Downregulated'].append(len(down))
                results['Total DEGs'].append(len(up) + len(down))
                
            except Exception as e:
                print(f"Error processing {file}: {e}")

# Create DataFrame
df_results = pd.DataFrame(results)
df_results = df_results.sort_values('Total DEGs', ascending=False)

# Prepare heatmap data
heatmap_data = df_results.set_index('Cell Type')[['Upregulated', 'Downregulated']]

# Plot
plt.figure(figsize=(8, 4))
sns.heatmap(
    heatmap_data.T,
    annot=True,
    fmt='d',
    cmap='RdYlBu_r',
    linewidths=0.5,
    cbar_kws={'label': 'Number of DEGs'}
)

plt.xticks(rotation=45, ha='right')
plt.ylabel('Regulation')
plt.title('DEG Summary per Cell Type')

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight')
plt.show()

print(f"Saved: {OUTPUT_FILE}")
