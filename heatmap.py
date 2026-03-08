"""
Step 4 — Heatmap
================
Heatmap of the top 50 DEGs with Z-score normalization.
Produces a JSON file for rendering as a Plotly heatmap in the report.

Input:
  output_GSE150404/df_exp.pkl
  output_GSE150404/stats.csv
  output_GSE150404/metadata.json

Output:
  output_GSE150404/heatmap.json
"""

import pandas as pd
import numpy as np
import json
import os
import pickle

OUTPUT_DIR = "output_GSE150404"

with open(os.path.join(OUTPUT_DIR, "df_exp.pkl"), 'rb') as f:
    df_exp = pickle.load(f)

stats = pd.read_csv(os.path.join(OUTPUT_DIR, "stats.csv"))

with open(os.path.join(OUTPUT_DIR, "metadata.json")) as f:
    meta = json.load(f)

cancer_cols = meta['cancer_cols']
normal_cols  = meta['normal_cols']

# ── Select Top 50 DEGs ───────────────────────────
print("Preparing top 50 DEGs for heatmap...")
top50 = stats[stats['DEG'] != 'NS'].sort_values('p_adj').head(50)

if top50.empty:
    print("  No significant DEGs found, using top 50 by p-value...")
    top50 = stats.sort_values('p_value').head(50)

# Subset expression data
df_hm = df_exp.loc[top50['probe_id']].copy()

# Use gene symbol as label if available
labels = top50.apply(
    lambda r: r['gene_symbol'] if pd.notna(r.get('gene_symbol')) else r['probe_id'],
    axis=1
).tolist()
df_hm.index = labels

# Column order: Cancer first, Normal second
col_order = cancer_cols + normal_cols
df_hm = df_hm[col_order]

# Z-score per gene
row_mean = df_hm.mean(axis=1)
row_std  = df_hm.std(axis=1).replace(0, 1)
df_hm_z  = df_hm.subtract(row_mean, axis=0).divide(row_std, axis=0)
df_hm_z  = df_hm_z.clip(-3, 3)

# ── Save JSON ────────────────────────────────────
result = {
    "z_matrix"  : df_hm_z.values.tolist(),
    "genes"     : df_hm_z.index.tolist(),
    "samples"   : col_order,
    "n_cancer"  : len(cancer_cols),
    "n_normal"  : len(normal_cols),
    "deg_status": top50['DEG'].tolist(),
}

out_path = os.path.join(OUTPUT_DIR, "heatmap.json")
with open(out_path, 'w') as f:
    json.dump(result, f)

print(f"  Genes : {len(labels)}")
print(f"  Saved : {out_path}")
