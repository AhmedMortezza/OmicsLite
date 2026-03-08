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

samples = meta['samples']
groups  = meta['groups']

# ── Top 12 DEGs ──────────────────────────────────
print("Preparing box plot data for top 12 DEGs...")
top12 = stats[stats['DEG'] != 'NS'].sort_values('p_adj').head(12)

if top12.empty:
    top12 = stats.sort_values('p_value').head(12)

genes_data = []
for _, row in top12.iterrows():
    label = row['gene_symbol'] if pd.notna(row.get('gene_symbol')) else row['probe_id']
    probe = row['probe_id']

    cancer_vals = df_exp.loc[probe, meta['cancer_cols']].tolist()
    normal_vals  = df_exp.loc[probe, meta['normal_cols']].tolist()

    genes_data.append({
        "gene"      : label,
        "probe_id"  : probe,
        "deg_status": row['DEG'],
        "logFC"     : round(float(row['logFC']), 4),
        "p_adj"     : float(row['p_adj']),
        "cancer"    : [round(v, 4) for v in cancer_vals],
        "normal"    : [round(v, 4) for v in normal_vals],
    })

result = {"genes": genes_data}
out_path = os.path.join(OUTPUT_DIR, "boxplot.json")
with open(out_path, 'w') as f:
    json.dump(result, f)

print(f"  Genes : {len(genes_data)}")
print(f"  Saved : {out_path}")
