import pandas as pd
import numpy as np
import json
import os
import pickle
import mygene
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

OUTPUT_DIR = "output_GSE150404"

# Load data
with open(os.path.join(OUTPUT_DIR, "df_exp.pkl"), 'rb') as f:
    df_exp = pickle.load(f)
with open(os.path.join(OUTPUT_DIR, "metadata.json")) as f:
    meta = json.load(f)

cancer_cols = meta['cancer_cols']
normal_cols  = meta['normal_cols']

# ── Gene Symbol Mapping ──────────────────────────
print("Mapping Probe IDs to Gene Symbols (MyGene API)...")
mg = mygene.MyGeneInfo()
mapping = {}
try:
    res = mg.querymany(df_exp.index.tolist(), scopes='reporter',
                       fields='symbol', species='human', as_dataframe=True)
    mapping = res['symbol'].dropna().to_dict()
    print(f"  Mapped: {len(mapping):,} probes")
except Exception as e:
    print(f"  Mapping error: {e}")

with open(os.path.join(OUTPUT_DIR, "mapping.json"), 'w') as f:
    json.dump(mapping, f)

# ── Statistics ───────────────────────────────────
print("Calculating logFC + t-test...")
mean_cancer = df_exp[cancer_cols].mean(axis=1)
mean_normal  = df_exp[normal_cols].mean(axis=1)
logFC = mean_cancer - mean_normal

t_stat, p_vals = ttest_ind(
    df_exp[cancer_cols].values,
    df_exp[normal_cols].values,
    axis=1
)
p_vals = np.where(np.isnan(p_vals), 1.0, p_vals)

# ── FDR Correction ───────────────────────────────
print("Applying FDR correction (Benjamini-Hochberg)...")
_, p_adj, _, _ = multipletests(p_vals, method='fdr_bh')

# ── Build Stats DataFrame ────────────────────────
stats = pd.DataFrame({
    'probe_id'        : df_exp.index,
    'gene_symbol'     : df_exp.index.map(mapping),
    'logFC'           : logFC.values,
    'p_value'         : p_vals,
    'p_adj'           : p_adj,
    'neg_log10_pv'    : -np.log10(np.clip(p_vals, 1e-300, 1)),
    'neg_log10_padj'  : -np.log10(np.clip(p_adj,  1e-300, 1)),
    'mean_expression' : df_exp.mean(axis=1).values,
    'mean_cancer'     : mean_cancer.values,
    'mean_normal'     : mean_normal.values,
})

def label_deg(row):
    if row['p_adj'] < 0.05 and row['logFC'] >  1.5: return 'Up'
    if row['p_adj'] < 0.05 and row['logFC'] < -1.5: return 'Down'
    return 'NS'

stats['DEG'] = stats.apply(label_deg, axis=1)

n_up   = (stats['DEG'] == 'Up').sum()
n_down = (stats['DEG'] == 'Down').sum()
print(f"  Up-regulated  : {n_up:,}")
print(f"  Down-regulated: {n_down:,}")

csv_path = os.path.join(OUTPUT_DIR, "stats.csv")
stats.to_csv(csv_path, index=False)
print(f"  Saved: {csv_path}")
