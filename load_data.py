"""
Step 1 — Load Data
==================
Download GSE150404 series matrix from NCBI GEO,
auto-detect Cancer/Normal groups from sample titles,
and save as pickle + metadata JSON.

Output:
  output_GSE150404/df_exp.pkl
  output_GSE150404/metadata.json
"""

import pandas as pd
import numpy as np
import urllib.request
import json
import os
import pickle
import gzip

OUTPUT_DIR = "output_GSE150404"
os.makedirs(OUTPUT_DIR, exist_ok=True)

GEO_ID   = "GSE150404"
URL      = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/{GEO_ID}/matrix/{GEO_ID}_series_matrix.txt.gz"
FILENAME = f"{GEO_ID}_series_matrix.txt.gz"

# ── 1. Download ──────────────────────────────────
if not os.path.exists(FILENAME):
    print(f"Downloading {GEO_ID}...")
    urllib.request.urlretrieve(URL, FILENAME)
    print("  Download complete.")
else:
    print(f"  File already exists: {FILENAME}")

# ── 2. Auto-detect Cancer/Normal from sample titles ──
print("Detecting groups from sample titles...")
titles = []
with gzip.open(FILENAME, 'rt') as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = [t.strip().strip('"') for t in line.strip().split("\t")[1:]]
            break

groups = []
for t in titles:
    t_lower = t.lower()
    if any(k in t_lower for k in ['cancer', 'tumor', 'tumour', 'carcinoma', 'case', 'hcc', 'met']):
        groups.append('Cancer')
    elif any(k in t_lower for k in ['normal', 'healthy', 'control', 'adjacent', 'non-tumor', 'nontumor']):
        groups.append('Normal')
    else:
        groups.append('Unknown')

# ── 3. Read expression matrix ────────────────────
print("Reading data matrix...")
df_exp = pd.read_csv(FILENAME, sep="\t", comment="!", index_col=0)
df_exp = df_exp.apply(pd.to_numeric, errors='coerce')
df_exp = df_exp.dropna(how='all')

samples = df_exp.columns.tolist()

# Print detection results
print("\nGroup detection results:")
for s, g in zip(samples, groups):
    print(f"  {g:10} | {s}")

# Warn if all samples are undetected
if all(g == 'Unknown' for g in groups):
    print("\n  [!] No keywords matched in sample titles.")
    print("  [!] Check sample names above and adjust keywords in the code.")

cancer_cols = [s for s, g in zip(samples, groups) if g == 'Cancer']
normal_cols  = [s for s, g in zip(samples, groups) if g == 'Normal']

print(f"\n  Probes  : {df_exp.shape[0]:,}")
print(f"  Samples : {len(samples)} (Cancer={len(cancer_cols)}, Normal={len(normal_cols)})")

# ── 4. Save output ───────────────────────────────
pkl_path = os.path.join(OUTPUT_DIR, "df_exp.pkl")
with open(pkl_path, 'wb') as f:
    pickle.dump(df_exp, f)

meta = {
    "geo_id"     : GEO_ID,
    "n_probes"   : df_exp.shape[0],
    "n_samples"  : len(samples),
    "n_cancer"   : len(cancer_cols),
    "n_normal"   : len(normal_cols),
    "samples"    : samples,
    "groups"     : groups,
    "cancer_cols": cancer_cols,
    "normal_cols": normal_cols,
}
with open(os.path.join(OUTPUT_DIR, "metadata.json"), 'w') as f:
    json.dump(meta, f, indent=2)

print(f"  Saved: {pkl_path}")
