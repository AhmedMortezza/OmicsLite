"""
Step 3 — PCA & UMAP
====================
Dimensionality reduction to visualize Cancer vs Normal separation.
Produces PCA coordinates (always) and UMAP coordinates (if available).

Input:
  output_GSE150404/df_exp.pkl
  output_GSE150404/metadata.json

Output:
  output_GSE150404/pca_umap.json  — Coordinates + variance explained
"""

import pandas as pd
import numpy as np
import json
import os
import pickle
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

OUTPUT_DIR = "output_GSE150404"

with open(os.path.join(OUTPUT_DIR, "df_exp.pkl"), 'rb') as f:
    df_exp = pickle.load(f)
with open(os.path.join(OUTPUT_DIR, "metadata.json")) as f:
    meta = json.load(f)

samples = meta['samples']
groups  = meta['groups']

# ── Preprocessing ────────────────────────────────
print("Preprocessing for dimensionality reduction...")
col_means = df_exp.mean(axis=1)
X = df_exp.T.copy()
X = X.fillna(col_means)

scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X)
print(f"  Matrix shape: {X_scaled.shape}")

# ── PCA ──────────────────────────────────────────
print("Computing PCA...")
pca = PCA(n_components=min(10, X_scaled.shape[0]-1))
pca_coords = pca.fit_transform(X_scaled)
var_exp = (pca.explained_variance_ratio_ * 100).tolist()
print(f"  PC1: {var_exp[0]:.1f}%  PC2: {var_exp[1]:.1f}%  PC3: {var_exp[2]:.1f}%")

result = {
    "pca": {
        "coords" : pca_coords[:, :3].tolist(),
        "var_exp": var_exp[:3],
        "samples": samples,
        "groups" : groups,
    }
}

# ── UMAP ─────────────────────────────────────────
try:
    from umap import UMAP
    print("Computing UMAP...")
    reducer = UMAP(n_components=2, random_state=42,
                   n_neighbors=min(15, len(samples) - 1))
    umap_coords = reducer.fit_transform(X_scaled)
    result["umap"] = {
        "coords" : umap_coords.tolist(),
        "samples": samples,
        "groups" : groups,
    }
    print("  UMAP complete.")
except ImportError:
    print("  UMAP not available (pip install umap-learn), skipping.")
    result["umap"] = None

out_path = os.path.join(OUTPUT_DIR, "pca_umap.json")
with open(out_path, 'w') as f:
    json.dump(result, f)
print(f"  Saved: {out_path}")
