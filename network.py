"""
Step 7 — STRING Protein Interaction Network
=============================================
Queries the STRING API for a protein-protein interaction network
built from the top differentially expressed genes.

Input:
  output_GSE150404/stats.csv

Output:
  output_GSE150404/network.json  — Nodes + edges for Plotly network
"""

import pandas as pd
import numpy as np
import json
import os
import requests

OUTPUT_DIR = "output_GSE150404"

stats = pd.read_csv(os.path.join(OUTPUT_DIR, "stats.csv"))

# Top 50 DEGs with a gene symbol
sig = stats[(stats['DEG'] != 'NS') & (stats['gene_symbol'].notna())]
query_genes = sig.sort_values('p_adj')['gene_symbol'].unique().tolist()[:50]
print(f"Querying STRING: {len(query_genes)} genes")

result = {"status": "skipped", "nodes": [], "edges": []}

if len(query_genes) < 3:
    print("  Too few genes, network skipped.")
else:
    try:
        deg_lookup   = stats.dropna(subset=['gene_symbol']).set_index('gene_symbol')['DEG'].to_dict()
        logfc_lookup = stats.dropna(subset=['gene_symbol']).set_index('gene_symbol')['logFC'].to_dict()

        r = requests.get(
            "https://string-db.org/api/json/network",
            params={
                'identifiers'    : '%0d'.join(query_genes),
                'species'        : 9606,
                'caller_identity': 'GSE150404_pipeline',
                'required_score' : 700,
            },
            timeout=25
        )

        if r.status_code != 200:
            raise ValueError(f"STRING API returned {r.status_code}")

        edges_raw = r.json()
        if not edges_raw:
            print("  No edges found in STRING.")
            result["status"] = "no_edges"
        else:
            # Collect unique nodes
            nodes_set = set()
            edge_list = []
            for e in edges_raw:
                a = e['preferredName_A']
                b = e['preferredName_B']
                nodes_set.add(a)
                nodes_set.add(b)
                edge_list.append({
                    "source": a,
                    "target": b,
                    "score" : round(float(e['score']), 3),
                })

            # Circular layout
            nodes_list = list(nodes_set)
            n = len(nodes_list)
            angles = np.linspace(0, 2 * np.pi, n, endpoint=False)

            nodes = []
            for i, gene in enumerate(nodes_list):
                nodes.append({
                    "id"   : gene,
                    "x"    : float(np.cos(angles[i])),
                    "y"    : float(np.sin(angles[i])),
                    "deg"  : deg_lookup.get(gene, 'NS'),
                    "logfc": round(float(logfc_lookup.get(gene, 0)), 3),
                })

            result = {
                "status" : "ok",
                "n_nodes": n,
                "n_edges": len(edge_list),
                "nodes"  : nodes,
                "edges"  : edge_list,
            }
            print(f"  Nodes: {n}  Edges: {len(edge_list)}")

    except requests.exceptions.ConnectionError:
        print("  Connection to STRING failed (offline?).")
        result["status"] = "connection_error"
    except Exception as e:
        print(f"  Network error: {e}")
        result["status"] = f"error: {e}"

out_path = os.path.join(OUTPUT_DIR, "network.json")
with open(out_path, 'w') as f:
    json.dump(result, f, indent=2)
print(f"  Saved: {out_path}")
