import pandas as pd
import numpy as np
import json
import os

OUTPUT_DIR = "output_GSE150404"

stats = pd.read_csv(os.path.join(OUTPUT_DIR, "stats.csv"))

# ── Prepare gene list ────────────────────────────
sig = stats[(stats['DEG'] != 'NS') & (stats['gene_symbol'].notna())]
sig_genes = sig.sort_values('p_adj')['gene_symbol'].unique().tolist()[:300]
print(f"Genes for enrichment: {len(sig_genes)}")

result = {"status": "skipped", "pathways": [], "table_rows": []}

if len(sig_genes) < 5:
    print("  Too few genes, enrichment skipped.")
else:
    try:
        import gseapy as gp
        print("Running Enrichr (KEGG + GO)...")

        enr = gp.enrichr(
            gene_list=sig_genes,
            gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023'],
            organism='human',
            outdir=None,
        )

        df_res = pd.DataFrame(enr.results)

        if df_res.empty:
            print("  No enrichment results returned.")
            result["status"] = "empty"
        else:
            # Filter FDR < 0.05, keep top 20
            df_sig = df_res[df_res['Adjusted P-value'] < 0.05]
            if df_sig.empty:
                print("  No pathways with FDR < 0.05, showing top 20 trends.")
                df_sig = df_res.head(20)
            else:
                df_sig = df_sig.head(20)

            print(f"  Pathways: {len(df_sig)}")

            # Data for bar chart
            pathways = []
            for _, r in df_sig.sort_values('Combined Score').iterrows():
                pathways.append({
                    "term"          : r['Term'],
                    "gene_set"      : r['Gene_set'],
                    "combined_score": round(float(r['Combined Score']), 2),
                    "p_adj"         : float(r['Adjusted P-value']),
                    "overlap"       : r['Overlap'],
                    "genes"         : r['Genes'],
                })

            # Data for table
            table_cols = ['Gene_set', 'Term', 'Overlap', 'Adjusted P-value', 'Genes']
            table_rows = df_sig[table_cols].fillna('—').to_dict(orient='records')

            result = {
                "status"    : "ok",
                "n_genes"   : len(sig_genes),
                "pathways"  : pathways,
                "table_rows": table_rows,
            }

    except ImportError:
        print("  gseapy not available (pip install gseapy)")
        result["status"] = "no_gseapy"
    except Exception as e:
        print(f"  Enrichment error: {e}")
        result["status"] = f"error: {e}"

out_path = os.path.join(OUTPUT_DIR, "enrichment.json")
with open(out_path, 'w') as f:
    json.dump(result, f, indent=2)
print(f"  Saved: {out_path}")
