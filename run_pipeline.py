import subprocess
import sys
import os
import time

STEPS = [
    ("1. Load Data",         "analysis/load_data.py"),
    ("2. DEG Analysis",      "analysis/deg_analysis.py"),
    ("3. PCA & UMAP",        "analysis/pca_umap.py"),
    ("4. Heatmap",           "analysis/heatmap.py"),
    ("5. Box Plot",          "analysis/boxplot.py"),
    ("6. Enrichment",        "analysis/enrichment.py"),
    ("7. Network (STRING)",  "analysis/network.py"),
    ("8. Build Report",      "report/build_report.py"),
]

print("=" * 60)
print("  GSE150404 PIPELINE — START")
print("=" * 60)

total_start = time.time()
for label, script in STEPS:
    print(f"\n▶  {label}")
    t = time.time()
    result = subprocess.run([sys.executable, script], capture_output=False)
    elapsed = time.time() - t
    if result.returncode != 0:
        print(f"   ❌ FAILED ({elapsed:.1f}s) — pipeline stopped.")
        sys.exit(1)
    print(f"   ✅ Done ({elapsed:.1f}s)")

print(f"\n{'='*60}")
print(f"  🎉 PIPELINE COMPLETE in {time.time()-total_start:.1f}s")
print(f"  📂 Output: output_GSE150404/")
print(f"  📄 Report: output_GSE150404/GSE150404_Report.html")
print(f"{'='*60}")
