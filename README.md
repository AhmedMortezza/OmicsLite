install all library
pip install pandas numpy scipy statsmodels scikit-learn umap-learn mygene gseapy plotly kaleido requests

GSE150404_project/
├── run_pipeline.py          # Main runner — executes all steps in order
├── requirements.txt         # Python dependencies
├── analysis/
│   ├── load_data.py         # Step 1: Download data, detect groups
│   ├── deg_analysis.py      # Step 2: logFC, t-test, FDR correction
│   ├── pca_umap.py          # Step 3: PCA and UMAP dimensionality reduction
│   ├── heatmap.py           # Step 4: Top 50 DEGs heatmap
│   ├── boxplot.py           # Step 5: Top 12 DEGs box plot
│   ├── enrichment.py        # Step 6: KEGG and GO enrichment via Enrichr
│   └── network.py           # Step 7: STRING protein interaction network
└── report/
    ├── build_report.py      # Step 8: Assemble HTML report + export JPGs
    └── style.css            # Report stylesheet

    
