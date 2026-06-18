# GSE_79 Visualization

This folder contains a Python script to produce sample-level visualizations (PCA, UMAP, hierarchical clustering dendrogram, and boxplots) before and after a simple normalization.

Files:
- `gse79_plot_samples.py`: main script
- `requirements.txt`: Python packages required

Usage:

1. Install dependencies (preferably in a virtual environment):

```bash
pip install -r requirements.txt
```

2. Run the script pointing to your counts matrix and metadata CSV. Example:

```bash
python gse79_plot_samples.py --counts path/to/GSE79_counts.csv --meta path/to/GSE79_meta.csv --outdir ./plots_GSE79
```

Notes:
- The script expects a counts CSV with genes as rows and samples as columns. If your file has samples as rows, the script will attempt to transpose.
- Metadata should have sample IDs as the index or include a `sample` column; if a `group` column exists it will be used for coloring PCA/UMAP plots.
- "Before" plots use raw counts (log2 for visualization). "After" plots use CPM normalization by default.
