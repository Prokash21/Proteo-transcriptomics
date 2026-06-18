# Six-Dataset Volcano Plots

This workflow creates volcano plots from the six mapped `resLFC` files in:

`June_17/venn_new/gene_id_mapping/mapped_resLFC`

## Thresholds

- `abs(log2FoldChange) >= 1`
- `padj < 0.05`

## Run

From the repository root:

```powershell
.\.venv\Scripts\python.exe June_17\volcano\code\make_volcano_plots.py
```

## Outputs

- Individual plots: `plots/volcano_*.png`, `.pdf`, `.svg`
- Collage: `plots/volcano_collage_6sets.png`, `.pdf`, `.svg`
- Plot data: `tables/volcano_data_*.csv`
- Summary counts: `tables/volcano_summary_counts.csv`
