# R Six-Set Venn Diagram

This workflow adapts:

`E:/DWCT/OLD/Others/PD_Project_2026/Venn/colon_lfc/python_runs/lfc_files_PD/venn_intersections_6sets/run_venn_6sets_R.R`

It uses the mapped MyGene files from:

`June_17/venn_new/gene_id_mapping/mapped_resLFC`

## Run

From the repository root:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' June_17\venn_new\venn_diagram_R\code\run_venn_6sets_R.R
```

## Outputs

Saved in `plots_and_tables/`:

- `venn_6sets_gene_id_R.png`
- `venn_6sets_gene_id_R.pdf`
- `all_possible_exact_intersections_gene_lists_R.csv`
- `exact_intersection_counts_R.csv`
- `genes_shared_by_all_6_R.csv`

SVG export is skipped unless the R package `svglite` is installed.
