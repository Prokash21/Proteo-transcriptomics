# Six-Set Venn Diagram

This folder creates plots and intersection tables from the mapped `resLFC` files in:

`June_17/venn_new/gene_id_mapping/mapped_resLFC`

Dataset labels are shown with accession-style names:

- `GSE180440`
- `GSE251845`
- `GSE50760`
- `GSE142279`
- `TCGA COAD`
- `TCGA READ`

## Run

From the repository root:

```powershell
.\.venv\Scripts\python.exe June_17\venn_new\venn_diagram\code\plot_venn_6sets.py
```

## Outputs

Saved in `plots_and_tables/`:

- `venn_6sets_gene_id.png`
- `venn_6sets_gene_id.pdf`
- `venn_6sets_gene_id.svg`
- `gene_set_sizes.png`
- `exact_intersection_counts.csv`
- `combination_overlap_counts.csv`
- `gene_membership_by_dataset.csv`
- `genes_shared_by_all_6.csv`

The Venn circles are schematic because six-set Venn areas cannot be made accurately proportional in a simple circle layout. Exact overlap values are saved in the CSV tables.
