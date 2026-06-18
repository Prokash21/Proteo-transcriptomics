# Shared-Gene Expression Matrices

This workflow extracts expression for genes shared by all six datasets from:

`June_17/venn_new/venn_diagram_R/plots_and_tables/genes_shared_by_all_6_R.csv`

Rows are mapped to gene symbols using the six mapped `resLFC` files in:

`June_17/venn_new/gene_id_mapping/mapped_resLFC`

## Outputs

- `matrices/merged_raw_shared_4663_genes_with_batch_row.csv`
- `matrices/merged_normalized_shared_4663_genes_with_batch_row.csv`
- `matrices/merged_raw_shared_4663_genes_expression_only.csv`
- `matrices/merged_normalized_shared_4663_genes_expression_only.csv`
- `metadata/batch_metadata.csv`
- `metadata/sample_metadata_all.csv`
- `reports/shared_gene_extraction_report.csv`

The `with_batch_row` files use the CSV header as sample names, the second row as batch numbers, and subsequent rows as genes.

## Batch Numbers

- `1`: GSE180440
- `2`: GSE251845
- `3`: GSE50760
- `4`: GSE142279
- `5`: TCGA COAD
- `6`: TCGA READ

Raw and normalized matrices are saved separately because their expression scales differ.
