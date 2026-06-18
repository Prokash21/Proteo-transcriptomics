# resLFC Gene ID Mapping

This folder maps the six `resLFC` files from `June_17/venn_new/data`.

## Code

- `code/map_reslfc_gene_ids.py`

The script detects whether the first column contains Entrez or Ensembl IDs, strips Ensembl version suffixes for querying, maps IDs with the MyGene API, and writes mapped CSV files.

## Run

From the repository root:

```powershell
.\.venv\Scripts\python.exe June_17\venn_new\gene_id_mapping\code\map_reslfc_gene_ids.py
```

## Outputs

- `mapped_resLFC/`: mapped CSV files
- `reports/id_mapping_summary.csv`: one-row summary per input file
- `reports/*_id_mapping_report.json`: detailed mapping report per input file

## Detected ID Types

| File | Detected ID type |
| --- | --- |
| `GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal.csv` | Entrez |
| `GSE_45_OR_02_resLFC_condition_tumor_vs_normal.csv` | Ensembl |
| `GSE_60_resLFC_raw_wolTreatment_tumor_vs_normal.csv` | Entrez |
| `GSE_79_02_resLFC_Treatment_tumor_vs_normal.csv` | Entrez |
| `tcga_coad_resLFC_Treatment_Tumor_vs_Normal.csv` | Ensembl |
| `tcga_read_resLFC_SampleType_Tumor_vs_Normal.csv` | Ensembl |
