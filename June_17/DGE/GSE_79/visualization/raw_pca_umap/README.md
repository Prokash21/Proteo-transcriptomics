# Raw PCA / UMAP

This folder contains a standalone Python script that generates publication-ready PCA and UMAP plots directly from the raw RNA-seq count matrix for GSE142279.

The script does not perform between-sample normalization. It uses a log2(x+1) transform only to stabilize the variance before embedding.

## Input files

- `../raw_data/count_data_raw_GSE142279.csv`
- `../raw_data/meta_data_GSE142279.csv`

## Usage

```bash
python raw_pca_umap.py \
  --counts ../raw_data/count_data_raw_GSE142279.csv \
  --meta ../raw_data/meta_data_GSE142279.csv \
  --outdir ./plots_raw_counts
```

## Output

- `raw_counts_pca.png` and `raw_counts_pca.pdf`
- `raw_counts_umap.png` and `raw_counts_umap.pdf`
- `raw_counts_pca_umap_coordinates.csv`

## Notes

- The metadata loader accepts a sample column named `sample`, `sample_name`, `sample _name`, or similar variants.
- The script keeps the top variable genes by default before PCA/UMAP. Adjust `--top-genes` if you want to retain more or fewer features.