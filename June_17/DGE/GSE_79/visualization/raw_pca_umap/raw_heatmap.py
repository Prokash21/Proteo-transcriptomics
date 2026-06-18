"""Create publication-ready hierarchical clustering heatmap from raw RNA-seq counts.

The script generates a heatmap with hierarchical clustering of both genes and samples.
A log2(x+1) transform is applied to the raw counts for visualization.

Example:
python raw_heatmap.py \
  --counts ../raw_data/count_data_raw_GSE142279.csv \
  --meta ../raw_data/meta_data_GSE142279.csv \
  --outdir ./plots_heatmap
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform


SAMPLE_COLUMN_CANDIDATES = [
    "sample",
    "sample_name",
    "sample name",
    "sample_id",
    "sample id",
    "sample_name ",
    "sample _name",
]

TISSUE_COLUMN_CANDIDATES = ["tissue", "group", "condition", "phenotype"]


def load_counts(path: Path) -> pd.DataFrame:
    counts = pd.read_csv(path, index_col=0)
    counts.index = counts.index.astype(str).str.strip()
    counts.columns = counts.columns.astype(str).str.strip()
    numeric = counts.apply(pd.to_numeric, errors="coerce")

    if numeric.shape[0] < numeric.shape[1]:
        numeric = numeric.transpose()

    numeric = numeric.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.dropna(axis=0, how="all").dropna(axis=1, how="all")
    numeric.columns = numeric.columns.astype(str).str.strip()
    return numeric


def load_meta(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path)
    columns_lower = {column.lower().strip(): column for column in meta.columns}

    sample_column = None
    for candidate in SAMPLE_COLUMN_CANDIDATES:
        if candidate in columns_lower:
            sample_column = columns_lower[candidate]
            break

    if sample_column is None:
        if meta.index.name is not None:
            meta.index = meta.index.astype(str).str.strip()
            return meta
        meta = pd.read_csv(path, index_col=0)
        meta.index = meta.index.astype(str).str.strip()
        return meta

    meta = meta.set_index(sample_column)
    meta.index = meta.index.astype(str).str.strip()
    meta.index.name = "sample"
    return meta


def align_samples(counts: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    common = [sample for sample in counts.columns if sample in meta.index]
    if not common:
        raise ValueError("No overlapping sample IDs found between counts and metadata.")
    return counts.loc[:, common], meta.loc[common].copy()


def log2_transform(counts: pd.DataFrame) -> pd.DataFrame:
    return np.log2(counts.astype(float) + 1.0)


def select_variable_genes(log_counts: pd.DataFrame, top_genes: int, min_detection: int) -> pd.DataFrame:
    detected = (log_counts > 0).sum(axis=1) >= min_detection
    filtered = log_counts.loc[detected].copy()
    if filtered.empty:
        filtered = log_counts.copy()

    variances = filtered.var(axis=1).sort_values(ascending=False)
    if top_genes and top_genes > 0:
        variances = variances.head(min(top_genes, len(variances)))

    selected = filtered.loc[variances.index]
    return selected


def infer_group_column(meta: pd.DataFrame) -> str | None:
    lowered = {column.lower().strip(): column for column in meta.columns}
    for candidate in TISSUE_COLUMN_CANDIDATES:
        if candidate in lowered:
            return lowered[candidate]
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot hierarchical clustering heatmap from raw RNA-seq counts.")
    parser.add_argument("--counts", required=True, type=Path, help="Raw count matrix CSV")
    parser.add_argument("--meta", required=True, type=Path, help="Sample metadata CSV")
    parser.add_argument("--outdir", default=Path("plots_heatmap"), type=Path, help="Output directory")
    parser.add_argument("--top-genes", type=int, default=500, help="Top variable genes to plot")
    parser.add_argument("--min-detection", type=int, default=3, help="Minimum number of samples with non-zero counts")
    parser.add_argument("--figsize-width", type=float, default=10, help="Figure width in inches")
    parser.add_argument("--figsize-height", type=float, default=12, help="Figure height in inches")
    parser.add_argument("--cmap", default="RdYlBu_r", help="Colormap name")
    args = parser.parse_args()

    sns.set_theme(style="white")

    counts = load_counts(args.counts)
    meta = load_meta(args.meta)
    counts, meta = align_samples(counts, meta)

    group_column = infer_group_column(meta)

    log_counts = log2_transform(counts)
    selected = select_variable_genes(log_counts, top_genes=args.top_genes, min_detection=args.min_detection)

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Create color map for sample groups
    group_colors = None
    if group_column and group_column in meta.columns:
        tissue_palette = {"tumor": "#B2182B", "normal": "#2166AC"}
        group_colors = pd.Series(
            [tissue_palette.get(str(val).strip().lower(), "#808080") for val in meta[group_column]],
            index=meta.index,
        )

    # Create heatmap with dendrogram
    g = sns.clustermap(
        selected.T,
        method="average",
        metric="euclidean",
        cmap=args.cmap,
        figsize=(args.figsize_width, args.figsize_height),
        col_colors=None,
        row_colors=group_colors,
        cbar_kws={"label": "log2(counts + 1)"},
        dendrogram_ratio=0.15,
        linewidths=0.5,
        linecolor="white",
    )

    g.ax_heatmap.set_xlabel("Genes", fontsize=12, weight="bold")
    g.ax_heatmap.set_ylabel("Samples", fontsize=12, weight="bold")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9)
    g.ax_heatmap.set_xticklabels([], fontsize=8, rotation=90)

    # Add legend for groups if available
    if group_column and group_column in meta.columns:
        tissue_palette = {"tumor": "#B2182B", "normal": "#2166AC"}
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=tissue_palette.get("tumor", "#333"), markersize=10, label="tumor"),
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=tissue_palette.get("normal", "#333"), markersize=10, label="normal"),
        ]
        g.ax_heatmap.legend(handles=handles, title=group_column.capitalize(), loc="upper left", bbox_to_anchor=(1.25, 1.0), frameon=True)

    g.fig.suptitle(f"Hierarchical Clustering Heatmap - Top {args.top_genes} Variable Genes", fontsize=14, weight="bold", y=0.98)

    g.fig.savefig(args.outdir / "raw_counts_heatmap.png", dpi=300, bbox_inches="tight")
    g.fig.savefig(args.outdir / "raw_counts_heatmap.pdf", bbox_inches="tight")
    plt.close(g.fig)

    print(f"Heatmap saved to {args.outdir}")


if __name__ == "__main__":
    main()
