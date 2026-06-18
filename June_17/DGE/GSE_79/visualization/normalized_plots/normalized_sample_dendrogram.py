"""Create a simple hierarchical clustering dendrogram of samples from normalized counts.

Matches the sample-only clustering logic used by gse79_plot_samples.py, while keeping
publication-style label coloring and clean figure styling.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


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


def infer_group_column(meta: pd.DataFrame) -> str | None:
    lowered = {column.lower().strip(): column for column in meta.columns}
    for candidate in TISSUE_COLUMN_CANDIDATES:
        if candidate in lowered:
            return lowered[candidate]
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot hierarchical dendrogram of samples from normalized counts.")
    parser.add_argument("--counts", required=True, type=Path, help="Normalized count matrix CSV")
    parser.add_argument("--meta", required=True, type=Path, help="Sample metadata CSV")
    parser.add_argument("--outdir", default=Path("plots_normalized_dendrogram"), type=Path, help="Output directory")
    parser.add_argument("--figsize-width", type=float, default=10, help="Figure width in inches")
    parser.add_argument("--figsize-height", type=float, default=6, help="Figure height in inches")
    args = parser.parse_args()

    sns.set_theme(style="white")

    counts = load_counts(args.counts)
    meta = load_meta(args.meta)
    counts, meta = align_samples(counts, meta)

    group_column = infer_group_column(meta)

    args.outdir.mkdir(parents=True, exist_ok=True)

    sample_dist = pdist(counts.T, metric="euclidean")
    Z = linkage(sample_dist, method="average")

    fig, ax = plt.subplots(figsize=(args.figsize_width, args.figsize_height))
    dendrogram(Z, labels=counts.columns, leaf_rotation=90, ax=ax, leaf_font_size=10)

    if group_column and group_column in meta.columns:
        tissue_colors = {"tumor": "#d62728", "normal": "#1f77b4"}
        for tick in ax.get_xmajorticklabels():
            sample = tick.get_text()
            tissue = str(meta.loc[sample, group_column]).strip().lower() if sample in meta.index else ""
            tick.set_color(tissue_colors.get(tissue, "#333333"))

        handles = [
            plt.Line2D([0], [0], color="#d62728", lw=3, label="tumor"),
            plt.Line2D([0], [0], color="#1f77b4", lw=3, label="normal"),
        ]
        ax.legend(handles=handles, title="Tissue", loc="upper right", frameon=True, edgecolor="black", fancybox=True, fontsize=11)

    ax.set_title("Hierarchical Clustering of Samples (Normalized Counts)", fontsize=14, weight="bold")
    ax.set_xlabel("Samples", fontsize=12, weight="bold")
    ax.set_ylabel("Distance", fontsize=12, weight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(args.outdir / "normalized_sample_dendrogram.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.outdir / "normalized_sample_dendrogram.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"Dendrogram saved to {args.outdir}")


if __name__ == "__main__":
    main()
