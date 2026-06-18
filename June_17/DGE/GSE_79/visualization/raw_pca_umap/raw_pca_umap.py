"""Create publication-ready PCA and UMAP plots from raw RNA-seq counts.

The script keeps the input counts unnormalized for plotting purposes. A log2(x+1)
transform is applied only to stabilize the variance before PCA/UMAP.

Example:
python raw_pca_umap.py \
  --counts ../raw_data/count_data_raw_GSE142279.csv \
  --meta ../raw_data/meta_data_GSE142279.csv \
  --outdir ./plots_raw_counts
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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap


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


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    scaler = StandardScaler(with_mean=True, with_std=True)
    scaled = scaler.fit_transform(matrix.T)
    return pd.DataFrame(scaled, index=matrix.columns, columns=matrix.index)


def infer_group_column(meta: pd.DataFrame) -> str | None:
    lowered = {column.lower().strip(): column for column in meta.columns}
    for candidate in TISSUE_COLUMN_CANDIDATES:
        if candidate in lowered:
            return lowered[candidate]
    return None


def style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)


def scatter_embedding(
    frame: pd.DataFrame,
    meta: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    xlabel: str,
    ylabel: str,
    out_base: Path,
    group_column: str | None,
    label_points: bool,
) -> None:
    palette = {"tumor": "#B2182B", "normal": "#2166AC"}

    fig, ax = plt.subplots(figsize=(6.2, 5.3), constrained_layout=True)
    if group_column and group_column in meta.columns:
        plot_data = frame.join(meta[[group_column]])
        sns.scatterplot(
            data=plot_data,
            x=x,
            y=y,
            hue=group_column,
            palette=palette,
            s=70,
            edgecolor="white",
            linewidth=0.8,
            ax=ax,
        )
        ax.legend(title=group_column.capitalize(), frameon=True, loc="best", edgecolor="black", fancybox=True)
    else:
        sns.scatterplot(
            data=frame,
            x=x,
            y=y,
            s=70,
            color="#333333",
            edgecolor="white",
            linewidth=0.8,
            ax=ax,
        )

    if label_points:
        for sample, row in frame.iterrows():
            ax.text(row[x], row[y], sample, fontsize=7, alpha=0.75, ha="left", va="bottom")

    ax.set_title(title, fontsize=14, weight="bold")
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    style_axes(ax)

    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot PCA and UMAP from raw RNA-seq counts.")
    parser.add_argument("--counts", required=True, type=Path, help="Raw count matrix CSV")
    parser.add_argument("--meta", required=True, type=Path, help="Sample metadata CSV")
    parser.add_argument("--outdir", default=Path("plots_raw_counts"), type=Path, help="Output directory")
    parser.add_argument("--top-genes", type=int, default=2000, help="Top variable genes to keep before embedding")
    parser.add_argument("--min-detection", type=int, default=3, help="Minimum number of samples with non-zero counts")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed for UMAP")
    parser.add_argument("--label-points", action="store_true", help="Annotate each sample point with its sample ID")
    args = parser.parse_args()

    sns.set_theme(style="white", context="paper", font_scale=1.2)

    counts = load_counts(args.counts)
    meta = load_meta(args.meta)
    counts, meta = align_samples(counts, meta)

    group_column = infer_group_column(meta)

    log_counts = log2_transform(counts)
    selected = select_variable_genes(log_counts, top_genes=args.top_genes, min_detection=args.min_detection)
    scaled = zscore_rows(selected)

    pca_n_components = min(10, scaled.shape[0], scaled.shape[1])
    if pca_n_components < 2:
        raise ValueError("Not enough samples or genes to compute a 2D PCA embedding.")

    pca = PCA(n_components=pca_n_components, random_state=args.random_state)
    pca_scores = pca.fit_transform(scaled)
    pca_frame = pd.DataFrame(pca_scores[:, :2], index=scaled.index, columns=["PC1", "PC2"])
    pca_frame["sample"] = pca_frame.index
    pca_var = pca.explained_variance_ratio_ * 100

    umap_input = pca_scores[:, : min(15, pca_scores.shape[1])]
    umap_model = umap.UMAP(
        n_neighbors=min(15, max(2, scaled.shape[0] - 1)),
        min_dist=0.25,
        metric="euclidean",
        random_state=args.random_state,
    )
    umap_scores = umap_model.fit_transform(umap_input)
    umap_frame = pd.DataFrame(umap_scores, index=scaled.index, columns=["UMAP1", "UMAP2"])
    umap_frame["sample"] = umap_frame.index

    args.outdir.mkdir(parents=True, exist_ok=True)

    scatter_embedding(
        pca_frame,
        meta,
        x="PC1",
        y="PC2",
        title="PCA of raw counts",
        xlabel=f"PC1 ({pca_var[0]:.1f}%)",
        ylabel=f"PC2 ({pca_var[1]:.1f}%)",
        out_base=args.outdir / "raw_counts_pca",
        group_column=group_column,
        label_points=args.label_points,
    )

    scatter_embedding(
        umap_frame,
        meta,
        x="UMAP1",
        y="UMAP2",
        title="UMAP of raw counts",
        xlabel="UMAP1",
        ylabel="UMAP2",
        out_base=args.outdir / "raw_counts_umap",
        group_column=group_column,
        label_points=args.label_points,
    )

    coords = pca_frame[["sample", "PC1", "PC2"]].join(umap_frame[["UMAP1", "UMAP2"]])
    coords = coords.join(meta)
    coords.to_csv(args.outdir / "raw_counts_pca_umap_coordinates.csv")


if __name__ == "__main__":
    main()