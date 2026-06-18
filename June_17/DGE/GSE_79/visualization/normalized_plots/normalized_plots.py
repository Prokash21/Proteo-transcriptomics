"""Create publication-ready PCA, UMAP, hierarchical tree, and boxplot from normalized counts.

This script uses the normalized expression matrix directly and writes all outputs into a
separate folder so raw and normalized figures stay organized.

Example:
python normalized_plots.py \
  --counts ../../Normalized_Count_Data_wol_03.csv \
  --meta ../raw_data/meta_data_GSE142279.csv \
  --outdir ./plots_normalized
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
import umap
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


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


def select_variable_genes(matrix: pd.DataFrame, top_genes: int, min_detection: int) -> pd.DataFrame:
    detected = (matrix > 0).sum(axis=1) >= min_detection
    filtered = matrix.loc[detected].copy()
    if filtered.empty:
        filtered = matrix.copy()

    variances = filtered.var(axis=1).sort_values(ascending=False)
    if top_genes and top_genes > 0:
        variances = variances.head(min(top_genes, len(variances)))

    return filtered.loc[variances.index]


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    scaler = StandardScaler(with_mean=True, with_std=True)
    scaled = scaler.fit_transform(matrix.T)
    return pd.DataFrame(scaled, index=matrix.columns, columns=matrix.index)


def make_scatter(frame: pd.DataFrame, meta: pd.DataFrame, x: str, y: str, out_base: Path, title: str, xlabel: str, ylabel: str, group_column: str | None, marker_size: int = 70) -> None:
    palette = {"tumor": "#B2182B", "normal": "#2166AC"}
    fig, ax = plt.subplots(figsize=(5.6, 4.8), constrained_layout=True)

    if group_column and group_column in meta.columns:
        plot_data = frame.join(meta[[group_column]])
        sns.scatterplot(
            data=plot_data,
            x=x,
            y=y,
            hue=group_column,
            palette=palette,
            s=marker_size,
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
            s=marker_size,
            color="#333333",
            edgecolor="white",
            linewidth=0.8,
            ax=ax,
        )

    ax.set_title(title, fontsize=13, weight="bold")
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)

    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_boxplot(counts: pd.DataFrame, meta: pd.DataFrame, out_base: Path, group_column: str | None) -> None:
    log_counts = np.log2(counts.astype(float) + 1.0)
    palette = {"tumor": "#B2182B", "normal": "#2166AC"}

    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    if group_column and group_column in meta.columns:
        sample_order = meta.sort_values(group_column).index.tolist()
    else:
        sample_order = log_counts.columns.tolist()

    boxplot_data_samples = []
    for sample in sample_order:
        sample_values = log_counts[sample].values
        for value in sample_values:
            group = meta.loc[sample, group_column] if group_column and group_column in meta.columns else "Unknown"
            boxplot_data_samples.append({
                "Sample": sample,
                "Value": value,
                group_column if group_column else "Group": group,
            })

    df_plot_samples = pd.DataFrame(boxplot_data_samples)

    sample_palette = {}
    for sample in sample_order:
        group = meta.loc[sample, group_column] if group_column and group_column in meta.columns else "Unknown"
        sample_palette[sample] = palette.get(str(group).strip().lower(), "#808080")

    sns.boxplot(
        data=df_plot_samples,
        x="Sample",
        y="Value",
        palette=sample_palette,
        ax=ax,
        width=0.7,
    )

    ax.set_xlabel("Samples", fontsize=11, weight="bold")
    ax.set_ylabel("log2(normalized + 1)", fontsize=11, weight="bold")
    ax.set_title("Expression Distribution per Sample", fontsize=12, weight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_dendrogram(counts: pd.DataFrame, meta: pd.DataFrame, out_base: Path, group_column: str | None, top_genes: int, min_detection: int) -> None:
    log_counts = np.log2(counts.astype(float) + 1.0)
    selected = select_variable_genes(log_counts, top_genes=top_genes, min_detection=min_detection)
    scaled = zscore_rows(selected)

    sample_dist = pdist(scaled.T, metric="euclidean")
    Z = linkage(sample_dist, method="average")

    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    dendrogram(Z, labels=scaled.columns, leaf_rotation=90, leaf_font_size=10, ax=ax)

    if group_column and group_column in meta.columns:
        tissue_palette = {"tumor": "#B2182B", "normal": "#2166AC"}
        for label in ax.get_xticklabels():
            sample_id = label.get_text()
            if sample_id in meta.index:
                tissue = str(meta.loc[sample_id, group_column]).strip().lower()
                label.set_color(tissue_palette.get(tissue, "#333333"))
                label.set_weight("bold")

        handles = [
            plt.Line2D([0], [0], color="#B2182B", lw=3, label="tumor"),
            plt.Line2D([0], [0], color="#2166AC", lw=3, label="normal"),
        ]
        ax.legend(handles=handles, title=group_column.capitalize(), loc="upper right", frameon=True, fontsize=11)

    ax.set_title("Hierarchical Clustering of Samples (Normalized Counts)", fontsize=14, weight="bold")
    ax.set_xlabel("Samples", fontsize=12, weight="bold")
    ax.set_ylabel("Distance", fontsize=12, weight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot PCA, UMAP, dendrogram, and boxplot from normalized counts.")
    parser.add_argument("--counts", required=True, type=Path, help="Normalized count matrix CSV")
    parser.add_argument("--meta", required=True, type=Path, help="Sample metadata CSV")
    parser.add_argument("--outdir", default=Path("plots_normalized"), type=Path, help="Output directory")
    parser.add_argument("--top-genes", type=int, default=2000, help="Top variable genes to keep before embedding")
    parser.add_argument("--min-detection", type=int, default=3, help="Minimum number of samples with non-zero values")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed for UMAP")
    args = parser.parse_args()

    sns.set_theme(style="white", context="paper", font_scale=1.1)

    counts = load_counts(args.counts)
    meta = load_meta(args.meta)
    counts, meta = align_samples(counts, meta)
    group_column = infer_group_column(meta)

    selected = select_variable_genes(counts, top_genes=args.top_genes, min_detection=args.min_detection)
    scaled = zscore_rows(selected)

    pca_n_components = min(10, scaled.shape[0], scaled.shape[1])
    if pca_n_components < 2:
        raise ValueError("Not enough samples or genes to compute a 2D PCA embedding.")

    pca = PCA(n_components=pca_n_components, random_state=args.random_state)
    pca_scores = pca.fit_transform(scaled)
    pca_frame = pd.DataFrame(pca_scores[:, :2], index=scaled.index, columns=["PC1", "PC2"])
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

    args.outdir.mkdir(parents=True, exist_ok=True)

    make_scatter(
        pca_frame,
        meta,
        x="PC1",
        y="PC2",
        out_base=args.outdir / "normalized_pca",
        title="PCA after normalization",
        xlabel=f"PC1 ({pca_var[0]:.1f}%)",
        ylabel=f"PC2 ({pca_var[1]:.1f}%)",
        group_column=group_column,
        marker_size=68,
    )
    make_scatter(
        umap_frame,
        meta,
        x="UMAP1",
        y="UMAP2",
        out_base=args.outdir / "normalized_umap",
        title="UMAP after normalization",
        xlabel="UMAP1",
        ylabel="UMAP2",
        group_column=group_column,
        marker_size=68,
    )
    make_dendrogram(
        counts,
        meta,
        out_base=args.outdir / "normalized_dendrogram",
        group_column=group_column,
        top_genes=args.top_genes,
        min_detection=args.min_detection,
    )
    make_boxplot(counts, meta, out_base=args.outdir / "normalized_boxplot_per_sample", group_column=group_column)

    coords = pca_frame.join(umap_frame)
    coords = coords.join(meta)
    coords.to_csv(args.outdir / "normalized_coordinates.csv")

    print(f"Normalized plots saved to {args.outdir}")


if __name__ == "__main__":
    main()
