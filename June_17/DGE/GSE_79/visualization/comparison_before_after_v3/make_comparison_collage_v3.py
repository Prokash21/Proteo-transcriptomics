"""Build a fully vector before/after collage for Illustrator.

This version does not stitch raster images. It draws the before/after panels directly
with matplotlib and saves the collage as PDF and SVG so the whole page stays editable
in Adobe Illustrator.

Layout:
UMAP before | UMAP after
Dendrogram before | Dendrogram after
Boxplot before | Boxplot after
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
PAGE_PRESETS = {
    "A4_landscape": (11.69, 8.27),
    "A4_portrait": (8.27, 11.69),
    "letter_landscape": (11.0, 8.5),
    "letter_portrait": (8.5, 11.0),
}


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


def style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)


def add_group_legend(ax: plt.Axes, title: str = "Tissue") -> None:
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#B2182B", markeredgecolor="white", markersize=8, label="tumor"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#2166AC", markeredgecolor="white", markersize=8, label="normal"),
    ]
    ax.legend(handles=handles, title=title, loc="best", frameon=True, edgecolor="black", fancybox=True, fontsize=8, title_fontsize=9)


def prepare_group_colors(meta: pd.DataFrame, group_column: str | None) -> pd.Series | None:
    if not group_column or group_column not in meta.columns:
        return None
    palette = {"tumor": "#B2182B", "normal": "#2166AC"}
    values = [palette.get(str(value).strip().lower(), "#808080") for value in meta[group_column]]
    return pd.Series(values, index=meta.index)


def variable_gene_matrix(mat: pd.DataFrame, top_genes: int = 2000, min_detection: int = 3) -> pd.DataFrame:
    log_mat = np.log2(mat.astype(float) + 1.0)
    detected = (log_mat > 0).sum(axis=1) >= min_detection
    filtered = log_mat.loc[detected].copy()
    if filtered.empty:
        filtered = log_mat.copy()
    variances = filtered.var(axis=1).sort_values(ascending=False)
    if top_genes and top_genes > 0:
        variances = variances.head(min(top_genes, len(variances)))
    return filtered.loc[variances.index]


def compute_umap_embedding(mat: pd.DataFrame, random_state: int) -> pd.DataFrame:
    selected = variable_gene_matrix(mat)
    scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(selected.T)
    pca_n_components = min(10, scaled.shape[0], scaled.shape[1])
    if pca_n_components < 2:
        raise ValueError("Not enough samples or genes to compute a 2D UMAP embedding.")
    pca = PCA(n_components=pca_n_components, random_state=random_state)
    pca_scores = pca.fit_transform(scaled)
    umap_input = pca_scores[:, : min(15, pca_scores.shape[1])]
    umap_model = umap.UMAP(
        n_neighbors=min(15, max(2, scaled.shape[0] - 1)),
        min_dist=0.25,
        metric="euclidean",
        random_state=random_state,
    )
    coords = umap_model.fit_transform(umap_input)
    return pd.DataFrame(coords, index=selected.columns, columns=["UMAP1", "UMAP2"])


def plot_umap_ax(mat: pd.DataFrame, meta: pd.DataFrame, ax: plt.Axes, title: str, random_state: int) -> None:
    emb = compute_umap_embedding(mat, random_state=random_state)
    group_column = infer_group_column(meta)

    if group_column and group_column in meta.columns:
        plot_data = emb.join(meta[[group_column]])
        sns.scatterplot(
            data=plot_data,
            x="UMAP1",
            y="UMAP2",
            hue=group_column,
            palette={"tumor": "#B2182B", "normal": "#2166AC"},
            s=48,
            edgecolor="white",
            linewidth=0.6,
            ax=ax,
        )
        add_group_legend(ax, title=group_column.capitalize())
    else:
        sns.scatterplot(
            data=emb,
            x="UMAP1",
            y="UMAP2",
            s=48,
            color="#333333",
            edgecolor="white",
            linewidth=0.6,
            ax=ax,
        )

    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_xlabel("UMAP1", fontsize=9)
    ax.set_ylabel("UMAP2", fontsize=9)
    style_axes(ax)


def plot_dendrogram_ax(mat: pd.DataFrame, meta: pd.DataFrame, ax: plt.Axes, title: str) -> None:
    group_column = infer_group_column(meta)
    sample_dist = pdist(mat.T, metric="euclidean")
    Z = linkage(sample_dist, method="average")
    dendrogram(Z, labels=mat.columns, leaf_rotation=90, leaf_font_size=7, ax=ax)

    if group_column and group_column in meta.columns:
        tissue_colors = {"tumor": "#d62728", "normal": "#1f77b4"}
        for tick in ax.get_xmajorticklabels():
            sample = tick.get_text()
            tissue = str(meta.loc[sample, group_column]).strip().lower() if sample in meta.index else ""
            tick.set_color(tissue_colors.get(tissue, "#333333"))
            tick.set_weight("bold")
        handles = [
            plt.Line2D([0], [0], color="#d62728", lw=2.5, label="tumor"),
            plt.Line2D([0], [0], color="#1f77b4", lw=2.5, label="normal"),
        ]
        ax.legend(handles=handles, title=group_column.capitalize(), loc="best", frameon=True, edgecolor="black", fancybox=True, fontsize=8, title_fontsize=9)

    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_xlabel("Samples", fontsize=9, weight="bold")
    ax.set_ylabel("Distance", fontsize=9, weight="bold")
    style_axes(ax)


def plot_boxplot_ax(mat: pd.DataFrame, meta: pd.DataFrame, ax: plt.Axes, title: str) -> None:
    group_column = infer_group_column(meta)
    logmat = np.log2(mat.astype(float) + 1.0)
    palette = {"tumor": "#B2182B", "normal": "#2166AC"}

    if group_column and group_column in meta.columns:
        sample_order = meta.sort_values(group_column).index.tolist()
    else:
        sample_order = logmat.columns.tolist()

    boxplot_data = []
    for sample in sample_order:
        for value in logmat[sample].values:
            group = meta.loc[sample, group_column] if group_column and group_column in meta.columns else "Unknown"
            boxplot_data.append({"Sample": sample, "Value": value, "Group": group})

    df_plot = pd.DataFrame(boxplot_data)
    sample_palette = {}
    for sample in sample_order:
        group = meta.loc[sample, group_column] if group_column and group_column in meta.columns else "Unknown"
        sample_palette[sample] = palette.get(str(group).strip().lower(), "#808080")

    sns.boxplot(
        data=df_plot,
        x="Sample",
        y="Value",
        palette=sample_palette,
        ax=ax,
        width=0.7,
        linewidth=0.8,
    )

    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_xlabel("Samples", fontsize=9, weight="bold")
    ax.set_ylabel("log2(expression + 1)", fontsize=9, weight="bold")
    ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    style_axes(ax)
    ax.grid(axis="y", alpha=0.25, linestyle="--")


def save_panel(fig: plt.Figure, outdir: Path, name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=300, bbox_inches="tight")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a vectorized before/after comparison collage.")
    parser.add_argument("--raw-counts", required=True, type=Path, help="Raw count matrix CSV")
    parser.add_argument("--normalized-counts", required=True, type=Path, help="Normalized count matrix CSV")
    parser.add_argument("--meta", required=True, type=Path, help="Sample metadata CSV")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--page", choices=list(PAGE_PRESETS.keys()), default="A4_landscape", help="Figure page size preset")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed for UMAP")
    args = parser.parse_args()

    sns.set_theme(style="white", context="paper", font_scale=1.0)

    raw_counts = load_counts(args.raw_counts)
    norm_counts = load_counts(args.normalized_counts)
    meta = load_meta(args.meta)
    raw_counts, meta = align_samples(raw_counts, meta)
    norm_counts, meta = align_samples(norm_counts, meta)

    before_dir = args.outdir / "before"
    after_dir = args.outdir / "after"
    args.outdir.mkdir(parents=True, exist_ok=True)
    before_dir.mkdir(parents=True, exist_ok=True)
    after_dir.mkdir(parents=True, exist_ok=True)

    width, height = PAGE_PRESETS[args.page]
    fig, axes = plt.subplots(3, 2, figsize=(width, height), constrained_layout=True)

    plot_umap_ax(raw_counts, meta, axes[0, 0], "UMAP before normalization", random_state=args.random_state)
    plot_umap_ax(norm_counts, meta, axes[0, 1], "UMAP after normalization", random_state=args.random_state)

    plot_dendrogram_ax(raw_counts, meta, axes[1, 0], "Dendrogram before normalization")
    plot_dendrogram_ax(norm_counts, meta, axes[1, 1], "Dendrogram after normalization")

    plot_boxplot_ax(raw_counts, meta, axes[2, 0], "Boxplot before normalization")
    plot_boxplot_ax(norm_counts, meta, axes[2, 1], "Boxplot after normalization")

    fig.suptitle("Raw vs Normalized Sample-Level Plots", fontsize=14, weight="bold")
    fig.savefig(args.outdir / "comparison_collage_v3.pdf", bbox_inches="tight")
    fig.savefig(args.outdir / "comparison_collage_v3.svg", bbox_inches="tight")
    fig.savefig(args.outdir / "comparison_collage_v3.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save individual vector panels too.
    panel_specs = [
        (raw_counts, "before", "umap_before", plot_umap_ax, "UMAP before normalization", {"random_state": args.random_state}),
        (norm_counts, "after", "umap_after", plot_umap_ax, "UMAP after normalization", {"random_state": args.random_state}),
        (raw_counts, "before", "dendrogram_before", plot_dendrogram_ax, "Dendrogram before normalization", {}),
        (norm_counts, "after", "dendrogram_after", plot_dendrogram_ax, "Dendrogram after normalization", {}),
        (raw_counts, "before", "boxplot_before", plot_boxplot_ax, "Boxplot before normalization", {}),
        (norm_counts, "after", "boxplot_after", plot_boxplot_ax, "Boxplot after normalization", {}),
    ]
    for mat, folder, name, plot_func, title, kwargs in panel_specs:
        fig_panel, ax_panel = plt.subplots(figsize=(5.0, 4.0), constrained_layout=True)
        plot_func(mat, meta, ax_panel, title, **kwargs)
        save_panel(fig_panel, args.outdir / folder, name)
        plt.close(fig_panel)

    print(f"Vector collage saved to {args.outdir}")


if __name__ == "__main__":
    main()
