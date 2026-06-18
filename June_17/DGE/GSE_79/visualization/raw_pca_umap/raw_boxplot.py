"""Create publication-ready boxplot of expression distribution from raw counts.

Shows the distribution of log-transformed counts across samples, grouped by tissue type.

Example:
python raw_boxplot.py \
  --counts ../raw_data/count_data_raw_GSE142279.csv \
  --meta ../raw_data/meta_data_GSE142279.csv \
  --outdir ./plots_boxplot
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


def infer_group_column(meta: pd.DataFrame) -> str | None:
    lowered = {column.lower().strip(): column for column in meta.columns}
    for candidate in TISSUE_COLUMN_CANDIDATES:
        if candidate in lowered:
            return lowered[candidate]
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot boxplot of expression distribution from raw RNA-seq counts.")
    parser.add_argument("--counts", required=True, type=Path, help="Raw count matrix CSV")
    parser.add_argument("--meta", required=True, type=Path, help="Sample metadata CSV")
    parser.add_argument("--outdir", default=Path("plots_boxplot"), type=Path, help="Output directory")
    parser.add_argument("--figsize-width", type=float, default=5.5, help="Figure width in inches")
    parser.add_argument("--figsize-height", type=float, default=4.5, help="Figure height in inches")
    args = parser.parse_args()

    sns.set_theme(style="white", context="paper", font_scale=1.0)

    counts = load_counts(args.counts)
    meta = load_meta(args.meta)
    counts, meta = align_samples(counts, meta)

    group_column = infer_group_column(meta)

    # Log-transform counts
    log_counts = log2_transform(counts)

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Prepare data for boxplot: long format
    boxplot_data = []
    for sample in log_counts.columns:
        sample_values = log_counts[sample].values
        for value in sample_values:
            group = meta.loc[sample, group_column] if group_column and group_column in meta.columns else "Unknown"
            boxplot_data.append({
                "Sample": sample,
                "Value": value,
                group_column if group_column else "Group": group,
            })

    df_plot = pd.DataFrame(boxplot_data)

    # Create boxplot
    fig, ax = plt.subplots(figsize=(args.figsize_width, args.figsize_height), constrained_layout=True)

    palette = {"tumor": "#B2182B", "normal": "#2166AC"}

    if group_column and group_column in meta.columns:
        sns.boxplot(
            data=df_plot,
            x=group_column,
            y="Value",
            palette=palette,
            ax=ax,
            width=0.6,
        )
    else:
        sns.boxplot(
            data=df_plot,
            y="Value",
            ax=ax,
            color="#2166AC",
            width=0.3,
        )

    ax.set_xlabel(group_column.capitalize() if group_column else "Samples", fontsize=11, weight="bold")
    ax.set_ylabel("log2(counts + 1)", fontsize=11, weight="bold")
    ax.set_title("Expression Distribution by Sample Group", fontsize=12, weight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    fig.savefig(args.outdir / "boxplot_by_group.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.outdir / "boxplot_by_group.pdf", bbox_inches="tight")

    # Create second boxplot: individual samples
    fig, ax = plt.subplots(figsize=(args.figsize_width + 1.5, args.figsize_height), constrained_layout=True)

    # Reorder by group for better visualization
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

    # Create hue palette for individual samples
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
    ax.set_ylabel("log2(counts + 1)", fontsize=11, weight="bold")
    ax.set_title("Expression Distribution per Sample", fontsize=12, weight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    fig.savefig(args.outdir / "boxplot_per_sample.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.outdir / "boxplot_per_sample.pdf", bbox_inches="tight")

    plt.close("all")
    print(f"Boxplots saved to {args.outdir}")


if __name__ == "__main__":
    main()
