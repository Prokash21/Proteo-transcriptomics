from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LFC_THRESHOLD = 1.0
PADJ_THRESHOLD = 0.05

DATASETS = {
    "GSE180440": "mapped_GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal_mygene.csv",
    "GSE251845": "mapped_GSE_45_OR_02_resLFC_condition_tumor_vs_normal_mygene.csv",
    "GSE50760": "mapped_GSE_60_resLFC_raw_wolTreatment_tumor_vs_normal_mygene.csv",
    "GSE142279": "mapped_GSE_79_02_resLFC_Treatment_tumor_vs_normal_mygene.csv",
    "TCGA COAD": "mapped_tcga_coad_resLFC_Treatment_Tumor_vs_Normal_mygene.csv",
    "TCGA READ": "mapped_tcga_read_resLFC_SampleType_Tumor_vs_Normal_mygene.csv",
}

COLORS = {
    "Not significant": "#B8B8B8",
    "Up": "#D84A4A",
    "Down": "#2F6B9A",
}


def prepare_volcano_data(path: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    df = pd.read_csv(path)
    required_columns = {"gene_id", "log2FoldChange", "padj"}
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")

    df = df.copy()
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df = df.dropna(subset=["log2FoldChange", "padj"])
    df = df[df["padj"] > 0]
    df["neg_log10_padj"] = -np.log10(df["padj"])

    df["status"] = "Not significant"
    df.loc[
        (df["log2FoldChange"] >= LFC_THRESHOLD) & (df["padj"] < PADJ_THRESHOLD),
        "status",
    ] = "Up"
    df.loc[
        (df["log2FoldChange"] <= -LFC_THRESHOLD) & (df["padj"] < PADJ_THRESHOLD),
        "status",
    ] = "Down"

    counts = {
        "total_tested": int(len(df)),
        "up": int((df["status"] == "Up").sum()),
        "down": int((df["status"] == "Down").sum()),
        "significant": int((df["status"].isin(["Up", "Down"])).sum()),
        "not_significant": int((df["status"] == "Not significant").sum()),
    }
    return df, counts


def top_genes_for_labels(df: pd.DataFrame, n: int = 8) -> pd.DataFrame:
    significant = df[df["status"].isin(["Up", "Down"])].copy()
    significant = significant.dropna(subset=["gene_id"])
    significant["gene_id"] = significant["gene_id"].astype(str).str.strip()
    significant = significant[significant["gene_id"] != ""]
    significant = significant.sort_values(["padj", "log2FoldChange"], ascending=[True, False])
    return significant.head(n)


def plot_volcano(
    df: pd.DataFrame,
    dataset_name: str,
    counts: dict[str, int],
    output_path: Path,
    annotate: bool = True,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#FBFBFA")

    for status in ["Not significant", "Down", "Up"]:
        subset = df[df["status"] == status]
        ax.scatter(
            subset["log2FoldChange"],
            subset["neg_log10_padj"],
            s=8 if status == "Not significant" else 13,
            c=COLORS[status],
            alpha=0.38 if status == "Not significant" else 0.72,
            linewidths=0,
            label=status,
            rasterized=True,
        )

    ax.axvline(-LFC_THRESHOLD, color="#444444", linestyle="--", linewidth=1)
    ax.axvline(LFC_THRESHOLD, color="#444444", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(PADJ_THRESHOLD), color="#444444", linestyle="--", linewidth=1)

    if annotate:
        y_limit = df["neg_log10_padj"].quantile(0.995)
        for _, row in top_genes_for_labels(df).iterrows():
            if row["neg_log10_padj"] <= y_limit * 0.25:
                continue
            ax.text(
                row["log2FoldChange"],
                row["neg_log10_padj"],
                str(row["gene_id"])[:18],
                fontsize=7,
                ha="center",
                va="bottom",
                color="#222222",
            )

    x_abs = float(np.nanpercentile(np.abs(df["log2FoldChange"]), 99.5))
    x_limit = max(3.0, min(max(x_abs * 1.1, LFC_THRESHOLD + 0.6), 12.0))
    ax.set_xlim(-x_limit, x_limit)
    ax.set_ylim(0, max(3.0, float(df["neg_log10_padj"].quantile(0.995)) * 1.08))

    ax.grid(alpha=0.18, linewidth=0.7)
    ax.set_xlabel("log2 Fold Change", fontsize=10)
    ax.set_ylabel("-log10 adjusted p-value", fontsize=10)
    ax.set_title(dataset_name, fontsize=14, fontweight="bold", pad=10)
    ax.text(
        0.02,
        0.98,
        f"Up: {counts['up']:,}\nDown: {counts['down']:,}\nThreshold: |LFC| ≥ {LFC_THRESHOLD:g}, padj < {PADJ_THRESHOLD:g}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#D0D0D0", "alpha": 0.94},
    )
    ax.legend(loc="upper right", frameon=True, fontsize=8, markerscale=1.5)

    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.12, top=0.88)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def save_collage(plot_data: list[tuple[str, pd.DataFrame, dict[str, int]]], output_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.patch.set_facecolor("white")

    for ax, (dataset_name, df, counts) in zip(axes.ravel(), plot_data):
        ax.set_facecolor("#FBFBFA")
        for status in ["Not significant", "Down", "Up"]:
            subset = df[df["status"] == status]
            ax.scatter(
                subset["log2FoldChange"],
                subset["neg_log10_padj"],
                s=5 if status == "Not significant" else 9,
                c=COLORS[status],
                alpha=0.30 if status == "Not significant" else 0.70,
                linewidths=0,
                rasterized=True,
            )

        ax.axvline(-LFC_THRESHOLD, color="#444444", linestyle="--", linewidth=0.9)
        ax.axvline(LFC_THRESHOLD, color="#444444", linestyle="--", linewidth=0.9)
        ax.axhline(-np.log10(PADJ_THRESHOLD), color="#444444", linestyle="--", linewidth=0.9)

        x_abs = float(np.nanpercentile(np.abs(df["log2FoldChange"]), 99.5))
        x_limit = max(3.0, min(max(x_abs * 1.1, LFC_THRESHOLD + 0.6), 12.0))
        ax.set_xlim(-x_limit, x_limit)
        ax.set_ylim(0, max(3.0, float(df["neg_log10_padj"].quantile(0.995)) * 1.08))
        ax.grid(alpha=0.16, linewidth=0.6)
        ax.set_title(dataset_name, fontsize=13, fontweight="bold")
        ax.text(
            0.02,
            0.98,
            f"Up {counts['up']:,} | Down {counts['down']:,}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#D0D0D0", "alpha": 0.93},
        )

    for ax in axes[:, 0]:
        ax.set_ylabel("-log10 adjusted p-value", fontsize=10)
    for ax in axes[1, :]:
        ax.set_xlabel("log2 Fold Change", fontsize=10)

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORS["Up"], label="Up", markersize=7),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORS["Down"], label="Down", markersize=7),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=COLORS["Not significant"],
            label="Not significant",
            markersize=7,
        ),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(
        "Volcano Plots Across Six Datasets",
        fontsize=18,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.965,
        f"Threshold: |log2FC| ≥ {LFC_THRESHOLD:g} and adjusted p-value < {PADJ_THRESHOLD:g}",
        ha="center",
        va="top",
        fontsize=10,
        color="#444444",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(output_dir / f"volcano_collage_6sets.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    volcano_dir = Path(__file__).resolve().parents[1]
    repo_root = Path(__file__).resolve().parents[2]
    input_dir = repo_root / "venn_new" / "gene_id_mapping" / "mapped_resLFC"
    plot_dir = volcano_dir / "plots"
    table_dir = volcano_dir / "tables"
    plot_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    plot_data = []
    for dataset_name, filename in DATASETS.items():
        df, counts = prepare_volcano_data(input_dir / filename)
        safe_name = dataset_name.replace(" ", "_").replace("-", "_")
        plot_volcano(df, dataset_name, counts, plot_dir / f"volcano_{safe_name}.png")
        df.to_csv(table_dir / f"volcano_data_{safe_name}.csv", index=False)

        summaries.append({"dataset": dataset_name, **counts})
        plot_data.append((dataset_name, df, counts))
        print(f"{dataset_name}: up={counts['up']:,}, down={counts['down']:,}, tested={counts['total_tested']:,}")

    save_collage(plot_data, plot_dir)
    pd.DataFrame(summaries).to_csv(table_dir / "volcano_summary_counts.csv", index=False)
    print(f"Saved volcano outputs to: {volcano_dir}")


if __name__ == "__main__":
    main()
