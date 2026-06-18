from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle


DATASETS = {
    "GSE180440": "mapped_GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal_mygene.csv",
    "GSE251845": "mapped_GSE_45_OR_02_resLFC_condition_tumor_vs_normal_mygene.csv",
    "GSE50760": "mapped_GSE_60_resLFC_raw_wolTreatment_tumor_vs_normal_mygene.csv",
    "GSE142279": "mapped_GSE_79_02_resLFC_Treatment_tumor_vs_normal_mygene.csv",
    "TCGA-COAD": "mapped_tcga_coad_resLFC_Treatment_Tumor_vs_Normal_mygene.csv",
    "TCGA-READ": "mapped_tcga_read_resLFC_SampleType_Tumor_vs_Normal_mygene.csv",
}

COLORS = {
    "GSE180440": "#2F6B9A",
    "GSE251845": "#E98B2A",
    "GSE50760": "#3E9B55",
    "GSE142279": "#D84A4A",
    "TCGA-COAD": "#955C94",
    "TCGA-READ": "#45A8A2",
}

DISPLAY_NAMES = {
    "GSE180440": "GSE180440",
    "GSE251845": "GSE251845",
    "GSE50760": "GSE50760",
    "GSE142279": "GSE142279",
    "TCGA-COAD": "TCGA COAD",
    "TCGA-READ": "TCGA READ",
}


def load_gene_sets(input_dir: Path) -> dict[str, set[str]]:
    gene_sets: dict[str, set[str]] = {}
    for dataset, filename in DATASETS.items():
        path = input_dir / filename
        if not path.exists():
            raise FileNotFoundError(f"Missing mapped input file: {path}")

        df = pd.read_csv(path, usecols=["gene_id"])
        genes = df["gene_id"].dropna().astype(str).str.strip()
        genes = genes[(genes != "") & (genes.str.upper() != "NA")]
        gene_sets[dataset] = set(genes)

    return gene_sets


def exact_intersection_rows(gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    names = list(gene_sets)
    all_genes = sorted(set.union(*gene_sets.values()))
    rows = []

    for gene in all_genes:
        present = [name for name in names if gene in gene_sets[name]]
        rows.append(
            {
                "gene_id": gene,
                "included_sets": "|".join(present),
                "n_sets": len(present),
                **{name: name in present for name in names},
            }
        )

    return pd.DataFrame(rows)


def exact_intersection_counts(membership: pd.DataFrame) -> pd.DataFrame:
    counts = (
        membership.groupby(["included_sets", "n_sets"], as_index=False)
        .agg(count=("gene_id", "size"))
        .sort_values(["n_sets", "count", "included_sets"], ascending=[False, False, True])
        .reset_index(drop=True)
    )
    counts.insert(0, "intersection", counts["included_sets"].str.replace("|", " & ", regex=False))
    return counts


def combination_overlap_counts(gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    names = list(gene_sets)
    rows = []
    for size in range(1, len(names) + 1):
        for subset in combinations(names, size):
            overlap = set.intersection(*(gene_sets[name] for name in subset))
            rows.append(
                {
                    "intersection": " & ".join(subset),
                    "included_sets": "|".join(subset),
                    "n_sets": size,
                    "count": len(overlap),
                }
            )
    return pd.DataFrame(rows).sort_values(["n_sets", "count"], ascending=[False, False])


def draw_six_set_venn_like_plot(
    gene_sets: dict[str, set[str]],
    exact_counts: pd.DataFrame,
    output_dir: Path,
) -> None:
    names = list(gene_sets)
    all_six_key = "|".join(names)
    all_six_count = int(exact_counts.loc[exact_counts["included_sets"] == all_six_key, "count"].sum())
    union_count = len(set.union(*gene_sets.values()))

    fig, ax = plt.subplots(figsize=(15.5, 10.5))
    fig.patch.set_facecolor("#FBFBFA")
    ax.set_facecolor("#FBFBFA")
    ax.set_aspect("equal")
    ax.axis("off")

    circle_specs = {
        "GSE180440": (-1.35, 0.78),
        "GSE251845": (0.00, 1.30),
        "GSE50760": (1.35, 0.78),
        "GSE142279": (-1.35, -0.72),
        "TCGA-COAD": (0.00, -1.24),
        "TCGA-READ": (1.35, -0.72),
    }
    label_specs = {
        "GSE180440": (-3.05, 1.75),
        "GSE251845": (0.00, 2.92),
        "GSE50760": (3.05, 1.75),
        "GSE142279": (-3.10, -2.02),
        "TCGA-COAD": (0.00, -3.28),
        "TCGA-READ": (3.10, -2.02),
    }

    for name in names:
        circle = Circle(
            circle_specs[name],
            radius=1.82,
            facecolor=COLORS[name],
            edgecolor=COLORS[name],
            alpha=0.28,
            linewidth=2.6,
        )
        ax.add_patch(circle)
        highlight = Circle(
            (circle_specs[name][0] - 0.18, circle_specs[name][1] + 0.18),
            radius=1.82,
            facecolor="none",
            edgecolor="white",
            alpha=0.35,
            linewidth=1.2,
        )
        ax.add_patch(highlight)

    for name in names:
        x, y = label_specs[name]
        ax.text(
            x,
            y,
            f"{DISPLAY_NAMES[name]}\n{len(gene_sets[name]):,}",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color=COLORS[name],
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": COLORS[name], "alpha": 0.85},
        )

    ax.text(
        0,
        0,
        f"Shared by all 6\n{all_six_count:,}",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        color="#1F1F1F",
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "white", "edgecolor": "#3D3D3D", "alpha": 0.95},
    )

    top_exact = exact_counts.head(10).copy()
    table_text = "\n".join(
        f"{row.included_sets.replace('|', ' + ')}: {int(row.count):,}" for row in top_exact.itertuples()
    )
    ax.text(
        3.80,
        0.88,
        "Top Exact Intersections\n" + table_text,
        ha="left",
        va="center",
        fontsize=8.8,
        linespacing=1.25,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "white", "edgecolor": "#D0D0D0", "alpha": 0.97},
    )

    set_text = "\n".join(f"{DISPLAY_NAMES[name]}: {len(gene_sets[name]):,}" for name in names)
    ax.text(
        3.80,
        -1.92,
        "Mapped Gene Sets\n" + set_text,
        ha="left",
        va="center",
        fontsize=9.2,
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "#F3F6F8", "edgecolor": "#D0D7DE", "alpha": 0.98},
    )

    ax.text(
        0,
        -3.95,
        f"Union: {union_count:,} mapped gene symbols. Circle areas are schematic; exact counts are saved as CSV.",
        ha="center",
        va="center",
        fontsize=10,
        color="#444444",
    )

    ax.set_xlim(-4.35, 7.15)
    ax.set_ylim(-4.35, 3.65)
    ax.set_title("Six-Dataset Gene Overlap Venn Diagram", fontsize=19, fontweight="bold", pad=22, color="#111111")

    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(output_dir / f"venn_6sets_gene_id.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_set_size_plot(gene_sets: dict[str, set[str]], output_dir: Path) -> None:
    sizes = pd.Series({name: len(genes) for name, genes in gene_sets.items()}).sort_values(ascending=True)
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.barh(sizes.index, sizes.values, color=[COLORS[name] for name in sizes.index])
    ax.set_xlabel("Mapped unique gene symbols")
    ax.set_title("Mapped Gene Set Sizes", fontsize=14, fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    for index, value in enumerate(sizes.values):
        ax.text(value, index, f" {value:,}", va="center", fontsize=9)
    fig.tight_layout()
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(output_dir / f"gene_set_sizes.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> tuple[Path, Path]:
    workflow_dir = Path(__file__).resolve().parents[1]
    venn_dir = Path(__file__).resolve().parents[2]
    input_dir = venn_dir / "gene_id_mapping" / "mapped_resLFC"
    output_dir = workflow_dir / "plots_and_tables"
    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir


def main() -> None:
    input_dir, output_dir = parse_args()
    gene_sets = load_gene_sets(input_dir)

    membership = exact_intersection_rows(gene_sets)
    exact_counts = exact_intersection_counts(membership)
    combination_counts = combination_overlap_counts(gene_sets)

    membership.to_csv(output_dir / "gene_membership_by_dataset.csv", index=False)
    exact_counts.to_csv(output_dir / "exact_intersection_counts.csv", index=False)
    combination_counts.to_csv(output_dir / "combination_overlap_counts.csv", index=False)

    all_six_genes = membership[membership["n_sets"] == len(gene_sets)][["gene_id"]]
    all_six_genes.to_csv(output_dir / "genes_shared_by_all_6.csv", index=False)

    draw_six_set_venn_like_plot(gene_sets, exact_counts, output_dir)
    save_set_size_plot(gene_sets, output_dir)

    print("Loaded mapped gene sets:")
    for name, genes in gene_sets.items():
        print(f"  {name}: {len(genes):,}")
    print(f"Saved outputs to: {output_dir}")


if __name__ == "__main__":
    main()
