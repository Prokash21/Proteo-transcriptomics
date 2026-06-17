from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from venn import venn


BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR.parent.parent / "converted_gene_id"
OUT_DIR = BASE_DIR

INPUT_FILES = {
    "GSE_40": "mapped_GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal_biomart.csv",
    "GSE_45": "mapped_GSE_45_OR_02_resLFC_condition_tumor_vs_normal_biomart.csv",
    "GSE_60": "mapped_GSE_60_resLFC_raw_wol_Treatment_tumor_vs_normal_biomart.csv",
    "GSE_79": "mapped_GSE_79_02_resLFC_Treatment_tumor_vs_normal_biomart.csv",
    "TCGA_COAD": "mapped_tcga_resLFC_Treatment_Tumor_vs_Normal_biomart.csv",
    "TCGA_READ": "mapped_tcga_read_resLFC_Treatment_Tumor_vs_Normal_biomart.csv",
}


def load_gene_sets():
    gene_sets = {}
    for label, filename in INPUT_FILES.items():
        path = INPUT_DIR / filename
        df = pd.read_csv(path)
        genes = set(df["gene_id"].dropna().astype(str).str.strip())
        genes = {g for g in genes if g and g.upper() != "NA"}
        gene_sets[label] = genes
    return gene_sets


def save_python_venn(gene_sets):
    plt.figure(figsize=(15, 12))
    venn(gene_sets, cmap="Set2", legend_loc="upper right", fontsize=10)
    plt.title("Python: venn package (6 sets)", fontsize=14, fontweight="bold")
    plt.tight_layout()
    png = OUT_DIR / "python_venn_package_6sets.png"
    pdf = OUT_DIR / "python_venn_package_6sets.pdf"
    plt.savefig(png, dpi=300, bbox_inches="tight")
    plt.savefig(pdf, bbox_inches="tight")
    plt.close()
    return png, pdf


def save_supervenn(gene_sets):
    from supervenn import supervenn

    labels = list(gene_sets.keys())
    sets = [gene_sets[k] for k in labels]

    fig = plt.figure(figsize=(16, 8))
    supervenn(
        sets,
        labels,
        min_width_for_annotation=1,
        side_plots=True,
        widths_minmax_ratio=0.02,
    )
    plt.title("Python: supervenn (6 sets)", fontsize=14, fontweight="bold")
    plt.tight_layout()
    png = OUT_DIR / "python_supervenn_6sets.png"
    pdf = OUT_DIR / "python_supervenn_6sets.pdf"
    plt.savefig(png, dpi=300, bbox_inches="tight")
    plt.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def save_matplotlib_venn_pairwise(gene_sets):
    from itertools import combinations
    from matplotlib_venn import venn2

    names = list(gene_sets.keys())
    pairs = list(combinations(names, 2))

    fig, axes = plt.subplots(3, 5, figsize=(24, 14))
    axes = axes.flatten()

    for idx, (a, b) in enumerate(pairs):
        ax = axes[idx]
        venn2([gene_sets[a], gene_sets[b]], set_labels=(a, b), ax=ax, alpha=0.65)
        ax.set_title(f"{a} vs {b}", fontsize=10)

    for idx in range(len(pairs), len(axes)):
        axes[idx].axis("off")

    fig.suptitle("Python: matplotlib-venn pairwise view (6 sets)", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    png = OUT_DIR / "python_matplotlib_venn_pairwise_6sets.png"
    pdf = OUT_DIR / "python_matplotlib_venn_pairwise_6sets.pdf"
    plt.savefig(png, dpi=300, bbox_inches="tight")
    plt.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def write_summary(gene_sets, outputs):
    summary_path = OUT_DIR / "python_plot_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Python 6-set visibility comparison\n")
        f.write("=================================\n\n")
        f.write("Input set sizes:\n")
        for k, v in gene_sets.items():
            f.write(f"- {k}: {len(v)}\n")
        f.write("\nGenerated plots:\n")
        for name, files in outputs.items():
            f.write(f"- {name}:\n")
            f.write(f"  PNG: {files[0].name}\n")
            f.write(f"  PDF: {files[1].name}\n")
        f.write("\nNotes:\n")
        f.write("- venn: classical overlap-style map, visually dense for 6 sets.\n")
        f.write("- supervenn: cleaner block matrix, usually better readability for many sets.\n")
        f.write("- matplotlib-venn pairwise: easier per-dataset comparison, but not one 6-set map.\n")


if __name__ == "__main__":
    gene_sets = load_gene_sets()
    outputs = {}
    outputs["venn"] = save_python_venn(gene_sets)
    outputs["supervenn"] = save_supervenn(gene_sets)
    outputs["matplotlib_venn_pairwise"] = save_matplotlib_venn_pairwise(gene_sets)
    write_summary(gene_sets, outputs)

    print("Python multi-library plots generated:")
    for name, files in outputs.items():
        print(f"  {name}: {files[0].name}, {files[1].name}")
    print("  summary: python_plot_summary.txt")
