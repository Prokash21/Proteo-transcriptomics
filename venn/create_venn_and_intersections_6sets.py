from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from venn import venn


BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "converted_gene_id"
OUTPUT_DIR = BASE_DIR / "venn_intersections_6sets"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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
    for dataset, filename in INPUT_FILES.items():
        path = INPUT_DIR / filename
        if not path.exists():
            raise FileNotFoundError(f"Missing input file: {path}")

        df = pd.read_csv(path)
        if "gene_id" not in df.columns:
            raise ValueError(f"Column 'gene_id' not found in {path}")

        genes = set(df["gene_id"].dropna().astype(str).str.strip())
        genes = {g for g in genes if g and g.upper() != "NA"}
        gene_sets[dataset] = genes

    return gene_sets


def exact_region_gene_sets(gene_sets):
    names = list(gene_sets.keys())
    all_rows = []

    for r in range(1, len(names) + 1):
        for subset in combinations(names, r):
            included = set(subset)
            excluded = [n for n in names if n not in included]

            overlap = set.intersection(*(gene_sets[n] for n in included))
            if excluded:
                overlap -= set.union(*(gene_sets[n] for n in excluded))

            key = " & ".join(subset)
            if overlap:
                for gene in sorted(overlap):
                    all_rows.append(
                        {
                            "intersection": key,
                            "included_sets": "|".join(subset),
                            "n_sets": len(subset),
                            "gene_id": gene,
                        }
                    )
            else:
                all_rows.append(
                    {
                        "intersection": key,
                        "included_sets": "|".join(subset),
                        "n_sets": len(subset),
                        "gene_id": "",
                    }
                )

    return pd.DataFrame(all_rows)


def intersection_counts(gene_sets):
    names = list(gene_sets.keys())
    rows = []

    for r in range(1, len(names) + 1):
        for subset in combinations(names, r):
            included = set(subset)
            excluded = [n for n in names if n not in included]
            overlap = set.intersection(*(gene_sets[n] for n in included))
            if excluded:
                overlap -= set.union(*(gene_sets[n] for n in excluded))

            rows.append(
                {
                    "intersection": " & ".join(subset),
                    "included_sets": "|".join(subset),
                    "n_sets": len(subset),
                    "count": len(overlap),
                }
            )

    return pd.DataFrame(rows).sort_values(["n_sets", "intersection"]).reset_index(drop=True)


def save_venn(gene_sets):
    plt.figure(figsize=(15, 12))
    venn(gene_sets, cmap="Pastel2", legend_loc="upper right", fontsize=10)
    plt.title("Gene Overlap Across 6 Datasets (mapped gene_id)", fontsize=14, fontweight="bold", pad=18)
    plt.tight_layout()

    png_path = OUTPUT_DIR / "venn_6sets_gene_id_python.png"
    pdf_path = OUTPUT_DIR / "venn_6sets_gene_id_python.pdf"

    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    return png_path, pdf_path


def main():
    gene_sets = load_gene_sets()

    print("Loaded gene_id sets:")
    for name, genes in gene_sets.items():
        print(f"  {name}: {len(genes)}")

    png_path, pdf_path = save_venn(gene_sets)
    print(f"Saved Python Venn PNG: {png_path}")
    print(f"Saved Python Venn PDF: {pdf_path}")

    intersections_df = exact_region_gene_sets(gene_sets)
    csv_path = OUTPUT_DIR / "all_possible_exact_intersections_gene_lists_python.csv"
    intersections_df.to_csv(csv_path, index=False)
    print(f"Saved intersections CSV (Python): {csv_path}")

    counts_df = intersection_counts(gene_sets)
    counts_path = OUTPUT_DIR / "exact_intersection_counts_python.csv"
    counts_df.to_csv(counts_path, index=False)
    print(f"Saved counts CSV (Python): {counts_path}")


if __name__ == "__main__":
    main()
