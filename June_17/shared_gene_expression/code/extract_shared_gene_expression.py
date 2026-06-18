from __future__ import annotations

from pathlib import Path

import pandas as pd


DATASETS = [
    {
        "batch": 1,
        "dataset": "GSE180440",
        "mapped_lfc": "mapped_GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal_mygene.csv",
        "raw": "GSE_40/raw_data/GSE180440_count_data_prepared.csv",
        "normalized": "GSE_40/Normalized_Count_Data_GSE_40.csv",
    },
    {
        "batch": 2,
        "dataset": "GSE251845",
        "mapped_lfc": "mapped_GSE_45_OR_02_resLFC_condition_tumor_vs_normal_mygene.csv",
        "raw": "GSE_45/raw_data/GSE251845_colon_raw_counts_prepared.csv",
        "normalized": "GSE_45/Normalized_Count_Data_outliers_removed_17.csv",
    },
    {
        "batch": 3,
        "dataset": "GSE50760",
        "mapped_lfc": "mapped_GSE_60_resLFC_raw_wolTreatment_tumor_vs_normal_mygene.csv",
        "raw": "GSE_60/raw_data/GSE50760_count_data_prepared_raw.csv",
        "normalized": "GSE_60/GSE_60_Normalized_Count_Data_wol.csv",
    },
    {
        "batch": 4,
        "dataset": "GSE142279",
        "mapped_lfc": "mapped_GSE_79_02_resLFC_Treatment_tumor_vs_normal_mygene.csv",
        "raw": "GSE_79/raw_data/count_data_raw_GSE142279.csv",
        "normalized": "GSE_79/Normalized_Count_Data_wol_03.csv",
    },
    {
        "batch": 5,
        "dataset": "TCGA COAD",
        "mapped_lfc": "mapped_tcga_coad_resLFC_Treatment_Tumor_vs_Normal_mygene.csv",
        "raw": "TCGA_coad/raw_data/TCGA_COAD_count_data_prepared.csv",
        "normalized": "TCGA_coad/Normalized_Count_Data_tcga.csv",
    },
    {
        "batch": 6,
        "dataset": "TCGA READ",
        "mapped_lfc": "mapped_tcga_read_resLFC_SampleType_Tumor_vs_Normal_mygene.csv",
        "raw": "TCGA_read/TCGA_read_count_data_prepared.csv",
        "normalized": "TCGA_read/Normalized_Count_Data_tcga_read.csv",
    },
]


def clean_id(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace('"', "").strip()


def strip_ensembl_version(value: str) -> str:
    if value.startswith("ENS") and "." in value:
        return value.rsplit(".", 1)[0]
    return value


def load_shared_genes(path: Path) -> list[str]:
    df = pd.read_csv(path)
    gene_col = "gene_id" if "gene_id" in df.columns else df.columns[0]
    genes = df[gene_col].dropna().astype(str).str.strip()
    genes = genes[(genes != "") & (genes.str.upper() != "NA")]
    return sorted(genes.unique())


def load_mapping(path: Path, shared_genes: set[str]) -> pd.DataFrame:
    mapping = pd.read_csv(
        path,
        usecols=["source_id", "query_source_id", "gene_id"],
        dtype=str,
    )
    mapping["source_id"] = mapping["source_id"].map(clean_id)
    mapping["query_source_id"] = mapping["query_source_id"].map(clean_id)
    mapping["gene_id"] = mapping["gene_id"].map(clean_id)
    mapping = mapping[mapping["gene_id"].isin(shared_genes)]

    source_map = mapping[["source_id", "gene_id"]].rename(columns={"source_id": "matrix_gene_id"})
    query_map = mapping[["query_source_id", "gene_id"]].rename(columns={"query_source_id": "matrix_gene_id"})
    combined = pd.concat([source_map, query_map], ignore_index=True)
    combined["matrix_gene_id"] = combined["matrix_gene_id"].map(strip_ensembl_version)
    combined = combined[(combined["matrix_gene_id"] != "") & (combined["gene_id"] != "")]
    combined = combined.drop_duplicates(subset=["matrix_gene_id"], keep="first")
    return combined


def read_and_map_matrix(
    matrix_path: Path,
    mapping: pd.DataFrame,
    shared_genes: list[str],
) -> tuple[pd.DataFrame, dict[str, int]]:
    df = pd.read_csv(matrix_path)
    original_gene_col = df.columns[0]
    sample_cols = list(df.columns[1:])

    df = df.rename(columns={original_gene_col: "matrix_gene_id"})
    df["matrix_gene_id"] = df["matrix_gene_id"].map(clean_id).map(strip_ensembl_version)
    df = df.merge(mapping, on="matrix_gene_id", how="inner")

    numeric = df[["gene_id"] + sample_cols].copy()
    for sample_col in sample_cols:
        numeric[sample_col] = pd.to_numeric(numeric[sample_col], errors="coerce")

    expression = numeric.groupby("gene_id", as_index=True)[sample_cols].sum(min_count=1)
    expression = expression.reindex(shared_genes)

    report = {
        "input_rows": int(len(df)),
        "input_samples": int(len(sample_cols)),
        "genes_found": int(expression.dropna(how="all").shape[0]),
        "genes_missing": int(expression.isna().all(axis=1).sum()),
    }
    return expression, report


def add_dataset_prefix(expression: pd.DataFrame, dataset_name: str) -> tuple[pd.DataFrame, dict[str, str]]:
    rename_map = {sample: f"{dataset_name.replace(' ', '_')}__{sample}" for sample in expression.columns}
    return expression.rename(columns=rename_map), rename_map


def write_matrix_with_batch_row(expression: pd.DataFrame, sample_metadata: pd.DataFrame, output_path: Path) -> None:
    batch_by_sample = dict(zip(sample_metadata["merged_sample"], sample_metadata["batch"]))
    output = expression.copy().reset_index().rename(columns={"gene_id": "Gene"})
    batch_row = {"Gene": "batch"}
    batch_row.update({sample: str(int(batch_by_sample[sample])) for sample in expression.columns})
    output = pd.concat([pd.DataFrame([batch_row]), output], ignore_index=True)
    output.to_csv(output_path, index=False)


def build_matrix(
    june_17_dir: Path,
    mapped_lfc_dir: Path,
    shared_genes: list[str],
    matrix_type: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    shared_gene_set = set(shared_genes)
    expressions = []
    sample_metadata_rows = []
    extraction_reports = []

    for dataset in DATASETS:
        mapping = load_mapping(mapped_lfc_dir / dataset["mapped_lfc"], shared_gene_set)
        expression, report = read_and_map_matrix(june_17_dir / dataset[matrix_type], mapping, shared_genes)
        expression, rename_map = add_dataset_prefix(expression, dataset["dataset"])
        expressions.append(expression)

        for original_sample, merged_sample in rename_map.items():
            sample_metadata_rows.append(
                {
                    "matrix_type": matrix_type,
                    "merged_sample": merged_sample,
                    "original_sample": original_sample,
                    "batch": dataset["batch"],
                    "dataset": dataset["dataset"],
                    "source_file": dataset[matrix_type],
                }
            )

        extraction_reports.append(
            {
                "matrix_type": matrix_type,
                "batch": dataset["batch"],
                "dataset": dataset["dataset"],
                "source_file": dataset[matrix_type],
                **report,
            }
        )

    merged = pd.concat(expressions, axis=1)
    merged.index.name = "gene_id"
    sample_metadata = pd.DataFrame(sample_metadata_rows)
    return merged, sample_metadata, extraction_reports


def main() -> None:
    june_17_dir = Path(__file__).resolve().parents[2]
    output_dir = june_17_dir / "shared_gene_expression"
    mapped_lfc_dir = june_17_dir / "venn_new" / "gene_id_mapping" / "mapped_resLFC"
    shared_genes_path = (
        june_17_dir
        / "venn_new"
        / "venn_diagram_R"
        / "plots_and_tables"
        / "genes_shared_by_all_6_R.csv"
    )

    matrix_dir = output_dir / "matrices"
    metadata_dir = output_dir / "metadata"
    report_dir = output_dir / "reports"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    shared_genes = load_shared_genes(shared_genes_path)
    all_sample_metadata = []
    all_reports = []

    for matrix_type in ["raw", "normalized"]:
        merged, sample_metadata, reports = build_matrix(june_17_dir, mapped_lfc_dir, shared_genes, matrix_type)

        write_matrix_with_batch_row(
            merged,
            sample_metadata,
            matrix_dir / f"merged_{matrix_type}_shared_4663_genes_with_batch_row.csv",
        )
        merged.to_csv(matrix_dir / f"merged_{matrix_type}_shared_4663_genes_expression_only.csv")
        sample_metadata.to_csv(metadata_dir / f"sample_metadata_{matrix_type}.csv", index=False)

        all_sample_metadata.append(sample_metadata)
        all_reports.extend(reports)

        print(
            f"{matrix_type}: genes={merged.shape[0]}, samples={merged.shape[1]}, "
            f"missing genes total={int(merged.isna().all(axis=1).sum())}"
        )

    sample_metadata_all = pd.concat(all_sample_metadata, ignore_index=True)
    sample_counts = (
        sample_metadata_all.groupby(["batch", "dataset", "matrix_type"], as_index=False)
        .agg(n_samples=("merged_sample", "size"))
        .pivot(index=["batch", "dataset"], columns="matrix_type", values="n_samples")
        .reset_index()
        .rename(columns={"raw": "raw_samples", "normalized": "normalized_samples"})
    )
    batch_metadata = pd.DataFrame(
        [
            {
                "batch": dataset["batch"],
                "dataset": dataset["dataset"],
                "raw_file": dataset["raw"],
                "normalized_file": dataset["normalized"],
                "mapped_lfc_file": dataset["mapped_lfc"],
            }
            for dataset in DATASETS
        ]
    ).merge(sample_counts, on=["batch", "dataset"], how="left")
    batch_metadata.to_csv(metadata_dir / "batch_metadata.csv", index=False)
    sample_metadata_all.to_csv(metadata_dir / "sample_metadata_all.csv", index=False)
    pd.DataFrame(all_reports).to_csv(report_dir / "shared_gene_extraction_report.csv", index=False)
    pd.DataFrame({"gene_id": shared_genes}).to_csv(metadata_dir / "shared_genes_used.csv", index=False)

    print(f"Saved shared-gene expression outputs to: {output_dir}")


if __name__ == "__main__":
    main()
