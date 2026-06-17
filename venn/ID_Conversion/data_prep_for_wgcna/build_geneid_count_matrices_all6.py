import csv
from collections import defaultdict
from pathlib import Path


def clean(value: str) -> str:
    if value is None:
        return ""
    v = value.strip()
    if v.startswith('"') and v.endswith('"') and len(v) >= 2:
        v = v[1:-1].strip()
    return v


def parse_count(value: str) -> int:
    v = clean(value)
    if v == "":
        return 0
    return int(float(v))


def read_target_genes(gene_file: Path):
    genes = []
    seen = set()
    with gene_file.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        _ = next(reader, None)
        for row in reader:
            if not row:
                continue
            gene_id = clean(row[0])
            if not gene_id or gene_id.upper() == "NA":
                continue
            if gene_id not in seen:
                seen.add(gene_id)
                genes.append(gene_id)
    return genes


def read_mapping(mapping_file: Path):
    gene_to_sources = defaultdict(list)
    with mapping_file.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id = clean(row.get("gene_id", ""))
            source_id = clean(row.get("source_id", ""))
            if not gene_id or not source_id:
                continue
            if source_id not in gene_to_sources[gene_id]:
                gene_to_sources[gene_id].append(source_id)
    return gene_to_sources


def read_counts(count_file: Path):
    with count_file.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        sample_ids = [clean(x) for x in header[1:]]

        counts_by_source = {}
        for row in reader:
            if not row:
                continue
            source_id = clean(row[0])
            if not source_id:
                continue
            vals = [parse_count(x) for x in row[1:]]

            if source_id in counts_by_source:
                counts_by_source[source_id] = [a + b for a, b in zip(counts_by_source[source_id], vals)]
            else:
                counts_by_source[source_id] = vals

    return sample_ids, counts_by_source


def build_dataset_matrix(dataset_name, genes, gene_to_sources, sample_ids, counts_by_source, out_dir: Path):
    out_matrix = out_dir / f"{dataset_name}_count_matrix_gene_id_from_6way_intersection.csv"
    out_missing = out_dir / f"{dataset_name}_gene_id_missing_or_unmapped.csv"
    out_summary = out_dir / f"{dataset_name}_gene_id_count_matrix_summary.txt"

    written_rows = 0
    missing = []
    multi_source_genes = 0

    with out_matrix.open("w", newline="", encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["gene_id", *sample_ids])

        for gene_id in genes:
            source_ids = gene_to_sources.get(gene_id, [])
            if len(source_ids) > 1:
                multi_source_genes += 1

            agg = [0] * len(sample_ids)
            used = 0

            for source_id in source_ids:
                vals = counts_by_source.get(source_id)
                if vals is None:
                    continue
                agg = [a + b for a, b in zip(agg, vals)]
                used += 1

            if used == 0:
                reason = "no_mapping" if not source_ids else "no_source_in_counts"
                missing.append((gene_id, "|".join(source_ids), reason))
                continue

            writer.writerow([gene_id, *agg])
            written_rows += 1

    with out_missing.open("w", newline="", encoding="utf-8") as f_miss:
        writer = csv.writer(f_miss)
        writer.writerow(["gene_id", "mapped_source_ids", "reason"])
        writer.writerows(missing)

    with out_summary.open("w", encoding="utf-8") as f_sum:
        f_sum.write(f"Dataset: {dataset_name}\n")
        f_sum.write(f"Input genes requested: {len(genes)}\n")
        f_sum.write(f"Genes written to output matrix: {written_rows}\n")
        f_sum.write(f"Genes missing mapping or count rows: {len(missing)}\n")
        f_sum.write(f"Genes mapped to multiple source_id values (aggregated by sum): {multi_source_genes}\n")
        f_sum.write(f"Output matrix: {out_matrix.name}\n")
        f_sum.write(f"Missing report: {out_missing.name}\n")

    return {
        "dataset": dataset_name,
        "input_genes": len(genes),
        "written": written_rows,
        "missing": len(missing),
        "multi_source": multi_source_genes,
        "matrix_file": out_matrix.name,
        "missing_file": out_missing.name,
        "summary_file": out_summary.name,
    }


def main():
    base_dir = Path(__file__).resolve().parent

    gene_file = base_dir / "all_6way_intersection_gene_id_4470.csv"
    out_dir = base_dir / "gene_id_count_matrices"
    out_dir.mkdir(parents=True, exist_ok=True)

    datasets = [
        {
            "name": "GSE_45",
            "map": base_dir / "mapped_lfc" / "mapped_GSE_45_OR_02_resLFC_condition_tumor_vs_normal_biomart.csv",
            "counts": base_dir / "counts" / "GSE_45_count_data_outliers_removed.csv",
        },
        {
            "name": "GSE_60",
            "map": base_dir / "mapped_lfc" / "mapped_GSE_60_resLFC_raw_wol_Treatment_tumor_vs_normal_biomart.csv",
            "counts": base_dir / "counts" / "GSE_60_count_data_raw_without_outlier.csv",
        },
        {
            "name": "GSE_79",
            "map": base_dir / "mapped_lfc" / "mapped_GSE_79_02_resLFC_Treatment_tumor_vs_normal_biomart.csv",
            "counts": base_dir / "counts" / "GSE_79_count_data.csv",
        },
        {
            "name": "TCGA-COAD",
            "map": base_dir / "mapped_lfc" / "mapped_TCGA_COAD_resLFC_Treatment_Tumor_vs_Normal_biomart.csv",
            "counts": base_dir / "counts" / "TCGA-COAD_count_data.csv",
        },
        {
            "name": "TCGA-READ",
            "map": base_dir / "mapped_lfc" / "mapped_TCGA_READ_resLFC_Treatment_Tumor_vs_Normal_biomart.csv",
            "counts": base_dir / "counts" / "TCGA-READ_count_data.csv",
        },
    ]

    genes = read_target_genes(gene_file)

    summary_rows = []
    for ds in datasets:
        gene_to_sources = read_mapping(ds["map"])
        sample_ids, counts_by_source = read_counts(ds["counts"])

        result = build_dataset_matrix(
            dataset_name=ds["name"],
            genes=genes,
            gene_to_sources=gene_to_sources,
            sample_ids=sample_ids,
            counts_by_source=counts_by_source,
            out_dir=out_dir,
        )
        summary_rows.append(result)
        print(f"Done: {ds['name']} -> {result['written']} genes written, missing={result['missing']}")

    combined_summary = out_dir / "all6_gene_id_count_matrices_summary.csv"
    with combined_summary.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "dataset",
                "input_genes",
                "written",
                "missing",
                "multi_source",
                "matrix_file",
                "missing_file",
                "summary_file",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Combined summary: {combined_summary}")


if __name__ == "__main__":
    main()
