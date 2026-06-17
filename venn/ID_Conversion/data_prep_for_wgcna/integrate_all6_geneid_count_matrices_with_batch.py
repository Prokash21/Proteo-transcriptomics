import csv
from pathlib import Path


def clean(value: str) -> str:
    if value is None:
        return ""
    v = value.strip()
    if v.startswith('"') and v.endswith('"') and len(v) >= 2:
        v = v[1:-1].strip()
    return v


def load_gene_order(gene_file: Path):
    genes = []
    seen = set()
    with gene_file.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        _ = next(reader, None)
        for row in reader:
            if not row:
                continue
            g = clean(row[0])
            if not g or g.upper() == "NA":
                continue
            if g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


def make_unique_sample_name(base_name: str, used_names: set) -> str:
    if base_name not in used_names:
        used_names.add(base_name)
        return base_name

    i = 2
    while True:
        candidate = f"{base_name}__dup{i}"
        if candidate not in used_names:
            used_names.add(candidate)
            return candidate
        i += 1


def load_dataset_matrix(path: Path, dataset: str, batch: str, used_sample_names: set):
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        raw_samples = [clean(x) for x in header[1:]]

        integrated_samples = []
        sample_metadata = []
        for s in raw_samples:
            prefixed = f"{dataset}__{s}"
            unique_name = make_unique_sample_name(prefixed, used_sample_names)
            integrated_samples.append(unique_name)
            sample_metadata.append(
                {
                    "integrated_sample_id": unique_name,
                    "original_sample_id": s,
                    "dataset": dataset,
                    "batch": batch,
                }
            )

        values_by_gene = {}
        for row in reader:
            if not row:
                continue
            gene_id = clean(row[0])
            if not gene_id:
                continue
            values_by_gene[gene_id] = [clean(x) for x in row[1:]]

    return integrated_samples, sample_metadata, values_by_gene


def main():
    base_dir = Path(__file__).resolve().parent
    input_dir = base_dir / "gene_id_count_matrices"
    out_dir = base_dir / "integrated_gene_id_count_matrices"
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_file = base_dir / "all_6way_intersection_gene_id_4470.csv"
    genes = load_gene_order(gene_file)

    datasets = [
        ("GSE_40", "batch1", input_dir / "GSE_40_count_matrix_gene_id_from_6way_intersection.csv"),
        ("GSE_45", "batch2", input_dir / "GSE_45_count_matrix_gene_id_from_6way_intersection.csv"),
        ("GSE_60", "batch3", input_dir / "GSE_60_count_matrix_gene_id_from_6way_intersection.csv"),
        ("GSE_79", "batch4", input_dir / "GSE_79_count_matrix_gene_id_from_6way_intersection.csv"),
        ("TCGA-COAD", "batch5", input_dir / "TCGA-COAD_count_matrix_gene_id_from_6way_intersection.csv"),
        ("TCGA-READ", "batch6", input_dir / "TCGA-READ_count_matrix_gene_id_from_6way_intersection.csv"),
    ]

    used_sample_names = set()
    merged_samples = []
    merged_metadata = []
    dataset_gene_values = {}
    samples_per_dataset = {}

    for dataset, batch, path in datasets:
        integrated_samples, sample_metadata, values_by_gene = load_dataset_matrix(
            path=path,
            dataset=dataset,
            batch=batch,
            used_sample_names=used_sample_names,
        )

        merged_samples.extend(integrated_samples)
        merged_metadata.extend(sample_metadata)
        dataset_gene_values[dataset] = values_by_gene
        samples_per_dataset[dataset] = len(integrated_samples)

    out_matrix = out_dir / "integrated_all6_count_matrix_4470_gene_id.csv"
    out_meta = out_dir / "integrated_all6_sample_batch_metadata.csv"
    out_summary = out_dir / "integrated_all6_summary.txt"

    missing_by_dataset = {d: 0 for d, _, _ in datasets}

    with out_matrix.open("w", newline="", encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["gene_id", *merged_samples])

        for gene_id in genes:
            row_values = [gene_id]
            for dataset, _, _ in datasets:
                vals = dataset_gene_values[dataset].get(gene_id)
                if vals is None:
                    vals = ["0"] * samples_per_dataset[dataset]
                    missing_by_dataset[dataset] += 1
                row_values.extend(vals)
            writer.writerow(row_values)

    with out_meta.open("w", newline="", encoding="utf-8") as f_meta:
        writer = csv.DictWriter(
            f_meta,
            fieldnames=["integrated_sample_id", "original_sample_id", "dataset", "batch"],
        )
        writer.writeheader()
        writer.writerows(merged_metadata)

    with out_summary.open("w", encoding="utf-8") as f_sum:
        f_sum.write("Integrated all 6 count matrices based on 4470 common gene_id\n")
        f_sum.write(f"Input gene count target: {len(genes)}\n")
        f_sum.write(f"Output matrix rows (excluding header): {len(genes)}\n")
        f_sum.write(f"Output matrix samples (all datasets): {len(merged_samples)}\n")
        f_sum.write("Batch map:\n")
        f_sum.write("  GSE_40 -> batch1\n")
        f_sum.write("  GSE_45 -> batch2\n")
        f_sum.write("  GSE_60 -> batch3\n")
        f_sum.write("  GSE_79 -> batch4\n")
        f_sum.write("  TCGA-COAD -> batch5\n")
        f_sum.write("  TCGA-READ -> batch6\n")
        f_sum.write("Missing genes per dataset (filled as 0 if any):\n")
        for dataset, _, _ in datasets:
            f_sum.write(f"  {dataset}: {missing_by_dataset[dataset]}\n")

    print(f"Wrote integrated matrix: {out_matrix}")
    print(f"Wrote sample metadata: {out_meta}")
    print(f"Wrote summary: {out_summary}")


if __name__ == "__main__":
    main()
