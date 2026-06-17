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


base_dir = Path(__file__).resolve().parent
file_genes = base_dir / "all_6way_intersection_gene_id_4470.csv"
file_map = base_dir / "mapped_lfc" / "mapped_GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal_biomart.csv"
file_counts = base_dir / "counts" / "GSE_40_count_data_removed.csv"

out_dir = base_dir / "gene_id_count_matrices"
out_dir.mkdir(parents=True, exist_ok=True)

out_matrix = out_dir / "GSE_40_count_matrix_gene_id_from_6way_intersection.csv"
out_missing = out_dir / "GSE_40_gene_id_missing_or_unmapped.csv"
out_summary = out_dir / "GSE_40_gene_id_count_matrix_summary.txt"

# 1) Read target genes, preserving order and removing empty/NA values.
genes = []
seen_gene = set()
with file_genes.open("r", newline="", encoding="utf-8-sig") as f:
    reader = csv.reader(f)
    _ = next(reader, None)  # header
    for row in reader:
        if not row:
            continue
        gene_id = clean(row[0])
        if not gene_id:
            continue
        if gene_id.upper() == "NA":
            continue
        if gene_id not in seen_gene:
            seen_gene.add(gene_id)
            genes.append(gene_id)

# 2) Read mapping: gene_id -> source_id (entrez ids as strings).
gene_to_sources = defaultdict(list)
with file_map.open("r", newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene_id = clean(row.get("gene_id", ""))
        source_id = clean(row.get("source_id", ""))
        if not gene_id or not source_id:
            continue
        if source_id not in gene_to_sources[gene_id]:
            gene_to_sources[gene_id].append(source_id)

# 3) Read count matrix keyed by source_id (first column).
with file_counts.open("r", newline="", encoding="utf-8-sig") as f:
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
        values = []
        for x in row[1:]:
            x = clean(x)
            if x == "":
                values.append(0)
            else:
                values.append(int(float(x)))
        counts_by_source[source_id] = values

# 4) Build gene_id-based matrix.
written_rows = 0
missing_map_or_counts = []
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
        missing_sources = []

        for source_id in source_ids:
            vals = counts_by_source.get(source_id)
            if vals is None:
                missing_sources.append(source_id)
                continue
            agg = [a + b for a, b in zip(agg, vals)]
            used += 1

        if used == 0:
            missing_map_or_counts.append((gene_id, "|".join(source_ids) if source_ids else "", "no_source_in_counts" if source_ids else "no_mapping"))
            continue

        writer.writerow([gene_id, *agg])
        written_rows += 1

# 5) Write missing/unmapped report and summary.
with out_missing.open("w", newline="", encoding="utf-8") as f_miss:
    writer = csv.writer(f_miss)
    writer.writerow(["gene_id", "mapped_source_ids", "reason"])
    for row in missing_map_or_counts:
        writer.writerow(row)

with out_summary.open("w", encoding="utf-8") as f_sum:
    f_sum.write(f"Input genes requested: {len(genes)}\n")
    f_sum.write(f"Genes written to output matrix: {written_rows}\n")
    f_sum.write(f"Genes missing mapping or count rows: {len(missing_map_or_counts)}\n")
    f_sum.write(f"Genes mapped to multiple source_id values (aggregated by sum): {multi_source_genes}\n")
    f_sum.write(f"Output matrix: {out_matrix.name}\n")
    f_sum.write(f"Missing report: {out_missing.name}\n")

print(f"Done. Wrote {written_rows} genes to: {out_matrix}")
print(f"Summary: {out_summary}")
