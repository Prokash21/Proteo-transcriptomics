from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


MYGENE_QUERY_URL = "https://mygene.info/v3/query"
ENSEMBL_RE = re.compile(r"^ENS[A-Z0-9]*G[0-9]+(?:\.[0-9]+)?$")
ENTREZ_RE = re.compile(r"^[0-9]+$")


def chunked(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def normalize_source_id(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace('"', "").strip()


def query_id_for_mapping(source_id: str, source_id_type: str) -> str:
    if source_id_type == "ensembl":
        return re.sub(r"\.[0-9]+$", "", source_id)
    return source_id


def detect_id_type(ids: pd.Series) -> tuple[str, dict[str, int]]:
    cleaned = ids.map(normalize_source_id)
    ensembl_count = int(cleaned.map(lambda value: bool(ENSEMBL_RE.match(value))).sum())
    entrez_count = int(cleaned.map(lambda value: bool(ENTREZ_RE.match(value))).sum())
    blank_count = int((cleaned == "").sum())
    other_count = int(len(cleaned) - ensembl_count - entrez_count - blank_count)

    if ensembl_count >= entrez_count:
        detected = "ensembl"
    else:
        detected = "entrez"

    return detected, {
        "total_rows": int(len(cleaned)),
        "ensembl_like": ensembl_count,
        "entrez_like": entrez_count,
        "blank": blank_count,
        "other": other_count,
    }


def extract_ensembl_gene(item: dict) -> str | None:
    ensembl = item.get("ensembl")
    if isinstance(ensembl, dict):
        return ensembl.get("gene")
    if isinstance(ensembl, list):
        for row in ensembl:
            if isinstance(row, dict) and row.get("gene"):
                return row.get("gene")
    return None


def fetch_mygene_mapping(
    ids: list[str],
    source_id_type: str,
    chunk_size: int,
    sleep_seconds: float,
    retries: int,
) -> pd.DataFrame:
    scopes = "ensembl.gene" if source_id_type == "ensembl" else "entrezgene"
    records: list[dict] = []

    for subset_number, subset in enumerate(chunked(ids, chunk_size), start=1):
        encoded = urlencode(
            {
                "q": ",".join(subset),
                "scopes": scopes,
                "fields": "symbol,name,ensembl.gene,taxid",
                "species": "human",
            }
        ).encode("utf-8")

        for attempt in range(1, retries + 1):
            try:
                request = Request(
                    MYGENE_QUERY_URL,
                    data=encoded,
                    headers={"Content-Type": "application/x-www-form-urlencoded"},
                    method="POST",
                )
                with urlopen(request, timeout=120) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                break
            except (HTTPError, URLError, TimeoutError) as error:
                if attempt == retries:
                    raise RuntimeError(
                        f"MyGene request failed for chunk {subset_number} after {retries} attempts: {error}"
                    ) from error
                time.sleep(sleep_seconds * attempt)

        if not isinstance(payload, list):
            continue

        for item in payload:
            source_id = normalize_source_id(item.get("query"))
            if not source_id:
                continue

            not_found = bool(item.get("notfound", False))
            records.append(
                {
                    "query_source_id": source_id,
                    "source_id_type": source_id_type,
                    "gene_id": item.get("symbol") if not not_found else None,
                    "entrez_id": str(item.get("_id")) if not not_found and item.get("_id") is not None else None,
                    "ensembl_gene_id": extract_ensembl_gene(item) if not not_found else None,
                    "gene_name": item.get("name") if not not_found else None,
                    "taxid": item.get("taxid") if not not_found else None,
                    "mygene_not_found": not_found,
                }
            )

        if sleep_seconds:
            time.sleep(sleep_seconds)

    columns = [
        "query_source_id",
        "source_id_type",
        "gene_id",
        "entrez_id",
        "ensembl_gene_id",
        "gene_name",
        "taxid",
        "mygene_not_found",
    ]
    if not records:
        return pd.DataFrame(columns=columns)

    mapping = pd.DataFrame.from_records(records)
    mapping = mapping.drop_duplicates(subset=["query_source_id", "source_id_type"], keep="first")
    return mapping[columns]


def valid_query_ids(ids: pd.Series, source_id_type: str) -> list[str]:
    cleaned = ids.map(lambda value: query_id_for_mapping(normalize_source_id(value), source_id_type))
    if source_id_type == "ensembl":
        return sorted({value for value in cleaned if ENSEMBL_RE.match(value)})
    return sorted({value for value in cleaned if ENTREZ_RE.match(value)})


def convert_file(input_file: Path, mapped_dir: Path, report_dir: Path, args: argparse.Namespace) -> dict:
    df = pd.read_csv(input_file)
    original_first_col = df.columns[0]
    df = df.rename(columns={original_first_col: "source_id"})
    df["source_id"] = df["source_id"].map(normalize_source_id)

    source_id_type, id_counts = detect_id_type(df["source_id"])
    df["source_id_type"] = source_id_type
    df["query_source_id"] = df["source_id"].map(lambda value: query_id_for_mapping(value, source_id_type))

    ids_to_query = valid_query_ids(df["source_id"], source_id_type)
    mapping = fetch_mygene_mapping(
        ids_to_query,
        source_id_type=source_id_type,
        chunk_size=args.chunk_size,
        sleep_seconds=args.sleep_seconds,
        retries=args.retries,
    )

    out = df.merge(mapping, on=["query_source_id", "source_id_type"], how="left")
    out.insert(1, "id_status", f"converted_from_{source_id_type}")
    out.loc[out["gene_id"].isna(), "id_status"] = f"{source_id_type}_unmapped"

    ordered_columns = [
        "source_id",
        "source_id_type",
        "query_source_id",
        "gene_id",
        "entrez_id",
        "ensembl_gene_id",
        "gene_name",
        "taxid",
        "id_status",
        "mygene_not_found",
    ]
    remaining_columns = [column for column in out.columns if column not in ordered_columns]
    out = out[ordered_columns + remaining_columns]

    mapped_file = mapped_dir / f"mapped_{input_file.stem}_mygene.csv"
    out.to_csv(mapped_file, index=False)

    report = {
        "input_file": input_file.name,
        "mapped_file": str(mapped_file),
        "detected_source_id_type": source_id_type,
        "id_pattern_counts": id_counts,
        "unique_ids_queried": len(ids_to_query),
        "mapped_rows": int(out["gene_id"].notna().sum()),
        "unmapped_rows": int(out["gene_id"].isna().sum()),
        "mapped_unique_gene_ids": int(out["gene_id"].dropna().nunique()),
    }
    report_file = report_dir / f"{input_file.stem}_id_mapping_report.json"
    report_file.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(
        f"{input_file.name}: {source_id_type}; "
        f"mapped {report['mapped_rows']}/{id_counts['total_rows']} rows; "
        f"saved {mapped_file.name}"
    )
    return report


def write_summary(reports: list[dict], report_dir: Path) -> None:
    rows = []
    for report in reports:
        rows.append(
            {
                "input_file": report["input_file"],
                "detected_source_id_type": report["detected_source_id_type"],
                "total_rows": report["id_pattern_counts"]["total_rows"],
                "ensembl_like": report["id_pattern_counts"]["ensembl_like"],
                "entrez_like": report["id_pattern_counts"]["entrez_like"],
                "other": report["id_pattern_counts"]["other"],
                "unique_ids_queried": report["unique_ids_queried"],
                "mapped_rows": report["mapped_rows"],
                "unmapped_rows": report["unmapped_rows"],
                "mapped_unique_gene_ids": report["mapped_unique_gene_ids"],
            }
        )
    pd.DataFrame(rows).to_csv(report_dir / "id_mapping_summary.csv", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Map resLFC Entrez/Ensembl IDs to human gene symbols with MyGene.")
    workflow_dir = Path(__file__).resolve().parents[1]
    venn_dir = Path(__file__).resolve().parents[2]
    parser.add_argument("--input-dir", type=Path, default=venn_dir / "data")
    parser.add_argument("--mapped-dir", type=Path, default=workflow_dir / "mapped_resLFC")
    parser.add_argument("--report-dir", type=Path, default=workflow_dir / "reports")
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--retries", type=int, default=3)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.mapped_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    input_files = sorted(args.input_dir.glob("*.csv"))
    if not input_files:
        raise FileNotFoundError(f"No CSV files found in {args.input_dir}")

    reports = [convert_file(input_file, args.mapped_dir, args.report_dir, args) for input_file in input_files]
    write_summary(reports, args.report_dir)
    print(f"Summary saved to {args.report_dir / 'id_mapping_summary.csv'}")


if __name__ == "__main__":
    main()
