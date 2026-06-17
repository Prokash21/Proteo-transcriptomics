#!/usr/bin/env python3
"""
Download PRIDE/ProteomeXchange dataset PXD001676.

- Fetches official file list from PRIDE API
- Converts FTP links to HTTPS for robust downloading
- Skips files that are already complete
- Saves a manifest with URLs and sizes
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterable, Optional

import requests

PROJECT = "PXD001676"
PRIDE_FILES_API = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{PROJECT}/files"
ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT_DIR / PROJECT
MANIFEST_PATH = OUTPUT_DIR / f"{PROJECT}_manifest.json"
TIMEOUT = 60
CHUNK_SIZE = 1024 * 1024  # 1 MiB


def ftp_to_https(url: str) -> str:
    if url.startswith("ftp://ftp.pride.ebi.ac.uk/"):
        return "https://ftp.pride.ebi.ac.uk/" + url[len("ftp://ftp.pride.ebi.ac.uk/") :]
    return url


def choose_download_url(public_locations: Iterable[dict]) -> Optional[str]:
    for entry in public_locations:
        value = entry.get("value", "")
        if value.startswith("ftp://ftp.pride.ebi.ac.uk/"):
            return ftp_to_https(value)
    return None


def fetch_manifest() -> list[dict]:
    resp = requests.get(PRIDE_FILES_API, timeout=TIMEOUT)
    resp.raise_for_status()
    payload = resp.json()
    if not isinstance(payload, list):
        raise RuntimeError("Unexpected PRIDE API response format (expected list)")
    return payload


def build_file_records(payload: list[dict]) -> list[dict]:
    records: list[dict] = []
    for item in payload:
        file_name = item.get("fileName")
        if not file_name:
            continue

        size = int(item.get("fileSizeBytes") or 0)
        category = ((item.get("fileCategory") or {}).get("value") or "UNKNOWN").upper()
        url = choose_download_url(item.get("publicFileLocations") or [])
        if not url:
            continue

        records.append(
            {
                "fileName": file_name,
                "fileSizeBytes": size,
                "category": category,
                "url": url,
            }
        )

    records.sort(key=lambda x: x["fileName"].lower())
    return records


def save_manifest(records: list[dict]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with MANIFEST_PATH.open("w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)


def download_file(url: str, dest: Path, expected_size: int) -> None:
    if dest.exists() and expected_size > 0 and dest.stat().st_size == expected_size:
        print(f"[skip] {dest.name} (already complete)")
        return

    mode = "wb"
    headers = {}
    resume_from = 0

    if dest.exists():
        resume_from = dest.stat().st_size
        if expected_size > 0 and resume_from < expected_size:
            mode = "ab"
            headers["Range"] = f"bytes={resume_from}-"
        elif expected_size > 0 and resume_from > expected_size:
            dest.unlink()
            resume_from = 0

    print(f"[get ] {dest.name}")
    with requests.get(url, stream=True, timeout=TIMEOUT, headers=headers) as resp:
        if resp.status_code == 416 and expected_size > 0 and dest.exists() and dest.stat().st_size == expected_size:
            print(f"[skip] {dest.name} (already complete)")
            return
        resp.raise_for_status()

        with dest.open(mode) as f:
            for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    f.write(chunk)

    if expected_size > 0:
        actual = dest.stat().st_size
        if actual != expected_size:
            raise RuntimeError(
                f"Size mismatch for {dest.name}: expected {expected_size}, got {actual}"
            )

    print(f"[done] {dest.name}")


def summarize(records: list[dict]) -> None:
    by_category: dict[str, int] = {}
    total_bytes = 0
    for rec in records:
        total_bytes += rec["fileSizeBytes"]
        by_category[rec["category"]] = by_category.get(rec["category"], 0) + 1

    gib = total_bytes / (1024 ** 3)
    print(f"Project: {PROJECT}")
    print(f"Files  : {len(records)}")
    print(f"Size   : {gib:.2f} GiB")
    print("Categories:")
    for cat, count in sorted(by_category.items()):
        print(f"  - {cat}: {count}")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Fetching manifest from {PRIDE_FILES_API} ...")
    payload = fetch_manifest()
    records = build_file_records(payload)

    if not records:
        print("No downloadable files found.")
        return 1

    save_manifest(records)
    summarize(records)
    print(f"Manifest saved to: {MANIFEST_PATH}")

    for rec in records:
        dest = OUTPUT_DIR / rec["fileName"]
        download_file(rec["url"], dest, rec["fileSizeBytes"])

    print("All files downloaded successfully.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted by user.")
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
