#!/usr/bin/env python3
"""
Download only the XML search result files for PXD001676 using the existing manifest.
Saves files to PXD001676/ and skips files that already exist and match expected size.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List

import requests

ROOT = Path(__file__).resolve().parents[1]
PXD_DIR = ROOT / "PXD001676"
def find_manifest(pxd_dir: Path) -> Path:
    # Prefer top-level manifest, otherwise search subdirectories
    candidate = pxd_dir / "PXD001676_manifest.json"
    if candidate.exists():
        return candidate
    for sub in pxd_dir.rglob('PXD001676_manifest.json'):
        return sub
    raise FileNotFoundError(f"Manifest not found in {pxd_dir}")

MANIFEST = find_manifest(PXD_DIR)
TIMEOUT = 60
CHUNK = 1024 * 1024


def load_manifest(path: Path) -> List[dict]:
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    with path.open('r', encoding='utf-8') as f:
        return json.load(f)


def is_xml_record(rec: dict) -> bool:
    fname = rec.get('fileName','').lower()
    cat = rec.get('category','').upper()
    return fname.endswith('_xml.xml') or cat == 'SEARCH'


def ftp_to_https(url: str) -> str:
    if url.startswith('ftp://ftp.pride.ebi.ac.uk/'):
        return 'https://ftp.pride.ebi.ac.uk/' + url[len('ftp://ftp.pride.ebi.ac.uk/'):]
    return url


def download(url: str, dest: Path, expected_size: int=0) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and expected_size>0 and dest.stat().st_size == expected_size:
        print(f"[skip] {dest.name} (already present)")
        return
    # Always start clean when the current file size does not match the manifest.
    # Some PRIDE endpoints do not honor Range reliably, so appending can corrupt XMLs.
    if dest.exists():
        dest.unlink()
    headers = {}
    mode = 'wb'
    print(f"[get ] {dest.name}")
    with requests.get(url, stream=True, timeout=TIMEOUT, headers=headers) as r:
        if r.status_code == 416 and expected_size>0 and dest.exists() and dest.stat().st_size==expected_size:
            print(f"[skip] {dest.name} (already complete)")
            return
        r.raise_for_status()
        with dest.open(mode) as f:
            for chunk in r.iter_content(chunk_size=CHUNK):
                if chunk:
                    f.write(chunk)
    if expected_size>0:
        actual = dest.stat().st_size
        if actual != expected_size:
            print(f"Warning: size mismatch for {dest.name}: expected {expected_size}, got {actual}")


def main():
    print('Loading manifest...')
    try:
        records = load_manifest(MANIFEST)
    except Exception as e:
        print('ERROR:', e)
        return 1

    xml_records = [r for r in records if is_xml_record(r)]
    if not xml_records:
        print('No XML/SEARCH records found in manifest.')
        return 1

    print(f'Found {len(xml_records)} XML records in manifest. Downloading to {PXD_DIR}...')
    for rec in xml_records:
        fname = rec.get('fileName')
        size = int(rec.get('fileSizeBytes') or 0)
        url = rec.get('url')
        if not url:
            # try to pick publicFileLocations style
            public = rec.get('publicFileLocations') or []
            if public:
                url = public[0].get('value')
        if not url:
            print(f"Skipping {fname}: no URL")
            continue
        url = ftp_to_https(url)
        dest = PXD_DIR / fname
        try:
            download(url, dest, size)
        except Exception as e:
            print(f"Failed: {fname} -> {e}")
    print('Done')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
