#!/usr/bin/env python3
"""Obtain the author-licensed inputs for v0.3.0 and check their checksums.

The repository does not contain the training matrix, the serialized module-clock
pipelines, or gene-level coefficients. This script checks the files you already
have, or downloads them, and stops if a checksum differs from Supplementary Table S20.

    python scripts/setup_v030.py --verify-only

Paths already on disk are given with environment variables (see the manifest below).
--verify-only does not download and does not write extracted clocks.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Checksums are the ones recorded in Supplementary Table S20.
MANIFEST = [
    {
        "name": "rodent expression matrix",
        "source": "Zenodo record 18763485, Expression_data_relative_rodents_Scaled.csv",
        "md5": "c283869cf69e35afa200c6a6fc281519",
        "dest": "data/rodent/Expression_data_relative_rodents_Scaled.csv",
        "env": "MAT_EXPRESSION",
    },
    {
        "name": "serialized rodent module clocks",
        "source": "Zenodo record 22166800, Rodent module clocks",
        "md5": "be9234547aacfbb2e83f0f1c5e44b3bd",
        "dest": "data/clocks/Rodent_module_clocks.zip",
        "env": "MAT_CLOCKS",
    },
    {
        "name": "published Table 5B",
        "source": "MOESM7 Table 5B sheet",
        "md5": "990a4d6cc6a9aaa892606f47fa2e1731",
        "dest": "data/supp/MOESM7_Table5B.xlsx",
        "env": "MAT_TABLE5B",
    },
]

# Public cohorts named in S20. They are downloaded from the provider, not from this repo.
ACCESSIONS = [
    ("GEO GSE136831", "primary lung scRNA-seq, Adams"),
    ("GEO GSE135893", "replication lung scRNA-seq, Habermann"),
    ("GEO GSE134692", "bulk lung RNA-seq"),
    ("GEO GSE124872", "Angelidis mouse lung ageing"),
    ("CELLxGENE census 2025-11-08", "HLCA 9f222629-9e39-47d0-b83f-e08d610c7479"),
    ("CELLxGENE census 2025-11-08", "Tabula Muris Senis lung 48b37086-25f7-4ecd-be66-f5bb378e3aea"),
    ("GTEx v8", "lung_counts.pkl, lung_meta.csv, lung_tech.tsv, SubjectPhenotypes.txt"),
]


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def locate(item) -> Path | None:
    env = os.environ.get(item["env"])
    if env:
        return Path(env)
    candidate = ROOT / item["dest"]
    return candidate if candidate.exists() else None


def verify_only() -> int:
    failed = False
    for item in MANIFEST:
        path = locate(item)
        if path is None or not path.exists():
            print(f"MISSING  {item['name']}: set {item['env']} to the file already on disk")
            print(f"         source: {item['source']}")
            failed = True
            continue
        digest = md5_file(path)
        if digest != item["md5"]:
            print(f"MISMATCH {item['name']}: {path}")
            print(f"         expected {item['md5']}")
            print(f"         got      {digest}")
            failed = True
            continue
        print(f"OK       {item['name']}: {path}")
    print("accessions (not checksummed here):")
    for acc, note in ACCESSIONS:
        print(f"         {acc} — {note}")
    if failed:
        print("stopped: a checksummed input is missing or does not match S20")
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify-only", action="store_true",
                    help="checksum files already on disk; do not download or extract")
    args = ap.parse_args()
    if args.verify_only:
        return verify_only()
    print("Downloading is not run from this entry point.")
    print("Place the files listed in the manifest and run: python scripts/setup_v030.py --verify-only")
    print("v0.2.0 retraining of module clocks (scripts/train_module_clocks.py) is withdrawn.")
    return verify_only()


if __name__ == "__main__":
    sys.exit(main())
