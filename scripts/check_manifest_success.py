#!/usr/bin/env python3
"""Check that every row in a sweep or node manifest has a success marker."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--field", choices=("output_dir", "coefficient_dir"),
        help="directory column to check; inferred when omitted",
    )
    args = parser.parse_args()

    with args.manifest.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit(f"manifest has no rows: {args.manifest}")
    field = args.field or ("output_dir" if "output_dir" in rows[0] else "coefficient_dir")
    missing = [row[field] for row in rows if not (Path(row[field]) / "_SUCCESS").is_file()]
    print(f"{len(rows) - len(missing)}/{len(rows)} success markers present ({field})")
    if missing:
        print("First missing paths:")
        for path in missing[:20]:
            print(f"  {path}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
