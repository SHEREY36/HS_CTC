#!/usr/bin/env python3
"""Run a closure sweep manifest sequentially for local testing."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--executable", type=Path, default=ROOT / "build" / "SphCyl")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if not args.executable.is_file():
        raise SystemExit(f"missing executable: {args.executable}; build with make -C build clean all")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    failed = []
    with args.manifest.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if args.limit is not None:
        rows = rows[: args.limit]
    for index, row in enumerate(rows, 1):
        output_dir = ROOT / row["output_dir"]
        if (output_dir / "_SUCCESS").is_file():
            print(f"[{index}/{len(rows)}] skip {row['node_key']} shard={row['shard_id']}")
            continue
        output_dir.mkdir(parents=True, exist_ok=True)
        command = [
            str(args.executable), row["alpha"], row["theta"], "1.0", row["aspect_ratio"],
            str(output_dir), row["seed"], row["nsamples"], "closure",
        ]
        print(f"[{index}/{len(rows)}] {row['node_key']} shard={row['shard_id']}", flush=True)
        result = subprocess.run(command, cwd=ROOT, env=env)
        if result.returncode == 0:
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "validate_closure_run.py"), str(output_dir), "--mark-success"],
                cwd=ROOT,
            )
        if result.returncode != 0:
            failed.append(row["node_key"])
    if failed:
        print("Failed nodes:", *failed, sep="\n  ", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
