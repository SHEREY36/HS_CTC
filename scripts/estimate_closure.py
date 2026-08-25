#!/usr/bin/env python3
"""Estimate the standalone 16-coefficient closure at one grid node."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from closure_common import estimate_from_runs, load_run, write_json


def discover_shards(node_dir: Path) -> list[Path]:
    shards = sorted(path for path in node_dir.glob("shard_*") if path.is_dir())
    return [path for path in shards if (path / "_SUCCESS").is_file()]


def write_long_csv(path: Path, estimate: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "alpha", "theta", "aspect_ratio", "n_events", "name", "group",
        "estimate", "standard_error", "ci_low", "ci_high", "ci_half_width",
        "precision_threshold", "precision_pass", "tail_pass",
        "top_0p1pct_absolute_contribution_fraction",
    )
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for coefficient in estimate["coefficients"]:
            row = {key: estimate[key] for key in ("alpha", "theta", "aspect_ratio", "n_events")}
            row.update(coefficient)
            writer.writerow({key: row[key] for key in fields})


def main() -> int:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--node-dir", type=Path)
    source.add_argument("--run-dir", type=Path, nargs="+")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=128)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260825)
    args = parser.parse_args()

    shard_dirs = discover_shards(args.node_dir) if args.node_dir else args.run_dir
    if not shard_dirs:
        raise SystemExit("no validated shard directories were found")
    runs = [load_run(path) for path in shard_dirs]
    estimate = estimate_from_runs(
        runs,
        n_blocks=args.blocks,
        n_bootstrap=args.bootstrap,
        bootstrap_seed=args.bootstrap_seed,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "estimate.json", estimate)
    write_long_csv(args.output_dir / "coefficients.csv", estimate)
    (args.output_dir / "_SUCCESS").write_text("estimated\n")
    print(json.dumps({
        "alpha": estimate["alpha"],
        "theta": estimate["theta"],
        "aspect_ratio": estimate["aspect_ratio"],
        "n_events": estimate["n_events"],
        "continuation_required": estimate["continuation_required"],
        "output_dir": str(args.output_dir),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
