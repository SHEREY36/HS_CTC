#!/usr/bin/env python3
"""Generate deterministic HS_CTC sweep and node manifests."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ALPHAS = tuple([round(0.50 + 0.05 * index, 3) for index in range(10)] + [0.975, 0.990])
THETAS = tuple(round(0.1 * index, 2) for index in range(1, 21))
ASPECT_RATIOS = (1.0, 1.5, 2.0, 3.0, 4.0)


def node_key(alpha: float, theta: float, aspect_ratio: float) -> str:
    return f"alpha_{alpha:.3f}_theta_{theta:.2f}_AR_{aspect_ratio:.2f}"


def seed_for(theta_index: int, ar_index: int, shard_id: int, base_seed: int = 12345) -> int:
    # Deliberately independent of alpha: this supplies common random numbers
    # along each alpha line while retaining independent shards.
    return base_seed + 1_000_003 * shard_id + 10_009 * theta_index + 1_009 * ar_index


def read_only_nodes(path: Path | None) -> set[tuple[float, float, float]] | None:
    if path is None:
        return None
    selected = set()
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            selected.add((float(row["alpha"]), float(row["theta"]), float(row["aspect_ratio"])))
    return selected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("pilot", "production", "continuation"), required=True)
    parser.add_argument("--samples", type=int)
    parser.add_argument("--shard-id", type=int)
    parser.add_argument("--base-seed", type=int, default=12345)
    parser.add_argument("--output-root", type=Path, default=Path("results"))
    parser.add_argument("--only-nodes", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nodes-output", type=Path)
    args = parser.parse_args()

    defaults = {
        "pilot": (20_000, 0),
        "production": (80_000, 1),
        "continuation": (100_000, 2),
    }
    default_samples, default_shard = defaults[args.stage]
    samples = args.samples if args.samples is not None else default_samples
    shard_id = args.shard_id if args.shard_id is not None else default_shard
    if samples <= 0 or shard_id < 0:
        raise SystemExit("samples must be positive and shard-id non-negative")
    selected = read_only_nodes(args.only_nodes)

    rows = []
    nodes = []
    node_index = 0
    for alpha in ALPHAS:
        for theta_index, theta in enumerate(THETAS):
            for ar_index, aspect_ratio in enumerate(ASPECT_RATIOS):
                triple = (alpha, theta, aspect_ratio)
                if selected is not None and triple not in selected:
                    continue
                key = node_key(*triple)
                node_dir = args.output_root / key
                run_dir = node_dir / f"shard_{shard_id:03d}"
                nodes.append({
                    "node_index": node_index,
                    "alpha": f"{alpha:.3f}",
                    "theta": f"{theta:.2f}",
                    "aspect_ratio": f"{aspect_ratio:.2f}",
                    "node_key": key,
                    "node_dir": str(node_dir),
                    "coefficient_dir": str(Path("coefficients") / key),
                })
                rows.append({
                    "task_index": len(rows),
                    "stage": args.stage,
                    "alpha": f"{alpha:.3f}",
                    "theta": f"{theta:.2f}",
                    "aspect_ratio": f"{aspect_ratio:.2f}",
                    "seed": seed_for(theta_index, ar_index, shard_id, args.base_seed),
                    "nsamples": samples,
                    "shard_id": shard_id,
                    "node_key": key,
                    "output_dir": str(run_dir),
                })
                node_index += 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, lineterminator="\n", fieldnames=rows[0].keys() if rows else (
            "task_index", "stage", "alpha", "theta", "aspect_ratio", "seed",
            "nsamples", "shard_id", "node_key", "output_dir",
        ))
        writer.writeheader()
        writer.writerows(rows)

    nodes_output = args.nodes_output or args.output.with_suffix(".nodes.csv")
    with nodes_output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, lineterminator="\n", fieldnames=(
            "node_index", "alpha", "theta", "aspect_ratio", "node_key", "node_dir", "coefficient_dir"
        ))
        writer.writeheader()
        seen = set()
        for row in nodes:
            if row["node_key"] not in seen:
                writer.writerow(row)
                seen.add(row["node_key"])
    print(f"wrote {len(rows)} tasks to {args.output}")
    print(f"wrote {len(seen)} nodes to {nodes_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
