#!/usr/bin/env python3
"""Aggregate node estimates, create alpha=1 limits, and plan continuation shards."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from closure_common import COEFFICIENTS
from make_sweep_manifest import ALPHAS, ASPECT_RATIOS, THETAS, node_key, seed_for


ANCHOR_ALPHAS = (0.90, 0.95, 0.975, 0.990)
FINAL_ALPHAS = tuple(round(0.50 + 0.05 * index, 2) for index in range(11))


def weighted_limit(
    x: np.ndarray, y: np.ndarray, se: np.ndarray
) -> tuple[float, float, float, float, float, float]:
    variance = np.maximum(se * se, 1.0e-24)
    weights = 1.0 / variance
    linear = np.column_stack((np.ones_like(x), x))
    normal = linear.T @ (weights[:, None] * linear)
    rhs = linear.T @ (weights * y)
    beta_linear = np.linalg.solve(normal, rhs)
    statistical_se = math.sqrt(float(np.linalg.inv(normal)[0, 0]))

    quadratic = np.column_stack((np.ones_like(x), x, x * x))
    normal_q = quadratic.T @ (weights[:, None] * quadratic)
    rhs_q = quadratic.T @ (weights * y)
    beta_quadratic = np.linalg.solve(normal_q, rhs_q)
    systematic = abs(float(beta_linear[0] - beta_quadratic[0]))
    total_se = math.hypot(statistical_se, systematic)
    estimate = float(beta_linear[0])
    return (
        estimate, total_se, estimate - 1.96 * total_se, estimate + 1.96 * total_se,
        statistical_se, systematic,
    )


def load_estimates(root: Path) -> list[dict]:
    estimates = []
    for path in sorted(root.rglob("estimate.json")):
        if not (path.parent / "_SUCCESS").is_file():
            continue
        value = json.loads(path.read_text())
        value["_path"] = str(path)
        estimates.append(value)
    return estimates


def metric_row(estimate: dict, name: str) -> tuple[float, float, float, float, str]:
    if name == "f_M":
        baseline = estimate["baseline"]
        return (
            baseline["f_M"], baseline["f_M_standard_error"],
            baseline["f_M_ci_low"], baseline["f_M_ci_high"], "baseline",
        )
    if name == "C_M":
        baseline = estimate["baseline"]
        return (
            baseline["C_M"], baseline["C_M_standard_error"],
            baseline["C_M_ci_low"], baseline["C_M_ci_high"], "baseline",
        )
    coefficient = next(item for item in estimate["coefficients"] if item["name"] == name)
    return (
        coefficient["estimate"], coefficient["standard_error"],
        coefficient["ci_low"], coefficient["ci_high"], coefficient["group"],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coefficients-root", type=Path, default=Path("coefficients"))
    parser.add_argument("--results-root", type=Path, default=Path("results"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--continuation-samples", type=int, default=100_000)
    parser.add_argument("--max-events", type=int, default=1_000_000)
    parser.add_argument("--base-seed", type=int, default=12345)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    estimates = load_estimates(args.coefficients_root)
    by_node = {}
    for item in estimates:
        key = (
            round(float(item["alpha"]), 3), round(float(item["theta"]), 2),
            round(float(item["aspect_ratio"]), 2),
        )
        if key in by_node:
            raise SystemExit(
                f"duplicate estimates for node {key}: {by_node[key]['_path']} and {item['_path']}"
            )
        by_node[key] = item
    expected = {(a, t, ar) for a in ALPHAS for t in THETAS for ar in ASPECT_RATIOS}
    missing = sorted(expected - set(by_node))
    if args.require_complete and missing:
        raise SystemExit(f"missing {len(missing)} simulated-node estimates; first missing node is {missing[0]}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metrics = ("f_M", "C_M") + COEFFICIENTS
    long_rows = []
    qa_rows = []
    for key in sorted(by_node):
        estimate = by_node[key]
        for metric in metrics:
            value, se, ci_low, ci_high, group = metric_row(estimate, metric)
            long_rows.append({
                "alpha": key[0], "theta": key[1], "aspect_ratio": key[2],
                "name": metric, "group": group, "estimate": value,
                "standard_error": se, "ci_low": ci_low, "ci_high": ci_high,
                "statistical_standard_error": se, "systematic_uncertainty": 0.0,
                "n_events": estimate["n_events"], "source": "monte_carlo_projection",
            })
        qa_rows.append({
            "alpha": key[0], "theta": key[1], "aspect_ratio": key[2],
            "n_events": estimate["n_events"], "n_shards": estimate["n_shards"],
            "continuation_required": estimate["continuation_required"],
            "baseline_precision_pass": estimate["baseline"]["precision_pass"],
            "coefficient_precision_pass": all(
                row["precision_pass"] for row in estimate["coefficients"]
            ),
            "projection_agreement_pass": estimate.get("projection_agreement_pass", False),
            "tail_diagnostic_pass": all(row["tail_pass"] for row in estimate["coefficients"]),
            "max_top_0p1pct_contribution_fraction": max(
                row["top_0p1pct_absolute_contribution_fraction"]
                for row in estimate["coefficients"]
            ),
            "leave_one_shard_out_max_absolute_change": estimate["leave_one_shard_out_max_absolute_change"],
            "leave_one_shard_out_f_M_max_absolute_change": estimate.get(
                "leave_one_shard_out_f_M_max_absolute_change"
            ),
            "max_elastic_relative_error": max(
                row["max_elastic_relative_error"] for row in estimate["validation"]
            ),
            "estimate_path": estimate["_path"],
        })

    # Extrapolate the physically defined alpha->1 limit at each theta/AR node.
    extrapolated: dict[
        tuple[float, float, str], tuple[float, float, float, float, float, float]
    ] = {}
    for theta in THETAS:
        for aspect_ratio in ASPECT_RATIOS:
            anchors = []
            for alpha in ANCHOR_ALPHAS:
                item = by_node.get((round(alpha, 3), theta, aspect_ratio))
                if item is not None:
                    anchors.append(item)
            if len(anchors) != len(ANCHOR_ALPHAS):
                continue
            x = np.array([1.0 - float(item["alpha"]) for item in anchors])
            for metric in metrics:
                rows = [metric_row(item, metric) for item in anchors]
                y = np.array([row[0] for row in rows])
                se = np.array([row[1] for row in rows])
                limit = weighted_limit(x, y, se)
                extrapolated[(theta, aspect_ratio, metric)] = limit
                long_rows.append({
                    "alpha": 1.0, "theta": theta, "aspect_ratio": aspect_ratio,
                    "name": metric, "group": rows[0][4], "estimate": limit[0],
                    "standard_error": limit[1], "ci_low": limit[2], "ci_high": limit[3],
                    "statistical_standard_error": limit[4],
                    "systematic_uncertainty": limit[5],
                    "n_events": sum(int(item["n_events"]) for item in anchors),
                    "source": "extrapolated_alpha_limit",
                })

    long_path = args.output_dir / "closure_coefficients_long.csv"
    long_fields = (
        "alpha", "theta", "aspect_ratio", "name", "group", "estimate",
        "standard_error", "statistical_standard_error", "systematic_uncertainty",
        "ci_low", "ci_high", "n_events", "source",
    )
    with long_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=long_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(long_rows)

    with (args.output_dir / "qa_summary.csv").open("w", newline="") as stream:
        fields = tuple(qa_rows[0].keys()) if qa_rows else (
            "alpha", "theta", "aspect_ratio", "n_events", "n_shards",
            "continuation_required", "baseline_precision_pass",
            "coefficient_precision_pass", "projection_agreement_pass",
            "tail_diagnostic_pass", "max_top_0p1pct_contribution_fraction",
            "leave_one_shard_out_max_absolute_change",
            "leave_one_shard_out_f_M_max_absolute_change", "max_elastic_relative_error",
            "estimate_path",
        )
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(qa_rows)

    continuation_rows = []
    for qa in qa_rows:
        if not qa["continuation_required"] or int(qa["n_events"]) >= args.max_events:
            continue
        alpha = float(qa["alpha"]); theta = float(qa["theta"]); ar = float(qa["aspect_ratio"])
        node_dir = args.results_root / node_key(alpha, theta, ar)
        existing_ids = []
        for path in node_dir.glob("shard_*"):
            try:
                existing_ids.append(int(path.name.split("_")[-1]))
            except ValueError:
                pass
        shard_id = max(existing_ids, default=1) + 1
        theta_index = THETAS.index(round(theta, 2))
        ar_index = ASPECT_RATIOS.index(round(ar, 2))
        continuation_rows.append({
            "task_index": len(continuation_rows), "stage": "continuation",
            "alpha": f"{alpha:.3f}", "theta": f"{theta:.2f}", "aspect_ratio": f"{ar:.2f}",
            "seed": seed_for(theta_index, ar_index, shard_id, args.base_seed),
            "nsamples": args.continuation_samples, "shard_id": shard_id,
            "node_key": node_key(alpha, theta, ar),
            "output_dir": str(node_dir / f"shard_{shard_id:03d}"),
        })
    continuation_path = args.output_dir / "continuation_manifest.csv"
    continuation_fields = (
        "task_index", "stage", "alpha", "theta", "aspect_ratio", "seed",
        "nsamples", "shard_id", "node_key", "output_dir",
    )
    with continuation_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=continuation_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(continuation_rows)

    continuation_nodes_path = args.output_dir / "continuation_nodes.csv"
    node_fields = (
        "node_index", "alpha", "theta", "aspect_ratio", "node_key", "node_dir", "coefficient_dir",
    )
    with continuation_nodes_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=node_fields, lineterminator="\n")
        writer.writeheader()
        for node_index, row in enumerate(continuation_rows):
            node_dir = args.results_root / row["node_key"]
            writer.writerow({
                "node_index": node_index, "alpha": row["alpha"], "theta": row["theta"],
                "aspect_ratio": row["aspect_ratio"], "node_key": row["node_key"],
                "node_dir": node_dir, "coefficient_dir": args.coefficients_root / row["node_key"],
            })

    # Dense final table: original 0.50..0.95 nodes plus extrapolated alpha=1.
    shape = (len(FINAL_ALPHAS), len(THETAS), len(ASPECT_RATIOS))
    cm = np.full(shape, np.nan)
    fm = np.full(shape, np.nan)
    cm_se = np.full(shape, np.nan)
    fm_se = np.full(shape, np.nan)
    cm_systematic = np.full(shape, np.nan)
    fm_systematic = np.full(shape, np.nan)
    beta = np.full(shape + (len(COEFFICIENTS),), np.nan)
    beta_se = np.full_like(beta, np.nan)
    beta_systematic = np.full_like(beta, np.nan)
    beta_low = np.full_like(beta, np.nan)
    beta_high = np.full_like(beta, np.nan)
    for ia, alpha in enumerate(FINAL_ALPHAS):
        for it, theta in enumerate(THETAS):
            for ir, ar in enumerate(ASPECT_RATIOS):
                if alpha < 1.0:
                    item = by_node.get((round(alpha, 3), theta, ar))
                    if item is None:
                        continue
                    fm[ia, it, ir] = item["baseline"]["f_M"]
                    cm[ia, it, ir] = item["baseline"]["C_M"]
                    fm_se[ia, it, ir] = item["baseline"]["f_M_standard_error"]
                    cm_se[ia, it, ir] = item["baseline"]["C_M_standard_error"]
                    fm_systematic[ia, it, ir] = 0.0
                    cm_systematic[ia, it, ir] = 0.0
                    for ic, name in enumerate(COEFFICIENTS):
                        row = next(value for value in item["coefficients"] if value["name"] == name)
                        beta[ia, it, ir, ic] = row["estimate"]
                        beta_se[ia, it, ir, ic] = row["standard_error"]
                        beta_systematic[ia, it, ir, ic] = 0.0
                        beta_low[ia, it, ir, ic] = row["ci_low"]
                        beta_high[ia, it, ir, ic] = row["ci_high"]
                else:
                    fm_limit = extrapolated.get((theta, ar, "f_M"))
                    cm_limit = extrapolated.get((theta, ar, "C_M"))
                    if fm_limit is None or cm_limit is None:
                        continue
                    fm[ia, it, ir] = fm_limit[0]
                    cm[ia, it, ir] = cm_limit[0]
                    fm_se[ia, it, ir] = fm_limit[1]
                    cm_se[ia, it, ir] = cm_limit[1]
                    fm_systematic[ia, it, ir] = fm_limit[5]
                    cm_systematic[ia, it, ir] = cm_limit[5]
                    for ic, name in enumerate(COEFFICIENTS):
                        limit = extrapolated[(theta, ar, name)]
                        beta[ia, it, ir, ic] = limit[0]
                        beta_se[ia, it, ir, ic] = limit[1]
                        beta_systematic[ia, it, ir, ic] = limit[5]
                        beta_low[ia, it, ir, ic] = limit[2]
                        beta_high[ia, it, ir, ic] = limit[3]
    np.savez_compressed(
        args.output_dir / "closure_grid.npz",
        alpha=np.asarray(FINAL_ALPHAS), theta=np.asarray(THETAS),
        aspect_ratio=np.asarray(ASPECT_RATIOS), coefficient_names=np.asarray(COEFFICIENTS),
        f_M=fm, f_M_standard_error=fm_se, f_M_systematic_uncertainty=fm_systematic,
        C_M=cm, C_M_standard_error=cm_se, C_M_systematic_uncertainty=cm_systematic,
        beta=beta, beta_standard_error=beta_se,
        beta_systematic_uncertainty=beta_systematic,
        beta_ci_low=beta_low, beta_ci_high=beta_high,
    )
    print(f"loaded {len(estimates)} node estimates; missing {len(missing)}")
    print(f"wrote {len(continuation_rows)} continuation tasks to {continuation_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
