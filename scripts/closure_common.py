#!/usr/bin/env python3
"""Shared schema, validation, and projection math for the closure pipeline."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


SCHEMA_VERSION = 1
COLUMNS = (
    "event_id",
    "c1_x", "c1_y", "c1_z",
    "c2_x", "c2_y", "c2_z",
    "omega1_x", "omega1_y", "omega1_z",
    "omega2_x", "omega2_y", "omega2_z",
    "u1_x", "u1_y", "u1_z",
    "u2_x", "u2_y", "u2_z",
    "delta_tr", "delta_rot", "delta_total",
    "et_elastic", "er_elastic", "et_inelastic", "er_inelastic",
    "e_initial", "elastic_rel_error", "n_contact", "impact_parameter",
    "contact_n_x", "contact_n_y", "contact_n_z",
    "centerline_x", "centerline_y", "centerline_z",
    "contact_lambda", "contact_mu",
)
INDEX = {name: idx for idx, name in enumerate(COLUMNS)}

COEFFICIENTS = (
    "beta_a2_tr",
    "beta_a2_rot",
    "beta_a11",
    "beta_PiPi",
    "beta_RR",
    "beta_QQ",
    "beta_PiR",
    "beta_PiQ",
    "beta_RQ",
    "beta_qtrqtr",
    "beta_qrotqrot",
    "beta_qtrqrot",
    "beta_a3_tr",
    "beta_a3_rot",
    "beta_a21",
    "beta_a12",
)
CORE_COEFFICIENTS = frozenset(COEFFICIENTS[:12])


def parse_metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{line_number}: expected key=value")
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def metadata_number(metadata: dict[str, str], key: str, cast=float):
    if key not in metadata:
        raise ValueError(f"metadata is missing {key!r}")
    return cast(metadata[key])


@dataclass(frozen=True)
class RunData:
    directory: Path
    metadata: dict[str, str]
    events: np.ndarray


def load_run(directory: Path | str, copy: bool = False) -> RunData:
    directory = Path(directory)
    metadata_path = directory / "metadata.txt"
    event_path = directory / "closure_events.bin"
    if not metadata_path.is_file():
        raise FileNotFoundError(metadata_path)
    if not event_path.is_file():
        raise FileNotFoundError(event_path)
    metadata = parse_metadata(metadata_path)
    version = metadata_number(metadata, "schema_version", int)
    n_columns = metadata_number(metadata, "n_columns", int)
    if version != SCHEMA_VERSION:
        raise ValueError(f"unsupported schema_version={version}; expected {SCHEMA_VERSION}")
    if n_columns != len(COLUMNS):
        raise ValueError(f"n_columns={n_columns}; expected {len(COLUMNS)}")
    listed_columns = tuple(metadata.get("columns", "").split(","))
    if listed_columns != COLUMNS:
        raise ValueError("metadata column order does not match schema version 1")
    item_count = event_path.stat().st_size // np.dtype("<f8").itemsize
    if event_path.stat().st_size % np.dtype("<f8").itemsize or item_count % n_columns:
        raise ValueError(f"{event_path} is truncated or has a non-integral record count")
    events = np.memmap(event_path, dtype="<f8", mode="r", shape=(item_count // n_columns, n_columns))
    if copy:
        events = np.asarray(events).copy()
    return RunData(directory=directory, metadata=metadata, events=events)


def validate_run(run: RunData, elastic_tolerance: float = 5.0e-3) -> dict:
    data = np.asarray(run.events)
    expected = metadata_number(run.metadata, "nsamples", int)
    errors: list[str] = []
    warnings: list[str] = []
    if len(data) != expected:
        errors.append(f"record count {len(data)} != nsamples {expected}")
    if data.size and not np.isfinite(data).all():
        errors.append("event data contain NaN or Inf")
    if not data.size:
        errors.append("event file is empty")
        return {"status": "fail", "errors": errors, "warnings": warnings, "n_events": 0}

    event_ids = data[:, INDEX["event_id"]]
    if not np.all(event_ids == np.floor(event_ids)):
        errors.append("event IDs are not integral")
    if len(np.unique(event_ids)) != len(event_ids):
        errors.append("event IDs are not unique within the shard")
    event_id_offset = metadata_number(run.metadata, "event_id_offset", int) \
        if "event_id_offset" in run.metadata else 0
    if len(data) == expected and not np.array_equal(
        np.sort(event_ids.astype(np.int64)),
        event_id_offset + np.arange(1, expected + 1, dtype=np.int64),
    ):
        errors.append("event IDs do not match event_id_offset + 1..nsamples")

    u1 = data[:, 13:16]
    u2 = data[:, 16:19]
    n1 = np.linalg.norm(u1, axis=1)
    n2 = np.linalg.norm(u2, axis=1)
    unit_error = float(max(np.max(np.abs(n1 - 1.0)), np.max(np.abs(n2 - 1.0))))
    if unit_error > 2.0e-6:
        errors.append(f"particle-axis unit-norm error {unit_error:.3e} exceeds 2e-6")

    omega1 = data[:, 7:10]
    omega2 = data[:, 10:13]
    perpendicular_error = float(max(
        np.max(np.abs(np.einsum("ij,ij->i", omega1, u1)) / np.maximum(1.0, np.linalg.norm(omega1, axis=1))),
        np.max(np.abs(np.einsum("ij,ij->i", omega2, u2)) / np.maximum(1.0, np.linalg.norm(omega2, axis=1))),
    ))
    if perpendicular_error > 2.0e-6:
        errors.append(f"omega-u perpendicularity error {perpendicular_error:.3e} exceeds 2e-6")

    delta_tr = data[:, INDEX["delta_tr"]]
    delta_rot = data[:, INDEX["delta_rot"]]
    delta = data[:, INDEX["delta_total"]]
    et_elastic = data[:, INDEX["et_elastic"]]
    er_elastic = data[:, INDEX["er_elastic"]]
    et_inelastic = data[:, INDEX["et_inelastic"]]
    er_inelastic = data[:, INDEX["er_inelastic"]]
    e_initial = data[:, INDEX["e_initial"]]
    energy_scale = np.maximum(1.0, np.abs(e_initial))
    identity_error = np.maximum.reduce((
        np.abs(delta - delta_tr - delta_rot),
        np.abs(delta_tr - (et_elastic - et_inelastic)),
        np.abs(delta_rot - (er_elastic - er_inelastic)),
    )) / energy_scale
    max_identity_error = float(np.max(identity_error))
    if max_identity_error > 5.0e-12:
        errors.append(f"dissipation identity error {max_identity_error:.3e} exceeds 5e-12")
    positive_floor = 1.0e-12 * np.maximum(1.0, np.abs(e_initial))
    bad_dissipation = int(np.count_nonzero(delta <= positive_floor))
    alpha = metadata_number(run.metadata, "alpha")
    if alpha < 1.0 and bad_dissipation:
        errors.append(f"{bad_dissipation} events have non-positive/roundoff-level dissipation")

    stored_elastic_error = data[:, INDEX["elastic_rel_error"]]
    calculated_elastic_error = (et_elastic + er_elastic - e_initial) / np.maximum(
        np.abs(e_initial), 1.0e-30
    )
    recorded_error_mismatch = float(np.max(np.abs(
        stored_elastic_error - calculated_elastic_error
    )))
    if recorded_error_mismatch > 5.0e-12:
        errors.append(
            f"stored elastic-error identity mismatch {recorded_error_mismatch:.3e} exceeds 5e-12"
        )
    if np.any(np.column_stack((et_elastic, er_elastic, et_inelastic, er_inelastic)) < -1.0e-12):
        errors.append("one or more stored energy components are negative")
    elastic_error = np.abs(stored_elastic_error)
    max_elastic_error = float(np.max(elastic_error))
    p999_elastic_error = float(np.quantile(elastic_error, 0.999))
    if max_elastic_error > elastic_tolerance:
        errors.append(
            f"elastic replay relative error {max_elastic_error:.3e} exceeds {elastic_tolerance:.3e}"
        )

    contact_normal = data[:, 30:33]
    contact_norm_error = float(np.max(np.abs(np.linalg.norm(contact_normal, axis=1) - 1.0)))
    if contact_norm_error > 2.0e-6:
        errors.append(f"contact-normal unit error {contact_norm_error:.3e} exceeds 2e-6")
    centerline = data[:, 33:36]
    centerline_norm_error = float(np.max(np.abs(np.linalg.norm(centerline, axis=1) - 1.0)))
    if centerline_norm_error > 2.0e-6:
        errors.append(f"centerline unit error {centerline_norm_error:.3e} exceeds 2e-6")
    contact_count = data[:, INDEX["n_contact"]]
    if np.any(contact_count != np.floor(contact_count)) or np.any(contact_count < 1.0):
        errors.append("contact counts must be positive integers")

    return {
        "status": "pass" if not errors else "fail",
        "errors": errors,
        "warnings": warnings,
        "n_events": int(len(data)),
        "max_axis_norm_error": unit_error,
        "max_omega_axis_dot": perpendicular_error,
        "max_dissipation_identity_error": max_identity_error,
        "nonpositive_dissipation_events": bad_dissipation,
        "max_elastic_relative_error": max_elastic_error,
        "p99p9_elastic_relative_error": p999_elastic_error,
        "max_recorded_elastic_error_mismatch": recorded_error_mismatch,
        "max_contact_normal_error": contact_norm_error,
        "max_centerline_norm_error": centerline_norm_error,
    }


def _tensor_basis() -> np.ndarray:
    basis = np.zeros((5, 3, 3), dtype=float)
    basis[0] = np.diag([1.0, -1.0, 0.0]) / math.sqrt(2.0)
    basis[1] = np.diag([1.0, 1.0, -2.0]) / math.sqrt(6.0)
    basis[2, 0, 1] = basis[2, 1, 0] = 1.0 / math.sqrt(2.0)
    basis[3, 0, 2] = basis[3, 2, 0] = 1.0 / math.sqrt(2.0)
    basis[4, 1, 2] = basis[4, 2, 1] = 1.0 / math.sqrt(2.0)
    return basis


TENSOR_BASIS = _tensor_basis()


def _quadratic_scores(vectors: np.ndarray) -> np.ndarray:
    return np.einsum("ni,kij,nj->nk", vectors, TENSOR_BASIS, vectors, optimize=True)


def projection_kernel(events: np.ndarray, metadata: dict[str, str]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Return the 16 event kernels whose residual-weighted means are beta."""
    velocity_scale = metadata_number(metadata, "velocity_scale")
    omega_scale = metadata_number(metadata, "omega_scale")
    if velocity_scale <= 0.0 or omega_scale <= 0.0:
        raise ValueError("velocity_scale and omega_scale must be positive")

    c1 = events[:, 1:4] / velocity_scale
    c2 = events[:, 4:7] / velocity_scale
    omega1 = events[:, 7:10] / omega_scale
    omega2 = events[:, 10:13] / omega_scale
    u1 = events[:, 13:16]
    u2 = events[:, 16:19]
    x1 = np.einsum("ij,ij->i", c1, c1)
    x2 = np.einsum("ij,ij->i", c2, c2)
    y1 = np.einsum("ij,ij->i", omega1, omega1)
    y2 = np.einsum("ij,ij->i", omega2, omega2)

    def radial(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, ...]:
        s1 = 1.5 - x
        s2 = 15.0 / 8.0 - 2.5 * x + 0.5 * x * x
        s3 = 35.0 / 16.0 - 35.0 / 8.0 * x + 7.0 / 4.0 * x * x - x * x * x / 6.0
        l1 = 1.0 - y
        l2 = 1.0 - 2.0 * y + 0.5 * y * y
        l3 = 1.0 - 3.0 * y + 1.5 * y * y - y * y * y / 6.0
        return s1, s2, s3, l1, l2, l3

    s11, s21, s31, l11, l21, l31 = radial(x1, y1)
    s12, s22, s32, l12, l22, l32 = radial(x2, y2)
    scalar = np.column_stack((
        s21 + s22,
        l21 + l22,
        s11 * l11 + s12 * l12,
        -s31 - s32,
        l31 + l32,
        s21 * l11 + s22 * l12,
        s11 * l21 + s12 * l22,
    ))

    pi1 = _quadratic_scores(c1)
    pi2 = _quadratic_scores(c2)
    wow1 = _quadratic_scores(omega1)
    wow2 = _quadratic_scores(omega2)
    uau1 = _quadratic_scores(u1)
    uau2 = _quadratic_scores(u2)
    # Dual scores: their one-particle response is dR/dlambda=A,dQ=0
    # and dR=0,dQ/dlambda=A, respectively.
    r1 = (10.0 / 7.0) * wow1 + (5.0 / 7.0) * uau1
    r2 = (10.0 / 7.0) * wow2 + (5.0 / 7.0) * uau2
    q1 = (10.0 / 7.0) * wow1 + (40.0 / 7.0) * uau1
    q2 = (10.0 / 7.0) * wow2 + (40.0 / 7.0) * uau2

    tensor_components = {
        "beta_PiPi": 8.0 * pi1 * pi2,
        "beta_RR": 8.0 * r1 * r2,
        "beta_QQ": 8.0 * q1 * q2,
        "beta_PiR": 4.0 * (pi1 * r2 + r1 * pi2),
        "beta_PiQ": 4.0 * (pi1 * q2 + q1 * pi2),
        "beta_RQ": 4.0 * (r1 * q2 + q1 * r2),
    }

    qtr1 = c1 * (x1[:, None] - 2.5)
    qtr2 = c2 * (x2[:, None] - 2.5)
    qrot1 = c1 * (y1[:, None] - 1.0)
    qrot2 = c2 * (y2[:, None] - 1.0)
    vector_components = {
        "beta_qtrqtr": qtr1 * qtr2,
        "beta_qrotqrot": qrot1 * qrot2,
        "beta_qtrqrot": qtr1 * qrot2 + qrot1 * qtr2,
    }

    kernel = np.empty((len(events), len(COEFFICIENTS)), dtype=float)
    kernel[:, 0:3] = scalar[:, 0:3]
    for offset, name in enumerate(COEFFICIENTS[3:9], 3):
        kernel[:, offset] = tensor_components[name].mean(axis=1)
    for offset, name in enumerate(COEFFICIENTS[9:12], 9):
        kernel[:, offset] = vector_components[name].mean(axis=1)
    kernel[:, 12:16] = scalar[:, 3:7]

    directional = {**tensor_components, **vector_components}
    return kernel, directional


def _splitmix64(values: np.ndarray) -> np.ndarray:
    mask = np.uint64(0xFFFFFFFFFFFFFFFF)
    z = (values + np.uint64(0x9E3779B97F4A7C15)) & mask
    z = ((z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)) & mask
    z = ((z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)) & mask
    return z ^ (z >> np.uint64(31))


def deterministic_blocks(event_ids: np.ndarray, run_seed: int, n_blocks: int) -> np.ndarray:
    ids = event_ids.astype(np.uint64, copy=False)
    seed_mix = np.uint64((int(run_seed) * 0xD6E8FEB86659FD93) & 0xFFFFFFFFFFFFFFFF)
    mixed = ids ^ seed_mix
    return (_splitmix64(mixed) % np.uint64(n_blocks)).astype(np.int64)


def estimate_from_runs(
    runs: Sequence[RunData],
    n_blocks: int = 128,
    n_bootstrap: int = 2000,
    bootstrap_seed: int = 20260825,
) -> dict:
    if not runs:
        raise ValueError("at least one run is required")
    if n_blocks <= 0 or n_bootstrap <= 1:
        raise ValueError("n_blocks must be positive and n_bootstrap must exceed one")
    reference = runs[0].metadata
    for key in ("alpha", "theta", "aspect_ratio", "velocity_scale", "omega_scale"):
        expected = metadata_number(reference, key)
        for run in runs[1:]:
            actual = metadata_number(run.metadata, key)
            if not math.isclose(actual, expected, rel_tol=2e-12, abs_tol=2e-12):
                raise ValueError(f"shard metadata mismatch for {key}: {actual} != {expected}")
    alpha = metadata_number(reference, "alpha")
    if alpha >= 1.0:
        raise ValueError("alpha=1 has zero dissipation; estimate its limit by aggregation")

    validations = [validate_run(run) for run in runs]
    failed = [v for v in validations if v["status"] != "pass"]
    if failed:
        raise ValueError("one or more shards failed validation: " + json.dumps(failed, sort_keys=True))

    arrays = [np.asarray(run.events) for run in runs]
    events = np.concatenate(arrays, axis=0) if len(arrays) > 1 else arrays[0]
    if len(np.unique(events[:, INDEX["event_id"]])) != len(events):
        raise ValueError("event IDs are not unique across the selected shards")
    # Tiny smoke datasets cannot populate 128 hashed blocks reliably. Production
    # nodes retain the requested 128; local tests use at least ten events/block.
    requested_n_blocks = n_blocks
    n_blocks = min(n_blocks, max(1, len(events) // 10))
    block_population_pass = n_blocks == requested_n_blocks
    kernel, directional = projection_kernel(events, reference)
    delta_tr = events[:, INDEX["delta_tr"]]
    delta = events[:, INDEX["delta_total"]]
    numerator = float(delta_tr.sum(dtype=np.float64))
    denominator = float(delta.sum(dtype=np.float64))
    f_m = numerator / denominator
    theta = metadata_number(reference, "theta")
    baseline_factor = 3.0 * theta / (3.0 * theta + 2.0)
    c_m = f_m / baseline_factor
    residual = delta_tr - f_m * delta
    coefficients = np.einsum("n,na->a", residual, kernel, optimize=True) / numerator

    # Build deterministic block sufficient statistics. Shard seeds are part of
    # the hash, so changing shard order does not change block membership.
    block_ids_parts = []
    for run, array in zip(runs, arrays):
        block_ids_parts.append(deterministic_blocks(
            array[:, INDEX["event_id"]], metadata_number(run.metadata, "seed", int), n_blocks
        ))
    block_ids = np.concatenate(block_ids_parts)
    n_by_block = np.bincount(block_ids, weights=delta_tr, minlength=n_blocks)
    d_by_block = np.bincount(block_ids, weights=delta, minlength=n_blocks)
    nk_by_block = np.empty((n_blocks, len(COEFFICIENTS)))
    dk_by_block = np.empty_like(nk_by_block)
    for idx in range(len(COEFFICIENTS)):
        nk_by_block[:, idx] = np.bincount(
            block_ids, weights=delta_tr * kernel[:, idx], minlength=n_blocks
        )
        dk_by_block[:, idx] = np.bincount(
            block_ids, weights=delta * kernel[:, idx], minlength=n_blocks
        )

    rng = np.random.default_rng(bootstrap_seed)
    counts = rng.multinomial(n_blocks, np.full(n_blocks, 1.0 / n_blocks), size=n_bootstrap)
    boot_n = counts @ n_by_block
    boot_d = counts @ d_by_block
    boot_f = boot_n / boot_d
    boot_nk = counts @ nk_by_block
    boot_dk = counts @ dk_by_block
    boot_beta = (boot_nk - boot_f[:, None] * boot_dk) / boot_n[:, None]
    boot_cm = boot_f / baseline_factor

    beta_se = boot_beta.std(axis=0, ddof=1)
    beta_ci = np.quantile(boot_beta, (0.025, 0.975), axis=0)
    f_se = float(boot_f.std(ddof=1))
    f_ci = np.quantile(boot_f, (0.025, 0.975))
    cm_se = float(boot_cm.std(ddof=1))
    cm_ci = np.quantile(boot_cm, (0.025, 0.975))

    coefficient_rows = []
    for idx, name in enumerate(COEFFICIENTS):
        absolute_contribution = np.abs(residual * kernel[:, idx])
        tail_count = max(1, int(math.ceil(0.001 * len(events))))
        if absolute_contribution.sum() > 0.0:
            largest = np.partition(absolute_contribution, len(events) - tail_count)[-tail_count:]
            tail_fraction = float(largest.sum() / absolute_contribution.sum())
        else:
            tail_fraction = 0.0
        half_width = 0.5 * float(beta_ci[1, idx] - beta_ci[0, idx])
        threshold = max(0.05, 0.15 * abs(float(coefficients[idx]))) if name in CORE_COEFFICIENTS else max(
            0.10, 0.25 * abs(float(coefficients[idx]))
        )
        coefficient_rows.append({
            "name": name,
            "group": "core" if name in CORE_COEFFICIENTS else "higher_order_candidate",
            "estimate": float(coefficients[idx]),
            "standard_error": float(beta_se[idx]),
            "ci_low": float(beta_ci[0, idx]),
            "ci_high": float(beta_ci[1, idx]),
            "ci_half_width": half_width,
            "precision_threshold": threshold,
            "precision_pass": bool(block_population_pass and half_width <= threshold),
            "top_0p1pct_absolute_contribution_fraction": tail_fraction,
            "tail_pass": bool(tail_fraction <= 0.10),
        })

    projection_diagnostics = {}
    for name, values in directional.items():
        estimates = np.einsum("n,nk->k", residual, values, optimize=True) / numerator

        directional_n = np.empty((n_blocks, values.shape[1]))
        directional_d = np.empty_like(directional_n)
        for direction in range(values.shape[1]):
            directional_n[:, direction] = np.bincount(
                block_ids, weights=delta_tr * values[:, direction], minlength=n_blocks
            )
            directional_d[:, direction] = np.bincount(
                block_ids, weights=delta * values[:, direction], minlength=n_blocks
            )
        boot_directional = (
            counts @ directional_n - boot_f[:, None] * (counts @ directional_d)
        ) / boot_n[:, None]
        boot_contrast = boot_directional - boot_directional.mean(axis=1, keepdims=True)
        contrast_ci = np.quantile(boot_contrast, (0.025, 0.975), axis=0)
        agreement = np.logical_and(contrast_ci[0] <= 0.0, contrast_ci[1] >= 0.0)
        projection_diagnostics[name] = {
            "direction_estimates": [float(value) for value in estimates],
            "direction_standard_errors": [
                float(value) for value in boot_directional.std(axis=0, ddof=1)
            ],
            "standard_deviation": float(np.std(estimates, ddof=1)),
            "centered_difference_ci_low": [float(value) for value in contrast_ci[0]],
            "centered_difference_ci_high": [float(value) for value in contrast_ci[1]],
            "agreement_pass": bool(np.all(agreement)),
        }

    loo_max = None
    loo_f_m_max = None
    if len(runs) > 1:
        loo_max = 0.0
        loo_f_m_max = 0.0
        boundaries = np.cumsum([0] + [len(array) for array in arrays])
        for shard_index in range(len(runs)):
            mask = np.ones(len(events), dtype=bool)
            mask[boundaries[shard_index]:boundaries[shard_index + 1]] = False
            n_loo = float(delta_tr[mask].sum())
            d_loo = float(delta[mask].sum())
            f_loo = n_loo / d_loo
            beta_loo = np.einsum(
                "n,na->a", delta_tr[mask] - f_loo * delta[mask], kernel[mask], optimize=True
            ) / n_loo
            loo_max = max(loo_max, float(np.max(np.abs(beta_loo - coefficients))))
            loo_f_m_max = max(loo_f_m_max, abs(f_loo - f_m))

    baseline_half_width = 0.5 * float(f_ci[1] - f_ci[0])
    baseline_precision_pass = block_population_pass and baseline_half_width <= 0.01 * abs(f_m)
    continuation_required = (
        not baseline_precision_pass
        or any(not row["precision_pass"] for row in coefficient_rows)
    )
    return {
        "schema_version": 1,
        "alpha": alpha,
        "theta": theta,
        "aspect_ratio": metadata_number(reference, "aspect_ratio"),
        "n_events": int(len(events)),
        "n_shards": len(runs),
        "baseline": {
            "B_theta": baseline_factor,
            "f_M": f_m,
            "f_M_standard_error": f_se,
            "f_M_ci_low": float(f_ci[0]),
            "f_M_ci_high": float(f_ci[1]),
            "C_M": c_m,
            "C_M_standard_error": cm_se,
            "C_M_ci_low": float(cm_ci[0]),
            "C_M_ci_high": float(cm_ci[1]),
            "precision_pass": bool(baseline_precision_pass),
        },
        "coefficients": coefficient_rows,
        "projection_diagnostics": projection_diagnostics,
        "projection_agreement_pass": bool(all(
            item["agreement_pass"] for item in projection_diagnostics.values()
        )),
        "leave_one_shard_out_max_absolute_change": loo_max,
        "leave_one_shard_out_f_M_max_absolute_change": loo_f_m_max,
        "continuation_required": bool(continuation_required),
        "validation": validations,
        "bootstrap": {
            "blocks": n_blocks,
            "requested_blocks": requested_n_blocks,
            "block_population_pass": block_population_pass,
            "replicates": n_bootstrap,
            "seed": bootstrap_seed,
        },
        "source_directories": [str(run.directory) for run in runs],
    }


def write_json(path: Path | str, value: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
