from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from closure_common import (  # noqa: E402
    COLUMNS,
    COEFFICIENTS,
    TENSOR_BASIS,
    estimate_from_runs,
    load_run,
    validate_run,
)
from make_sweep_manifest import ALPHAS, ASPECT_RATIOS, THETAS, seed_for  # noqa: E402
from aggregate_closure import weighted_limit  # noqa: E402


def base_states(rng: np.random.Generator, n: int):
    c = rng.normal(scale=math.sqrt(0.5), size=(n, 3))
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    omega = rng.normal(scale=math.sqrt(0.5), size=(n, 3))
    omega -= np.einsum("ij,ij->i", omega, u)[:, None] * u
    return c, omega, u


class ScoreNormalizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rng = np.random.default_rng(4102)
        cls.c, cls.omega, cls.u = base_states(cls.rng, 300_000)

    def test_sonine_scores_are_dual_to_named_cumulants(self):
        def gamma_expectation(poly, shape):
            return sum(
                coefficient * math.gamma(shape + power) / math.gamma(shape)
                for power, coefficient in enumerate(poly)
            )

        def norm(poly, shape):
            return gamma_expectation(np.polynomial.polynomial.polymul(poly, poly), shape)

        s1 = np.array([1.5, -1.0])
        s2 = np.array([15.0 / 8.0, -2.5, 0.5])
        s3 = np.array([35.0 / 16.0, -35.0 / 8.0, 7.0 / 4.0, -1.0 / 6.0])
        l1 = np.array([1.0, -1.0])
        l2 = np.array([1.0, -2.0, 0.5])
        l3 = np.array([1.0, -3.0, 1.5, -1.0 / 6.0])
        self.assertAlmostEqual(norm(s1, 1.5), 1.5, places=12)
        self.assertAlmostEqual(norm(s2, 1.5), 15.0 / 8.0, places=12)
        self.assertAlmostEqual(norm(s3, 1.5), 35.0 / 16.0, places=11)
        self.assertAlmostEqual(norm(l1, 1.0), 1.0, places=12)
        self.assertAlmostEqual(norm(l2, 1.0), 1.0, places=12)
        self.assertAlmostEqual(norm(l3, 1.0), 1.0, places=11)
        self.assertAlmostEqual(norm(s2, 1.5) * norm(l1, 1.0) / (15.0 / 8.0), 1.0, places=12)
        self.assertAlmostEqual(norm(s1, 1.5) * norm(l2, 1.0) / 1.5, 1.0, places=12)

    def test_spin_alignment_dual_scores(self):
        for a in TENSOR_BASIS:
            cac = np.einsum("ni,ij,nj->n", self.c, a, self.c)
            wow = np.einsum("ni,ij,nj->n", self.omega, a, self.omega)
            uau = np.einsum("ni,ij,nj->n", self.u, a, self.u)
            score_r = (10.0 / 7.0) * wow + (5.0 / 7.0) * uau
            score_q = (10.0 / 7.0) * wow + (40.0 / 7.0) * uau

            dpi = 2.0 * np.einsum("ni,nj,n->ij", self.c, self.c, cac) / len(cac)

            def responses(score):
                dr = 3.0 * np.einsum("ni,nj,n->ij", self.omega, self.omega, score) / len(score)
                dq = 1.5 * np.einsum("ni,nj,n->ij", self.u, self.u, score) / len(score)
                return dr, dq

            dr_r, dq_r = responses(score_r)
            dr_q, dq_q = responses(score_q)
            self.assertLess(np.linalg.norm(dpi - a), 0.035)
            self.assertLess(np.linalg.norm(dr_r - a), 0.04)
            self.assertLess(np.linalg.norm(dq_r), 0.04)
            self.assertLess(np.linalg.norm(dr_q), 0.065)
            self.assertLess(np.linalg.norm(dq_q - a), 0.065)

    def test_heat_flux_scores_are_dual_in_all_directions(self):
        x = np.einsum("ij,ij->i", self.c, self.c)
        y = np.einsum("ij,ij->i", self.omega, self.omega)
        score_tr = self.c * (x[:, None] - 2.5)
        score_rot = self.c * (y[:, None] - 1.0)
        response_tr = 0.8 * score_tr.T @ score_tr / len(score_tr)
        response_rot = 2.0 * score_rot.T @ score_rot / len(score_rot)
        cross = score_tr.T @ score_rot / len(score_tr)
        np.testing.assert_allclose(response_tr, np.eye(3), atol=0.035)
        np.testing.assert_allclose(response_rot, np.eye(3), atol=0.035)
        np.testing.assert_allclose(cross, np.zeros((3, 3)), atol=0.025)


class PipelineTests(unittest.TestCase):
    def make_run(self, root: Path, name: str, seed: int, n: int = 1500):
        rng = np.random.default_rng(seed)
        c1, omega1, u1 = base_states(rng, n)
        c2, omega2, u2 = base_states(rng, n)
        delta = rng.uniform(0.05, 0.25, n)
        delta_tr = delta * (0.58 + 0.03 * (c1[:, 0] + c2[:, 0]))
        delta_rot = delta - delta_tr
        e0 = np.full(n, 10.0)
        et_el = np.full(n, 6.0)
        er_el = np.full(n, 4.0)
        normal = rng.normal(size=(n, 3)); normal /= np.linalg.norm(normal, axis=1)[:, None]
        centre = rng.normal(size=(n, 3)); centre /= np.linalg.norm(centre, axis=1)[:, None]
        data = np.zeros((n, len(COLUMNS)))
        event_id_offset = seed * 10_000_000
        data[:, 0] = event_id_offset + np.arange(1, n + 1)
        data[:, 1:4] = c1; data[:, 4:7] = c2
        data[:, 7:10] = omega1; data[:, 10:13] = omega2
        data[:, 13:16] = u1; data[:, 16:19] = u2
        data[:, 19] = delta_tr; data[:, 20] = delta_rot; data[:, 21] = delta
        data[:, 22] = et_el; data[:, 23] = er_el
        data[:, 24] = et_el - delta_tr; data[:, 25] = er_el - delta_rot
        data[:, 26] = e0; data[:, 27] = 0.0; data[:, 28] = 1.0
        data[:, 29] = rng.uniform(0.0, 1.0, n)
        data[:, 30:33] = normal; data[:, 33:36] = centre
        directory = root / name; directory.mkdir()
        data.astype("<f8").tofile(directory / "closure_events.bin")
        metadata = {
            "schema_version": "1", "dtype": "<f8", "n_columns": str(len(COLUMNS)),
            "columns": ",".join(COLUMNS), "alpha": "0.8", "theta": "1.0",
            "aspect_ratio": "2.0", "velocity_scale": "1.0", "omega_scale": "1.0",
            "mass": "1.0", "moi_perpendicular": "1.0", "nsamples": str(n),
            "seed": str(seed), "event_id_offset": str(event_id_offset),
            "output_mode": "closure",
        }
        (directory / "metadata.txt").write_text("".join(f"{k}={v}\n" for k, v in metadata.items()))
        return load_run(directory)

    def test_schema_validation_and_estimation_are_shard_order_invariant(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run1 = self.make_run(root, "shard_000", 51)
            run2 = self.make_run(root, "shard_001", 92)
            self.assertEqual(validate_run(run1)["status"], "pass")
            first = estimate_from_runs([run1, run2], n_blocks=32, n_bootstrap=80, bootstrap_seed=7)
            second = estimate_from_runs([run2, run1], n_blocks=32, n_bootstrap=80, bootstrap_seed=7)
            self.assertAlmostEqual(first["baseline"]["f_M"], second["baseline"]["f_M"], places=13)
            b1 = np.array([row["estimate"] for row in first["coefficients"]])
            b2 = np.array([row["estimate"] for row in second["coefficients"]])
            np.testing.assert_allclose(b1, b2, rtol=2e-12, atol=2e-12)
            self.assertEqual(len(first["coefficients"]), len(COEFFICIENTS))
            self.assertIn("projection_agreement_pass", first)

    def test_corrupt_and_alpha_one_runs_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run = self.make_run(root, "valid", 63, n=100)
            run.metadata["alpha"] = "1.0"
            with self.assertRaisesRegex(ValueError, "alpha=1"):
                estimate_from_runs([run], n_blocks=8, n_bootstrap=10)

            corrupt = root / "corrupt"
            corrupt.mkdir()
            (corrupt / "metadata.txt").write_text((run.directory / "metadata.txt").read_text())
            payload = (run.directory / "closure_events.bin").read_bytes()
            (corrupt / "closure_events.bin").write_bytes(payload[:-1])
            with self.assertRaisesRegex(ValueError, "truncated"):
                load_run(corrupt)

    def test_grid_and_common_random_seed_contract(self):
        self.assertEqual(len(ALPHAS) * len(THETAS) * len(ASPECT_RATIOS), 1200)
        self.assertEqual(seed_for(2, 3, 4), seed_for(2, 3, 4))
        self.assertNotEqual(seed_for(2, 3, 4), seed_for(2, 3, 5))

    def test_alpha_limit_recovers_a_linear_near_elastic_trend(self):
        x = 1.0 - np.array([0.90, 0.95, 0.975, 0.99])
        y = 2.75 - 1.2 * x
        estimate, standard_error, low, high, statistical, systematic = weighted_limit(
            x, y, np.full(4, 0.02)
        )
        self.assertAlmostEqual(estimate, 2.75, places=12)
        self.assertLess(low, estimate)
        self.assertGreater(high, estimate)
        self.assertGreater(standard_error, 0.0)
        self.assertGreater(statistical, 0.0)
        self.assertAlmostEqual(systematic, 0.0, places=12)


if __name__ == "__main__":
    unittest.main()
