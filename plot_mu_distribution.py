#!/usr/bin/env python3
"""Plot the mu_in distribution from an HS_CTC run.

Reads column 4 of chi.txt, overlays the theoretical isotropic PDF p(mu)=2mu
on [0, 1], and also plots the empirical CDF against F(mu)=mu^2.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

# Matplotlib needs a writable config/cache directory on this machine.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot mu_in distribution from chi.txt (column 4)."
    )
    parser.add_argument(
        "run_dir",
        help="Run directory that contains chi.txt, e.g. mu_check_AR2_alpha1_N10k",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=30,
        help="Number of histogram bins (default: 30)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively in addition to saving it",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    chi_path = run_dir / "chi.txt"
    out_path = run_dir / "mu_distribution.png"

    if not chi_path.exists():
        raise SystemExit(f"Missing file: {chi_path}")

    data = np.loadtxt(chi_path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 4:
        raise SystemExit(
            f"{chi_path} has {data.shape[1]} columns; expected at least 4 with mu_in in column 4."
        )

    mu = data[:, 3]
    if not np.all(np.isfinite(mu)):
        raise SystemExit("mu_in contains NaN or Inf values.")
    if np.any((mu < 0.0) | (mu > 1.0)):
        raise SystemExit("mu_in contains values outside [0, 1].")

    mu_sorted = np.sort(mu)
    n = mu.size
    mu_grid = np.linspace(0.0, 1.0, 400)
    expected_pdf = 2.0 * mu_grid
    expected_cdf = mu_grid**2
    empirical_cdf = np.arange(1, n + 1) / n

    mean_mu = mu.mean()
    std_mu = mu.std(ddof=1) if n > 1 else 0.0

    fig, (ax_pdf, ax_cdf) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax_pdf.hist(
        mu,
        bins=args.bins,
        range=(0.0, 1.0),
        density=True,
        alpha=0.65,
        color="#4C72B0",
        edgecolor="black",
        linewidth=0.6,
        label="CTC histogram",
    )
    ax_pdf.plot(mu_grid, expected_pdf, color="#C44E52", linewidth=2.2, label=r"Theory: $p(\mu)=2\mu$")
    ax_pdf.set_xlabel(r"$\mu$")
    ax_pdf.set_ylabel("Density")
    ax_pdf.set_title("PDF Check")
    ax_pdf.set_xlim(0.0, 1.0)
    ax_pdf.legend(frameon=False)

    ax_cdf.plot(mu_sorted, empirical_cdf, color="#4C72B0", linewidth=2.0, label="Empirical CDF")
    ax_cdf.plot(mu_grid, expected_cdf, color="#C44E52", linewidth=2.2, linestyle="--", label=r"Theory: $F(\mu)=\mu^2$")
    ax_cdf.set_xlabel(r"$\mu$")
    ax_cdf.set_ylabel("CDF")
    ax_cdf.set_title("CDF Check")
    ax_cdf.set_xlim(0.0, 1.0)
    ax_cdf.set_ylim(0.0, 1.0)
    ax_cdf.legend(frameon=False)

    fig.suptitle(
        f"mu_in sanity check: N={n}, mean={mean_mu:.4f}, expected mean=0.6667",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")

    print(f"Saved plot to: {out_path}")
    print(f"N = {n}")
    print(f"mean(mu) = {mean_mu:.6f}   expected = 0.666667")
    print(f"std(mu)  = {std_mu:.6f}")
    print(f"min(mu)  = {mu.min():.6f}")
    print(f"max(mu)  = {mu.max():.6f}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
