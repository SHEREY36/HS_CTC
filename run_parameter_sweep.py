#!/usr/bin/env python3
"""
Parameter sweep orchestration for HS_CTC simulations.

Sweeps over:
  - Coefficient of restitution (alpha_pp)
  - Temperature ratio r = kTm/kTI  (kTI fixed at 1.0, kTm varies)
  - Aspect ratio (AR = (LCYL + DIA) / DIA)

Temperature parameterisation:
  kTI = 1.0  (fixed, avoids singularity at kTI=0)
  kTm = r    where r in [0.1, 0.2, ..., 2.0] (step 0.1)
  => ratio r = kTm/kTI = kTm
"""

import subprocess
import itertools
import os
from pathlib import Path
import numpy as np
import argparse

# -----------------------------------------------------------------------
# Default parameter ranges
# -----------------------------------------------------------------------
ALPHA_PP_VALUES = np.round(np.linspace(0.5, 1.0, 11), 3)   # 0.5, 0.55, ..., 1.0  (11 values)
KTM_VALUES      = np.round(np.arange(0.1, 2.01, 0.1), 2)   # 0.1, 0.2, ..., 2.0   (20 values)
KTI             = 1.0                                        # fixed
AR_VALUES       = [1.0, 1.5, 2.0, 3.0, 4.0]                # aspect ratios to sweep

# Configuration
EXECUTABLE      = "./SphCyl"
BASE_OUTPUT_DIR = "results"


def run_simulation(alpha, ktm, ar, output_dir, num_threads=20):
    """Run one simulation; kTI is always 1.0."""
    os.makedirs(output_dir, exist_ok=True)

    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(num_threads)

    # SphCyl <alpha> <kTm> <kTI> <AR> <output_dir>
    cmd = [EXECUTABLE,
           f"{alpha:.4f}",
           f"{ktm:.4f}",
           f"{KTI:.4f}",
           f"{ar:.4f}",
           output_dir]

    print(f"  alpha={alpha:.3f}  kTm={ktm:.2f}  kTI={KTI:.1f}  AR={ar:.2f}  "
          f"r={ktm/KTI:.2f}  -> {output_dir}")

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)

    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.decode().strip()}")
        return False

    return True


def main():
    parser = argparse.ArgumentParser(description='Run HS_CTC parameter sweep')
    parser.add_argument('--threads',    type=int,   default=20,              help='OpenMP threads per run')
    parser.add_argument('--output-dir', type=str,   default=BASE_OUTPUT_DIR, help='Base output directory')
    parser.add_argument('--alpha',      type=float, nargs='+', default=None,
                        help='Override alpha values (space separated)')
    parser.add_argument('--ar',         type=float, nargs='+', default=None,
                        help='Override AR values (space separated)')
    args = parser.parse_args()

    alpha_vals = args.alpha if args.alpha else ALPHA_PP_VALUES
    ar_vals    = args.ar    if args.ar    else AR_VALUES

    Path(args.output_dir).mkdir(exist_ok=True)

    combos = list(itertools.product(alpha_vals, KTM_VALUES, ar_vals))
    total  = len(combos)

    print(f"Parameter sweep: {len(alpha_vals)} alpha × {len(KTM_VALUES)} kTm × {len(ar_vals)} AR = {total} runs")
    print(f"kTI fixed at {KTI}  (temperature ratio r = kTm/kTI)\n")

    failed = []
    for i, (alpha, ktm, ar) in enumerate(combos, 1):
        print(f"[{i:4d}/{total}]", end="  ")
        out = f"{args.output_dir}/alpha_{alpha:.3f}_r{ktm/KTI:.2f}_AR{ar:.2f}"
        ok = run_simulation(alpha, ktm, ar, out, args.threads)
        if not ok:
            failed.append((alpha, ktm, ar))

    print(f"\nDone.  {total - len(failed)}/{total} succeeded.")
    if failed:
        print("Failed runs:")
        for p in failed:
            print(f"  alpha={p[0]:.3f}  kTm={p[1]:.2f}  AR={p[2]:.2f}")
    else:
        print(f"Results in: {args.output_dir}/")


if __name__ == "__main__":
    main()
