#!/usr/bin/env python3
"""
Parameter sweep orchestration for HS_CTC simulations.
Iterates over coefficient of restitution (alpha) and temperature ratios (kTm/kTI).
"""

import subprocess
import itertools
import os
from pathlib import Path
import numpy as np
import argparse

# Default parameter ranges
ALPHA_PP_VALUES = np.linspace(0.5, 1.0, 11)  # 0.5, 0.55, ..., 1.0
KTM_VALUES = [1.0]  # Fixed translational temperature
KTI_VALUES = [0.0, 0.5, 1.0, 1.5, 2.0]  # Temperature ratios

# Configuration
EXECUTABLE = "./SphCyl"
BASE_OUTPUT_DIR = "results"

def run_simulation(alpha, ktm, kti, output_dir, num_threads=20):
    """Run single simulation with given parameters"""
    os.makedirs(output_dir, exist_ok=True)

    # Set environment
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(num_threads)

    # Run with command-line args: SphCyl <alpha> <ktm> <kti> <output_dir>
    cmd = [EXECUTABLE, str(alpha), str(ktm), str(kti), output_dir]

    print(f"Running: alpha={alpha:.3f}, kTm={ktm}, kTI={kti}, threads={num_threads}")
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env
    )

    if result.returncode != 0:
        print(f"ERROR: {result.stderr.decode()}")
        return False

    print(f"  Completed: {output_dir}")
    return True

def main():
    parser = argparse.ArgumentParser(description='Run HS_CTC parameter sweep')
    parser.add_argument('--threads', type=int, default=20, help='Number of OpenMP threads')
    parser.add_argument('--output-dir', type=str, default=BASE_OUTPUT_DIR, help='Base output directory')
    args = parser.parse_args()

    # Create base output directory
    Path(args.output_dir).mkdir(exist_ok=True)

    # Generate all parameter combinations
    param_combinations = list(itertools.product(
        ALPHA_PP_VALUES,
        KTM_VALUES,
        KTI_VALUES
    ))

    print(f"Total simulations: {len(param_combinations)}")
    print(f"OpenMP threads per simulation: {args.threads}")

    for i, (alpha, ktm, kti) in enumerate(param_combinations, 1):
        output_dir = f"{args.output_dir}/alpha_{alpha:.3f}_ktm_{ktm:.2f}_kti_{kti:.2f}"

        print(f"\n[{i}/{len(param_combinations)}]")
        success = run_simulation(alpha, ktm, kti, output_dir, args.threads)

        if not success:
            print(f"Failed: {output_dir}")
            continue

    print("\n\nAll simulations complete!")
    print(f"Results in: {args.output_dir}/")

if __name__ == "__main__":
    main()
