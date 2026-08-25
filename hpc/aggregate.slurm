#!/bin/bash
#SBATCH --job-name=hs_ctc_aggregate
#SBATCH --account=morri353
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --error=logs/aggregate_%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"

args=(
    --coefficients-root "${COEFFICIENTS_ROOT:-coefficients}"
    --results-root "${RESULTS_ROOT:-results}"
    --output-dir "${AGGREGATE_OUTPUT_DIR:-closure_output}"
)
if [[ "${REQUIRE_COMPLETE:-0}" == 1 ]]; then
    args+=(--require-complete)
fi
python3 scripts/aggregate_closure.py "${args[@]}"
