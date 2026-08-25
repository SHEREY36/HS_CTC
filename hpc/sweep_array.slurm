#!/bin/bash
#SBATCH --job-name=hs_ctc_sweep
#SBATCH --account=morri353
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=06:00:00
#SBATCH --output=logs/sweep_%A_%a.out
#SBATCH --error=logs/sweep_%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
: "${MANIFEST:?submit through hpc/submit_manifest.sh or export MANIFEST}"
TASK_OFFSET=${TASK_OFFSET:-0}
GLOBAL_TASK_ID=$((TASK_OFFSET + SLURM_ARRAY_TASK_ID))
ROW_NUMBER=$((GLOBAL_TASK_ID + 2))

row=$(awk -F, -v row_number="$ROW_NUMBER" 'NR == row_number {print; exit}' "$MANIFEST")
row=${row%$'\r'}
if [[ -z "$row" ]]; then
    echo "No manifest row for global task $GLOBAL_TASK_ID" >&2
    exit 2
fi
IFS=, read -r task_index stage alpha theta aspect_ratio seed nsamples shard_id node_key output_dir <<<"$row"
if [[ "$task_index" != "$GLOBAL_TASK_ID" ]]; then
    echo "Manifest task_index=$task_index does not match global task $GLOBAL_TASK_ID" >&2
    exit 2
fi

if [[ -f "$output_dir/_SUCCESS" ]]; then
    echo "Already validated; skipping $output_dir"
    exit 0
fi
mkdir -p "$output_dir"

command -v gfortran >/dev/null
command -v python3 >/dev/null
[[ -x build/SphCyl ]] || { echo "build/SphCyl is missing; run make -C build clean all" >&2; exit 2; }
[[ -f system_input.dat ]] || { echo "system_input.dat is missing" >&2; exit 2; }

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

echo "node=$node_key shard=$shard_id alpha=$alpha theta=$theta AR=$aspect_ratio samples=$nsamples seed=$seed"
srun ./build/SphCyl "$alpha" "$theta" 1.0 "$aspect_ratio" "$output_dir" "$seed" "$nsamples" closure
python3 scripts/validate_closure_run.py "$output_dir" --mark-success
