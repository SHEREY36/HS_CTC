#!/bin/bash
#SBATCH --job-name=hs_ctc_estimate
#SBATCH --account=morri353
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/estimate_%A_%a.out
#SBATCH --error=logs/estimate_%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
: "${MANIFEST:?submit through hpc/submit_manifest.sh or export MANIFEST}"
TASK_OFFSET=${TASK_OFFSET:-0}
GLOBAL_TASK_ID=$((TASK_OFFSET + SLURM_ARRAY_TASK_ID))
ROW_NUMBER=$((GLOBAL_TASK_ID + 2))

row=$(awk -F, -v row_number="$ROW_NUMBER" 'NR == row_number {print; exit}' "$MANIFEST")
row=${row%$'\r'}
if [[ -z "$row" ]]; then
    echo "No node-manifest row for global task $GLOBAL_TASK_ID" >&2
    exit 2
fi
IFS=, read -r node_index alpha theta aspect_ratio node_key node_dir coefficient_dir <<<"$row"
if [[ "$node_index" != "$GLOBAL_TASK_ID" ]]; then
    echo "Node manifest index mismatch" >&2
    exit 2
fi

mkdir -p "$coefficient_dir"
rm -f "$coefficient_dir/_SUCCESS"
python3 scripts/estimate_closure.py \
    --node-dir "$node_dir" \
    --output-dir "$coefficient_dir" \
    --blocks 128 --bootstrap 2000 --bootstrap-seed 20260825
