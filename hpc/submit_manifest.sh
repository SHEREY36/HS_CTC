#!/bin/bash
# Submit a CSV manifest in Slurm-compatible array chunks.

set -euo pipefail

usage() {
    echo "Usage: $0 MANIFEST [SLURM_SCRIPT] [CONCURRENCY] [DEPENDENCY_JOB_ID]" >&2
}

[[ $# -ge 1 ]] || { usage; exit 2; }
MANIFEST_INPUT=$1
SLURM_SCRIPT=${2:-hpc/sweep_array.slurm}
CONCURRENCY=${3:-50}
DEPENDENCY_JOB_ID=${4:-}

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT_DIR"
mkdir -p logs manifests

if [[ "$MANIFEST_INPUT" = /* ]]; then
    MANIFEST=$MANIFEST_INPUT
else
    MANIFEST=$ROOT_DIR/$MANIFEST_INPUT
fi
[[ -f "$MANIFEST" ]] || { echo "Manifest not found: $MANIFEST" >&2; exit 2; }
[[ -f "$SLURM_SCRIPT" ]] || { echo "Slurm script not found: $SLURM_SCRIPT" >&2; exit 2; }

TASK_COUNT=$(( $(wc -l < "$MANIFEST") - 1 ))
(( TASK_COUNT > 0 )) || { echo "Manifest has no tasks: $MANIFEST" >&2; exit 2; }

MAX_ARRAY_SIZE=$(scontrol show config 2>/dev/null | awk -F= '/MaxArraySize/ {gsub(/[[:space:]]/, "", $2); print $2; exit}')
MAX_ARRAY_SIZE=${MAX_ARRAY_SIZE:-1001}
CHUNK_SIZE=$MAX_ARRAY_SIZE
JOB_LOG="manifests/submitted_jobs_$(date +%Y%m%d_%H%M%S).txt"

offset=0
while (( offset < TASK_COUNT )); do
    remaining=$((TASK_COUNT - offset))
    count=$CHUNK_SIZE
    (( remaining < count )) && count=$remaining
    last=$((count - 1))
    args=(--parsable --array="0-${last}%${CONCURRENCY}" \
        --export="ALL,MANIFEST=${MANIFEST},TASK_OFFSET=${offset}")
    if [[ -n "$DEPENDENCY_JOB_ID" ]]; then
        args+=(--dependency="afterok:${DEPENDENCY_JOB_ID}")
    fi
    job_id=$(sbatch "${args[@]}" "$SLURM_SCRIPT")
    echo "$job_id offset=$offset count=$count manifest=$MANIFEST script=$SLURM_SCRIPT" | tee -a "$JOB_LOG"
    offset=$((offset + count))
done

echo "Submission record: $JOB_LOG"
