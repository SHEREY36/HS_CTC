#!/bin/bash
# Compatibility entry point for the standalone closure sweep.
set -euo pipefail

stage=${1:-pilot}
case "$stage" in
    pilot) samples=20000; shard=0 ;;
    production) samples=80000; shard=1 ;;
    *) echo "Usage: bash submit_array.sh [pilot|production]" >&2; exit 2 ;;
esac

mkdir -p manifests logs results coefficients
manifest="manifests/${stage}.csv"
python3 scripts/make_sweep_manifest.py \
    --stage "$stage" --samples "$samples" --shard-id "$shard" --output "$manifest"
exec bash hpc/submit_manifest.sh "$manifest" hpc/sweep_array.slurm 50
