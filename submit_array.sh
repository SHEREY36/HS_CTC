#!/bin/bash
#
# SLURM job array for HS_CTC parameter sweep.
#
# Parameter space:
#   11 alpha  x  20 kTm  x  5 AR  =  1100 jobs
#   kTI fixed at 1.0;  r = kTm/kTI  in [0.1 .. 2.0]
#
#SBATCH --job-name=hs_ctc_sweep
#SBATCH --array=0-1099
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=6:00:00
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH -A morri353

mkdir -p logs

# -----------------------------------------------------------------------
# Parameter arrays  (must match run_parameter_sweep.py)
# -----------------------------------------------------------------------
ALPHA_VALUES=(0.500 0.550 0.600 0.650 0.700 0.750 0.800 0.850 0.900 0.950 1.000)
KTM_VALUES=(0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 \
            1.10 1.20 1.30 1.40 1.50 1.60 1.70 1.80 1.90 2.00)
KTI=1.0
AR_VALUES=(1.0 1.5 2.0 3.0 4.0)

N_ALPHA=${#ALPHA_VALUES[@]}   # 11
N_KTM=${#KTM_VALUES[@]}       # 20
N_AR=${#AR_VALUES[@]}         # 5

# 3D index mapping:
#   task_id = ALPHA_IDX * (N_KTM * N_AR) + KTM_IDX * N_AR + AR_IDX
ALPHA_IDX=$(( SLURM_ARRAY_TASK_ID / (N_KTM * N_AR) ))
REMAINDER=$(( SLURM_ARRAY_TASK_ID % (N_KTM * N_AR) ))
KTM_IDX=$(( REMAINDER / N_AR ))
AR_IDX=$(( REMAINDER % N_AR ))

ALPHA=${ALPHA_VALUES[$ALPHA_IDX]}
KTM=${KTM_VALUES[$KTM_IDX]}
AR=${AR_VALUES[$AR_IDX]}

# Temperature ratio r = kTm/kTI
R=$(echo "$KTM / $KTI" | bc -l | xargs printf "%.2f")

# -----------------------------------------------------------------------
OUTPUT_DIR="results/alpha_${ALPHA}_r${R}_AR${AR}"
mkdir -p "$OUTPUT_DIR"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Task ${SLURM_ARRAY_TASK_ID}: alpha=${ALPHA}  kTm=${KTM}  kTI=${KTI}  r=${R}  AR=${AR}"
echo "Threads: $OMP_NUM_THREADS   Output: $OUTPUT_DIR"

# SphCyl <alpha> <kTm> <kTI> <AR> <output_dir>
./SphCyl "$ALPHA" "$KTM" "$KTI" "$AR" "$OUTPUT_DIR"

echo "Completed: $OUTPUT_DIR"
