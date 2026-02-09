#!/bin/bash
#SBATCH --job-name=hs_ctc_sweep
#SBATCH --array=0-54  # 11 alpha * 5 kti = 55 jobs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH -A morri353

# Create logs directory
mkdir -p logs

# Parameter arrays (must match run_parameter_sweep.py)
ALPHA_VALUES=(0.500 0.550 0.600 0.650 0.700 0.750 0.800 0.850 0.900 0.950 1.000)
KTM=1.0
KTI_VALUES=(0.0 0.5 1.0 1.5 2.0)

# Calculate which parameters to use for this array task
N_ALPHA=${#ALPHA_VALUES[@]}
N_KTI=${#KTI_VALUES[@]}

ALPHA_IDX=$((SLURM_ARRAY_TASK_ID / N_KTI))
KTI_IDX=$((SLURM_ARRAY_TASK_ID % N_KTI))

ALPHA=${ALPHA_VALUES[$ALPHA_IDX]}
KTI=${KTI_VALUES[$KTI_IDX]}

# Create output directory
OUTPUT_DIR="results/alpha_${ALPHA}_ktm_${KTM}_kti_${KTI}"
mkdir -p $OUTPUT_DIR

# Set OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Job array task: $SLURM_ARRAY_TASK_ID"
echo "Parameters: alpha=$ALPHA, kTm=$KTM, kTI=$KTI"
echo "Threads: $OMP_NUM_THREADS"
echo "Output: $OUTPUT_DIR"

# Run simulation
./SphCyl $ALPHA $KTM $KTI $OUTPUT_DIR

echo "Completed: $OUTPUT_DIR"
