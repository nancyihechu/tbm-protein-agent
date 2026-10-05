#!/bin/bash
#SBATCH --job-name=af3_test_target
#SBATCH --partition=dkhf3-lab-gpu
#SBATCH --account=dkhf3-lab
#SBATCH --gres=gpu:L40S:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/cluster/pixstor/dkhf3-lab/users/nancy/TBM_AGENT/examples/test_target/logs/slurm-%j.out

set -eo pipefail

EX=/cluster/pixstor/dkhf3-lab/users/nancy/TBM_AGENT/examples/test_target
SIF=/cluster/software/src/SINGULARITY_IMAGES/alphafold3/NEW_VERSION_NO_MEMORY_ISSUES/alphafold3.sif
INPUT=$EX/af3/test_target.json
OUT=$EX/af3/output

mkdir -p "$OUT"
export TMPDIR=$EX/af3/tmp_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"

echo "[$(date)] Host: $(hostname)  GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader)"

# Official DeepMind AlphaFold3 container (3.0.0); MSA search + inference
singularity exec --nv -B /cluster:/cluster --env TMPDIR=$TMPDIR "$SIF" \
    python /app/alphafold/run_alphafold.py \
        --json_path="$INPUT" \
        --model_dir=/cluster/VAST/cietest/alphafold3_model_parameters \
        --db_dir=/cluster/VAST/cietest/alphafold3_database \
        --output_dir="$OUT" \
        --jackhmmer_n_cpu=8 \
        --nhmmer_n_cpu=8

echo "[$(date)] AlphaFold3 returned $?"
rm -rf "$TMPDIR"
ls -lR "$OUT"
