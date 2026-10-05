#!/bin/bash
# AlphaFold3 for one target: sbatch run_af3_target.slurm <TARGET>
# Input:  targets/<TARGET>/af3_model/<TARGET>.json
# Output: targets/<TARGET>/af3_model/
# Uses the official DeepMind AlphaFold3 3.0.1 module. If the first attempt runs
# out of GPU memory, it retries once with run_script_a3_unified_mem.
#SBATCH --job-name=af3
#SBATCH --partition=dkhf3-lab-gpu
#SBATCH --account=dkhf3-lab
#SBATCH --gres=gpu:L40S:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00

set -o pipefail
shopt -s expand_aliases

TARGET=${1:?usage: sbatch run_af3_target.slurm <TARGET>}
TBM=/cluster/pixstor/dkhf3-lab/users/nancy/TBM_AGENT
OUT=$TBM/targets/$TARGET/af3_model
INPUT=$OUT/$TARGET.json
LOG=$OUT/af3_run.log
[[ -f "$INPUT" ]] || { echo "Missing $INPUT"; exit 1; }

export TMPDIR=$OUT/tmp_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"

module load alphafold3/alphafold3_v301_deepmind

ARGS=(--json_path="$INPUT"
      --model_dir=/cluster/VAST/cietest/alphafold3_model_parameters
      --db_dir=/cluster/VAST/cietest/alphafold3_database
      --output_dir="$OUT"
      --jackhmmer_n_cpu=8
      --nhmmer_n_cpu=8)

echo "[$(date)] Target: $TARGET  Host: $(hostname)  GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader)"
echo "[$(date)] Attempt 1: run_script_a3 (AlphaFold3 3.0.1)"
run_script_a3 "${ARGS[@]}" 2>&1 | tee "$LOG"
rc=$?

if [[ $rc -ne 0 ]] && grep -qiE "RESOURCE_EXHAUSTED|out of memory|OOM" "$LOG"; then
    echo "[$(date)] Out of GPU memory; attempt 2: run_script_a3_unified_mem"
    # Reuse the MSA/templates from attempt 1 if its data JSON was written
    DATA=$(ls "$OUT"/*/*_data.json 2>/dev/null | head -1)
    if [[ -n "$DATA" ]]; then
        ARGS=("${ARGS[@]/--json_path=*/--json_path=$DATA}" --norun_data_pipeline)
    fi
    run_script_a3_unified_mem "${ARGS[@]}" 2>&1 | tee -a "$LOG"
    rc=$?
fi

echo "[$(date)] AlphaFold3 returned $rc"
rm -rf "$TMPDIR"
exit $rc
