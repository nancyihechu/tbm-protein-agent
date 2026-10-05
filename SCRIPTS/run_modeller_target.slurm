#!/bin/bash
# MODELLER for one target:
#   sbatch run_modeller_target.slurm <TARGET> <START> <END> <TEMPLATE> [TEMPLATE ...]
# Uses targets/<TARGET>/alignment/<TARGET>.pir and targets/<TARGET>/templates/,
# writes 5 models (DOPE + GA341) to targets/<TARGET>/tbm_models/.
#SBATCH --job-name=modeller
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00

set -euo pipefail
TARGET=$1; START=$2; END=$3; shift 3
TBM=/cluster/pixstor/dkhf3-lab/users/nancy/TBM_AGENT
D=$TBM/targets/$TARGET

module load miniconda3/v26.7.1-1_py314
set +u
eval "$(conda shell.bash hook)"
conda activate /cluster/pixstor/dkhf3-lab/users/nancy/envs/tbm_tools
set -u

echo "[$(date)] $TARGET residues $START-$END templates: $*"
python $TBM/SCRIPTS/run_modeller.py \
    --alignment "$D/alignment/$TARGET.pir" \
    --template "$@" \
    --template-dir "$D/templates" \
    --query "$TARGET" --start "$START" --end "$END" \
    --n-models 5 \
    --out-dir "$D/tbm_models" > "$D/tbm_models.log" 2>&1
tail -12 "$D/tbm_models.log"
