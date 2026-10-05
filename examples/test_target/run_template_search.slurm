#!/bin/bash
#SBATCH --job-name=template_search
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/pixstor/dkhf3-lab/users/nancy/TBM_AGENT/examples/test_target/logs/slurm-%j.out

set -euo pipefail

EX=/cluster/pixstor/dkhf3-lab/users/nancy/TBM_AGENT/examples/test_target
OUT=$EX/hhsearch_output
UNIREF=/cluster/VAST/cietest/rosettafold_databases/UniRef30_2020_06/UniRef30_2020_06
PDB70=/cluster/VAST/cietest/alphafold_database/pdb70/pdb70
CPUS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUT"

# Load the RoseTTAFold conda env
module load rosettafold/v1.1.0
# conda activate hooks reference unset vars, so relax nounset around them
set +u
source /cluster/software/conda-envs-global/rosetta_1_2025/miniconda3/etc/profile.d/conda.sh
conda activate RoseTTAFold
set -u

echo "[$(date)] Host: $(hostname)  CPUs: $CPUS"
echo "hhblits:  $(which hhblits)"
echo "hhsearch: $(which hhsearch)"

# Step 1: build MSA with hhblits against UniRef30
echo "[$(date)] Step 1: hhblits"
hhblits -i "$EX/test_target.fasta" \
        -d "$UNIREF" \
        -oa3m "$OUT/test_target.a3m" \
        -o /dev/null \
        -n 3 -e 0.001 -cpu "$CPUS"

# Step 2: search templates with hhsearch against PDB70
echo "[$(date)] Step 2: hhsearch"
hhsearch -i "$OUT/test_target.a3m" \
         -d "$PDB70" \
         -o "$OUT/test_target.hhr" \
         -cpu "$CPUS"

echo "[$(date)] Done."
ls -lh "$OUT"
