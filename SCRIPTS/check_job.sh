#!/bin/bash
# Show your queued/running jobs and tail the most recent SLURM log in the repo.
# Usage: check_job.sh [N_LINES]
TBM=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
N=${1:-30}

squeue -u "$USER"
echo
LATEST=$(find "$TBM" -name "slurm-*.out" -printf "%T@ %p\n" 2>/dev/null | sort -n | tail -1 | cut -d" " -f2-)
if [[ -n "$LATEST" ]]; then
    echo "=== $LATEST (last $N lines) ==="
    tail -n "$N" "$LATEST"
else
    echo "No slurm-*.out files under $TBM yet."
fi
