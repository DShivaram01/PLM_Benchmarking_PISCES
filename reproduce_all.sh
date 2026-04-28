#!/bin/bash
set -e

echo "========================================================"
echo "  PLM Benchmarking PISCES — full reproduction pipeline"
echo "========================================================"
echo ""

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

# Use the plm conda environment's Python directly so the script works
# whether or not the environment is already activated.
PYTHON=/home/aubsxd002/.conda/envs/plm/bin/python

echo "Step 1/6: Overall correlations"
$PYTHON scripts/01_overall_correlations.py
echo ""

echo "Step 2/6: Identity-bin correlations"
$PYTHON scripts/02_identity_bin_correlations.py
echo ""

echo "Step 3/6: Length-ratio correlations"
$PYTHON scripts/03_length_ratio_correlations.py
echo ""

echo "Step 4/6: TM-threshold correlations"
$PYTHON scripts/04_tm_threshold_correlations.py
echo ""

echo "Step 5/6: BLAST missing-hit analysis"
$PYTHON scripts/05_blast_missing_analysis.py
echo ""

echo "Step 6/6: Summary leaderboard"
$PYTHON scripts/06_summary_leaderboard.py
echo ""

echo "========================================================"
echo "  All done.  Results written to: results/"
echo "========================================================"
