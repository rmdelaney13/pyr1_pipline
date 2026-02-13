#!/bin/bash
#SBATCH --job-name=cluster_docks
#SBATCH --output=clustering_%j.out
#SBATCH --error=clustering_%j.err
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb472_asc2
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# ============================================================================
# Clustering-Only Script
# ============================================================================
#
# Run this after all array tasks complete to cluster the final results.
#
# Usage:
#   sbatch run_clustering_only.sh config.txt
#   # or directly:
#   bash run_clustering_only.sh config.txt
# ============================================================================

if [ -z "$1" ]; then
    echo "ERROR: No config file provided"
    echo "Usage: bash run_clustering_only.sh config.txt"
    exit 1
fi

CONFIG_FILE="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_SCRIPT="${SCRIPT_DIR}/run_docking_workflow.py"

# Determine which Python to use
if [ -n "$CONDA_PREFIX" ]; then
    PYTHON_CMD="python"
elif [ -n "$VIRTUAL_ENV" ]; then
    PYTHON_CMD="python"
else
    PYTHON_CMD="python3"
fi

echo "=========================================="
echo "Running Clustering Step"
echo "=========================================="
echo "Config: $CONFIG_FILE"
echo "Start time: $(date)"
echo "=========================================="

$PYTHON_CMD "$WORKFLOW_SCRIPT" "$CONFIG_FILE" \
    --skip-create-table \
    --skip-docking

EXIT_CODE=$?

echo "=========================================="
echo "End time: $(date)"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Clustering completed successfully"
else
    echo "Clustering FAILED with exit code $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE
