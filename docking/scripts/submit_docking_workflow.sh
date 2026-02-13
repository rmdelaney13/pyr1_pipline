#!/bin/bash
#SBATCH --job-name=docking_workflow
#SBATCH --output=docking_%A_%a.out
#SBATCH --error=docking_%A_%a.err
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb472_asc2
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# ============================================================================
# SLURM Submission Script for Docking Workflow
# ============================================================================
#
# This script orchestrates the complete docking pipeline:
#   1. create_table.py (run once by array task 0)
#   2. grade_conformers_glycine_shaved_docking_multiple_slurm.py (each array task)
#   3. cluster_docked_post_array.py (run once after all arrays complete)
#
# Usage:
#   sbatch --array=0-N submit_docking_workflow.sh config.txt
#
# Where N is (ArrayTaskCount - 1) from your config file.
# For example, if ArrayTaskCount=10, submit with --array=0-9
#
# The array tasks will run in parallel, then you can manually run clustering:
#   python run_docking_workflow.py config.txt --skip-create-table --skip-docking
#
# Or use the two-stage approach below.
# ============================================================================

# Check if config file is provided
if [ -z "$1" ]; then
    echo "ERROR: No config file provided"
    echo "Usage: sbatch --array=0-N submit_docking_workflow.sh config.txt"
    exit 1
fi

CONFIG_FILE="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_SCRIPT="${SCRIPT_DIR}/run_docking_workflow.py"

# Determine which Python to use (prefer conda environment)
if [ -n "$CONDA_PREFIX" ]; then
    PYTHON_CMD="python"
    echo "Using Python from conda environment: $CONDA_PREFIX"
elif [ -n "$VIRTUAL_ENV" ]; then
    PYTHON_CMD="python"
    echo "Using Python from virtual environment: $VIRTUAL_ENV"
else
    PYTHON_CMD="python3"
    echo "Using system Python: $(which $PYTHON_CMD)"
fi

# Get the array index (SLURM_ARRAY_TASK_ID)
ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

echo "=========================================="
echo "Docking Workflow - Array Task $ARRAY_INDEX"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Config: $CONFIG_FILE"
echo "Workflow Script: $WORKFLOW_SCRIPT"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "=========================================="

# Task 0 creates the table; all tasks run docking; no task runs clustering yet
if [ "$ARRAY_INDEX" -eq 0 ]; then
    # First array task: create table and run docking
    echo "Array task 0: Creating table and running docking"
    $PYTHON_CMD "$WORKFLOW_SCRIPT" "$CONFIG_FILE" \
        --array-index "$ARRAY_INDEX" \
        --skip-clustering
else
    # Other array tasks: skip table creation, run docking only
    echo "Array task $ARRAY_INDEX: Running docking only"
    $PYTHON_CMD "$WORKFLOW_SCRIPT" "$CONFIG_FILE" \
        --array-index "$ARRAY_INDEX" \
        --skip-create-table \
        --skip-clustering
fi

EXIT_CODE=$?

echo "=========================================="
echo "End time: $(date)"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Array task $ARRAY_INDEX completed successfully"
else
    echo "Array task $ARRAY_INDEX FAILED with exit code $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE
