#!/bin/bash

# ============================================================================
# Complete Docking Workflow - Single Command Submission
# ============================================================================
#
# This script submits the entire workflow with automatic job dependencies:
#   1. Array job (parallel docking tasks)
#   2. Clustering job (runs automatically after all array tasks complete)
#
# Usage:
#   bash submit_complete_workflow.sh config.txt [array_count]
#
# Example:
#   bash submit_complete_workflow.sh config.txt 10
#   # Submits array tasks 0-9, then auto-runs clustering when done
#
# ============================================================================

if [ -z "$1" ]; then
    echo "ERROR: No config file provided"
    echo ""
    echo "Usage: bash submit_complete_workflow.sh config.txt [array_count]"
    echo ""
    echo "Arguments:"
    echo "  config.txt    - Path to your config file (required)"
    echo "  array_count   - Number of array tasks (optional, reads from config if not provided)"
    echo ""
    echo "Example:"
    echo "  bash submit_complete_workflow.sh myconfig.txt 10"
    exit 1
fi

CONFIG_FILE="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Determine array count
if [ -n "$2" ]; then
    ARRAY_COUNT="$2"
else
    # Try to read from config file
    ARRAY_COUNT=$(grep -E "^\s*ArrayTaskCount\s*=" "$CONFIG_FILE" | head -1 | sed 's/.*=\s*//' | sed 's/\s*#.*//' | tr -d ' ')
    if [ -z "$ARRAY_COUNT" ]; then
        echo "ERROR: Could not determine ArrayTaskCount from config file"
        echo "Please specify array count as second argument or set ArrayTaskCount in config"
        exit 1
    fi
fi

# Validate array count
if ! [[ "$ARRAY_COUNT" =~ ^[0-9]+$ ]] || [ "$ARRAY_COUNT" -lt 1 ]; then
    echo "ERROR: Invalid array count: $ARRAY_COUNT"
    echo "Must be a positive integer"
    exit 1
fi

ARRAY_MAX=$((ARRAY_COUNT - 1))

echo "============================================================================"
echo "Submitting Complete Docking Workflow"
echo "============================================================================"
echo "Config file: $CONFIG_FILE"
echo "Array tasks: 0-$ARRAY_MAX (total: $ARRAY_COUNT)"
echo "Script directory: $SCRIPT_DIR"
echo ""

# Submit array job
echo "Step 1: Submitting array job..."
ARRAY_JOB_ID=$(sbatch --parsable --array=0-$ARRAY_MAX "$SCRIPT_DIR/submit_docking_workflow.sh" "$CONFIG_FILE")

if [ -z "$ARRAY_JOB_ID" ]; then
    echo "ERROR: Failed to submit array job"
    exit 1
fi

echo "  ✓ Array job submitted: $ARRAY_JOB_ID"
echo "  ✓ Array tasks: $ARRAY_JOB_ID"_"0 to $ARRAY_JOB_ID"_"$ARRAY_MAX"
echo ""

# Submit clustering job with dependency on array job completion
echo "Step 2: Submitting clustering job (will start after array completes)..."
CLUSTER_JOB_ID=$(sbatch --parsable --dependency=afterok:$ARRAY_JOB_ID "$SCRIPT_DIR/run_clustering_only.sh" "$CONFIG_FILE")

if [ -z "$CLUSTER_JOB_ID" ]; then
    echo "ERROR: Failed to submit clustering job"
    echo "Array job $ARRAY_JOB_ID is still running"
    exit 1
fi

echo "  ✓ Clustering job submitted: $CLUSTER_JOB_ID"
echo "  ✓ Dependency: Will start after job $ARRAY_JOB_ID completes"
echo ""

echo "============================================================================"
echo "Workflow Submitted Successfully!"
echo "============================================================================"
echo ""
echo "Job IDs:"
echo "  Array job:      $ARRAY_JOB_ID (tasks 0-$ARRAY_MAX)"
echo "  Clustering job: $CLUSTER_JOB_ID (will run after array completes)"
echo ""
echo "Monitor progress:"
echo "  squeue -u \$USER"
echo "  squeue -j $ARRAY_JOB_ID,$CLUSTER_JOB_ID"
echo ""
echo "View outputs:"
echo "  tail -f docking_${ARRAY_JOB_ID}_*.out"
echo "  tail -f clustering_${CLUSTER_JOB_ID}.out"
echo ""
echo "Cancel jobs:"
echo "  scancel $ARRAY_JOB_ID $CLUSTER_JOB_ID"
echo ""
echo "============================================================================"
echo ""
echo "The workflow will complete automatically. No further action needed!"
echo "Check back later for results in your OutputDir/clustered_final/"
echo "============================================================================"
