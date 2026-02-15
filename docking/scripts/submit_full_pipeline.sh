#!/bin/bash
#SBATCH --job-name=pipeline_orchestrator
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb472_asc2
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# ============================================================================
# Pipeline Orchestrator (Phase 1): SDF → Docking → Design → AF3 Prep
# ============================================================================
#
# Submits a lightweight SLURM orchestrator job (1 CPU, 4GB, 8hr) that chains
# docking through AF3 JSON prep. Survives SSH disconnects — submit and log out.
# The real compute happens in array jobs it spawns (MPNN, Rosetta, docking).
#
# Recommended 3-phase workflow:
#
#   Phase 1 (this script, ~4 hours, submit-and-forget):
#     bash submit_full_pipeline.sh config.txt \
#       --design-args "--skip-af3-submit --skip-af3-analyze"
#
#   Phase 2 (seconds, from login node, after Phase 1 finishes):
#     python run_design_pipeline.py config.txt --af3-submit-only
#
#   Phase 3 (from login node, after AF3 GPU jobs finish):
#     python run_design_pipeline.py config.txt --af3-analyze-only
#
# Options:
#   --design-only              Skip docking, run design pipeline only
#   --docking-only             Run docking only, skip design pipeline
#   --design-args "ARGS"       Extra args for design pipeline
#
# Monitor progress:
#   tail -f /scratch/.../logs/pipeline_orchestrator_<JOBID>.out
#   cat /scratch/.../logs/pipeline_status.log
#
# ============================================================================

set -euo pipefail

# ── Parse arguments ──────────────────────────────────────────────────────────
# When called directly (not via sbatch), this script submits itself.
# When running as the orchestrator, _PIPELINE_ORCHESTRATOR_RUNNING is set.

CONFIG_FILE=""
DESIGN_ONLY=false
DOCKING_ONLY=false
EXTRA_DESIGN_ARGS=""

# Parse all arguments (works both in submission wrapper and inside SLURM)
args=("$@")
positional=()
for ((i=0; i<${#args[@]}; i++)); do
    case "${args[$i]}" in
        --design-only)  DESIGN_ONLY=true ;;
        --docking-only) DOCKING_ONLY=true ;;
        --design-args)
            i=$((i+1))
            EXTRA_DESIGN_ARGS="${args[$i]}"
            ;;
        *)  positional+=("${args[$i]}") ;;
    esac
done

if [ ${#positional[@]} -lt 1 ]; then
    echo "ERROR: No config file provided"
    echo ""
    echo "Usage: bash submit_full_pipeline.sh config.txt [options]"
    echo ""
    echo "Options:"
    echo "  --design-only              Skip docking, run design pipeline only"
    echo "  --docking-only             Run docking only, skip design pipeline"
    echo "  --design-args \"ARGS\"       Extra args for design pipeline (e.g. \"--rosetta-to-af3\")"
    echo ""
    echo "Examples:"
    echo "  bash submit_full_pipeline.sh config.txt"
    echo "  bash submit_full_pipeline.sh config.txt --design-only"
    echo "  bash submit_full_pipeline.sh config.txt --design-args \"--af3-analyze-only\""
    exit 1
fi

CONFIG_FILE="${positional[0]}"

# Convert config to absolute path
if [[ "$CONFIG_FILE" != /* ]]; then
    CONFIG_FILE="$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")"
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Determine pipeline root (this script lives in docking/scripts/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Read SCRATCH_ROOT from config for log directory
SCRATCH_ROOT=$(grep -E "^\s*SCRATCH_ROOT\s*=" "$CONFIG_FILE" | head -1 | sed 's/.*=\s*//' | sed 's/\s*#.*//' | tr -d ' ')
if [ -z "$SCRATCH_ROOT" ]; then
    echo "ERROR: Could not read SCRATCH_ROOT from config file"
    exit 1
fi
LOG_DIR="${SCRATCH_ROOT}/logs"
STATUS_LOG="${LOG_DIR}/pipeline_status.log"

# Read CAMPAIGN_ROOT for working directory
CAMPAIGN_ROOT=$(grep -E "^\s*CAMPAIGN_ROOT\s*=" "$CONFIG_FILE" | head -1 | sed 's/.*=\s*//' | sed 's/\s*#.*//' | tr -d ' ')
if [ -z "$CAMPAIGN_ROOT" ]; then
    CAMPAIGN_ROOT="$(dirname "$CONFIG_FILE")"
fi

# ── Self-submit if not already the orchestrator job ──────────────────────────
# We use a custom env var instead of SLURM_JOB_ID because the user may
# already be inside an interactive SLURM session (sinteractive / salloc).
if [ -z "${_PIPELINE_ORCHESTRATOR_RUNNING:-}" ]; then
    mkdir -p "$LOG_DIR"

    echo "============================================================================"
    echo "Submitting Full Pipeline Orchestrator"
    echo "============================================================================"
    echo "Config:     $CONFIG_FILE"
    echo "Pipe root:  $PIPE_ROOT"
    echo "Log dir:    $LOG_DIR"
    echo "Mode:       $(if $DESIGN_ONLY; then echo 'design-only'; elif $DOCKING_ONLY; then echo 'docking-only'; else echo 'full (docking + design)'; fi)"
    if [ -n "$EXTRA_DESIGN_ARGS" ]; then
        echo "Design args: $EXTRA_DESIGN_ARGS"
    fi
    echo ""

    # Build the sbatch command, forwarding all original arguments
    JOB_ID=$(sbatch --parsable \
        --output="${LOG_DIR}/pipeline_orchestrator_%j.out" \
        --error="${LOG_DIR}/pipeline_orchestrator_%j.err" \
        --export="ALL,PIPELINE_SCRIPT_DIR=${SCRIPT_DIR},_PIPELINE_ORCHESTRATOR_RUNNING=1" \
        "$SCRIPT_DIR/submit_full_pipeline.sh" "$@")

    echo "Orchestrator job submitted: $JOB_ID"
    echo ""
    echo "Monitor progress:"
    echo "  tail -f ${LOG_DIR}/pipeline_orchestrator_${JOB_ID}.out"
    echo "  cat ${LOG_DIR}/pipeline_status.log"
    echo ""
    echo "Check queue:"
    echo "  squeue -u \$USER"
    echo ""
    echo "Cancel everything:"
    echo "  scancel $JOB_ID    # orchestrator (child jobs will also stop)"
    echo "============================================================================"
    exit 0
fi

# ── Running inside SLURM — execute the pipeline ─────────────────────────────
mkdir -p "$LOG_DIR"

# Status logging helper
log_status() {
    local stage="$1"
    local status="$2"
    local detail="${3:-}"
    local timestamp
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
    local line="[$timestamp] $stage: $status"
    if [ -n "$detail" ]; then
        line="$line — $detail"
    fi
    echo "$line" | tee -a "$STATUS_LOG"
}

# Write header to status log
cat > "$STATUS_LOG" <<HEADER
============================================================================
Pipeline Status Log
============================================================================
Job ID:      $SLURM_JOB_ID
Node:        $(hostname)
Config:      $CONFIG_FILE
Started:     $(date)
Mode:        $(if $DESIGN_ONLY; then echo 'design-only'; elif $DOCKING_ONLY; then echo 'docking-only'; else echo 'full (docking + design)'; fi)
============================================================================

HEADER

echo "============================================================================"
echo "Pipeline Orchestrator Running"
echo "============================================================================"
echo "Job ID:      $SLURM_JOB_ID"
echo "Node:        $(hostname)"
echo "Config:      $CONFIG_FILE"
echo "Pipe root:   $PIPE_ROOT"
echo "Status log:  $STATUS_LOG"
echo "Started:     $(date)"
echo "============================================================================"
echo ""

cd "$CAMPAIGN_ROOT"

# Determine which Python to use
if [ -n "${CONDA_PREFIX:-}" ]; then
    PYTHON_CMD="python"
    log_status "SETUP" "Using conda Python" "$CONDA_PREFIX"
elif [ -n "${VIRTUAL_ENV:-}" ]; then
    PYTHON_CMD="python"
    log_status "SETUP" "Using venv Python" "$VIRTUAL_ENV"
else
    PYTHON_CMD="python3"
    log_status "SETUP" "Using system Python" "$(which $PYTHON_CMD)"
fi

DOCKING_SCRIPT="$PIPE_ROOT/docking/scripts/run_docking_workflow.py"
DESIGN_SCRIPT="$PIPE_ROOT/design/scripts/run_design_pipeline.py"

EXIT_CODE=0

# ── Stage 1: Docking ────────────────────────────────────────────────────────
if ! $DESIGN_ONLY; then
    log_status "DOCKING" "STARTED" "submitting array jobs with --slurm --wait"
    if $PYTHON_CMD "$DOCKING_SCRIPT" "$CONFIG_FILE" --slurm --wait; then
        log_status "DOCKING" "COMPLETED"
    else
        log_status "DOCKING" "FAILED" "exit code $?"
        log_status "PIPELINE" "ABORTED" "docking failed, design not started"
        exit 1
    fi
fi

# ── Stage 2: Design ─────────────────────────────────────────────────────────
if ! $DOCKING_ONLY; then
    log_status "DESIGN" "STARTED" "running with --wait $EXTRA_DESIGN_ARGS"

    # Build design command
    DESIGN_CMD="$PYTHON_CMD $DESIGN_SCRIPT $CONFIG_FILE --wait"
    if [ -n "$EXTRA_DESIGN_ARGS" ]; then
        DESIGN_CMD="$DESIGN_CMD $EXTRA_DESIGN_ARGS"
    fi

    if eval $DESIGN_CMD; then
        log_status "DESIGN" "COMPLETED"
    else
        log_status "DESIGN" "FAILED" "exit code $?"
        EXIT_CODE=1
    fi
fi

# ── Done ─────────────────────────────────────────────────────────────────────
echo ""
echo "============================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    log_status "PIPELINE" "COMPLETED" "all stages finished successfully"
    echo "PIPELINE COMPLETE"
else
    log_status "PIPELINE" "FINISHED WITH ERRORS" "check logs above"
    echo "PIPELINE FINISHED WITH ERRORS"
fi
echo "Ended: $(date)"
echo "============================================================================"

exit $EXIT_CODE
