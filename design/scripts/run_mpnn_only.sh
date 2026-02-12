#!/bin/bash
# ============================================================================
# Run LigandMPNN Design ONLY (Standalone for Troubleshooting)
# ============================================================================
#
# This script runs ONLY the LigandMPNN design step on a directory of PDBs.
# Use this for troubleshooting or when you need to re-run just MPNN.
#
# Usage:
#   bash run_mpnn_only.sh <input_pdb_dir> <output_dir> [array_count]
#
# Example:
#   bash run_mpnn_only.sh /path/to/docked/pdbs /path/to/mpnn/output 50
#
# Arguments:
#   input_pdb_dir  - Directory containing .pdb files to design
#   output_dir     - Where to save MPNN output
#   array_count    - (Optional) Number of parallel array tasks (default: auto-detect from PDB count)
#
# ============================================================================

if [ -z "$1" ] || [ -z "$2" ]; then
    echo "ERROR: Missing required arguments"
    echo ""
    echo "Usage: bash run_mpnn_only.sh <input_pdb_dir> <output_dir> [array_count]"
    echo ""
    echo "Example:"
    echo "  bash run_mpnn_only.sh ./clustered_final ./mpnn_output 50"
    exit 1
fi

INPUT_PDB_DIR="$(realpath "$1")"
OUTPUT_DIR="$(realpath "$2")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Validate input directory
if [ ! -d "$INPUT_PDB_DIR" ]; then
    echo "ERROR: Input PDB directory not found: $INPUT_PDB_DIR"
    exit 1
fi

# Count PDB files
PDB_COUNT=$(ls "$INPUT_PDB_DIR"/*.pdb 2>/dev/null | wc -l)
if [ "$PDB_COUNT" -eq 0 ]; then
    echo "ERROR: No PDB files found in $INPUT_PDB_DIR"
    exit 1
fi

# Determine array count
if [ -n "$3" ]; then
    ARRAY_COUNT="$3"
else
    ARRAY_COUNT="$PDB_COUNT"
fi

ARRAY_MAX=$((ARRAY_COUNT - 1))

echo "============================================================================"
echo "Running LigandMPNN Design (Standalone)"
echo "============================================================================"
echo "Input PDBs: $INPUT_PDB_DIR"
echo "  Found $PDB_COUNT PDB files"
echo "Output dir: $OUTPUT_DIR"
echo "Array tasks: 0-$ARRAY_MAX (total: $ARRAY_COUNT)"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Path to MPNN template script
MPNN_TEMPLATE="$PIPELINE_ROOT/design/instructions/ligand_alignment_mpnni_grouped.sh"

if [ ! -f "$MPNN_TEMPLATE" ]; then
    echo "ERROR: MPNN script template not found at: $MPNN_TEMPLATE"
    exit 1
fi

# Create customized MPNN script
MPNN_CUSTOM="$OUTPUT_DIR/submit_mpnn_custom.sh"
cp "$MPNN_TEMPLATE" "$MPNN_CUSTOM"

# Update paths in the script
sed -i "s|^PDB_DIR=.*|PDB_DIR=\"$INPUT_PDB_DIR\"|" "$MPNN_CUSTOM"
sed -i "s|^OUTPUT_BASE=.*|OUTPUT_BASE=\"$OUTPUT_DIR\"|" "$MPNN_CUSTOM"
sed -i "s|#SBATCH --array=1-.*|#SBATCH --array=1-$ARRAY_COUNT|" "$MPNN_CUSTOM"

echo "✓ Created custom MPNN script: $MPNN_CUSTOM"
echo ""
echo "Submitting to SLURM..."

# Submit to SLURM
JOB_ID=$(sbatch --parsable "$MPNN_CUSTOM")

if [ -z "$JOB_ID" ]; then
    echo "ERROR: Failed to submit MPNN job"
    exit 1
fi

echo "✓ MPNN job submitted: $JOB_ID"
echo ""
echo "Monitor progress:"
echo "  squeue -j $JOB_ID"
echo "  tail -f LigandMPNN_batch_${JOB_ID}_*.out"
echo ""
echo "Cancel job:"
echo "  scancel $JOB_ID"
echo ""
echo "Output will be saved to: $OUTPUT_DIR"
echo "============================================================================"
