#!/bin/bash
# ============================================================================
# Run Rosetta Relax ONLY (Standalone for Troubleshooting)
# ============================================================================
#
# This script runs ONLY the Rosetta relax step on MPNN output.
# Use this for troubleshooting or when you need to re-run just Rosetta.
#
# Usage:
#   bash run_rosetta_only.sh <template_pdb_dir> <mpnn_output_dir> <output_dir> <ligand_params> [array_count]
#
# Example:
#   bash run_rosetta_only.sh ./clustered_final ./mpnn_output ./rosetta_output /path/to/ligand.params 500
#
# Arguments:
#   template_pdb_dir  - Directory with original/parent PDB files
#   mpnn_output_dir   - Directory with MPNN output (contains .fa files)
#   output_dir        - Where to save Rosetta output
#   ligand_params     - Path to ligand .params file
#   array_count       - (Optional) Number of parallel array tasks (default: auto-detect from FASTA count)
#
# ============================================================================

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    echo "ERROR: Missing required arguments"
    echo ""
    echo "Usage: bash run_rosetta_only.sh <template_pdb_dir> <mpnn_output_dir> <output_dir> <ligand_params> [array_count]"
    echo ""
    echo "Example:"
    echo "  bash run_rosetta_only.sh ./clustered_final ./mpnn_output ./rosetta_output /path/to/ligand.params 500"
    exit 1
fi

TEMPLATE_DIR="$(realpath "$1")"
MPNN_OUTPUT_DIR="$(realpath "$2")"
OUTPUT_DIR="$(realpath "$3")"
LIGAND_PARAMS="$(realpath "$4")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Validate directories and files
if [ ! -d "$TEMPLATE_DIR" ]; then
    echo "ERROR: Template PDB directory not found: $TEMPLATE_DIR"
    exit 1
fi

if [ ! -d "$MPNN_OUTPUT_DIR" ]; then
    echo "ERROR: MPNN output directory not found: $MPNN_OUTPUT_DIR"
    exit 1
fi

if [ ! -f "$LIGAND_PARAMS" ]; then
    echo "ERROR: Ligand params file not found: $LIGAND_PARAMS"
    exit 1
fi

# Count FASTA files
FASTA_COUNT=$(find "$MPNN_OUTPUT_DIR" -type f -name "*.fa" | wc -l)
if [ "$FASTA_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTA files found in $MPNN_OUTPUT_DIR"
    exit 1
fi

# Determine array count
if [ -n "$5" ]; then
    ARRAY_COUNT="$5"
else
    ARRAY_COUNT="$FASTA_COUNT"
fi

ARRAY_MAX=$((ARRAY_COUNT - 1))

echo "============================================================================"
echo "Running Rosetta Relax (Standalone)"
echo "============================================================================"
echo "Template PDBs: $TEMPLATE_DIR"
echo "MPNN output: $MPNN_OUTPUT_DIR"
echo "  Found $FASTA_COUNT FASTA files"
echo "Ligand params: $LIGAND_PARAMS"
echo "Output dir: $OUTPUT_DIR"
echo "Array tasks: 0-$ARRAY_MAX (total: $ARRAY_COUNT)"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Path to Rosetta template script
ROSETTA_TEMPLATE="$PIPELINE_ROOT/design/instructions/submit_pyrosetta_general_threading_relax.sh"

if [ ! -f "$ROSETTA_TEMPLATE" ]; then
    echo "ERROR: Rosetta script template not found at: $ROSETTA_TEMPLATE"
    exit 1
fi

# Create customized Rosetta script
ROSETTA_CUSTOM="$OUTPUT_DIR/submit_rosetta_custom.sh"
cp "$ROSETTA_TEMPLATE" "$ROSETTA_CUSTOM"

# Update paths in the script
sed -i "s|^TEMPLATE_DIR=.*|TEMPLATE_DIR=\"$TEMPLATE_DIR\"|" "$ROSETTA_CUSTOM"
sed -i "s|^MPNN_OUTPUT_BASE=.*|MPNN_OUTPUT_BASE=\"$MPNN_OUTPUT_DIR\"|" "$ROSETTA_CUSTOM"
sed -i "s|^OUTPUT_DIR=.*|OUTPUT_DIR=\"$OUTPUT_DIR\"|" "$ROSETTA_CUSTOM"
sed -i "s|^LIGAND_PARAMS=.*|LIGAND_PARAMS=\"$LIGAND_PARAMS\"|" "$ROSETTA_CUSTOM"
sed -i "s|#SBATCH --array=1-.*|#SBATCH --array=1-$ARRAY_COUNT|" "$ROSETTA_CUSTOM"

echo "✓ Created custom Rosetta script: $ROSETTA_CUSTOM"
echo ""
echo "Submitting to SLURM..."

# Submit to SLURM
JOB_ID=$(sbatch --parsable "$ROSETTA_CUSTOM")

if [ -z "$JOB_ID" ]; then
    echo "ERROR: Failed to submit Rosetta job"
    exit 1
fi

echo "✓ Rosetta job submitted: $JOB_ID"
echo ""
echo "Monitor progress:"
echo "  squeue -j $JOB_ID"
echo "  tail -f thread_relax_${JOB_ID}_*.out"
echo ""
echo "Cancel job:"
echo "  scancel $JOB_ID"
echo ""
echo "Output will be saved to: $OUTPUT_DIR"
echo "============================================================================"
