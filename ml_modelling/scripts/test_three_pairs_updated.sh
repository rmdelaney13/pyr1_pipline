#!/bin/bash
#
# Test script for running 3 example pairs through the ML docking pipeline
# Examples: abscisic_Acid (wt), nitazene (PYR1^nitav2), Lithocholic_Acid (seq16_designed)
#
# Usage:
#   bash ml_modelling/scripts/test_three_pairs_updated.sh
#
# This uses the NEW batch processing workflow with prepare_pairs_dataset.py
#

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$SCRIPT_DIR/../.."

echo "=========================================="
echo "ML Pipeline Test: 3 Example Pairs"
echo "=========================================="
echo ""
echo "Pairs:"
echo "  1. abscisic_Acid + wt (wildtype)"
echo "  2. nitazene + PYR1^nitav2"
echo "  3. Lithocholic_Acid + seq16_designed"
echo ""

# ===== CONFIGURATION =====
PROJECT_ROOT="$ROOT_DIR"
CACHE_DIR="${ROOT_DIR}/ml_modelling/test_cache"  # Local cache for testing
TEMPLATE_PDB="${ROOT_DIR}/docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb"

# Input CSV with 3 test pairs (already exists)
INPUT_CSV="${ROOT_DIR}/ml_modelling/data/test_three_pairs.csv"

# Prepared CSV (output of prepare_pairs_dataset.py)
PREPARED_CSV="${ROOT_DIR}/ml_modelling/test_three_pairs_prepared.csv"

echo "Configuration:"
echo "  Project root: $PROJECT_ROOT"
echo "  Input CSV: $INPUT_CSV"
echo "  Template PDB: $TEMPLATE_PDB"
echo "  Cache dir: $CACHE_DIR"
echo ""

# Check if input CSV exists
if [ ! -f "$INPUT_CSV" ]; then
    echo "ERROR: Test CSV not found: $INPUT_CSV"
    exit 1
fi

# Check if template PDB exists
if [ ! -f "$TEMPLATE_PDB" ]; then
    echo "ERROR: Template PDB not found: $TEMPLATE_PDB"
    echo "Please check the path!"
    exit 1
fi

# Create cache directory
mkdir -p "$CACHE_DIR"

# ===== STEP 1: PREPARE PAIRS CSV =====
echo "[1/2] Preparing pairs CSV (adding pair_id and label_source)..."
python "${PROJECT_ROOT}/ml_modelling/scripts/prepare_pairs_dataset.py" \
    --pairs-csv "$INPUT_CSV" \
    --output "$PREPARED_CSV" \
    --label 1 \
    --label-tier positive \
    --label-source test_pairs

echo "✓ Prepared CSV: $PREPARED_CSV"
echo ""

# Show what we're about to process
echo "Prepared pairs:"
head -4 "$PREPARED_CSV"
echo ""

# ===== STEP 2: RUN ORCHESTRATOR =====
echo "[2/2] Running pipeline orchestration (local, no SLURM)..."
echo ""

# Run locally with reduced docking repeats for quick testing
python "${PROJECT_ROOT}/ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py" \
    --pairs-csv "$PREPARED_CSV" \
    --cache-dir "$CACHE_DIR" \
    --template-pdb "$TEMPLATE_PDB" \
    --docking-repeats 10 \
    --max-pairs 3

# For SLURM (uncomment below, comment above):
# python "${PROJECT_ROOT}/ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py" \
#     --pairs-csv "$PREPARED_CSV" \
#     --cache-dir "$CACHE_DIR" \
#     --template-pdb "$TEMPLATE_PDB" \
#     --docking-repeats 50 \
#     --max-pairs 3 \
#     --use-slurm

echo ""
echo "=========================================="
echo "Pipeline orchestration complete!"
echo "=========================================="
echo ""
echo "Check results in: $CACHE_DIR"
echo ""
echo "Expected outputs per pair:"
echo "  • {pair_id}/conformers/conformers_final.sdf"
echo "  • {pair_id}/mutant.pdb"
echo "  • {pair_id}/docking/*.pdb"
echo "  • {pair_id}/docking/hbond_geometry_summary.csv"
echo "  • {pair_id}/clustered/*.pdb"
echo "  • {pair_id}/metadata.json"
echo ""
echo "Next steps:"
echo "  1. Verify conformer generation:"
echo "     ls $CACHE_DIR/*/conformers/"
echo ""
echo "  2. Check docking outputs:"
echo "     ls $CACHE_DIR/*/docking/"
echo ""
echo "  3. View H-bond geometry:"
echo "     cat $CACHE_DIR/*/docking/hbond_geometry_summary.csv"
echo ""
echo "  4. View processing summary:"
echo "     cat $CACHE_DIR/processing_summary.csv"
echo ""
echo "  5. Check metadata for each pair:"
echo "     cat $CACHE_DIR/*/metadata.json"
echo ""

# Show summary if available
if [ -f "$CACHE_DIR/processing_summary.csv" ]; then
    echo "Processing summary:"
    cat "$CACHE_DIR/processing_summary.csv"
    echo ""
fi
