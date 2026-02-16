#!/bin/bash
# Quick test of batch processing workflow
# Tests with 3 pairs from existing data

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$SCRIPT_DIR/../.."

echo "═══════════════════════════════════════════════════════════"
echo "PYR1 ML Batch Processing - Quick Test"
echo "═══════════════════════════════════════════════════════════"
echo ""

# 1. Prepare test dataset (first 3 pairs from existing data)
echo "[1/3] Preparing test pairs CSV (3 pairs)..."
python "$SCRIPT_DIR/prepare_pairs_dataset.py" \
    --existing-csv "$ROOT_DIR/ml_modelling/data/ligand_smiles_signature.csv" \
    --output "$ROOT_DIR/ml_modelling/test_pairs.csv" \
    --label 1 \
    --label-tier positive \
    --label-source validated_binder

# Limit to 3 pairs for quick test
head -4 "$ROOT_DIR/ml_modelling/test_pairs.csv" > "$ROOT_DIR/ml_modelling/test_3_pairs.csv"

echo "✓ Test pairs CSV created: ml_modelling/test_3_pairs.csv"
echo ""

# 2. Show what will be processed
echo "[2/3] Test pairs summary:"
cat "$ROOT_DIR/ml_modelling/test_3_pairs.csv"
echo ""

# 3. Show next command to run orchestrator
echo "[3/3] Ready to run orchestrator!"
echo ""
echo "To process these 3 pairs, run:"
echo ""
echo "python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \\"
echo "    --pairs-csv ml_modelling/test_3_pairs.csv \\"
echo "    --cache-dir /scratch/test_ml_cache \\"
echo "    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \\"
echo "    --docking-repeats 10 \\"
echo "    --max-pairs 3"
echo ""
echo "This will:"
echo "  • Generate conformers for 3 ligands"
echo "  • Thread 3 variant structures"
echo "  • Dock with 10 repeats (fast test)"
echo "  • Cluster and analyze results"
echo ""
echo "Expected runtime: ~10-15 minutes (local, no SLURM)"
echo ""
echo "═══════════════════════════════════════════════════════════"
