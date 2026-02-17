#!/bin/bash
#SBATCH --job-name=ml_relax
#SBATCH --output=ml_relax_%A_%a.out
#SBATCH --error=ml_relax_%A_%a.err
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb472_asc2
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# SLURM array wrapper for ML pipeline relax stage.
# Reads a TSV manifest file and runs relax_general_universal.py
# for the row corresponding to SLURM_ARRAY_TASK_ID.
#
# Usage: sbatch --array=1-N submit_relax_ml.sh /path/to/relax_manifest.tsv
#
# Manifest format (tab-separated, one row per structure):
#   input_pdb  output_pdb  ligand_params  xml_path  ligand_chain  water_chain

# Change to the directory where sbatch was called (project root)
# This ensures relative paths in the Python script resolve correctly
cd "$SLURM_SUBMIT_DIR" || exit 1

MANIFEST="$1"

if [ -z "$MANIFEST" ] || [ ! -f "$MANIFEST" ]; then
    echo "ERROR: Manifest file not found: $MANIFEST"
    exit 1
fi

# Read the line for this array task (1-indexed)
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")

if [ -z "$LINE" ]; then
    echo "ERROR: No entry for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r INPUT_PDB OUTPUT_PDB LIGAND_PARAMS XML_PATH LIG_CHAIN WAT_CHAIN <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: Relaxing ${INPUT_PDB}"
echo "  Output: ${OUTPUT_PDB}"
echo "  Params: ${LIGAND_PARAMS}"

python design/rosetta/relax_general_universal.py \
    "$INPUT_PDB" "$OUTPUT_PDB" "$LIGAND_PARAMS" \
    --xml_path "$XML_PATH" \
    --ligand_chain "$LIG_CHAIN" \
    --water_chain "$WAT_CHAIN"

echo "Task ${SLURM_ARRAY_TASK_ID}: Done (exit code $?)"
