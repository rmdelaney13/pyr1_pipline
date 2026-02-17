#!/bin/bash
#SBATCH --job-name=cst_relax
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb472_asc2
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Backbone-constrained FastRelax for threaded mutant structures.
# Single job (not array) - one mutant PDB per pair.
#
# Usage:
#   sbatch submit_constrained_relax.sh <input_pdb> <output_pdb> [params_file]
#
# Arguments:
#   $1 - Input PDB (threaded mutant)
#   $2 - Output PDB (relaxed mutant)
#   $3 - Optional: extra_res_fa params file

# SLURM starts in $HOME by default - change to submit directory
cd "$SLURM_SUBMIT_DIR" || exit 1

INPUT_PDB="$1"
OUTPUT_PDB="$2"
PARAMS_FILE="$3"

if [ -z "$INPUT_PDB" ] || [ -z "$OUTPUT_PDB" ]; then
    echo "ERROR: Usage: submit_constrained_relax.sh <input_pdb> <output_pdb> [params_file]"
    exit 1
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

echo "Constrained Relax"
echo "  Input:  $INPUT_PDB"
echo "  Output: $OUTPUT_PDB"
echo "  Params: ${PARAMS_FILE:-none}"
echo "  Job ID: ${SLURM_JOB_ID:-local}"

CMD="python ml_modelling/scripts/constrained_relax.py \"$INPUT_PDB\" \"$OUTPUT_PDB\""
if [ -n "$PARAMS_FILE" ]; then
    CMD="$CMD --params \"$PARAMS_FILE\""
fi

eval $CMD

echo "Done (exit code $?)"
