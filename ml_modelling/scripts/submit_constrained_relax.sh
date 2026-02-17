#!/bin/bash
#SBATCH --job-name=repack_mut
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb472_asc2
#SBATCH --time=00:05:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Sidechain repack around mutation sites (PackRotamersMover, not FastRelax).
# Single job (not array) - one mutant PDB per pair.
#
# Usage:
#   sbatch submit_constrained_relax.sh <input_pdb> <output_pdb> <mutations> [params_file]
#
# Arguments:
#   $1 - Input PDB (threaded mutant)
#   $2 - Output PDB (relaxed mutant)
#   $3 - Mutation signature (e.g., "59K;120A;160G")
#   $4 - Optional: extra_res_fa params file

# SLURM starts in $HOME by default - change to submit directory
cd "$SLURM_SUBMIT_DIR" || exit 1

INPUT_PDB="$1"
OUTPUT_PDB="$2"
MUTATIONS="$3"
PARAMS_FILE="$4"

if [ -z "$INPUT_PDB" ] || [ -z "$OUTPUT_PDB" ] || [ -z "$MUTATIONS" ]; then
    echo "ERROR: Usage: submit_constrained_relax.sh <input_pdb> <output_pdb> <mutations> [params_file]"
    exit 1
fi

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

echo "Sidechain Repack"
echo "  Input:     $INPUT_PDB"
echo "  Output:    $OUTPUT_PDB"
echo "  Mutations: $MUTATIONS"
echo "  Params:    ${PARAMS_FILE:-none}"
echo "  Job ID:    ${SLURM_JOB_ID:-local}"

CMD="python ml_modelling/scripts/constrained_relax.py \"$INPUT_PDB\" \"$OUTPUT_PDB\" --mutations \"$MUTATIONS\""
if [ -n "$PARAMS_FILE" ]; then
    CMD="$CMD --params \"$PARAMS_FILE\""
fi

eval $CMD

echo "Done (exit code $?)"
