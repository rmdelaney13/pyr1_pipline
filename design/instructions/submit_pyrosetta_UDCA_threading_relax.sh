#!/bin/bash
#SBATCH --job-name=Thread_Relax
#SBATCH --output=thread_relax_%A_%a.out
#SBATCH --error=thread_relax_%A_%a.err
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb472_asc2
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-912   # <--- UPDATE THIS NUMBER to match the count from the 'find' command below

# === Configuration ===

# 1. WHERE ARE THE ORIGINAL INPUT PDBS? (The Parents)
# Based on your logs, they are likely here:
TEMPLATE_DIR="/scratch/alpine/ryde3462/UDCA_design_normal/filtered_1"

# 2. WHERE ARE THE MPNN OUTPUTS? (The folder containing arrayXXX_mpnn folders)
MPNN_OUTPUT_BASE="/scratch/alpine/ryde3462/UDCA_design_normal/mpnn_output_2"

# 3. WHERE DO YOU WANT THE RESULTS?
OUTPUT_DIR="/scratch/alpine/ryde3462/UDCA_design_normal/relax_2"
mkdir -p "$OUTPUT_DIR"

# 4. PATHS TO SCRIPTS/PARAMS
PYTHON_SCRIPT="/projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/UDCA_relax_threaded.py"
LIGAND_PARAMS="/projects/ryde3462/bile_acids/UDCA/conformers_deprot/0/0.params"

# === Processing Logic ===

# Find all FASTA files in the mpnn_output_1 folder
# This looks recursively for any file ending in .fa
mapfile -t FASTA_FILES < <(find "$MPNN_OUTPUT_BASE" -type f -name "*.fa" | sort)
TOTAL=${#FASTA_FILES[@]}

# Select the specific file for this Slurm Array ID
INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [ $INDEX -ge $TOTAL ]; then
    echo "Task ID $SLURM_ARRAY_TASK_ID exceeds total files ($TOTAL). Exiting."
    exit 0
fi

FASTA_FILE="${FASTA_FILES[$INDEX]}"
BASENAME=$(basename "$FASTA_FILE" .fa)

# Locate the matching Parent PDB in the Template Directory
TEMPLATE_PDB="${TEMPLATE_DIR}/${BASENAME}.pdb"

if [[ ! -f "$TEMPLATE_PDB" ]]; then
    echo "CRITICAL ERROR: Parent PDB not found!"
    echo "Looking for: $TEMPLATE_PDB"
    echo "Please check if TEMPLATE_DIR is correct."
    exit 1
fi

echo "------------------------------------------------"
echo "Task ID:    $SLURM_ARRAY_TASK_ID"
echo "Design:     $BASENAME"
echo "Template:   $TEMPLATE_PDB"
echo "Sequences:  $FASTA_FILE"
echo "------------------------------------------------"

# Run the Python script
python "$PYTHON_SCRIPT" "$TEMPLATE_PDB" "$FASTA_FILE" "$LIGAND_PARAMS" "$OUTPUT_DIR/${BASENAME}"

echo "Done."
