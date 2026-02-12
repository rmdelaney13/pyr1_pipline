#!/bin/bash
#SBATCH --job-name=LigandMPNN_batch            # Job name
#SBATCH --output=LigandMPNN_batch_%A_%a.out      # Standard output and error log
#SBATCH --error=LigandMPNN_batch_%A_%a.err
#SBATCH --partition=amilan
#SBATCH --account=ucb472_asc2
#SBATCH --ntasks=1                             # Run a single task
#SBATCH --qos=normal
#SBATCH --mem=8G                              # Total memory per task
#SBATCH --time=00:20:00                        # Time limit hrs:min:secu
#SBATCH --array=1-34                    # Array range: 900 jobs total

# ------------------------------
# Load Required Modules
# -----------------------------
module load anaconda   # Adjust Python version as needed
# module load other_modules_if_required
conda activate ligandmpnn_env


# ------------------------------
# Define Variables
# ------------------------------
PDB_DIR="/scratch/alpine/ryde3462/CA_design_flipped/filtered_1" # Directory with .pdb files
OUTPUT_BASE="/scratch/alpine/ryde3462/CA_design_flipped/mpnn_output_2"              # Base directory for LigandMPNN outputs
MODEL_SCRIPT="/projects/ryde3462/software/LigandMPNN/run.py"     # Path to LigandMPNN run script

# ------------------------------
# Generate List of PDB Files
# ------------------------------
mapfile -t pdb_files < <(ls "${PDB_DIR}"/*.pdb | sort)
TOTAL_FILES=${#pdb_files[@]}
echo "Total number of PDB files: ${TOTAL_FILES}"

# ------------------------------
# Determine Batch Size per Job
# ------------------------------
# For example, if you want to process 7 files per job:
GROUP_SIZE=1

# Calculate the start and end indices for this job
# (Bash arrays are 0-indexed, but SLURM_ARRAY_TASK_ID starts at 1)
START_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) * GROUP_SIZE ))
END_INDEX=$(( SLURM_ARRAY_TASK_ID * GROUP_SIZE - 1 ))
if [ $END_INDEX -ge $TOTAL_FILES ]; then
    END_INDEX=$(( TOTAL_FILES - 1 ))
fi

echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing files from index ${START_INDEX} to ${END_INDEX}"

# ------------------------------
# Loop over the Selected Files
# ------------------------------
for (( i=START_INDEX; i<=END_INDEX; i++ )); do
    PDB_FILE="${pdb_files[$i]}"
    PDB_BASENAME=$(basename "${PDB_FILE}" .pdb)
    OUT_FOLDER="${OUTPUT_BASE}/${PDB_BASENAME}_mpnn"
    
    echo "Processing file: ${PDB_FILE}"
    mkdir -p "${OUT_FOLDER}"

    # ------------------------------
    # Run LigandMPNN for this file
    # ------------------------------
    cd /projects/ryde3462/software/LigandMPNN
    python "${MODEL_SCRIPT}" \
        --seed 111 \
        --model_type "ligand_mpnn" \
        --pdb_path "${PDB_FILE}" \
        --redesigned_residues "A59 A79 A81 A90 A92 A106 A108 A115 A118 A120 A139 A157 A158 A161 A162 A165" \
        --out_folder "${OUT_FOLDER}" \
        --number_of_batches 1 \
        --batch_size 20 \
        --temperature 0.3 \
        --omit_AA_per_residue /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/LCA_omit_design_59.json \
	--bias_AA_per_residue /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/LCA_bias_AA_per_residue_flipped.json \
        --pack_side_chains 0 \
        --checkpoint_ligand_mpnn "/projects/ryde3462/software/LigandMPNN/model_params/ligandmpnn_v_32_020_25.pt" \
        --pack_with_ligand_context 1 

    echo "LigandMPNN processing completed for ${PDB_FILE}. Output saved to ${OUT_FOLDER}."
done

echo "All files in batch (job ${SLURM_ARRAY_TASK_ID}) have been processed."

