#!/bin/bash
#SBATCH --job-name=ligand_docking
#SBATCH --output=output_%A_%a.out
#SBATCH --error=error_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --qos=normal
#SBATCH --account=ucb-general
#SBATCH --partition=amilan

CONFIG_PATH="${1:-}"
if [[ -z "$CONFIG_PATH" ]]; then
  echo "Usage: sbatch --array=1-900 $0 /path/to/config_multiple.txt"
  exit 1
fi

HOME_DIR="$(cd "$(dirname "$CONFIG_PATH")" && pwd)"
mkdir -p "${HOME_DIR}/logs"
cd "${HOME_DIR}"

PIPE_ROOT="/projects/ryde3462/software/pyr1_pipeline"
PYTHON_BIN="/projects/ryde3462/software/projects/ryde3462/conda/envs/ligand_alignment/bin/python"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Python executable not found at: ${PYTHON_BIN}"
  echo "Update PYTHON_BIN in ${0} to your cluster env path."
  exit 1
fi

"${PYTHON_BIN}" "${PIPE_ROOT}/docking/scripts/run_docking_from_sdf.py" \
  "${CONFIG_PATH}" --mode sequence --array-index "${SLURM_ARRAY_TASK_ID}"
