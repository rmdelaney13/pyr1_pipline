# Individual Step Execution Guide

This guide shows you how to run each step of the design pipeline individually for troubleshooting.

## Overview

You have **two workflow options**:

1. **Complete Automated Workflow** - Runs everything from docking → MPNN → Rosetta → AF3 prep
2. **Individual Steps** - Run each step separately for troubleshooting

---

## Complete Automated Workflow

### Option 1: Full Pipeline (Docking + Design)

```bash
# 1. Run docking workflow
bash docking/scripts/submit_complete_workflow.sh config.txt

# 2. Wait for docking to complete, then run design
python design/scripts/run_design_pipeline.py config.txt
```

### Option 2: Design Only (Docking Already Complete)

```bash
# If docking is already done, run just the design pipeline
python design/scripts/run_design_pipeline.py config.txt
```

---

## Individual Steps for Troubleshooting

### Step 1: Docking (Individual Components)

#### 1a. Create Table Only
```bash
python docking/scripts/create_table.py config.txt
```

#### 1b. Run Docking Only (Single Array Task)
```bash
python docking/scripts/grade_conformers_glycine_shaved_docking_multiple_slurm.py config.txt 0
```

#### 1c. Cluster Results Only
```bash
python docking/scripts/cluster_docked_post_array.py config.txt
```

Or use the clustering-only SLURM script:
```bash
sbatch docking/scripts/run_clustering_only.sh config.txt
```

---

### Step 2: LigandMPNN Design (Individual)

#### Quick Method (Wrapper Script)
```bash
bash design/scripts/run_mpnn_only.sh <input_pdb_dir> <output_dir> [array_count]
```

**Example:**
```bash
bash design/scripts/run_mpnn_only.sh \
    /scratch/alpine/user/output/clustered_final \
    /scratch/alpine/user/design/mpnn_output \
    50
```

#### Manual Method (Direct Script)
```bash
# 1. Edit the MPNN script with your paths
nano design/instructions/ligand_alignment_mpnni_grouped.sh

# 2. Update these variables:
#    PDB_DIR="/path/to/your/pdbs"
#    OUTPUT_BASE="/path/to/output"

# 3. Submit to SLURM
sbatch --array=1-50 design/instructions/ligand_alignment_mpnni_grouped.sh
```

#### Via Python Orchestrator (Skip Later Steps)
```bash
python design/scripts/run_design_pipeline.py config.txt \
    --skip-rosetta \
    --skip-aggregate \
    --skip-filter \
    --skip-af3-prep
```

---

### Step 3: Rosetta Relax (Individual)

#### Quick Method (Wrapper Script)
```bash
bash design/scripts/run_rosetta_only.sh \
    <template_pdb_dir> \
    <mpnn_output_dir> \
    <output_dir> \
    <ligand_params> \
    [array_count]
```

**Example:**
```bash
bash design/scripts/run_rosetta_only.sh \
    /scratch/alpine/user/output/clustered_final \
    /scratch/alpine/user/design/mpnn_output \
    /scratch/alpine/user/design/rosetta_output \
    /projects/user/ligand.params \
    500
```

#### Manual Method (Direct Script)
```bash
# 1. Edit the Rosetta script with your paths
nano design/instructions/submit_pyrosetta_general_threading_relax.sh

# 2. Update these variables:
#    TEMPLATE_DIR="/path/to/parent/pdbs"
#    MPNN_OUTPUT_BASE="/path/to/mpnn/output"
#    OUTPUT_DIR="/path/to/rosetta/output"
#    LIGAND_PARAMS="/path/to/ligand.params"

# 3. Submit to SLURM
sbatch --array=1-500 design/instructions/submit_pyrosetta_general_threading_relax.sh
```

#### Via Python Orchestrator (Skip Earlier Steps)
```bash
python design/scripts/run_design_pipeline.py config.txt \
    --skip-mpnn \
    --skip-af3-prep
```

---

### Step 4: Score Aggregation (Individual)

```bash
python design/instructions/aggregate_scores.py \
    <rosetta_output_dir> \
    --output scores.csv
```

**Example:**
```bash
python design/instructions/aggregate_scores.py \
    /scratch/alpine/user/design/rosetta_output \
    --output /scratch/alpine/user/design/scores/all_scores.csv
```

---

### Step 5: Filtering (Individual)

```bash
python design/instructions/relax_2_filter__allpolar_unsats.py \
    <scores_csv> \
    <rosetta_output_dir> \
    <filtered_output_dir> \
    --target_n 1000 \
    --max_unsat 1 \
    --max_per_parent 20
```

**Example:**
```bash
python design/instructions/relax_2_filter__allpolar_unsats.py \
    /scratch/alpine/user/design/scores/all_scores.csv \
    /scratch/alpine/user/design/rosetta_output \
    /scratch/alpine/user/design/filtered \
    --target_n 1000 \
    --max_unsat 1 \
    --max_per_parent 20
```

---

### Step 6: AF3 Preparation (Individual)

#### Only Run AF3 Prep (Skip Everything Else)
```bash
python design/scripts/run_design_pipeline.py config.txt --af3-prep-only
```

This will:
1. Generate FASTA from filtered PDBs
2. Extract SMILES from SDF
3. Create AF3 JSON inputs (binary and ternary)

---

## Common Troubleshooting Scenarios

### Scenario 1: Docking Failed, Need to Re-run Just One Array Task
```bash
# Re-run array index 5
python docking/scripts/grade_conformers_glycine_shaved_docking_multiple_slurm.py config.txt 5
```

### Scenario 2: MPNN Completed, But Rosetta Failed
```bash
# Just run Rosetta (using completed MPNN output)
bash design/scripts/run_rosetta_only.sh \
    /scratch/alpine/user/clustered_final \
    /scratch/alpine/user/mpnn_output \
    /scratch/alpine/user/rosetta_output_retry \
    /projects/user/ligand.params
```

### Scenario 3: Want to Test Different Filter Settings
```bash
# Re-run filtering with different parameters
python design/instructions/relax_2_filter__allpolar_unsats.py \
    scores.csv \
    rosetta_output \
    filtered_strict \
    --target_n 500 \
    --max_unsat 0 \
    --max_per_parent 10
```

### Scenario 4: Need to Re-generate AF3 Inputs with Different SMILES
```bash
# 1. Edit the AF3 template JSON manually with new SMILES
nano design/templates/pyr1_binary_template.json

# 2. Re-run AF3 prep only
python design/scripts/run_design_pipeline.py config.txt --af3-prep-only
```

---

## Quick Reference Table

| Step | Wrapper Script | Direct Script | Python Orchestrator |
|------|---------------|---------------|---------------------|
| **Docking: Create Table** | - | `create_table.py` | `run_docking_workflow.py --prepare-only` |
| **Docking: Run** | - | `grade_conformers_glycine_shaved_docking_multiple_slurm.py` | `run_docking_workflow.py --skip-clustering` |
| **Docking: Cluster** | - | `cluster_docked_post_array.py` | `run_docking_workflow.py --skip-create-table --skip-docking` |
| **MPNN Design** | `run_mpnn_only.sh` | `ligand_alignment_mpnni_grouped.sh` | `run_design_pipeline.py --skip-rosetta --skip-aggregate --skip-filter --skip-af3-prep` |
| **Rosetta Relax** | `run_rosetta_only.sh` | `submit_pyrosetta_general_threading_relax.sh` | `run_design_pipeline.py --skip-mpnn` |
| **Score Aggregation** | - | `aggregate_scores.py` | `run_design_pipeline.py --skip-mpnn --skip-rosetta` |
| **Filtering** | - | `relax_2_filter__allpolar_unsats.py` | `run_design_pipeline.py --skip-mpnn --skip-rosetta --skip-aggregate` |
| **AF3 Prep** | - | Manual (extract_smiles.py + make_af3_jsons.py) | `run_design_pipeline.py --af3-prep-only` |

---

## Directory Structure

After running the full pipeline, you'll have:

```
your_output/
├── docking_output/           # Docking results
│   ├── docked_*.pdb
│   └── clustered_final/      # Input to MPNN
│       └── *.pdb
├── design/
│   ├── iteration_1/
│   │   ├── mpnn_output/      # LigandMPNN designs
│   │   │   └── *_mpnn/
│   │   │       └── *.fa
│   │   ├── rosetta_output/   # Rosetta relax results
│   │   │   └── *.pdb
│   │   ├── scores/           # Aggregated scores
│   │   │   └── iteration_1_scores.csv
│   │   └── filtered/         # Filtered designs
│   │       └── *.pdb
│   └── af3_inputs/           # AF3 preparation
│       ├── binary/
│       │   └── *.json
│       └── ternary/
│           └── *.json
```

---

## Tips

1. **Always check job status** before running the next step:
   ```bash
   squeue -u $USER
   ```

2. **Monitor output files** for errors:
   ```bash
   tail -f *.out
   tail -f *.err
   ```

3. **Count outputs** to verify completion:
   ```bash
   # Count PDBs in a directory
   ls -1 *.pdb | wc -l

   # Count FASTA files
   find . -name "*.fa" | wc -l
   ```

4. **Use dry-run mode** to test without submitting:
   ```bash
   python design/scripts/run_design_pipeline.py config.txt --dry-run
   ```

5. **Check disk space** before large runs:
   ```bash
   df -h /scratch/alpine/$USER
   ```
