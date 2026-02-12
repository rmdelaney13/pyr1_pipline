# Unified Docking Workflow

This directory contains a unified workflow for running the complete docking pipeline from SDF input to clustered results.

## Overview

The workflow consists of three main steps:

1. **create_table.py** - Generate params/PDBs and alignment table from SDF conformers
2. **grade_conformers_glycine_shaved_docking_multiple_slurm.py** - Dock conformers with filtering
3. **cluster_docked_post_array.py** - Cluster final results across all array tasks

## Quick Start

### Option 1: Single Local Run (No SLURM)

```bash
# Complete workflow in one command
python scripts/run_docking_workflow.py config.txt
```

### Option 2: SLURM Array Job - Complete Automation (Recommended for Large Campaigns)

```bash
# ONE COMMAND - submits everything with automatic dependencies!
bash scripts/submit_complete_workflow.sh config.txt
```

This automatically:
- Reads `ArrayTaskCount` from your config.txt
- Submits array jobs 0-(N-1) for parallel docking
- Submits clustering job with dependency (runs after arrays complete)
- No manual intervention needed!

**Alternative manual two-step approach**:
```bash
# Step 1: Submit array job (replace 9 with ArrayTaskCount-1)
sbatch --array=0-9 scripts/submit_docking_workflow.sh config.txt

# Step 2: After all array tasks complete, run clustering
sbatch scripts/run_clustering_only.sh config.txt
# Or run directly:
python scripts/run_docking_workflow.py config.txt --skip-create-table --skip-docking
```

### Option 3: Test Locally with Multiple Arrays

```bash
# Simulate 4 array tasks locally (useful for debugging)
python scripts/run_docking_workflow.py config.txt --local-arrays 4
```

## Configuration

Edit your `config.txt` file to set key parameters:

### Essential Settings

```ini
[DEFAULT]
# Root directories
PIPE_ROOT = /path/to/pyr1_pipeline
CAMPAIGN_ROOT = /path/to/your/campaign
SCRATCH_ROOT = /scratch/path/for/outputs

# Input structure templates
PrePDBFileName = %(PIPE_ROOT)s/docking/ligand_alignment/files_for_PYR1_docking/3QN1_H2O.pdb
PostPDBFileName = %(PIPE_ROOT)s/docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb

[create_table]
# Input ligand conformers (supports wildcards)
MoleculeSDFs = %(CAMPAIGN_ROOT)s/conformers/*.sdf

# Dynamic alignment settings
DynamicAcceptorAlignment = True
MaxDynamicAlignments = 20
TargetAtomTriplets = O2-C11-C9; O2-C9-C11

[grade_conformers]
# SLURM array configuration
ArrayTaskCount = 10  # Set to 1 for single run, >1 for array jobs

# Output location
OutputDir = %(SCRATCH_ROOT)s/docked

# Clustering settings
EnablePoseClusteringInArrayTask = False  # False for array jobs
ClusterRMSDCutoff = 0.75

# Docking parameters
MaxScore = -300
EnableHBondGeometryFilter = True
MaxPerturbTries = 5
```

## Detailed Usage

### run_docking_workflow.py

The master orchestration script that coordinates all steps.

```bash
# Basic usage
python scripts/run_docking_workflow.py config.txt

# Advanced options
python scripts/run_docking_workflow.py config.txt \
    --mode glycine \                    # Docking mode (auto|glycine|sequence)
    --skip-create-table \              # Skip table creation if already done
    --skip-clustering \                # Skip final clustering
    --array-index 0 \                  # Run specific array index
    --local-arrays 10 \                # Test with 10 local array tasks
    --prepare-only                     # Only run create_table and exit
```

### Workflow Stages

#### Stage 1: Table Creation (create_table.py)

This step:
- Reads SDF conformer files
- Generates Rosetta params and PDB files
- Detects H-bond acceptor atoms using RDKit
- Creates alignment CSV with atom triplets
- Outputs: CSV table, PKL file, conformer PDBs, params files

**Skip this step** if you already have the CSV/params from a previous run:
```bash
python scripts/run_docking_workflow.py config.txt --skip-create-table
```

#### Stage 2: Docking (grade_conformers_glycine_shaved_docking_multiple_slurm.py)

This step:
- Loads conformers from the alignment table
- Perturbs and docks conformers into the binding pocket
- Applies H-bond geometry filters
- Scores poses with Rosetta
- Outputs: Passing docked PDBs, H-bond geometry CSV

**For SLURM arrays**, each array task processes a different subset of conformers.

#### Stage 3: Clustering (cluster_docked_post_array.py)

This step:
- Collects all passing PDBs from all array tasks
- Clusters based on ligand heavy-atom RMSD
- Keeps best-scoring representative per cluster
- Outputs: Clustered PDBs, cluster_summary.csv

**Run clustering separately** after all array tasks complete:
```bash
bash scripts/run_clustering_only.sh config.txt
```

## SLURM Array Jobs

### Why Use Arrays?

- **Parallelization**: Process thousands of conformers in parallel
- **Efficiency**: Each array task handles a subset of conformers
- **Fault Tolerance**: Individual task failures don't affect others
- **Resource Management**: Better cluster utilization

### Array Job Workflow

1. **Set ArrayTaskCount** in config.txt:
   ```ini
   [grade_conformers]
   ArrayTaskCount = 10  # Will split work across 10 tasks
   ```

2. **Submit array job**:
   ```bash
   # If ArrayTaskCount=10, submit with --array=0-9
   sbatch --array=0-9 scripts/submit_docking_workflow.sh config.txt
   ```

3. **Monitor progress**:
   ```bash
   squeue -u $USER
   tail -f docking_*.out
   ```

4. **After all tasks complete, run clustering**:
   ```bash
   sbatch scripts/run_clustering_only.sh config.txt
   ```

### Array Task Distribution

With `ArrayTaskCount=10`, conformers are distributed:
- Array 0: conformers 0, 10, 20, 30, ...
- Array 1: conformers 1, 11, 21, 31, ...
- Array 2: conformers 2, 12, 22, 32, ...
- ...
- Array 9: conformers 9, 19, 29, 39, ...

## Output Files

### After Docking (Per Array Task)

```
$OutputDir/
├── a0000_rep_0001.pdb          # Passing docked pose from array 0
├── a0000_rep_0002.pdb
├── a0001_rep_0001.pdb          # Passing docked pose from array 1
├── hbond_geometry_summary_array0000.csv  # H-bond metrics for array 0
└── hbond_geometry_summary_array0001.csv
```

### After Clustering

```
$OutputDir/clustered_final/
├── cluster_0001_a0003_rep_0042.pdb  # Best pose from cluster 1
├── cluster_0002_a0007_rep_0123.pdb  # Best pose from cluster 2
├── cluster_summary.csv              # Cluster representatives with scores
└── ...
```

The `cluster_summary.csv` contains:
- `cluster_id`: Cluster number
- `source_pdb`: Original PDB file name
- `score`: Rosetta score
- `output_pdb`: Path to clustered representative

## Troubleshooting

### "No candidate PDB files found"

**Cause**: No conformers passed docking filters
**Solution**:
- Check H-bond filter settings (relax `EnableHBondGeometryFilter`)
- Increase `MaxPerturbTries` to allow more sampling
- Relax `MaxScore` cutoff
- Review `hbond_geometry_summary*.csv` to see why conformers failed

### "create_table.py failed"

**Cause**: Missing RDKit or invalid SDF files
**Solution**:
- Ensure RDKit is installed: `conda install -c conda-forge rdkit`
- Validate SDF files: Open in a molecular viewer
- Check `MoleculeSDFs` path in config.txt

### Array job only runs 1 task

**Cause**: Forgot `--array` flag in sbatch
**Solution**:
```bash
# Wrong:
sbatch scripts/submit_docking_workflow.sh config.txt

# Correct:
sbatch --array=0-9 scripts/submit_docking_workflow.sh config.txt
```

### Out of memory errors

**Cause**: Large protein or many waters
**Solution**:
- Increase `--mem` in SLURM script (edit `submit_docking_workflow.sh`)
- Remove distant waters from PostPDBFileName
- Reduce `MaxPerturbTries`

## Examples

### Example 1: Small Test Run

```bash
# Edit config to use only 1 SDF file
# Set ArrayTaskCount = 1

python scripts/run_docking_workflow.py config.txt
```

### Example 2: Large Campaign with 100 Conformers

```bash
# Edit config:
# ArrayTaskCount = 10
# MoleculeSDFs = path/to/*.sdf (100 SDF files)

# Submit array job
sbatch --array=0-9 scripts/submit_docking_workflow.sh config.txt

# Wait for completion, then cluster
sbatch scripts/run_clustering_only.sh config.txt
```

### Example 3: Re-run Clustering with Different RMSD

```bash
# Edit config:
# ClusterRMSDCutoff = 1.0  (was 0.75)

# Re-run clustering only
python scripts/run_docking_workflow.py config.txt \
    --skip-create-table \
    --skip-docking
```

## Advanced: Custom Pipelines

### Run Only Table Creation

```bash
python scripts/run_docking_workflow.py config.txt --prepare-only
```

### Run Docking for Specific Array Index

```bash
# Useful for re-running failed array task
python scripts/run_docking_workflow.py config.txt \
    --array-index 5 \
    --skip-create-table \
    --skip-clustering
```

### Chain with Other Scripts

```bash
# Generate conformers with RDKit
python generate_conformers.py input.sdf -o conformers/

# Run docking workflow
python scripts/run_docking_workflow.py config.txt

# Analyze results
python analyze_clusters.py $OutputDir/clustered_final/cluster_summary.csv
```

## Related Scripts

- **run_docking_from_sdf.py** - Simpler workflow (table + docking, no clustering)
- **create_table.py** - Standalone table creation
- **grade_conformers_glycine_shaved.py** - Non-array docking script
- **cluster_docked_post_array.py** - Standalone clustering script

## Contact

For questions or issues, check the main pipeline README or contact the pipeline maintainer.
