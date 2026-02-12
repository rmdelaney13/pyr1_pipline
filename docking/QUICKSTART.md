# Docking Workflow - Quick Start Guide

## One-Command Execution Options

### 🚀 Option 1: Single Local Run

```bash
python scripts/run_docking_workflow.py config.txt
```

**When to use**: Small test runs, debugging, or single-molecule docking

---

### 🚀 Option 2: SLURM Array Job - TRUE One Command! (Recommended)

```bash
# Single command submits everything with automatic dependencies
bash scripts/submit_complete_workflow.sh config.txt
```

**When to use**: Large campaigns with many conformers (recommended for production)

**What it does**:
- Submits array job (reads `ArrayTaskCount` from config.txt)
- Auto-submits clustering job that runs after all arrays complete
- Uses SLURM job dependencies - no manual intervention needed!

**Alternative (manual two-step)**:
```bash
# Step 1: Submit array job
sbatch --array=0-9 scripts/submit_docking_workflow.sh config.txt
# Step 2: After completion, submit clustering
sbatch scripts/run_clustering_only.sh config.txt
```

---

### 🚀 Option 3: Local Array Test

```bash
python scripts/run_docking_workflow.py config.txt --local-arrays 4
```

**When to use**: Testing array distribution locally without SLURM

---

## What Gets Run

Your single command executes these three steps automatically:

1. **create_table.py** - Generates params/PDBs and alignment table from SDF
2. **grade_conformers_glycine_shaved_docking_multiple_slurm.py** - Docks conformers with filtering (runs N times for array jobs)
3. **cluster_docked_post_array.py** - Clusters final results

---

## Before You Run

### 1. Configure `config.txt`

Edit these key sections:

```ini
[DEFAULT]
CAMPAIGN_ROOT = /path/to/your/campaign
SCRATCH_ROOT = /scratch/path/for/outputs

[create_table]
MoleculeSDFs = %(CAMPAIGN_ROOT)s/conformers/*.sdf
TargetAtomTriplets = O2-C11-C9; O2-C9-C11

[grade_conformers]
ArrayTaskCount = 10  # Set to 1 for single run, 10+ for arrays
OutputDir = %(SCRATCH_ROOT)s/docked
```

### 2. Prepare Input Files

- ✅ SDF conformer file(s) in `MoleculeSDFs` location
- ✅ Pre/Post PDB templates (usually from pipeline repo)
- ✅ Params files if your protein has non-standard residues

### 3. Activate Environment

```bash
conda activate ligand_alignment  # or your PyRosetta environment
```

---

## After Running

### Check Results

```bash
# View final clustered results
ls $SCRATCH_ROOT/docked/clustered_final/*.pdb

# Review cluster summary
cat $SCRATCH_ROOT/docked/clustered_final/cluster_summary.csv

# Analyze sampling adequacy
python scripts/analyze_sampling_adequacy.py config.txt --ligand-heavy-atoms 15
```

### Typical Output Structure

```
$SCRATCH_ROOT/docked/
├── a0000_rep_0001.pdb              # Docked poses from array 0
├── a0001_rep_0001.pdb              # Docked poses from array 1
├── hbond_geometry_summary_array0000.csv
├── hbond_geometry_summary_array0001.csv
└── clustered_final/
    ├── cluster_0001_*.pdb          # Best representative from cluster 1
    ├── cluster_0002_*.pdb          # Best representative from cluster 2
    └── cluster_summary.csv         # Scores and cluster info
```

---

## Common Commands

### Re-run clustering with different RMSD cutoff

```bash
# Edit config.txt: ClusterRMSDCutoff = 1.0
python scripts/run_docking_workflow.py config.txt --skip-create-table --skip-docking
```

### Re-run a specific array task that failed

```bash
python scripts/run_docking_workflow.py config.txt --array-index 5 --skip-create-table --skip-clustering
```

### Only prepare table (don't dock)

```bash
python scripts/run_docking_workflow.py config.txt --prepare-only
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "No candidate PDB files found" | Relax filters in config (MaxScore, HBond params) |
| "RDKit not found" | `conda install -c conda-forge rdkit` |
| Array job only runs 1 task | Add `--array=0-N` to sbatch command |
| Out of memory | Increase `--mem` in submit_docking_workflow.sh |

---

## Need More Details?

See [WORKFLOW_README.md](WORKFLOW_README.md) for comprehensive documentation.
