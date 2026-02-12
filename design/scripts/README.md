# Design Pipeline Scripts

This directory contains the automated design pipeline orchestrator and helper scripts.

---

## Main Script

### **run_design_pipeline.py**

Master orchestrator for the MPNN → Rosetta → Filter → AF3 prep workflow.

**Usage:**
```bash
python run_design_pipeline.py config.txt
```

**What it does:**
1. Reads config to determine iteration count, paths, and parameters
2. For each iteration:
   - Generates custom MPNN submit script with correct paths/array count
   - Submits MPNN SLURM job
   - Generates custom Rosetta submit script with correct paths/array count
   - Submits Rosetta SLURM job
   - Aggregates Rosetta scores to CSV
   - Filters designs by metrics
3. After all iterations:
   - Generates FASTA from final filtered PDBs
   - Extracts SMILES from ligand SDF
   - Updates AF3 templates with SMILES
   - Generates binary and ternary AF3 JSONs

**Key features:**
- Automatic path management (no manual script editing!)
- Array count calculation based on input files
- Iteration support (1-2 rounds)
- SMILES auto-extraction and insertion
- Stage-by-stage execution with skip flags

**Options:**
- `--iteration N` - Run specific iteration only
- `--skip-mpnn` - Skip LigandMPNN stage
- `--skip-rosetta` - Skip Rosetta stage
- `--skip-aggregate` - Skip score aggregation
- `--skip-filter` - Skip filtering
- `--skip-af3-prep` - Skip AF3 preparation
- `--af3-prep-only` - Only prepare AF3 inputs
- `--dry-run` - Generate scripts but don't submit jobs
- `--wait` - Wait for SLURM jobs to complete

---

## Helper Scripts

### **extract_smiles.py**

Extract SMILES string from SDF file.

**Usage:**
```bash
python extract_smiles.py ligand.sdf
# Output: C1=CC2=C(C(=C1)O)NC(=CC2=O)C(=O)O

python extract_smiles.py ligand.sdf --output smiles.txt
```

**Requirements:**
- RDKit: `conda install -c conda-forge rdkit`

---

### **update_template_smiles.py**

Update SMILES in AF3 JSON template.

**Usage:**
```bash
# From SDF
python update_template_smiles.py template.json \
    --sdf ligand.sdf \
    --output updated_template.json

# From SMILES string
python update_template_smiles.py template.json \
    --smiles "C1=CC2=C(C(=C1)O)NC(=CC2=O)C(=O)O" \
    --in-place
```

**What it does:**
- Finds `sequences[N].ligand.smiles` in JSON
- Replaces with new SMILES
- Validates update was successful

---

## Examples

### Full Pipeline
```bash
python run_design_pipeline.py /scratch/alpine/ryde3462/xan_design/config.txt
```

### Two Iterations
```bash
# In config.txt: DesignIterationRounds = 2
python run_design_pipeline.py config.txt
```

### Run Only Iteration 2
```bash
python run_design_pipeline.py config.txt --iteration 2
```

### Dry Run (Generate Scripts, Don't Submit)
```bash
python run_design_pipeline.py config.txt --dry-run
```

### Only AF3 Prep (Manual MPNN/Rosetta Already Done)
```bash
python run_design_pipeline.py config.txt --af3-prep-only
```

---

## How It Works

### 1. Configuration Loading

```python
cfg = DesignPipelineConfig(config_file)
# Reads [design] section
# Sets up all paths and parameters
```

### 2. Iteration Loop

For each iteration (1 to N):

**MPNN Stage:**
- Reads template MPNN script
- Replaces `PDB_DIR`, `OUTPUT_BASE` with iteration-specific paths
- Calculates array count from input PDB count
- Writes custom submit script
- Submits to SLURM

**Rosetta Stage:**
- Reads template Rosetta script
- Replaces `TEMPLATE_DIR`, `MPNN_OUTPUT_BASE`, `OUTPUT_DIR`
- Calculates array count from FASTA count
- Writes custom submit script
- Submits to SLURM

**Aggregation:**
- Scans Rosetta output for `.sc` files
- Aggregates to single CSV

**Filtering:**
- Filters by dG_sep, unsats, charge, etc.
- Limits per parent for diversity
- Copies top N PDBs to filtered directory

### 3. AF3 Preparation

After all iterations:
- Generate FASTA from final filtered PDBs
- Extract SMILES from SDF (if provided)
- Update templates with SMILES
- Generate binary and ternary JSONs

---

## Directory Management

The pipeline automatically creates and manages:

```
$SCRATCH_ROOT/design/
├── iteration_1/
│   ├── mpnn_output/
│   ├── rosetta_output/
│   ├── scores/
│   └── filtered/
├── iteration_2/
│   └── ...
└── af3_inputs/
    ├── binary/
    └── ternary/
```

**No manual directory creation needed!**

---

## Error Handling

The pipeline validates:
- Config file exists and is readable
- Input directories exist and contain expected files
- Required scripts exist
- SLURM submissions succeed

**Exit codes:**
- 0 = Success
- 1 = Error occurred

---

## Requirements

- Python 3.6+
- Access to SLURM cluster
- Conda environments:
  - `ligandmpnn_env` (for MPNN)
  - `ligand_alignment` (for Rosetta)
- RDKit (for SMILES extraction)
- Existing scripts:
  - LigandMPNN submit template
  - Rosetta submit template
  - aggregate_scores.py
  - filter script
  - split_and_mutate_to_fasta.py
  - make_af3_jsons.py

---

## Extending the Pipeline

### Add Custom Filter

1. Create new filter script (e.g., `custom_filter.py`)
2. Update config:
   ```ini
   [design]
   FilterScript = %(PIPE_ROOT)s/design/instructions/custom_filter.py
   ```
3. Ensure it accepts same arguments as `relax_2_filter__allpolar_unsats.py`

### Add Pre/Post-Processing Steps

Modify `run_design_pipeline.py`:

```python
# In run_full_pipeline():

# Before MPNN
if not args.skip_mpnn:
    preprocess_pdbs(cfg, iter_num)  # Your function
    run_mpnn_design(cfg, iter_num)

# After filtering
if not args.skip_filter:
    filter_designs(cfg, iter_num, scores_csv)
    postprocess_filtered(cfg, iter_num)  # Your function
```

---

## Debugging

### Enable Verbose Output

```python
# In run_design_pipeline.py, add:
import logging
logging.basicConfig(level=logging.DEBUG)
```

### Check Generated Scripts

After dry run:
```bash
# View generated MPNN script
cat $SCRATCH_ROOT/design/iteration_1/mpnn_output/submit_mpnn.sh

# View generated Rosetta script
cat $SCRATCH_ROOT/design/iteration_1/rosetta_output/submit_rosetta.sh
```

### Manual Submission

If auto-submission fails:
```bash
# Submit manually
sbatch $SCRATCH_ROOT/design/iteration_1/mpnn_output/submit_mpnn.sh
sbatch $SCRATCH_ROOT/design/iteration_1/rosetta_output/submit_rosetta.sh
```

---

## See Also

- [../DESIGN_PIPELINE_README.md](../DESIGN_PIPELINE_README.md) - Complete documentation
- [../QUICKSTART.md](../QUICKSTART.md) - Quick reference
- [../CONFIG_TEMPLATE.txt](../CONFIG_TEMPLATE.txt) - Config template
- [../SETUP_CHECKLIST.md](../SETUP_CHECKLIST.md) - Setup guide
