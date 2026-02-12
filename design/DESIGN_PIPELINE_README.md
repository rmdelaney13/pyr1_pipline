# Design Pipeline: Automated MPNN → Rosetta → AF3 Prep

This pipeline automates the design workflow from docked structures to AF3-ready sequences.

## Overview

```
Docked PDBs → LigandMPNN → Rosetta Relax → Filter → (Optional: Iterate) → FASTA → AF3 JSONs
```

**What this automates:**
- ✅ LigandMPNN sequence design (no manual directory editing!)
- ✅ Rosetta threading and relaxation (no manual directory editing!)
- ✅ Score aggregation
- ✅ Filtering by Rosetta metrics
- ✅ Optional iteration (1-2 rounds)
- ✅ FASTA generation
- ✅ AF3 JSON preparation with automatic SMILES insertion

**What remains manual:**
- ⚠️ AF3 execution (GPU cluster)
- ⚠️ AF3 analysis

---

## Quick Start

### 1. Setup Config File

Add the `[design]` section to your config.txt. See [CONFIG_TEMPLATE.txt](CONFIG_TEMPLATE.txt) for the full template.

**Minimum required settings:**
```ini
[design]
DesignRoot = design
DesignIterationRounds = 1
DesignResidues = 59 79 81 90 92 106 108 115 118 120 139 157 158 161 162 165
LigandParams = %(CAMPAIGN_ROOT)s/conformers/0/0.params
LigandSDF = %(CAMPAIGN_ROOT)s/conformers/0/0.sdf
```

### 2. Ensure Docking is Complete

The design pipeline starts from clustered docked PDBs:

```bash
# Check that docking output exists
ls $SCRATCH_ROOT/docked/clustered_final/*.pdb
```

If not, run docking first:
```bash
python docking/scripts/run_docking_workflow.py config.txt
```

### 3. Run Design Pipeline

**One-command execution:**
```bash
# Full pipeline: MPNN → Rosetta → Filter → AF3 prep
python design/scripts/run_design_pipeline.py config.txt
```

**This will:**
1. Generate custom MPNN submit script with correct paths
2. Submit MPNN array job to SLURM
3. Generate custom Rosetta submit script with correct paths
4. Submit Rosetta array job to SLURM
5. Aggregate Rosetta scores to CSV
6. Filter designs by metrics (dG_sep, unsats, charge, etc.)
7. Generate FASTA from filtered PDBs
8. Create AF3 binary + ternary JSON inputs with SMILES

---

## Output Structure

```
$SCRATCH_ROOT/design/
├── iteration_1/
│   ├── mpnn_output/
│   │   ├── submit_mpnn.sh              # Auto-generated MPNN script
│   │   └── array001_mpnn/*.fa          # MPNN FASTA outputs
│   │
│   ├── rosetta_output/
│   │   ├── submit_rosetta.sh           # Auto-generated Rosetta script
│   │   ├── array001_design_0.pdb       # Relaxed PDBs
│   │   └── array001_design_0.sc        # Score files
│   │
│   ├── scores/
│   │   └── iteration_1_scores.csv      # Aggregated scores
│   │
│   └── filtered/
│       ├── filtered.csv                # Filtered designs CSV
│       ├── filtered.fasta              # Final sequences
│       └── *.pdb                       # Filtered PDBs
│
├── iteration_2/                        # If DesignIterationRounds >= 2
│   └── [same structure]
│
└── af3_inputs/
    ├── binary/
    │   ├── template_with_smiles.json   # Template with SMILES
    │   └── *.json                      # AF3 input JSONs
    │
    └── ternary/
        ├── template_with_smiles.json
        └── *.json
```

---

## Stage-by-Stage Execution

### Run Specific Iteration

```bash
# Only run iteration 2 (assumes iteration 1 is complete)
python design/scripts/run_design_pipeline.py config.txt --iteration 2
```

### Skip Stages

```bash
# Skip MPNN (assumes already run)
python design/scripts/run_design_pipeline.py config.txt --skip-mpnn

# Skip Rosetta
python design/scripts/run_design_pipeline.py config.txt --skip-rosetta

# Only prepare AF3 inputs (assumes everything else done)
python design/scripts/run_design_pipeline.py config.txt --af3-prep-only
```

### Dry Run (Generate Scripts Without Submitting)

```bash
# See what would be run without actually submitting SLURM jobs
python design/scripts/run_design_pipeline.py config.txt --dry-run
```

This creates the custom submit scripts in the output directories, which you can inspect and manually submit if desired.

---

## Iteration Workflow

### Single Iteration (Default)

```
Iteration 1:
  Input:  Clustered docked PDBs
  Output: Filtered designs (e.g., top 1000)
```

### Two Iterations

```
Iteration 1:
  Input:  Clustered docked PDBs
  MPNN → Rosetta → Filter (top 1000)

Iteration 2:
  Input:  Filtered PDBs from iteration 1
  MPNN → Rosetta → Filter (top 1000 from this round)
```

**Why iterate?**
- Refine sequences based on Rosetta energy landscape
- Potentially find better designs than single-pass
- Typically diminishing returns after 2 iterations

**To enable:**
```ini
[design]
DesignIterationRounds = 2
```

---

## Filtering Parameters

Designs are filtered based on Rosetta metrics. Adjust in config.txt:

```ini
[design]
# Keep top 1000 designs
FilterTargetN = 1000

# Maximum buried unsatisfied polars (lower = stricter)
FilterMaxUnsats = 1

# Maximum designs per parent dock (for diversity)
FilterMaxPerParent = 20

# Optional: Relax specific requirements
FilterIgnoreO1 = False      # Require O1 polar contact
FilterIgnoreO2 = False      # Require O2 polar contact
FilterIgnoreCharge = False  # Require charge satisfaction
```

**Common adjustments:**
- **Too few designs pass:** Increase `FilterMaxUnsats`, decrease `FilterTargetN`
- **Too many similar designs:** Decrease `FilterMaxPerParent`
- **Ligand-specific filters failing:** Set `FilterIgnoreO1=True` or `FilterIgnoreO2=True`

---

## SMILES Handling

### Automatic SMILES Extraction (Recommended)

Specify your ligand SDF in config:
```ini
[design]
LigandSDF = %(CAMPAIGN_ROOT)s/conformers/0/0.sdf
```

The pipeline will:
1. Extract SMILES from SDF using RDKit
2. Update AF3 templates with the SMILES
3. Generate JSONs with correct ligand structure

### Manual SMILES Update

If you need to manually update templates:

```bash
# Extract SMILES from SDF
python design/scripts/extract_smiles.py conformers/0/0.sdf
# Output: C1=CC2=C(C(=C1)O)NC(=CC2=O)C(=O)O

# Update template JSON
python design/scripts/update_template_smiles.py \
  design/templates/pyr1_binary_template.json \
  --smiles "C1=CC2=C(C(=C1)O)NC(=CC2=O)C(=O)O" \
  --output updated_template.json
```

---

## Monitoring Progress

### Check SLURM Jobs

```bash
# View running jobs
squeue -u $USER

# Monitor specific job
squeue -j <job_id>

# Check output logs
tail -f *.out
tail -f *.err
```

### Check Outputs

```bash
# Count MPNN FASTA outputs
find $SCRATCH_ROOT/design/iteration_1/mpnn_output -name "*.fa" | wc -l

# Count Rosetta PDB outputs
find $SCRATCH_ROOT/design/iteration_1/rosetta_output -name "*.pdb" | wc -l

# View filtered designs summary
head $SCRATCH_ROOT/design/iteration_1/filtered/filtered.csv
```

### Common Issues

**"No PDB files found in input directory"**
- Check that docking completed successfully
- Verify `OutputDir` and `ClusteredOutputDir` in config.txt

**"No FASTA files found in MPNN output"**
- MPNN jobs may still be running or failed
- Check SLURM logs in `mpnn_output/` directory

**"No .sc files found"**
- Rosetta jobs may still be running or failed
- Check SLURM logs in `rosetta_output/` directory

**"Too few designs passed filters"**
- Relax filtering criteria (increase `FilterMaxUnsats`)
- Check score distribution in `scores/iteration_1_scores.csv`

---

## Manual Workflow (Pre-Automation)

For reference, here's what you used to do manually:

### Old Way ❌
```bash
# 1. Edit ligand_alignment_mpnni_grouped.sh
#    - Change PDB_DIR="/scratch/.../filtered_1"
#    - Change OUTPUT_BASE="/scratch/.../mpnn_output_1"
#    - Update --array count
sbatch ligand_alignment_mpnni_grouped.sh

# 2. Edit submit_pyrosetta_general_threading_relax.sh
#    - Change TEMPLATE_DIR="/scratch/.../filtered_1"
#    - Change MPNN_OUTPUT_BASE="/scratch/.../mpnn_output_1"
#    - Change OUTPUT_DIR="/scratch/.../relax_1"
#    - Update --array count
sbatch submit_pyrosetta_general_threading_relax.sh

# 3. Manually aggregate scores
python aggregate_scores.py /scratch/.../relax_1 --output relax_1.csv

# 4. Manually filter
python relax_2_filter__allpolar_unsats.py \
  relax_1.csv /scratch/.../relax_1 /scratch/.../filtered_1 \
  --target_n 1000 --max_unsat 1 --max_per_parent 20

# 5. Manually generate FASTA
python split_and_mutate_to_fasta.py \
  /scratch/.../filtered_1 /scratch/.../filtered_1/filtered.fasta

# 6. Manually create AF3 JSONs (forgetting to update SMILES!)
python make_af3_jsons.py --template old_template.json --fasta filtered.fasta --outdir af3_input

# 7. Repeat for iteration 2...
```

### New Way ✅
```bash
# Everything in one command:
python design/scripts/run_design_pipeline.py config.txt
```

---

## Advanced: Custom Filtering

If you need custom filtering logic beyond the standard script:

1. Create a new filter script (e.g., `custom_filter.py`)
2. Update config:
   ```ini
   [design]
   FilterScript = %(PIPE_ROOT)s/design/instructions/custom_filter.py
   ```
3. Ensure it accepts the same arguments as `relax_2_filter__allpolar_unsats.py`

---

## Integration with Docking Pipeline

### Full Workflow: SDF → AF3-Ready

```bash
# Step 1: Docking (automated)
python docking/scripts/run_docking_workflow.py config.txt

# Step 2: Design (automated)
python design/scripts/run_design_pipeline.py config.txt

# Step 3: AF3 (manual - GPU)
# Submit JSONs from $SCRATCH_ROOT/design/af3_inputs/binary to AF3 cluster
# Submit JSONs from $SCRATCH_ROOT/design/af3_inputs/ternary to AF3 cluster

# Step 4: AF3 Analysis (manual)
# Use existing analysis scripts on AF3 outputs
```

### SLURM Job Dependencies (Advanced)

To chain docking → design automatically:

```bash
# Submit docking
DOCK_JOB=$(sbatch --parsable docking/scripts/submit_complete_workflow.sh config.txt)

# Submit design after docking completes
sbatch --dependency=afterok:$DOCK_JOB design/scripts/submit_design_pipeline.sh config.txt
```

(Note: `submit_design_pipeline.sh` would need to be created as a SLURM wrapper)

---

## Troubleshooting

### RDKit Import Errors

If SMILES extraction fails:
```bash
conda install -c conda-forge rdkit
```

### Path Issues

If scripts can't find files:
```bash
# Check config paths
grep -E "PIPE_ROOT|CAMPAIGN_ROOT|SCRATCH_ROOT" config.txt

# Verify docking outputs exist
ls $(grep OutputDir config.txt | cut -d'=' -f2 | xargs)/clustered_final/
```

### SLURM Job Failures

Check logs in iteration output directories:
```bash
cd $SCRATCH_ROOT/design/iteration_1/mpnn_output
ls *.out *.err

cd $SCRATCH_ROOT/design/iteration_1/rosetta_output
ls *.out *.err
```

### Memory/Time Issues

Edit generated submit scripts and resubmit:
```bash
# Increase memory/time in auto-generated script
vim $SCRATCH_ROOT/design/iteration_1/mpnn_output/submit_mpnn.sh
# Change #SBATCH --mem=8G to #SBATCH --mem=16G

# Manually submit
sbatch $SCRATCH_ROOT/design/iteration_1/mpnn_output/submit_mpnn.sh
```

---

## Next Steps

After the design pipeline completes:

1. **Review filtered designs:**
   ```bash
   # View top designs
   head -20 $SCRATCH_ROOT/design/iteration_1/filtered/filtered.csv

   # Visualize in PyMOL
   pymol $SCRATCH_ROOT/design/iteration_1/filtered/*.pdb
   ```

2. **Submit to AF3:**
   - Copy JSONs to GPU cluster
   - Submit binary predictions (no HAB)
   - Submit ternary predictions (with HAB)

3. **Analyze AF3 results:**
   - Use existing `binary_analysis.py` and `ternary_analysis.py` scripts
   - Filter by pLDDT, ipTM, RMSD
   - Select top candidates for experimental testing

4. **Optional: Another iteration:**
   ```bash
   # Run iteration 2 on best AF3 results
   # (requires manual selection of top AF3 designs as input)
   ```

---

## Reference: Config Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DesignRoot` | `design` | Output directory name |
| `DesignIterationRounds` | `1` | Number of MPNN→Rosetta cycles (1-2) |
| `DesignResidues` | (see template) | Residues to design (space-separated) |
| `LigandParams` | *(required)* | Ligand params file path |
| `LigandSDF` | *(optional)* | SDF for SMILES extraction |
| `FilterTargetN` | `1000` | Target number of designs after filtering |
| `FilterMaxUnsats` | `1` | Max buried unsatisfied polars |
| `FilterMaxPerParent` | `20` | Max designs per parent (diversity) |
| `MPNNBatchSize` | `40` | Sequences per parent |
| `MPNNTemperature` | `0.3` | MPNN sampling temperature |
| `AF3BinaryTemplate` | (see template) | Binary complex template JSON |
| `AF3TernaryTemplate` | (see template) | Ternary complex template JSON |

---

## Contact

For issues or questions:
- Check existing documentation in `docking/` directory
- Review SLURM logs in output directories
- Verify config.txt settings match your cluster paths
