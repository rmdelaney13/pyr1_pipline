# Design Pipeline Setup Checklist

Follow these steps to set up the automated design pipeline on your cluster.

---

## ☐ 1. Sync Pipeline Files to Cluster

Copy these new files from your local machine to the cluster:

```bash
# On your local machine:
cd "C:\Users\rmdel\OneDrive - UCB-O365\Whitehead Lab\pyr1_pipeline"

# Copy to cluster
scp -r design/scripts/ youruser@login.rc.colorado.edu:/projects/youruser/pyr1_pipeline/design/
scp design/CONFIG_TEMPLATE.txt youruser@login.rc.colorado.edu:/projects/youruser/pyr1_pipeline/design/
scp design/DESIGN_PIPELINE_README.md youruser@login.rc.colorado.edu:/projects/youruser/pyr1_pipeline/design/
scp design/QUICKSTART.md youruser@login.rc.colorado.edu:/projects/youruser/pyr1_pipeline/design/
```

Or use rsync for easier syncing:
```bash
rsync -avz design/ youruser@login.rc.colorado.edu:/projects/youruser/pyr1_pipeline/design/
```

---

## ☐ 2. Verify Helper Scripts Exist

Check that these scripts are available (they should be in your LigandMPNN repo):

```bash
# On cluster:
ls /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/split_and_mutate_to_fasta.py
ls /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/make_af3_jsons.py
```

If missing, you may need to sync them from the LigandMPNN repository.

---

## ☐ 3. Create Config Section

Add the `[design]` section to your campaign config file. Use [CONFIG_TEMPLATE.txt](CONFIG_TEMPLATE.txt) as a guide.

**Example for xanthurenic acid:**

```bash
# On cluster, edit your config
vim /scratch/alpine/ryde3462/xan_design/config.txt
```

Add:
```ini
[design]
DesignRoot = design
DesignIterationRounds = 1

# MPNN settings
MPNNScript = %(PIPE_ROOT)s/design/instructions/ligand_alignment_mpnni_grouped.sh
MPNNBatchSize = 40
MPNNTemperature = 0.3
DesignResidues = 59 79 81 90 92 106 108 115 118 120 139 157 158 161 162 165
MPNNOmitDesignFile = /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/LCA_omit_design_59.json
MPNNBiasFile = /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/instructions/LCA_bias_AA_per_residue.json

# Rosetta settings
RosettaScript = %(PIPE_ROOT)s/design/instructions/submit_pyrosetta_general_threading_relax.sh
RosettaRelaxPy = %(PIPE_ROOT)s/design/instructions/general_relax.py
LigandParams = /projects/ryde3462/xanthurenic_acid/conformers/0/0.params

# Filtering
FilterScript = %(PIPE_ROOT)s/design/instructions/relax_2_filter__allpolar_unsats.py
FilterTargetN = 1000
FilterMaxUnsats = 1
FilterMaxPerParent = 20

# AF3
AF3BinaryTemplate = %(PIPE_ROOT)s/design/templates/pyr1_binary_template.json
AF3TernaryTemplate = %(PIPE_ROOT)s/design/templates/pyr1_ternary_template.json
LigandSDF = /projects/ryde3462/xanthurenic_acid/conformers/0/0.sdf

# Helper scripts
AggregateScoresScript = %(PIPE_ROOT)s/design/instructions/aggregate_scores.py
SplitToFastaScript = /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/split_and_mutate_to_fasta.py
MakeAF3JSONsScript = /projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/scripts/make_af3_jsons.py
```

---

## ☐ 4. Verify Docking Output Exists

```bash
# Check clustered docked PDBs
ls /scratch/alpine/ryde3462/xan_design/docked/clustered_final/*.pdb | wc -l
```

If empty, run docking first:
```bash
python /projects/ryde3462/pyr1_pipeline/docking/scripts/run_docking_workflow.py config.txt
```

---

## ☐ 5. Test Pipeline (Dry Run)

```bash
# Activate environment
conda activate ligand_alignment

# Test without submitting jobs
python /projects/ryde3462/pyr1_pipeline/design/scripts/run_design_pipeline.py \
    /scratch/alpine/ryde3462/xan_design/config.txt \
    --dry-run
```

**What to check:**
- ✓ Config loads without errors
- ✓ Paths are correct
- ✓ Custom scripts generated in output directories
- ✓ No missing file errors

---

## ☐ 6. Run First Iteration

```bash
# Full pipeline execution
python /projects/ryde3462/pyr1_pipeline/design/scripts/run_design_pipeline.py \
    /scratch/alpine/ryde3462/xan_design/config.txt
```

**Monitor progress:**
```bash
# Watch SLURM queue
squeue -u ryde3462

# Check MPNN outputs
watch -n 10 "find /scratch/alpine/ryde3462/xan_design/design/iteration_1/mpnn_output -name '*.fa' | wc -l"

# Check Rosetta outputs
watch -n 10 "find /scratch/alpine/ryde3462/xan_design/design/iteration_1/rosetta_output -name '*.pdb' | wc -l"
```

---

## ☐ 7. Verify Outputs

After pipeline completes:

```bash
# Check filtered designs
wc -l /scratch/alpine/ryde3462/xan_design/design/iteration_1/filtered/filtered.csv
ls /scratch/alpine/ryde3462/xan_design/design/iteration_1/filtered/*.pdb | wc -l

# Check FASTA
wc -l /scratch/alpine/ryde3462/xan_design/design/iteration_1/filtered/filtered.fasta

# Check AF3 JSONs
ls /scratch/alpine/ryde3462/xan_design/design/af3_inputs/binary/*.json | wc -l
ls /scratch/alpine/ryde3462/xan_design/design/af3_inputs/ternary/*.json | wc -l
```

**Expected outputs:**
- CSV with ~FilterTargetN designs (e.g., 1000)
- Same number of PDB files
- FASTA with sequences
- Binary and ternary JSON files for AF3

---

## ☐ 8. Verify SMILES in Templates

Check that SMILES was correctly inserted:

```bash
# Should show xanthurenic acid SMILES
grep -A 2 '"ligand"' /scratch/alpine/ryde3462/xan_design/design/af3_inputs/binary/template_with_smiles.json | grep smiles
```

Expected: `"smiles": "C1=CC2=C(C(=C1)O)NC(=CC2=O)C(=O)O"`

If wrong, manually update:
```bash
python /projects/ryde3462/pyr1_pipeline/design/scripts/update_template_smiles.py \
    /scratch/alpine/ryde3462/xan_design/design/af3_inputs/binary/template_with_smiles.json \
    --sdf /projects/ryde3462/xanthurenic_acid/conformers/0/0.sdf \
    --in-place
```

---

## ☐ 9. Optional: Test Second Iteration

If you want to test iteration:

```bash
# Edit config
vim /scratch/alpine/ryde3462/xan_design/config.txt
# Change: DesignIterationRounds = 2

# Run again
python /projects/ryde3462/pyr1_pipeline/design/scripts/run_design_pipeline.py \
    /scratch/alpine/ryde3462/xan_design/config.txt
```

Or run just iteration 2:
```bash
python /projects/ryde3462/pyr1_pipeline/design/scripts/run_design_pipeline.py \
    /scratch/alpine/ryde3462/xan_design/config.txt \
    --iteration 2
```

---

## ☐ 10. Document for Lab Members

Create a project-specific README:

```bash
# In your campaign directory
cat > /scratch/alpine/ryde3462/xan_design/DESIGN_README.md << 'EOF'
# Xanthurenic Acid Design Pipeline

## Quick Run

```bash
conda activate ligand_alignment
python /projects/ryde3462/pyr1_pipeline/design/scripts/run_design_pipeline.py \
    /scratch/alpine/ryde3462/xan_design/config.txt
```

## Outputs

- Filtered designs: `/scratch/alpine/ryde3462/xan_design/design/iteration_1/filtered/`
- AF3 inputs: `/scratch/alpine/ryde3462/xan_design/design/af3_inputs/`

## Next Steps

1. Submit AF3 binary jobs from `af3_inputs/binary/`
2. Submit AF3 ternary jobs from `af3_inputs/ternary/`
3. Analyze results with existing AF3 analysis scripts

See full docs: `/projects/ryde3462/pyr1_pipeline/design/DESIGN_PIPELINE_README.md`
EOF
```

---

## Troubleshooting

### Import Errors

```bash
# If docking_pipeline_utils not found:
cd /projects/ryde3462/pyr1_pipeline/design/scripts
export PYTHONPATH=/projects/ryde3462/pyr1_pipeline/docking/scripts:$PYTHONPATH
```

### Path Not Found Errors

Double-check all paths in config match your cluster:
```bash
grep -E "PIPE_ROOT|CAMPAIGN_ROOT|SCRATCH_ROOT" /scratch/alpine/ryde3462/xan_design/config.txt
```

### RDKit Errors

```bash
conda activate ligand_alignment
conda install -c conda-forge rdkit
```

### SLURM Job Failures

Check logs in output directories:
```bash
ls /scratch/alpine/ryde3462/xan_design/design/iteration_1/mpnn_output/*.err
ls /scratch/alpine/ryde3462/xan_design/design/iteration_1/rosetta_output/*.err
```

---

## Success Criteria

✅ **Pipeline is working correctly if:**
1. Dry run completes without errors
2. SLURM jobs submit successfully
3. MPNN generates FASTA files
4. Rosetta generates PDB and .sc files
5. Filtering produces expected number of designs
6. FASTA is generated
7. AF3 JSONs are created with correct SMILES
8. No manual directory editing required!

---

## Support

- Documentation: `/projects/ryde3462/pyr1_pipeline/design/DESIGN_PIPELINE_README.md`
- Quick ref: `/projects/ryde3462/pyr1_pipeline/design/QUICKSTART.md`
- Config template: `/projects/ryde3462/pyr1_pipeline/design/CONFIG_TEMPLATE.txt`
