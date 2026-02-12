# Design Pipeline - Quick Start

## One-Command Workflow

```bash
# After docking is complete, run:
python design/scripts/run_design_pipeline.py config.txt
```

This single command executes:
1. ✅ LigandMPNN design
2. ✅ Rosetta relax
3. ✅ Score aggregation
4. ✅ Filtering
5. ✅ FASTA generation
6. ✅ AF3 JSON preparation

---

## Setup (One-Time)

### 1. Add Config Section

Copy from [CONFIG_TEMPLATE.txt](CONFIG_TEMPLATE.txt) to your `config.txt`:

```ini
[design]
DesignRoot = design
DesignIterationRounds = 1
DesignResidues = 59 79 81 90 92 106 108 115 118 120 139 157 158 161 162 165
LigandParams = %(CAMPAIGN_ROOT)s/conformers/0/0.params
LigandSDF = %(CAMPAIGN_ROOT)s/conformers/0/0.sdf
FilterTargetN = 1000
FilterMaxUnsats = 1
FilterMaxPerParent = 20
```

### 2. Verify Docking Complete

```bash
ls $SCRATCH_ROOT/docked/clustered_final/*.pdb
```

---

## Common Commands

```bash
# Full pipeline
python design/scripts/run_design_pipeline.py config.txt

# Specific iteration only
python design/scripts/run_design_pipeline.py config.txt --iteration 2

# Dry run (generate scripts, don't submit)
python design/scripts/run_design_pipeline.py config.txt --dry-run

# Only AF3 prep (after manual MPNN/Rosetta)
python design/scripts/run_design_pipeline.py config.txt --af3-prep-only

# Skip stages
python design/scripts/run_design_pipeline.py config.txt --skip-mpnn
python design/scripts/run_design_pipeline.py config.txt --skip-rosetta
```

---

## Output Locations

```bash
# Check MPNN outputs
ls $SCRATCH_ROOT/design/iteration_1/mpnn_output/

# Check Rosetta outputs
ls $SCRATCH_ROOT/design/iteration_1/rosetta_output/

# Check filtered designs
cat $SCRATCH_ROOT/design/iteration_1/filtered/filtered.csv
cat $SCRATCH_ROOT/design/iteration_1/filtered/filtered.fasta

# Check AF3 inputs
ls $SCRATCH_ROOT/design/af3_inputs/binary/*.json
ls $SCRATCH_ROOT/design/af3_inputs/ternary/*.json
```

---

## Monitoring

```bash
# Check SLURM jobs
squeue -u $USER

# Count outputs
find $SCRATCH_ROOT/design/iteration_1/mpnn_output -name "*.fa" | wc -l
find $SCRATCH_ROOT/design/iteration_1/rosetta_output -name "*.pdb" | wc -l
find $SCRATCH_ROOT/design/iteration_1/filtered -name "*.pdb" | wc -l
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "No PDB files found" | Check docking completed: `ls $SCRATCH_ROOT/docked/clustered_final/` |
| "No FASTA files found" | MPNN still running or failed. Check logs in `mpnn_output/*.err` |
| Too few designs pass | Increase `FilterMaxUnsats` or decrease `FilterTargetN` |
| SMILES not updating | Verify `LigandSDF` path in config, install RDKit |

---

## Full Documentation

See [DESIGN_PIPELINE_README.md](DESIGN_PIPELINE_README.md) for complete details.

---

## Old vs New Workflow

### Before (Manual) ❌
```bash
# Edit ligand_alignment_mpnni_grouped.sh (directories, array count)
sbatch ligand_alignment_mpnni_grouped.sh

# Edit submit_pyrosetta_general_threading_relax.sh (directories, array count)
sbatch submit_pyrosetta_general_threading_relax.sh

python aggregate_scores.py /scratch/.../relax_1 --output relax_1.csv
python relax_2_filter__allpolar_unsats.py relax_1.csv ...
python split_and_mutate_to_fasta.py /scratch/.../filtered_1 filtered.fasta
python make_af3_jsons.py --template ... --fasta ... --outdir ...

# Manually edit SMILES in templates
# Repeat for iteration 2...
```

### After (Automated) ✅
```bash
python design/scripts/run_design_pipeline.py config.txt
```

---

## Next Steps After Pipeline

1. **Review designs:**
   ```bash
   head -20 $SCRATCH_ROOT/design/iteration_1/filtered/filtered.csv
   ```

2. **Submit to AF3 (GPU cluster):**
   - Binary: `$SCRATCH_ROOT/design/af3_inputs/binary/*.json`
   - Ternary: `$SCRATCH_ROOT/design/af3_inputs/ternary/*.json`

3. **Analyze AF3 results** using existing scripts

---

## Config Quick Reference

Essential settings:
- `DesignResidues` - Which residues to design
- `LigandParams` - Path to ligand params file
- `FilterTargetN` - How many designs to keep (default: 1000)
- `FilterMaxUnsats` - Strictness (lower = stricter, default: 1)

Iteration:
- `DesignIterationRounds = 1` - Single pass
- `DesignIterationRounds = 2` - Two rounds of refinement
