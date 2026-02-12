# Workflow Options: Complete vs Individual Steps

This document provides a quick overview of how to run the pyr1_pipeline.

## Two Workflow Modes

### 🚀 Mode 1: Complete Automated Workflow
**Use when:** Running production pipelines, everything is working correctly

```bash
# Full pipeline: Docking → MPNN → Rosetta → AF3 Prep
bash docking/scripts/submit_complete_workflow.sh config.txt
# Wait for docking to complete...
python design/scripts/run_design_pipeline.py config.txt
```

**What it does:**
- Docking: Creates params, runs docking arrays, clusters results
- Design: Runs MPNN, Rosetta relax, filters, prepares AF3 inputs
- Everything is automated with SLURM job dependencies

---

### 🔧 Mode 2: Individual Steps
**Use when:** Troubleshooting, testing parameters, re-running specific steps

#### Docking Steps (Individual)
```bash
# Step 1: Create table only
python docking/scripts/create_table.py config.txt

# Step 2: Run docking (single array task for testing)
python docking/scripts/grade_conformers_glycine_shaved_docking_multiple_slurm.py config.txt 0

# Step 3: Cluster results
python docking/scripts/cluster_docked_post_array.py config.txt
```

#### Design Steps (Individual)
```bash
# Step 1: LigandMPNN only
bash design/scripts/run_mpnn_only.sh \
    /path/to/docked_pdbs \
    /path/to/mpnn_output \
    50

# Step 2: Rosetta only
bash design/scripts/run_rosetta_only.sh \
    /path/to/template_pdbs \
    /path/to/mpnn_output \
    /path/to/rosetta_output \
    /path/to/ligand.params \
    500

# Step 3: Aggregate scores
python design/instructions/aggregate_scores.py rosetta_output --output scores.csv

# Step 4: Filter designs
python design/instructions/relax_2_filter__allpolar_unsats.py \
    scores.csv \
    rosetta_output \
    filtered_output \
    --target_n 1000 \
    --max_unsat 1

# Step 5: AF3 prep only
python design/scripts/run_design_pipeline.py config.txt --af3-prep-only
```

---

## Quick Decision Tree

```
Do you need to troubleshoot or test specific steps?
│
├─ NO → Use Complete Automated Workflow
│        ├─ Docking: submit_complete_workflow.sh
│        └─ Design: run_design_pipeline.py
│
└─ YES → Use Individual Steps
         ├─ Need to re-run just one docking array? → grade_conformers_*.py <index>
         ├─ MPNN failed? → run_mpnn_only.sh
         ├─ Rosetta failed? → run_rosetta_only.sh
         ├─ Want different filter settings? → relax_2_filter__allpolar_unsats.py
         └─ Need to regenerate AF3 inputs? → run_design_pipeline.py --af3-prep-only
```

---

## Common Use Cases

### Use Case 1: First Time Running Pipeline
**Recommendation:** Start with individual steps to verify each works

```bash
# 1. Test docking on single array task
python docking/scripts/run_docking_workflow.py config.txt --array-index 0 --skip-clustering

# 2. If successful, run full docking
bash docking/scripts/submit_complete_workflow.sh config.txt

# 3. After docking completes, test MPNN on small subset
bash design/scripts/run_mpnn_only.sh clustered_final mpnn_test 5

# 4. If successful, run full design pipeline
python design/scripts/run_design_pipeline.py config.txt
```

### Use Case 2: Production Run (Everything Works)
**Recommendation:** Use complete automated workflow

```bash
# Submit docking
bash docking/scripts/submit_complete_workflow.sh config.txt

# Monitor: squeue -u $USER
# When docking completes, submit design
python design/scripts/run_design_pipeline.py config.txt
```

### Use Case 3: Re-running After Failure
**Recommendation:** Use individual steps to resume

```bash
# Example: Rosetta failed due to memory issue

# 1. Don't re-run MPNN (already done)
# 2. Just run Rosetta with more memory
bash design/scripts/run_rosetta_only.sh \
    clustered_final \
    mpnn_output \
    rosetta_retry \
    ligand.params

# 3. Continue from aggregation
python design/instructions/aggregate_scores.py rosetta_retry --output scores.csv
python design/instructions/relax_2_filter__allpolar_unsats.py ...
python design/scripts/run_design_pipeline.py config.txt --af3-prep-only
```

### Use Case 4: Testing Different Parameters
**Recommendation:** Use individual steps

```bash
# Test different MPNN temperatures
bash design/scripts/run_mpnn_only.sh pdbs mpnn_temp_0.1 10
# Edit script to change temperature, then:
bash design/scripts/run_mpnn_only.sh pdbs mpnn_temp_0.5 10

# Test different filter stringency
python design/instructions/relax_2_filter__allpolar_unsats.py scores.csv rosetta filtered_strict --max_unsat 0
python design/instructions/relax_2_filter__allpolar_unsats.py scores.csv rosetta filtered_loose --max_unsat 2
```

---

## Documentation

- **Complete Workflow Details:** See `INTEGRATED_PIPELINE_SUMMARY.md`
- **Individual Steps Guide:** See `design/INDIVIDUAL_STEPS_GUIDE.md`
- **Docking Workflow:** See `docking/WORKFLOW_README.md`
- **Quick Start:** See `QUICK_START.md`

---

## Summary Table

| Workflow Type | Pros | Cons | When to Use |
|---------------|------|------|-------------|
| **Complete Automated** | ✅ Fully automated<br>✅ Handles dependencies<br>✅ Minimal intervention | ❌ Hard to debug failures<br>❌ All-or-nothing | Production runs<br>Validated pipelines |
| **Individual Steps** | ✅ Fine-grained control<br>✅ Easy to troubleshoot<br>✅ Resume from failures<br>✅ Test parameters | ❌ More manual work<br>❌ Must track state | Development<br>Troubleshooting<br>Parameter testing<br>Partial re-runs |

---

## Getting Help

- Check log files: `*.out` and `*.err` in your working directory
- Monitor jobs: `squeue -u $USER`
- Check disk space: `df -h /scratch/alpine/$USER`
- Review config: Ensure all paths are correct in `config.txt`

For detailed step-by-step instructions, see `design/INDIVIDUAL_STEPS_GUIDE.md`
