# Cluster Setup & Testing Guide

This guide shows you how to set up the pipeline on the cluster and test with kyna_test.

## 📁 Directory Structure on Cluster

```
/projects/ryde3462/
├── software/
│   └── pyr1_pipeline/              # ← Your GitHub clone (pipeline code)
│       ├── docking/
│       ├── design/
│       └── templates/
│
├── kyna_test/                      # ← Your test project
│   ├── config.txt                  # ← Copy from pipeline/kyna_test_config.txt
│   ├── conformers/
│   │   └── kyna_conf.sdf          # ← Already there
│   └── (outputs will be created)
│
└── (other projects...)

/scratch/alpine/ryde3462/
└── kyna_test/                      # ← Large outputs go here
    ├── docked/
    └── design/
```

---

## 🚀 Step-by-Step Setup

### 1. Copy Config to Your Project Directory

```bash
# On the cluster
cd /projects/ryde3462/kyna_test

# Copy the config from the pipeline repo
cp /projects/ryde3462/software/pyr1_pipeline/kyna_test_config.txt config.txt

# Verify it's there
ls -lh config.txt
```

### 2. Verify Your SDF File

```bash
cd /projects/ryde3462/kyna_test/conformers
ls -lh kyna_conf.sdf

# Should see your kyna_conf.sdf file
```

### 3. Test Step 1: Create Table (Table Generation)

This generates params and PDB files from your SDF.

```bash
cd /projects/ryde3462/kyna_test

# Run create_table.py
python /projects/ryde3462/software/pyr1_pipeline/docking/scripts/create_table.py config.txt
```

**Expected output:**
- `kyna_ligands.csv` - Alignment table
- `kyna_ligands.pkl` - Pickle file
- `conformers/0/` - Directory with:
  - `0.pdb` - Ligand structure
  - `0.params` - Rosetta params
  - `0.sdf` - Individual conformer

**Verify:**
```bash
# Check table was created
cat kyna_ligands.csv

# Check conformer files
ls -lh conformers/0/
```

### 4. Test Step 2: Run Single Docking Task (Local Test)

Test docking locally before submitting SLURM jobs.

```bash
cd /projects/ryde3462/kyna_test

# Run array index 0 locally (no SLURM)
python /projects/ryde3462/software/pyr1_pipeline/docking/scripts/grade_conformers_glycine_shaved.py config.txt 0
```

**Expected output:**
- `/scratch/alpine/ryde3462/kyna_test/docked/array_0/` - Docked PDBs
- Log messages showing docking progress

**This will take a few minutes.** If it works, you're ready for SLURM!

### 5. Test Step 3: Submit Full Docking via SLURM

Now submit the real workflow with SLURM arrays.

```bash
cd /projects/ryde3462/kyna_test

# Submit complete docking workflow (2 array tasks as configured)
bash /projects/ryde3462/software/pyr1_pipeline/docking/scripts/submit_complete_workflow.sh config.txt
```

**What happens:**
1. Submits array job (tasks 0-1)
2. Automatically submits clustering job (runs after arrays complete)
3. You get job IDs

**Monitor:**
```bash
# Check job status
squeue -u ryde3462

# Watch array job output
tail -f docking_*.out

# Watch clustering output (after arrays finish)
tail -f clustering_*.out
```

**Results will be in:**
```bash
/scratch/alpine/ryde3462/kyna_test/docked/clustered_final/
```

### 6. Test Step 4: Run Design Pipeline

After docking completes, run the design pipeline.

```bash
cd /projects/ryde3462/kyna_test

# Full design pipeline
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt
```

**Or test individual design steps:**

#### Test MPNN Only:
```bash
bash /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_mpnn_only.sh \
    /scratch/alpine/ryde3462/kyna_test/docked/clustered_final \
    /scratch/alpine/ryde3462/kyna_test/design/mpnn_output \
    5
```

#### Test Rosetta Only (after MPNN completes):
```bash
bash /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_rosetta_only.sh \
    /scratch/alpine/ryde3462/kyna_test/docked/clustered_final \
    /scratch/alpine/ryde3462/kyna_test/design/mpnn_output \
    /scratch/alpine/ryde3462/kyna_test/design/rosetta_output \
    /projects/ryde3462/kyna_test/conformers/0/0.params \
    50
```

---

## 🔧 Troubleshooting

### Issue: "Module not found" errors

**Solution:** Activate the correct conda environment

```bash
# For docking (needs PyRosetta)
module load anaconda
conda activate pyrosetta

# For MPNN
conda activate ligandmpnn_env
```

### Issue: "File not found" errors

**Check these paths in your config.txt:**

```bash
# Verify pipeline root exists
ls /projects/ryde3462/software/pyr1_pipeline

# Verify template PDB files exist
ls /projects/ryde3462/software/pyr1_pipeline/docking/ligand_alignment/files_for_PYR1_docking/

# Verify your SDF exists
ls /projects/ryde3462/kyna_test/conformers/kyna_conf.sdf
```

### Issue: SLURM jobs fail immediately

**Check logs:**
```bash
# Find your job ID
squeue -u ryde3462

# Check error log
cat docking_JOBID_0.err
cat docking_JOBID_0.out
```

### Issue: Out of disk space

**Check scratch usage:**
```bash
df -h /scratch/alpine/ryde3462
```

If full, clean up old outputs:
```bash
rm -rf /scratch/alpine/ryde3462/kyna_test/docked/array_*
# (Keep clustered_final though!)
```

---

## 📊 Expected File Sizes

After running on kyna_test, you should see approximately:

```
kyna_test/
├── config.txt                      (~8 KB)
├── kyna_ligands.csv               (~1-10 KB)
├── conformers/
│   ├── kyna_conf.sdf              (~varies)
│   └── 0/
│       ├── 0.pdb                  (~5 KB)
│       ├── 0.params               (~5 KB)
│       └── 0.sdf                  (~5 KB)

/scratch/alpine/ryde3462/kyna_test/
├── docked/
│   ├── array_0/                   (~100s of MB, can delete after clustering)
│   ├── array_1/                   (~100s of MB, can delete after clustering)
│   └── clustered_final/           (~10-50 MB) ← KEEP THIS
│       ├── *.pdb                  (clustered docked poses)
│       └── cluster_summary.csv
│
└── design/
    ├── mpnn_output/               (~100s of MB)
    ├── rosetta_output/            (~GBs)
    └── filtered/                  (~100s of MB) ← Final designs
```

---

## ✅ Success Checklist

After each step, verify:

- [ ] **Step 1 (create_table):** `kyna_ligands.csv` exists and has rows
- [ ] **Step 2 (docking):** `/scratch/.../docked/clustered_final/` has PDB files
- [ ] **Step 3 (MPNN):** `/scratch/.../design/mpnn_output/` has `.fa` files
- [ ] **Step 4 (Rosetta):** `/scratch/.../design/rosetta_output/` has `.pdb` files
- [ ] **Step 5 (filtering):** `/scratch/.../design/filtered/` has filtered PDBs
- [ ] **Step 6 (AF3 prep):** `/scratch/.../design/af3_inputs/` has `.json` files

---

## 🎯 Quick Command Reference

```bash
# Test locally (no SLURM)
cd /projects/ryde3462/kyna_test
python /projects/ryde3462/software/pyr1_pipeline/docking/scripts/create_table.py config.txt
python /projects/ryde3462/software/pyr1_pipeline/docking/scripts/grade_conformers_glycine_shaved.py config.txt 0

# Submit full workflow via SLURM
bash /projects/ryde3462/software/pyr1_pipeline/docking/scripts/submit_complete_workflow.sh config.txt

# Monitor jobs
squeue -u ryde3462
tail -f *.out

# Run design after docking
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt
```

---

## 📖 Next Steps

Once kyna_test works:

1. **Create a new project directory** for your real ligand
2. **Copy the config template** and customize paths
3. **Run the full pipeline** with confidence!

Example for a new project:
```bash
mkdir /projects/ryde3462/my_new_ligand
cd /projects/ryde3462/my_new_ligand
mkdir conformers
cp /projects/ryde3462/software/pyr1_pipeline/templates/unified_config_template.txt config.txt
nano config.txt  # Edit CAMPAIGN_ROOT, SCRATCH_ROOT, etc.
```

Happy docking! 🚀
