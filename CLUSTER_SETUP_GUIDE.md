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

# Submit complete docking workflow (3 array tasks as configured)
python /projects/ryde3462/software/pyr1_pipeline/docking/scripts/run_docking_workflow.py config.txt --slurm
```

**What happens:**
1. Creates alignment table (if not already done)
2. Submits array job (tasks 0-2)
3. Automatically submits clustering job (runs after arrays complete)
4. You get job IDs

**Monitor:**
```bash
# Check job status
squeue -u ryde3462

# Watch array job output
tail -f /scratch/alpine/ryde3462/kyna_test/logs/docking_*.out

# Watch clustering output (after arrays finish)
tail -f /scratch/alpine/ryde3462/kyna_test/logs/clustering_*.out
```

**Results will be in:**
```bash
/scratch/alpine/ryde3462/kyna_test/docked/clustered_final/
```

### 6. Run Everything: SDF → Docking → Design → AF3 (Recommended)

The pipeline runs in 3 phases. Phase 1 uses a SLURM orchestrator job so it
survives SSH disconnects. Phases 2 and 3 are quick commands you run from a
login node after each phase completes.

#### Phase 1: Docking → AF3 Prep (~4 hours, submit-and-forget)

This submits a lightweight orchestrator SLURM job (1 CPU, 4GB, 8hr wall time)
that automatically chains: docking arrays → clustering → MPNN → Rosetta →
filter → FASTA → AF3 JSON prep. It stops before submitting AF3 GPU jobs so it
doesn't hold a CPU node while waiting on GPU queue.

```bash
cd /projects/ryde3462/kyna_test
PIPE=/projects/ryde3462/software/pyr1_pipeline

# Full pipeline through AF3 prep (stops before GPU submission)
bash $PIPE/docking/scripts/submit_full_pipeline.sh config.txt \
  --design-args "--skip-af3-submit --skip-af3-analyze"
```

You can log out after this. The orchestrator runs inside SLURM.

**Monitor progress:**
```bash
# Watch the orchestrator output (live)
tail -f /scratch/alpine/ryde3462/kyna_test/logs/pipeline_orchestrator_<JOBID>.out

# Quick status check (stage-by-stage timestamps)
cat /scratch/alpine/ryde3462/kyna_test/logs/pipeline_status.log

# Check all your jobs
squeue -u ryde3462
```

#### Phase 2: Submit AF3 GPU Jobs (seconds, from login node)

Once Phase 1 finishes (check `pipeline_status.log` or `squeue`), submit the
AF3 GPU jobs. This takes seconds — no orchestrator needed:

```bash
cd /projects/ryde3462/kyna_test
python $PIPE/design/scripts/run_design_pipeline.py config.txt --af3-submit-only
```

GPU jobs go into the SLURM queue on their own. Wait time depends on GPU
availability — check `squeue -u ryde3462` periodically.

#### Phase 3: Analyze AF3 Results (after GPU jobs finish)

Once all AF3 GPU jobs have completed:

```bash
cd /projects/ryde3462/kyna_test
python $PIPE/design/scripts/run_design_pipeline.py config.txt --af3-analyze-only
```

This extracts pLDDT, ipTM, ligand RMSD, and filters by quality cutoffs.
Final results land in `/scratch/.../design/af3_analysis/filtered_metrics.csv`.

#### Other Orchestrator Options

```bash
# Design only (docking already done), through AF3 prep
bash $PIPE/docking/scripts/submit_full_pipeline.sh config.txt \
  --design-only --design-args "--skip-af3-submit --skip-af3-analyze"

# Docking only (skip design entirely)
bash $PIPE/docking/scripts/submit_full_pipeline.sh config.txt --docking-only

# Jump from existing Rosetta output to AF3 prep
bash $PIPE/docking/scripts/submit_full_pipeline.sh config.txt \
  --design-only --design-args "--rosetta-to-af3 --skip-af3-submit --skip-af3-analyze"
```

The `--wait` flag tells the pipeline to poll SLURM and wait for each batch of jobs
(MPNN, Rosetta, AF3) to finish before continuing to the next stage.

**What happens automatically (10 stages):**

| Stage | What | Output |
|-------|------|--------|
| 1 | LigandMPNN design on docked PDBs | `mpnn_output/*.fa` |
| 2 | Rosetta relax of MPNN sequences | `rosetta_output/*.pdb` + `.sc` |
| 3 | Aggregate Rosetta scores | `scores/iteration_1_scores.csv` |
| 4 | Filter by Rosetta metrics | `filtered/filtered.csv` + PDBs |
| 5 | Generate FASTA from filtered PDBs | `filtered/filtered.fasta` |
| 6 | Generate AF3 JSON inputs (binary + ternary) | `af3_inputs/binary/*.json`, `af3_inputs/ternary/*.json` |
| 7-8 | Batch JSONs + submit AF3 GPU jobs | `af3_output/binary/`, `af3_output/ternary/` |
| 9 | Analyze AF3 results (pLDDT, ipTM, RMSD) | `af3_analysis/master_metrics.csv` |
| 10 | Filter AF3 results by quality cutoffs | `af3_analysis/filtered_metrics.csv` |

**Monitor jobs during the run:**
```bash
# Check SLURM queue
squeue -u ryde3462

# Watch pipeline output (in another terminal)
# Pipeline prints progress to stdout as each stage completes
```

#### Running Stages Individually

If you want to run stages separately (e.g., after a failure or to re-run a step):

```bash
cd /projects/ryde3462/kyna_test

# Run MPNN + Rosetta + filter, stop before AF3 submission
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt --skip-af3-submit --skip-af3-analyze --wait

# Skip MPNN/Rosetta, jump from existing Rosetta output to AF3 prep + submission
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt --rosetta-to-af3 --wait

# Only prepare AF3 JSONs (no submission)
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt --af3-prep-only

# Only batch and submit AF3 GPU jobs (JSONs already prepared)
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt --af3-submit-only

# Only run AF3 analysis + filtering (after GPU jobs complete)
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt --af3-analyze-only

# Dry run - generate all scripts but don't submit to SLURM
python /projects/ryde3462/software/pyr1_pipeline/design/scripts/run_design_pipeline.py config.txt --dry-run
```

#### All Design Pipeline Flags

| Flag | Effect |
|------|--------|
| `--wait` | Wait for each SLURM job batch to finish before continuing |
| `--dry-run` | Generate scripts but don't submit to SLURM |
| `--iteration N` | Run only iteration N |
| `--skip-mpnn` | Skip LigandMPNN stage |
| `--skip-rosetta` | Skip Rosetta relax stage |
| `--skip-aggregate` | Skip score aggregation |
| `--skip-filter` | Skip Rosetta filtering |
| `--skip-af3-prep` | Skip AF3 JSON generation |
| `--skip-af3-submit` | Stop before submitting AF3 GPU jobs |
| `--skip-af3-analyze` | Skip AF3 analysis and filtering |
| `--af3-prep-only` | Only do FASTA + AF3 JSON prep (skip MPNN/Rosetta) |
| `--af3-submit-only` | Only batch and submit AF3 GPU jobs |
| `--af3-analyze-only` | Only run AF3 analysis + quality filtering |
| `--rosetta-to-af3` | Skip MPNN/Rosetta, aggregate → filter → AF3 |

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

## 📊 Expected Output Structure

After running the full pipeline on kyna_test, you should see:

```
kyna_test/                              (CAMPAIGN_ROOT = /projects/ryde3462/kyna_test)
├── config.txt
├── kyna_ligands.csv
├── kyna_ligands.pkl
├── conformers/
│   ├── kyna_conf.sdf
│   └── 0/
│       ├── 0.pdb
│       ├── 0.params
│       └── 0.sdf

/scratch/alpine/ryde3462/kyna_test/     (SCRATCH_ROOT)
├── docked/
│   ├── array_0/                        (can delete after clustering)
│   ├── array_1/
│   └── clustered_final/               ← Input to design pipeline
│       ├── *.pdb
│       └── cluster_summary.csv
│
└── design/
    ├── iteration_1/
    │   ├── mpnn_output/
    │   │   ├── submit_mpnn.sh          (auto-generated SLURM script)
    │   │   └── array001_mpnn/*.fa
    │   │
    │   ├── rosetta_output/
    │   │   ├── submit_rosetta.sh       (auto-generated SLURM script)
    │   │   ├── *.pdb                   (relaxed structures)
    │   │   └── *.sc                    (score files)
    │   │
    │   ├── scores/
    │   │   └── iteration_1_scores.csv
    │   │
    │   └── filtered/
    │       ├── filtered.csv            (top N designs)
    │       ├── filtered.fasta          (sequences for AF3)
    │       └── *.pdb
    │
    ├── af3_inputs/
    │   ├── binary/
    │   │   ├── *.json                  (one per design)
    │   │   ├── batch_01/               (batched for GPU submission)
    │   │   └── batch_02/
    │   └── ternary/
    │       ├── *.json
    │       ├── batch_01/
    │       └── batch_02/
    │
    ├── af3_output/
    │   ├── binary/
    │   │   ├── submit_af3_binary.sh    (auto-generated SLURM script)
    │   │   ├── logs/
    │   │   └── */                      (AF3 prediction outputs)
    │   └── ternary/
    │       ├── submit_af3_ternary.sh
    │       ├── logs/
    │       └── */
    │
    └── af3_analysis/
        ├── binary_metrics.csv
        ├── ternary_metrics.csv
        ├── master_metrics.csv          ← Combined with ligand RMSD
        └── filtered_metrics.csv        ← Final quality-filtered results
```

---

## ✅ Success Checklist

After each step, verify:

- [ ] **Step 1 (create_table):** `kyna_ligands.csv` exists and has rows
- [ ] **Step 2 (docking):** `/scratch/.../docked/clustered_final/` has PDB files
- [ ] **Step 3 (MPNN):** `/scratch/.../design/iteration_1/mpnn_output/` has `.fa` files
- [ ] **Step 4 (Rosetta):** `/scratch/.../design/iteration_1/rosetta_output/` has `.pdb` and `.sc` files
- [ ] **Step 5 (scores):** `/scratch/.../design/iteration_1/scores/iteration_1_scores.csv` exists
- [ ] **Step 6 (filtering):** `/scratch/.../design/iteration_1/filtered/` has filtered PDBs + `filtered.csv`
- [ ] **Step 7 (FASTA):** `/scratch/.../design/iteration_1/filtered/filtered.fasta` exists
- [ ] **Step 8 (AF3 prep):** `/scratch/.../design/af3_inputs/binary/` and `ternary/` have `.json` files
- [ ] **Step 9 (AF3 submit):** Jobs appear in `squeue -u ryde3462` (or already completed)
- [ ] **Step 10 (AF3 analysis):** `/scratch/.../design/af3_analysis/master_metrics.csv` exists
- [ ] **Step 11 (AF3 filter):** `/scratch/.../design/af3_analysis/filtered_metrics.csv` has passing designs

---

## 🎯 Quick Command Reference

```bash
cd /projects/ryde3462/kyna_test
PIPE=/projects/ryde3462/software/pyr1_pipeline

# --- PHASE 1: Docking → AF3 Prep (submit-and-forget, ~4 hours) ---

bash $PIPE/docking/scripts/submit_full_pipeline.sh config.txt \
  --design-args "--skip-af3-submit --skip-af3-analyze"

# Design only (docking already done)
bash $PIPE/docking/scripts/submit_full_pipeline.sh config.txt \
  --design-only --design-args "--skip-af3-submit --skip-af3-analyze"

# Docking only
bash $PIPE/docking/scripts/submit_full_pipeline.sh config.txt --docking-only

# --- PHASE 2: Submit AF3 GPU Jobs (from login node, takes seconds) ---

python $PIPE/design/scripts/run_design_pipeline.py config.txt --af3-submit-only

# --- PHASE 3: Analyze AF3 Results (after GPU jobs finish) ---

python $PIPE/design/scripts/run_design_pipeline.py config.txt --af3-analyze-only

# --- MONITORING ---

tail -f /scratch/alpine/ryde3462/kyna_test/logs/pipeline_orchestrator_*.out
cat /scratch/alpine/ryde3462/kyna_test/logs/pipeline_status.log
squeue -u ryde3462

# --- DOCKING ONLY ---

# Submit docking via SLURM (returns immediately)
python $PIPE/docking/scripts/run_docking_workflow.py config.txt --slurm

# Submit docking via SLURM and wait for completion
python $PIPE/docking/scripts/run_docking_workflow.py config.txt --slurm --wait

# Test locally (no SLURM)
python $PIPE/docking/scripts/run_docking_workflow.py config.txt --local-arrays 2

# --- DESIGN PIPELINE ONLY (after docking completes) ---

# Run everything: MPNN → Rosetta → Filter → AF3 prep → AF3 submit → AF3 analyze
python $PIPE/design/scripts/run_design_pipeline.py config.txt --wait

# Dry run (see what would happen without submitting jobs)
python $PIPE/design/scripts/run_design_pipeline.py config.txt --dry-run

# --- PARTIAL DESIGN RUNS ---

# MPNN + Rosetta + filter only (stop before AF3)
python $PIPE/design/scripts/run_design_pipeline.py config.txt --skip-af3-submit --skip-af3-analyze --wait

# Jump from existing Rosetta output → AF3
python $PIPE/design/scripts/run_design_pipeline.py config.txt --rosetta-to-af3 --wait

# Submit AF3 jobs only (JSONs already prepared)
python $PIPE/design/scripts/run_design_pipeline.py config.txt --af3-submit-only

# Analyze AF3 results only (after GPU jobs complete)
python $PIPE/design/scripts/run_design_pipeline.py config.txt --af3-analyze-only

# --- MONITORING ---
squeue -u ryde3462
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
