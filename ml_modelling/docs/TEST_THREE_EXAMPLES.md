# Testing ML Pipeline with 3 Example Pairs

This guide walks you through testing the ML docking pipeline with three specific examples:
1. **Win wt** - abscisic Acid with wild-type PYR1
2. **PYR1nitav2** - nitazene with PYR1^nitav2 variant
3. **lca seq16** - Lithocholic Acid with seq16_designed variant

## Quick Start (On Your Cluster)

### 1. Verify Test Data

The test CSV file has been created at:
```bash
ml_modelling/data/test_three_pairs.csv
```

View the file:
```bash
cat ml_modelling/data/test_three_pairs.csv
```

Expected content:
- 3 rows (plus header)
- Columns: pair_id, ligand_name, ligand_smiles, variant_name, variant_signature, label, label_tier

### 2. Update Cluster Paths

Edit the test script to match your cluster setup:
```bash
nano ml_modelling/scripts/test_three_pairs.sh
```

**Update these variables:**
```bash
PROJECT_ROOT="/path/to/pyr1_pipeline"          # Your pyr1_pipeline location
CACHE_DIR="${SCRATCH}/ml_dataset_test"         # Output directory (use $SCRATCH or fast storage)
TEMPLATE_PDB="${PROJECT_ROOT}/templates/PYR1_WT_clean.pdb"  # WT PYR1 template PDB
```

### 3. Run the Test

**Option A: Interactive (Local) Run**
Good for debugging, runs on login node:
```bash
cd /path/to/pyr1_pipeline
bash ml_modelling/scripts/test_three_pairs.sh
```

**Option B: SLURM Submission** (Recommended)
Edit `test_three_pairs.sh` to uncomment the `--use-slurm` version, then:
```bash
bash ml_modelling/scripts/test_three_pairs.sh
```

Monitor jobs:
```bash
squeue -u $USER
watch -n 10 squeue -u $USER  # Auto-refresh every 10 seconds
```

### 4. Check Outputs

After the pipeline runs, check the cache directory:

```bash
CACHE_DIR="${SCRATCH}/ml_dataset_test"  # Or whatever you set above

# View processing summary
cat $CACHE_DIR/processing_summary.csv

# List all pair directories
ls -lh $CACHE_DIR/

# Check a specific pair's outputs
PAIR_ID="test_001"  # Or test_002, test_003
ls -lhR $CACHE_DIR/$PAIR_ID/
```

Expected directory structure per pair:
```
$CACHE_DIR/
├── test_001/
│   ├── metadata.json           # Stage completion tracking
│   ├── conformers/
│   │   └── conformers_final.sdf    # 10 ligand conformers
│   ├── mutant.pdb              # Threaded mutant structure
│   ├── docking/
│   │   ├── docking_0.log
│   │   └── docked_poses/       # Raw docking outputs
│   ├── clustered/
│   │   ├── clustering_stats.json
│   │   ├── clustering_stats.csv
│   │   └── best_pose.pdb
│   ├── relax/                  # (if implemented)
│   └── af3_binary/             # (if implemented)
├── test_002/
│   └── ...
├── test_003/
│   └── ...
└── processing_summary.csv
```

## Pipeline Stages

The ML pipeline runs these stages for each pair:

### Stage 1: Conformer Generation
- Uses RDKit ETKDG to generate 30 conformers
- MMFF energy minimization
- Clusters to 10 final diverse conformers
- **Output:** `conformers/conformers_final.sdf`

### Stage 2: Mutation Threading
- Threads variant mutations onto WT PYR1 template
- Uses PyRosetta for side-chain placement
- **Output:** `mutant.pdb`

### Stage 3: Docking to Mutant
- Docks all 10 conformers to the mutant pocket
- 50 docking attempts per conformer (500 total poses)
- Direct docking to mutant (no glycine shaving!)
- **Output:** `docking/docked_poses/`

### Stage 4: Clustering & Statistics
- Clusters 500 poses by RMSD (2.0 Å cutoff)
- Calculates convergence ratio (top cluster size / total)
- Detects clashes and binding failures
- **Output:** `clustered/clustering_stats.json`

### Stage 5: Rosetta Relax
- Relaxes best pose with Rosetta
- Skipped if severe clash detected (score > 5)
- **Output:** `relax/relaxed.pdb` (NOT YET IMPLEMENTED)

### Stage 6: AF3 Predictions
- Binary: PYR1 + ligand
- Ternary: PYR1 + ABA + ligand
- Uses Rosetta-relaxed pose as template
- **Output:** `af3_binary/` and `af3_ternary/` (NOT YET IMPLEMENTED)

## Interpreting Results

### Check Metadata
```bash
cat $CACHE_DIR/test_001/metadata.json
```

Look for:
- `"status": "complete"` for each stage
- No `"errors"` entries

### Check Clustering Statistics
```bash
cat $CACHE_DIR/test_001/clustered/clustering_stats.json | python -m json.tool
```

**Key metrics:**
- **convergence_ratio**: 0.5+ is good (50%+ poses in top cluster)
- **best_overall_score**: <0 is good, >5 suggests clash
- **top_cluster_size**: Larger = more consistent binding
- **clash_count**: Should be 0

Example good result:
```json
{
  "global_stats": {
    "convergence_ratio": 0.68,
    "best_overall_score": -2.3,
    "top_cluster_size": 340,
    "total_poses": 500,
    "num_clusters": 8,
    "clash_count": 0
  }
}
```

Example bad result (clash):
```json
{
  "global_stats": {
    "convergence_ratio": 0.12,
    "best_overall_score": 8.7,
    "top_cluster_size": 60,
    "total_poses": 500,
    "num_clusters": 45,
    "clash_count": 12
  }
}
```

### Visualize Docked Poses

Download the best pose and view in PyMOL:
```bash
# On cluster
scp user@cluster:$CACHE_DIR/test_001/clustered/best_pose.pdb .

# On local machine
pymol best_pose.pdb
```

In PyMOL:
```python
# Show binding pocket residues
select pocket, resi 59+81+83+92+94+108+110+117+120+122+141+159+160+163+164+167
show sticks, pocket
show spheres, organic  # Ligand
zoom organic
```

## Troubleshooting

### Pipeline Fails at Conformer Generation
**Error:** `ModuleNotFoundError: No module named 'ligand_conformers'`

**Fix:** Make sure you're running from the project root:
```bash
cd /path/to/pyr1_pipeline
python -m ligand_conformers --help  # Should work
```

### Pipeline Fails at Threading
**Error:** `FileNotFoundError: template PDB not found`

**Fix:** Update `TEMPLATE_PDB` path in `test_three_pairs.sh`
```bash
# Find your template PDB
find . -name "*PYR1*clean*.pdb"
# Update TEMPLATE_PDB in test_three_pairs.sh
```

### SLURM Jobs Fail Immediately
**Error:** SLURM jobs exit with non-zero status

**Fix:** Check the log files:
```bash
cat $CACHE_DIR/test_001/docking_0.log
# Look for Python errors, missing modules, path issues
```

Common issues:
- Wrong conda environment (missing PyRosetta, RDKit)
- Module load errors (need to load alphafold3, rosetta, etc.)
- Path issues (relative vs absolute paths)

### Empty Docking Directory
**Error:** `docking/` directory exists but is empty

**Fix:** Check if docking script path is correct:
```bash
ls pyr1_pipeline/docking/scripts/grade_conformers_mutant_docking.py
# If missing, update path in orchestrate_ml_dataset_pipeline.py
```

## Manual Run (Alternative)

If the orchestrator script doesn't work, you can run stages manually:

### Step 1: Generate Conformers
```bash
python -m ligand_conformers \
    --input "CC1=CC(=O)CC(C1/C=C/C(=C/C(=O)O)/C)(C)C" \
    --input-type smiles \
    --output $CACHE_DIR/test_001/conformers \
    --num-confs 30 \
    --k-final 10 \
    --cluster-rmsd-cutoff 1.0
```

### Step 2: Thread Mutations
```bash
python scripts/thread_variant_to_pdb.py \
    --template templates/PYR1_WT_clean.pdb \
    --signature "" \
    --output $CACHE_DIR/test_001/mutant.pdb \
    --chain A
```

### Step 3: Run Docking
```bash
python docking/scripts/grade_conformers_mutant_docking.py \
    --mutant-pdb $CACHE_DIR/test_001/mutant.pdb \
    --ligand-sdf $CACHE_DIR/test_001/conformers/conformers_final.sdf \
    --output-dir $CACHE_DIR/test_001/docking \
    --repeats 50
```

### Step 4: Cluster Poses
```bash
python docking/scripts/cluster_docked_with_stats.py \
    --input-dir $CACHE_DIR/test_001/docking \
    --output-dir $CACHE_DIR/test_001/clustered \
    --rmsd-cutoff 2.0 \
    --stats-csv clustering_stats.csv
```

## Next Steps

After successfully running the 3 test pairs:

1. **Verify all 3 pairs completed successfully**
   - Check `processing_summary.csv`
   - Inspect clustering statistics for each pair

2. **Visualize binding modes**
   - Download `best_pose.pdb` files
   - Compare binding modes in PyMOL

3. **Run on larger dataset**
   - Create a CSV with more pairs (e.g., 10-30 pairs)
   - Test with SLURM submission
   - Monitor resource usage (time, memory)

4. **Implement remaining stages**
   - Add Rosetta relax (Stage 5)
   - Add AF3 predictions (Stage 6)
   - Test feature aggregation

5. **Generate ML features**
   - Run `aggregate_ml_features.py` when stages 5-6 are implemented
   - Train baseline ML models

## Expected Runtime

For 3 pairs (local run, no SLURM):
- Conformer generation: ~2-5 min per pair
- Threading: <1 min per pair
- Docking: ~10-20 min per pair (50 repeats)
- Clustering: ~1 min per pair
- **Total: ~30-60 minutes**

For 3 pairs (SLURM, parallel):
- **Total wall time: ~15-30 minutes**

## Questions?

See full documentation:
- [README.md](../README.md) - Main ML dataset overview
- [IMPLEMENTATION_GUIDE.md](IMPLEMENTATION_GUIDE.md) - Complete implementation guide
- [../docking/WORKFLOW_README.md](../../docking/WORKFLOW_README.md) - Docking details

---

**Status:** Test ready - waiting for cluster run
**Last Updated:** 2026-02-16
