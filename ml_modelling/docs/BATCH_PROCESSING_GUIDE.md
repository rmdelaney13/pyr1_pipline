# Batch Processing Guide: Lists of Ligands × Variants → ML Dataset

This guide shows how to take **lists of ligands and sequences** and generate the complete ML dataset for PYR1 biosensor modeling.

## Overview

The batch processing workflow:
```
Lists (ligands + variants)
  ↓
1. Prepare pairs CSV ← YOU ARE HERE
  ↓
2. Run orchestrator (conformers → docking → relax → AF3)
  ↓
3. Aggregate features
  ↓
4. Train ML models
```

---

## Step 1: Prepare Your Input Data

You have **three options** for input:

### Option A: Use Existing Positive Pairs (Recommended for getting started)

You already have 288 validated binder pairs in:
```
ml_modelling/data/ligand_smiles_signature.csv
```

**Columns:**
- `ligand_name`: Ligand name
- `ligand_smiles_or_ligand_ID`: SMILES string
- `PYR1_variant_name`: Variant name
- `PYR1_variant_signature`: Mutations (e.g., "59K;120A;160G")

### Option B: Generate Cartesian Product (All Ligands × All Variants)

**Ligands CSV** (`ligands.csv`):
```csv
ligand_name,ligand_smiles
WIN-55212-2,CCOc1ccc(cc1)C(=O)Nc2ccc(cc2)S(=O)(=O)N
JWH-018,C1=CC=C2C(=C1)C(=CN2CCCCCC)C(=O)C3=CC=CC=C3
diazinon,CCOP(=S)(OCC)OC1=CC(=NC(=C1)C(C)C)C(C)C
```

**Variants CSV** (`variants.csv`):
```csv
variant_name,variant_signature,label,label_tier
PYR1^4F,83K;115Q;120G;159V;160G,1,positive
PYR1_neg_001,59R;81M;120A,0,T1
PYR1_neg_002,59K;81V;120A,0,T2
```

**Note:** `label` (0/1) and `label_tier` are optional in variants CSV. Defaults: `label=1`, `label_tier=positive`

### Option C: Custom Pairs (Specific Combinations Only)

**Pairs CSV** (`custom_pairs.csv`):
```csv
ligand_name,ligand_smiles,variant_name,variant_signature,label,label_tier
WIN-55212-2,SMILES_STRING,PYR1^4F,83K;115Q;120G;159V;160G,1,positive
WIN-55212-2,SMILES_STRING,PYR1_neg_001,59R;81M;120A,0,T1
```

---

## Step 2: Convert to Orchestrator Format

Use the `prepare_pairs_dataset.py` script to create a standardized CSV.

### Example 1: Convert Existing Positives

```bash
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --existing-csv ml_modelling/data/ligand_smiles_signature.csv \
    --output ml_modelling/pairs_dataset_positives.csv \
    --label 1 \
    --label-tier positive \
    --label-source validated_binder
```

**Output:** `pairs_dataset_positives.csv` with columns:
- `pair_id`: Unique ID (auto-generated)
- `ligand_name`, `ligand_smiles`
- `variant_name`, `variant_signature`
- `label` (1), `label_tier` (positive), `label_source` (validated_binder)

### Example 2: Generate Cartesian Product (10 Ligands × 100 Variants = 1000 Pairs)

```bash
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --ligands-csv ligands.csv \
    --variants-csv variants.csv \
    --output ml_modelling/pairs_dataset_cartesian.csv \
    --cartesian
```

**Result:** All possible combinations of ligands and variants.

### Example 3: Merge Positives + Negatives

```bash
# First, prepare positives
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --existing-csv ml_modelling/data/ligand_smiles_signature.csv \
    --output ml_modelling/pairs_positives.csv \
    --label 1 \
    --label-tier positive

# Then, prepare negatives (from custom CSV)
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --pairs-csv negatives_curated.csv \
    --output ml_modelling/pairs_negatives.csv

# Finally, merge them
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --existing-csv ml_modelling/pairs_positives.csv \
    --additional-csv ml_modelling/pairs_negatives.csv \
    --output ml_modelling/pairs_dataset_full.csv \
    --remove-duplicates
```

---

## Step 3: Run the Orchestrator

Once you have a prepared pairs CSV, run the orchestrator to process all pairs:

```bash
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/pairs_dataset_positives.csv \
    --cache-dir /scratch/ml_dataset_cache \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --docking-repeats 50 \
    --use-slurm \
    --max-pairs 10  # Optional: test with 10 pairs first
```

### What the Orchestrator Does (Per Pair):

1. **Generate conformers** (RDKit ETKDG + MMFF)
   - Output: `{cache_dir}/{pair_id}/conformers/conformers_final.sdf`

2. **Thread variant** (PyRosetta)
   - Output: `{cache_dir}/{pair_id}/mutant.pdb`

3. **Dock conformers** (Rosetta with water constraints)
   - Output: `{cache_dir}/{pair_id}/docking/*.pdb`
   - Geometry CSV: `{cache_dir}/{pair_id}/docking/hbond_geometry_summary.csv`

4. **Cluster poses** (RMSD-based)
   - Output: `{cache_dir}/{pair_id}/clustered/*.pdb`
   - Stats: `{cache_dir}/{pair_id}/clustered/clustering_stats.json`

5. **Relax best pose** (FastRelax)
   - Output: `{cache_dir}/{pair_id}/relax/relaxed_pose.pdb`

6. **AF3 predictions** (binary + ternary)
   - Output: `{cache_dir}/{pair_id}/af3_binary/`, `{cache_dir}/{pair_id}/af3_ternary/`

### Caching & Resumability

Each pair gets a `metadata.json` tracking completion:
```json
{
  "pair_id": "WIN-55212-2_PYR1^4F_abc123",
  "stages": {
    "conformers": {"status": "complete"},
    "threading": {"status": "complete"},
    "docking": {"status": "complete"},
    "clustering": {"status": "complete"},
    "relax": {"status": "complete"}
  }
}
```

**Re-running the orchestrator skips completed pairs automatically.**

---

## Step 4: Aggregate Features (TODO - coming soon)

```bash
python ml_modelling/scripts/aggregate_ml_features.py \
    --cache-dir /scratch/ml_dataset_cache \
    --pairs-csv ml_modelling/pairs_dataset_positives.csv \
    --output ml_modelling/features_table.csv
```

**Output columns:**
- Rosetta: `dG_sep`, `buried_unsats`, `sasa_interface`, `hbonds_interface`, `total_score`
- AF3: `ipTM`, `pLDDT_protein`, `pLDDT_ligand`, `interface_PAE`, `ligand_RMSD`
- Conformer: `min_energy`, `rmsd_range`, `conformer_count`
- Labels: `pair_id`, `ligand_name`, `variant_name`, `label`, `label_tier`

---

## Step 5: Train ML Models (TODO - coming soon)

```bash
python ml_modelling/scripts/baseline_models.py \
    --features ml_modelling/features_table.csv \
    --output ml_modelling/baseline_results.csv \
    --splits 0.6 0.2 0.2  # Train/val/test
```

---

## Complete Example: 10 Ligands × 30 Variants = 300 Pairs

### Input Files

**`ligands_pilot.csv`** (10 cannabinoids):
```csv
ligand_name,ligand_smiles
WIN-55212-2,CCOc1ccc(cc1)C(=O)Nc2ccc(cc2)S(=O)(=O)N
JWH-018,C1=CC=C2C(=C1)C(=CN2CCCCCC)C(=O)C3=CC=CC=C3
4F-MDMB-BUTINACA,O=C(N[C@H](C(OC)=O)C(C)(C)C)C1=NN(CCCCF)C2=C1C=CC=C2
JWH-167,C1=CC=C2C(=C1)C(=CN2CCCCCC)C(=O)C3=CC=CC=C3
...
```

**`variants_pilot.csv`** (10 positives + 10 T1 + 10 T2):
```csv
variant_name,variant_signature,label,label_tier,label_source
PYR1^4F,83K;115Q;120G;159V;160G,1,positive,validated
PYR1^JWH,83F;115Q;120G;159L;160G,1,positive,validated
...
PYR1_neg_T1_001,59R;81M;120A,0,T1,historical_screen
PYR1_neg_T1_002,59Q;92M;108V,0,T1,historical_screen
...
PYR1_neg_T2_001,83K;120G,0,T2,near_neighbor
PYR1_neg_T2_002,115Q;160G,0,T2,near_neighbor
...
```

### Workflow

```bash
# 1. Generate cartesian product (10 × 30 = 300 pairs)
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --ligands-csv ligands_pilot.csv \
    --variants-csv variants_pilot.csv \
    --cartesian \
    --output ml_modelling/pilot_300_pairs.csv

# 2. Run pilot (test with 10 pairs first)
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/pilot_300_pairs.csv \
    --cache-dir /scratch/pilot_cache \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --docking-repeats 30 \
    --max-pairs 10

# 3. If successful, run full pilot (300 pairs)
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/pilot_300_pairs.csv \
    --cache-dir /scratch/pilot_cache \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --docking-repeats 50 \
    --use-slurm

# 4. Aggregate features (after all pairs complete)
python ml_modelling/scripts/aggregate_ml_features.py \
    --cache-dir /scratch/pilot_cache \
    --pairs-csv ml_modelling/pilot_300_pairs.csv \
    --output ml_modelling/pilot_features.csv

# 5. Train baseline models
python ml_modelling/scripts/baseline_models.py \
    --features ml_modelling/pilot_features.csv \
    --output ml_modelling/pilot_results.csv
```

---

## Expected Runtimes (Phase 1 Pilot: 300 Pairs)

| Stage | Per Pair | Total (Serial) | Total (Parallel, 300 cores) |
|-------|----------|----------------|----------------------------|
| Conformers | 2 min | 10 hours | 2 min |
| Threading | 1 min | 5 hours | 1 min |
| Docking (50 repeats) | 5 min | 25 hours | 5 min |
| Clustering | 2 min | 10 hours | 2 min |
| Relax | 15 min | 75 hours | 15 min |
| **Total (CPU)** | **~25 min** | **125 hours** | **~25 min** |
| AF3 (8 A100s) | 30 min | 150 hours | ~19 hours |

**Wall time estimate (300 cores + 8 A100s):** ~20 hours

---

## Troubleshooting

### Issue: "Missing required columns"
**Solution:** Make sure your input CSV has the required columns. Use `--existing-csv` for the old format, or `--pairs-csv` for custom format.

### Issue: "Pair already complete, skipping"
**Solution:** This is normal caching behavior. To re-run a pair, delete its cache directory:
```bash
rm -rf /scratch/ml_dataset_cache/{pair_id}
```

### Issue: "Docking failed with severe clash"
**Solution:** Check the geometry diagnostics CSV:
```bash
cat /scratch/ml_dataset_cache/{pair_id}/docking/hbond_geometry_summary.csv
```
If all poses have poor geometry, the ligand may be incompatible with this variant.

### Issue: "AF3 predictions not implemented"
**Solution:** AF3 integration is planned but not yet implemented in the orchestrator. You can run AF3 manually on the relaxed poses.

---

## Quick Reference: Common Use Cases

### Use Case 1: Process Existing 288 Positive Pairs

```bash
# Prepare CSV
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --existing-csv ml_modelling/data/ligand_smiles_signature.csv \
    --output ml_modelling/pairs_288_positives.csv

# Run orchestrator
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/pairs_288_positives.csv \
    --cache-dir /scratch/positives_cache \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --use-slurm
```

### Use Case 2: Test 1 Ligand with 5 New Variants

**`test_variants.csv`:**
```csv
variant_name,variant_signature,label,label_tier
PYR1_test_001,59K;120A;160G,0,T2
PYR1_test_002,59Q;120A;160V,0,T2
PYR1_test_003,81Y;120H;160G,0,T2
PYR1_test_004,92M;122Q;160T,0,T3
PYR1_test_005,83W;110Y;164Y,0,T3
```

**`test_ligand.csv`:**
```csv
ligand_name,ligand_smiles
WIN-55212-2,CCOc1ccc(cc1)C(=O)Nc2ccc(cc2)S(=O)(=O)N
```

```bash
# Generate 1 × 5 = 5 pairs
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --ligands-csv test_ligand.csv \
    --variants-csv test_variants.csv \
    --cartesian \
    --output ml_modelling/test_5_pairs.csv

# Run locally (no SLURM)
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/test_5_pairs.csv \
    --cache-dir /scratch/test_cache \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --docking-repeats 30
```

---

## Summary

✅ **You can now:**
1. Take lists of ligands (SMILES) and variants (signatures)
2. Generate the necessary pairs CSV (cartesian product or custom)
3. Run the orchestrator to dock all pairs and generate features
4. Build ML models to score and predict binding

✅ **Key Scripts:**
- `prepare_pairs_dataset.py` - Convert lists → pairs CSV
- `orchestrate_ml_dataset_pipeline.py` - Run full pipeline on all pairs
- `aggregate_ml_features.py` - Extract features from outputs (TODO)
- `baseline_models.py` - Train ML models (TODO)

✅ **Next Steps:**
- Test with 10 pairs to validate the workflow
- Run pilot with 300 pairs (10 ligands × 30 variants)
- Curate negative dataset (Tier 1, 2, 3)
- Scale to 2000+ pairs for production dataset

---

**End of Batch Processing Guide**
