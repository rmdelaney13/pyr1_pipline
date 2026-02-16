# PYR1 ML Dataset Generation: Complete End-to-End Workflow

**📖 Related Guides:**
- **[Batch Processing Guide](BATCH_PROCESSING_GUIDE.md)** - Process lists of ligands × variants
- **[Project Plan](PYR1_ML_DATASET_PROJECT_PLAN.md)** - Full project scope and timeline

## Overview

This document describes the complete workflow for generating ML datasets from ligand SMILES and PYR1 variant sequences to water-constrained docked poses.

**Pipeline Stages**:
0. **Prepare Pairs Dataset** (lists → standardized pairs CSV) - **NEW for batch processing**
1. **Ligand Conformer Generation** (SMILES → 3D conformers)
2. **Alignment Table Creation** (auto-detect H-bond acceptor atoms)
3. **Variant Threading** (sequence → mutant PDB)
4. **Water-Constrained Docking** (conformers + mutant → high-quality poses)
5. **Feature Extraction** (poses → ML features)

---

## Prerequisites

### Required Software

```bash
# Python packages
pip install rdkit pandas numpy biopython

# PyRosetta (for docking and threading)
# Follow installation instructions at: https://www.pyrosetta.org/

# AlphaFold3 (optional, for AF3 predictions)
# See AF3 documentation for setup
```

### Required Files

- **Template PDB**: WT PYR1 structure with crystallographic waters
  - Example: `3QN1_nolig_H2O.pdb`
- **Ligand Library**: CSV with SMILES and metadata
  - Columns: `ligand_name`, `smiles`, `rotatable_bonds` (optional)
- **Variant Library**: CSV with variant signatures
  - Columns: `variant_name`, `signature` (e.g., "59K;120A;160G")

---

## Stage 0: Prepare Pairs Dataset (Batch Processing)

**⚠️ NEW:** For batch processing multiple ligands and variants, see the [Batch Processing Guide](BATCH_PROCESSING_GUIDE.md).

### Input
- List of ligands (CSV with `ligand_name`, `ligand_smiles`)
- List of variants (CSV with `variant_name`, `variant_signature`)
- OR: Existing pairs CSV (e.g., `ml_modelling/data/ligand_smiles_signature.csv`)

### Output
- `pairs_dataset.csv`: Standardized pairs for orchestrator
  - Columns: `pair_id`, `ligand_name`, `ligand_smiles`, `variant_name`, `variant_signature`, `label`, `label_tier`, `label_source`

### Command

```bash
# Option 1: Convert existing positives (288 pairs)
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --existing-csv ml_modelling/data/ligand_smiles_signature.csv \
    --output ml_modelling/pairs_dataset.csv

# Option 2: Generate cartesian product (all ligands × all variants)
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --ligands-csv ligands.csv \
    --variants-csv variants.csv \
    --cartesian \
    --output ml_modelling/pairs_dataset.csv

# Option 3: Merge positives + negatives
python ml_modelling/scripts/prepare_pairs_dataset.py \
    --existing-csv positives.csv \
    --additional-csv negatives.csv \
    --output ml_modelling/pairs_dataset.csv
```

### Batch Orchestration

Once you have a pairs CSV, run the orchestrator to process all pairs automatically:

```bash
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/pairs_dataset.csv \
    --cache-dir /scratch/ml_dataset_cache \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --docking-repeats 50 \
    --use-slurm \
    --max-pairs 10  # Test with 10 pairs first
```

**The orchestrator handles Stages 1-5 automatically for each pair with caching and resumability.**

---

## Stage 1: Ligand Conformer Generation

### Input
- SMILES string

### Output
- `conformers_final.sdf`: Clustered conformers (typically 10-50 per ligand)
- `conformers_all.sdf`: All generated conformers (before clustering)
- Metadata JSON with RMSD matrix

### Command

```bash
# Single ligand
python -m ligand_conformers \
    --input "CCOc1ccc(cc1)C(=O)Nc2ccc(cc2)S(=O)(=O)N" \
    --input-type smiles \
    --output campaign_root/conformers_deprot/WIN_55212 \
    --num-confs 150 \
    --k-final 10 \
    --cluster-rmsd-cutoff 1.0

# Batch processing from CSV
python ml_modelling/scripts/process_ligand_smiles.py \
    --input ligands.csv \
    --output-dir campaign_root/conformers_deprot
```

### Key Parameters

- `--num-confs`: Number of initial conformers to generate (3× final)
- `--k-final`: Number of conformers to keep after clustering
- `--cluster-rmsd-cutoff`: RMSD threshold for clustering (Å)
- `--refine`: Enable OpenMM energy minimization (slower, higher quality)

### Expected Runtime
- ~1-5 minutes per ligand (depends on size and rotatable bonds)
- Rigid molecules (<5 rotatable bonds): 50 conformers
- Flexible molecules (>10 rotatable bonds): 200+ conformers

---

## Stage 2: Alignment Table Creation

### Input
- Conformer SDF files (from Stage 1)
- Template PDB with target residue

### Output
- `ligands_alignment.csv`: Table with alignment atom triplets
  - Columns: `Molecule Name`, `Molecule ID`, `Conformer Range`, `Molecule Atoms`, `Target Atoms`

### Command

```bash
python docking/scripts/create_table.py ml_modelling/config/mutant_docking_example.conf
```

### Configuration (in config file)

```ini
[create_table]
MoleculeSDFs = campaign_root/conformers_deprot/*/*.sdf
CSVFileName = campaign_root/ligands_alignment.csv

# Auto-generate alignment triplets from ligand acceptor atoms
DynamicAcceptorAlignment = True
MaxDynamicAlignments = 20
IncludeReverseNeighborOrder = False

# Target atoms in template ligand (water acceptor + neighbors)
TargetAtomTriplets = O2-C11-C9

# Auto-detect acceptors using RDKit
AcceptorMode = auto
```

### What It Does

1. Scans each ligand for H-bond acceptor atoms (O, N)
2. For each acceptor, finds 2 neighbors to define orientation plane
3. Creates alignment triplet: `acceptor-neighbor1-neighbor2`
4. Maps SDF atom indices → Rosetta PDB atom names
5. Writes CSV table with all valid alignment combinations

### Expected Output

```csv
Molecule Name,Molecule ID,Conformer Range,Molecule Atoms,Target Atoms
WIN-55212,WIN_55212,1_10_,O15-C14-C13,O2-C11-C9
WIN-55212,WIN_55212,1_10_,O15-C13-C14,O2-C11-C9
WIN-55212,WIN_55212,1_10_,N19-C18-C20,O2-C11-C9
```

---

## Stage 3: Variant Threading

### Input
- Template PDB (WT structure)
- Variant signature (e.g., "59K;120A;160G")

### Output
- Mutant PDB with threaded mutations

### Command

```bash
# Single variant
python scripts/thread_variant_to_pdb.py \
    --template docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --signature "59K;120A;160G" \
    --output campaign_root/threaded_variants/mutant_59K_120A_160G.pdb \
    --chain A

# Batch processing from CSV
python scripts/thread_variant_to_pdb.py \
    --template 3QN1_nolig_H2O.pdb \
    --csv variants.csv \
    --output-dir campaign_root/threaded_variants/
```

### Variant Signature Formats

All formats are supported:
- Semicolon: `59K;120A;160G` (position + target AA)
- Underscore: `K59Q_Y120A_A160G` (WT + position + target AA)
- Space-separated: `59K 120A 160G`

### What It Does

1. Loads template PDB into PyRosetta
2. Parses variant signature → list of mutations
3. Applies mutations using `mutate_residue()` (preserves backbone)
4. Packs sidechains around mutations
5. Writes mutant PDB

### Expected Runtime
- ~30 seconds per variant
- Batch processing: parallel via array jobs

---

## Stage 4: Water-Constrained Docking

### Input
- Mutant PDB (from Stage 3)
- Conformer alignment table CSV (from Stage 2)
- Conformer files (from Stage 1)

### Output
- Docked pose PDBs (`rep_1.pdb`, `rep_2.pdb`, ...)
- Geometry CSV (`hbond_geometry_summary.csv`)

### Command

```bash
# Single job (all conformers)
python docking/scripts/grade_conformers_mutant_docking.py \
    ml_modelling/config/mutant_docking_example.conf

# SLURM array job (parallel processing)
#SBATCH --array=0-9
python docking/scripts/grade_conformers_mutant_docking.py \
    ml_modelling/config/mutant_docking_example.conf \
    $SLURM_ARRAY_TASK_ID
```

### Configuration (in config file)

```ini
[mutant_docking]
MutantPDB = campaign_root/threaded_variants/mutant_59K_120A_160G.pdb
OutputDir = scratch/docked_mutant_59K_120A_160G

# H-bond filter (CRITICAL for ML quality)
EnableHBondGeometryFilter = True
UseClosestLigandAcceptor = True
HBondDistanceIdeal = 2.8
HBondDistanceIdealBuffer = 0.8
HBondConstraintWeight = 4.0
EnforceFinalIdealGeometry = True

# Sampling
MaxPerturbTries = 30
Rotation = 25.0
Translation = 0.5

# Clustering
EnablePoseClusteringInArrayTask = False  # Disable for array jobs
ClusterRMSDCutoff = 0.75
```

### What It Does

1. **Load mutant structure** (with real sidechains)
2. **For each conformer**:
   - Align to pocket using SVD (triplet atoms)
   - Graft into mutant pose
3. **Perturbation loop** (30 tries):
   - Random rotation/translation
   - Add H-bond constraints (ligand → nearest water)
   - Add water network constraints
   - Minimize with constraints
   - Evaluate H-bond geometry (distance, donor/acceptor angles)
   - **Accept if passes H-bond filter**
4. **Pack sidechains** around accepted pose
5. **Post-pack validation** (CRITICAL):
   - Re-evaluate H-bond geometry
   - Reject if geometry degrades during packing
6. **Score** and optionally **cluster** by RMSD
7. **Output**:
   - Accepted PDB files
   - Geometry CSV with detailed diagnostics

### H-Bond Filter Details

**Pre-minimization check**:
- Distance: 2.0-3.6 Å (ideal ± buffer)
- Donor angle: >120° (O-H...Acceptor)
- Quality score: >0.25 (composite)

**Post-pack validation** (stricter):
- Enforces ideal geometry windows
- Rejects poses where packing breaks H-bond
- Ensures ML dataset quality

### Expected Runtime
- ~0.5-2 seconds per conformer (depends on MaxPerturbTries)
- 1000 conformers: 10-30 minutes (single job)
- Array jobs: 2-5 minutes per task (100 conformers/task)

### Output Files

**Docked poses**:
```
scratch/docked_mutant_59K_120A_160G/
├── rep_1.pdb          # First accepted pose
├── rep_2.pdb          # Second accepted pose
├── ...
└── hbond_geometry_summary.csv
```

**Geometry CSV columns**:
- `conf_idx`: Conformer index
- `accepted`: Passed perturbation loop?
- `passed_score`: Score < threshold?
- `saved_cluster`: Unique by RMSD?
- `score`: Rosetta energy (REU)
- `distance`: H-bond distance (Å)
- `donor_angle`: O-H...Acc angle (°)
- `acceptor_angle`: H-bond surface angle (°)
- `quality`: Composite score (0-1)
- `water_residue`: Which water forms H-bond
- `postpack_geometry_passed`: Post-pack validation

---

## Stage 5: Post-Processing and Clustering

### For Array Jobs

After all array tasks complete, aggregate outputs and cluster globally:

```bash
python docking/scripts/cluster_docked_post_array.py \
    ml_modelling/config/mutant_docking_example.conf
```

### What It Does

1. Loads all PDBs from array task outputs
2. Extracts ligand heavy atom coordinates
3. Clusters by RMSD (default 0.75 Å cutoff)
4. Keeps best-scoring pose per cluster
5. Writes clustered poses to final output directory

### Expected Output

```
scratch/docked_mutant_59K_120A_160G/clustered_final/
├── cluster_001_score_-450.2.pdb
├── cluster_002_score_-428.7.pdb
├── ...
└── clustering_summary.csv
```

---

## Stage 6: Feature Extraction for ML

### Input
- Docked poses (from Stage 4/5)
- Mutant PDB
- Variant signature

### Output
- ML feature table (CSV)

### Command

```bash
python ml_modelling/scripts/aggregate_ml_features.py \
    --input-dir scratch/docked_mutant_59K_120A_160G/clustered_final \
    --variant-signature "59K;120A;160G" \
    --ligand-name "WIN-55212" \
    --output ml_features.csv
```

### Features Extracted

**Docking features**:
- `rosetta_score`: Binding energy (REU)
- `hbond_distance`: Water H-bond distance (Å)
- `hbond_donor_angle`: H-bond linearity (°)
- `hbond_quality`: Composite H-bond score

**Geometric features**:
- `ligand_rmsd_to_template`: RMSD vs. reference pose
- `pocket_volume`: Pocket cavity volume (Å³)
- `sasa_ligand`: Ligand solvent-accessible surface area

**Variant features**:
- `num_mutations`: Number of mutations vs. WT
- `mutation_positions`: Comma-separated list
- `mutation_types`: Substitution types (polar→nonpolar, etc.)

---

## Complete Example: Single Ligand + Variant Pair

### Step-by-Step

```bash
# 1. Generate conformers from SMILES
python -m ligand_conformers \
    --input "CCOc1ccc(cc1)C(=O)Nc2ccc(cc2)S(=O)(=O)N" \
    --input-type smiles \
    --output campaign_root/conformers_deprot/WIN_55212 \
    --num-confs 150 \
    --k-final 10

# 2. Create alignment table
python docking/scripts/create_table.py mutant_docking_example.conf

# 3. Thread variant to PDB
python scripts/thread_variant_to_pdb.py \
    --template 3QN1_nolig_H2O.pdb \
    --signature "59K;120A;160G" \
    --output campaign_root/threaded_variants/mutant_59K_120A_160G.pdb

# 4. Dock with water constraints
python docking/scripts/grade_conformers_mutant_docking.py mutant_docking_example.conf

# 5. Extract ML features
python ml_modelling/scripts/aggregate_ml_features.py \
    --input-dir scratch/docked_mutant_59K_120A_160G \
    --variant-signature "59K;120A;160G" \
    --ligand-name "WIN-55212" \
    --output ml_features.csv
```

### Expected Total Runtime
- Conformer generation: 2 min
- Alignment table: 10 sec
- Variant threading: 30 sec
- Docking (10 conformers): 30 sec
- Feature extraction: 5 sec
- **Total: ~4 minutes**

---

## Batch Processing with Orchestrator

For large-scale datasets (many ligand-variant pairs):

```bash
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --ligands ligands.csv \
    --variants variants.csv \
    --output-dir campaign_root \
    --template 3QN1_nolig_H2O.pdb \
    --slurm  # Optional: submit SLURM jobs
```

### What It Does

1. **Enumerate pairs**: Cross-product of ligands × variants
2. **Generate conformers**: For each unique ligand (cached)
3. **Thread variants**: For each unique variant (cached)
4. **Dock**: For each (ligand, variant) pair
5. **Track progress**: Skip completed pairs (metadata.json)
6. **Aggregate features**: Combine all outputs → final ML table

### Expected Runtime (Example)
- 100 ligands × 50 variants = 5,000 pairs
- Conformers: 100 × 2 min = 200 min (3.3 hrs, cached)
- Threading: 50 × 30 sec = 25 min (cached)
- Docking: 5,000 × 30 sec = 42 hrs (parallel: 4 hrs on 10 nodes)
- **Total: ~8 hours wall time** (with caching and parallelization)

---

## Quality Control

### Check H-Bond Geometry

```bash
# View geometry CSV
python -c "
import pandas as pd
df = pd.read_csv('scratch/docked_mutant_59K_120A_160G/hbond_geometry_summary.csv')

# Filter to accepted poses
accepted = df[df['accepted'] == True]
print(f'Acceptance rate: {len(accepted)/len(df):.1%}')

# H-bond statistics
print(f'Distance: {accepted["distance"].mean():.2f} ± {accepted["distance"].std():.2f} Å')
print(f'Donor angle: {accepted["donor_angle"].mean():.1f} ± {accepted["donor_angle"].std():.1f}°')
print(f'Quality: {accepted["quality"].mean():.3f} ± {accepted["quality"].std():.3f}')
"
```

### Visual Inspection

```bash
# Load in PyMOL
pymol campaign_root/threaded_variants/mutant_59K_120A_160G.pdb \
      scratch/docked_mutant_59K_120A_160G/rep_1.pdb

# Check water-ligand H-bond
select water, resn TP3 and (resn TP3 around 3.5 and organic)
distance hbond, water and name O, organic and elem O+N
```

---

## Troubleshooting

### Issue: Low acceptance rate (<10%)

**Causes**:
- Conformers incompatible with pocket shape
- H-bond filter too strict

**Solutions**:
- Increase `MaxPerturbTries` (30 → 50)
- Loosen filter: `HBondQualityMin = 0.1`
- Check water placement in mutant PDB

### Issue: Slow docking (>5 sec/conformer)

**Causes**:
- Large pocket (many waters to check)
- Slow minimization (steric clashes)

**Solutions**:
- Reduce `MaxPerturbTries` (30 → 20)
- Disable water network constraints: `auto_setup_water_constraints()` line 717

### Issue: No poses pass post-pack validation

**Causes**:
- Packing breaks H-bond geometry
- Ideal window too narrow

**Solutions**:
- Disable strict validation: `EnforceFinalIdealGeometry = False`
- Widen buffer: `HBondDistanceIdealBuffer = 1.0`

---

## References

- **PyRosetta Documentation**: https://www.pyrosetta.org/
- **RDKit Conformer Generation**: https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
- **Rosetta H-Bond Geometry**: Gray et al., J Mol Biol (2003)
- **ETKDG Algorithm**: Riniker & Landrum, JCIM (2015)

---

## Contact

For questions or issues:
- Open an issue on GitHub: https://github.com/whitehead-lab/pyr1_pipeline
- Email: <lab-email>

---

## Version History

- **2026-02-16**: Initial version with water-constrained docking
- **Future**: AlphaFold3 integration, LigandMPNN design mode
