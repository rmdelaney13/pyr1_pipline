# Mutant Docking Configuration - CORRECTED Solution

## Problem Re-Analysis ✅

The mutant docking workflow needs **TWO separate PDB files**:

1. **Template PDB (NO ligand)** - For mutation threading
   - Example: `3QN1_nolig_H2O.pdb`
   - Chains: A, B, D (protein + waters)
   - Used to create mutant.pdb via threading

2. **Reference PDB (WITH ligand)** - For SVD alignment
   - Example: `3QN1_H2O.pdb`
   - Chains: A, B, D, **X** (protein + waters + ligand)
   - Used to extract template ligand coordinates (O2, C11, C9)

## Why Both Are Needed

### The Alignment Workflow

```
1. create_table.py generates alignment CSV:
   ↓
   Molecule Atoms: [O1, C5, C7]  (atoms in NEW ligand)
   Target Atoms:   [O2, C11, C9] (atoms in TEMPLATE ligand)

2. Docking script performs SVD alignment:
   ↓
   Get xyz of [O1, C5, C7] from conformer
   Get xyz of [O2, C11, C9] from target_res  ← NEEDS template ligand!
   ↓
   Compute rotation + translation matrix
   ↓
   Transform conformer to align with template

3. Graft aligned conformer into mutant_pose
   ↓
   Perturb, minimize, dock
```

**Critical point**: The `conf.align_to_target(target_res)` call at [line 659](../../docking/scripts/grade_conformers_mutant_docking.py#L659) expects `target_res` to be the **template ligand** with atoms O2, C11, C9 - NOT a protein residue!

## The Corrected Fix

### Changes Made

1. **[orchestrate_ml_dataset_pipeline.py](../scripts/orchestrate_ml_dataset_pipeline.py)**
   - Added `--reference-pdb` CLI argument (points to 3QN1_H2O.pdb WITH ligand)
   - Modified `run_docking()` to include `PrePDBFileName` in config
   - Updated defaults:
     - `reference_chain = 'X'` (ligand chain)
     - `reference_residue = 1` (ligand residue number)

2. **[grade_conformers_mutant_docking.py](../../docking/scripts/grade_conformers_mutant_docking.py)**
   - Modified to load `PrePDBFileName` as a separate reference structure
   - Extracts `target_res` from reference PDB (template ligand)
   - Uses this for SVD alignment
   - Still docks to `mutant_pose` (no ligand)

### New Config Format

```ini
[DEFAULT]
CSVFileName = alignment_table.csv
PathToConformers = conformers_params
ChainLetter = X                    # Template ligand chain
ResidueNumber = 1                  # Template ligand residue
LigandResidueNumber = 1
AutoGenerateAlignment = False
PrePDBFileName = 3QN1_H2O.pdb     # ← NEW: Reference PDB with ligand

[mutant_docking]
MutantPDB = mutant.pdb             # Mutant structure (NO ligand)
LigandSDF = conformers_final.sdf
OutputDir = docking
DockingRepeats = 50
ArrayTaskCount = 10
```

## Updated Usage

### Command Line

```bash
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/data/test_three_examples.csv \
    --cache-dir /scratch/alpine/$USER/test_three_cache \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --reference-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_H2O.pdb \
    --docking-repeats 50 \
    --use-slurm \
    --reference-chain X \
    --reference-residue 1
```

**Key difference**: Now you pass BOTH PDBs:
- `--template-pdb`: NO ligand (for threading)
- `--reference-pdb`: WITH ligand (for alignment)

## Workflow Diagram

```
Input SMILES
   ↓
Conformer Generation (RDKit)
   ↓
Alignment Table (using template ligand atom names O2-C11-C9)
   ↓
┌─────────────────────────────────────────────────────────┐
│ Mutation Threading                                      │
│   template_pdb: 3QN1_nolig_H2O.pdb (NO ligand)         │
│   ↓                                                      │
│   mutant.pdb (chains A, B, D)                           │
└─────────────────────────────────────────────────────────┘
   ↓
┌─────────────────────────────────────────────────────────┐
│ Docking                                                  │
│   reference_pdb: 3QN1_H2O.pdb (WITH ligand chain X)    │
│   ↓                                                      │
│   Extract target_res = chain X, residue 1               │
│   ↓                                                      │
│   For each conformer:                                   │
│     1. Align to target_res (using O2, C11, C9 coords)  │
│     2. Graft into mutant.pdb                            │
│     3. Perturb, minimize, pack                          │
│     4. Validate H-bond geometry                         │
│     5. Score and cluster                                │
└─────────────────────────────────────────────────────────┘
   ↓
Docked poses (rep_*.pdb)
```

## File Comparison

| File | Chains | Purpose | Usage |
|------|--------|---------|-------|
| **3QN1_H2O.pdb** | A, B, D, X | Reference WITH ligand | `--reference-pdb` for alignment |
| **3QN1_nolig_H2O.pdb** | A, B, D | Template NO ligand | `--template-pdb` for threading |
| **mutant.pdb** | A, B, D | Generated mutant | Created by threading |

## Template Ligand Atoms

From 3QN1_H2O.pdb chain X:
```
HETATM 1424  C9  A8T X   1      -0.162  21.525  30.699
HETATM 1426  O2  A8T X   1      -1.932  22.648  29.642
HETATM 1427  C11 A8T X   1      -0.615  23.909  31.306
```

These O2, C11, C9 coordinates are used as the alignment target for all conformers.

## Testing

### Verify Both PDBs Exist
```bash
ls -lh docking/ligand_alignment/files_for_PYR1_docking/3QN1_H2O.pdb
ls -lh docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb
```

### Check Chains
```bash
# Reference PDB (should have chain X)
grep "^HETATM" docking/ligand_alignment/files_for_PYR1_docking/3QN1_H2O.pdb | \
    awk '{print $5}' | sort -u
# Output: A B D X

# Template PDB (should NOT have chain X)
grep "^ATOM\|^HETATM" docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb | \
    awk '{print $5}' | sort -u
# Output: A B D
```

### Run Test
```bash
python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \
    --pairs-csv ml_modelling/data/test_three_examples.csv \
    --cache-dir ml_modelling/cache/test_corrected \
    --template-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_nolig_H2O.pdb \
    --reference-pdb docking/ligand_alignment/files_for_PYR1_docking/3QN1_H2O.pdb \
    --docking-repeats 10 \
    --max-pairs 1
```

### Validate Config
```bash
cat ml_modelling/cache/test_corrected/test_001/docking/docking_config.txt
```

Should show:
```ini
[DEFAULT]
...
ChainLetter = X                                          # ✅ Ligand chain
ResidueNumber = 1                                        # ✅ Ligand residue
PrePDBFileName = .../3QN1_H2O.pdb                       # ✅ Reference with ligand

[mutant_docking]
MutantPDB = .../mutant.pdb                              # ✅ Mutant without ligand
...
```

## Key Differences from Original Fix

| Aspect | Original (WRONG) | Corrected (RIGHT) |
|--------|------------------|-------------------|
| ChainLetter | A (protein) | X (ligand) |
| ResidueNumber | 120 (pocket residue) | 1 (ligand) |
| target_res source | mutant_pose | reference_pose |
| target_res type | Protein residue | Template ligand |
| Alignment atoms | Won't match | O2, C11, C9 ✅ |

## Why the Original Fix Was Wrong

The original fix set `ChainLetter=A, ResidueNumber=120` (a protein residue). This would have failed at alignment because:

1. The alignment table specifies Target Atoms = `[O2, C11, C9]`
2. A protein residue doesn't have atoms named O2, C11, C9
3. `conf.align_to_target(target_res)` would fail to find these atoms
4. SVD alignment would crash

## Summary

✅ **Correct approach**: Use TWO PDBs
- Template PDB (no ligand) for threading → mutant.pdb
- Reference PDB (with ligand) for alignment → target_res

✅ **ChainLetter = X, ResidueNumber = 1** (template ligand)
✅ **PrePDBFileName = 3QN1_H2O.pdb** (with ligand chain X)

This ensures the SVD alignment can find the target atom coordinates (O2, C11, C9) in the template ligand structure.

---

**Status**: ✅ Corrected and ready for testing
**Files Modified**: 2 (orchestrator + mutant docking script)
**Breaking Changes**: Requires new `--reference-pdb` argument
