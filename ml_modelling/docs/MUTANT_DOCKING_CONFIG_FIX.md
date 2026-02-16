# Mutant Docking Configuration Fix

## Problem Summary

The mutant docking workflow is failing because the orchestrator generates a config that inherits `ChainLetter = X` (ligand chain) from the template config, but the mutant PDB has **no ligand chain** - only protein chains (A, B, D).

## Root Cause Analysis

### What's Happening

1. **Orchestrator** ([orchestrate_ml_dataset_pipeline.py:304-316](../scripts/orchestrate_ml_dataset_pipeline.py#L304-L316)) generates this config:
   ```ini
   [DEFAULT]
   CSVFileName = {alignment_csv}
   PathToConformers = {conformers_dir}

   [mutant_docking]
   MutantPDB = {mutant_pdb}
   LigandSDF = {conformers_sdf}
   OutputDir = {output_dir}
   DockingRepeats = {docking_repeats}
   ArrayTaskCount = {array_tasks}
   ```

2. **Mutant docking script** ([grade_conformers_mutant_docking.py:935-938](../../docking/scripts/grade_conformers_mutant_docking.py#L935-L938)) reads:
   ```python
   chain_letter = _cfg_str(def_section, 'ChainLetter', 'A')
   residue_number = _cfg_int(def_section, 'ResidueNumber', 1)
   ```

   Since orchestrator doesn't set these, it inherits from whatever global `[DEFAULT]` section exists.

3. **If the global config has `ChainLetter = X`** (from template PDB with ligand), the script tries to find chain X in the mutant PDB at [line 1020](../../docking/scripts/grade_conformers_mutant_docking.py#L1020):
   ```python
   target_idx = mutant_pose.pdb_info().pdb2pose(chain_letter, residue_number)
   ```

   **This fails** because the mutant PDB has NO chain X!

### Why This Matters

The `target_res` extracted from the mutant PDB is used for:
1. **Alignment atom definitions** ([line 1032](../../docking/scripts/grade_conformers_mutant_docking.py#L1032)) - if `AutoGenerateAlignment=True`
2. **Conformer alignment** ([line 659](../../docking/scripts/grade_conformers_mutant_docking.py#L659)) - each conformer is aligned to target_res

## Solution

The orchestrator must **explicitly set ChainLetter and ResidueNumber** to point to a **protein residue that exists in the mutant PDB**.

### Recommended Fix

Modify `run_docking()` in [orchestrate_ml_dataset_pipeline.py](../scripts/orchestrate_ml_dataset_pipeline.py#L273-L373) to add these parameters:

```python
def run_docking(
    mutant_pdb: Path,
    conformers_sdf: Path,
    alignment_csv: Path,
    conformers_dir: Path,
    output_dir: Path,
    docking_repeats: int = 50,
    use_slurm: bool = False,
    array_tasks: int = 10,
    reference_chain: str = 'A',        # NEW
    reference_residue: int = 120       # NEW (pick a pocket residue)
) -> Optional[str]:
    """
    Run docking to mutant pocket.

    Args:
        ...
        reference_chain: Chain in mutant PDB to use for alignment reference
        reference_residue: Residue number for alignment reference (PDB numbering)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  Docking to mutant pocket ({docking_repeats} repeats)...")

    # Create config file
    config_path = output_dir / 'docking_config.txt'
    with open(config_path, 'w') as f:
        f.write(f"""[DEFAULT]
CSVFileName = {alignment_csv}
PathToConformers = {conformers_dir}
ChainLetter = {reference_chain}          # ADDED
ResidueNumber = {reference_residue}      # ADDED
LigandResidueNumber = 1                  # ADDED (for clarity)
AutoGenerateAlignment = False            # ADDED (table already exists)

[mutant_docking]
MutantPDB = {mutant_pdb}
LigandSDF = {conformers_sdf}
OutputDir = {output_dir}
DockingRepeats = {docking_repeats}
ArrayTaskCount = {array_tasks}
""")

    # ... rest of function unchanged
```

### Alternative: Load Reference Ligand from Template PDB

If you need alignment to a **template ligand structure** (not a protein residue), you'd need to:

1. Add a `PrePDBFileName` parameter pointing to the template PDB with ligand
2. Modify the mutant docking script to:
   - Load the reference PDB
   - Extract target_res from the reference PDB (chain X, residue 1)
   - Use this reference_target_res for alignment
   - But still dock to the mutant_pose

This is more complex and likely not necessary since:
- The alignment table already has alignment atoms defined
- You just need a spatial reference point for SVD alignment

## Which Residue to Use?

Choose a **stable pocket residue** near the binding site. Good candidates:
- A conserved residue in the binding pocket
- A residue involved in ligand recognition (e.g., Tyr 120, Lys 59)
- NOT a residue that's being mutated frequently

Example for PYR1:
```python
reference_chain = 'A'
reference_residue = 120  # Tyr120 gate residue (conserved)
```

## Testing the Fix

After implementing:

1. Check that the config has correct `ChainLetter` and `ResidueNumber`
   ```bash
   cat {cache_dir}/test_001/docking/docking_config.txt
   ```

2. Verify the mutant PDB has that residue:
   ```bash
   grep "^ATOM" mutant.pdb | grep " A " | awk '{print $5}' | sort -u | grep " 120"
   ```

3. Run a single test case locally (not SLURM) to verify alignment works

## Impact Assessment

This is a **critical bug** that blocks the entire ML dataset generation pipeline. The fix is:
- Low risk (just adding explicit config parameters)
- High impact (unblocks docking → enables relax → enables AF3 → enables ML dataset)
- Backward compatible (doesn't change the docking algorithm)

## Files to Modify

1. **[ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py](../scripts/orchestrate_ml_dataset_pipeline.py)**
   - Function: `run_docking()` (lines 273-373)
   - Function: `process_single_pair()` (line 589 - add reference residue params to call)

2. **Optional**: Add to orchestrator argument parser:
   ```python
   parser.add_argument('--reference-chain', default='A', help='Reference chain for alignment')
   parser.add_argument('--reference-residue', type=int, default=120, help='Reference residue for alignment (PDB numbering)')
   ```

---

**Status**: Documentation complete, ready for implementation
**Priority**: P0 (blocking)
**Estimated effort**: 15 minutes
