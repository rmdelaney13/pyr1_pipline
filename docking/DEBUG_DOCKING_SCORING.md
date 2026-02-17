# Debugging Docking & Scoring — Unrealistic Scores

## Quick Reference: What the Scores Mean

```
score           = sf_all(protein + ligand + waters)   # absolute Rosetta energy
baseline_score  = sf_all(protein + waters, NO ligand) # computed once at start
relative_score  = score - baseline_score              # energy CHANGE from adding ligand
```

**Expected ranges:**
- `baseline_score`: -200 to -1500 REU (depends on protein size and relaxation)
- `relative_score`: -100 to +50 REU for good binders; +50 to +100 for marginal
- `relative_score` > 100: rejected by `MaxScore = 100` threshold

**If you see `relative_score` in the thousands:** the ligand is in catastrophic clashes.

---

## Symptom → Cause Lookup

| Symptom | Likely Cause | Debug Step |
|---------|-------------|------------|
| `relative_score` in thousands (2000+) | Severe VdW clashes between ligand and protein | Steps 1, 2, 3 |
| `quality = 0.0` for most rows | No water within H-bond distance of ligand acceptor | Step 5 |
| `distance/donor_angle/acceptor_angle` all empty | `evaluate_hbond_geometry` found no nearby water oxygen | Step 5 |
| `postpack_geometry_passed=True, postpack_ideal_passed=False` | Water exists but H-bond geometry outside ideal window | Step 6 |
| `accepted=True` but `passed_score=False` for ALL rows | Every pose clashes — alignment or params problem | Steps 1, 2, 4 |
| `baseline_score` very positive (> 0) | Mutant PDB has internal clashes (needs relaxation) | Step 7 |

---

## Step-by-Step Debug Protocol

### Step 1: Visualize a Docked Pose

Even rejected poses are saved to `rotated/` directories. Open one in PyMOL:

```bash
# Find rotated (pre-pack) PDBs
ls <pair_cache>/docking/rotated_*/rotated_*.pdb

# Or find passed PDBs (if any exist)
ls <pair_cache>/docking/*rep_*.pdb
```

In PyMOL:
```
load rotated_0.pdb
show sticks, organic         # show ligand
show surface, polymer        # show protein surface
# Look for: ligand buried inside protein, ligand clashing with sidechains
```

**What to look for:**
- Ligand atoms overlapping protein atoms → VdW clash (fa_rep issue)
- Ligand completely outside the binding pocket → alignment problem
- Ligand in pocket but buried too deep → perturbation too aggressive
- Waters displaced or overlapping ligand → water scoring issue

---

### Step 2: Score Breakdown by Energy Term

Add this diagnostic snippet to `grade_conformers_mutant_docking.py` AFTER line 862
(`score = sf_all(copy_pose)`):

```python
# === DIAGNOSTIC: per-term energy breakdown ===
if global_conf_idx <= 3:  # Only for first 3 conformers
    energies = copy_pose.energies()
    sf_all(copy_pose)  # ensure energies are populated
    from pyrosetta.rosetta.core.scoring import ScoreType
    terms_to_check = [
        ScoreType.fa_rep,    # VdW repulsion (clashes)
        ScoreType.fa_atr,    # VdW attraction
        ScoreType.fa_elec,   # Electrostatics
        ScoreType.fa_sol,    # Solvation
        ScoreType.hbond_sc,  # Sidechain H-bonds
        ScoreType.ref,       # Reference energies
    ]
    logger.info("=== SCORE BREAKDOWN for conformer %d ===", conf_idx)
    for term in terms_to_check:
        logger.info("  %-15s = %.2f", term.name, sf_all.score_by_scoretype(copy_pose, term))
    logger.info("  TOTAL           = %.2f", score)
    logger.info("  BASELINE        = %.2f", baseline_score)
    logger.info("  RELATIVE        = %.2f", relative_score)
    logger.info("=== END BREAKDOWN ===")
# === END DIAGNOSTIC ===
```

**What to look for:**
- `fa_rep` >> 1000: severe atom overlaps (clashes)
- `fa_elec` >> 500: charge-charge repulsion (bad atom types in params)
- `fa_sol` >> 500: solvation penalty (ligand buried improperly)

---

### Step 3: Score With vs Without Waters

Add this after the score breakdown (also diagnostic):

```python
# === DIAGNOSTIC: score without waters ===
if global_conf_idx <= 3:
    no_water_pose = copy_pose.clone()
    # Delete waters from highest index to lowest
    water_res = [i for i in range(1, no_water_pose.total_residue()+1)
                 if no_water_pose.residue(i).is_water()]
    for wr in reversed(water_res):
        no_water_pose.delete_residue_slow(wr)
    score_no_water = sf_all(no_water_pose)
    logger.info("Score without waters: %.2f (vs with waters: %.2f, diff: %.2f)",
                score_no_water, score, score - score_no_water)
# === END DIAGNOSTIC ===
```

If removing waters dramatically reduces the score → water placement/params are the issue.

---

### Step 4: Check Ligand Params File

The conformer params files are generated from SDF by `molfile_to_params.py`.
Verify the params look reasonable:

```bash
# Check the params file for the ligand being docked
cat <conformers_dir>/<mol_id>/<mol_id>.params

# Look for:
# - ATOM lines with reasonable types (aromatic C should be "aroC", not "CH1")
# - BOND lines that match the chemistry
# - CHI lines for rotatable bonds
# - NBR_ATOM and NBR_RADIUS
```

**Common params issues:**
- Wrong atom types → extreme `fa_rep` or `fa_elec`
- Missing or wrong `NBR_ATOM` → scoring doesn't apply neighbor corrections
- Hydrogen count mismatch → extra repulsion

Quick validation:
```python
# In a PyRosetta session:
pose = Pose()
res_set = pyrosetta.generate_nonstandard_residue_set(pose, ['path/to/ligand.params'])
pyrosetta.pose_from_file(pose, res_set, 'path/to/ligand_0001.pdb')
print("Ligand-only score:", pyrosetta.get_fa_scorefxn()(pose))
# Should be roughly -10 to +50 for a drug-like molecule
# If it's > 500 → params file is broken
```

---

### Step 5: Check Water Presence and Proximity

The H-bond geometry filter requires a water molecule near the ligand acceptor atom.
If waters are missing or displaced, `quality = 0.0`.

```python
# === DIAGNOSTIC: water proximity check ===
# Add after ligand grafting (after line 694)
if global_conf_idx <= 3:
    lig = grafted_pose.residue(lig_idx)
    logger.info("Ligand residue %d: %s (%d atoms)", lig_idx, lig.name3(), lig.natoms())

    # Count waters in the pose
    water_count = sum(1 for i in range(1, grafted_pose.total_residue()+1)
                      if grafted_pose.residue(i).is_water())
    logger.info("Waters in pose: %d", water_count)

    # Find closest water to each acceptor
    for ai in range(1, lig.natoms()+1):
        try:
            if not lig.heavyatom_is_an_acceptor(ai):
                continue
        except:
            continue
        acc_xyz = lig.xyz(ai)
        best_dist = None
        for ri in range(1, grafted_pose.total_residue()+1):
            w = grafted_pose.residue(ri)
            if not w.is_water():
                continue
            for oi in range(1, w.natoms()+1):
                if w.atom_type(oi).element() == "O":
                    d = acc_xyz.distance(w.xyz(oi))
                    if best_dist is None or d < best_dist:
                        best_dist = d
        logger.info("  Acceptor %s: closest water O = %.2f Å",
                    lig.atom_name(ai).strip(), best_dist if best_dist else -1)
# === END DIAGNOSTIC ===
```

**What to look for:**
- `Waters in pose: 0` → waters weren't loaded (check PDB format / TP3 naming)
- All water distances > 5.0 Å → ligand placed too far from water network
- Water distance 1.0-2.0 Å → ligand overlapping with water → severe clash

---

### Step 6: Check Collision Grid Coverage

The backbone-only collision grid (`include_sc=False`) misses sidechain clashes.
This is intentional (ligand should fit between sidechains after repacking), but
if the pocket has large sidechains (W, F, Y, R, K), clashes may persist after repacking.

To test if sidechain clashes are the problem, temporarily switch to full-atom collision:

```python
# In grade_conformers_mutant_docking.py, change line 622:
#   include_sc=False,     # Original
    include_sc=True,      # Test: include sidechains
```

**If this dramatically reduces accepted conformers:** sidechains are blocking the pocket.
Consider:
- Pre-relaxing the mutant structure before docking
- Using a smaller `vdw_modifier` (e.g., 0.5 instead of 0.7) for more tolerance
- Checking if mutations introduce bulky residues into the binding pocket

---

### Step 7: Check Baseline Score of Mutant PDB

If the mutant PDB itself has bad geometry, baseline will be extreme:

```python
# Quick check script:
import pyrosetta
pyrosetta.init('-extra_res_fa A8T.params -mute all')
from pyrosetta import Pose
pose = Pose()
pyrosetta.pose_from_file(pose, 'path/to/mutant.pdb')
sf = pyrosetta.get_fa_scorefxn()
print(f"Baseline score: {sf(pose):.2f} REU")
print(f"Residues: {pose.total_residue()}")
print(f"Score per residue: {sf(pose)/pose.total_residue():.2f}")
# Reasonable: -2 to -5 REU/residue
# Bad: > 0 REU/residue (clashes in structure)
```

If baseline is positive → the mutant PDB needs relaxation before docking.

---

### Step 8: Verify Ligand Alignment

The ligand is aligned to a template ligand (A8T in 3QN1_H2O.pdb) via SVD.
If the alignment atoms are wrong, the ligand could end up in the wrong place entirely.

```bash
# Check the alignment table CSV
cat <pair_cache>/docking/alignment_table.csv
# Look for Molecule Atoms and Target Atoms columns
# These should be atom names that exist in both the conformer and template
```

```python
# Verify alignment visually:
# Load both reference (with A8T) and a docked pose
# In PyMOL:
load 3QN1_H2O.pdb, reference
load rotated_0.pdb, docked
align reference and polymer, docked and polymer
# The docked ligand should overlap the template ligand position
```

---

## Config Parameters That Affect Scoring

From `docking_config.txt` (generated by orchestrator):

```ini
[mutant_docking]
MaxScore = 100              # relative_score threshold (REU)
EnablePocketProximityFilter = True
PocketMaxDistance = 8.0      # max COM distance from pocket center (Å)

# Not in orchestrator config (using defaults):
# Rotation = 25.0           # degrees of random rotation per try
# Translation = 0.5         # Å of random translation per try
# VDW_Modifier = 0.7        # collision grid leniency (lower = more tolerant)
# MaxPerturbTries = 30       # attempts before giving up on a conformer
# HBondConstraintWeight = 4.0  # constraint strength during minimization
```

**Tuning tips for hard-to-dock ligands:**
- Reduce `Rotation` to 10-15 and `Translation` to 0.25 for finer sampling near the template position
- Reduce `VDW_Modifier` to 0.5 to allow tighter approach before rejection
- Increase `MaxPerturbTries` to 50-100 for more sampling
- Set `HBondQualityMin = 0.1` to reject poses with no real H-bond interaction
- Set `MaxScore = 200` temporarily to see if any poses are close to acceptable

---

## Common Root Causes Summary

### 1. Ligand-Sidechain VdW Clashes (Most Likely)
**Evidence:** `fa_rep` term >> 1000 in score breakdown
**Fix:** Either include sidechains in collision grid, or pre-relax the mutant structure

### 2. Water Scoring Artifacts
**Evidence:** Scores drop dramatically when waters are removed (Step 3)
**Fix:** Check water residue naming (must be TP3/HOH), ensure waters aren't overlapping ligand

### 3. Bad Ligand Params
**Evidence:** Ligand-only score >> 100 (Step 4)
**Fix:** Regenerate params with `molfile_to_params.py -n <name> -p <name> --conformers-in-one-file`

### 4. Ligand Misplacement
**Evidence:** Ligand in wrong position when visualized (Step 1)
**Fix:** Check alignment atoms in CSV, verify SVD alignment by overlaying with template

### 5. Unrelaxed Mutant Structure
**Evidence:** Baseline score positive or very high for the protein size
**Fix:** Run FastRelax on the mutant PDB before docking

---

## CSV Column Reference

| Column | Source | Meaning |
|--------|--------|---------|
| `conf_idx` | Loop counter | Conformer index (1-based) |
| `conf_num` | Conformer file | Original conformer number from SDF |
| `accepted` | Pre-pack geometry | Passed initial H-bond + collision check |
| `accepted_try` | Perturbation loop | Which try (1-30) succeeded |
| `passed_score` | Score filter | `relative_score < MaxScore` AND post-pack geometry OK |
| `saved_cluster` | RMSD clustering | New cluster representative (not a duplicate) |
| `score` | `sf_all(pose)` | Absolute Rosetta energy (protein + ligand + water) |
| `relative_score` | `score - baseline` | Energy change from adding ligand |
| `distance` | Post-pack H-bond | Acceptor-water O distance (Å) |
| `donor_angle` | Post-pack H-bond | D-H...O angle (degrees, ideal ≈ 180) |
| `acceptor_angle` | Post-pack H-bond | H...O-acceptor angle (degrees) |
| `quality` | Post-pack H-bond | Composite quality (0.0 = no interaction, 1.0 = perfect) |
| `water_residue` | Post-pack H-bond | Rosetta residue index of bridging water |
| `strict_window_passed` | Ideal window | Distance within ±0.8 Å of 2.8 Å AND donor angle > 140° |
| `postpack_geometry_passed` | Post-pack check | Basic H-bond geometry criteria met |
| `postpack_ideal_passed` | Post-pack ideal | Strict ideal geometry window met |
| `pocket_distance` | Pocket filter | Ligand COM distance from pocket center |
| `pocket_passed` | Pocket filter | `pocket_distance ≤ PocketMaxDistance` |
| `output_pdb` | File output | Path to saved PDB (only for cluster reps) |

---

## Pipeline Flow Diagram

```
mutant.pdb (no ligand) + conformers SDF
    │
    ▼
[1] Load mutant, compute baseline_score = sf_all(mutant)
    │
    ▼
[2] Build backbone collision grid (exclude waters, backbone only)
    │
    ▼
[3] For each conformer × docking_repeats:
    │
    ├─ [3a] Align ligand to template via SVD
    │
    ├─ [3b] Graft ligand onto mutant via append_residue_by_jump
    │
    ├─ [3c] Perturbation loop (up to MaxPerturbTries):
    │   ├─ Random rigid-body perturbation (Rotation°, Translation Å)
    │   ├─ Backbone collision check → skip if clashing
    │   ├─ Add H-bond constraints (to water or target residue)
    │   ├─ Minimize with constraints (ligand jump only)
    │   └─ Evaluate H-bond geometry → accept if quality ≥ threshold
    │
    ├─ [3d] Pack sidechains (RestrictToRepacking)
    │
    ├─ [3e] score = sf_all(pose), relative_score = score - baseline
    │
    ├─ [3f] Post-pack H-bond geometry check
    │   └─ Must pass both basic AND ideal window criteria
    │
    ├─ [3g] Score threshold: relative_score < MaxScore (100)
    │
    ├─ [3h] Pocket proximity check: COM within PocketMaxDistance (8.0 Å)
    │
    └─ [3i] RMSD clustering (keep unique representatives)
            └─ Save PDB for new cluster reps
```
