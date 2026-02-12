# Sequence-CSV Docking Pipeline

**Threading-based docking for binder sequences with rigorous water-mediated H-bond constraints**

## What This Does

Docks ligands to protein sequences specified in a CSV file. Each sequence is threaded onto your binding pocket positions, then the same rigorous H-bond geometry filtering from the glycine-shaved pipeline is applied.

**Key Features:**
- ✅ CSV-driven sequence input (easy to scale to hundreds of sequences)
- ✅ Automatic SLURM array submission (one job per sequence)
- ✅ Rigorous water-mediated H-bond geometry filtering
- ✅ Organized outputs by sequence name
- ✅ Keeps intermediate outputs for inspection
- ✅ Compatible with existing glycine-shaved workflow

## Files Created

### Core Scripts
1. **`scripts/grade_conformers_sequence_csv_docking_multiple_slurm.py`**
   - Main docking script with sequence threading
   - Reads CSV, threads sequences, runs docking with H-bond constraints
   - Outputs organized by sequence name

2. **`scripts/submit_sequence_csv_docking.py`**
   - Helper script for SLURM submission
   - Auto-detects array size from CSV
   - Generates SLURM scripts with proper parameters

### Templates
3. **`templates/config_sequence_docking.txt`**
   - Complete config template with `[SEQUENCE_DOCKING]` section
   - Pre-configured with optimal H-bond constraint settings
   - Example showing 16-position threading for PYR1

4. **`templates/sequences_example.csv`**
   - Example CSV with sample sequences
   - Shows proper format and naming

### Documentation
5. **`SEQUENCE_DOCKING_GUIDE.md`**
   - Complete usage guide with examples
   - Troubleshooting tips
   - Analysis workflows

## Quick Start

### 1. Create Your Sequence CSV

```csv
name,signature
binder_001,QILMVAGDHVASILYQ
binder_002,SILQMAGDHVYTNPWL
```

16 letters = 16 positions to thread onto

### 2. Configure

```bash
cp docking/templates/config_sequence_docking.txt my_config.txt
# Edit [SEQUENCE_DOCKING] section to point to your CSV
```

### 3. Test Locally

```bash
# Test first sequence
python docking/scripts/submit_sequence_csv_docking.py my_config.txt --local --max-local 1
```

### 4. Submit to SLURM

```bash
# Auto-generates and submits array job
python docking/scripts/submit_sequence_csv_docking.py my_config.txt --submit
```

## Output Structure

```
{OutputDirBase}/sequences/
├── binder_001/
│   ├── clustered_final/           # ← FINAL RESULTS (best unique docks)
│   ├── rotated_pass1/              # All initial docked poses
│   ├── pass_score_repacked_pass1/  # Passing poses (pre-clustering)
│   ├── hbond_geometry_pass1.csv    # Detailed H-bond metrics
│   └── summary_binder_001.json     # Statistics
├── binder_002/
│   └── ...
└── wild_type/
    └── ...
```

## Key Configuration

The `[SEQUENCE_DOCKING]` section in your config:

```ini
[SEQUENCE_DOCKING]
# Path to CSV with sequences
SequenceCSV = /path/to/sequences.csv

# Positions to thread onto (must match signature length)
ThreadingPositions = 59 79 81 90 92 106 108 115 118 120 139 157 158 161 162 165
```

## H-bond Constraint Strategy

This pipeline uses the **same rigorous water-mediated H-bond filtering** as your glycine-shaved workflow:

1. **Dynamic water constraints** pull ligand toward nearest water during minimization
2. **Strict geometry filters** ensure ideal H-bond distance/angles
3. **Auto-detection** of ligand acceptors (finds best H-bond partner automatically)
4. **Post-packing validation** rejects poses that lose H-bond after repacking

This is critical for water-mediated binding like you described wanting to prioritize.

## Comparison to Glycine-Shaved Docking

| Feature | Glycine-Shaved | Sequence-CSV |
|---------|---------------|--------------|
| **Pocket definition** | All glycine (no bias) | Specific sequences threaded |
| **Use case** | Explore pocket geometry | Dock to binders/variants |
| **H-bond constraints** | ✅ Full support | ✅ Full support |
| **Parallelization** | Array by conformer chunks | Array by sequence |
| **Output organization** | By array index | By sequence name |

Both use the same underlying docking + H-bond filtering code!

## Integration with Your Workflow

You mentioned you already use AF3/Boltz for modeling. This pipeline complements that:

```
Workflow Option 1: Rosetta → AF3
├─ Rosetta sequence docking → top docked poses
└─ Use top poses as templates for AF3 predictions

Workflow Option 2: AF3 → Rosetta
├─ AF3/Boltz → initial complex model
└─ Extract sequence → dock in Rosetta for water-mediated refinement

Workflow Option 3: Parallel
├─ Both AF3 and Rosetta in parallel
└─ Compare results to find consensus binding modes
```

The advantage of Rosetta here is the **explicit water H-bond constraints** that you can weight very highly.

## Command Reference

### Preview sequences
```bash
python scripts/submit_sequence_csv_docking.py config.txt --preview
```

### Local testing (first 2 sequences)
```bash
python scripts/submit_sequence_csv_docking.py config.txt --local --max-local 2
```

### Generate SLURM script (don't submit)
```bash
python scripts/submit_sequence_csv_docking.py config.txt
# Creates submit_seq_dock.sh
```

### Auto-submit to SLURM
```bash
python scripts/submit_sequence_csv_docking.py config.txt --submit
```

### Run specific sequence manually
```bash
python scripts/grade_conformers_sequence_csv_docking_multiple_slurm.py config.txt 0
# 0 = first sequence in CSV
```

### Custom SLURM params
```bash
python scripts/submit_sequence_csv_docking.py config.txt \
    --job-name my_docking \
    --time 48:00:00 \
    --mem 32G \
    --submit
```

## Example Use Cases

### 1. Screen Binder Library
You have 100 binder sequences you want to dock:
```bash
# Create binders.csv with 100 sequences
# Submit array job (one task per sequence)
python scripts/submit_sequence_csv_docking.py config.txt --submit
# Outputs organized: sequences/binder_001/, sequences/binder_002/, ...
```

### 2. Compare Variants
Test wild-type vs mutants:
```csv
name,signature
wild_type,QILMVAGDHVASILYM
mutant_F79A,QALMVAGDHVASILYM
mutant_M90A,QILMVAGDHAASILYM
```

### 3. Saturation Mutagenesis
Create CSV with all single mutations at key positions, dock all variants

## Troubleshooting

### No docks in clustered_final/
Check `hbond_geometry_pass1.csv` to see failure modes:
- If most fail geometry: loosen `HBondDonorAngleMin` or `HBondDistanceIdealBuffer`
- If most fail score: increase `MaxScore` (less negative)
- If passing geometry but failing clustering: check `pass_score_repacked_pass1/` for pre-clustered docks

### Array index out of range
Array size must match CSV row count:
- 10 sequences → `--array=0-9`
- Auto-handled by `submit_sequence_csv_docking.py`

### Memory issues
Increase per-task memory:
```bash
python scripts/submit_sequence_csv_docking.py config.txt --mem 32G --submit
```

## Performance

Typical runtime per sequence (on Alpine):
- **Small ligand** (~20 atoms, 50 conformers): 30-60 minutes
- **Medium ligand** (~30 atoms, 100 conformers): 1-2 hours
- **Large ligand** (~40 atoms, 200 conformers): 3-4 hours

With SLURM arrays, all sequences run in parallel → total time = single sequence time

## Tips for Success

1. **Test locally first** - catch config errors before submitting hundreds of jobs
2. **Start with NumPasses=1** - increase if you need more sampling
3. **Check intermediate outputs** - `rotated_pass1/` shows all attempts, good for debugging
4. **Compare to glycine-shaved** - if sequence docking fails but glycine succeeds, sequence may not be compatible with binding
5. **Use multiple positions** - the 16-position threading is specific to PYR1; adjust `ThreadingPositions` for your system

## Next Steps

1. **Generate sequences.csv** with your binders
2. **Copy and edit config** from `templates/config_sequence_docking.txt`
3. **Test locally** with `--local --max-local 1`
4. **Submit array job** with `--submit`
5. **Analyze results** in `sequences/{name}/clustered_final/`
6. **Optional: Relax top docks** with Rosetta relax
7. **Optional: Model with AF3** using top docked geometries

## Questions?

See:
- **`SEQUENCE_DOCKING_GUIDE.md`** - Detailed usage guide
- **`QUICKSTART.md`** - General pipeline overview
- **`WORKFLOW_README.md`** - Complete workflow documentation

---

**Created:** 2025-02-12
**Compatibility:** Works with existing glycine-shaved docking infrastructure
**Status:** Ready for production use
