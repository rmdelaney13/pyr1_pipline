# Sequence-CSV Docking Workflow Guide

This guide explains how to dock ligands to protein sequences specified in a CSV file, with each sequence threaded onto specific pocket positions.

## Overview

The sequence-CSV docking pipeline allows you to:
- **Thread protein sequences** from a CSV onto specific positions (e.g., the 16 positions defining your binding pocket)
- **Run one array job per sequence** for efficient parallelization
- **Use the same rigorous H-bond geometry constraints** as the glycine-shaved pipeline
- **Organize outputs by sequence name** for easy comparison

## Quick Start

### 1. Prepare Your Sequence CSV

Create a CSV file with two columns: `name` and `signature`

```csv
name,signature
binder_001,QILMVAGDHVASILYQ
binder_002,SILQMAGDHVYTNPWL
wild_type,QILMVAGDHVASILYM
```

**Requirements:**
- **16 amino acids** (matching the 16 pocket positions)
- **Valid amino acid codes**: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
- **Unique names** for each sequence (used for output directories)

### 2. Configure Your Docking

Copy the template config and edit it:

```bash
cp docking/templates/config_sequence_docking.txt my_sequence_docking.txt
```

Edit the `[SEQUENCE_DOCKING]` section:

```ini
[SEQUENCE_DOCKING]
# Path to your CSV
SequenceCSV = /path/to/your/sequences.csv

# Threading positions (16 positions for PYR1 pocket)
ThreadingPositions = 59 79 81 90 92 106 108 115 118 120 139 157 158 161 162 165
```

### 3. Test Locally (Recommended)

Before submitting to SLURM, test with the first sequence:

```bash
# Test first sequence only
python docking/scripts/submit_sequence_csv_docking.py my_sequence_docking.txt --local --max-local 1

# Preview all sequences without running
python docking/scripts/submit_sequence_csv_docking.py my_sequence_docking.txt --preview
```

### 4. Submit to SLURM

Once testing looks good, submit the full array:

```bash
# Generate submission script and review
python docking/scripts/submit_sequence_csv_docking.py my_sequence_docking.txt

# Auto-submit
python docking/scripts/submit_sequence_csv_docking.py my_sequence_docking.txt --submit
```

This will create one SLURM array job with one task per sequence.

## Output Organization

Outputs are organized by sequence name (Option C):

```
{OutputDirBase}/
└── sequences/
    ├── binder_001/
    │   ├── post_mutate_binder_001.pdb      # Threaded structure
    │   ├── rotated_pass1/                   # Initial docked poses
    │   ├── pass_score_repacked_pass1/       # Passing poses (pre-clustering)
    │   ├── clustered_final/                 # ✓ FINAL CLUSTERED RESULTS
    │   │   ├── pass1_repacked_cluster0001.pdb
    │   │   ├── pass1_repacked_cluster0002.pdb
    │   │   └── ...
    │   ├── hbond_geometry_pass1.csv         # H-bond geometry metrics
    │   └── summary_binder_001.json          # Docking statistics
    ├── binder_002/
    │   └── ...
    └── wild_type/
        └── ...
```

### Key Output Files

1. **`clustered_final/`**: Best unique docked poses (use these for modeling)
2. **`rotated_pass1/`**: All initial docked poses (for inspection)
3. **`pass_score_repacked_pass1/`**: All poses passing score cutoff (before clustering)
4. **`hbond_geometry_pass*.csv`**: Detailed H-bond geometry for each conformer
5. **`summary_{seq_name}.json`**: Statistics summary

## Understanding the Pipeline

### Sequence Threading

For each sequence, the pipeline:
1. Loads the 16-letter signature from CSV
2. Threads each amino acid onto the corresponding position:
   ```
   Position:  59 79 81 90 92 106 108 115 118 120 139 157 158 161 162 165
   Signature: Q  I  L  M  V  A   G   D   H   V   A   S   I   L   Y   Q
   ```
3. Saves the threaded structure as `post_mutate_{seq_name}.pdb`

### Docking with H-bond Constraints

The docking uses the same rigorous water-mediated H-bond filtering as the glycine-shaved pipeline:

1. **Initial sampling**: Rigid-body perturbation with collision checking
2. **Water constraints**: Automatically detects nearby waters and adds H-bond constraints
3. **Minimization**: Pulls ligand toward water with strong constraint weights
4. **Geometry filtering**: Accepts only poses with ideal H-bond geometry
5. **Score filtering**: Keeps poses below energy threshold
6. **Clustering**: Removes redundant poses by RMSD

### Configuration Parameters

Key parameters in `[grade_conformers]`:

```ini
# Sampling
MaxPerturbTries = 30              # Attempts per conformer
Rotation = 25                      # Degrees rotation
Translation = 0.5                  # Angstroms translation

# H-bond filtering (CRITICAL for water-mediated binding)
EnableHBondGeometryFilter = True
HBondConstraintWeight = 4.0        # Strong pull during minimization
HBondDistanceIdeal = 2.8           # Ideal H-bond distance (Å)
HBondDonorAngleMin = 120           # Minimum donor angle
EnforceFinalIdealGeometry = True   # Strict final check

# Scoring and clustering
MaxScore = -300                    # Energy cutoff
ClusterRMSDCutoff = 0.75          # RMSD threshold (Å)
```

## Analyzing Results

### 1. Check Summary Statistics

```bash
# View summary for a specific sequence
cat {OutputDirBase}/sequences/binder_001/summary_binder_001.json
```

Key metrics:
- `total_conformers`: Total conformers attempted
- `accepted_conformers`: Passed H-bond geometry
- `pass_score_conformers`: Passed energy cutoff
- `final_unique_clustered_docks`: Final unique poses in `clustered_final/`

### 2. Inspect H-bond Geometry

```bash
# Load in pandas or Excel
import pandas as pd
df = pd.read_csv('hbond_geometry_pass1.csv')

# Filter for accepted poses
accepted = df[df['saved_cluster'] == True]
print(accepted[['conf_num', 'score', 'distance', 'donor_angle', 'quality']])
```

### 3. Visualize Docked Poses

Load the final clustered poses in PyMOL or ChimeraX:

```python
# PyMOL
load {OutputDirBase}/sequences/binder_001/clustered_final/*.pdb

# Color by sequence
color marine, binder_001*
color orange, binder_002*
```

### 4. Compare Across Sequences

To find the best overall docks across all sequences:

```bash
# Collect all summaries
python << EOF
import json
import glob

summaries = []
for f in glob.glob('{OutputDirBase}/sequences/*/summary_*.json'):
    with open(f) as h:
        data = json.load(h)
        summaries.append({
            'name': data['sequence_name'],
            'total_docks': data['totals']['final_unique_clustered_docks'],
            'pass_score': data['totals']['pass_score_conformers']
        })

for s in sorted(summaries, key=lambda x: -x['total_docks']):
    print(f"{s['name']:20s} {s['total_docks']:4d} docks  ({s['pass_score']:4d} passed score)")
EOF
```

## Advanced Usage

### Running Multiple Passes

Increase sampling by running multiple independent passes:

```ini
[MULTIPLE_PASSES]
NumPasses = 3
```

Each pass uses different random perturbations. Final outputs aggregate all passes.

### Custom SLURM Parameters

```bash
python submit_sequence_csv_docking.py config.txt \
    --job-name my_docking \
    --time 48:00:00 \
    --partition amilan \
    --mem 32G \
    --submit
```

### Running Specific Sequences Manually

```bash
# Run sequence at index 0 (first row in CSV)
python grade_conformers_sequence_csv_docking_multiple_slurm.py config.txt 0

# Run sequence at index 5
python grade_conformers_sequence_csv_docking_multiple_slurm.py config.txt 5
```

## Troubleshooting

### "Signature length does not match positions count"

Make sure your signature has exactly 16 amino acids (or however many positions you specified).

### "Array index out of range"

Check that your SLURM array size matches the CSV:
```bash
# If you have 8 sequences, array should be 0-7
#SBATCH --array=0-7
```

### No docks in clustered_final/

Common causes:
1. **Too strict geometry filter**: Lower `HBondDonorAngleMin` or increase `HBondDistanceIdealBuffer`
2. **Score cutoff too low**: Increase `MaxScore` (less negative)
3. **Not enough conformers**: Check your input SDF has diverse conformers

Check `hbond_geometry_pass1.csv` to see why conformers failed.

### Out of memory

Increase SLURM memory:
```bash
python submit_sequence_csv_docking.py config.txt --mem 32G --submit
```

## Tips for Success

1. **Start small**: Test with 1-2 sequences locally before submitting full array
2. **Check intermediate outputs**: The `rotated_pass*/` dirs show all attempts, good for debugging
3. **Adjust H-bond parameters**: If you're not getting docks, the geometry filter might be too strict
4. **Use multiple passes**: More sampling = better coverage of conformational space
5. **Compare to AF3/Boltz**: You mentioned you already have AF3 models - compare Rosetta docks to those for validation

## Integration with Existing Workflow

This sequence-CSV pipeline complements your existing tools:

| Tool | Use Case |
|------|----------|
| **Glycine-shaved docking** | Explore pocket without sequence bias |
| **Sequence-CSV docking** | Dock to specific sequences (binders, variants) |
| **AF3/Boltz** | Full complex prediction |
| **Rosetta (this pipeline)** | Fine-grained water-mediated interactions |

You can use Rosetta docks as starting structures for AF3, or vice versa!

## Next Steps

After docking:
1. **Relax structures**: Use Rosetta relax with constraints
2. **Score with AF3**: Use docked poses as templates for AF3 modeling
3. **Analyze contacts**: Identify key ligand-water-protein interactions
4. **Design improvements**: Use insights to design better binders

## Questions?

Check the main docking guides:
- `docking/QUICKSTART.md` - General pipeline overview
- `docking/WORKFLOW_README.md` - Detailed workflow steps
- `docking/SAMPLING_GUIDE.md` - Sampling and clustering strategies
