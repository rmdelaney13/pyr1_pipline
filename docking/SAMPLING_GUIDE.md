# Docking Sampling and Clustering Guide

## How to Know if You Have Enough Docks

This is indeed a judgment call, but here are practical guidelines to assess sampling adequacy.

---

## Quick Assessment Tool

Run this after your docking completes:

```bash
python scripts/analyze_sampling_adequacy.py config.txt --ligand-heavy-atoms 15
```

This will analyze your results and give recommendations for:
- Whether you have sufficient sampling
- If your RMSD cutoff is appropriate
- Score distribution and convergence

---

## 1. RMSD Cutoff Selection

### Rule of Thumb by Ligand Size

| Ligand Heavy Atoms | Recommended RMSD | Rationale |
|-------------------|------------------|-----------|
| **< 10 atoms** | 0.5 Å | Small, rigid ligands - tight clustering needed |
| **10-20 atoms** | **0.75 Å** | Standard medium-sized ligands (your default) |
| **20-30 atoms** | 1.0 Å | Larger ligands with more conformational freedom |
| **> 30 atoms** | 1.5 Å | Very large/flexible ligands |

### Adjust Based on Results

**INCREASE RMSD cutoff if:**
- You get >200 clusters (over-fragmented)
- Visual inspection shows clusters contain near-identical poses
- Many clusters have only 1-2 members
- **Example**: Change from 0.75 → 1.0 Å

**DECREASE RMSD cutoff if:**
- You get <10 clusters (under-fragmented)
- Visual inspection shows clearly different poses in the same cluster
- You want more granular discrimination
- **Example**: Change from 0.75 → 0.5 Å

### How to Re-cluster with Different RMSD

```bash
# Edit config.txt
# ClusterRMSDCutoff = 1.0  # was 0.75

# Re-run clustering only
python scripts/run_docking_workflow.py config.txt --skip-create-table --skip-docking
```

---

## 2. Assessing Sampling Adequacy

### Metric 1: Cluster Count

| Cluster Count | Assessment | Action |
|--------------|------------|--------|
| **< 10** | ⚠️ Insufficient sampling | Run more conformers OR relax filters |
| **10-50** | ✓ Reasonable | Good for focused binding site |
| **50-200** | ✓ Excellent | Very thorough sampling |
| **> 200** | ⚠️ Possible over-clustering | Consider increasing RMSD cutoff |

### Metric 2: Cluster Saturation

**Goal**: See diminishing returns as you add conformers

```python
# Good saturation example:
# After 1000 conformers: 50 clusters
# After 2000 conformers: 65 clusters  ← Saturation
# After 3000 conformers: 70 clusters  ← Saturated!

# Under-sampled example:
# After 1000 conformers: 20 clusters
# After 2000 conformers: 45 clusters  ← Still growing
# After 3000 conformers: 75 clusters  ← Need more!
```

**How to check**:
1. Run initial campaign (e.g., 1000 conformers)
2. Count clusters: `ls $OutputDir/clustered_final/*.pdb | wc -l`
3. Run additional conformers (e.g., 1000 more)
4. Re-cluster and compare cluster count
5. If cluster count plateaus → you're saturated! ✓

### Metric 3: Pass Rate Analysis

Check your H-bond geometry summary CSVs:

```bash
# Count total conformers
grep -h "," hbond_geometry_summary*.csv | wc -l

# Count passing conformers
grep -h ",True,True," hbond_geometry_summary*.csv | wc -l
```

| Pass Rate | Assessment | Action |
|-----------|------------|--------|
| **< 0.1%** | ⚠️ Too strict | Relax filters (MaxScore, H-bond params) |
| **0.1-1%** | ✓ Selective | Good for finding best binders |
| **1-10%** | ✓ Reasonable | Balanced sampling |
| **> 10%** | ⚠️ Too permissive | May want stricter filters |

### Metric 4: Score Distribution

**Ideal score distribution**:
- Broad distribution initially
- Convergence to a "funnel" at lower scores
- Top clusters should be in the bottom 10% of scores

**Red flags**:
- All scores very similar → May be stuck in local minimum
- Wide spread with no clustering → Under-sampled
- Bimodal distribution → May have multiple binding modes (good!)

---

## 3. Practical Examples

### Example 1: Small Rigid Ligand (10 heavy atoms)

**Settings**:
```ini
[grade_conformers]
ClusterRMSDCutoff = 0.5
ArrayTaskCount = 5
```

**Expected results**:
- ~20-50 clusters if well-sampled
- High pass rates (>1%) due to rigidity
- Tight score convergence

**Assessment**: If you get <10 clusters → run more arrays

---

### Example 2: Medium Flexible Ligand (18 heavy atoms)

**Settings**:
```ini
[grade_conformers]
ClusterRMSDCutoff = 0.75
ArrayTaskCount = 10
```

**Expected results**:
- ~50-100 clusters if well-sampled
- Moderate pass rates (0.5-2%)
- Broader score distribution

**Assessment**: Check redundancy (poses/cluster). If >10 → likely saturated

---

### Example 3: Large Flexible Ligand (30 heavy atoms)

**Settings**:
```ini
[grade_conformers]
ClusterRMSDCutoff = 1.0
ArrayTaskCount = 20
```

**Expected results**:
- ~100-200 clusters if well-sampled
- Lower pass rates (0.1-1%)
- Wide score distribution

**Assessment**: May need iterative refinement of binding region

---

## 4. Workflow for Determining Adequacy

### Step 1: Initial Run
```bash
# Start with ArrayTaskCount = 5-10
bash scripts/submit_complete_workflow.sh config.txt
```

### Step 2: Analyze Results
```bash
# Run analysis
python scripts/analyze_sampling_adequacy.py config.txt --ligand-heavy-atoms 15

# Count clusters
ls $OutputDir/clustered_final/*.pdb | wc -l

# Check top scores
sort -t',' -k4 -n $OutputDir/clustered_final/cluster_summary.csv | head -20
```

### Step 3: Assess and Iterate

**If under-sampled** (<20 clusters):
```bash
# Increase ArrayTaskCount in config.txt
# ArrayTaskCount = 20  # was 10

# Run more arrays (they'll add to existing results)
bash scripts/submit_complete_workflow.sh config.txt
```

**If RMSD seems wrong**:
```bash
# Edit config.txt: adjust ClusterRMSDCutoff
# Re-cluster existing results
python scripts/run_docking_workflow.py config.txt --skip-create-table --skip-docking
```

**If filters too strict** (pass rate <0.1%):
```ini
# Edit config.txt
[grade_conformers]
MaxScore = -250  # was -300 (less strict)
HBondDistanceIdealBuffer = 1.0  # was 0.8 (more permissive)
MaxPerturbTries = 10  # was 5 (more attempts)
```

### Step 4: Visual Inspection

**Critical final step**: Load top 10-20 clusters in PyMOL/Chimera

```pymol
# In PyMOL
cd $OutputDir/clustered_final/
load cluster_0001_*.pdb
load cluster_0002_*.pdb
...
# Examine if they represent distinct binding modes
```

**Good signs**:
- ✓ Distinct binding modes with different key interactions
- ✓ Chemically reasonable poses
- ✓ H-bonds to expected residues

**Red flags**:
- ⚠️ All poses look nearly identical → Over-sampled or too tight RMSD
- ⚠️ Poses look random/unrealistic → Filters too loose
- ⚠️ Missing expected interactions → Under-sampled

---

## 5. RMSD Cutoff Philosophy

### Conceptual Framework

**RMSD cutoff represents**: "How different do two poses need to be to represent distinct binding hypotheses?"

**Consider**:
- **Binding site volume**: Larger sites → larger RMSD cutoff
- **Ligand flexibility**: More rotatable bonds → larger RMSD cutoff
- **Experimental context**:
  - Lead optimization (subtle differences matter) → smaller cutoff
  - Hit discovery (broad diversity) → larger cutoff
- **Downstream analysis**:
  - MD simulations planned → keep more clusters (smaller cutoff)
  - Quick screening → fewer clusters (larger cutoff)

### Literature Values

- **Fragment screening**: 0.5-1.0 Å (small, rigid)
- **Drug-like molecules**: 0.75-1.5 Å (standard)
- **Peptides/large ligands**: 1.5-2.5 Å (flexible)
- **Covalent docking**: 0.5-0.75 Å (precise geometry)

---

## 6. Common Scenarios

### Scenario: "I have 500 clusters!"

**Diagnosis**: RMSD cutoff too tight OR extremely diverse sampling

**Solutions**:
1. Try ClusterRMSDCutoff = 1.5 Å (from 0.75)
2. Take top 50 by score and re-cluster those separately
3. This might actually be good! Check if binding site is large/promiscuous

### Scenario: "I only have 5 clusters"

**Diagnosis**: Under-sampled OR RMSD cutoff too loose OR very restrictive filters

**Solutions**:
1. Check pass rates - if <0.1%, relax filters
2. Run more array tasks (ArrayTaskCount = 20+)
3. Try ClusterRMSDCutoff = 0.5 Å (from 0.75)
4. Check input SDF has diverse conformers

### Scenario: "Top 10 clusters all have same score (~within 2 REU)"

**Diagnosis**: Good sampling convergence! Multiple degenerate solutions.

**Action**:
- ✓ This is actually ideal
- Visually inspect to see if they're truly different binding modes
- Consider all of them for downstream analysis

### Scenario: "Scores range from -400 to -250 REU across clusters"

**Diagnosis**: Very diverse sampling, possibly multiple binding sites/modes

**Action**:
- Visually inspect top vs bottom clusters
- May indicate promiscuous binding
- Focus analysis on top scoring clusters (bottom 20%)

---

## 7. Quick Decision Tree

```
START: How many clusters do you have?
│
├─ < 10 clusters
│   ├─ Pass rate < 0.1%? → RELAX FILTERS
│   └─ Pass rate OK? → RUN MORE CONFORMERS
│
├─ 10-50 clusters
│   ├─ Visually distinct? → ✓ GOOD SAMPLING
│   └─ Look similar? → INCREASE RMSD to 1.0-1.5
│
├─ 50-200 clusters
│   └─ ✓ EXCELLENT SAMPLING
│       └─ For simpler analysis: INCREASE RMSD to 1.0-1.5
│
└─ > 200 clusters
    ├─ Many clusters with 1 member? → INCREASE RMSD
    └─ Highly diverse binding? → ✓ GOOD, but consider focusing on top 50
```

---

## 8. Summary Recommendations

### For Your Situation (KYNA-like, 0.75 Å RMSD)

**0.75 Å is a good starting point for medium-sized ligands!**

**Adjust based on**:
1. **Ligand size**: If >20 heavy atoms → try 1.0 Å
2. **Results**: If >100 clusters → try 1.0 Å; if <10 clusters → try 0.5 Å
3. **Visual inspection**: Most important! Do clusters look distinct?

**Target cluster count**: 20-100 clusters is ideal for most applications

**Run the analysis tool**:
```bash
python scripts/analyze_sampling_adequacy.py config.txt --ligand-heavy-atoms <your_atom_count>
```

This will give you specific recommendations for YOUR system!

---

## Further Reading

- Rosetta Ligand Docking documentation
- RMSD clustering in drug discovery (J. Chem. Inf. Model.)
- Your own visual inspection (most valuable!)
