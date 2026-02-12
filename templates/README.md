# Pipeline Config Templates

This directory contains configuration file templates for the pyr1_pipeline.

---

## 🌟 Recommended: Unified Config Template

**File:** [unified_config_template.txt](unified_config_template.txt)

**Use this for:** Complete docking + design workflows

**What it includes:**
- ✅ All docking parameters (`[DEFAULT]`, `[create_table]`, `[grade_conformers]`)
- ✅ All design parameters (`[design]`)
- ✅ Fully commented with examples
- ✅ Ready to copy and edit

**Quick start:**
```bash
# Copy to your campaign directory
cp templates/unified_config_template.txt /scratch/youruser/your_ligand/config.txt

# Edit the key sections
vim /scratch/youruser/your_ligand/config.txt

# Run pipelines
python docking/scripts/run_docking_workflow.py config.txt
python design/scripts/run_design_pipeline.py config.txt
```

**See also:** [CONFIG_GUIDE.md](CONFIG_GUIDE.md) for detailed parameter explanations

---

## 📁 Other Templates

### **docking/templates/config.txt**

**Use this for:** Docking-only workflows (no design)

Includes:
- `[DEFAULT]` - Root paths
- `[create_table]` - Ligand input
- `[grade_conformers]` - Docking parameters

### **design/CONFIG_TEMPLATE.txt**

**Use this for:** Adding design to existing docking config

Shows only the `[design]` section for reference. Recommended to use unified template instead.

---

## 📚 Documentation

| File | Description |
|------|-------------|
| [unified_config_template.txt](unified_config_template.txt) | Complete config template (docking + design) |
| [CONFIG_GUIDE.md](CONFIG_GUIDE.md) | Parameter explanations and examples |
| [../docking/WORKFLOW_README.md](../docking/WORKFLOW_README.md) | Docking pipeline documentation |
| [../design/DESIGN_PIPELINE_README.md](../design/DESIGN_PIPELINE_README.md) | Design pipeline documentation |

---

## 🎯 Choosing the Right Template

```
┌─────────────────────────────────────────────────┐
│  What do you want to do?                       │
└─────────────────────────────────────────────────┘
                    │
        ┌───────────┴───────────┐
        │                       │
        ▼                       ▼
   Docking + Design        Docking Only
        │                       │
        ▼                       ▼
unified_config_template    docking/templates/config.txt
```

**Most users:** Use `unified_config_template.txt`

---

## ⚙️ Config File Structure

A complete config file has these sections:

```ini
[DEFAULT]
# Root paths and protein templates

[MULTIPLE_PASSES]
# Multi-pass docking control (optional)

[create_table]
# Ligand input files and alignment

[grade_conformers]
# Docking parameters and SLURM settings

[design]
# Design pipeline: MPNN, Rosetta, filtering, AF3
```

---

## 🔧 Customization

### For a New Ligand

1. Copy `unified_config_template.txt`
2. Edit these sections:
   - `CAMPAIGN_ROOT` and `SCRATCH_ROOT` in `[DEFAULT]`
   - `MoleculeSDFs` in `[create_table]`
   - `TargetAtomTriplets` in `[create_table]`
   - `LigandParams` and `LigandSDF` in `[design]`

### For a New Protein Target

1. Copy `unified_config_template.txt`
2. Edit these additional sections:
   - `PrePDBFileName` and `PostPDBFileName` in `[DEFAULT]`
   - `DesignResidues` in `[design]`
   - `GlycineShavePositions` in `[grade_conformers]`
   - Optionally create new MPNN bias/omit files

---

## 📋 Quick Reference

### Essential Parameters

| Section | Parameter | What It Does |
|---------|-----------|--------------|
| `[DEFAULT]` | `CAMPAIGN_ROOT` | Where input files live |
| `[DEFAULT]` | `SCRATCH_ROOT` | Where outputs go |
| `[create_table]` | `MoleculeSDFs` | Input ligand conformers |
| `[grade_conformers]` | `ArrayTaskCount` | SLURM parallelization |
| `[grade_conformers]` | `MaxScore` | Energy cutoff |
| `[design]` | `DesignResidues` | Residues to redesign |
| `[design]` | `FilterTargetN` | How many designs to keep |

### Common Adjustments

**More permissive docking:**
```ini
[grade_conformers]
MaxScore = -250
```

**More designs in output:**
```ini
[design]
FilterTargetN = 2000
FilterMaxUnsats = 2
```

**Two design iterations:**
```ini
[design]
DesignIterationRounds = 2
```

---

## 🧪 Testing Your Config

Before running the full pipeline:

```bash
# Test config parsing
python -c "
from configparser import ConfigParser
cfg = ConfigParser()
cfg.read('config.txt')
print('Sections:', cfg.sections())
print('OutputDir:', cfg.get('grade_conformers', 'OutputDir'))
"

# Test docking dry-run
python docking/scripts/run_docking_workflow.py config.txt --prepare-only

# Test design dry-run
python design/scripts/run_design_pipeline.py config.txt --dry-run
```

---

## 📞 Support

- **Config questions:** See [CONFIG_GUIDE.md](CONFIG_GUIDE.md)
- **Docking issues:** See [../docking/WORKFLOW_README.md](../docking/WORKFLOW_README.md)
- **Design issues:** See [../design/DESIGN_PIPELINE_README.md](../design/DESIGN_PIPELINE_README.md)
- **Setup help:** See [../design/SETUP_CHECKLIST.md](../design/SETUP_CHECKLIST.md)
