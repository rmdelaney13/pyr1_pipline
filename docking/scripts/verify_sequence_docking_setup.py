#!/usr/bin/env python
"""
Verification script for sequence-CSV docking setup.

Checks that all required files are present and configuration is valid.

Usage:
    python verify_sequence_docking_setup.py config.txt
"""

import argparse
import os
import sys
from configparser import ConfigParser

import pandas as pd


def check_file_exists(path, description):
    """Check if file exists and print status."""
    exists = os.path.exists(path)
    status = "✓" if exists else "✗"
    print(f"  {status} {description}: {path}")
    return exists


def verify_config(config_path):
    """Verify configuration file."""
    print("\n" + "="*80)
    print("CONFIGURATION VERIFICATION")
    print("="*80)

    if not check_file_exists(config_path, "Config file"):
        return False

    config = ConfigParser()
    try:
        with open(config_path, "r", encoding="utf-8-sig") as f:
            config.read_file(f)
        print("  ✓ Config file parseable")
    except Exception as e:
        print(f"  ✗ Config file parse error: {e}")
        return False

    # Check required sections
    required_sections = ["DEFAULT", "SEQUENCE_DOCKING", "grade_conformers"]
    all_present = True
    for section in required_sections:
        if section in config:
            print(f"  ✓ Section [{section}]")
        else:
            print(f"  ✗ Missing section [{section}]")
            all_present = False

    if not all_present:
        return False

    # Check SEQUENCE_DOCKING keys
    seq_section = config["SEQUENCE_DOCKING"]
    csv_path = seq_section.get("SequenceCSV", None)
    positions = seq_section.get("ThreadingPositions", None)

    if csv_path:
        print(f"  ✓ SequenceCSV configured: {csv_path}")
    else:
        print(f"  ✗ SequenceCSV not configured")
        return False

    if positions:
        pos_list = positions.split()
        print(f"  ✓ ThreadingPositions configured ({len(pos_list)} positions)")
    else:
        print(f"  ✗ ThreadingPositions not configured")
        return False

    return True


def verify_csv(csv_path):
    """Verify sequence CSV."""
    print("\n" + "="*80)
    print("SEQUENCE CSV VERIFICATION")
    print("="*80)

    if not check_file_exists(csv_path, "Sequence CSV"):
        return False

    try:
        df = pd.read_csv(csv_path)
        print(f"  ✓ CSV readable ({len(df)} sequences)")
    except Exception as e:
        print(f"  ✗ CSV read error: {e}")
        return False

    # Check columns
    required_cols = {"name", "signature"}
    if required_cols.issubset(set(df.columns)):
        print(f"  ✓ Required columns present: {required_cols}")
    else:
        print(f"  ✗ Missing columns. Found: {df.columns.tolist()}, need: {required_cols}")
        return False

    # Check for duplicates
    if df["name"].duplicated().any():
        dups = df[df["name"].duplicated()]["name"].tolist()
        print(f"  ✗ Duplicate sequence names: {dups}")
        return False
    else:
        print(f"  ✓ No duplicate names")

    # Check signature lengths
    sig_lengths = df["signature"].str.len()
    if sig_lengths.nunique() == 1:
        print(f"  ✓ All signatures same length ({sig_lengths.iloc[0]} AA)")
    else:
        print(f"  ✗ Inconsistent signature lengths: {sig_lengths.unique()}")
        return False

    # Check amino acid validity
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_seqs = []
    for idx, row in df.iterrows():
        sig = str(row["signature"]).upper()
        invalid_aa = set(sig) - valid_aa
        if invalid_aa:
            invalid_seqs.append(f"{row['name']} (invalid: {invalid_aa})")

    if invalid_seqs:
        print(f"  ✗ Invalid amino acids found:")
        for s in invalid_seqs[:5]:  # Show first 5
            print(f"      {s}")
        return False
    else:
        print(f"  ✓ All amino acids valid")

    # Preview sequences
    print(f"\n  Sequence preview:")
    for idx, row in df.head(3).iterrows():
        print(f"    [{idx}] {row['name']:20s} {row['signature']}")
    if len(df) > 3:
        print(f"    ... and {len(df) - 3} more")

    return True


def verify_scripts():
    """Verify required scripts are present."""
    print("\n" + "="*80)
    print("SCRIPT VERIFICATION")
    print("="*80)

    script_dir = os.path.dirname(os.path.abspath(__file__))

    scripts = {
        "Docking script": "grade_conformers_sequence_csv_docking_multiple_slurm.py",
        "Submission helper": "submit_sequence_csv_docking.py",
        "Utility functions": "docking_pipeline_utils.py",
    }

    all_present = True
    for description, filename in scripts.items():
        path = os.path.join(script_dir, filename)
        if not check_file_exists(path, description):
            all_present = False

    return all_present


def verify_legacy_modules():
    """Verify legacy modules are accessible."""
    print("\n" + "="*80)
    print("LEGACY MODULE VERIFICATION")
    print("="*80)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    legacy_dir = os.path.normpath(os.path.join(script_dir, "..", "legacy"))

    legacy_modules = [
        "alignment.py",
        "conformer_prep.py",
        "collision_check.py",
    ]

    all_present = True
    for module in legacy_modules:
        path = os.path.join(legacy_dir, module)
        if not check_file_exists(path, f"Legacy module: {module}"):
            all_present = False

    return all_present


def verify_templates():
    """Verify template files are present."""
    print("\n" + "="*80)
    print("TEMPLATE VERIFICATION")
    print("="*80)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    template_dir = os.path.normpath(os.path.join(script_dir, "..", "templates"))

    templates = {
        "Config template": "config_sequence_docking.txt",
        "Example CSV": "sequences_example.csv",
    }

    all_present = True
    for description, filename in templates.items():
        path = os.path.join(template_dir, filename)
        if not check_file_exists(path, description):
            all_present = False

    return all_present


def verify_pdb_files(config):
    """Verify PDB template files."""
    print("\n" + "="*80)
    print("PDB TEMPLATE VERIFICATION")
    print("="*80)

    default = config["DEFAULT"]
    pre_pdb = default.get("PrePDBFileName", "")
    post_pdb = default.get("PostPDBFileName", "")

    all_present = True
    if not check_file_exists(pre_pdb, "Pre-PDB template"):
        all_present = False
    if not check_file_exists(post_pdb, "Post-PDB template"):
        all_present = False

    return all_present


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument(
        "config_file",
        help="Path to config file to verify"
    )
    parser.add_argument(
        "--skip-csv",
        action="store_true",
        help="Skip CSV verification (useful if CSV doesn't exist yet)"
    )

    args = parser.parse_args()

    print("="*80)
    print("SEQUENCE-CSV DOCKING SETUP VERIFICATION")
    print("="*80)
    print(f"Config: {args.config_file}\n")

    results = {}

    # Verify config
    results["config"] = verify_config(args.config_file)

    if results["config"]:
        config = ConfigParser()
        with open(args.config_file, "r", encoding="utf-8-sig") as f:
            config.read_file(f)

        # Verify CSV
        if not args.skip_csv:
            seq_section = config["SEQUENCE_DOCKING"]
            csv_path = seq_section.get("SequenceCSV", None)
            if csv_path:
                results["csv"] = verify_csv(csv_path)
            else:
                results["csv"] = False
        else:
            print("\n" + "="*80)
            print("SEQUENCE CSV VERIFICATION")
            print("="*80)
            print("  ⊘ Skipped (--skip-csv)")
            results["csv"] = True

        # Verify PDB files
        results["pdbs"] = verify_pdb_files(config)
    else:
        results["csv"] = False
        results["pdbs"] = False

    # Verify scripts
    results["scripts"] = verify_scripts()

    # Verify legacy modules
    results["legacy"] = verify_legacy_modules()

    # Verify templates
    results["templates"] = verify_templates()

    # Summary
    print("\n" + "="*80)
    print("VERIFICATION SUMMARY")
    print("="*80)

    all_passed = all(results.values())

    for component, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status:8s} {component.upper()}")

    print("="*80)

    if all_passed:
        print("\n✓ ALL CHECKS PASSED")
        print("\nYou're ready to run sequence-CSV docking!")
        print("\nNext steps:")
        print("  1. Test locally: python scripts/submit_sequence_csv_docking.py config.txt --local --max-local 1")
        print("  2. Submit to SLURM: python scripts/submit_sequence_csv_docking.py config.txt --submit")
        return 0
    else:
        print("\n✗ SOME CHECKS FAILED")
        print("\nPlease fix the issues above before running docking.")
        print("\nCommon fixes:")
        print("  - Missing CSV: Create a CSV with 'name' and 'signature' columns")
        print("  - Wrong paths: Update paths in config to match your system")
        print("  - Missing templates: Re-copy from docking/templates/")
        return 1


if __name__ == "__main__":
    sys.exit(main())
