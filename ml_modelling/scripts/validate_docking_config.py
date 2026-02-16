#!/usr/bin/env python
"""
Validate mutant docking configuration before running pipeline.

Usage:
  python ml_modelling/scripts/validate_docking_config.py \
      --config cache/test_001/docking/docking_config.txt \
      --mutant-pdb cache/test_001/mutant.pdb

Checks:
  1. Config has ChainLetter and ResidueNumber
  2. ChainLetter is NOT 'X' (ligand chain)
  3. Mutant PDB contains the specified chain and residue
  4. Mutant PDB has no ligand (expected for mutant docking)
  5. Mutant PDB has waters (needed for H-bond filtering)

Author: Claude Code (Whitehead Lab PYR1 Pipeline)
Date: 2026-02-16
"""

import argparse
import sys
from pathlib import Path
from configparser import ConfigParser


def parse_pdb_residues(pdb_path):
    """
    Extract unique (chain, resnum) pairs from PDB.

    Returns:
        dict: {chain: set(residue_numbers)}
    """
    chains = {}
    waters = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21:22].strip()
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue

                resname = line[17:20].strip()

                if resname in {'WAT', 'HOH', 'TP3'}:
                    waters.append((chain, resnum))
                else:
                    if chain not in chains:
                        chains[chain] = set()
                    chains[chain].add(resnum)

    return chains, waters


def validate_config(config_path, mutant_pdb_path):
    """
    Validate docking configuration.

    Returns:
        (bool, list): (is_valid, list_of_errors)
    """
    errors = []
    warnings = []

    # Load config
    config = ConfigParser()
    config.read(config_path)

    if 'DEFAULT' not in config:
        errors.append("Config missing [DEFAULT] section")
        return False, errors

    default = config['DEFAULT']

    # Check ChainLetter
    if 'ChainLetter' not in default:
        errors.append("Config missing ChainLetter parameter")
    else:
        chain_letter = default['ChainLetter'].strip()
        if chain_letter == 'X':
            errors.append(
                f"ChainLetter = {chain_letter} (ligand chain) - this will fail! "
                "Mutant PDB has no ligand. Use a protein chain like 'A'."
            )
        print(f"✓ ChainLetter = {chain_letter}")

    # Check ResidueNumber
    if 'ResidueNumber' not in default:
        errors.append("Config missing ResidueNumber parameter")
    else:
        try:
            residue_number = int(default['ResidueNumber'])
            print(f"✓ ResidueNumber = {residue_number}")
        except ValueError:
            errors.append(f"Invalid ResidueNumber: {default['ResidueNumber']}")
            residue_number = None

    # Check LigandResidueNumber
    if 'LigandResidueNumber' not in default:
        warnings.append("Config missing LigandResidueNumber (recommended: 1)")
    else:
        lig_res_num = default['LigandResidueNumber']
        print(f"✓ LigandResidueNumber = {lig_res_num}")

    # Check AutoGenerateAlignment
    if 'AutoGenerateAlignment' not in default:
        warnings.append("Config missing AutoGenerateAlignment (recommended: False)")
    else:
        auto_align = default['AutoGenerateAlignment'].strip().lower()
        if auto_align in {'true', '1', 'yes'}:
            warnings.append(
                "AutoGenerateAlignment = True - this may be slow. "
                "If alignment table already exists, set to False."
            )
        print(f"✓ AutoGenerateAlignment = {auto_align}")

    # Validate against mutant PDB
    if not Path(mutant_pdb_path).exists():
        errors.append(f"Mutant PDB not found: {mutant_pdb_path}")
        return (len(errors) == 0), errors

    chains, waters = parse_pdb_residues(mutant_pdb_path)

    print(f"\nMutant PDB structure:")
    print(f"  Chains: {sorted(chains.keys())}")
    for chain in sorted(chains.keys()):
        resnums = sorted(chains[chain])
        print(f"    Chain {chain}: {len(resnums)} residues (range {min(resnums)}-{max(resnums)})")
    print(f"  Waters: {len(waters)}")

    # Check for ligand (should NOT exist)
    ligand_chains = [c for c in chains if c in {'X', 'Y', 'Z', 'L'}]
    if ligand_chains:
        warnings.append(
            f"Found potential ligand chain(s): {ligand_chains}. "
            "For mutant docking, ligand should be absent."
        )
    else:
        print(f"✓ No ligand chain detected (expected for mutant docking)")

    # Check reference residue exists
    if 'ChainLetter' in default and residue_number is not None:
        chain_letter = default['ChainLetter'].strip()

        if chain_letter not in chains:
            errors.append(
                f"Chain {chain_letter} not found in mutant PDB! "
                f"Available chains: {sorted(chains.keys())}"
            )
        elif residue_number not in chains[chain_letter]:
            errors.append(
                f"Residue {chain_letter}:{residue_number} not found in mutant PDB! "
                f"Chain {chain_letter} has residues: {sorted(chains[chain_letter])[:10]}..."
            )
        else:
            print(f"✓ Reference residue {chain_letter}:{residue_number} exists in mutant PDB")

    # Check waters (needed for H-bond filtering)
    if len(waters) == 0:
        warnings.append(
            "No waters found in mutant PDB. "
            "H-bond geometry filtering requires waters!"
        )
    else:
        print(f"✓ {len(waters)} waters detected (needed for H-bond filtering)")

    # Print warnings
    if warnings:
        print("\n⚠ WARNINGS:")
        for w in warnings:
            print(f"  - {w}")

    # Return validation result
    return (len(errors) == 0), errors


def main():
    parser = argparse.ArgumentParser(
        description='Validate mutant docking configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--config', required=True, help='Path to docking_config.txt')
    parser.add_argument('--mutant-pdb', required=True, help='Path to mutant.pdb')

    args = parser.parse_args()

    print("="*60)
    print("VALIDATING MUTANT DOCKING CONFIGURATION")
    print("="*60)
    print(f"Config: {args.config}")
    print(f"Mutant PDB: {args.mutant_pdb}")
    print()

    is_valid, errors = validate_config(args.config, args.mutant_pdb)

    print("\n" + "="*60)
    if is_valid:
        print("✅ VALIDATION PASSED")
        print("="*60)
        print("Configuration is valid. Ready to run mutant docking.")
        return 0
    else:
        print("❌ VALIDATION FAILED")
        print("="*60)
        print("ERRORS:")
        for e in errors:
            print(f"  ✗ {e}")
        print()
        print("Fix these errors before running docking.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
