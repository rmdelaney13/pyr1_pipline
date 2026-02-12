#!/usr/bin/env python3
"""
Update SMILES in AF3 JSON templates.

This script finds the ligand SMILES field in AF3 JSON templates and updates it.

Usage:
    # Update from SDF
    python update_template_smiles.py template.json --sdf ligand.sdf --output updated_template.json

    # Update from SMILES string
    python update_template_smiles.py template.json --smiles "C1=CC2=C(C(=C1)O)NC(=CC2=O)C(=O)O"

    # Update in-place
    python update_template_smiles.py template.json --sdf ligand.sdf --in-place
"""

import argparse
import json
import sys
from pathlib import Path


def extract_smiles_from_sdf(sdf_path):
    """Extract SMILES from SDF using RDKit"""
    try:
        from rdkit import Chem
        supplier = Chem.SDMolSupplier(str(sdf_path))
        mol = next(iter(supplier))
        if mol:
            return Chem.MolToSmiles(mol)
        return None
    except ImportError:
        print("ERROR: RDKit not installed. Install with: conda install -c conda-forge rdkit",
              file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to extract SMILES from {sdf_path}: {e}", file=sys.stderr)
        return None


def update_smiles_in_template(template_path, smiles, output_path=None):
    """
    Update SMILES in AF3 JSON template.

    Args:
        template_path: Path to template JSON
        smiles: SMILES string to insert
        output_path: Output path (None = print to stdout)

    Returns:
        True if SMILES was updated, False otherwise
    """
    with open(template_path, 'r') as f:
        data = json.load(f)

    # Find ligand SMILES in sequences
    updated = False
    locations = []

    if 'sequences' in data:
        for i, seq in enumerate(data['sequences']):
            if 'ligand' in seq:
                if 'smiles' in seq['ligand']:
                    old_smiles = seq['ligand']['smiles']
                    seq['ligand']['smiles'] = smiles
                    updated = True
                    locations.append(f"sequences[{i}].ligand.smiles")
                    print(f"Updated sequences[{i}].ligand.smiles")
                    print(f"  Old: {old_smiles}")
                    print(f"  New: {smiles}")

    if not updated:
        print("WARNING: No ligand SMILES field found in template", file=sys.stderr)
        print("Template structure:", file=sys.stderr)
        print(json.dumps(data, indent=2)[:500] + "...", file=sys.stderr)
        return False

    # Write output
    if output_path:
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"\n✓ Updated template saved to: {output_path}")
    else:
        print(json.dumps(data, indent=2))

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Update SMILES in AF3 JSON template",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        'template',
        help='Input template JSON file'
    )

    smiles_group = parser.add_mutually_exclusive_group(required=True)
    smiles_group.add_argument(
        '--smiles',
        help='SMILES string to use'
    )
    smiles_group.add_argument(
        '--sdf',
        help='Extract SMILES from SDF file'
    )

    parser.add_argument(
        '--output', '-o',
        help='Output file (default: print to stdout)'
    )
    parser.add_argument(
        '--in-place', '-i',
        action='store_true',
        help='Update template in-place'
    )

    args = parser.parse_args()

    # Get SMILES
    if args.sdf:
        print(f"Extracting SMILES from {args.sdf}...")
        smiles = extract_smiles_from_sdf(args.sdf)
        if not smiles:
            print("ERROR: Failed to extract SMILES", file=sys.stderr)
            sys.exit(1)
        print(f"  SMILES: {smiles}\n")
    else:
        smiles = args.smiles

    # Determine output path
    if args.in_place:
        output_path = args.template
    elif args.output:
        output_path = args.output
    else:
        output_path = None

    # Update template
    success = update_smiles_in_template(args.template, smiles, output_path)

    if not success:
        sys.exit(1)


if __name__ == '__main__':
    main()
