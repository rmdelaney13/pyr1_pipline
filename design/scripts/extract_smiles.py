#!/usr/bin/env python3
"""
Extract SMILES string from SDF file using RDKit.

Usage:
    python extract_smiles.py input.sdf
    python extract_smiles.py input.sdf --output smiles.txt
"""

import argparse
import sys

try:
    from rdkit import Chem
except ImportError:
    print("ERROR: RDKit is not installed", file=sys.stderr)
    print("Install with: conda install -c conda-forge rdkit", file=sys.stderr)
    sys.exit(1)


def extract_smiles(sdf_path, canonical=True):
    """
    Extract SMILES from SDF file.

    Args:
        sdf_path: Path to SDF file
        canonical: If True, return canonical SMILES

    Returns:
        SMILES string or None if extraction failed
    """
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path))
        mol = next(iter(supplier))

        if mol is None:
            print(f"ERROR: Could not read molecule from {sdf_path}", file=sys.stderr)
            return None

        if canonical:
            smiles = Chem.MolToSmiles(mol)
        else:
            smiles = Chem.MolToSmiles(mol, canonical=False)

        return smiles

    except Exception as e:
        print(f"ERROR: Failed to extract SMILES: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Extract SMILES from SDF file"
    )
    parser.add_argument(
        'sdf_file',
        help='Input SDF file'
    )
    parser.add_argument(
        '--output', '-o',
        help='Output file (default: print to stdout)'
    )
    parser.add_argument(
        '--non-canonical',
        action='store_true',
        help='Use non-canonical SMILES'
    )

    args = parser.parse_args()

    smiles = extract_smiles(args.sdf_file, canonical=not args.non_canonical)

    if smiles is None:
        sys.exit(1)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(smiles + '\n')
        print(f"SMILES written to {args.output}")
    else:
        print(smiles)


if __name__ == '__main__':
    main()
