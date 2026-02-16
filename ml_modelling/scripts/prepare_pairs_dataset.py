#!/usr/bin/env python
"""
Prepare pairs dataset for ML pipeline orchestration.

This script takes input CSV(s) with ligands and variants and produces
a standardized pairs CSV for the orchestrator.

Input formats:
  1. Existing format: ligand_smiles_signature.csv (positives only)
  2. Separate lists: ligands.csv + variants.csv (cartesian product)
  3. Custom pairs: user-defined combinations

Output format (for orchestrator):
  pair_id, ligand_name, ligand_smiles, variant_name, variant_signature, label, label_tier, label_source

Usage:
  # Convert existing positives
  python prepare_pairs_dataset.py \
      --existing-csv ml_modelling/data/ligand_smiles_signature.csv \
      --output ml_modelling/pairs_dataset.csv

  # Generate cartesian product (all ligands × all variants)
  python prepare_pairs_dataset.py \
      --ligands-csv ligands.csv \
      --variants-csv variants.csv \
      --output ml_modelling/pairs_dataset.csv \
      --cartesian

  # Generate specific pairs (read from a mapping file)
  python prepare_pairs_dataset.py \
      --pairs-csv custom_pairs.csv \
      --output ml_modelling/pairs_dataset.csv

Author: Claude Code (Whitehead Lab PYR1 Pipeline)
Date: 2026-02-16
"""

import argparse
import hashlib
import logging
from pathlib import Path
from typing import List, Dict, Optional

import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def generate_pair_id(ligand_name: str, variant_name: str) -> str:
    """
    Generate unique pair ID from ligand and variant names.

    Format: {ligand_name}_{variant_name}_hash
    Hash ensures uniqueness even if names collide.
    """
    combined = f"{ligand_name}_{variant_name}"
    # Sanitize for filesystem (remove special characters)
    sanitized = combined.replace('/', '_').replace('\\', '_').replace(' ', '_')
    # Add short hash for uniqueness
    hash_suffix = hashlib.md5(combined.encode()).hexdigest()[:8]
    return f"{sanitized}_{hash_suffix}"


def normalize_variant_signature(signature: str) -> str:
    """
    Normalize variant signature to standard format.

    Handles two formats:
      1. Semicolon-separated: "59K;120A;160G"
      2. Underscore-separated: "K59Q_Y120A_A160G"

    Output: Always semicolon-separated, sorted by position: "59K;120A;160G"
    """
    # Handle NaN/None/empty
    if pd.isna(signature) or not signature:
        return ""

    signature = str(signature).strip()

    # Handle underscore format (e.g., "K59Q_Y120A_A160G")
    if '_' in signature:
        mutations = []
        for mut in signature.split('_'):
            mut = mut.strip()
            if not mut:
                continue
            # Extract position and new AA
            # Format: {WT_AA}{position}{new_AA}
            if len(mut) >= 3 and mut[0].isalpha() and mut[-1].isalpha():
                position = ''.join(c for c in mut[1:-1] if c.isdigit())
                new_aa = mut[-1]
                if position:
                    mutations.append(f"{position}{new_aa}")
                else:
                    logger.warning(f"Could not parse mutation (no position): {mut}")
            else:
                logger.warning(f"Could not parse mutation (invalid format): {mut}")
        signature = ';'.join(mutations)

    # Sort by position (numeric)
    mutations = [m.strip() for m in signature.split(';') if m.strip()]

    # Filter out invalid mutations (no position)
    valid_mutations = []
    for m in mutations:
        position_str = ''.join(c for c in m if c.isdigit())
        if position_str:
            valid_mutations.append(m)
        else:
            logger.warning(f"Skipping mutation with no position: {m}")

    if not valid_mutations:
        return ""

    # Sort by position
    sorted_mutations = sorted(valid_mutations, key=lambda m: int(''.join(c for c in m if c.isdigit())))
    return ';'.join(sorted_mutations)


def validate_smiles(smiles: str) -> bool:
    """Basic SMILES validation (checks for non-empty and reasonable length)."""
    if not smiles or len(smiles) < 3:
        return False
    # Additional validation could use RDKit if available
    return True


# ═══════════════════════════════════════════════════════════════════════
# CONVERSION FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

def convert_existing_csv(input_csv: Path, label: int = 1, label_tier: str = 'positive',
                        label_source: str = 'validated_binder') -> pd.DataFrame:
    """
    Convert existing ligand_smiles_signature.csv to orchestrator format.

    Expected input columns:
      - ligand_name
      - ligand_smiles_or_ligand_ID
      - PYR1_variant_name
      - PYR1_variant_signature

    Output columns:
      - pair_id
      - ligand_name
      - ligand_smiles
      - variant_name
      - variant_signature
      - label (1=binder, 0=non-binder)
      - label_tier (positive, T1, T2, T3)
      - label_source (provenance)
    """
    logger.info(f"Loading existing CSV: {input_csv}")
    df = pd.read_csv(input_csv)

    # Check required columns
    required = ['ligand_name', 'ligand_smiles_or_ligand_ID', 'PYR1_variant_name', 'PYR1_variant_signature']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    logger.info(f"Loaded {len(df)} rows")

    # Filter out rows with missing required data
    df = df.dropna(subset=['ligand_name', 'ligand_smiles_or_ligand_ID', 'PYR1_variant_name'])
    logger.info(f"After filtering NaN: {len(df)} rows")

    # Rename columns to standard format
    df_out = pd.DataFrame()
    df_out['ligand_name'] = df['ligand_name']
    df_out['ligand_smiles'] = df['ligand_smiles_or_ligand_ID']
    df_out['variant_name'] = df['PYR1_variant_name']
    df_out['variant_signature'] = df['PYR1_variant_signature'].apply(normalize_variant_signature)

    # Generate pair IDs
    df_out['pair_id'] = df_out.apply(
        lambda row: generate_pair_id(row['ligand_name'], row['variant_name']), axis=1
    )

    # Add labels
    df_out['label'] = label
    df_out['label_tier'] = label_tier
    df_out['label_source'] = label_source

    # Validate SMILES
    invalid_smiles = df_out[~df_out['ligand_smiles'].apply(validate_smiles)]
    if len(invalid_smiles) > 0:
        logger.warning(f"Found {len(invalid_smiles)} rows with invalid SMILES (will be included but may fail later)")

    # Reorder columns
    df_out = df_out[['pair_id', 'ligand_name', 'ligand_smiles', 'variant_name', 'variant_signature',
                     'label', 'label_tier', 'label_source']]

    logger.info(f"✓ Converted {len(df_out)} pairs")

    return df_out


def generate_cartesian_product(ligands_csv: Path, variants_csv: Path,
                               label: int = 1, label_tier: str = 'generated',
                               label_source: str = 'cartesian_product') -> pd.DataFrame:
    """
    Generate all combinations of ligands × variants (cartesian product).

    Ligands CSV format:
      - ligand_name
      - ligand_smiles

    Variants CSV format:
      - variant_name
      - variant_signature
      - (optional) label, label_tier
    """
    logger.info(f"Loading ligands from: {ligands_csv}")
    ligands = pd.read_csv(ligands_csv)

    logger.info(f"Loading variants from: {variants_csv}")
    variants = pd.read_csv(variants_csv)

    # Validate columns
    if not {'ligand_name', 'ligand_smiles'}.issubset(ligands.columns):
        raise ValueError("Ligands CSV must have columns: ligand_name, ligand_smiles")

    if not {'variant_name', 'variant_signature'}.issubset(variants.columns):
        raise ValueError("Variants CSV must have columns: variant_name, variant_signature")

    logger.info(f"Generating cartesian product: {len(ligands)} ligands × {len(variants)} variants = {len(ligands) * len(variants)} pairs")

    # Create cartesian product
    ligands['key'] = 1
    variants['key'] = 1
    df_out = ligands.merge(variants, on='key').drop('key', axis=1)

    # Normalize variant signatures
    df_out['variant_signature'] = df_out['variant_signature'].apply(normalize_variant_signature)

    # Generate pair IDs
    df_out['pair_id'] = df_out.apply(
        lambda row: generate_pair_id(row['ligand_name'], row['variant_name']), axis=1
    )

    # Add labels (use variant-specific labels if present, otherwise default)
    if 'label' not in df_out.columns:
        df_out['label'] = label
    if 'label_tier' not in df_out.columns:
        df_out['label_tier'] = label_tier
    if 'label_source' not in df_out.columns:
        df_out['label_source'] = label_source

    # Reorder columns
    df_out = df_out[['pair_id', 'ligand_name', 'ligand_smiles', 'variant_name', 'variant_signature',
                     'label', 'label_tier', 'label_source']]

    logger.info(f"✓ Generated {len(df_out)} pairs")

    return df_out


def load_custom_pairs(pairs_csv: Path) -> pd.DataFrame:
    """
    Load custom pairs CSV (already in orchestrator format or easily mappable).

    Required columns:
      - ligand_name, ligand_smiles, variant_name, variant_signature
    Optional columns:
      - pair_id (auto-generated if missing)
      - label, label_tier, label_source (defaults applied if missing)
    """
    logger.info(f"Loading custom pairs from: {pairs_csv}")
    df = pd.read_csv(pairs_csv)

    # Validate required columns
    required = ['ligand_name', 'ligand_smiles', 'variant_name', 'variant_signature']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Normalize variant signatures
    df['variant_signature'] = df['variant_signature'].apply(normalize_variant_signature)

    # Generate pair IDs if missing
    if 'pair_id' not in df.columns:
        df['pair_id'] = df.apply(
            lambda row: generate_pair_id(row['ligand_name'], row['variant_name']), axis=1
        )

    # Add default labels if missing
    if 'label' not in df.columns:
        df['label'] = 1
    if 'label_tier' not in df.columns:
        df['label_tier'] = 'unknown'
    if 'label_source' not in df.columns:
        df['label_source'] = 'custom_input'

    # Reorder columns
    df = df[['pair_id', 'ligand_name', 'ligand_smiles', 'variant_name', 'variant_signature',
             'label', 'label_tier', 'label_source']]

    logger.info(f"✓ Loaded {len(df)} pairs")

    return df


def merge_multiple_datasets(datasets: List[pd.DataFrame], remove_duplicates: bool = True) -> pd.DataFrame:
    """
    Merge multiple datasets (e.g., positives + negatives).

    Args:
        datasets: List of DataFrames in orchestrator format
        remove_duplicates: Remove duplicate pair_ids (keeps first occurrence)

    Returns:
        Merged DataFrame
    """
    logger.info(f"Merging {len(datasets)} datasets")

    df_merged = pd.concat(datasets, ignore_index=True)

    if remove_duplicates:
        before = len(df_merged)
        df_merged = df_merged.drop_duplicates(subset='pair_id', keep='first')
        after = len(df_merged)
        if before > after:
            logger.warning(f"Removed {before - after} duplicate pair_ids")

    logger.info(f"✓ Final dataset: {len(df_merged)} pairs")

    # Summary statistics
    logger.info(f"  Unique ligands: {df_merged['ligand_name'].nunique()}")
    logger.info(f"  Unique variants: {df_merged['variant_name'].nunique()}")
    logger.info(f"  Label distribution:")
    for label, count in df_merged['label'].value_counts().items():
        logger.info(f"    label={label}: {count} ({count/len(df_merged)*100:.1f}%)")

    if 'label_tier' in df_merged.columns:
        logger.info(f"  Tier distribution:")
        for tier, count in df_merged['label_tier'].value_counts().items():
            logger.info(f"    {tier}: {count} ({count/len(df_merged)*100:.1f}%)")

    return df_merged


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Prepare pairs dataset for ML pipeline orchestration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert existing positives
  python prepare_pairs_dataset.py \\
      --existing-csv ml_modelling/data/ligand_smiles_signature.csv \\
      --output ml_modelling/pairs_dataset.csv

  # Generate cartesian product
  python prepare_pairs_dataset.py \\
      --ligands-csv ligands.csv \\
      --variants-csv variants.csv \\
      --output ml_modelling/pairs_dataset.csv \\
      --cartesian

  # Merge positives + negatives
  python prepare_pairs_dataset.py \\
      --existing-csv positives.csv \\
      --additional-csv negatives.csv \\
      --output ml_modelling/pairs_dataset.csv
"""
    )

    # Input options (mutually exclusive groups)
    input_group = parser.add_argument_group('Input Options')
    input_group.add_argument('--existing-csv', type=Path,
                            help='Existing ligand_smiles_signature.csv (positives)')
    input_group.add_argument('--ligands-csv', type=Path,
                            help='Ligands CSV (for cartesian product)')
    input_group.add_argument('--variants-csv', type=Path,
                            help='Variants CSV (for cartesian product)')
    input_group.add_argument('--pairs-csv', type=Path,
                            help='Custom pairs CSV (already mapped)')
    input_group.add_argument('--additional-csv', type=Path, action='append',
                            help='Additional CSVs to merge (can be specified multiple times)')

    # Options
    parser.add_argument('--cartesian', action='store_true',
                       help='Generate cartesian product (requires --ligands-csv and --variants-csv)')
    parser.add_argument('--output', type=Path, required=True,
                       help='Output pairs CSV for orchestrator')
    parser.add_argument('--label', type=int, default=1,
                       help='Default label for pairs (1=binder, 0=non-binder)')
    parser.add_argument('--label-tier', type=str, default='positive',
                       help='Default label tier (positive, T1, T2, T3)')
    parser.add_argument('--label-source', type=str, default='validated_binder',
                       help='Default label source (provenance)')
    parser.add_argument('--remove-duplicates', action='store_true', default=True,
                       help='Remove duplicate pair_ids (default: True)')

    args = parser.parse_args()

    # Validate input options
    datasets = []

    if args.existing_csv:
        df = convert_existing_csv(
            args.existing_csv,
            label=args.label,
            label_tier=args.label_tier,
            label_source=args.label_source
        )
        datasets.append(df)

    if args.cartesian:
        if not (args.ligands_csv and args.variants_csv):
            parser.error("--cartesian requires both --ligands-csv and --variants-csv")
        df = generate_cartesian_product(
            args.ligands_csv,
            args.variants_csv,
            label=args.label,
            label_tier=args.label_tier,
            label_source=args.label_source
        )
        datasets.append(df)

    if args.pairs_csv:
        df = load_custom_pairs(args.pairs_csv)
        datasets.append(df)

    if args.additional_csv:
        for csv_path in args.additional_csv:
            logger.info(f"Loading additional CSV: {csv_path}")
            df = load_custom_pairs(csv_path)
            datasets.append(df)

    if not datasets:
        parser.error("No input data specified. Use --existing-csv, --cartesian, or --pairs-csv")

    # Merge datasets
    if len(datasets) > 1:
        df_final = merge_multiple_datasets(datasets, remove_duplicates=args.remove_duplicates)
    else:
        df_final = datasets[0]

    # Save output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df_final.to_csv(args.output, index=False)
    logger.info(f"✓ Saved {len(df_final)} pairs to: {args.output}")

    # Print first few rows
    print("\nFirst 5 rows:")
    print(df_final.head())

    print(f"\n[SUCCESS] Ready to run orchestrator:")
    print(f"  python ml_modelling/scripts/orchestrate_ml_dataset_pipeline.py \\")
    print(f"      --pairs-csv {args.output} \\")
    print(f"      --cache-dir <cache_dir> \\")
    print(f"      --template-pdb <template_pdb>")


if __name__ == '__main__':
    main()
