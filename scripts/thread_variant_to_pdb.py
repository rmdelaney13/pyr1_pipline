#!/usr/bin/env python3
"""
Thread PYR1 variant signature onto template PDB structure.

This script parses variant signatures in multiple formats and applies mutations
using PyRosetta to generate mutant structures for docking.

Supported formats:
    - Semicolon: "59K;120A;160G" (position + target AA)
    - Underscore: "K59Q_Y120A_A160G" (WT + position + target AA)
    - Space-separated: "59K 120A 160G"
    - Mixed format with underscores and numbers

Usage:
    python thread_variant_to_pdb.py \
        --template 3QN1_nolig_H2O.pdb \
        --signature "59K;120A;160G" \
        --output mutant_59K_120A_160G.pdb \
        --chain A

    # Batch processing from CSV
    python thread_variant_to_pdb.py \
        --template 3QN1_nolig_H2O.pdb \
        --csv variants.csv \
        --output-dir mutants/

Author: Claude Code (Whitehead Lab PYR1 Pipeline)
Date: 2026-02-16
"""

import argparse
import sys
import os
import re
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple, List
import logging

try:
    import pyrosetta
    from pyrosetta import pose_from_pdb
    from pyrosetta.toolbox import mutate_residue
    PYROSETTA_AVAILABLE = True
except ImportError:
    PYROSETTA_AVAILABLE = False
    print("WARNING: PyRosetta not available. Install to enable mutation threading.")


# Logging setup
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Amino acid mapping
AA_1_TO_3 = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

AA_3_TO_1 = {v: k for k, v in AA_1_TO_3.items()}


def parse_variant_signature(signature: str) -> Dict[int, str]:
    """
    Parse variant signature into {position: target_amino_acid} dict.

    Supports multiple formats:
        "59K;120A;160G" → {59: 'K', 120: 'A', 160: 'G'}
        "K59Q;Y120A;A160G" → {59: 'Q', 120: 'A', 160: 'G'}
        "K59Q_Y120A_A160G" → {59: 'Q', 120: 'A', 160: 'G'}
        "59K 120A 160G" → {59: 'K', 120: 'A', 160: 'G'}

    Args:
        signature: Variant signature string

    Returns:
        Dictionary mapping position (int) to target amino acid (single letter)

    Raises:
        ValueError: If signature format is unrecognized
    """
    if not signature or pd.isna(signature) or signature.strip() == '':
        return {}

    mutations = {}

    # Normalize: replace underscores and multiple spaces with single semicolon
    normalized = signature.replace('_', ';').replace(' ', ';')
    normalized = re.sub(r';+', ';', normalized)  # Remove duplicate semicolons

    for mut in normalized.split(';'):
        mut = mut.strip()
        if not mut:
            continue

        # Try format: "K59Q" or "59Q" (WT + position + target OR position + target)
        match = re.match(r'^([A-Z])?(\d+)([A-Z])$', mut)
        if match:
            wt_aa, pos, target_aa = match.groups()
            mutations[int(pos)] = target_aa
            continue

        # Try format: "K-59-Q" (hyphen-separated)
        match = re.match(r'^([A-Z])-(\d+)-([A-Z])$', mut)
        if match:
            wt_aa, pos, target_aa = match.groups()
            mutations[int(pos)] = target_aa
            continue

        # If no match, log warning and skip
        logger.warning(f"Could not parse mutation: '{mut}' in signature '{signature}'")

    if not mutations:
        raise ValueError(f"No valid mutations parsed from signature: '{signature}'")

    return mutations


def convert_wt_to_deletion_numbering(wt_position: int, deletion_start: int = 69, deletion_length: int = 2) -> int:
    """
    Convert WT PYR1 sequence numbering to deletion PDB numbering.

    NOTE: The 3QN1_nolig_H2O.pdb template PRESERVES original WT numbering
    (with a gap at 69-70), so this conversion is NOT needed for that template.
    Use --no-deletion when threading onto 3QN1_nolig_H2O.pdb.

    This function is only needed for PDBs that have been renumbered
    continuously after deletion (no gap in residue numbers).

    The 3QN1 structure has residues 69-70 missing (disordered loop).

    Conversion rule (for renumbered PDBs only):
        - Positions before deletion_start: unchanged
        - Positions within deletion: DO NOT EXIST (will raise error)
        - Positions after deletion: subtract deletion_length

    Args:
        wt_position: Position in WT PYR1 sequence (variant signature numbering)
        deletion_start: First deleted residue (default 67)
        deletion_length: Number of deleted residues (default 2)

    Returns:
        Corresponding position in deletion PDB

    Raises:
        ValueError: If position falls within deleted region (67-68)
    """
    deletion_end = deletion_start + deletion_length - 1  # 68

    if wt_position < deletion_start:
        # Before deletion: no change
        return wt_position
    elif deletion_start <= wt_position <= deletion_end:
        # Within deleted region: ERROR
        raise ValueError(
            f"Position {wt_position} is in the deleted region (residues {deletion_start}-{deletion_end} removed from PDB). "
            f"Cannot thread mutations at deleted positions!"
        )
    else:
        # After deletion: shift by deletion length
        return wt_position - deletion_length


def apply_mutations(
    pose,
    mutations: Dict[int, str],
    chain: str = 'A',
    validate: bool = True,
    has_deletion: bool = True,
    deletion_start: int = 69,
    deletion_length: int = 2
) -> Tuple[object, List[str]]:
    """
    Apply mutations to a PyRosetta pose.

    IMPORTANT: This function handles PDB numbering offsets due to deletions.
    The 3QN1 PYR1 structure has residues 67-68 DELETED, so positions after 68
    must be shifted by -2.

    Args:
        pose: PyRosetta Pose object
        mutations: {wt_position: target_aa} dictionary (WT sequence numbering!)
        chain: Chain ID (default 'A')
        validate: If True, check that positions exist and report current AAs
        has_deletion: If True, convert WT numbering to deletion PDB numbering
        deletion_start: First deleted residue (default 67)
        deletion_length: Number of deleted residues (default 2)

    Returns:
        (mutated_pose, log_messages): Tuple of mutated pose and list of log messages

    Raises:
        ValueError: If PyRosetta is not available or position in deleted region
        RuntimeError: If critical mutation fails
    """
    if not PYROSETTA_AVAILABLE:
        raise ValueError("PyRosetta is required for mutation threading")

    pdb_info = pose.pdb_info()
    log_messages = []

    for wt_position, target_aa in sorted(mutations.items()):
        # Convert WT numbering to deletion PDB numbering
        if has_deletion:
            try:
                pdb_position = convert_wt_to_deletion_numbering(
                    wt_position, deletion_start, deletion_length
                )
                conversion_msg = f"Position {wt_position} (WT) → {pdb_position} (PDB with Δ{deletion_start}-{deletion_start+deletion_length-1})"
                logger.info(conversion_msg)
                log_messages.append(conversion_msg)
            except ValueError as e:
                error_msg = f"ERROR: {e}"
                logger.error(error_msg)
                log_messages.append(error_msg)
                if validate:
                    raise
                continue
        else:
            pdb_position = wt_position

        # Convert PDB numbering to Rosetta pose numbering
        rosetta_position = pdb_info.pdb2pose(chain, pdb_position)

        if rosetta_position == 0:
            msg = f"WARNING: Position {pdb_position} (chain {chain}) not found in structure"
            logger.warning(msg)
            log_messages.append(msg)
            continue

        # Get current amino acid
        current_residue = pose.residue(rosetta_position)
        current_aa = current_residue.name1()

        if current_aa == target_aa:
            msg = f"Position {pdb_position} (Rosetta {rosetta_position}): already {target_aa} (no change)"
            logger.info(msg)
            log_messages.append(msg)
        else:
            msg = f"Position {pdb_position} (Rosetta {rosetta_position}): {current_aa} → {target_aa}"
            logger.info(msg)
            log_messages.append(msg)

            try:
                mutate_residue(pose, rosetta_position, target_aa)
            except Exception as e:
                error_msg = f"ERROR: Failed to mutate position {pdb_position}: {e}"
                logger.error(error_msg)
                log_messages.append(error_msg)
                if validate:
                    raise RuntimeError(error_msg)

    return pose, log_messages


def thread_variant(
    template_pdb: str,
    signature: str,
    output_pdb: str,
    chain: str = 'A',
    validate: bool = True,
    log_file: str = None,
    has_deletion: bool = True,
    deletion_start: int = 69,
    deletion_length: int = 2
) -> bool:
    """
    Complete workflow: parse signature, load template, apply mutations, save output.

    IMPORTANT: Default behavior assumes 3QN1 PYR1 structure with Δ67-68 deletion.
    If using a different template, set has_deletion=False.

    Args:
        template_pdb: Path to template PDB file
        signature: Variant signature string (uses WT sequence numbering!)
        output_pdb: Path to output mutated PDB file
        chain: Chain ID to mutate (default 'A')
        validate: If True, fail on critical errors
        log_file: Optional path to save mutation log
        has_deletion: If True, convert WT numbering to account for deletions
        deletion_start: First deleted residue (default 67)
        deletion_length: Number of deleted residues (default 2)

    Returns:
        True if successful, False otherwise
    """
    if not PYROSETTA_AVAILABLE:
        logger.error("PyRosetta is not available. Cannot perform threading.")
        return False

    try:
        # Parse signature
        logger.info(f"Parsing signature: {signature}")
        mutations = parse_variant_signature(signature)
        logger.info(f"Parsed {len(mutations)} mutations: {mutations}")

        if has_deletion:
            logger.info(f"Template has deletion: Δ{deletion_start}-{deletion_start+deletion_length-1}")
            logger.info(f"WT positions ≥{deletion_start+deletion_length} will be shifted by -{deletion_length}")

        # Initialize PyRosetta (quietly)
        pyrosetta.init('-mute all')

        # Load template
        logger.info(f"Loading template: {template_pdb}")
        if not os.path.exists(template_pdb):
            raise FileNotFoundError(f"Template PDB not found: {template_pdb}")

        pose = pose_from_pdb(template_pdb)
        logger.info(f"Loaded structure with {pose.total_residue()} residues")

        # Apply mutations (with deletion numbering conversion)
        logger.info(f"Applying mutations to chain {chain}...")
        pose, log_messages = apply_mutations(
            pose, mutations, chain, validate,
            has_deletion=has_deletion,
            deletion_start=deletion_start,
            deletion_length=deletion_length
        )

        # Repack sidechains to resolve clashes from mutations
        logger.info("Repacking sidechains around mutated positions...")
        scorefxn = pyrosetta.get_fa_scorefxn()
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
        packer.task_factory(tf)
        score_before = scorefxn(pose)
        packer.apply(pose)
        score_after = scorefxn(pose)
        logger.info(f"Repack complete: score {score_before:.1f} → {score_after:.1f}")

        # Save output
        logger.info(f"Saving mutated structure: {output_pdb}")
        os.makedirs(os.path.dirname(output_pdb) or '.', exist_ok=True)
        pose.dump_pdb(output_pdb)

        # Save log if requested
        if log_file:
            with open(log_file, 'w') as f:
                f.write(f"Template: {template_pdb}\n")
                f.write(f"Signature: {signature}\n")
                f.write(f"Chain: {chain}\n")
                f.write(f"Has deletion: {has_deletion}\n")
                if has_deletion:
                    f.write(f"Deletion: Δ{deletion_start}-{deletion_start+deletion_length-1}\n")
                f.write(f"Output: {output_pdb}\n\n")
                f.write("Mutations applied:\n")
                for msg in log_messages:
                    f.write(f"  {msg}\n")

        logger.info(f"✓ Successfully created mutant: {output_pdb}")
        return True

    except Exception as e:
        logger.error(f"✗ Threading failed: {e}")
        if validate:
            raise
        return False


def batch_thread_from_csv(
    csv_path: str,
    template_pdb: str,
    output_dir: str,
    chain: str = 'A',
    validate: bool = False
) -> pd.DataFrame:
    """
    Batch process variants from CSV file.

    CSV format:
        variant_name,variant_signature
        PYR1^4F,59K;120A;160G
        PYR1^WIN,83F;115Q;120G

    Args:
        csv_path: Path to CSV file
        template_pdb: Path to template PDB
        output_dir: Directory to save mutated structures
        chain: Chain ID to mutate
        validate: If True, stop on first error

    Returns:
        DataFrame with columns: variant_name, signature, output_pdb, status, log_file
    """
    logger.info(f"Batch processing from CSV: {csv_path}")

    # Load CSV
    df = pd.read_csv(csv_path)
    required_cols = ['variant_name', 'variant_signature']

    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"CSV missing required column: '{col}'")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Process each variant
    results = []

    for idx, row in df.iterrows():
        variant_name = row['variant_name']
        signature = row['variant_signature']

        # Generate output paths
        safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', variant_name)
        output_pdb = os.path.join(output_dir, f"{safe_name}.pdb")
        log_file = os.path.join(output_dir, f"{safe_name}_log.txt")

        logger.info(f"\n[{idx+1}/{len(df)}] Processing {variant_name}...")

        # Thread variant
        success = thread_variant(
            template_pdb=template_pdb,
            signature=signature,
            output_pdb=output_pdb,
            chain=chain,
            validate=validate,
            log_file=log_file
        )

        results.append({
            'variant_name': variant_name,
            'signature': signature,
            'output_pdb': output_pdb if success else None,
            'status': 'SUCCESS' if success else 'FAILED',
            'log_file': log_file if success else None
        })

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    # Save summary
    summary_path = os.path.join(output_dir, 'threading_summary.csv')
    results_df.to_csv(summary_path, index=False)
    logger.info(f"\n✓ Batch processing complete. Summary: {summary_path}")
    logger.info(f"  Success: {(results_df['status'] == 'SUCCESS').sum()}/{len(results_df)}")

    return results_df


def main():
    parser = argparse.ArgumentParser(
        description='Thread PYR1 variant signature onto template PDB',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Single variant
    python thread_variant_to_pdb.py \\
        --template 3QN1_nolig_H2O.pdb \\
        --signature "59K;120A;160G" \\
        --output mutant.pdb

    # Batch processing from CSV
    python thread_variant_to_pdb.py \\
        --template 3QN1_nolig_H2O.pdb \\
        --csv variants.csv \\
        --output-dir mutants/

    # With validation and logging
    python thread_variant_to_pdb.py \\
        --template 3QN1_nolig_H2O.pdb \\
        --signature "K59Q;Y120A;A160G" \\
        --output mutant.pdb \\
        --log mutation_log.txt \\
        --validate
        """
    )

    # Single variant mode
    parser.add_argument('--template', required=True, help='Template PDB file')
    parser.add_argument('--signature', help='Variant signature (e.g., "59K;120A;160G")')
    parser.add_argument('--output', help='Output mutated PDB file')
    parser.add_argument('--log', help='Optional mutation log file')

    # Batch mode
    parser.add_argument('--csv', help='CSV file with variant_name, variant_signature columns')
    parser.add_argument('--output-dir', help='Output directory for batch mode')

    # Options
    parser.add_argument('--chain', default='A', help='Chain ID to mutate (default: A)')
    parser.add_argument('--validate', action='store_true', help='Fail on critical errors')
    parser.add_argument('--test', action='store_true', help='Test signature parsing only (no PyRosetta)')

    # Deletion handling (3QN1 has residues 67-68 deleted)
    parser.add_argument('--no-deletion', action='store_true',
                       help='Template does NOT have deletion (use WT numbering directly)')
    parser.add_argument('--deletion-start', type=int, default=69,
                       help='First deleted residue position (default: 69)')
    parser.add_argument('--deletion-length', type=int, default=2,
                       help='Number of deleted residues (default: 2)')

    args = parser.parse_args()

    # Check PyRosetta availability
    if not args.test and not PYROSETTA_AVAILABLE:
        logger.error("PyRosetta is not installed. Install with: conda install -c conda-forge pyrosetta")
        sys.exit(1)

    # Test mode: just parse signature
    if args.test:
        if not args.signature:
            logger.error("--signature required for test mode")
            sys.exit(1)

        logger.info(f"Testing signature parsing: {args.signature}")
        try:
            mutations = parse_variant_signature(args.signature)
            logger.info(f"✓ Parsed {len(mutations)} mutations:")
            for pos, aa in sorted(mutations.items()):
                logger.info(f"  Position {pos} → {aa}")
            sys.exit(0)
        except Exception as e:
            logger.error(f"✗ Parsing failed: {e}")
            sys.exit(1)

    # Batch mode
    if args.csv:
        if not args.output_dir:
            logger.error("--output-dir required for batch mode")
            sys.exit(1)

        results_df = batch_thread_from_csv(
            csv_path=args.csv,
            template_pdb=args.template,
            output_dir=args.output_dir,
            chain=args.chain,
            validate=args.validate
        )

        # Exit with error if any failed
        if (results_df['status'] == 'FAILED').any():
            sys.exit(1)
        else:
            sys.exit(0)

    # Single variant mode
    else:
        if not args.signature or not args.output:
            logger.error("--signature and --output required for single variant mode")
            parser.print_help()
            sys.exit(1)

        success = thread_variant(
            template_pdb=args.template,
            signature=args.signature,
            output_pdb=args.output,
            chain=args.chain,
            validate=args.validate,
            log_file=args.log,
            has_deletion=not args.no_deletion,  # Default: True (has deletion)
            deletion_start=args.deletion_start,
            deletion_length=args.deletion_length
        )

        sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
