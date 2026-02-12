"""
UNIVERSAL design filtering script.

Auto-detects:
  - All polar contact columns (O1_polar_contact, N1_polar_contact, S1_polar_contact, etc.)
  - All charge satisfaction columns
  - Handles variable numbers of polar atoms per ligand

Usage:
    python relax_filter_universal.py input.csv relax_dir output_dir \\
        --target_n 500 \\
        --max_unsat 5 \\
        --require_all_polar_satisfied \\
        --require_charge_satisfied
"""
import argparse
import os
import re
import shutil
from typing import Tuple, List

import pandas as pd


def detect_polar_contact_columns(df: pd.DataFrame) -> List[str]:
    """
    Auto-detect all polar contact columns in the DataFrame.

    Returns:
        List of column names matching pattern: {atom}_polar_contact
    """
    pattern = re.compile(r'^[A-Z][0-9]*_polar_contact$', re.IGNORECASE)
    polar_cols = [col for col in df.columns if pattern.match(col)]
    return sorted(polar_cols)


def detect_charge_satisfaction_columns(df: pd.DataFrame) -> List[str]:
    """
    Auto-detect charge satisfaction columns.

    Returns:
        List of column names for charge satisfaction
    """
    charge_patterns = [
        'charge_satisfied',
        'carboxylate_satisfied',
        'amine_satisfied',
        'sulfonate_satisfied',
        'phosphate_satisfied'
    ]
    found_cols = [col for col in charge_patterns if col in df.columns]
    return found_cols


def yes_mask(series: pd.Series) -> pd.Series:
    """Matches 'yes', 'true', '1', '1.0' case-insensitive."""
    return series.astype(str).str.strip().str.lower().isin(["yes", "true", "1", "1.0"])


def pdb_name_from_csv_filename(csv_filename: str) -> str:
    """Extract PDB filename from CSV row filename."""
    base = os.path.basename(str(csv_filename)).strip()
    stem, _ext = os.path.splitext(base)
    # Remove common suffixes
    if stem.endswith("_relaxed_score"):
        stem = stem[:-len("_relaxed_score")]
    if stem.endswith("_score"):
        stem = stem[:-len("_score")]
    return f"{stem}.pdb"


def get_parent_dock(filename: str) -> str:
    """
    Extract parent dock name from filename.
    Assumes format: 'arrayXXX_passY_repackedZ_design_N.sc'
    """
    return str(filename).split("_design_")[0]


def filter_designs(
    df: pd.DataFrame,
    target_n: int,
    max_unsat: int,
    max_per_parent: int,
    polar_cols_to_check: List[str],
    charge_cols_to_check: List[str],
    require_any_polar: bool = False
) -> pd.DataFrame:
    """
    Universal filtering with auto-detected columns.

    Args:
        df: Input dataframe
        target_n: Final number of designs to return
        max_unsat: Maximum buried unsatisfied polars
        max_per_parent: Maximum designs per parent dock
        polar_cols_to_check: List of polar contact columns to require
        charge_cols_to_check: List of charge satisfaction columns to require
        require_any_polar: If True, require ANY polar contact (not all)

    Returns:
        Filtered dataframe
    """
    d = df.copy()

    # Clean and convert
    d["dG_sep"] = pd.to_numeric(d["dG_sep"], errors="coerce")
    d["buried_unsatisfied_polars"] = pd.to_numeric(d["buried_unsatisfied_polars"], errors="coerce")
    d = d.dropna(subset=["dG_sep", "buried_unsatisfied_polars"]).copy()
    d = d.drop_duplicates(subset=["filename"], keep="first")

    # Start with all True
    mask_combined = pd.Series([True] * len(d), index=d.index)

    # Apply polar contact filters
    if polar_cols_to_check:
        if require_any_polar:
            # Require at least ONE polar contact
            any_polar_mask = pd.Series([False] * len(d), index=d.index)
            for col in polar_cols_to_check:
                if col in d.columns:
                    any_polar_mask |= yes_mask(d[col])
            mask_combined &= any_polar_mask
        else:
            # Require ALL specified polar contacts
            for col in polar_cols_to_check:
                if col in d.columns:
                    mask_combined &= yes_mask(d[col])
                else:
                    print(f"Warning: Polar column '{col}' not found in data, skipping")

    # Apply charge satisfaction filters
    for col in charge_cols_to_check:
        if col in d.columns:
            mask_combined &= yes_mask(d[col])
        else:
            print(f"Warning: Charge column '{col}' not found in data, skipping")

    # Apply unsats filter
    mask_unsats = d["buried_unsatisfied_polars"] <= max_unsat

    valid_candidates = d[mask_combined & mask_unsats].copy()

    # Sort by score
    valid_candidates = valid_candidates.sort_values("dG_sep", ascending=True)

    # Limit per parent dock
    valid_candidates["parent_dock"] = valid_candidates["filename"].apply(get_parent_dock)
    balanced_candidates = valid_candidates.groupby("parent_dock").head(max_per_parent)

    # Final cutoff
    final_set = balanced_candidates.sort_values("dG_sep", ascending=True).head(target_n)

    return final_set.reset_index(drop=True)


def copy_pdbs(df: pd.DataFrame, src_dir: str, out_dir: str) -> Tuple[int, int]:
    """Copy PDB files to output directory."""
    copied, missing = 0, 0
    for _, row in df.iterrows():
        pdb_name = pdb_name_from_csv_filename(row["filename"])
        src = os.path.join(src_dir, pdb_name)
        dst = os.path.join(out_dir, pdb_name)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            copied += 1
        else:
            missing += 1
    return copied, missing


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Universal design filter with auto-detected polar contacts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Require all polar atoms satisfied + charge satisfied:
  python %(prog)s input.csv relax_dir output_dir --require_all_polar --require_charge

  # Require only specific atoms:
  python %(prog)s input.csv relax_dir output_dir --require_polar O1,O2 --require_charge

  # Require ANY polar contact (not all):
  python %(prog)s input.csv relax_dir output_dir --require_any_polar --require_charge

  # No polar requirements, only charge:
  python %(prog)s input.csv relax_dir output_dir --require_charge
        """
    )

    parser.add_argument("input_csv", help="Input relax CSV")
    parser.add_argument("relax_dir", help="Directory containing relaxed PDBs")
    parser.add_argument("output_dir", help="Output directory for filtered results")

    # Basic filters
    parser.add_argument("--target_n", type=int, default=500,
                       help="Final number of designs to output (default: 500)")
    parser.add_argument("--max_unsat", type=int, default=5,
                       help="Maximum allowed buried unsatisfied polars (default: 5)")
    parser.add_argument("--max_per_parent", type=int, default=5,
                       help="Maximum designs per parent dock (default: 5)")

    # Polar contact filters (mutually exclusive groups)
    polar_group = parser.add_mutually_exclusive_group()
    polar_group.add_argument("--require_all_polar", action="store_true",
                            help="Require ALL detected polar atoms to have contacts")
    polar_group.add_argument("--require_any_polar", action="store_true",
                            help="Require at least ONE polar atom to have contact")
    polar_group.add_argument("--require_polar", type=str,
                            help="Comma-separated list of specific polar atoms to require (e.g., 'O1,O2,N1')")
    polar_group.add_argument("--ignore_polar", action="store_true",
                            help="Ignore all polar contact requirements")

    # Charge satisfaction filters
    charge_group = parser.add_mutually_exclusive_group()
    charge_group.add_argument("--require_charge", action="store_true",
                             help="Require overall charge_satisfied=yes (default behavior)")
    charge_group.add_argument("--require_all_charges", action="store_true",
                             help="Require ALL charge types satisfied (carboxylate, amine, etc.)")
    charge_group.add_argument("--ignore_charge", action="store_true",
                             help="Ignore charge satisfaction")

    # Output options
    parser.add_argument("--output_csv_name", default="filtered.csv",
                       help="Output CSV filename (default: filtered.csv)")
    parser.add_argument("--no_copy_pdbs", action="store_true",
                       help="Do not copy PDB files")

    args = parser.parse_args()

    if args.target_n <= 0:
        raise ValueError("--target_n must be > 0")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Read CSV
    df = pd.read_csv(args.input_csv)

    # Auto-detect available columns
    detected_polar_cols = detect_polar_contact_columns(df)
    detected_charge_cols = detect_charge_satisfaction_columns(df)

    print(f"Auto-detected {len(detected_polar_cols)} polar contact columns:")
    for col in detected_polar_cols:
        print(f"  - {col}")

    print(f"\nAuto-detected {len(detected_charge_cols)} charge satisfaction columns:")
    for col in detected_charge_cols:
        print(f"  - {col}")

    # Determine which polar columns to check
    polar_cols_to_check = []
    require_any = False

    if args.ignore_polar:
        print("\nIgnoring all polar contact requirements")
    elif args.require_all_polar:
        polar_cols_to_check = detected_polar_cols
        print(f"\nRequiring ALL {len(polar_cols_to_check)} polar atoms to have contacts")
    elif args.require_any_polar:
        polar_cols_to_check = detected_polar_cols
        require_any = True
        print(f"\nRequiring at least ONE of {len(polar_cols_to_check)} polar atoms to have contact")
    elif args.require_polar:
        # Parse comma-separated list
        requested = [x.strip() for x in args.require_polar.split(',')]
        for atom in requested:
            col_name = f"{atom}_polar_contact"
            if col_name in detected_polar_cols:
                polar_cols_to_check.append(col_name)
            else:
                print(f"Warning: Requested polar atom '{atom}' not found in data")
        print(f"\nRequiring specific polar atoms: {polar_cols_to_check}")
    else:
        # Default: require all polar atoms
        polar_cols_to_check = detected_polar_cols
        print(f"\nDefault: Requiring ALL {len(polar_cols_to_check)} polar atoms")

    # Determine which charge columns to check
    charge_cols_to_check = []

    if args.ignore_charge:
        print("\nIgnoring all charge requirements")
    elif args.require_all_charges:
        charge_cols_to_check = detected_charge_cols
        print(f"\nRequiring ALL {len(charge_cols_to_check)} charge types satisfied")
    elif args.require_charge or (not args.ignore_charge and not args.require_all_charges):
        # Default: use 'charge_satisfied' if present
        if 'charge_satisfied' in detected_charge_cols:
            charge_cols_to_check = ['charge_satisfied']
            print("\nRequiring overall charge_satisfied")
        else:
            print("\nWarning: charge_satisfied column not found, ignoring charge requirements")

    # Filter designs
    print(f"\nFiltering with:")
    print(f"  Max unsats: {args.max_unsat}")
    print(f"  Max per parent: {args.max_per_parent}")
    print(f"  Target designs: {args.target_n}")

    selected = filter_designs(
        df,
        target_n=args.target_n,
        max_unsat=args.max_unsat,
        max_per_parent=args.max_per_parent,
        polar_cols_to_check=polar_cols_to_check,
        charge_cols_to_check=charge_cols_to_check,
        require_any_polar=require_any
    )

    # Write output
    out_csv = os.path.join(args.output_dir, args.output_csv_name)
    selected.to_csv(out_csv, index=False)

    # Stats
    num_designs = len(selected)
    unique_docks = selected["parent_dock"].nunique() if "parent_dock" in selected else 0

    print(f"\n{'='*60}")
    print(f"✓ Wrote filtered CSV to: {out_csv}")
    print(f"✓ Total designs selected: {num_designs}")
    print(f"✓ Unique parent docks: {unique_docks}")

    if num_designs > 0:
        print(f"\nScore statistics:")
        print(f"  Best dG_sep: {selected['dG_sep'].min():.2f}")
        print(f"  Worst dG_sep: {selected['dG_sep'].max():.2f}")
        print(f"  Mean dG_sep: {selected['dG_sep'].mean():.2f}")

    # Copy PDBs
    if args.no_copy_pdbs:
        print("\nSkipping PDB copying (--no_copy_pdbs set)")
        return

    print("\nCopying PDB files...")
    copied, missing = copy_pdbs(selected, args.relax_dir, args.output_dir)
    print(f"✓ PDBs copied: {copied}")
    if missing:
        print(f"⚠ PDBs missing: {missing}")

    print(f"\n{'='*60}")
    print("Done!")


if __name__ == "__main__":
    main()
