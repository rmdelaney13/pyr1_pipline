import argparse
import os
import shutil
from typing import Tuple

import pandas as pd

# Updated requirements to match your actual score file columns
REQUIRED_COLS = [
    "filename",
    "dG_sep",
    "buried_unsatisfied_polars",
    "O1_polar_contact",
    "O2_polar_contact",
    "charge_satisfied"
]

def require_cols(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")

def yes_mask(series: pd.Series) -> pd.Series:
    """Matches 'yes', 'true', '1', '1.0' case-insensitive."""
    return series.astype(str).str.strip().str.lower().isin(["yes", "true", "1", "1.0"])

def pdb_name_from_csv_filename(csv_filename: str) -> str:
    base = os.path.basename(str(csv_filename)).strip()
    stem, _ext = os.path.splitext(base)
    # Remove common suffixes if present
    if stem.endswith("_relaxed_score"):
        stem = stem[:-len("_relaxed_score")]
    if stem.endswith("_score"):
        stem = stem[:-len("_score")]
    return f"{stem}.pdb"

def get_parent_dock(filename: str) -> str:
    """
    Extracts the parent dock name from the filename.
    Assumes format: 'arrayXXX_passY_repackedZ_design_N.sc'
    Splits by '_design_' and takes the prefix.
    """
    return str(filename).split("_design_")[0]

def filter_designs(df: pd.DataFrame, target_n: int, max_unsat: int, max_per_parent: int, 
                   check_o1: bool, check_o2: bool, check_charge: bool) -> pd.DataFrame:
    require_cols(df)
    d = df.copy()

    # 1. Clean and Convert Data Types
    d["dG_sep"] = pd.to_numeric(d["dG_sep"], errors="coerce")
    d["buried_unsatisfied_polars"] = pd.to_numeric(d["buried_unsatisfied_polars"], errors="coerce")
    d = d.dropna(subset=["dG_sep", "buried_unsatisfied_polars"]).copy()
    d = d.drop_duplicates(subset=["filename"], keep="first")

    # 2. Apply Dynamic Filters
    # Start with a mask of all True
    mask_combined = pd.Series([True] * len(d), index=d.index)

    if check_o1:
        mask_combined &= yes_mask(d["O1_polar_contact"])
    
    if check_o2:
        mask_combined &= yes_mask(d["O2_polar_contact"])
        
    if check_charge:
        mask_combined &= yes_mask(d["charge_satisfied"])

    mask_unsats = d["buried_unsatisfied_polars"] <= max_unsat

    valid_candidates = d[mask_combined & mask_unsats].copy()

    # 3. Sort by Score (Best dG_sep at top)
    valid_candidates = valid_candidates.sort_values("dG_sep", ascending=True)

    # 4. Identify Parent Dock and Limit per Parent
    valid_candidates["parent_dock"] = valid_candidates["filename"].apply(get_parent_dock)

    # Group by parent and take the top N best scoring for that parent
    balanced_candidates = valid_candidates.groupby("parent_dock").head(max_per_parent)

    # 5. Final Sort (Global score) and Cutoff
    final_set = balanced_candidates.sort_values("dG_sep", ascending=True).head(target_n)

    return final_set.reset_index(drop=True)

def copy_pdbs(df: pd.DataFrame, src_dir: str, out_dir: str) -> Tuple[int, int]:
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
    parser = argparse.ArgumentParser(description="Filter designs with configurable polar contact requirements.")
    parser.add_argument("input_csv", help="Input relax CSV")
    parser.add_argument("relax_dir", help="Directory containing relaxed PDBs")
    parser.add_argument("output_dir", help="Output directory for filtered results")
    parser.add_argument("--target_n", type=int, default=500, help="Final number of designs to output")
    parser.add_argument("--max_unsat", type=int, default=5, help="Maximum allowed buried unsatisfied polars")
    parser.add_argument("--max_per_parent", type=int, default=5, help="Maximum designs allowed per parent dock structure")
    parser.add_argument("--output_csv_name", default="filtered.csv")
    parser.add_argument("--no_copy_pdbs", action="store_true", help="Do not copy PDBs")
    
    # New Filter Flags
    parser.add_argument("--ignore_o1", action="store_true", help="If set, O1 polar contact is NOT required")
    parser.add_argument("--ignore_o2", action="store_true", help="If set, O2 polar contact is NOT required")
    parser.add_argument("--ignore_charge", action="store_true", help="If set, charge satisfaction is NOT required")

    args = parser.parse_args()

    if args.target_n <= 0:
        raise ValueError("--target_n must be > 0")

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    
    # Determine requirements based on flags
    req_o1 = not args.ignore_o1
    req_o2 = not args.ignore_o2
    req_charge = not args.ignore_charge

    print(f"Filtering: max {args.max_per_parent} per parent, max {args.max_unsat} unsats...")
    print(f"Criteria Active:")
    print(f"  - O1 Required: {req_o1}")
    print(f"  - O2 Required: {req_o2}")
    print(f"  - Charge Required: {req_charge}")

    selected = filter_designs(df,
                              target_n=args.target_n,
                              max_unsat=args.max_unsat,
                              max_per_parent=args.max_per_parent,
                              check_o1=req_o1,
                              check_o2=req_o2,
                              check_charge=req_charge)

    out_csv = os.path.join(args.output_dir, args.output_csv_name)
    selected.to_csv(out_csv, index=False)

    # --- Stats ---
    num_designs = len(selected)
    unique_docks = selected["parent_dock"].nunique() if "parent_dock" in selected else 0

    print(f"\nWrote filtered CSV to: {out_csv}")
    print(f"Total designs selected: {num_designs}")
    print(f"Unique parent docks represented: {unique_docks}")

    if args.no_copy_pdbs:
        print("Skipping PDB copying (--no_copy_pdbs set)")
        return

    copied, missing = copy_pdbs(selected, args.relax_dir, args.output_dir)
    print(f"PDBs copied: {copied}")
    if missing:
        print(f"PDBs missing (not found in relax_dir): {missing}")

if __name__ == "__main__":
    main()