#!/usr/bin/env python3
"""
relax_1_filter_buckets_with_unsats.py

Combines the logic of "Buckets" (Priority vs Fill) with "Hard Filters" (Max Unsats).

Selection Logic:
  1. Global Filter: discard any design where buried_unsatisfied_polars > max_unsat.
  2. Bucket 1 (Priority): all_polar_and_charge == yes.
  3. Bucket 2 (Fill): (O1_polar_contact == yes OR O2_polar_contact == yes).
     Used only if Bucket 1 < target_n.

Deduplicates by filename (best dG_sep wins).

Usage:
  python relax_1_filter_buckets_with_unsats.py INPUT.csv RELAX_DIR OUTPUT_DIR --target_n 1000 --max_unsat 5
"""

import argparse
import os
import shutil
from typing import Dict, Tuple

import pandas as pd


REQUIRED_COLS = [
    "filename",
    "dG_sep",
    "all_polar_and_charge",
    "O1_polar_contact",
    "O2_polar_contact",
    "buried_unsatisfied_polars",  # Added this requirement
]


def require_cols(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")


def yes_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(["yes", "true", "1", "1.0"])


def pdb_name_from_csv_filename(csv_filename: str) -> str:
    base = os.path.basename(str(csv_filename)).strip()
    stem, _ext = os.path.splitext(base)

    if stem.endswith("_relaxed_score"):
        stem = stem[: -len("_relaxed_score")]
    if stem.endswith("_score"):
        stem = stem[: -len("_score")]

    return f"{stem}.pdb"


def select_designs(
    df: pd.DataFrame,
    target_n: int,
    max_unsat: int,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    
    require_cols(df)

    d = df.copy()
    
    # 1. CLEAN DATA
    d["dG_sep"] = pd.to_numeric(d["dG_sep"], errors="coerce")
    d["buried_unsatisfied_polars"] = pd.to_numeric(d["buried_unsatisfied_polars"], errors="coerce")
    
    # Drop rows with broken scores
    d = d.dropna(subset=["dG_sep", "buried_unsatisfied_polars"]).copy()

    # Sort by Score (best first) and Deduplicate (keep best score per filename)
    d = d.sort_values("dG_sep", ascending=True)
    d = d.drop_duplicates(subset=["filename"], keep="first").reset_index(drop=True)

    # 2. GLOBAL FILTER: Max Unsats
    # We discard anything with too many unsats immediately, for BOTH buckets.
    d = d[d["buried_unsatisfied_polars"] <= max_unsat].copy()

    # 3. DEFINE MASKS
    mask_allpolar = yes_mask(d["all_polar_and_charge"])
    mask_o1o2 = yes_mask(d["O1_polar_contact"]) | yes_mask(d["O2_polar_contact"])

    # 4. BUCKET 1: Perfect Matches (All Polar + Charge)
    bucket1 = d[mask_allpolar].copy().sort_values("dG_sep", ascending=True)

    # If Bucket 1 fills the quota, we stop here.
    if len(bucket1) >= target_n:
        final_df = bucket1.head(target_n).reset_index(drop=True)
        stats = {
            "n_selected": int(len(final_df)),
            "n_bucket1_allpolar": int(len(final_df)),
            "n_bucket2_o1o2_fill": 0,
            "max_unsat_used": max_unsat
        }
        return final_df, stats

    # 5. BUCKET 2: Fill with "Partial" Matches (O1 OR O2)
    # We exclude designs already in Bucket 1 to avoid duplication logic issues,
    # though strictly speaking mask_allpolar is a subset of mask_o1o2 usually.
    kept_filenames = set(bucket1["filename"].tolist())
    remaining = d[~d["filename"].isin(kept_filenames)].copy()

    # Apply the "O1 OR O2" filter to the remaining
    bucket2_candidates = remaining[mask_o1o2.loc[remaining.index]].copy()
    bucket2_candidates = bucket2_candidates.sort_values("dG_sep", ascending=True)

    # Take only what we need
    n_need = target_n - len(bucket1)
    bucket2 = bucket2_candidates.head(n_need).copy()

    # 6. COMBINE
    final_df = pd.concat([bucket1, bucket2], ignore_index=True)
    final_df = final_df.sort_values("dG_sep", ascending=True).reset_index(drop=True)

    stats = {
        "n_selected": int(len(final_df)),
        "n_bucket1_allpolar": int(len(bucket1)),
        "n_bucket2_o1o2_fill": int(len(bucket2)),
        "max_unsat_used": max_unsat
    }
    return final_df, stats


def copy_pdbs(final_df: pd.DataFrame, relax_dir: str, output_dir: str) -> Tuple[int, int]:
    copied = 0
    missing = 0

    for _, row in final_df.iterrows():
        pdb_name = pdb_name_from_csv_filename(row["filename"])
        src = os.path.join(relax_dir, pdb_name)
        dst = os.path.join(output_dir, pdb_name)

        if os.path.exists(src):
            shutil.copy2(src, dst)
            copied += 1
        else:
            missing += 1

    return copied, missing


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter: 1. Max Unsat check. 2. Prioritize All-Polar. 3. Fill with O1-or-O2."
    )
    parser.add_argument("input_csv", help="Input relax CSV")
    parser.add_argument("relax_dir", help="Directory containing relaxed PDBs")
    parser.add_argument("output_dir", help="Output directory for filtered results")
    parser.add_argument("--target_n", type=int, default=500, help="Final number of designs to output")
    parser.add_argument("--max_unsat", type=int, default=5, help="Maximum allowed buried unsatisfied polars (default: 5)")
    parser.add_argument("--output_csv_name", default="filtered.csv")
    parser.add_argument("--no_copy_pdbs", action="store_true", help="Do not copy PDBs")
    args = parser.parse_args()

    if args.target_n <= 0:
        raise ValueError("--target_n must be > 0")

    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    final_df, stats = select_designs(df, target_n=args.target_n, max_unsat=args.max_unsat)

    out_csv = os.path.join(args.output_dir, args.output_csv_name)
    final_df.to_csv(out_csv, index=False)

    # Verify composition in final set
    allpolar_ct = int(yes_mask(final_df["all_polar_and_charge"]).sum())
    o1_ct = int(yes_mask(final_df["O1_polar_contact"]).sum())
    o2_ct = int(yes_mask(final_df["O2_polar_contact"]).sum())
    o1o2_ct = int((yes_mask(final_df["O1_polar_contact"]) | yes_mask(final_df["O2_polar_contact"])).sum())

    print("\n==== FILTER SUMMARY (BUCKETS + UNSAT) ====")
    print(f"Max Unsats Allowed: {stats['max_unsat_used']}")
    print(f"Selected total:     {stats['n_selected']}")
    print("")
    print("Bucket Logic:")
    print(f"  Bucket 1 (All Polar & Charge): {stats['n_bucket1_allpolar']}")
    print(f"  Bucket 2 (Fill w/ O1 or O2):   {stats['n_bucket2_o1o2_fill']}")
    print("")
    print("Properties in Final Set:")
    print(f"  all_polar_and_charge: {allpolar_ct}")
    print(f"  Any O1 or O2 contact: {o1o2_ct}")
    print("========================================\n")

    if args.no_copy_pdbs:
        print(f"Filtered CSV written: {out_csv}")
        print("PDB copying: skipped (--no_copy_pdbs set)")
        return

    copied, missing = copy_pdbs(final_df, args.relax_dir, args.output_dir)
    print(f"Filtered CSV written: {out_csv}")
    print(f"PDBs copied: {copied}")
    if missing:
        print(f"PDBs missing (not found in relax_dir): {missing}")


if __name__ == "__main__":
    main()
