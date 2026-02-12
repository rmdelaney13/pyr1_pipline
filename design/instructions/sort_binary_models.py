#!/usr/bin/env python3
"""
sort_binary_models.py

Sorts AF3 models by:
 1. Ligand pLDDT (must be >= plddt_min)
 2. Ligand iPTM (must be >= iptm_min)
 3. Distance (must be <= max_dist)

Copies/Hardlinks selected models and outputs a filtered CSV.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import os
from pathlib import Path
from typing import Optional, List


def safe_float(v: Optional[str]) -> Optional[float]:
    if v is None:
        return None
    v = str(v).strip()
    if v == "" or v.lower() in {"nan", "none", "null", "na"}:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description="Filter AF3 models by pLDDT, iPTM, and Max Distance.")
    ap.add_argument("--run_dir", required=True, help="Directory containing per-seq subfolders.")
    ap.add_argument("--csv", required=True, help="Input CSV file path.")

    # ID Column
    ap.add_argument("--id_col", default="target", help="Column with seq IDs.")

    # 1. pLDDT Filter (Higher is better)
    ap.add_argument("--plddt_col", default="ligand_plddt", help="Column name for pLDDT.")
    ap.add_argument("--plddt_min", type=float, default=65.0, help="Minimum pLDDT (default: 65.0).")

    # 2. iPTM Filter (Higher is better)
    ap.add_argument("--iptm_col", default="ligand_iptm", help="Column name for iPTM.")
    ap.add_argument("--iptm_min", type=float, default=0.70, help="Minimum iPTM (default: 0.7).")

    # 3. Distance Filter (Lower is better)
    ap.add_argument("--dist_col", default="min_dist_to_ligand_O_aligned", help="Column name for distance.")
    ap.add_argument("--max_dist", type=float, default=3.0, help="Maximum allowed distance (default: 3.0).")

    # Output arguments
    ap.add_argument("--out_dir", required=True, help="Directory to put selected model files.")
    ap.add_argument("--out_csv", help="Optional: Path to write the filtered CSV summary.")
    ap.add_argument("--model_suffix", default="_model.cif", help="Model filename suffix.")

    ap.add_argument("--copy_mode", choices=["copy", "hardlink", "symlink"], default="copy", help="File copy method.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing files in out_dir.")
    ap.add_argument("--dry_run", action="store_true", help="Print actions without writing.")

    args = ap.parse_args()

    # Resolve paths
    run_dir = Path(args.run_dir).resolve()
    csv_path = Path(args.csv).resolve() if Path(args.csv).is_absolute() else (Path.cwd() / args.csv).resolve()
    out_dir = Path(args.out_dir).resolve() if Path(args.out_dir).is_absolute() else (Path.cwd() / args.out_dir).resolve()

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    if not run_dir.exists():
        raise FileNotFoundError(f"run_dir not found: {run_dir}")

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    kept = 0
    seen = 0
    missing = 0
    badrows = 0

    selected_rows = []
    
    # Track stats of selected models
    sel_plddts: List[float] = []
    sel_iptms: List[float] = []
    sel_dists: List[float] = []

    print(f"--- Filtering Criteria ---")
    print(f"pLDDT >= {args.plddt_min}")
    print(f"iPTM  >= {args.iptm_min}")
    print(f"Dist  <= {args.max_dist}")
    print(f"--------------------------")

    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        if fieldnames is None:
            raise ValueError("CSV appears to have no header.")

        # Check columns exist
        required = {args.id_col, args.plddt_col, args.iptm_col, args.dist_col}
        missing_cols = required - set(fieldnames)
        if missing_cols:
            raise ValueError(f"CSV missing required columns: {sorted(missing_cols)}")

        for row in reader:
            seen += 1
            seq_id = (row.get(args.id_col) or "").strip()

            # Parse values
            plddt = safe_float(row.get(args.plddt_col))
            iptm  = safe_float(row.get(args.iptm_col))
            dist  = safe_float(row.get(args.dist_col))

            # Skip bad data
            if not seq_id or plddt is None or iptm is None or dist is None:
                badrows += 1
                continue

            # --- FILTERING LOGIC ---
            if plddt < args.plddt_min:
                continue

            if iptm < args.iptm_min:
                continue

            if dist > args.max_dist:
                continue

            # Check if file exists
            src = run_dir / seq_id / f"{seq_id}{args.model_suffix}"
            if not src.exists():
                if args.dry_run: print(f"Missing file: {src}")
                missing += 1
                continue

            # Keep it
            selected_rows.append(row)
            sel_plddts.append(plddt)
            sel_iptms.append(iptm)
            sel_dists.append(dist)
            kept += 1

            dst = out_dir / src.name

            if dst.exists() and not args.overwrite:
                continue

            if args.dry_run:
                continue

            if dst.exists() and args.overwrite:
                dst.unlink()

            if args.copy_mode == "copy":
                shutil.copy2(src, dst)
            elif args.copy_mode == "hardlink":
                try:
                    os.link(src, dst)
                except FileExistsError:
                    pass
            elif args.copy_mode == "symlink":
                dst.symlink_to(src)

    # Write output CSV if requested
    if args.out_csv and not args.dry_run and selected_rows:
        out_csv_path = Path(args.out_csv).resolve() if Path(args.out_csv).is_absolute() else (Path.cwd() / args.out_csv).resolve()
        out_csv_path.parent.mkdir(parents=True, exist_ok=True)

        with out_csv_path.open("w", newline="") as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(selected_rows)
        print(f"Filtered CSV written to: {out_csv_path}")

    print("\n=== SORT SUMMARY ===")
    print(f"Total Rows Checked: {seen}")
    print(f"Models Selected:    {kept}")
    if kept > 0:
        print(f"  > pLDDT range:    {min(sel_plddts):.1f} - {max(sel_plddts):.1f}")
        print(f"  > iPTM range:     {min(sel_iptms):.2f} - {max(sel_iptms):.2f}")
        print(f"  > Dist range:     {min(sel_dists):.2f} - {max(sel_dists):.2f}")
    print(f"Bad Data Rows:      {badrows}")
    print(f"Missing Files:      {missing}")

if __name__ == "__main__":
    main()
