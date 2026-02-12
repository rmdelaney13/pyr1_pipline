#!/usr/bin/env python3
"""
sort_and_stage_low_rmsd.py

1. Reads the Master Summary CSV.
2. Filters designs where ligand_rmsd <= Threshold.
3. Sorts them (Lowest RMSD first).
4. Copies/Links the corresponding Ternary and Binary CIFs into an output folder.
"""

import argparse
import csv
import sys
import os
import shutil
from pathlib import Path

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
def safe_float(val):
    if not val: return None
    try:
        if str(val).lower() in ["na", "nan", "none", ""]: return None
        return float(val)
    except ValueError: return None

def find_cif(base_dir, target_name):
    """
    Finds .cif file. Checks:
    1. Flat file: base_dir/target_model.cif
    2. Subfolder: base_dir/target/target_model.cif
    """
    base = Path(base_dir)
    # Check flat
    f = base / f"{target_name}_model.cif"
    if f.exists(): return f
    # Check subfolder
    f = base / target_name / f"{target_name}_model.cif"
    if f.exists(): return f
    # Fallback: check just .cif
    f = base / f"{target_name}.cif"
    if f.exists(): return f
    return None

def stage_file(src, dst, mode='hardlink'):
    if dst.exists(): return
    if mode == 'hardlink':
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)
    elif mode == 'symlink':
        dst.symlink_to(src)
    else:
        shutil.copy2(src, dst)

# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Filter by RMSD and Stage Files")
    # Inputs
    parser.add_argument("--input_csv", required=True, help="Path to master_summary.csv")
    parser.add_argument("--ternary_dir", required=True, help="Path to original Ternary files")
    parser.add_argument("--binary_dir", required=True, help="Path to original Binary files")
    
    # Outputs
    parser.add_argument("--output_dir", required=True, help="Directory to create for results")
    parser.add_argument("--threshold", type=float, default=2.0, help="Max RMSD to keep (default: 2.0)")
    parser.add_argument("--copy_mode", choices=['copy', 'hardlink', 'symlink'], default='hardlink')
    
    args = parser.parse_args()

    # Setup Output Directories
    out_path = Path(args.output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    models_dir = out_path / "models"
    models_dir.mkdir(exist_ok=True)
    
    csv_out_path = out_path / "low_rmsd_sorted.csv"

    print(f"Reading {args.input_csv}...")
    
    kept_rows = []
    
    # 1. Read and Filter
    try:
        with open(args.input_csv, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            
            for row in reader:
                rmsd_val = safe_float(row.get("ligand_rmsd"))
                if rmsd_val is not None and rmsd_val <= args.threshold:
                    kept_rows.append((rmsd_val, row))
                    
    except FileNotFoundError:
        print(f"Error: Could not find file {args.input_csv}")
        sys.exit(1)

    # 2. Sort (Lowest RMSD first)
    kept_rows.sort(key=lambda x: x[0])
    
    print(f"Found {len(kept_rows)} designs with RMSD <= {args.threshold} A.")
    print(f"Staging files to: {models_dir}")

    # 3. Stage Files & Write CSV
    with open(csv_out_path, 'w', newline='') as f:
        # Add a new column for the new filename for easy reference
        if "staged_filename_base" not in fieldnames:
            fieldnames.append("staged_filename_base")
            
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for rank, (rmsd, row) in enumerate(kept_rows, 1):
            t_target = row.get("ternary_target")
            d_id = row.get("design_id")
            
            # Find Source Files
            t_src = find_cif(args.ternary_dir, t_target)
            b_src = find_cif(args.binary_dir, d_id)
            
            # Define Dest Names
            # Format: rank_01_designID_ternary.cif
            base_name = f"rank_{rank:03d}_{d_id}"
            t_dst = models_dir / f"{base_name}_ternary.cif"
            b_dst = models_dir / f"{base_name}_binary.cif"
            
            # Copy/Link Ternary
            if t_src:
                stage_file(t_src, t_dst, args.copy_mode)
            else:
                print(f"Warning: Missing Ternary file for {t_target}")

            # Copy/Link Binary
            if b_src:
                stage_file(b_src, b_dst, args.copy_mode)
            else:
                print(f"Warning: Missing Binary file for {d_id}")

            # Update CSV Row
            row["staged_filename_base"] = base_name
            writer.writerow(row)

    print(f"\nDone.")
    print(f"  - CSV: {csv_out_path}")
    print(f"  - Models: {models_dir}")

if __name__ == "__main__":
    main()
