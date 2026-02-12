import argparse
import os
import shutil
from typing import Tuple
import pandas as pd

REQUIRED_COLS = ["filename", "dG_sep", "interface_unsats", "lig_unsat_polars", "charge_satisfied"]

def require_cols(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")

def yes_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(["yes", "true", "1", "1.0"])

def pdb_name_from_csv_filename(csv_filename: str) -> str:
    base = os.path.basename(str(csv_filename)).strip()
    if base.endswith(".sc"): base = base[:-3]
    stem, _ext = os.path.splitext(base)
    return f"{stem}.pdb"

def get_parent_dock(filename: str) -> str:
    return str(filename).split("_design_")[0]

def copy_pdbs(df: pd.DataFrame, src_dir: str, out_dir: str) -> Tuple[int, int]:
    copied, missing = 0, 0
    for _, row in df.iterrows():
        pdb_name = pdb_name_from_csv_filename(row["filename"])
        src = os.path.join(src_dir, pdb_name)
        dst = os.path.join(out_dir, pdb_name)
        if os.path.exists(src):
            shutil.copy2(src, dst)
            copied += 1
        else: missing += 1
    return copied, missing

def filter_designs(df: pd.DataFrame, target_n: int, max_interface_unsat: float, 
                    max_lig_unsat: float, max_per_parent: int, 
                    req_charge: bool) -> pd.DataFrame:
    require_cols(df)
    d = df.copy()

    # Numeric conversion
    for col in ["dG_sep", "interface_unsats", "lig_unsat_polars"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d = d.dropna(subset=["dG_sep", "interface_unsats", "lig_unsat_polars"]).drop_duplicates(subset=["filename"])

    total_count = len(d)
    
    # Track failures for troubleshooting
    fail_interface = d[d["interface_unsats"] > max_interface_unsat]
    fail_ligand = d[d["lig_unsat_polars"] > max_lig_unsat]
    fail_charge = d[~yes_mask(d["charge_satisfied"])] if req_charge else pd.DataFrame()

    mask_combined = (d["interface_unsats"] <= max_interface_unsat) & (d["lig_unsat_polars"] <= max_lig_unsat)
    if req_charge:
        mask_combined &= yes_mask(d["charge_satisfied"])

    valid_candidates = d[mask_combined].copy()
    
    print(f"\n--- Troubleshooting Diagnostics ---")
    print(f"Total designs processed: {total_count}")
    print(f"Failed max_interface_unsat (>{max_interface_unsat}): {len(fail_interface)}")
    print(f"Failed max_lig_unsat (>{max_lig_unsat}): {len(fail_ligand)}")
    if req_charge:
        print(f"Failed charge_satisfied requirement: {len(fail_charge)}")
    
    if total_count > 0 and len(valid_candidates) == 0:
        print(f"\nTIP: Your best (lowest) lig_unsat_polars is {d['lig_unsat_polars'].min()}.")
        print(f"TIP: Your best (lowest) interface_unsats is {d['interface_unsats'].min()}.")

    valid_candidates = valid_candidates.sort_values("dG_sep", ascending=True)
    valid_candidates["parent_dock"] = valid_candidates["filename"].apply(get_parent_dock)
    return valid_candidates.groupby("parent_dock").head(max_per_parent).sort_values("dG_sep", ascending=True).head(target_n)

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    parser.add_argument("relax_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--target_n", type=int, default=500)
    parser.add_argument("--max_interface_unsat", type=float, default=3.0)
    parser.add_argument("--max_lig_unsat", type=float, default=0.0)
    parser.add_argument("--max_per_parent", type=int, default=5)
    parser.add_argument("--require_charge", action="store_true")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    selected = filter_designs(df, args.target_n, args.max_interface_unsat, 
                              args.max_lig_unsat, args.max_per_parent, args.require_charge)

    selected.to_csv(os.path.join(args.output_dir, "filtered_results.csv"), index=False)
    print(f"\n--- Final Results ---")
    print(f"Designs passing filters: {len(selected)}")
    if not selected.empty:
        copied, missing = copy_pdbs(selected, args.relax_dir, args.output_dir)
        print(f"PDBs copied: {copied}")

if __name__ == "__main__":
    main()
