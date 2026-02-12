#!/usr/bin/env python3
"""
aggregate_all_metrics.py

Combines Binary and Ternary metrics into a single CSV and calculates 
Ligand RMSD (Binary vs Ternary) on the fly.
"""

import argparse
import csv
import numpy as np
import sys
import os
from pathlib import Path
from Bio.PDB import MMCIFParser, Superimposer

# -----------------------------
# 1. Structure / RMSD Logic
# -----------------------------
def get_structure_and_atoms(cif_path, protein_chain="A", ligand_chain="B"):
    """Parses CIF and extracts CA atoms (for alignment) and Ligand atoms (for RMSD)."""
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("model", cif_path)
    except Exception:
        return None, None, None

    ca_atoms = []
    ligand_atoms = []

    for model in structure:
        for chain in model:
            if chain.id == protein_chain:
                for res in chain:
                    for atom in res:
                        if atom.get_name() == "CA":
                            ca_atoms.append(atom)
            
            if chain.id == ligand_chain:
                for res in chain:
                    for atom in res:
                        # Exclude Hydrogens for cleaner RMSD
                        if atom.element != "H":
                            ligand_atoms.append(atom)
        break # Only process model 0

    # Sort ligand atoms by name to ensure 1-to-1 mapping
    if ligand_atoms:
        ligand_atoms.sort(key=lambda a: a.get_name())
    
    return structure, ca_atoms, ligand_atoms

def calculate_pair_rmsd(ternary_path, binary_path):
    """Aligns Binary(A) to Ternary(A), then measures RMSD of Ligand(B)."""
    # Load Ternary (Reference)
    t_struct, t_ca, t_lig = get_structure_and_atoms(ternary_path)
    if not t_struct or not t_ca or not t_lig:
        return None

    # Load Binary (Mobile)
    b_struct, b_ca, b_lig = get_structure_and_atoms(binary_path)
    if not b_struct or not b_ca or not b_lig:
        return None

    # Quality Checks
    if len(t_ca) != len(b_ca):
        return None # Protein length mismatch
    if len(t_lig) != len(b_lig):
        return None # Ligand atom count mismatch

    # Superimpose Proteins
    si = Superimposer()
    si.set_atoms(t_ca, b_ca)
    si.apply(b_struct.get_atoms())

    # Calculate Ligand RMSD
    diffs = []
    for a_t, a_b in zip(t_lig, b_lig):
        d = a_t.get_coord() - a_b.get_coord()
        diffs.append(np.sum(d * d))

    return np.sqrt(np.average(diffs))

# -----------------------------
# 2. File Finding Helpers
# -----------------------------
def find_cif(base_dir, target_name):
    """
    Finds .cif file for a target. 
    Checks:
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
    
    return None

def normalize_ternary_target(t_name):
    """Extracts 'design_id' from ternary target name (removes _model suffix if present)."""
    t = t_name.strip()
    if t.endswith("_model"):
        return t[:-6]
    return t

# -----------------------------
# 3. Data Processing
# -----------------------------
def load_csv_map(csv_path):
    """Reads CSV into a dict keyed by 'target'."""
    data = {}
    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            tid = row.get("target", "").strip()
            if tid:
                data[tid] = row
    return data

def get_val(row, keys, type_func=float):
    """Helper to safely extract a value from multiple possible column names."""
    for k in keys:
        if k in row and row[k] and row[k].lower() not in ["na", "nan", "none", ""]:
            try:
                return type_func(row[k])
            except:
                pass
    return None

def main():
    parser = argparse.ArgumentParser(description="Aggregate Metrics + Calculate RMSD")
    parser.add_argument("--ternary_csv", required=True)
    parser.add_argument("--binary_csv", required=True)
    parser.add_argument("--ternary_dir", required=True, help="Folder with Ternary CIFs")
    parser.add_argument("--binary_dir", required=True, help="Folder with Binary CIFs")
    parser.add_argument("--output_csv", required=True)
    
    args = parser.parse_args()

    print("Loading CSVs...")
    t_data = load_csv_map(args.ternary_csv)
    b_data = load_csv_map(args.binary_csv) # Keyed by design_id

    # Column mappings (New vs Old names)
    dist_cols = ["min_dist_to_ligand_O_aligned", "O3_O1_dist_aligned", "Dist", "dist"]
    plddt_cols = ["ligand_plddt", "plddt"]
    iptm_cols  = ["ligand_iptm", "chainB_iptm", "iptm"]

    master_rows = []
    
    print(f"Processing {len(t_data)} ternary targets...")

    for i, (t_target, t_row) in enumerate(t_data.items(), 1):
        design_id = normalize_ternary_target(t_target)
        
        # 1. Get Ternary Metrics
        t_plddt = get_val(t_row, plddt_cols)
        t_iptm  = get_val(t_row, iptm_cols)
        t_dist  = get_val(t_row, dist_cols)

        # 2. Get Binary Metrics
        b_row = b_data.get(design_id)
        b_plddt, b_iptm, b_dist = None, None, None
        
        if b_row:
            b_plddt = get_val(b_row, plddt_cols)
            b_iptm  = get_val(b_row, iptm_cols)
            b_dist  = get_val(b_row, dist_cols)

        # 3. Calculate RMSD
        rmsd_val = None
        
        t_file = find_cif(args.ternary_dir, t_target)
        b_file = find_cif(args.binary_dir, design_id)

        if t_file and b_file:
            rmsd_val = calculate_pair_rmsd(t_file, b_file)

        # 4. Build Output Row
        master_rows.append({
            "design_id": design_id,
            "ternary_target": t_target,
            # Ternary
            "ternary_plddt": f"{t_plddt:.2f}" if t_plddt else "NA",
            "ternary_iptm":  f"{t_iptm:.2f}"  if t_iptm  else "NA",
            "ternary_dist":  f"{t_dist:.3f}"  if t_dist  else "NA",
            # Binary
            "binary_plddt":  f"{b_plddt:.2f}" if b_plddt else "NA",
            "binary_iptm":   f"{b_iptm:.2f}"  if b_iptm  else "NA",
            "binary_dist":   f"{b_dist:.3f}"  if b_dist  else "NA",
            # Comparison
            "ligand_rmsd":   f"{rmsd_val:.3f}" if rmsd_val is not None else "NA"
        })

        if i % 50 == 0:
            print(f"  Processed {i}...")

    # Write Result
    headers = [
        "design_id", "ternary_target", 
        "ternary_plddt", "ternary_iptm", "ternary_dist",
        "binary_plddt", "binary_iptm", "binary_dist",
        "ligand_rmsd"
    ]

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)
    
    with open(args.output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(master_rows)

    print(f"\nDone! Master CSV written to: {args.output_csv}")

if __name__ == "__main__":
    main()
