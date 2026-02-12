#!/usr/bin/env python3

"""
no_hab_analysis_aggregate_all.py

Aggregates all AlphaFold3 inference results into a single CSV with:
  - ligand chain B average pLDDT
  - O3–O1 aligned distance
No filtering, sorting, or copying is performed.
"""

import os
import json
import glob
import csv
import argparse
import numpy as np
from pathlib import Path
from Bio.PDB import MMCIFParser, Superimposer

# ----------------------------
# AF3‐confidence & pLDDT parsing
# ----------------------------
def find_confidence_jsons(base_dir):
    design_map = {}
    for p in Path(base_dir).rglob("*_confidences.json"):
        if "_summary" in p.name: 
            continue
        name = p.stem.rsplit("_confidences", 1)[0]
        design_map[name] = p
    for p in Path(base_dir).rglob("*_data.json"):
        if "_summary" in p.name: 
            continue
        name = p.stem.rsplit("_data", 1)[0]
        if name not in design_map:
            design_map[name] = p
    return design_map

def extract_plddt(inference_dir):
    jmap = find_confidence_jsons(inference_dir)
    out = {}
    for name, js_path in jmap.items():
        try:
            data = json.loads(js_path.read_text())
            plddts = np.array(data.get("atom_plddts", []), dtype=float)
            chains = np.array(data.get("atom_chain_ids", []), dtype=str)
            if plddts.size == 0 or chains.size == 0:
                out[name] = None
                continue
            mask = (chains == "B")
            out[name] = float(plddts[mask].mean()) if mask.any() else None
        except Exception as e:
            print(f"Warning: Failed to extract pLDDT from {js_path}: {e}")
            out[name] = None
    return out

# ----------------------------
# Geometry parsing
# ----------------------------
def get_full_structure_and_ca_atoms(file_path):
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure('model', file_path)
    except Exception as e:
        print(f"Error parsing CIF file: {e}")
        return None, None

    ca_atoms = []
    for model in structure:
        all_ca_atoms_unsorted = []
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == 'CA':
                        all_ca_atoms_unsorted.append(atom)
        ca_atoms = sorted(all_ca_atoms_unsorted, key=lambda a: (a.get_parent().get_parent().id, a.get_parent().id[1]))
        break
    return structure, ca_atoms

def get_specific_atom(structure, chain_id, res_name, res_id, atom_name):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.get_resname().strip() == res_name and residue.get_id()[1] == res_id:
                        for atom in residue:
                            if atom.get_name() == atom_name:
                                return atom
    return None

def extract_o3_o1_distance(ref_model_path, af_cif_path):
    ref_structure, ref_ca_atoms = get_full_structure_and_ca_atoms(ref_model_path)
    af_structure, af_ca_atoms = get_full_structure_and_ca_atoms(af_cif_path)

    if ref_structure is None or af_structure is None or not ref_ca_atoms or not af_ca_atoms:
        return None

    if len(ref_ca_atoms) != len(af_ca_atoms):
        return None

    ref_o1_atom = get_specific_atom(ref_structure, 'A', 'WI5', 201, 'O1')
    af_o3_atom = get_specific_atom(af_structure, 'B', 'LIG_B', 1, 'O3')

    if ref_o1_atom is None or af_o3_atom is None:
        return None

    superimposer = Superimposer()
    superimposer.set_atoms(ref_ca_atoms, af_ca_atoms)
    superimposer.apply([af_o3_atom])

    return np.linalg.norm(af_o3_atom.get_coord() - ref_o1_atom.get_coord())

# ----------------------------
# Main
# ----------------------------
def main(args):
    plddt_map = extract_plddt(args.inference_dir)
    all_targets = list(plddt_map.keys())

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    with open(args.output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['target', 'chainB_avg_plddt', 'o3_o1_dist_aligned'])

        for tgt in all_targets:
            cif_path = glob.glob(os.path.join(args.inference_dir, tgt, f"{tgt}_model.cif"))
            if not cif_path:
                print(f"Warning: No CIF file found for target {tgt}.")
                continue

            cif_path = cif_path[0]
            plddtB = plddt_map.get(tgt, None)
            o3_o1_dist = extract_o3_o1_distance(args.ref_model, cif_path)

            writer.writerow([
                tgt,
                f"{plddtB:.2f}" if plddtB is not None else 'NA',
                f"{o3_o1_dist:.3f}" if o3_o1_dist is not None else 'NA'
            ])

    print(f"\nAggregated {len(all_targets)} targets.")
    print(f"Results written to {args.output_csv}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Aggregate AF3 results without filtering or sorting.")
    parser.add_argument('--inference_dir', required=True, help='Root directory containing AF3 JSONs and *_model.cif files')
    parser.add_argument('--ref_model', required=True, help='Reference CIF template path')
    parser.add_argument('--output_csv', required=True, help='Path for combined CSV output')
    args = parser.parse_args()
    main(args)

