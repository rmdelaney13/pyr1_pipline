#!/usr/bin/env python3

"""
no_hab_analysis_no_sort_v2_chainiptm.py

Aggregates all AlphaFold3 inference results into a single CSV with:
  - ligand chain B chain_iptm (from *_summary_confidences.json)
  - O3–O1 aligned distance (same as before)

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
# JSON discovery
# ----------------------------
def find_summary_confidence_jsons(base_dir):
    """
    Map {target_name: path_to_*_summary_confidences.json}
    """
    out = {}
    for p in Path(base_dir).rglob("*_summary_confidences.json"):
        name = p.stem.rsplit("_summary_confidences", 1)[0]
        out[name] = p
    return out

def find_full_confidence_json(base_dir, target_name):
    """
    Return path to {target}_confidences.json if present, else None.
    """
    p = Path(base_dir) / target_name / f"{target_name}_confidences.json"
    return p if p.exists() else None

def chain_order_from_full_confidences(full_conf_path):
    """
    Best-effort: infer chain order used by AF3 summary arrays.
    AF3 outputs often include token_chain_ids / atom_chain_ids; we take unique order of appearance.
    Returns list like ['A','B',...], or None if not possible.
    """
    try:
        data = json.loads(Path(full_conf_path).read_text())
    except Exception:
        return None

    for key in ("token_chain_ids", "atom_chain_ids", "chain_ids"):
        if key in data and isinstance(data[key], list) and len(data[key]) > 0:
            seq = data[key]
            order = []
            seen = set()
            for c in seq:
                if c not in seen:
                    seen.add(c)
                    order.append(str(c))
            return order if order else None
    return None

def extract_chainB_iptm(summary_path, full_conf_path=None, ligand_chain_id="B"):
    """
    Extract chain_iptm value corresponding to ligand_chain_id (default 'B').
    Uses chain order inferred from full confidences if provided; else falls back to ['A','B'] if 2 chains.
    Returns float or None.
    """
    try:
        sdata = json.loads(Path(summary_path).read_text())
        chain_iptm = sdata.get("chain_iptm", None)
        if chain_iptm is None:
            return None
        # ensure list-like
        if not isinstance(chain_iptm, list) or len(chain_iptm) == 0:
            return None

        # Try to infer chain order robustly
        chain_order = chain_order_from_full_confidences(full_conf_path) if full_conf_path else None

        # Fallback: common case is protein A + ligand B
        if chain_order is None:
            if len(chain_iptm) == 2:
                chain_order = ["A", "B"]
            else:
                # cannot safely map
                return None

        if ligand_chain_id not in chain_order:
            return None

        idx = chain_order.index(ligand_chain_id)
        if idx >= len(chain_iptm):
            return None

        return float(chain_iptm[idx])
    except Exception as e:
        print(f"Warning: Failed to extract chain_iptm from {summary_path}: {e}")
        return None

# ----------------------------
# Geometry parsing (unchanged)
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
        ca_atoms = sorted(
            all_ca_atoms_unsorted,
            key=lambda a: (a.get_parent().get_parent().id, a.get_parent().id[1])
        )
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
    af_o3_atom  = get_specific_atom(af_structure,  'B', 'LIG_B', 1, 'O3')

    if ref_o1_atom is None or af_o3_atom is None:
        return None

    superimposer = Superimposer()
    superimposer.set_atoms(ref_ca_atoms, af_ca_atoms)
    superimposer.apply([af_o3_atom])

    return float(np.linalg.norm(af_o3_atom.get_coord() - ref_o1_atom.get_coord()))

# ----------------------------
# Main
# ----------------------------
def main(args):
    summary_map = find_summary_confidence_jsons(args.inference_dir)
    all_targets = sorted(summary_map.keys())

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    with open(args.output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['target', 'chainB_iptm', 'o3_o1_dist_aligned'])

        for tgt in all_targets:
            # CIF
            cif_path = glob.glob(os.path.join(args.inference_dir, tgt, f"{tgt}_model.cif"))
            if not cif_path:
                print(f"Warning: No CIF file found for target {tgt}.")
                continue
            cif_path = cif_path[0]

            # chain_iptm from summary + (optional) mapping from full confidences
            summary_path = summary_map[tgt]
            full_conf_path = find_full_confidence_json(args.inference_dir, tgt)
            chainB_iptm = extract_chainB_iptm(summary_path, full_conf_path, ligand_chain_id=args.ligand_chain)

            # geometry
            o3_o1_dist = extract_o3_o1_distance(args.ref_model, cif_path)

            writer.writerow([
                tgt,
                f"{chainB_iptm:.3f}" if chainB_iptm is not None else 'NA',
                f"{o3_o1_dist:.3f}" if o3_o1_dist is not None else 'NA'
            ])

    print(f"\nAggregated {len(all_targets)} targets.")
    print(f"Results written to {args.output_csv}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Aggregate AF3 results using chain_iptm (ligand chain) + O3–O1 distance.")
    parser.add_argument('--inference_dir', required=True, help='Root directory containing AF3 outputs (per-target subdirs)')
    parser.add_argument('--ref_model', required=True, help='Reference CIF template path')
    parser.add_argument('--output_csv', required=True, help='Path for combined CSV output')
    parser.add_argument('--ligand_chain', default='B', help='Ligand chain ID for chain_iptm extraction (default: B)')
    args = parser.parse_args()
    main(args)

