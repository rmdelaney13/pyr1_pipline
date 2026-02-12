#!/usr/bin/env python3

"""
binary_analysis.py

Aggregates AlphaFold3 results by recursively searching for JSONs.
Outputs:
  1. Ligand Chain pLDDT
  2. Ligand Chain iPTM
  3. <ligand_oxygen_atom>–O1 aligned distance

Distance:
  AF3: ligand_chain, residue 1, atom <ligand_oxygen_atom> (default O3)
  REF: chain A, residue 201, atom O1
"""

import os
import json
import csv
import argparse
import numpy as np
from pathlib import Path
from Bio.PDB import MMCIFParser, Superimposer


# ----------------------------
# 1. File Discovery
# ----------------------------
def find_af3_files(base_dir):
    targets = {}
    base = Path(base_dir)

    print(f"Scanning {base_dir} recursively...")

    for summary_path in base.rglob("*_summary_confidences.json"):
        target_name = summary_path.name.replace("_summary_confidences.json", "")
        parent_dir = summary_path.parent

        full_conf_path = parent_dir / f"{target_name}_confidences.json"
        cif_path = parent_dir / f"{target_name}_model.cif"

        targets[target_name] = {
            "summary": str(summary_path),
            "full": str(full_conf_path) if full_conf_path.exists() else None,
            "cif": str(cif_path) if cif_path.exists() else None,
        }

    return targets


# ----------------------------
# 2. JSON Parsing
# ----------------------------
def get_chain_order(full_conf_path):
    if not full_conf_path:
        return None

    try:
        data = json.loads(Path(full_conf_path).read_text())
        for key in ["token_chain_ids", "atom_chain_ids", "chain_ids"]:
            if key in data and isinstance(data[key], list):
                unique = []
                seen = set()
                for c in data[key]:
                    c = str(c)
                    if c not in seen:
                        unique.append(c)
                        seen.add(c)
                return unique
    except Exception:
        pass

    return ["A", "B"]


def extract_json_metrics(file_map, ligand_chain="B"):
    plddt = None
    iptm = None

    # pLDDT
    if file_map.get("full"):
        try:
            data = json.loads(Path(file_map["full"]).read_text())
            plddts = np.array(data.get("atom_plddts", []), dtype=float)
            chains = np.array(data.get("atom_chain_ids", []), dtype=str)

            mask = (chains == ligand_chain)
            if mask.any():
                plddt = float(plddts[mask].mean())
        except Exception as e:
            print(f"Warning reading pLDDT for {file_map['full']}: {e}")

    # iPTM
    if file_map.get("summary"):
        try:
            data = json.loads(Path(file_map["summary"]).read_text())
            chain_iptms = data.get("chain_iptm", [])

            if chain_iptms:
                chain_order = get_chain_order(file_map.get("full"))
                if chain_order and ligand_chain in chain_order:
                    idx = chain_order.index(ligand_chain)
                    if idx < len(chain_iptms):
                        iptm = float(chain_iptms[idx])
                elif len(chain_iptms) == 2 and ligand_chain == "B":
                    iptm = float(chain_iptms[1])
        except Exception as e:
            print(f"Warning reading iPTM for {file_map['summary']}: {e}")

    return plddt, iptm


# ----------------------------
# 3. Geometry
# ----------------------------
def get_full_structure_and_ca_atoms(file_path):
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("model", file_path)
    except Exception:
        return None, None

    ca_atoms = []
    for model in structure:
        all_ca = []
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "CA":
                        all_ca.append(atom)

        ca_atoms = sorted(
            all_ca,
            key=lambda a: (a.get_parent().get_parent().id, a.get_parent().get_id()[1]),
        )
        break

    return structure, ca_atoms


def get_specific_atom(structure, chain_id, res_id, atom_name):
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if residue.get_id()[1] != res_id:
                    continue
                for atom in residue:
                    if atom.get_name() == atom_name:
                        return atom
    return None


def extract_ligox_o1_distance(
    ref_model_path,
    af_cif_path,
    ligand_chain="B",
    ligand_oxygen_atom="O3",
    ligand_res_id=1,
    ref_chain="A",
    ref_res_id=201,
    ref_atom="O1",
):
    if not af_cif_path:
        return None

    ref_struct, ref_ca = get_full_structure_and_ca_atoms(ref_model_path)
    af_struct, af_ca = get_full_structure_and_ca_atoms(af_cif_path)

    if not ref_struct or not af_struct or not ref_ca or not af_ca:
        return None
    if len(ref_ca) != len(af_ca):
        return None

    ref_o1 = get_specific_atom(ref_struct, ref_chain, ref_res_id, ref_atom)
    af_ox = get_specific_atom(af_struct, ligand_chain, ligand_res_id, ligand_oxygen_atom)

    if not ref_o1 or not af_ox:
        return None

    superimposer = Superimposer()
    superimposer.set_atoms(ref_ca, af_ca)
    superimposer.apply([af_ox])

    return float(np.linalg.norm(af_ox.get_coord() - ref_o1.get_coord()))


# ----------------------------
# 4. Main
# ----------------------------
def main(args):
    target_map = find_af3_files(args.inference_dir)
    all_targets = sorted(target_map.keys())
    print(f"Found {len(all_targets)} targets.")

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    dist_col = f"{args.ligand_oxygen_atom}_O1_dist_aligned"

    with open(args.output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["target", "ligand_plddt", "ligand_iptm", dist_col])

        for i, tgt in enumerate(all_targets, 1):
            files = target_map[tgt]

            plddt, iptm = extract_json_metrics(files, ligand_chain=args.ligand_chain)
            dist = extract_ligox_o1_distance(
                args.ref_model,
                files.get("cif"),
                ligand_chain=args.ligand_chain,
                ligand_oxygen_atom=args.ligand_oxygen_atom,
                ligand_res_id=args.ligand_res_id,
                ref_chain=args.ref_chain,
                ref_res_id=args.ref_res_id,
                ref_atom=args.ref_atom,
            )

            writer.writerow([
                tgt,
                f"{plddt:.2f}" if plddt is not None else "NA",
                f"{iptm:.2f}" if iptm is not None else "NA",
                f"{dist:.3f}" if dist is not None else "NA",
            ])

            if i % 100 == 0:
                print(f"Processed {i} targets...")

    print(f"\nDone. Results written to {args.output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate AF3 metrics using recursive search.")
    parser.add_argument("--inference_dir", required=True)
    parser.add_argument("--ref_model", required=True)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--ligand_chain", default="B")

    parser.add_argument("--ligand_oxygen_atom", default="O3",
                        help="Ligand atom to measure from (default O3)")
    parser.add_argument("--ligand_res_id", type=int, default=1)
    parser.add_argument("--ref_chain", default="A")
    parser.add_argument("--ref_res_id", type=int, default=201)
    parser.add_argument("--ref_atom", default="O1")

    args = parser.parse_args()
    main(args)

