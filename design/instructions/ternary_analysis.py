#!/usr/bin/env python3

"""
ternary_analysis.py (Updated)
Aggregates AlphaFold3 results by recursively searching for JSONs.
Calculates the minimum distance from the Reference O1 atom to the CLOSEST Oxygen in the aligned ligand.
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
# 3. Geometry (UPDATED)
# ----------------------------
def get_full_structure_and_ca_atoms(file_path, chain_filter=None):
    """
    Parses CIF and returns structure and CA atoms for alignment.
    """
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("model", file_path)
    except Exception:
        return None, None

    ca_atoms = []
    for model in structure:
        for chain in model:
            # FILTER: Skip chains that don't match the requested chain_filter
            if chain_filter and chain.id != chain_filter:
                continue

            for residue in chain:
                for atom in residue:
                    if atom.get_name() == "CA":
                        ca_atoms.append(atom)
        
        # Sort CA atoms to ensure consistent alignment order
        if ca_atoms:
            ca_atoms = sorted(
                ca_atoms,
                key=lambda a: (a.get_parent().get_parent().id, a.get_parent().get_id()[1]),
            )
        break

    return structure, ca_atoms


def get_specific_atom(structure, chain_id, res_id, atom_name):
    """Retrieves a single specific atom (e.g., O1 on the reference)."""
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


def get_ligand_oxygen_atoms(structure, chain_id, res_id):
    """Retrieves ALL oxygen atoms in the specified ligand residue."""
    oxygens = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if residue.get_id()[1] != res_id:
                    continue
                for atom in residue:
                    # Check element (O) or fallback to name starts with O
                    # BioPython usually sets element, but name check is a safe fallback for ligands
                    if atom.element == "O" or atom.get_name().startswith("O"):
                        oxygens.append(atom)
    return oxygens


def extract_min_dist_to_ligand_O_aligned(
    ref_model_path,
    af_cif_path,
    ligand_chain="B",
    ligand_res_id=1,
    ref_chain="A",
    ref_res_id=201,
    ref_atom="O1",
):
    """
    1. Aligns AF model (Chain A) to Ref model (Chain A).
    2. Transforms AF ligand (Chain B).
    3. Calculates distances from Ref Atom (O1) to ALL Oxygen atoms in AF ligand.
    4. Returns the MINIMUM distance.
    """
    if not af_cif_path:
        return None

    # Load Structures & CA atoms for alignment
    ref_struct, ref_ca = get_full_structure_and_ca_atoms(ref_model_path, chain_filter=ref_chain)
    af_struct, af_ca = get_full_structure_and_ca_atoms(af_cif_path, chain_filter=ref_chain)

    if not ref_struct or not af_struct or not ref_ca or not af_ca:
        return None

    # Check length mismatch
    if len(ref_ca) != len(af_ca):
        return None

    # 1. Get the Reference Atom (The stationary target)
    ref_target_atom = get_specific_atom(ref_struct, ref_chain, ref_res_id, ref_atom)
    if not ref_target_atom:
        return None

    # 2. Get ALL atoms in the AF ligand to apply the transformation
    #    (We transform the whole ligand, then check distances)
    af_ligand_residue = None
    af_ligand_atoms = []
    
    for model in af_struct:
        for chain in model:
            if chain.id == ligand_chain:
                for residue in chain:
                    if residue.get_id()[1] == ligand_res_id:
                        af_ligand_residue = residue
                        af_ligand_atoms = list(residue.get_atoms())
                        break
    
    if not af_ligand_atoms:
        return None

    # 3. Align AF Protein -> Ref Protein
    superimposer = Superimposer()
    superimposer.set_atoms(ref_ca, af_ca)
    
    # Apply rotation/translation to the ligand atoms
    superimposer.apply(af_ligand_atoms)

    # 4. Find all Oxygens in the (now aligned) ligand
    aligned_oxygens = [
        a for a in af_ligand_atoms 
        if a.element == "O" or a.get_name().startswith("O")
    ]

    if not aligned_oxygens:
        return None

    # 5. Calculate min distance
    min_dist = float('inf')
    ref_coord = ref_target_atom.get_coord()

    for ox in aligned_oxygens:
        dist = float(np.linalg.norm(ox.get_coord() - ref_coord))
        if dist < min_dist:
            min_dist = dist

    return min_dist


# ----------------------------
# 4. Main
# ----------------------------
def main(args):
    target_map = find_af3_files(args.inference_dir)
    all_targets = sorted(target_map.keys())
    print(f"Found {len(all_targets)} targets.")

    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    # Updated Column Name
    dist_col = "min_dist_to_ligand_O_aligned"

    with open(args.output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["target", "ligand_plddt", "ligand_iptm", dist_col])

        for i, tgt in enumerate(all_targets, 1):
            files = target_map[tgt]

            plddt, iptm = extract_json_metrics(files, ligand_chain=args.ligand_chain)
            
            # Calculate Min Distance to any Oxygen
            dist = extract_min_dist_to_ligand_O_aligned(
                args.ref_model,
                files.get("cif"),
                ligand_chain=args.ligand_chain,
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
    parser.add_argument("--ligand_res_id", type=int, default=1)
    
    parser.add_argument("--ref_chain", default="A")
    parser.add_argument("--ref_res_id", type=int, default=201)
    parser.add_argument("--ref_atom", default="O1", help="Reference atom name (default O1)")

    args = parser.parse_args()
    main(args)
