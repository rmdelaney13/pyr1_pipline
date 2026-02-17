# -*- coding: utf-8 -*-
"""
UNIVERSAL Relax and score script for protein-ligand-water systems.

Auto-detects:
  - All polar atoms (N, O, S) in the ligand
  - Charged groups (carboxylates, amines, sulfonates, phosphates)
  - Generates dynamic scoring columns for each polar atom found

Usage:
    python relax_general_universal.py input.pdb output.pdb ligand.params [--config config.txt]
"""
import sys
import os
import math
import argparse
from collections import defaultdict
import pyrosetta
from pyrosetta import init, pose_from_file, create_score_function
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, IncludeCurrent
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, NeighborhoodResidueSelector, AndResidueSelector,
    NotResidueSelector, ResidueIndexSelector, OrResidueSelector
)
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, angle_constraint
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint, AngleConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
from pyrosetta.rosetta.protocols import jd2


# =============================================================================
# UNIVERSAL LIGAND ANALYSIS FUNCTIONS
# =============================================================================

def get_ligand_polar_atoms(pose, ligand_res_idx, include_sulfur=True):
    """
    Auto-detect all polar heavy atoms in the ligand.

    Returns:
        dict: {atom_name: element} for all polar atoms
    """
    polar_elements = {"N", "O"}
    if include_sulfur:
        polar_elements.add("S")

    lig_res = pose.residue(ligand_res_idx)
    polar_atoms = {}

    for i in range(1, lig_res.natoms() + 1):
        atom_name = lig_res.atom_name(i).strip()
        element = lig_res.atom_type(i).element()

        # Skip hydrogens
        if element == "H":
            continue

        if element in polar_elements:
            polar_atoms[atom_name] = element

    return polar_atoms


def detect_charged_groups(pose, ligand_res_idx):
    """
    Auto-detect charged groups in the ligand by analyzing bonding patterns.

    Returns:
        dict: {
            'carboxylates': [(C_atom, [O_atoms])],
            'amines': [N_atoms],
            'sulfonates': [(S_atom, [O_atoms])],
            'phosphates': [(P_atom, [O_atoms])]
        }
    """
    lig_res = pose.residue(ligand_res_idx)
    result = {
        'carboxylates': [],
        'amines': [],
        'sulfonates': [],
        'phosphates': []
    }

    # Build atom connectivity map from residue type bond graph
    # (Residue.is_bonded checks inter-residue bonds, not intra-residue atoms)
    atom_neighbors = defaultdict(list)
    res_type = lig_res.type()
    try:
        for i in range(1, lig_res.natoms() + 1):
            for j in res_type.bonded_neighbor(i):
                atom_neighbors[i].append(j)
    except Exception:
        # Fallback: distance-based bond detection
        for i in range(1, lig_res.natoms() + 1):
            for j in range(i + 1, lig_res.natoms() + 1):
                if lig_res.xyz(i).distance(lig_res.xyz(j)) < 1.9:
                    atom_neighbors[i].append(j)
                    atom_neighbors[j].append(i)

    # Detect carboxylates: C bonded to 2+ oxygens
    for i in range(1, lig_res.natoms() + 1):
        if lig_res.atom_type(i).element() == "C":
            bonded_oxygens = []
            for j in atom_neighbors[i]:
                if lig_res.atom_type(j).element() == "O":
                    bonded_oxygens.append(lig_res.atom_name(j).strip())

            if len(bonded_oxygens) >= 2:
                c_name = lig_res.atom_name(i).strip()
                result['carboxylates'].append((c_name, bonded_oxygens))

    # Detect amines: N atoms (can be charged or uncharged)
    for i in range(1, lig_res.natoms() + 1):
        if lig_res.atom_type(i).element() == "N":
            # Check if bonded to hydrogens (primary/secondary/tertiary amine)
            bonded_h = sum(1 for j in atom_neighbors[i]
                          if lig_res.atom_type(j).element() == "H")
            if bonded_h >= 1:  # Has at least one H
                result['amines'].append(lig_res.atom_name(i).strip())

    # Detect sulfonates: S bonded to 3+ oxygens
    for i in range(1, lig_res.natoms() + 1):
        if lig_res.atom_type(i).element() == "S":
            bonded_oxygens = []
            for j in atom_neighbors[i]:
                if lig_res.atom_type(j).element() == "O":
                    bonded_oxygens.append(lig_res.atom_name(j).strip())

            if len(bonded_oxygens) >= 3:
                s_name = lig_res.atom_name(i).strip()
                result['sulfonates'].append((s_name, bonded_oxygens))

    # Detect phosphates: P bonded to 3+ oxygens
    for i in range(1, lig_res.natoms() + 1):
        if lig_res.atom_type(i).element() == "P":
            bonded_oxygens = []
            for j in atom_neighbors[i]:
                if lig_res.atom_type(j).element() == "O":
                    bonded_oxygens.append(lig_res.atom_name(j).strip())

            if len(bonded_oxygens) >= 3:
                p_name = lig_res.atom_name(i).strip()
                result['phosphates'].append((p_name, bonded_oxygens))

    return result


def check_polar_contact(pose, ligand_res_idx, lig_atom_name, max_dist=3.5):
    """
    Check if a specific ligand atom has polar contacts.
    Uses element-based detection (more reliable than is_acceptor).

    Returns:
        bool: True if polar contact found within max_dist
    """
    polar_elements = {"N", "O", "S"}

    try:
        lig_res = pose.residue(ligand_res_idx)
        if not lig_res.has(lig_atom_name):
            return False

        lig_xyz = lig_res.xyz(lig_atom_name)

        for i in range(1, pose.total_residue() + 1):
            if i == ligand_res_idx:
                continue

            res = pose.residue(i)
            for j in range(1, res.natoms() + 1):
                element = res.atom_type(j).element()
                if element in polar_elements:
                    if lig_xyz.distance(res.xyz(j)) <= max_dist:
                        return True

        return False
    except Exception as e:
        print(f"Warning: Could not check polar contact for {lig_atom_name}: {e}")
        return False


def check_charge_satisfaction(pose, ligand_res_idx, charged_groups,
                              hbond_dist=3.0, hbond_angle=120,
                              salt_bridge_dist=3.5):
    """
    Universal charge satisfaction checker.

    Checks:
      - Carboxylates: Need 2 H-bonds OR 1 salt bridge to Arg/Lys
      - Amines: Need 1 H-bond OR 1 salt bridge to Asp/Glu
      - Sulfonates: Need 2 H-bonds OR 1 salt bridge to Arg/Lys
      - Phosphates: Need 2 H-bonds OR 1 salt bridge to Arg/Lys

    Returns:
        dict: {
            'carboxylate_satisfied': bool,
            'amine_satisfied': bool,
            'sulfonate_satisfied': bool,
            'phosphate_satisfied': bool,
            'overall_satisfied': bool
        }
    """
    pdb_info = pose.pdb_info()
    lig_res = pose.residue(ligand_res_idx)

    # Get H-bond set
    hb_set = HBondSet()
    fill_hbond_set(pose, False, hb_set)

    result = {
        'carboxylate_satisfied': True,  # Default true if no groups
        'amine_satisfied': True,
        'sulfonate_satisfied': True,
        'phosphate_satisfied': True
    }

    # Check carboxylates (anionic)
    for c_atom, o_atoms in charged_groups.get('carboxylates', []):
        satisfied = False

        # Count H-bonds to oxygens
        hbond_count = 0
        for o_name in o_atoms:
            if not lig_res.has(o_name):
                continue
            o_idx = lig_res.atom_index(o_name)

            for i in range(1, hb_set.nhbonds() + 1):
                h = hb_set.hbond(i)
                if h.acc_res() == ligand_res_idx and h.acc_atm() == o_idx:
                    d = h.get_HAdist(pose)
                    a = h.get_AHDangle(pose)
                    if d < hbond_dist and a > hbond_angle:
                        hbond_count += 1

        # Check salt bridge to Arg/Lys
        salt_bridge = False
        pos_res = {"ARG": ["NH1", "NH2"], "LYS": ["NZ"]}
        for o_name in o_atoms:
            if not lig_res.has(o_name):
                continue
            o_xyz = lig_res.xyz(o_name)

            for i in range(1, pose.total_residue() + 1):
                if i == ligand_res_idx:
                    continue
                r = pose.residue(i)
                if r.name3() in pos_res:
                    for at in pos_res[r.name3()]:
                        if r.has(at):
                            if o_xyz.distance(r.xyz(r.atom_index(at))) <= salt_bridge_dist:
                                salt_bridge = True
                                break

        satisfied = (hbond_count >= 2) or salt_bridge
        result['carboxylate_satisfied'] = result['carboxylate_satisfied'] and satisfied

    # Check amines (cationic)
    for n_atom in charged_groups.get('amines', []):
        satisfied = False

        if not lig_res.has(n_atom):
            continue
        n_xyz = lig_res.xyz(n_atom)

        # Count H-bonds (amine as donor)
        hbond_count = 0
        n_idx = lig_res.atom_index(n_atom)
        for i in range(1, hb_set.nhbonds() + 1):
            h = hb_set.hbond(i)
            if h.don_res() == ligand_res_idx:
                # Check if donor heavy atom is our nitrogen
                don_heavy = h.don_hatm() - 1  # Approximation
                if don_heavy == n_idx or abs(don_heavy - n_idx) <= 2:
                    d = h.get_HAdist(pose)
                    a = h.get_AHDangle(pose)
                    if d < hbond_dist and a > hbond_angle:
                        hbond_count += 1

        # Check salt bridge to Asp/Glu
        salt_bridge = False
        neg_res = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
        for i in range(1, pose.total_residue() + 1):
            if i == ligand_res_idx:
                continue
            r = pose.residue(i)
            if r.name3() in neg_res:
                for at in neg_res[r.name3()]:
                    if r.has(at):
                        if n_xyz.distance(r.xyz(r.atom_index(at))) <= salt_bridge_dist:
                            salt_bridge = True
                            break

        satisfied = (hbond_count >= 1) or salt_bridge
        result['amine_satisfied'] = result['amine_satisfied'] and satisfied

    # Check sulfonates and phosphates (similar to carboxylates)
    for group_type in ['sulfonates', 'phosphates']:
        for central_atom, o_atoms in charged_groups.get(group_type, []):
            satisfied = False

            hbond_count = 0
            for o_name in o_atoms:
                if not lig_res.has(o_name):
                    continue
                o_idx = lig_res.atom_index(o_name)

                for i in range(1, hb_set.nhbonds() + 1):
                    h = hb_set.hbond(i)
                    if h.acc_res() == ligand_res_idx and h.acc_atm() == o_idx:
                        d = h.get_HAdist(pose)
                        a = h.get_AHDangle(pose)
                        if d < hbond_dist and a > hbond_angle:
                            hbond_count += 1

            salt_bridge = False
            pos_res = {"ARG": ["NH1", "NH2"], "LYS": ["NZ"]}
            for o_name in o_atoms:
                if not lig_res.has(o_name):
                    continue
                o_xyz = lig_res.xyz(o_name)

                for i in range(1, pose.total_residue() + 1):
                    if i == ligand_res_idx:
                        continue
                    r = pose.residue(i)
                    if r.name3() in pos_res:
                        for at in pos_res[r.name3()]:
                            if r.has(at):
                                if o_xyz.distance(r.xyz(r.atom_index(at))) <= salt_bridge_dist:
                                    salt_bridge = True
                                    break

            satisfied = (hbond_count >= 2) or salt_bridge
            key = group_type[:-1] + '_satisfied'  # Remove 's' from plural
            result[key] = result[key] and satisfied

    # Overall satisfaction
    result['overall_satisfied'] = all([
        result['carboxylate_satisfied'],
        result['amine_satisfied'],
        result['sulfonate_satisfied'],
        result['phosphate_satisfied']
    ])

    return result


# =============================================================================
# SYSTEM-SPECIFIC FUNCTIONS (Keep for backwards compatibility)
# =============================================================================

def constrain_water_network(pose, water_res, ligand_res, ligand_o_atom="O3"):
    """
    Optional: System-specific water network constraints.
    Set ligand_o_atom to match your ligand's water-interacting oxygen.
    """
    pdb_info = pose.pdb_info()

    def add_distance_constraint(r1, a1, r2, a2, dist=1.9, sd=0.2):
        if not pose.residue(r1).has(a1) or not pose.residue(r2).has(a2):
            return
        aid1 = AtomID(pose.residue(r1).atom_index(a1), r1)
        aid2 = AtomID(pose.residue(r2).atom_index(a2), r2)
        pose.add_constraint(AtomPairConstraint(aid1, aid2, HarmonicFunc(dist, sd)))

    def find_residue(chain, number):
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.chain(i) == chain and pdb_info.number(i) == number:
                return i
        return None

    # PYR1-specific residues (modify for your system)
    a86_res = find_residue("A", 86)
    a114_res = find_residue("A", 114)
    c360_res = find_residue("C", 360)

    if a86_res and a114_res and c360_res:
        add_distance_constraint(water_res, "H1", a86_res, "O")
        add_distance_constraint(water_res, "O", a114_res, "H", dist=1.9, sd=0.1)
        add_distance_constraint(water_res, "O", c360_res, "HE1")
        add_distance_constraint(water_res, "H2", ligand_res, ligand_o_atom)


def apply_interface_relax(pose, scorefxn, water_chain="D"):
    """
    Apply FastRelax to interface region.
    Water is kept FIXED (no bb/chi/jump movement) to preserve pocket geometry.
    """
    ligand_sel = ChainSelector("B")

    # Build water selector: try chain first, then residue name
    water_names = {'TP3', 'HOH', 'WAT', 'TIP', 'TIP3'}
    water_residues = set()
    pdb_info = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        if pdb_info.chain(i) == water_chain:
            water_residues.add(i)
        elif pose.residue(i).name3().strip() in water_names:
            water_residues.add(i)

    near10 = NeighborhoodResidueSelector(ligand_sel, 10.0)
    nearA10 = AndResidueSelector(near10, ChainSelector("A"))
    nearC10 = AndResidueSelector(near10, ChainSelector("C"))

    mm = MoveMap()
    total = pose.total_residue()
    for i in range(1, total + 1):
        # Water is FIXED - do not allow movement
        if i in water_residues:
            mm.set_bb(i, False)
            mm.set_chi(i, False)
        elif nearA10.apply(pose)[i] or nearC10.apply(pose)[i]:
            mm.set_bb(i, True)
            mm.set_chi(i, True)
        else:
            mm.set_bb(i, False)
            mm.set_chi(i, False)

    ft = pose.fold_tree()
    for j in range(1, ft.num_jump() + 1):
        down = ft.downstream_jump_residue(j)
        # Only allow ligand jump to move, NOT water
        if ligand_sel.apply(pose)[down] and down not in water_residues:
            mm.set_jump(j, True)

    tf = TaskFactory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(IncludeCurrent())
    packer = PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(pose)

    relax = FastRelax(scorefxn)
    relax.set_movemap(mm)
    relax.apply(pose)


def compute_total_hb_to_lig(pose, ligand_chain="B"):
    """Count total H-bonds to ligand."""
    hb_set = HBondSet()
    fill_hbond_set(pose, False, hb_set)
    pdb_info = pose.pdb_info()
    return sum(1 for i in range(1, hb_set.nhbonds() + 1)
               if ligand_chain in {pdb_info.chain(hb_set.hbond(i).don_res()),
                                   pdb_info.chain(hb_set.hbond(i).acc_res())})


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Universal relax and score for protein-ligand complexes"
    )
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("output_pdb", help="Output PDB file")
    parser.add_argument("ligand_params", help="Ligand params file")
    parser.add_argument("--xml_path",
                       default="/projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/rosetta/interface_scoring.xml",
                       help="Path to RosettaScripts XML for interface scoring")
    parser.add_argument("--ligand_chain", default="B", help="Ligand chain ID")
    parser.add_argument("--water_chain", default="D", help="Water chain ID (if present)")
    parser.add_argument("--skip_water_constraints", action="store_true",
                       help="Skip system-specific water constraints")

    args = parser.parse_args()

    # Initialize PyRosetta
    init(f"-extra_res_fa {args.ligand_params} -ex1 -ex2aro -use_input_sc "
         f"-relax:fast -relax:default_repeats 5 "
         f"-corrections::beta_nov16 true "
         f"-score:weights beta_nov16 "
         f"-mute all")

    pose = pose_from_file(args.input_pdb)
    pdb_info = pose.pdb_info()

    # Find ligand residue
    lig_res_idx = None
    for i in range(1, pose.total_residue() + 1):
        if pdb_info.chain(i) == args.ligand_chain:
            lig_res_idx = i
            break

    if lig_res_idx is None:
        print(f"ERROR: Could not find ligand on chain {args.ligand_chain}")
        sys.exit(1)

    print(f"Found ligand at residue {lig_res_idx}")

    # Auto-detect ligand properties
    polar_atoms = get_ligand_polar_atoms(pose, lig_res_idx)
    charged_groups = detect_charged_groups(pose, lig_res_idx)

    print(f"\nDetected {len(polar_atoms)} polar atoms: {list(polar_atoms.keys())}")
    print(f"Detected charged groups:")
    print(f"  Carboxylates: {len(charged_groups['carboxylates'])}")
    print(f"  Amines: {len(charged_groups['amines'])}")
    print(f"  Sulfonates: {len(charged_groups['sulfonates'])}")
    print(f"  Phosphates: {len(charged_groups['phosphates'])}")

    # Setup scoring
    sf = create_score_function("beta_nov16")
    sf.set_weight(atom_pair_constraint, 1.0)
    sf.set_weight(angle_constraint, 1.0)

    # Optional: Water constraints (system-specific)
    if not args.skip_water_constraints:
        # Find water: first by chain, then by residue name (TP3, HOH, WAT)
        wat_res_idx = None
        water_names = {'TP3', 'HOH', 'WAT', 'TIP', 'TIP3'}

        # Method 1: search by chain letter
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.chain(i) == args.water_chain:
                wat_res_idx = i
                break

        # Method 2: search by residue name (TP3/HOH/WAT)
        if wat_res_idx is None:
            for i in range(1, pose.total_residue() + 1):
                if pose.residue(i).name3().strip() in water_names:
                    wat_res_idx = i
                    break

        if wat_res_idx is not None:
            water_res = pose.residue(wat_res_idx)
            print(f"Found water at residue {wat_res_idx} "
                  f"(chain={pdb_info.chain(wat_res_idx)}, name={water_res.name3().strip()})")
            # Try to find which oxygen interacts with water (heuristic)
            water_lig_atom = None
            for atom_name in polar_atoms.keys():
                if polar_atoms[atom_name] == "O":
                    water_lig_atom = atom_name
                    break
            if water_lig_atom:
                constrain_water_network(pose, wat_res_idx, lig_res_idx, water_lig_atom)
                print(f"Applied water constraints using ligand atom {water_lig_atom}")
        else:
            print("No water found (checked chain and residue names), skipping water constraints")

    # Relax
    print("\nRunning relax...")
    apply_interface_relax(pose, sf)
    post_relax_total_score = sf(pose)

    # Score with XML
    if os.path.exists(args.xml_path):
        print("Running interface scoring XML...")
        rsp = RosettaScriptsParser()
        parsed_protocol = rsp.generate_mover(args.xml_path)
        parsed_protocol.apply(pose)
        scores = pose.scores
    else:
        print(f"Warning: XML not found at {args.xml_path}, skipping")
        scores = {}

    # Evaluate all polar contacts (DYNAMIC)
    polar_contact_results = {}
    for atom_name, element in polar_atoms.items():
        has_contact = check_polar_contact(pose, lig_res_idx, atom_name)
        polar_contact_results[f"{atom_name}_polar_contact"] = "yes" if has_contact else "no"

    # Evaluate charge satisfaction
    charge_results = check_charge_satisfaction(pose, lig_res_idx, charged_groups)

    # Count H-bonds
    tot_hb = compute_total_hb_to_lig(pose, args.ligand_chain)

    # Compile scores
    def _fmt(val, decimals=2):
        """Format a score value: numeric → rounded string, else 'N/A'."""
        if isinstance(val, (int, float)):
            return f"{val:.{decimals}f}"
        return "N/A"

    score_data = {
        "post_relax_total_score": f"{post_relax_total_score:.2f}",
        "total_hbonds_to_ligand": tot_hb,
        "dG_sep": _fmt(scores.get('dG_separated')),
        "buried_unsatisfied_polars": _fmt(scores.get('delta_unsatHbonds')),
        "shape_complementarity": _fmt(scores.get('sc_value')),
        "dsasa_int": _fmt(scores.get('dSASA_int')),
        "packstat": _fmt(scores.get('packstat'), 4),
        "dSASA_hphobic": _fmt(scores.get('dSASA_hphobic')),
        "dSASA_polar": _fmt(scores.get('dSASA_polar')),
        "nres_int": _fmt(scores.get('nres_int'), 0),
        "hbonds_int": _fmt(scores.get('hbonds_int'), 0),
        "per_residue_energy_int": _fmt(scores.get('per_residue_energy_int'), 3),
    }

    # Add dynamic polar contact columns
    score_data.update(polar_contact_results)

    # Add charge satisfaction
    score_data["charge_satisfied"] = 1 if charge_results['overall_satisfied'] else 0
    score_data["carboxylate_satisfied"] = 1 if charge_results['carboxylate_satisfied'] else 0
    score_data["amine_satisfied"] = 1 if charge_results['amine_satisfied'] else 0
    score_data["sulfonate_satisfied"] = 1 if charge_results['sulfonate_satisfied'] else 0
    score_data["phosphate_satisfied"] = 1 if charge_results['phosphate_satisfied'] else 0

    # Write output
    pose.dump_pdb(args.output_pdb)

    output_score_path = args.output_pdb.replace('.pdb', '_score.sc')
    with open(output_score_path, "w") as f:
        f.write("SCORES:\n")
        for key, value in score_data.items():
            f.write(f"{key}: {value}\n")

    print(f"\n✓ Relaxed structure written to: {args.output_pdb}")
    print(f"✓ Scores written to: {output_score_path}")
    print(f"\nPolar contact summary:")
    for key, value in polar_contact_results.items():
        print(f"  {key}: {value}")
    print(f"\nCharge satisfaction: {score_data['charge_satisfied']}")


if __name__ == "__main__":
    main()
