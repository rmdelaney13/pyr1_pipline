# -*- coding: utf-8 -*-
"""
Usage:
    python relax_threaded.py template.pdb sequences.fa ligand.params output_prefix
"""
import sys
import os
import math
import pyrosetta
from pyrosetta import init, pose_from_file, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (RestrictToRepacking, IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT )
from pyrosetta.rosetta.core.select.residue_selector import (ChainSelector, NeighborhoodResidueSelector, AndResidueSelector, NotResidueSelector, OrResidueSelector, ResidueIndexSelector)
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint, AngleConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols import jd2
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, angle_constraint, buried_unsatisfied_penalty, fa_elec, fa_intra_rep, fa_dun

# =============================================================================
# 1. FASTA & THREADING HELPERS
# =============================================================================

AA_1_TO_3 = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

def get_sequences_from_fasta(fasta_path):
    seqs = []
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    current_id = "native"
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if "id=" in line:
                parts = line.split(',')
                for p in parts:
                    if "id=" in p:
                        current_id = "design_" + p.split('=')[1].strip()
                        break
            else:
                current_id = "native"
        else:
            if line:
                seqs.append((current_id, line))
    return seqs

def thread_sequence(pose, full_sequence):
    seq_parts = full_sequence.split(':')
    pdb_info = pose.pdb_info()
    current_seq_part_idx = 0
    seen_chains = []
    
    for i in range(1, pose.total_residue() + 1):
        chain = pdb_info.chain(i)
        if chain not in seen_chains:
            seen_chains.append(chain)
            chain_res_indices = [j for j in range(1, pose.total_residue() + 1) if pdb_info.chain(j) == chain]
            
            if not pose.residue(chain_res_indices[0]).is_protein():
                continue

            if current_seq_part_idx < len(seq_parts):
                target_seq = seq_parts[current_seq_part_idx]
                if len(target_seq) == len(chain_res_indices):
                    for k, res_idx in enumerate(chain_res_indices):
                        target_aa = target_seq[k]
                        current_name1 = pose.residue(res_idx).name1()
                        if current_name1 != target_aa:
                            try:
                                new_name3 = AA_1_TO_3.get(target_aa, "ALA") 
                                MutateResidue(res_idx, new_name3).apply(pose)
                            except Exception as e:
                                print(f"Error mutating res {res_idx} to {target_aa}: {e}")
                    current_seq_part_idx += 1

# =============================================================================
# 2. SELECTION HELPERS (UPDATED)
# =============================================================================

def get_indices(pose):
    """
    Finds Ligand and Water indices dynamically based on names/chains.
    """
    ligand_idx = None
    water_idx = None
    
    # Iterate through all residues
    for i in range(1, pose.total_residue() + 1):
        name = pose.residue(i).name3()
        chain = pose.pdb_info().chain(i)
        
        # Check for Water (TP3 or HOH)
        if name in ["TP3", "HOH", "WAT"]:
            water_idx = i
        
        # Check for Ligand (Chain B, but NOT the water)
        elif chain == "B" and name not in ["TP3", "HOH", "WAT"]:
            ligand_idx = i
            
    return ligand_idx, water_idx

# =============================================================================
# 3. CONSTRAINT & RELAX HELPERS
# =============================================================================

def find_ligand_oxygen(pose, water_res, ligand_res):
    lig = pose.residue(ligand_res)
    if not lig.natoms(): return None
    wat = pose.residue(water_res)
    if not wat.has("H2"): return None
    wat_H2_xyz = wat.xyz("H2")
    best_atom = None
    best_dist = float('inf')
    for idx in range(1, lig.natoms() + 1):
        if lig.atom_type(idx).element() == "O":
            name = lig.atom_name(idx).strip()
            d = wat_H2_xyz.distance(lig.xyz(name))
            if d < best_dist:
                best_dist = d; best_atom = name
    return best_atom

def constrain_water_network(pose, water_res, ligand_res):
    pdb_info = pose.pdb_info()
    lig_o = find_ligand_oxygen(pose, water_res, ligand_res)
    if not lig_o: return

    def add_dist(r1, a1, r2, a2, dist=2.0, sd=0.5):
        aid1 = AtomID(pose.residue(r1).atom_index(a1), r1)
        aid2 = AtomID(pose.residue(r2).atom_index(a2), r2)
        pose.add_constraint(AtomPairConstraint(aid1, aid2, HarmonicFunc(dist, sd)))

    def add_ang(r1, a1, r2, a2, r3, a3, angle=180.0, sd=20.0):
        aid1 = AtomID(pose.residue(r1).atom_index(a1), r1)
        aid2 = AtomID(pose.residue(r2).atom_index(a2), r2)
        aid3 = AtomID(pose.residue(r3).atom_index(a3), r3)
        pose.add_constraint(AngleConstraint(aid1, aid2, aid3, HarmonicFunc(angle, sd)))

    def find_res(c, n):
        for i in range(1, pose.total_residue()+1):
            if pdb_info.chain(i) == c and pdb_info.number(i) == n: return i
        return None

    a86 = find_res("A", 86); a114 = find_res("A", 114); c360 = find_res("C", 360)

    if a86 and a114 and c360:
        add_dist(water_res, "H1", a86, "O")
        add_dist(water_res, "O", a114, "H")
        add_dist(water_res, "O", c360, "HE1")
        add_dist(water_res, "H2", ligand_res, lig_o, dist=2.0, sd=0.5)
        add_ang(a86, "O", water_res, "H1", water_res, "O")
        if pose.residue(a114).has("N"): add_ang(a114, "N", a114, "H", water_res, "O", 180, 40)
        if pose.residue(c360).has("NE1"): add_ang(c360, "NE1", c360, "HE1", water_res, "O")
        add_ang(ligand_res, lig_o, water_res, "H2", water_res, "O")
        add_ang(water_res, "H1", water_res, "O", water_res, "H2", 104.5, 5.0)

def apply_interface_relax(pose, scorefxn, ligand_idx, water_idx):
    """
    Relax using split shell sizes:
      - 11.0 A for Chain A
      - 6.0 A for Chain C
    """
    if ligand_idx is None:
        print("Error: No ligand found for relax selection!")
        return
        
    # 1. Select Ligand/Water
    ligand_sel = ResidueIndexSelector(ligand_idx)
    
    if water_idx:
        water_sel = ResidueIndexSelector(water_idx)
        lig_plus_water = OrResidueSelector(ligand_sel, water_sel)
    else:
        lig_plus_water = ligand_sel

    # 2. Define Shells by Chain
    # Shell for Chain A (Wide 11.0 A)
    # We select neighbors of ligand/water within 11.0 A, then RESTRICT to Chain A
    bubble_A = NeighborhoodResidueSelector(lig_plus_water, 11.0)
    shell_A = AndResidueSelector(bubble_A, ChainSelector("A"))

    # Shell for Chain C (Narrow 6.0 A)
    # We select neighbors of ligand/water within 6.0 A, then RESTRICT to Chain C
    bubble_C = NeighborhoodResidueSelector(lig_plus_water, 4.0)
    shell_C = AndResidueSelector(bubble_C, ChainSelector("C"))

    # Combine A and C shells
    interface_shell = OrResidueSelector(shell_A, shell_C)

    # Final Shell: Interface (A+C) + Ligand + Water themselves
    shell = OrResidueSelector(interface_shell, lig_plus_water)
    
    # 3. Standard Setup (MoveMap, FoldTree, TaskFactory)
    mm = MoveMap()
    subset = shell.apply(pose)
    mm.set_bb(subset)
    mm.set_chi(subset)
    
    ft = pose.fold_tree()
    pdb_info = pose.pdb_info()
    lig_chain = pdb_info.chain(ligand_idx)
    wat_chain = pdb_info.chain(water_idx) if water_idx else ""

    for j in range(1, ft.num_jump()+1):
        down = ft.downstream_jump_residue(j)
        c = pdb_info.chain(down)
        if c == lig_chain or (water_idx and c == wat_chain):
            mm.set_jump(j, True)

    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    tf.push_back(RestrictToRepacking())
    
    # Freeze everything NOT in our custom shell
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(shell)))

    relax = FastRelax(scorefxn)
    relax.set_task_factory(tf)
    relax.set_movemap(mm)
    relax.dualspace(False)
    relax.apply(pose)

# =============================================================================
# 4. SCORING HELPERS
# =============================================================================

def check_polar_contact_with_ligand_atom(pose, ligand_res, atom_name, max_dist=3.5):
    polar_atoms = {"N", "O"}
    if not pose.residue(ligand_res).has(atom_name): return False
    lig_xyz = pose.residue(ligand_res).xyz(atom_name)
    for i in range(1, pose.total_residue() + 1):
        if i == ligand_res: continue
        res = pose.residue(i)
        for j in range(1, res.natoms() + 1):
            if res.atom_type(j).element() in polar_atoms:
                if lig_xyz.distance(res.xyz(j)) <= max_dist: return True
    return False

def compute_total_hb_to_lig(pose, ligand_idx):
    if not ligand_idx: return 0
    hb_set = HBondSet()
    fill_hbond_set(pose, False, hb_set)
    count = 0
    for i in range(1, hb_set.nhbonds() + 1):
        h = hb_set.hbond(i)
        if h.don_res() == ligand_idx or h.acc_res() == ligand_idx:
            count += 1
    return count

def check_carboxylic_acid_satisfaction(pose, ligand_idx, receptor_chain="A",
                                     ox_names=["O5", "O4"], carbon_name="C24", pos_threshold=3.5):
    if not ligand_idx: return {"satisfied": False}
    
    ligand_res = pose.residue(ligand_idx)
    atoms = [ligand_res.atom_name(j).strip() for j in range(1, ligand_res.natoms() + 1)]
    if not (all(ox in atoms for ox in ox_names) and carbon_name in atoms):
        return {"satisfied": False}

    ox_idx = {ox: ligand_res.atom_index(ox) for ox in ox_names if ligand_res.has(ox)}
    hb_set = HBondSet()
    fill_hbond_set(pose, False, hb_set)
    hbcount = 0
    for i in range(1, hb_set.nhbonds() + 1):
        h = hb_set.hbond(i)
        if h.acc_res() == ligand_idx and h.acc_atm() in ox_idx.values():
            d, a = h.get_HAdist(pose), h.get_AHDangle(pose)
            if d < 3.0 and a > 120:
                hbcount += 1
                
    sat_sb = False
    pdb_info = pose.pdb_info()
    pos_res = {"ARG": ["NH1", "NH2"], "LYS": ["NZ"]}
    for ox, idx in ox_idx.items():
        pos = ligand_res.xyz(idx)
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.chain(i) == receptor_chain:
                r = pose.residue(i)
                nm = r.name3()
                if nm in pos_res:
                    for at in pos_res[nm]:
                        if r.has(at):
                            if pos.distance(r.xyz(r.atom_index(at))) <= pos_threshold:
                                sat_sb = True
    return {"satisfied": hbcount >= 2 or sat_sb}

def score_and_dump(pose, outp, sf, xml_mover, ligand_idx, water_idx):
    if xml_mover: 
        xml_mover.apply(pose)
    
    scores = pose.scores
    dG_sep = scores.get("dG_separated", "N/A")
    hbonds = scores.get("hbonds_int", "N/A")
    delta_unsat = scores.get("delta_unsatHbonds", "N/A")
    
    if ligand_idx:
        o1_contact = "yes" if check_polar_contact_with_ligand_atom(pose, ligand_idx, "O1") else "no"
        o2_contact = "yes" if check_polar_contact_with_ligand_atom(pose, ligand_idx, "O2") else "no"
        tot_hb = compute_total_hb_to_lig(pose, ligand_idx)
        charge_sat = 1 if check_carboxylic_acid_satisfaction(pose, ligand_idx).get("satisfied") else 0
    else:
        o1_contact = "no"; o2_contact = "no"; tot_hb = 0; charge_sat = 0

    score_data = {
        "O1_polar_contact": o1_contact,
        "O2_polar_contact": o2_contact,
        "dG_sep": f"{dG_sep:.2f}" if isinstance(dG_sep, float) else dG_sep,
        "total_hbonds_to_ligand": tot_hb,
        "hbond_int": hbonds,
        "charge_satisfied": charge_sat,
        "buried_unsatisfied_polars": f"{delta_unsat:.2f}" if isinstance(delta_unsat, float) else delta_unsat,
    }

    pose.dump_pdb(outp)
    output_score_path = outp.replace('.pdb', '.sc')
    with open(output_score_path, "w") as f:
        f.write("SCORES:\n")
        for key, value in score_data.items():
            f.write(f"{key}: {value}\n")

# =============================================================================
# 5. MAIN LOOP
# =============================================================================

def main():
    if len(sys.argv) != 5:
        print("Usage: python relax_threaded.py template.pdb sequences.fa ligand.params output_prefix")
        sys.exit(1)
        
    template_pdb, fasta_file, params_file, out_prefix = sys.argv[1:]
    
    init(f"-extra_res_fa {params_file} -ex1 -ex2aro -ex2 -use_input_sc -relax:fast -relax:default_repeats 5 -corrections::beta_nov16 true -mute all")

    template_pose = pose_from_file(template_pdb)
    
    sf = create_score_function("beta_nov16")
    
    # Weights for design/relax
    #sf.set_weight(fa_dun, 0.8)   # Default is often ~0.7, increasing to 0.8 helps
    # Boost Intra-residue Repulsion (punish self-clash)
    #sf.set_weight(fa_intra_rep, 0.02) # Default is often 0.002 or lower
    
    # Weights for design/relax
    sf.set_weight(atom_pair_constraint, 0.2)
    sf.set_weight(angle_constraint, 0.2)
    #sf.set_weight(buried_unsatisfied_penalty, 1.4)
    #sf.set_weight(fa_elec, 1.2)
    
    xml_path = "/projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/rosetta/interface_scoring.xml"
    rsp = RosettaScriptsParser()
    xml_mover = rsp.generate_mover(xml_path)

    sequences = get_sequences_from_fasta(fasta_file)
    print(f"Loaded {len(sequences)} sequences from {fasta_file}")

    for i, (design_id, seq) in enumerate(sequences):
        print(f"Processing {design_id}...")
        
        work_pose = template_pose.clone()
        thread_sequence(work_pose, seq)
        
        # FIND LIGAND AND WATER DYNAMICALLY
        lig_idx, wat_idx = get_indices(work_pose)
        
        if not lig_idx:
            print("CRITICAL WARNING: Ligand not found in Chain B! Check input.")
        
        # APPLY CONSTRAINTS (If water exists)
        if lig_idx and wat_idx:
            try:
                constrain_water_network(work_pose, wat_idx, lig_idx)
            except Exception as e:
                print(f"Warning: Failed to apply constraints: {e}")
        else:
            print("Warning: Skipping constraints (Ligand or Water missing).")
            
        # RELAX (Passing indices explicitly)
        apply_interface_relax(work_pose, sf, lig_idx, wat_idx)
        
        out_name = f"{out_prefix}_{design_id}.pdb"
        score_and_dump(work_pose, out_name, sf, xml_mover, lig_idx, wat_idx)

if __name__ == "__main__":
    main()