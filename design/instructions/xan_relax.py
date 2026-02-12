# -*- coding: utf-8 -*-
"""
Updated Relax Script for Xanthurenic Acid (XAN)
Fixed: Hydrogen Bond Geometry Reporting using xyz coordinates
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
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, angle_constraint

# =============================================================================
# CONFIGURATION
# =============================================================================

DESIGN_RESIDUES = [
    "A59", "A79", "A81", "A90", "A92", "A106", "A108", "A115", 
    "A118", "A120", "A139", "A141", "A157", "A158", "A161", "A162", "A165"
]

EXTRA_RELAX_RESIDUES = ["C360"] 

AA_1_TO_3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

# =============================================================================
# 1. THREADING & UTILITY
# =============================================================================

def get_sequences_from_fasta(fasta_path):
    seqs = []
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    current_id = "native"
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            current_id = "design_" + line.split('id=')[1].split(',')[0].strip() if "id=" in line else "native"
        elif line:
            seqs.append((current_id, line))
    return seqs

def thread_sequence(pose, full_sequence):
    seq_parts = full_sequence.split(':')
    pdb_info = pose.pdb_info()
    curr_part = 0
    seen_chains = []
    for i in range(1, pose.total_residue() + 1):
        chain = pdb_info.chain(i)
        if chain not in seen_chains:
            seen_chains.append(chain)
            indices = [j for j in range(1, pose.total_residue() + 1) if pdb_info.chain(j) == chain]
            if not pose.residue(indices[0]).is_protein(): continue
            if curr_part < len(seq_parts):
                target_seq = seq_parts[curr_part]
                for k, res_idx in enumerate(indices):
                    if k < len(target_seq) and pose.residue(res_idx).name1() != target_seq[k]:
                        MutateResidue(res_idx, AA_1_TO_3.get(target_seq[k], "ALA")).apply(pose)
                curr_part += 1

def get_indices(pose):
    ligand_idx = water_idx = None
    for i in range(1, pose.total_residue() + 1):
        name = pose.residue(i).name3()
        chain = pose.pdb_info().chain(i)
        if name in ["TP3", "HOH", "WAT"]: water_idx = i
        elif chain == "B" and pose.residue(i).is_ligand(): ligand_idx = i
    return ligand_idx, water_idx

# =============================================================================
# 2. CONSTRAINTS & RELAX
# =============================================================================

def find_ligand_oxygen(pose, water_res, ligand_res):
    lig = pose.residue(ligand_res)
    wat = pose.residue(water_res)
    wat_xyz = wat.xyz("O") 
    best_atom = None
    best_dist = 99.0
    for name in ["O1", "O2", "O3", "O4"]:
        if lig.has(name):
            d = wat_xyz.distance(lig.xyz(name))
            if d < best_dist:
                best_dist = d; best_atom = name
    return best_atom

def constrain_water_network(pose, water_res, ligand_res):
    pdb = pose.pdb_info()
    lig_o = find_ligand_oxygen(pose, water_res, ligand_res)
    if not lig_o: return

    def add_dist(r1, a1, r2, a2, dist=2.0, sd=0.4):
        aid1 = AtomID(pose.residue(r1).atom_index(a1), r1)
        aid2 = AtomID(pose.residue(r2).atom_index(a2), r2)
        pose.add_constraint(AtomPairConstraint(aid1, aid2, HarmonicFunc(dist, sd)))

    a86 = pdb.pdb2pose("A", 86); a114 = pdb.pdb2pose("A", 114); c360 = pdb.pdb2pose("C", 360) 

    if a86 and a114 and c360:
        add_dist(water_res, "H1", a86, "O")
        add_dist(water_res, "O", a114, "H")
        add_dist(water_res, "O", c360, "HE1")
        add_dist(water_res, "H2", ligand_res, lig_o, dist=2.0, sd=0.4)

def apply_interface_relax(pose, scorefxn, ligand_idx, water_idx, design_residues, extra_residues):
    if ligand_idx is None: return
    pdb_info = pose.pdb_info()
    lig_sel = ResidueIndexSelector(ligand_idx)
    water_sel = ResidueIndexSelector(water_idx) if water_idx else lig_sel
    lig_plus_water = OrResidueSelector(lig_sel, water_sel)
    ligand_neighbors = NeighborhoodResidueSelector(lig_plus_water, 10.0, False)
    
    design_indices = [pdb_info.pdb2pose(item[0], int(item[1:])) for item in design_residues if pdb_info.pdb2pose(item[0], int(item[1:])) != 0]
    
    if design_indices:
        design_sel = ResidueIndexSelector(",".join(map(str, design_indices)))
        design_bubble = NeighborhoodResidueSelector(design_sel, 6.0, True) 
        candidate_shell = OrResidueSelector(ligand_neighbors, design_bubble)
    else:
        candidate_shell = ligand_neighbors

    final_shell = AndResidueSelector(candidate_shell, ChainSelector("A"))
    final_shell = OrResidueSelector(final_shell, lig_plus_water)
    for item in extra_residues:
        idx = pdb_info.pdb2pose(item[0], int(item[1:]))
        if idx: final_shell = OrResidueSelector(final_shell, ResidueIndexSelector(idx))

    mm = MoveMap()
    subset = final_shell.apply(pose)
    mm.set_bb(subset); mm.set_chi(subset)
    for j in range(1, pose.fold_tree().num_jump()+1): mm.set_jump(j, True)

    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(final_shell)))

    relax = FastRelax(scorefxn, 5)
    relax.set_task_factory(tf)
    relax.set_movemap(mm)
    relax.apply(pose)

# =============================================================================
# 3. GEOMETRY & SCORING HELPERS
# =============================================================================

def get_ligand_hbond_geometry(pose, ligand_idx, atom_names=["O1", "O2", "O3", "O4"]):
    """Extracts D-H...A distance and angle using explicit xyz coordinates."""
    geo_results = {f"{n}_{k}": 0.0 for n in atom_names for k in ["dist", "angle"]}
    if not ligand_idx: return geo_results
    
    hbond_set = HBondSet()
    fill_hbond_set(pose, False, hbond_set)
    lig_res = pose.residue(ligand_idx)

    for i in range(1, hbond_set.nhbonds() + 1):
        hb = hbond_set.hbond(i)
        if hb.acc_res() == ligand_idx:
            acc_atm_idx = hb.acc_atm()
            acc_atm_name = lig_res.atom_name(acc_atm_idx).strip()
            
            if acc_atm_name in atom_names:
                # don_item() returns the index of the heavy atom donor
                # don_hatm() returns the index of the hydrogen atom
                d_idx = hb.don_item()
                h_idx = hb.don_hatm()
                
                d_xyz = pose.xyz(AtomID(d_idx, hb.don_res()))
                h_xyz = pose.xyz(AtomID(h_idx, hb.don_res()))
                a_xyz = pose.xyz(AtomID(acc_atm_idx, hb.acc_res()))

                dist = h_xyz.distance(a_xyz)
                # D-H...A angle calculation
                angle_rad = pyrosetta.rosetta.numeric.angle_radians(d_xyz, h_xyz, a_xyz)
                angle_deg = math.degrees(angle_rad)

                current_best = geo_results.get(f"{acc_atm_name}_dist", 99.0)
                if current_best == 0.0 or dist < current_best:
                    geo_results[f"{acc_atm_name}_dist"] = round(dist, 2)
                    geo_results[f"{acc_atm_name}_angle"] = round(angle_deg, 1)
    return geo_results

def check_polar_contact_with_ligand_atom(pose, ligand_res, atom_name, max_dist=3.5):
    if not pose.residue(ligand_res).has(atom_name): return False
    lig_xyz = pose.residue(ligand_res).xyz(atom_name)
    for i in range(1, pose.total_residue() + 1):
        if i == ligand_res: continue
        res = pose.residue(i)
        for j in range(1, res.natoms() + 1):
            a_type = res.atom_type(j)
            if a_type.is_acceptor() or a_type.is_donor():
                if lig_xyz.distance(res.xyz(j)) <= max_dist: return True
    return False

def check_carboxylate_satisfaction(pose, ligand_idx, receptor_chain="A", ox_names=["O3", "O4"]):
    if not ligand_idx: return False
    lig = pose.residue(ligand_idx)
    hb_set = HBondSet(); fill_hbond_set(pose, False, hb_set)
    ox_indices = [lig.atom_index(ox) for ox in ox_names if lig.has(ox)]
    
    hb_count = sum(1 for h in range(1, hb_set.nhbonds()+1) 
                   if hb_set.hbond(h).acc_res() == ligand_idx and hb_set.hbond(h).acc_atm() in ox_indices)
    
    sb_found = False
    pos_res = {"ARG": ["NH1", "NH2"], "LYS": ["NZ"]}
    for idx in ox_indices:
        oxy_xyz = lig.xyz(idx)
        for i in range(1, pose.total_residue() + 1):
            res = pose.residue(i)
            if res.name3() in pos_res and pose.pdb_info().chain(i) == receptor_chain:
                for at in pos_res[res.name3()]:
                    if res.has(at) and oxy_xyz.distance(res.xyz(at)) < 3.8:
                        sb_found = True; break
    return hb_count >= 2 or sb_found

def score_and_dump(pose, outp, sf, xml_mover, ligand_idx, water_idx):
    if xml_mover: xml_mover.apply(pose)
    scores = pose.scores
    
    contact_data = {}
    hbond_geo = {}
    if ligand_idx:
        for ox in ["O1", "O2", "O3", "O4"]:
            contact_data[f"{ox}_contact"] = "yes" if check_polar_contact_with_ligand_atom(pose, ligand_idx, ox) else "no"
        charge_sat = 1 if check_carboxylate_satisfaction(pose, ligand_idx) else 0
        hbond_geo = get_ligand_hbond_geometry(pose, ligand_idx)
    else:
        contact_data = {f"{ox}_contact": "no" for ox in ["O1", "O2", "O3", "O4"]}
        charge_sat = 0
        hbond_geo = {f"{ox}_{k}": 0.0 for ox in ["O1", "O2", "O3", "O4"] for k in ["dist", "angle"]}

    score_data = {
        **contact_data,
        **hbond_geo,
        "dG_sep": scores.get("dG_separated", "N/A"),
        "charge_satisfied": charge_sat,
    }
    
    pose.dump_pdb(outp)
    with open(outp.replace('.pdb', '.sc'), "w") as f:
        for k, v in score_data.items(): f.write(f"{k}: {v}\n")

# =============================================================================
# 4. MAIN
# =============================================================================

def main():
    if len(sys.argv) != 5:
        print("Usage: python xan_relax.py template.pdb sequences.fa ligand.params output_prefix")
        sys.exit(1)

    template_pdb, fasta_file, params_file, out_prefix = sys.argv[1:]
    init(f"-extra_res_fa {params_file} -ex1 -ex2 -use_input_sc -relax:fast -corrections::beta_nov16 true -mute all")

    template_pose = pose_from_file(template_pdb)
    sf = create_score_function("beta_nov16")
    sf.set_weight(atom_pair_constraint, 0.3)
    # Note: angle_constraint weight is only needed if you specifically add AngleConstraints 
    # but having it set doesn't hurt.

    xml_path = "/projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/rosetta/interface_scoring.xml"
    xml_mover = RosettaScriptsParser().generate_mover(xml_path) if os.path.exists(xml_path) else None

    sequences = get_sequences_from_fasta(fasta_file)
    for design_id, seq in sequences:
        print(f"Processing {design_id}...")
        work_pose = template_pose.clone()
        thread_sequence(work_pose, seq)
        lig_idx, wat_idx = get_indices(work_pose)
        if lig_idx and wat_idx:
            constrain_water_network(work_pose, wat_idx, lig_idx)
        apply_interface_relax(work_pose, sf, lig_idx, wat_idx, DESIGN_RESIDUES, EXTRA_RELAX_RESIDUES)
        score_and_dump(work_pose, f"{out_prefix}_{design_id}.pdb", sf, xml_mover, lig_idx, wat_idx)

if __name__ == "__main__":
    main()