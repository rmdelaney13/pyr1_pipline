# -*- coding: utf-8 -*-
import sys
import os
import pyrosetta
from pyrosetta import init, pose_from_file, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (RestrictToRepacking, IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT )
from pyrosetta.rosetta.core.select.residue_selector import (NeighborhoodResidueSelector, OrResidueSelector, ResidueIndexSelector, NotResidueSelector)
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet, fill_hbond_set
from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, angle_constraint

# =============================================================================
# CONFIGURATION & CONSTANTS
# =============================================================================

DESIGN_RESIDUES = [
    "A59", "A79", "A81", "A90", "A92", "A106", "A108", "A115", 
    "A118", "A120", "A139", "A141", "A157", "A158", "A161", "A162", "A165"
]
EXTRA_RELAX_RESIDUES = ["C360"] 

AA_1_TO_3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS',
             'L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

# =============================================================================
# CORE HELPERS
# =============================================================================
# -*- coding: utf-8 -*-
import sys
import os
import pyrosetta
from pyrosetta import init, pose_from_file, create_score_function
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (RestrictToRepacking, IncludeCurrent, OperateOnResidueSubset, PreventRepackingRLT )
# UPDATED IMPORTS HERE
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, 
    NeighborhoodResidueSelector, 
    AndResidueSelector, 
    NotResidueSelector, 
    OrResidueSelector, 
    ResidueIndexSelector
)
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

AA_1_TO_3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS',
             'L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

# =============================================================================
# HELPERS
# =============================================================================

def get_indices(pose):
    lig_idx, wat_idx = None, None
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i)
        if res.is_ligand():
            if res.name3() in ["HOH", "TP3", "WAT"]:
                wat_idx = i
            else:
                lig_idx = i
    return lig_idx, wat_idx

def thread_sequence(pose, full_sequence):
    seq_parts = full_sequence.split(':')
    pdb_info = pose.pdb_info()
    current_part = 0
    seen_chains = []
    for i in range(1, pose.total_residue() + 1):
        chain = pdb_info.chain(i)
        if chain not in seen_chains:
            seen_chains.append(chain)
            chain_res = [j for j in range(1, pose.total_residue()+1) if pdb_info.chain(j)==chain]
            if not pose.residue(chain_res[0]).is_protein(): continue
            if current_part < len(seq_parts):
                target = seq_parts[current_part]
                for k, res_idx in enumerate(chain_res):
                    if k < len(target) and pose.residue(res_idx).name1() != target[k]:
                        MutateResidue(res_idx, AA_1_TO_3.get(target[k], "ALA")).apply(pose)
                current_part += 1

def get_ligand_metrics(pose, ligand_idx):
    if not ligand_idx: return 0, "N/A"
    hb_set = HBondSet()
    fill_hbond_set(pose, False, hb_set)
    lig_res = pose.residue(ligand_idx)
    polar_indices = [i for i in range(1, lig_res.natoms() + 1) if lig_res.atom_type(i).element() in ["N", "O"]]
    satisfied_indices = set()
    for i in range(1, hb_set.nhbonds() + 1):
        hb = hb_set.hbond(i)
        if hb.don_res() == ligand_idx: satisfied_indices.add(hb.don_hatm())
        if hb.acc_res() == ligand_idx: satisfied_indices.add(hb.acc_atm())
    unsat_count = len([idx for idx in polar_indices if idx not in satisfied_indices])
    charge_satisfied = "Yes"
    charged_found = False
    for i in range(1, lig_res.natoms() + 1):
        q = lig_res.atomic_charge(i)
        if abs(q) < 0.5: continue
        charged_found = True
        counter_ion_present = False
        for r_idx in range(1, pose.total_residue() + 1):
            if r_idx == ligand_idx: continue
            rec = pose.residue(r_idx)
            if lig_res.xyz(i).distance(rec.nbr_atom_xyz()) < 4.5:
                if q < -0.5 and rec.name3() in ["ARG", "LYS", "HIS"]: counter_ion_present = True; break
                if q > 0.5 and rec.name3() in ["ASP", "GLU"]: counter_ion_present = True; break
        if not counter_ion_present: charge_satisfied = "No"
    if not charged_found: charge_satisfied = "Neutral"
    return unsat_count, charge_satisfied

# =============================================================================
# RELAX & CONSTRAINTS
# =============================================================================

def apply_water_constraints(pose, wat_idx, lig_idx):
    if not wat_idx or not lig_idx: return
    lig = pose.residue(lig_idx)
    wat_o_xyz = pose.residue(wat_idx).xyz("O")
    best_atom, min_d = None, 99.0
    for i in range(1, lig.natoms() + 1):
        if lig.atom_type(i).element() == "O":
            d = wat_o_xyz.distance(lig.xyz(i))
            if d < min_d: min_d = d; best_atom = lig.atom_name(i).strip()
    if not best_atom: return
    def add_dist_cst(r1, a1, r2, a2, dist):
        aid1 = AtomID(pose.residue(r1).atom_index(a1), r1)
        aid2 = AtomID(pose.residue(r2).atom_index(a2), r2)
        pose.add_constraint(AtomPairConstraint(aid1, aid2, HarmonicFunc(dist, 0.5)))
    p_info = pose.pdb_info()
    a86, a114, c360 = p_info.pdb2pose("A", 86), p_info.pdb2pose("A", 114), p_info.pdb2pose("C", 360)
    if all([a86, a114, c360]):
        add_dist_cst(wat_idx, "H1", a86, "O", 2.0)
        add_dist_cst(wat_idx, "O", a114, "H", 2.0)
        add_dist_cst(wat_idx, "O", c360, "HE1", 2.0)
        add_dist_cst(wat_idx, "H2", lig_idx, best_atom, 2.0)

def run_relax(pose, sf, lig_idx, wat_idx):
    mm = MoveMap()
    p_info = pose.pdb_info()
    
    # 1. Ligand + Water Selection
    lig_sel = ResidueIndexSelector(str(lig_idx))
    wat_sel = ResidueIndexSelector(str(wat_idx)) if wat_idx else None
    
    # 2. Neighborhood Shell (8.0A) restricted to Chain A
    core_sel = OrResidueSelector(lig_sel, wat_sel) if wat_sel else lig_sel
    nbr_sel = NeighborhoodResidueSelector(core_sel, 8.0, False)
    chain_a_sel = ChainSelector("A")
    shell_chain_a = AndResidueSelector(nbr_sel, chain_a_sel)
    
    # 3. Design Residues + C360
    design_sel = ResidueIndexSelector()
    for item in DESIGN_RESIDUES:
        idx = p_info.pdb2pose(item[0], int(item[1:]))
        if idx: design_sel.append_index(idx)
    
    c360_idx = p_info.pdb2pose("C", 360)
    c360_sel = ResidueIndexSelector(str(c360_idx)) if c360_idx else None

    # Final Whitelist
    final_sel = OrResidueSelector(shell_chain_a, design_sel)
    final_sel = OrResidueSelector(final_sel, core_sel)
    if c360_sel: final_sel = OrResidueSelector(final_sel, c360_sel)
        
    subset = final_sel.apply(pose)
    mm.set_bb(subset)
    mm.set_chi(subset)
    mm.set_jump(True)
    
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(final_sel)))
    
    relax = FastRelax(sf, 5)
    relax.set_task_factory(tf)
    relax.set_movemap(mm)
    relax.apply(pose)

# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) < 5:
        print("Usage: python relax.py template.pdb seqs.fa lig.params prefix")
        return

    init(f"-extra_res_fa {sys.argv[3]} -ex1 -ex2 -relax:fast -corrections::beta_nov16 true -mute all")
    template = pose_from_file(sys.argv[1])
    sf = create_score_function("beta_nov16")
    sf.set_weight(atom_pair_constraint, 0.25)
    sf.set_weight(angle_constraint, 0.25)

    xml_path = "/projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/rosetta/interface_scoring.xml"
    xml_mover = RosettaScriptsParser().generate_mover(xml_path)

    seqs = []
    with open(sys.argv[2], 'r') as f:
        curr_id = "native"
        for line in f:
            if line.startswith(">"):
                curr_id = line.split('id=')[1].split(',')[0].strip() if "id=" in line else "native"
            else:
                seqs.append((curr_id, line.strip()))

    for design_id, sequence in seqs:
        work_pose = template.clone()
        thread_sequence(work_pose, sequence)
        l_idx, w_idx = get_indices(work_pose)
        
        apply_water_constraints(work_pose, w_idx, l_idx)
        run_relax(work_pose, sf, l_idx, w_idx)
        
        xml_mover.apply(work_pose)
        unsat_polars, charge_stat = get_ligand_metrics(work_pose, l_idx)
        
        scores = {
            "dG_sep": work_pose.scores.get("dG_separated", "N/A"),
            "lig_unsat_polars": unsat_polars,
            "charge_satisfied": charge_stat,
            "interface_unsats": work_pose.scores.get("delta_unsatHbonds", "N/A")
        }
        
        out_name = f"{sys.argv[4]}_{design_id}"
        work_pose.dump_pdb(f"{out_name}.pdb")
        with open(f"{out_name}.sc", "w") as f:
            for k, v in scores.items(): f.write(f"{k}: {v}\n")

if __name__ == "__main__":
    main()
def get_indices(pose):
    lig_idx, wat_idx = None, None
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i)
        if res.is_ligand():
            if res.name3() in ["HOH", "TP3", "WAT"]:
                wat_idx = i
            else:
                lig_idx = i
    return lig_idx, wat_idx

def thread_sequence(pose, full_sequence):
    seq_parts = full_sequence.split(':')
    pdb_info = pose.pdb_info()
    current_part = 0
    seen_chains = []
    for i in range(1, pose.total_residue() + 1):
        chain = pdb_info.chain(i)
        if chain not in seen_chains:
            seen_chains.append(chain)
            chain_res = [j for j in range(1, pose.total_residue()+1) if pdb_info.chain(j)==chain]
            if not pose.residue(chain_res[0]).is_protein(): continue
            if current_part < len(seq_parts):
                target = seq_parts[current_part]
                for k, res_idx in enumerate(chain_res):
                    if k < len(target) and pose.residue(res_idx).name1() != target[k]:
                        MutateResidue(res_idx, AA_1_TO_3.get(target[k], "ALA")).apply(pose)
                current_part += 1

# =============================================================================
# LIGAND METRICS (GENERALIZED)
# =============================================================================

def get_ligand_metrics(pose, ligand_idx):
    if not ligand_idx: return 0, "N/A"
    
    hb_set = HBondSet()
    fill_hbond_set(pose, False, hb_set)
    lig_res = pose.residue(ligand_idx)
    
    # 1. Unsatisfied Polars based on Elements N/O
    polar_indices = [i for i in range(1, lig_res.natoms() + 1) if lig_res.atom_type(i).element() in ["N", "O"]]
    
    satisfied_indices = set()
    for i in range(1, hb_set.nhbonds() + 1):
        hb = hb_set.hbond(i)
        if hb.don_res() == ligand_idx: satisfied_indices.add(hb.don_hatm())
        if hb.acc_res() == ligand_idx: satisfied_indices.add(hb.acc_atm())
        
    unsat_count = len([idx for idx in polar_indices if idx not in satisfied_indices])

    # 2. Charge Satisfaction (Formal/Strong Partial Charges)
    # Threshold 0.5 picks up your OOC (-0.7) and COO (+0.68) atoms
    charge_satisfied = "Yes"
    charged_found = False
    for i in range(1, lig_res.natoms() + 1):
        q = lig_res.atomic_charge(i)
        if abs(q) < 0.5: continue
        charged_found = True
        
        counter_ion_present = False
        for r_idx in range(1, pose.total_residue() + 1):
            if r_idx == ligand_idx: continue
            rec = pose.residue(r_idx)
            if lig_res.xyz(i).distance(rec.nbr_atom_xyz()) < 4.5:
                # If ligand negative, look for ARG/LYS/HIS. If positive, ASP/GLU.
                if q < -0.5 and rec.name3() in ["ARG", "LYS", "HIS"]: counter_ion_present = True; break
                if q > 0.5 and rec.name3() in ["ASP", "GLU"]: counter_ion_present = True; break
        
        if not counter_ion_present: charge_satisfied = "No"
            
    if not charged_found: charge_satisfied = "Neutral"
    return unsat_count, charge_satisfied

# =============================================================================
# CONSTRAINTS & RELAX
# =============================================================================

def apply_water_constraints(pose, wat_idx, lig_idx):
    if not wat_idx or not lig_idx: return
    
    # Find nearest ligand oxygen to water O
    lig = pose.residue(lig_idx)
    wat_o_xyz = pose.residue(wat_idx).xyz("O")
    best_atom, min_d = None, 99.0
    for i in range(1, lig.natoms() + 1):
        if lig.atom_type(i).element() == "O":
            d = wat_o_xyz.distance(lig.xyz(i))
            if d < min_d:
                min_d = d; best_atom = lig.atom_name(i).strip()
    
    if not best_atom: return

    def add_dist_cst(r1, a1, r2, a2, dist):
        aid1 = AtomID(pose.residue(r1).atom_index(a1), r1)
        aid2 = AtomID(pose.residue(r2).atom_index(a2), r2)
        pose.add_constraint(AtomPairConstraint(aid1, aid2, HarmonicFunc(dist, 0.5)))

    p_info = pose.pdb_info()
    a86, a114, c360 = p_info.pdb2pose("A", 86), p_info.pdb2pose("A", 114), p_info.pdb2pose("C", 360)
    
    if all([a86, a114, c360]):
        add_dist_cst(wat_idx, "H1", a86, "O", 2.0)
        add_dist_cst(wat_idx, "O", a114, "H", 2.0)
        add_dist_cst(wat_idx, "O", c360, "HE1", 2.0)
        add_dist_cst(wat_idx, "H2", lig_idx, best_atom, 2.0)

def run_relax(pose, sf, lig_idx):
    mm = MoveMap()
    p_info = pose.pdb_info()
    
    # 1. Select around ligand (8.0A is usually the "sweet spot")
    lig_sel = ResidueIndexSelector(str(lig_idx))
    nbr_sel = NeighborhoodResidueSelector(lig_sel, 7.0, True)
    
    # 2. Filter for Chain A only
    chain_a_sel = ChainSelector("A")
    shell_chain_a = AndResidueSelector(nbr_sel, chain_a_sel)
    
    # 3. Include design residues (forced)
    design_sel = ResidueIndexSelector()
    for item in DESIGN_RESIDUES:
        idx = p_info.pdb2pose(item[0], int(item[1:]))
        if idx: design_sel.append_index(idx)
        
    # 4. Include specific Chain C residue (C360)
    c360_idx = p_info.pdb2pose("C", 360)
    c360_sel = ResidueIndexSelector(str(c360_idx)) if c360_idx else None
    
    # Combine all: (Shell in A) OR (Design residues) OR (Ligand) OR (C360)
    final_sel = OrResidueSelector(shell_chain_a, design_sel)
    final_sel = OrResidueSelector(final_sel, lig_sel)
    if c360_sel:
        final_sel = OrResidueSelector(final_sel, c360_sel)
        
    subset = final_sel.apply(pose)
    
    # Apply to MoveMap
    mm.set_bb(subset)
    mm.set_chi(subset)
    mm.set_jump(True)
    
    # TaskFactory to restrict packing to this same subset
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())
    tf.push_back(RestrictToRepacking())
    tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), NotResidueSelector(final_sel)))
    
    relax = FastRelax(sf, 5)
    relax.set_task_factory(tf)
    relax.set_movemap(mm)
    relax.apply(pose)
    
    
# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) < 5:
        print("Usage: python relax.py template.pdb seqs.fa lig.params prefix")
        return

    init(f"-extra_res_fa {sys.argv[3]} -ex1 -ex2 -relax:fast -corrections::beta_nov16 true -mute all")
    
    template = pose_from_file(sys.argv[1])
    sf = create_score_function("beta_nov16")
    sf.set_weight(atom_pair_constraint, 0.25)
    sf.set_weight(angle_constraint, 0.25)

    # XML for interface scoring
    xml_path = "/projects/ryde3462/software/LigandMPNN/ligand_alignment_mpnn/rosetta/interface_scoring.xml"
    parser = RosettaScriptsParser()
    xml_mover = parser.generate_mover(xml_path)

    # Get sequences from FASTA
    seqs = []
    with open(sys.argv[2], 'r') as f:
        curr_id = "native"
        for line in f:
            if line.startswith(">"):
                curr_id = line.split('id=')[1].split(',')[0].strip() if "id=" in line else "native"
            else:
                seqs.append((curr_id, line.strip()))

    for design_id, sequence in seqs:
        work_pose = template.clone()
        thread_sequence(work_pose, sequence)
        
        l_idx, w_idx = get_indices(work_pose)
        apply_water_constraints(work_pose, w_idx, l_idx)
        run_relax(work_pose, sf, l_idx)
        
        # Final Scoring
        xml_mover.apply(work_pose)
        unsat_polars, charge_stat = get_ligand_metrics(work_pose, l_idx)
        
        scores = {
            "dG_sep": work_pose.scores.get("dG_separated", "N/A"),
            "lig_unsat_polars": unsat_polars,
            "charge_satisfied": charge_stat,
            "interface_unsats": work_pose.scores.get("delta_unsatHbonds", "N/A")
        }
        
        out_pdb = f"{sys.argv[4]}_{design_id}.pdb"
        work_pose.dump_pdb(out_pdb)
        with open(out_pdb.replace(".pdb", ".sc"), "w") as f:
            for k, v in scores.items(): f.write(f"{k}: {v}\n")

if __name__ == "__main__":
    main()
