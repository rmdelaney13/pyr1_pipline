#!/usr/bin/env python
"""
Shell-restricted backbone-constrained FastRelax for threaded mutant structures.

Only relaxes sidechains within a 10 A shell around mutated residues (same
approach as relax_general_universal.py's interface relax). Backbone is
constrained with CoordinateConstraints; everything outside the shell is frozen.

Input:  Protein-only PDB (no ligand). May include water.
Output: Relaxed PDB with optimized sidechains near mutation sites.

Uses ref2015 score function (same as docking stage).

Usage:
    python constrained_relax.py input.pdb output.pdb --mutations "59K;120A;160G"
"""
import argparse
import re
import sys
from pathlib import Path


BACKBONE_ATOMS = ["N", "CA", "C", "O"]


def parse_mutation_positions(signature):
    """
    Parse variant signature like "59K;120A;160G" into PDB residue numbers.

    Returns list of ints, e.g. [59, 120, 160].
    """
    if not signature or signature.strip() == "":
        return []
    positions = []
    for token in signature.split(";"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r'^(\d+)', token)
        if match:
            positions.append(int(match.group(1)))
    return positions


def find_shell_residues(pose, mutation_pdb_numbers, shell_distance=10.0):
    """
    Find pose residue indices within shell_distance of mutated positions.

    Args:
        pose: PyRosetta Pose
        mutation_pdb_numbers: list of PDB residue numbers (chain A)
        shell_distance: neighborhood radius in Angstroms

    Returns:
        set of pose residue indices in the shell (including mutation sites)
    """
    from pyrosetta.rosetta.core.select.residue_selector import (
        ResidueIndexSelector, NeighborhoodResidueSelector
    )

    # Map PDB numbers to pose indices
    pdb_info = pose.pdb_info()
    mutation_pose_indices = []
    for pdb_num in mutation_pdb_numbers:
        pose_idx = pdb_info.pdb2pose('A', pdb_num)
        if pose_idx > 0:
            mutation_pose_indices.append(pose_idx)
        else:
            print(f"  Warning: PDB residue A:{pdb_num} not found in pose")

    if not mutation_pose_indices:
        return set()

    # Build selector: residues within shell_distance of any mutated residue
    index_str = ",".join(str(i) for i in mutation_pose_indices)
    mut_sel = ResidueIndexSelector(index_str)
    shell_sel = NeighborhoodResidueSelector(mut_sel, shell_distance, True)

    # Evaluate selector
    shell_vec = shell_sel.apply(pose)
    shell_set = set()
    for i in range(1, pose.total_residue() + 1):
        if shell_vec[i]:
            shell_set.add(i)

    return shell_set


def add_backbone_coordinate_constraints(pose, residue_indices, sd=0.5):
    """
    Add CoordinateConstraints on backbone heavy atoms (N, CA, C, O)
    for specified residues only.

    Requires that addVirtualResAsRoot(pose) has already been called.
    """
    from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.id import AtomID

    vrt_res = pose.total_residue()
    vrt_atom_id = AtomID(1, vrt_res)

    func = HarmonicFunc(0.0, sd)

    n_constrained = 0
    for i in residue_indices:
        residue = pose.residue(i)
        if not residue.is_protein():
            continue
        for atom_name in BACKBONE_ATOMS:
            if residue.has(atom_name):
                atom_idx = residue.atom_index(atom_name)
                atom_id = AtomID(atom_idx, i)
                xyz = residue.xyz(atom_name)
                constraint = CoordinateConstraint(
                    atom_id, vrt_atom_id, xyz, func
                )
                pose.add_constraint(constraint)
                n_constrained += 1

    print(f"Added {n_constrained} backbone CoordinateConstraints (sd={sd} A)")


def strip_vrt_from_pdb(pdb_path):
    """Remove VRT residue lines from a PDB file in-place."""
    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    filtered = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname in ("VRT", "VIRT", "XXX"):
                continue
        filtered.append(line)

    with open(pdb_path, 'w') as f:
        f.writelines(filtered)


def main():
    parser = argparse.ArgumentParser(
        description="Shell-restricted constrained FastRelax for threaded mutant PDBs"
    )
    parser.add_argument("input_pdb", help="Input PDB (threaded mutant)")
    parser.add_argument("output_pdb", help="Output relaxed PDB")
    parser.add_argument("--mutations", required=True,
                        help="Variant signature, e.g. '59K;120A;160G'")
    parser.add_argument("--params", default=None,
                        help="Optional extra_res_fa params file")
    parser.add_argument("--sd", type=float, default=0.5,
                        help="Harmonic constraint SD in Angstroms (default: 0.5)")
    parser.add_argument("--cycles", type=int, default=5,
                        help="FastRelax cycles (default: 5)")
    parser.add_argument("--shell", type=float, default=10.0,
                        help="Shell radius around mutations in Angstroms (default: 10.0)")
    args = parser.parse_args()

    # Parse mutation positions
    mutation_positions = parse_mutation_positions(args.mutations)
    if not mutation_positions:
        print("No mutations specified - copying input to output unchanged")
        import shutil
        Path(args.output_pdb).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.input_pdb, args.output_pdb)
        return

    print(f"Mutation positions: {mutation_positions}")

    # Ensure output directory exists
    Path(args.output_pdb).parent.mkdir(parents=True, exist_ok=True)

    # Build PyRosetta init flags
    init_flags = "-ex1 -ex2aro -use_input_sc -relax:fast -mute all"
    if args.params:
        init_flags = f"-extra_res_fa {args.params} {init_flags}"

    import pyrosetta
    from pyrosetta import init, pose_from_file
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.core.pose import addVirtualResAsRoot
    from pyrosetta.rosetta.core.scoring import coordinate_constraint

    init(init_flags)

    # Load pose
    print(f"Loading: {args.input_pdb}")
    pose = pose_from_file(args.input_pdb)
    print(f"Loaded {pose.total_residue()} residues")

    # Find shell residues around mutation sites
    shell_residues = find_shell_residues(pose, mutation_positions, args.shell)
    print(f"Shell: {len(shell_residues)} residues within {args.shell} A of mutations")

    if not shell_residues:
        print("WARNING: No shell residues found - copying input unchanged")
        import shutil
        shutil.copy2(args.input_pdb, args.output_pdb)
        return

    # Identify water residues (keep frozen)
    water_names = {'TP3', 'HOH', 'WAT', 'TIP', 'TIP3'}
    water_residues = set()
    for i in range(1, pose.total_residue() + 1):
        if pose.residue(i).name3().strip() in water_names:
            water_residues.add(i)

    # Add virtual residue as fold-tree root (anchor for CoordinateConstraints)
    addVirtualResAsRoot(pose)
    print(f"Added VRT root (now {pose.total_residue()} residues)")

    # Add backbone constraints only within the shell
    add_backbone_coordinate_constraints(pose, shell_residues, sd=args.sd)

    # Build MoveMap: chi-only for shell protein residues, everything else frozen
    # (Same pattern as relax_general_universal.py apply_interface_relax)
    mm = MoveMap()
    for i in range(1, pose.total_residue() + 1):
        if i in water_residues:
            mm.set_bb(i, False)
            mm.set_chi(i, False)
        elif i in shell_residues:
            mm.set_bb(i, False)
            mm.set_chi(i, True)
        else:
            mm.set_bb(i, False)
            mm.set_chi(i, False)

    # Score function: ref2015 + coordinate_constraint weight
    sfxn = pyrosetta.create_score_function("ref2015")
    sfxn.set_weight(coordinate_constraint, 1.0)

    # Run FastRelax with restricted MoveMap
    print(f"Running FastRelax ({args.cycles} cycles, {len(shell_residues)} mobile residues)...")
    relax = FastRelax(sfxn, args.cycles)
    relax.set_movemap(mm)
    relax.apply(pose)

    score_after = sfxn(pose)
    print(f"Post-relax score: {score_after:.2f}")

    # Save output
    pose.dump_pdb(args.output_pdb)
    strip_vrt_from_pdb(args.output_pdb)
    print(f"Saved relaxed structure: {args.output_pdb}")


if __name__ == "__main__":
    main()
