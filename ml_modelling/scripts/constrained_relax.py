#!/usr/bin/env python
"""
Backbone-constrained FastRelax for threaded mutant structures.

Applies CoordinateConstraints on backbone heavy atoms (N, CA, C, O) with
HarmonicFunc(0.0, 0.5 A) to relieve steric strain from mutation threading
while keeping backbone geometry close to the crystal template.

Input:  Protein-only PDB (no ligand). May include water.
Output: Relaxed PDB with optimized sidechains and near-native backbone.

Uses ref2015 score function (same as docking stage).

Usage:
    python constrained_relax.py input.pdb output.pdb [--params PARAMS] [--sd 0.5]
"""
import argparse
import sys
from pathlib import Path


BACKBONE_ATOMS = ["N", "CA", "C", "O"]


def add_backbone_coordinate_constraints(pose, sd=0.5):
    """
    Add CoordinateConstraints on backbone heavy atoms (N, CA, C, O)
    for all protein residues.

    Requires that addVirtualResAsRoot(pose) has already been called.
    The VRT residue (last residue) provides the fixed reference atom.

    Args:
        pose: PyRosetta Pose with VRT as root
        sd: Standard deviation for HarmonicFunc in Angstroms
    """
    from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.id import AtomID

    vrt_res = pose.total_residue()
    vrt_atom_id = AtomID(1, vrt_res)  # Atom 1 ("ORIG") of VRT residue

    func = HarmonicFunc(0.0, sd)

    n_constrained = 0
    for i in range(1, vrt_res):  # All residues except VRT
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
    """
    Remove VRT residue lines from a PDB file in-place.

    PyRosetta's dump_pdb includes the virtual residue added by
    addVirtualResAsRoot. This strips those lines so downstream
    tools see a clean protein PDB.
    """
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
        description="Backbone-constrained FastRelax for threaded mutant PDBs"
    )
    parser.add_argument("input_pdb", help="Input PDB (threaded mutant)")
    parser.add_argument("output_pdb", help="Output relaxed PDB")
    parser.add_argument("--params", default=None,
                        help="Optional extra_res_fa params file")
    parser.add_argument("--sd", type=float, default=0.5,
                        help="Harmonic constraint SD in Angstroms (default: 0.5)")
    parser.add_argument("--cycles", type=int, default=5,
                        help="FastRelax cycles (default: 5)")
    args = parser.parse_args()

    # Ensure output directory exists
    Path(args.output_pdb).parent.mkdir(parents=True, exist_ok=True)

    # Build PyRosetta init flags
    init_flags = "-ex1 -ex2aro -use_input_sc -relax:fast -mute all"
    if args.params:
        init_flags = f"-extra_res_fa {args.params} {init_flags}"

    import pyrosetta
    from pyrosetta import init, pose_from_file
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.pose import addVirtualResAsRoot
    from pyrosetta.rosetta.core.scoring import coordinate_constraint

    init(init_flags)

    # Load pose
    print(f"Loading: {args.input_pdb}")
    pose = pose_from_file(args.input_pdb)
    print(f"Loaded {pose.total_residue()} residues")

    # Add virtual residue as fold-tree root (anchor for CoordinateConstraints)
    addVirtualResAsRoot(pose)
    print(f"Added VRT root (now {pose.total_residue()} residues)")

    # Add backbone coordinate constraints
    add_backbone_coordinate_constraints(pose, sd=args.sd)

    # Score function: ref2015 + coordinate_constraint weight
    sfxn = pyrosetta.create_score_function("ref2015")
    sfxn.set_weight(coordinate_constraint, 1.0)

    # Run FastRelax
    print(f"Running FastRelax ({args.cycles} cycles)...")
    relax = FastRelax(sfxn, args.cycles)
    relax.apply(pose)

    score_after = sfxn(pose)
    print(f"Post-relax score: {score_after:.2f}")

    # Save output
    pose.dump_pdb(args.output_pdb)
    strip_vrt_from_pdb(args.output_pdb)
    print(f"Saved relaxed structure: {args.output_pdb}")


if __name__ == "__main__":
    main()
