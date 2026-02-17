#!/usr/bin/env python
"""
Sidechain repack around threaded mutation sites.

Repacks sidechains within a 10 A shell around mutated residues using
PackRotamersMover (same approach as the proven docking script
grade_conformers_docked_to_sequence_multiple_slurm1.py). No FastRelax,
no backbone movement, no coordinate constraints — just discrete rotamer
optimization. Completes in seconds.

Input:  Protein-only PDB (no ligand). May include water.
Output: Repacked PDB with optimized sidechains near mutation sites.

Uses ref2015 score function (same as docking stage).

Usage:
    python constrained_relax.py input.pdb output.pdb --mutations "59K;120A;160G"
"""
import argparse
import re
import sys
from pathlib import Path


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

    index_str = ",".join(str(i) for i in mutation_pose_indices)
    mut_sel = ResidueIndexSelector(index_str)
    shell_sel = NeighborhoodResidueSelector(mut_sel, shell_distance, True)

    shell_vec = shell_sel.apply(pose)
    shell_set = set()
    for i in range(1, pose.total_residue() + 1):
        if shell_vec[i]:
            shell_set.add(i)

    return shell_set


def main():
    parser = argparse.ArgumentParser(
        description="Sidechain repack around mutation sites for threaded mutant PDBs"
    )
    parser.add_argument("input_pdb", help="Input PDB (threaded mutant)")
    parser.add_argument("output_pdb", help="Output repacked PDB")
    parser.add_argument("--mutations", required=True,
                        help="Variant signature, e.g. '59K;120A;160G'")
    parser.add_argument("--params", default=None,
                        help="Optional extra_res_fa params file")
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
    init_flags = "-ex1 -ex2aro -use_input_sc -mute all"
    if args.params:
        init_flags = f"-extra_res_fa {args.params} {init_flags}"

    import pyrosetta
    from pyrosetta import init, pose_from_file
    from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
    from pyrosetta.rosetta.core.pack.task.operation import (
        RestrictToRepacking, IncludeCurrent, PreventRepacking
    )
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

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

    # Score function
    sfxn = pyrosetta.create_score_function("ref2015")

    score_before = sfxn(pose)
    print(f"Score before repack: {score_before:.2f}")

    # Build TaskFactory: repack shell residues only, prevent everything else
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(IncludeCurrent())
    tf.push_back(RestrictToRepacking())

    # Prevent repacking outside the shell
    prevent = PreventRepacking()
    for i in range(1, pose.total_residue() + 1):
        if i not in shell_residues:
            prevent.include_residue(i)
    tf.push_back(prevent)

    # Pack sidechains
    print(f"Repacking sidechains ({len(shell_residues)} residues)...")
    packer = PackRotamersMover(sfxn)
    packer.task_factory(tf)
    packer.apply(pose)

    score_after = sfxn(pose)
    print(f"Score after repack:  {score_after:.2f} (delta: {score_after - score_before:+.2f})")

    # Save output
    pose.dump_pdb(args.output_pdb)
    print(f"Saved repacked structure: {args.output_pdb}")


if __name__ == "__main__":
    main()
