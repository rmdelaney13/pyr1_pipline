#!/usr/bin/env python
"""
Diagnostic script for debugging unrealistic docking scores.

Runs a SINGLE conformer through the same pipeline as grade_conformers_mutant_docking.py
but with detailed energy breakdowns at every stage.

Usage:
    python debug_docking_scoring.py <config_file> [--max-conformers N] [--dump-pdbs]

    Uses the same config file as grade_conformers_mutant_docking.py.
    Runs from the docking/scripts/ directory (same as production scripts).

Example:
    cd docking/scripts
    python debug_docking_scoring.py /path/to/docking_config.txt --max-conformers 2 --dump-pdbs
"""

import argparse
import os
import sys
import time
import logging
from configparser import ConfigParser

import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("debug_docking")

# Path setup (same as production script)
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
sys.path.insert(0, os.path.join(script_dir, ".."))
sys.path.insert(0, os.path.join(script_dir, "..", "legacy"))

import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.scoring import ScoreType
import pyrosetta.rosetta.protocols.rigid as rigid_moves
import pyrosetta.rosetta.core.scoring as scoring

import docking_pipeline_utils as dpu
import alignment
import conformer_prep
import collision_check

# Re-use config helpers from production script
from grade_conformers_mutant_docking import (
    _cfg_clean, _cfg_float, _cfg_int, _cfg_bool, _cfg_str,
    setup_packer_task, prepare_dataframe,
    add_hbond_constraint, add_hbond_constraint_to_water, auto_setup_water_constraints,
)


# ═══════════════════════════════════════════════════════════════════
# DIAGNOSTIC HELPERS
# ═══════════════════════════════════════════════════════════════════

SCORE_TERMS = [
    ScoreType.fa_rep,          # VdW repulsion (clashes)
    ScoreType.fa_atr,          # VdW attraction
    ScoreType.fa_elec,         # Electrostatics
    ScoreType.fa_sol,          # Solvation / desolvation
    ScoreType.hbond_sc,        # Sidechain H-bonds
    ScoreType.hbond_bb_sc,     # BB-SC H-bonds
    ScoreType.hbond_lr_bb,     # Long-range BB H-bonds
    ScoreType.hbond_sr_bb,     # Short-range BB H-bonds
    ScoreType.ref,             # Reference energies
    ScoreType.rama_prepro,     # Ramachandran
    ScoreType.omega,           # Omega angle
    ScoreType.fa_dun,          # Dunbrack rotamer
    ScoreType.p_aa_pp,         # Probability of AA given phi/psi
    ScoreType.atom_pair_constraint,  # Constraints
]


def score_breakdown(sf, pose, label=""):
    """Score a pose and print per-term energy breakdown."""
    total = sf(pose)
    print(f"\n{'=' * 65}")
    print(f"  SCORE BREAKDOWN: {label}")
    print(f"{'=' * 65}")

    terms = []
    for term in SCORE_TERMS:
        try:
            val = sf.score_by_scoretype(pose, term)
        except Exception:
            val = 0.0
        if abs(val) > 0.01:
            terms.append((term.name, val))

    terms.sort(key=lambda x: -abs(x[1]))
    for name, val in terms:
        bar = "#" * min(int(abs(val) / 10), 50)
        sign = "+" if val > 0 else ""
        print(f"  {name:<25s} {sign}{val:>12.2f}  {bar}")

    accounted = sum(v for _, v in terms)
    other = total - accounted
    if abs(other) > 0.1:
        print(f"  {'(other terms)':<25s} {'+' if other > 0 else ''}{other:>12.2f}")

    print(f"  {'─' * 40}")
    print(f"  {'TOTAL':<25s} {'+' if total > 0 else ''}{total:>12.2f}")
    print(f"{'=' * 65}")
    return total


def per_residue_energy(sf, pose, residue_idx, label=""):
    """Show energy contribution of a single residue to the total."""
    sf(pose)
    energies = pose.energies()
    total_res = 0.0
    print(f"\n  Per-term energy for residue {residue_idx} "
          f"({pose.residue(residue_idx).name3()}): {label}")

    terms = []
    for term in SCORE_TERMS:
        try:
            wt = sf.get_weight(term)
            if abs(wt) < 1e-6:
                continue
            raw = energies.residue_total_energies(residue_idx)[term]
            weighted = raw * wt
        except Exception:
            continue
        if abs(weighted) > 0.01:
            terms.append((term.name, weighted))
            total_res += weighted

    terms.sort(key=lambda x: -abs(x[1]))
    for name, val in terms:
        print(f"    {name:<25s} {'+' if val > 0 else ''}{val:>10.2f}")
    print(f"    {'RESIDUE TOTAL':<25s} {'+' if total_res > 0 else ''}{total_res:>10.2f}")
    return total_res


def water_analysis(pose, lig_idx):
    """Analyze all waters in the pose relative to the ligand."""
    lig = pose.residue(lig_idx)
    waters = []
    for ri in range(1, pose.total_residue() + 1):
        res = pose.residue(ri)
        if not res.is_water():
            continue
        o_xyz = None
        for ai in range(1, res.natoms() + 1):
            if res.atom_type(ai).element() == "O":
                o_xyz = np.array([res.xyz(ai).x, res.xyz(ai).y, res.xyz(ai).z])
                break
        if o_xyz is None:
            continue

        # Distance to each ligand heavy atom
        min_dist = None
        closest_atom = None
        for ai in range(1, lig.natoms() + 1):
            if lig.atom_is_hydrogen(ai):
                continue
            a_xyz = np.array([lig.xyz(ai).x, lig.xyz(ai).y, lig.xyz(ai).z])
            d = float(np.linalg.norm(o_xyz - a_xyz))
            if min_dist is None or d < min_dist:
                min_dist = d
                closest_atom = lig.atom_name(ai).strip()

        waters.append({
            "residue_idx": ri,
            "name3": res.name3(),
            "min_dist_to_lig": min_dist,
            "closest_lig_atom": closest_atom,
        })

    print(f"\n  WATER ANALYSIS ({len(waters)} waters in pose)")
    if not waters:
        print("  *** NO WATERS FOUND — H-bond geometry will always fail ***")
        return waters

    waters.sort(key=lambda w: w["min_dist_to_lig"])
    print(f"  {'Res#':<6s} {'Name':<6s} {'Dist(Å)':<10s} {'Closest Lig Atom':<20s} {'Status'}")
    print(f"  {'─' * 55}")
    for w in waters[:15]:  # Show top 15 closest
        d = w["min_dist_to_lig"]
        if d < 1.5:
            status = "*** CLASH ***"
        elif d < 2.5:
            status = "H-bond range"
        elif d < 4.0:
            status = "Near"
        else:
            status = ""
        print(f"  {w['residue_idx']:<6d} {w['name3']:<6s} {d:<10.3f} {w['closest_lig_atom']:<20s} {status}")

    clashing = sum(1 for w in waters if w["min_dist_to_lig"] < 1.5)
    hbond_range = sum(1 for w in waters if 1.5 <= w["min_dist_to_lig"] <= 3.5)
    if clashing > 0:
        print(f"\n  *** {clashing} water(s) CLASH with ligand (< 1.5 Å) — major score penalty ***")
    print(f"  Waters in H-bond range (1.5-3.5 Å): {hbond_range}")
    return waters


def ligand_acceptor_analysis(pose, lig_idx):
    """Find all acceptor atoms on the ligand and their closest water."""
    lig = pose.residue(lig_idx)
    print(f"\n  LIGAND ACCEPTOR ANALYSIS (residue {lig_idx}: {lig.name3()})")
    print(f"  {'Atom':<8s} {'Element':<8s} {'IsAcceptor':<12s} {'Closest Water':<15s} {'Dist(Å)':<10s}")
    print(f"  {'─' * 55}")

    for ai in range(1, lig.natoms() + 1):
        if lig.atom_is_hydrogen(ai):
            continue
        name = lig.atom_name(ai).strip()
        elem = lig.atom_type(ai).element()
        try:
            is_acc = lig.heavyatom_is_an_acceptor(ai)
        except Exception:
            is_acc = False

        if not is_acc and elem not in {"O", "N"}:
            continue

        a_xyz = np.array([lig.xyz(ai).x, lig.xyz(ai).y, lig.xyz(ai).z])
        best_dist = None
        best_water = None
        for ri in range(1, pose.total_residue() + 1):
            w = pose.residue(ri)
            if not w.is_water():
                continue
            for oi in range(1, w.natoms() + 1):
                if w.atom_type(oi).element() != "O":
                    continue
                o_xyz = np.array([w.xyz(oi).x, w.xyz(oi).y, w.xyz(oi).z])
                d = float(np.linalg.norm(o_xyz - a_xyz))
                if best_dist is None or d < best_dist:
                    best_dist = d
                    best_water = ri

        print(f"  {name:<8s} {elem:<8s} {'YES' if is_acc else 'no':<12s} "
              f"{'res ' + str(best_water) if best_water else 'NONE':<15s} "
              f"{best_dist:.3f}" if best_dist else "N/A")


def collision_analysis(mutant_pose, conf, lig_idx):
    """Test collision detection with different settings."""
    water_indices = set(
        i for i in range(1, mutant_pose.total_residue() + 1)
        if mutant_pose.residue(i).is_water()
    )

    configs = [
        ("BB-only, vdw=0.7 (production)", False, 0.7, water_indices),
        ("BB-only, vdw=0.5 (lenient)",     False, 0.5, water_indices),
        ("BB+SC, vdw=0.7 (strict)",        True,  0.7, water_indices),
        ("BB+SC, vdw=0.5 (moderate)",      True,  0.5, water_indices),
        ("BB+SC+water, vdw=0.7 (everything)", True, 0.7, set()),
    ]

    print(f"\n  COLLISION ANALYSIS")
    print(f"  {'Configuration':<40s} {'Collides?':<12s}")
    print(f"  {'─' * 55}")

    for label, inc_sc, vdw_mod, excluded in configs:
        try:
            grid = collision_check.CollisionGrid(
                mutant_pose,
                bin_width=1.0,
                vdw_modifier=vdw_mod,
                include_sc=inc_sc,
                excluded_residues=excluded,
            )
            collides = conf.check_collision(grid)
        except Exception as e:
            collides = f"ERROR: {e}"
        print(f"  {label:<40s} {str(collides):<12s}")


def score_with_without_waters(sf, pose, lig_idx, label=""):
    """Compare scores with and without explicit waters."""
    score_with = sf(pose)

    no_water = pose.clone()
    water_res = [i for i in range(1, no_water.total_residue() + 1)
                 if no_water.residue(i).is_water()]
    # Track ligand index shift
    lig_shift = sum(1 for w in water_res if w < lig_idx)
    for wr in reversed(water_res):
        no_water.delete_residue_slow(wr)
    score_without = sf(no_water)
    water_contribution = score_with - score_without

    # Score protein-only (no ligand, no water)
    no_lig = no_water.clone()
    adjusted_lig = lig_idx - lig_shift
    if adjusted_lig <= no_lig.total_residue():
        no_lig.delete_residue_slow(adjusted_lig)
    score_protein_only_no_water = sf(no_lig)

    print(f"\n  WATER IMPACT ANALYSIS: {label}")
    print(f"  {'Score with waters + ligand:':<40s} {score_with:>12.2f}")
    print(f"  {'Score WITHOUT waters + ligand:':<40s} {score_without:>12.2f}")
    print(f"  {'Water contribution:':<40s} {water_contribution:>12.2f}")
    print(f"  {'Protein only (no lig, no water):':<40s} {score_protein_only_no_water:>12.2f}")
    print(f"  {'Ligand contribution (no water):':<40s} {score_without - score_protein_only_no_water:>12.2f}")

    return {
        "with_water": score_with,
        "without_water": score_without,
        "water_contribution": water_contribution,
        "protein_only_no_water": score_protein_only_no_water,
    }


def score_ligand_only(sf, pose, lig_idx):
    """Score the ligand residue in isolation."""
    lig_only = Pose()
    lig_only.append_residue_by_jump(pose.residue(lig_idx), 0)
    score = sf(lig_only)
    print(f"\n  LIGAND-ONLY SCORE: {score:.2f} REU")
    if score > 100:
        print(f"  *** WARNING: Ligand internal energy is very high ({score:.0f} REU) ***")
        print(f"  This suggests bad atom types or geometry in the params file.")
    elif score > 50:
        print(f"  NOTE: Moderate ligand strain ({score:.0f} REU). May need better conformers.")
    else:
        print(f"  Ligand internal energy looks reasonable.")
    return score


# ═══════════════════════════════════════════════════════════════════
# MAIN DEBUG LOOP
# ═══════════════════════════════════════════════════════════════════

def debug_docking(config_path, max_conformers=2, dump_pdbs=False):
    """Run diagnostic docking on a few conformers with full score breakdowns."""

    config = ConfigParser()
    config.read(config_path)

    if "mutant_docking" not in config:
        logger.error("Config file must have [mutant_docking] section")
        sys.exit(1)

    section = config["mutant_docking"]
    def_section = config["DEFAULT"]

    mutant_pdb = _cfg_str(section, "MutantPDB", "")
    csv_file = _cfg_str(def_section, "CSVFileName", "")
    path_to_conformers = _cfg_str(def_section, "PathToConformers", "")
    chain_letter = _cfg_str(def_section, "ChainLetter", "A")
    residue_number = _cfg_int(def_section, "ResidueNumber", 1)
    lig_res_num = _cfg_int(def_section, "LigandResidueNumber", 1)
    auto_align = _cfg_bool(def_section, "AutoGenerateAlignment", False)
    pre_pdb = _cfg_str(def_section, "PrePDBFileName", "")

    # Resolve relative paths: try CWD first, then project root (two levels up from
    # docking/scripts/), then relative to the config file directory.
    project_root = os.path.abspath(os.path.join(script_dir, "..", ".."))
    config_dir = os.path.dirname(os.path.abspath(config_path))
    search_bases = [os.getcwd(), project_root, config_dir]

    def resolve_path(p):
        if not p or os.path.isabs(p):
            return p
        for base in search_bases:
            candidate = os.path.join(base, p)
            if os.path.exists(candidate):
                return os.path.abspath(candidate)
        return p  # return as-is if not found anywhere

    mutant_pdb = resolve_path(mutant_pdb)
    csv_file = resolve_path(csv_file)
    path_to_conformers = resolve_path(path_to_conformers)
    pre_pdb = resolve_path(pre_pdb)
    # Resolve params list paths too (done below after parsing)

    hbond_ideal = _cfg_float(section, "HBondDistanceIdeal", 2.8)
    dist_buf = _cfg_float(section, "HBondDistanceIdealBuffer", 0.8)
    hbond_min = max(0.0, hbond_ideal - dist_buf)
    hbond_max = hbond_ideal + dist_buf
    max_pass_score = _cfg_float(section, "MaxScore", 100.0)

    params_list = [resolve_path(p) for p in _cfg_str(def_section, "ParamsList", "").split()]

    # ── Print config summary ──
    print("\n" + "█" * 65)
    print("  DOCKING SCORE DEBUGGER")
    print("█" * 65)
    print(f"  Config:         {config_path}")
    print(f"  Mutant PDB:     {mutant_pdb}")
    print(f"  Reference PDB:  {pre_pdb}")
    print(f"  Conformers CSV: {csv_file}")
    print(f"  Conformers dir: {path_to_conformers}")
    print(f"  Params:         {params_list}")
    print(f"  MaxScore:       {max_pass_score}")
    print(f"  HBond ideal:    {hbond_ideal} Å (window: {hbond_min:.1f}-{hbond_max:.1f})")
    print("█" * 65)

    # ── Validate files exist ──
    for label, path in [("Mutant PDB", mutant_pdb), ("CSV", csv_file)]:
        if not path or not os.path.exists(path):
            logger.error(f"{label} not found: {path}")
            sys.exit(1)

    # ── Init PyRosetta ──
    if params_list:
        pyrosetta.init(f"-extra_res_fa {' '.join(params_list)} -mute all")
    else:
        pyrosetta.init("-mute all")

    sf = pyrosetta.get_fa_scorefxn()

    # ── Load mutant pose ──
    mutant_pose = Pose()
    if params_list:
        res_set = pyrosetta.generate_nonstandard_residue_set(mutant_pose, params_list)
        pyrosetta.pose_from_file(mutant_pose, res_set, mutant_pdb)
    else:
        pyrosetta.pose_from_file(mutant_pose, mutant_pdb)

    n_res = mutant_pose.total_residue()
    n_water = sum(1 for i in range(1, n_res + 1) if mutant_pose.residue(i).is_water())
    n_protein = n_res - n_water
    print(f"\n  Mutant pose: {n_res} residues ({n_protein} protein, {n_water} water)")

    # ── DIAGNOSTIC 1: Baseline Score ──
    print("\n" + "▓" * 65)
    print("  DIAGNOSTIC 1: BASELINE SCORE (protein + waters, NO ligand)")
    print("▓" * 65)
    baseline = score_breakdown(sf, mutant_pose, "Mutant pose (baseline)")
    per_res = baseline / n_protein if n_protein > 0 else 0
    print(f"\n  Score per protein residue: {per_res:.2f} REU/res")
    if per_res > 0:
        print("  *** WARNING: Positive per-residue score suggests the mutant PDB")
        print("  *** has internal clashes and may need relaxation before docking.")
    elif per_res > -2:
        print("  NOTE: Marginal per-residue score. Structure may benefit from relaxation.")
    else:
        print("  Baseline quality looks reasonable.")

    # Show water contribution to baseline
    if n_water > 0:
        no_water_baseline = mutant_pose.clone()
        wrs = [i for i in range(1, no_water_baseline.total_residue() + 1)
               if no_water_baseline.residue(i).is_water()]
        for wr in reversed(wrs):
            no_water_baseline.delete_residue_slow(wr)
        base_no_water = sf(no_water_baseline)
        print(f"\n  Baseline without waters: {base_no_water:.2f}")
        print(f"  Water contribution to baseline: {baseline - base_no_water:.2f}")

    # ── Load reference pose and get target residue ──
    target_res = None
    target_res_idx = None
    pocket_center = None

    if pre_pdb and os.path.exists(pre_pdb):
        reference_pose = Pose()
        if params_list:
            res_set_ref = pyrosetta.generate_nonstandard_residue_set(reference_pose, params_list)
            pyrosetta.pose_from_file(reference_pose, res_set_ref, pre_pdb)
        else:
            pyrosetta.pose_from_file(reference_pose, pre_pdb)

        target_idx = reference_pose.pdb_info().pdb2pose(chain_letter, residue_number)
        if target_idx > 0:
            target_res = reference_pose.residue(target_idx)
            target_res_idx = target_idx
            template_coords = dpu.ligand_heavy_atom_coords(reference_pose, target_idx)
            pocket_center = dpu.ligand_com(template_coords)
            print(f"\n  Reference target: {chain_letter}{residue_number} -> "
                  f"pose idx {target_idx} ({target_res.name3()})")
            print(f"  Pocket center: [{pocket_center[0]:.2f}, {pocket_center[1]:.2f}, {pocket_center[2]:.2f}]")
        else:
            logger.error(f"Cannot resolve target residue {chain_letter}{residue_number} in reference")
            sys.exit(1)
    else:
        logger.error(f"Reference PDB not found: {pre_pdb}")
        sys.exit(1)

    # ── Load conformer table ──
    df = prepare_dataframe(csv_file, path_to_conformers, target_res, lig_res_num, auto_align)
    print(f"\n  Conformer table: {len(df)} rows")

    if df.empty:
        logger.error("No conformers in table")
        sys.exit(1)

    # ── Setup movers (same as production) ──
    packer = setup_packer_task()

    sf_cst = pyrosetta.get_fa_scorefxn()
    sf_cst.set_weight(scoring.atom_pair_constraint,
                      _cfg_float(section, "HBondConstraintWeight", 4.0))

    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mm = pyrosetta.MoveMap()
    mm.set_jump(False)
    min_mover.movemap(mm)
    min_mover.score_function(sf_cst)

    rotation = _cfg_float(section, "Rotation", 25.0)
    translation = _cfg_float(section, "Translation", 0.5)
    max_tries = _cfg_int(section, "MaxPerturbTries", 30)

    # Build collision grid (same as production)
    water_indices = set(
        i for i in range(1, mutant_pose.total_residue() + 1)
        if mutant_pose.residue(i).is_water()
    )
    backbone_grid = collision_check.CollisionGrid(
        mutant_pose,
        bin_width=_cfg_float(section, "BinWidth", 1.0),
        vdw_modifier=_cfg_float(section, "VDW_Modifier", 0.7),
        include_sc=False,
        excluded_residues=water_indices,
    )

    # ── DIAGNOSTIC 2+: Process conformers ──
    conf_count = 0
    for conf_idx, conf in enumerate(
        conformer_prep.yield_ligand_poses(df, path_to_conformers, False, lig_res_num),
        start=1,
    ):
        if not conf:
            continue
        if conf_count >= max_conformers:
            break
        conf_count += 1

        print("\n\n" + "▓" * 65)
        print(f"  DIAGNOSTIC: CONFORMER {conf_idx} (conf_num={getattr(conf, 'conf_num', 'NA')})")
        print("▓" * 65)

        # Get alignment atoms
        molecule_atoms = None
        target_atoms = None
        if hasattr(conf, "molecule_atoms") and hasattr(conf, "target_atoms"):
            molecule_atoms = conf.molecule_atoms
            target_atoms = conf.target_atoms
        elif hasattr(conf, "lig_atoms") and hasattr(conf, "target_atoms"):
            molecule_atoms = conf.lig_atoms
            target_atoms = conf.target_atoms
        elif hasattr(conf, "lig_aid") and hasattr(conf, "t_aid"):
            molecule_atoms = conf.lig_aid
            target_atoms = conf.t_aid

        if molecule_atoms is None or target_atoms is None:
            print("  *** SKIP: Missing alignment atoms ***")
            continue

        print(f"  Molecule atoms: {molecule_atoms}")
        print(f"  Target atoms:   {target_atoms}")

        # Align
        try:
            conf.align_to_target(target_res)
        except Exception as exc:
            print(f"  *** ALIGNMENT FAILED: {exc} ***")
            continue

        # Graft ligand
        try:
            grafted = mutant_pose.clone()
            anchor = len(grafted.chain_sequence(1))
            grafted.append_residue_by_jump(conf.pose.residue(1), anchor)
            lig_idx = grafted.total_residue()

            # Find ligand jump
            lig_jump_num = None
            ft = grafted.fold_tree()
            for jj in range(1, ft.num_jump() + 1):
                if ft.downstream_jump_residue(jj) == lig_idx:
                    lig_jump_num = jj
                    break
        except Exception as exc:
            print(f"  *** GRAFTING FAILED: {exc} ***")
            continue

        print(f"  Ligand grafted at residue {lig_idx}, jump {lig_jump_num}")
        print(f"  Ligand name: {grafted.residue(lig_idx).name3()}, "
              f"atoms: {grafted.residue(lig_idx).natoms()}")

        # ── Score right after grafting (no perturbation) ──
        print("\n" + "-" * 65)
        print("  STAGE: Immediately after grafting (no perturbation)")
        print("-" * 65)
        score_breakdown(sf, grafted, "Post-graft (no perturbation)")
        per_residue_energy(sf, grafted, lig_idx, "Ligand residue")
        score_ligand_only(sf, grafted, lig_idx)
        water_analysis(grafted, lig_idx)
        ligand_acceptor_analysis(grafted, lig_idx)

        # ── Try perturbations ──
        print("\n" + "-" * 65)
        print(f"  STAGE: Perturbation loop (rotation={rotation}°, translation={translation} Å)")
        print("-" * 65)

        best_pose = None
        best_quality = -1
        best_try = None

        for try_idx in range(1, min(max_tries, 10) + 1):
            copy_pose = grafted.clone()
            rigid_moves.RigidBodyPerturbMover(
                lig_jump_num, rotation, translation
            ).apply(copy_pose)

            # Sync coords for collision check
            for atom_id in range(1, copy_pose.residue(lig_idx).natoms() + 1):
                conf.pose.residue(1).set_xyz(
                    atom_id, copy_pose.residue(lig_idx).xyz(atom_id)
                )

            collides = conf.check_collision(backbone_grid)
            if collides:
                if try_idx <= 3:
                    print(f"  Try {try_idx}: backbone collision → skip")
                continue

            # Closest acceptor
            acceptor_name = molecule_atoms[0]
            neighbor_names = list(molecule_atoms[1:3])
            from grade_conformers_mutant_docking import _closest_ligand_acceptor_to_any_water
            closest_atom, _, _ = _closest_ligand_acceptor_to_any_water(copy_pose, lig_idx)
            if closest_atom is not None:
                acceptor_name = closest_atom
            neighbor_names = []

            # Add constraints
            add_hbond_constraint_to_water(
                copy_pose, lig_idx, acceptor_name,
                dist_ideal=hbond_ideal,
                dist_sd=_cfg_float(section, "HBondConstraintSD", 0.25),
                capture_max=_cfg_float(section, "HBondConstraintCaptureMax", 8.0),
            )
            auto_setup_water_constraints(copy_pose, lig_idx)

            # Minimize
            mm.set_jump(False)
            if lig_jump_num:
                mm.set_jump(lig_jump_num, True)
            min_mover.apply(copy_pose)

            # Evaluate geometry
            hbond_result = dpu.evaluate_hbond_geometry(
                copy_pose, lig_idx, acceptor_name, neighbor_names,
                hbond_min, hbond_max, hbond_ideal,
                _cfg_float(section, "HBondDonorAngleMin", 120.0),
                _cfg_float(section, "HBondAcceptorAngleMin", 90.0),
                0.0,
            )

            quality = hbond_result.get("quality", 0.0) or 0.0
            dist = hbond_result.get("distance")
            donor_a = hbond_result.get("donor_angle")

            status = "PASS" if hbond_result.get("passed") else "fail"
            dist_str = f"{dist:.2f}" if dist else "None"
            donor_str = f"{donor_a:.1f}" if donor_a else "None"

            # Pre-pack score
            pre_pack_score = sf(copy_pose)
            relative = pre_pack_score - baseline

            print(f"  Try {try_idx}: {status} | quality={quality:.3f} | "
                  f"dist={dist_str} | donor={donor_str} | "
                  f"score={pre_pack_score:.0f} (rel={relative:.0f})")

            if quality > best_quality:
                best_quality = quality
                best_pose = copy_pose.clone()
                best_try = try_idx

        if best_pose is None:
            print("\n  *** ALL TRIES COLLIDED — no pose to analyze ***")
            print("  This means ligand cannot fit in the pocket even with backbone-only collision.")
            collision_analysis(mutant_pose, conf, lig_idx)
            continue

        # ── Detailed analysis of best pre-pack pose ──
        print(f"\n  Best pre-pack pose: try {best_try}, quality={best_quality:.3f}")

        print("\n" + "-" * 65)
        print("  STAGE: Best pre-pack pose (before sidechain repacking)")
        print("-" * 65)
        pre_pack_total = score_breakdown(sf, best_pose, "Pre-pack best pose")
        per_residue_energy(sf, best_pose, lig_idx, "Ligand (pre-pack)")
        water_analysis(best_pose, lig_idx)

        # ── Pack sidechains ──
        print("\n" + "-" * 65)
        print("  STAGE: After sidechain repacking")
        print("-" * 65)
        packed_pose = best_pose.clone()
        packer.apply(packed_pose)
        post_pack_total = score_breakdown(sf, packed_pose, "Post-pack")
        per_residue_energy(sf, packed_pose, lig_idx, "Ligand (post-pack)")

        relative_score = post_pack_total - baseline
        print(f"\n  ┌─────────────────────────────────────────────┐")
        print(f"  │  FINAL RELATIVE SCORE: {relative_score:>10.2f} REU        │")
        print(f"  │  Threshold (MaxScore): {max_pass_score:>10.2f} REU        │")
        print(f"  │  Verdict: {'PASS ✓' if relative_score < max_pass_score else 'FAIL ✗':>10s}                       │")
        print(f"  └─────────────────────────────────────────────┘")

        # ── Water impact ──
        print("\n" + "-" * 65)
        print("  STAGE: Water impact analysis")
        print("-" * 65)
        score_with_without_waters(sf, packed_pose, lig_idx, "Post-pack")

        # ── Collision with sidechains ──
        print("\n" + "-" * 65)
        print("  STAGE: Collision grid analysis (post-pack)")
        print("-" * 65)
        # Re-sync conformer coords to packed pose
        for atom_id in range(1, packed_pose.residue(lig_idx).natoms() + 1):
            conf.pose.residue(1).set_xyz(
                atom_id, packed_pose.residue(lig_idx).xyz(atom_id)
            )
        collision_analysis(mutant_pose, conf, lig_idx)

        # ── Post-pack H-bond geometry ──
        print("\n" + "-" * 65)
        print("  STAGE: Post-pack H-bond geometry")
        print("-" * 65)
        postpack_hbond = dpu.evaluate_hbond_geometry(
            packed_pose, lig_idx,
            acceptor_name, [],
            hbond_min, hbond_max, hbond_ideal,
            _cfg_float(section, "HBondDonorAngleMin", 120.0),
            _cfg_float(section, "HBondAcceptorAngleMin", 90.0),
            0.0,
        )
        for k, v in sorted(postpack_hbond.items()):
            print(f"  {k}: {v}")

        # ── Dump PDBs ──
        if dump_pdbs:
            debug_dir = os.path.join(os.path.dirname(config_path), "debug_pdbs")
            os.makedirs(debug_dir, exist_ok=True)
            pre_path = os.path.join(debug_dir, f"debug_conf{conf_idx}_prepack.pdb")
            post_path = os.path.join(debug_dir, f"debug_conf{conf_idx}_postpack.pdb")
            graft_path = os.path.join(debug_dir, f"debug_conf{conf_idx}_graft_nopert.pdb")
            best_pose.dump_pdb(pre_path)
            packed_pose.dump_pdb(post_path)
            grafted.dump_pdb(graft_path)
            print(f"\n  Dumped PDBs:")
            print(f"    Pre-perturbation graft: {graft_path}")
            print(f"    Pre-pack (best try):    {pre_path}")
            print(f"    Post-pack:              {post_path}")

    # ── Summary ──
    print("\n\n" + "█" * 65)
    print("  SUMMARY")
    print("█" * 65)
    print(f"  Conformers analyzed: {conf_count}")
    print(f"  Baseline score: {baseline:.2f} REU ({per_res:.2f}/residue)")
    print(f"  MaxScore threshold: {max_pass_score:.2f} REU (relative)")
    print(f"\n  If all relative scores are in the thousands:")
    print(f"    1. Check SCORE BREAKDOWN — is fa_rep the dominant term?")
    print(f"       → YES: ligand has severe VdW clashes with protein sidechains")
    print(f"       → Consider: include_sc=True in collision grid, or relax mutant first")
    print(f"    2. Check WATER IMPACT — does removing water drop the score significantly?")
    print(f"       → YES: water placement is causing scoring artifacts")
    print(f"       → Consider: check water naming (TP3/HOH), check water overlap with ligand")
    print(f"    3. Check LIGAND-ONLY SCORE — is it > 100 REU?")
    print(f"       → YES: ligand params file has issues")
    print(f"       → Regenerate with molfile_to_params.py")
    print(f"    4. Check COLLISION ANALYSIS — does BB+SC collision detect clashes?")
    print(f"       → YES but production only checks BB: sidechains are the clash source")
    print(f"       → Consider: reducing perturbation magnitude or enabling SC collision")
    print("█" * 65)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Debug docking scoring issues",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("config", help="Docking config file (same as grade_conformers_mutant_docking.py)")
    parser.add_argument("--max-conformers", type=int, default=2,
                        help="Number of conformers to analyze (default: 2)")
    parser.add_argument("--dump-pdbs", action="store_true",
                        help="Save debug PDBs for visualization in PyMOL")
    args = parser.parse_args()

    debug_docking(args.config, args.max_conformers, args.dump_pdbs)
