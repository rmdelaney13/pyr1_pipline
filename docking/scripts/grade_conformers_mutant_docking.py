#!/usr/bin/env python
"""
Conformer docking to PRE-THREADED MUTANT PYR1 structures with water-constrained alignment.

This script is designed for ML dataset generation where we need to evaluate
how specific PYR1 variants bind ligands. Unlike the glycine-shaved workflow
(used for design), this script docks directly to the mutant pocket with real
sidechains, using water-mediated H-bond constraints for high-quality poses.

Workflow:
1. Load pre-threaded mutant PDB (from thread_variant_to_pdb.py)
2. Load conformer CSV table with alignment atoms (from create_table.py)
3. Align conformers to pocket using SVD
4. Perturb + minimize with water H-bond constraints
5. Filter by H-bond geometry (pre-pack and post-pack)
6. Score and cluster by RMSD
7. Output geometry diagnostics CSV

Usage:
  Single job:   python grade_conformers_mutant_docking.py config.txt
  Array task:   python grade_conformers_mutant_docking.py config.txt [array_index]
                (or set SLURM_ARRAY_TASK_ID environment variable)

Config file format:
  [mutant_docking]
  MutantPDB = /path/to/mutant_59K_120A_160G.pdb
  OutputDir = /path/to/output/
  EnableHBondGeometryFilter = True
  HBondDistanceIdeal = 2.8
  MaxPerturbTries = 30
  ArrayTaskCount = 10  # For SLURM array parallelization

Author: Claude Code (Whitehead Lab PYR1 Pipeline)
Date: 2026-02-16
"""

import argparse
import os
import sys
import time
import logging
from configparser import ConfigParser
from typing import List, Dict, Any, Optional

import pandas as pd
import numpy as np

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# Add script directory and legacy to path for utility imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
sys.path.insert(0, os.path.join(script_dir, ".."))
sys.path.insert(0, os.path.join(script_dir, "..", "legacy"))

# Import PyRosetta (will be initialized in main)
import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.scoring import ScoreFunction
import pyrosetta.rosetta.protocols.grafting as graft
import pyrosetta.rosetta.protocols.rigid as rigid_moves
import pyrosetta.rosetta.core.scoring as scoring

# Import pipeline utilities
import docking_pipeline_utils as dpu
import alignment
import conformer_prep


# ═══════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

def _cfg_clean(raw):
    """Remove inline comments from config values."""
    if raw is None:
        return ""
    return str(raw).split("#", 1)[0].strip()


def _cfg_float(section, key, default):
    """Parse float from config section."""
    raw = section.get(key, None)
    cleaned = _cfg_clean(raw)
    if cleaned == "":
        cleaned = str(default)
    return float(cleaned)


def _cfg_int(section, key, default):
    """Parse int from config section."""
    raw = section.get(key, None)
    cleaned = _cfg_clean(raw)
    if cleaned == "":
        cleaned = str(default)
    return int(float(cleaned))


def _cfg_bool(section, key, default):
    """Parse bool from config section."""
    raw = section.get(key, None)
    cleaned = _cfg_clean(raw).lower()
    if cleaned == "":
        return bool(default)
    return cleaned in {"1", "true", "yes", "on"}


def _cfg_str(section, key, default):
    """Parse string from config section."""
    raw = section.get(key, None)
    cleaned = _cfg_clean(raw)
    return cleaned if cleaned else str(default)


def _resolve_array_runtime(cli_array_index, section):
    """Determine array task index and count from CLI or environment."""
    env_idx = os.environ.get("SLURM_ARRAY_TASK_ID", "").strip()
    env_count = os.environ.get("SLURM_ARRAY_TASK_COUNT", "").strip()

    if cli_array_index is not None:
        array_index = int(cli_array_index)
    elif env_idx:
        array_index = int(env_idx)
    else:
        array_index = 0

    if env_count:
        array_count = max(1, int(env_count))
    else:
        array_count = max(1, _cfg_int(section, "ArrayTaskCount", 1))

    if array_index < 0:
        raise ValueError(f"Array index must be >= 0, got {array_index}")
    if array_index >= array_count:
        raise ValueError(
            f"Array index {array_index} is out of range for array count {array_count}"
        )
    return array_index, array_count


def _slice_dataframe_for_array(df, array_index, array_count):
    """Slice DataFrame for SLURM array task."""
    if array_count <= 1:
        return df.reset_index(drop=True)
    keep = (np.arange(len(df)) % array_count) == array_index
    return df.loc[keep].reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════
# ROSETTA SETUP HELPERS
# ═══════════════════════════════════════════════════════════════════

def setup_packer_task():
    """Configure packer task factory for sidechain repacking."""
    from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
    from pyrosetta.rosetta.protocols import minimization_packing as pack_min

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.IncludeCurrent())
    tf.push_back(operation.NoRepackDisulfides())
    tf.push_back(operation.RestrictToRepacking())
    packer = pack_min.PackRotamersMover()
    packer.task_factory(tf)
    return packer


# ═══════════════════════════════════════════════════════════════════
# CONSTRAINT HELPERS
# ═══════════════════════════════════════════════════════════════════

def add_hbond_constraint(pose, lig_res_idx, lig_atom_name, target_res_idx,
                        target_atom_name, dist_ideal=2.8, dist_sd=0.2):
    """Add direct harmonic constraint between ligand and target atoms."""
    from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.id import AtomID

    try:
        if not pose.residue(lig_res_idx).has(lig_atom_name):
            return
        if not pose.residue(target_res_idx).has(target_atom_name):
            return

        lig_atom_id = AtomID(
            pose.residue(lig_res_idx).atom_index(lig_atom_name),
            lig_res_idx
        )
        target_atom_id = AtomID(
            pose.residue(target_res_idx).atom_index(target_atom_name),
            target_res_idx
        )

        func = HarmonicFunc(dist_ideal, dist_sd)
        cst = AtomPairConstraint(lig_atom_id, target_atom_id, func)
        pose.add_constraint(cst)
    except Exception:
        pass


def add_hbond_constraint_to_water(
    pose, lig_res_idx, lig_atom_name, dist_ideal=2.8,
    dist_sd=0.3, capture_max=8.0
):
    """
    Add dynamic pull between ligand acceptor and nearest water oxygen.
    Helps minimization achieve ideal H-bonding distance.
    """
    from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.id import AtomID

    try:
        lig = pose.residue(lig_res_idx)
        if not lig.has(lig_atom_name):
            return
        lig_atom_idx = lig.atom_index(lig_atom_name)
        lig_xyz = lig.xyz(lig_atom_idx)

        # Find nearest water oxygen
        best = None
        best_dist = None
        for resid in range(1, pose.total_residue() + 1):
            water = pose.residue(resid)
            if not water.is_water():
                continue
            for o_idx in range(1, water.natoms() + 1):
                if water.atom_type(o_idx).element() != "O":
                    continue
                d = lig_xyz.distance(water.xyz(o_idx))
                if best_dist is None or d < best_dist:
                    best_dist = d
                    best = (resid, o_idx)

        if best is None or best_dist is None or best_dist > float(capture_max):
            return

        wat_res_idx, wat_o_idx = best
        lig_atom_id = AtomID(lig_atom_idx, lig_res_idx)
        wat_atom_id = AtomID(wat_o_idx, wat_res_idx)
        func = HarmonicFunc(float(dist_ideal), float(dist_sd))
        pose.add_constraint(AtomPairConstraint(lig_atom_id, wat_atom_id, func))
    except Exception:
        pass


def auto_setup_water_constraints(pose, lig_res_idx, dist_cutoff=3.5):
    """
    Automatically add constraints between nearby waters and potential H-bond partners.
    Stabilizes water network during minimization.
    """
    from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.id import AtomID

    lig_res = pose.residue(lig_res_idx)
    water_indices = [
        i for i in range(1, pose.total_residue() + 1)
        if pose.residue(i).is_water() and
        pose.residue(i).nbr_atom_xyz().distance(lig_res.nbr_atom_xyz()) < 6.0
    ]

    for wat_res_idx in water_indices:
        wat_res = pose.residue(wat_res_idx)
        potential_acceptors, potential_donors = [], []

        # Scan nearby residues for H-bond partners
        for i in range(1, pose.total_residue() + 1):
            if i == wat_res_idx:
                continue
            curr_res = pose.residue(i)
            if curr_res.nbr_atom_xyz().distance(wat_res.nbr_atom_xyz()) > (dist_cutoff + 5.0):
                continue

            for atom_i in range(1, curr_res.natoms() + 1):
                dist = curr_res.xyz(atom_i).distance(wat_res.xyz("O"))
                if dist > dist_cutoff:
                    continue

                atom_name = curr_res.atom_name(atom_i).strip()
                is_lig = (i == lig_res_idx)
                # Prioritize ligand atoms slightly
                priority_offset = -0.5 if is_lig else 0.0

                if curr_res.atom_is_polar_hydrogen(atom_i):
                    potential_donors.append((dist + priority_offset, i, atom_name))
                elif curr_res.heavyatom_is_an_acceptor(atom_i):
                    potential_acceptors.append((dist + priority_offset, i, atom_name))

        potential_donors.sort()
        potential_acceptors.sort()

        # Add constraint to best donor (H -> water O)
        if potential_donors:
            best = potential_donors[0]
            pose.add_constraint(
                AtomPairConstraint(
                    AtomID(wat_res.atom_index("O"), wat_res_idx),
                    AtomID(pose.residue(best[1]).atom_index(best[2]), best[1]),
                    HarmonicFunc(2.0, 0.2)
                )
            )

        # Add constraint to best acceptor (water H -> acceptor)
        if potential_acceptors:
            best = potential_acceptors[0]
            pose.add_constraint(
                AtomPairConstraint(
                    AtomID(wat_res.atom_index("H1"), wat_res_idx),
                    AtomID(pose.residue(best[1]).atom_index(best[2]), best[1]),
                    HarmonicFunc(2.0, 0.2)
                )
            )


# ═══════════════════════════════════════════════════════════════════
# GEOMETRY VALIDATION HELPERS
# ═══════════════════════════════════════════════════════════════════

def _passes_hbond_ideal_window(hbond_result, params):
    """
    Check if H-bond geometry falls within ideal windows (pre-minimization).
    Optional stricter gate controlled by params['enforce_ideal_window'].
    """
    if not params.get("enforce_ideal_window", False):
        return True

    dist = hbond_result.get("distance", None)
    donor = hbond_result.get("donor_angle", None)
    acceptor = hbond_result.get("acceptor_angle", None)

    if dist is None or donor is None:
        return False

    if abs(dist - params["hbond_ideal"]) > params["distance_ideal_buffer"]:
        return False

    donor_floor = params["donor_ideal"] - params["donor_ideal_buffer"]
    if donor < donor_floor:
        return False

    if params.get("require_acceptor_ideal_window", False):
        if acceptor is None:
            return False
        acceptor_floor = params["acceptor_ideal"] - params["acceptor_ideal_buffer"]
        if acceptor < acceptor_floor:
            return False

    return True


def _passes_hbond_ideal_window_explicit(hbond_result, params):
    """
    Always enforce ideal geometry checks (used for final post-pack filtering).
    This is the CRITICAL filter for ML dataset quality.
    """
    dist = hbond_result.get("distance", None)
    donor = hbond_result.get("donor_angle", None)
    acceptor = hbond_result.get("acceptor_angle", None)

    if dist is None or donor is None:
        return False

    if abs(dist - params["hbond_ideal"]) > params["distance_ideal_buffer"]:
        return False

    donor_floor = params["donor_ideal"] - params["donor_ideal_buffer"]
    if donor < donor_floor:
        return False

    if params.get("require_acceptor_ideal_window", False):
        if acceptor is None:
            return False
        acceptor_floor = params["acceptor_ideal"] - params["acceptor_ideal_buffer"]
        if acceptor < acceptor_floor:
            return False

    return True


def _closest_ligand_acceptor_to_any_water(pose, lig_res_idx):
    """
    Find ligand acceptor atom nearest to any water oxygen.
    Returns (atom_name, min_distance, water_residue_idx).
    """
    import numpy as np

    lig = pose.residue(lig_res_idx)
    acceptor_indices = []

    for i in range(1, lig.natoms() + 1):
        try:
            is_acceptor = lig.heavyatom_is_an_acceptor(i)
        except Exception:
            is_acceptor = False
        if not is_acceptor:
            continue
        if lig.atom_type(i).element() not in {"O", "N"}:
            continue
        acceptor_indices.append(i)

    if not acceptor_indices:
        return (None, None, None)

    best_atom = None
    best_dist = None
    best_water = None

    for acc_idx in acceptor_indices:
        acc_xyz = np.array([
            float(lig.xyz(acc_idx).x),
            float(lig.xyz(acc_idx).y),
            float(lig.xyz(acc_idx).z)
        ], dtype=float)

        for resid in range(1, pose.total_residue() + 1):
            water = pose.residue(resid)
            if not water.is_water():
                continue
            for o_idx in range(1, water.natoms() + 1):
                if water.atom_type(o_idx).element() != "O":
                    continue
                o_xyz = np.array([
                    float(water.xyz(o_idx).x),
                    float(water.xyz(o_idx).y),
                    float(water.xyz(o_idx).z)
                ], dtype=float)
                dist = float(np.linalg.norm(o_xyz - acc_xyz))
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best_atom = lig.atom_name(acc_idx).strip()
                    best_water = resid

    return (best_atom, best_dist, best_water)


def _water_diagnostics(pose, lig_res_idx, acceptor_atom_name):
    """
    Return water diagnostics for debugging H-bond failures.
    Returns (num_waters, nearest_water_O_distance_to_acceptor).
    """
    import numpy as np

    lig = pose.residue(lig_res_idx)
    acc_idx = dpu.atom_index_by_name(lig, acceptor_atom_name) if acceptor_atom_name else None

    if acc_idx is None:
        # Fallback: find first O or N atom
        for i in range(1, lig.natoms() + 1):
            if lig.atom_type(i).element() in {"O", "N"}:
                acc_idx = i
                break

    if acc_idx is None:
        return (0, None)

    acc_xyz = dpu._xyz_to_np(lig.xyz(acc_idx))
    num_waters = 0
    nearest = None

    for resid in range(1, pose.total_residue() + 1):
        water = pose.residue(resid)
        if not water.is_water():
            continue
        num_waters += 1
        for o_idx in range(1, water.natoms() + 1):
            if water.atom_type(o_idx).element() != "O":
                continue
            o_xyz = dpu._xyz_to_np(water.xyz(o_idx))
            dist = float(np.linalg.norm(o_xyz - acc_xyz))
            if nearest is None or dist < nearest:
                nearest = dist

    return (num_waters, nearest)


# ═══════════════════════════════════════════════════════════════════
# CLUSTERING HELPERS
# ═══════════════════════════════════════════════════════════════════

def _find_existing_cluster(coords, clusters, rmsd_cutoff):
    """Check if coords match any existing cluster within RMSD cutoff."""
    for idx, cluster in enumerate(clusters):
        if dpu.ligand_rmsd(coords, cluster["coords"]) <= rmsd_cutoff:
            return idx
    return None


# ═══════════════════════════════════════════════════════════════════
# TABLE PREPARATION (CSV LOADING)
# ═══════════════════════════════════════════════════════════════════

def prepare_dataframe(csv_file, path_to_conformers, target_res, lig_res_num, auto_align):
    """
    Load and prepare conformer table from CSV.
    Adapted from glycine-shaved script for mutant docking use case.

    Args:
        csv_file: Path to CSV table (from create_table.py)
        path_to_conformers: Directory containing conformer params/pdb files
        target_res: Rosetta Residue object for alignment target
        lig_res_num: Ligand residue number in pose (usually 1)
        auto_align: Whether to auto-generate alignment atoms (usually False)

    Returns:
        Prepared DataFrame with alignment atom columns
    """
    df = pd.read_csv(csv_file, index_col=False)
    df["Molecule File Stem"] = df["Molecule ID"].apply(lambda name: f"{name}/{name}")
    df["Conformer Range"] = df["Conformer Range"].apply(
        lambda v: tuple(str(v).split("_")[:2]) if "_" in str(v) else (0, 0)
    )

    if auto_align:
        for i, row in df.iterrows():
            if pd.notnull(row.get("Molecule Atoms")):
                continue
            mol_id = row["Molecule ID"]
            params = f"{path_to_conformers}/{mol_id}/{mol_id}.params"
            pdb = f"{path_to_conformers}/{mol_id}/{mol_id}_0001.pdb"
            if not os.path.exists(params):
                continue

            lig = Pose()
            pyrosetta.pose_from_file(
                lig,
                pyrosetta.generate_nonstandard_residue_set(lig, [params]),
                pdb
            )
            mol_ats, tar_ats = alignment.auto_align_residue_to_residue(
                lig, lig.residue(lig_res_num), target_res
            )
            df.at[i, "Molecule Atoms"] = "-".join(mol_ats)
            df.at[i, "Target Atoms"] = "-".join(tar_ats)

    df["Molecule Atoms"] = df["Molecule Atoms"].apply(
        lambda x: tuple(str(x).split("-")) if pd.notnull(x) else None
    )
    df["Target Atoms"] = df["Target Atoms"].apply(
        lambda x: tuple(str(x).split("-")) if pd.notnull(x) else None
    )

    return df.dropna(subset=["Molecule Atoms", "Target Atoms"]).copy()


# ═══════════════════════════════════════════════════════════════════
# MAIN DOCKING FUNCTION
# ═══════════════════════════════════════════════════════════════════

def align_and_dock_conformers(
    mutant_pose,
    target_res,
    target_res_idx,
    path_to_conformers,
    df,
    params,
    output_dir=None,
    rename_water_to_tp3=True,
    cluster_enabled=True,
    cluster_rmsd_cutoff=0.75,
    geometry_csv_path=None,
    output_name_prefix="",
):
    """
    Main docking loop: align conformers to mutant pocket and evaluate H-bonds.

    Args:
        mutant_pose: Pre-threaded mutant PYR1 structure
        target_res: Target residue for alignment (from reference structure)
        target_res_idx: Target residue index in mutant_pose (may be None if glycine-shaved)
        path_to_conformers: Directory with conformer params/pdb files
        df: DataFrame with conformer metadata (from prepare_dataframe)
        params: Dictionary of docking parameters (from config)
        output_dir: Output directory for docked PDBs
        rename_water_to_tp3: Whether to rename waters to TP3 for LigandMPNN
        cluster_enabled: Whether to perform RMSD clustering in-loop
        cluster_rmsd_cutoff: RMSD cutoff for clustering (Angstroms)
        geometry_csv_path: Path to write geometry diagnostics CSV
        output_name_prefix: Prefix for output PDB files

    Returns:
        Dictionary with statistics (total, passed_score, clustered)
    """
    # Initialize Rosetta movers and score functions
    packer = setup_packer_task()
    sf_all = pyrosetta.get_fa_scorefxn()
    sf_cst = pyrosetta.get_fa_scorefxn()
    sf_cst.set_weight(
        scoring.atom_pair_constraint,
        float(params.get("hbond_constraint_weight", 1.0))
    )

    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mm = pyrosetta.MoveMap()
    mm.set_jump(True)
    min_mover.movemap(mm)
    min_mover.score_function(sf_cst)

    # Statistics tracking
    stats = {"total": 0, "passed_score": 0, "clustered": 0}
    clusters = []
    geometry_rows = []

    max_tries = int(params.get("max_tries", 30))
    debug_every = int(params.get("debug_every", 10))

    # Main conformer loop
    for conf_idx, conf in enumerate(
        conformer_prep.yield_ligand_poses(df, path_to_conformers, False, params['lig_res_num']),
        start=1
    ):
        if not conf:
            continue

        if conf_idx == 1 or conf_idx % debug_every == 0:
            logger.info(
                "Processing conformer %d (conf_num=%s). Stats: seen=%d passed=%d clustered=%d",
                conf_idx, getattr(conf, "conf_num", "NA"),
                stats["total"], stats["passed_score"], stats["clustered"]
            )

        # Extract alignment atoms (support multiple attribute naming conventions)
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
        else:
            logger.error(f"Conformer missing alignment attributes. Found: {dir(conf)}")
            continue

        if molecule_atoms is None or target_atoms is None or \
           len(molecule_atoms) < 1 or len(target_atoms) < 1:
            logger.error("Conformer has empty alignment atom lists; skipping.")
            continue

        # SVD alignment to target residue
        try:
            conf.align_to_target(target_res)
        except Exception as exc:
            logger.error(f"Alignment failed for conformer {conf_idx}: {exc}")
            continue

        # Graft ligand into mutant pose
        try:
            grafted_pose = graft.insert_pose_into_pose(
                mutant_pose,
                conf.pose,
                len(mutant_pose.chain_sequence(1))
            )
            lig_idx = len(mutant_pose.chain_sequence(1)) + 1
        except Exception as exc:
            logger.error(f"Grafting failed for conformer {conf_idx}: {exc}")
            continue

        # Perturbation loop
        keep_candidate = False
        accepted_try_idx = None
        last_hbond_result = None
        strict_ok = False
        out_path = None
        score = None

        for try_idx in range(1, max_tries + 1):
            copy_pose = grafted_pose.clone()

            # Apply random perturbation
            rigid_moves.RigidBodyPerturbMover(
                params['jump_num'],
                params['rotation'],
                params['translation']
            ).apply(copy_pose)

            # Determine acceptor atom for H-bond constraint
            acceptor_name = molecule_atoms[0]
            neighbor_names = list(molecule_atoms[1:3])

            if params.get("use_closest_acceptor", True):
                closest_atom, _, _ = _closest_ligand_acceptor_to_any_water(copy_pose, lig_idx)
                if closest_atom is not None:
                    acceptor_name = closest_atom
                # No robust generic neighbor definition for dynamic atom
                neighbor_names = []

            # Add constraints
            if (not params.get("use_closest_acceptor", True)) and \
               target_res_idx is not None and target_res_idx > 0:
                # Fixed constraint to target residue atom
                add_hbond_constraint(
                    copy_pose, lig_idx, acceptor_name,
                    target_res_idx, target_atoms[0], params['hbond_ideal']
                )

            if params.get("use_closest_acceptor", True) and \
               params.get("dynamic_water_pull_constraint", True):
                # Dynamic constraint to nearest water
                add_hbond_constraint_to_water(
                    copy_pose, lig_idx, acceptor_name,
                    dist_ideal=params['hbond_ideal'],
                    dist_sd=params.get('hbond_constraint_sd', 0.3),
                    capture_max=params.get('hbond_constraint_capture_max', 8.0),
                )

            # Auto-constrain water network
            auto_setup_water_constraints(copy_pose, lig_idx)

            # Minimize with constraints
            t_min0 = time.time()
            try:
                min_mover.apply(copy_pose)
            except Exception as exc:
                logger.warning(f"Minimization failed for conformer {conf_idx}, try {try_idx}: {exc}")
                continue
            min_dt = time.time() - t_min0

            if min_dt > float(params.get("slow_min_threshold_sec", 5.0)):
                logger.info(
                    "Slow minimization: conformer %d (conf_num=%s), try %d/%d took %.2fs",
                    conf_idx, getattr(conf, "conf_num", "NA"), try_idx, max_tries, min_dt
                )

            # Evaluate H-bond geometry (pre-pack)
            try:
                hbond_result = dpu.evaluate_hbond_geometry(
                    copy_pose, lig_idx, acceptor_name, neighbor_names,
                    params['hbond_min'], params['hbond_max'], params['hbond_ideal'],
                    params['donor_angle'], params['acceptor_angle'], params['quality_min']
                )
            except Exception as exc:
                logger.warning(f"H-bond evaluation failed for conformer {conf_idx}, try {try_idx}: {exc}")
                continue

            last_hbond_result = hbond_result
            strict_ok = _passes_hbond_ideal_window(hbond_result, params)

            # Check if we accept this perturbation
            if ((hbond_result["passed"] and strict_ok) or
                (not params['use_hbond_filter'])):
                keep_candidate = True
                accepted_try_idx = try_idx
                break

            # Log failure on last try
            if try_idx == max_tries:
                num_waters, nearest_water_dist = _water_diagnostics(
                    copy_pose, lig_idx, molecule_atoms[0]
                )
                logger.info(
                    "Conformer %d (conf_num=%s) failed hbond filter after %d tries; "
                    "best quality=%.3f dist=%s donor=%s acceptor=%s water_res=%s "
                    "waters=%d nearest_O=%.3f strict=%s",
                    conf_idx, getattr(conf, "conf_num", "NA"), max_tries,
                    float(hbond_result.get("quality", 0.0)),
                    f"{hbond_result.get('distance'):.3f}" if hbond_result.get("distance") else "None",
                    f"{hbond_result.get('donor_angle'):.1f}" if hbond_result.get("donor_angle") else "None",
                    f"{hbond_result.get('acceptor_angle'):.1f}" if hbond_result.get("acceptor_angle") else "None",
                    hbond_result.get("water_residue", None),
                    num_waters,
                    nearest_water_dist if nearest_water_dist is not None else -1.0,
                    strict_ok,
                )

        # Record failure if no perturbation worked
        if not keep_candidate:
            geometry_rows.append({
                "conf_idx": conf_idx,
                "conf_num": getattr(conf, "conf_num", "NA"),
                "accepted": False,
                "accepted_try": accepted_try_idx,
                "passed_score": False,
                "saved_cluster": False,
                "score": None,
                "distance": None if last_hbond_result is None else last_hbond_result.get("distance"),
                "donor_angle": None if last_hbond_result is None else last_hbond_result.get("donor_angle"),
                "acceptor_angle": None if last_hbond_result is None else last_hbond_result.get("acceptor_angle"),
                "quality": None if last_hbond_result is None else last_hbond_result.get("quality"),
                "water_residue": None if last_hbond_result is None else last_hbond_result.get("water_residue"),
                "strict_window_passed": strict_ok,
                "output_pdb": None,
            })
            stats["total"] += 1
            continue

        # Pack sidechains
        try:
            packer.apply(copy_pose)
            score = sf_all(copy_pose)
        except Exception as exc:
            logger.warning(f"Packing failed for conformer {conf_idx}: {exc}")
            stats["total"] += 1
            continue

        # Post-pack H-bond validation (CRITICAL for ML dataset quality)
        try:
            postpack_hbond_result = dpu.evaluate_hbond_geometry(
                copy_pose, lig_idx, acceptor_name, neighbor_names,
                params['hbond_min'], params['hbond_max'], params['hbond_ideal'],
                params['donor_angle'], params['acceptor_angle'], params['quality_min'],
            )
        except Exception as exc:
            logger.warning(f"Post-pack H-bond evaluation failed for conformer {conf_idx}: {exc}")
            stats["total"] += 1
            continue

        saved_cluster = False
        postpack_pass = bool(postpack_hbond_result.get("passed", False))
        postpack_strict_ok = (
            _passes_hbond_ideal_window_explicit(postpack_hbond_result, params)
            if params.get("enforce_final_ideal_geometry", True)
            else True
        )

        # Final acceptance criteria
        if (not params['use_hbond_filter']) or (postpack_pass and postpack_strict_ok):
            if score < params['max_pass_score']:
                stats["passed_score"] += 1
                coords = dpu.ligand_heavy_atom_coords(copy_pose, lig_idx)

                # Optional RMSD clustering
                c_idx = None
                if cluster_enabled:
                    c_idx = _find_existing_cluster(coords, clusters, cluster_rmsd_cutoff)

                if c_idx is None:
                    # Save unique pose
                    out_dir = output_dir or "output"
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir, exist_ok=True)

                    out = os.path.join(out_dir, f"{output_name_prefix}rep_{stats['total']}.pdb")
                    dpu.dump_pose_pdb(copy_pose, out, rename_water=rename_water_to_tp3)

                    clusters.append({"coords": coords, "score": score, "path": out})
                    stats["clustered"] += 1
                    out_path = out
                    saved_cluster = True
        else:
            # Log post-pack geometry rejection
            logger.info(
                "Post-pack geometry reject: conformer %d (conf_num=%s) "
                "dist=%s donor=%s acceptor=%s passed=%s strict=%s",
                conf_idx, getattr(conf, "conf_num", "NA"),
                f"{postpack_hbond_result.get('distance'):.3f}" if postpack_hbond_result.get("distance") else "None",
                f"{postpack_hbond_result.get('donor_angle'):.1f}" if postpack_hbond_result.get("donor_angle") else "None",
                f"{postpack_hbond_result.get('acceptor_angle'):.1f}" if postpack_hbond_result.get("acceptor_angle") else "None",
                postpack_pass,
                postpack_strict_ok,
            )

        # Record geometry results
        geometry_rows.append({
            "conf_idx": conf_idx,
            "conf_num": getattr(conf, "conf_num", "NA"),
            "accepted": True,
            "accepted_try": accepted_try_idx,
            "passed_score": bool(score is not None and score < params['max_pass_score']),
            "saved_cluster": saved_cluster,
            "score": score,
            "distance": postpack_hbond_result.get("distance"),
            "donor_angle": postpack_hbond_result.get("donor_angle"),
            "acceptor_angle": postpack_hbond_result.get("acceptor_angle"),
            "quality": postpack_hbond_result.get("quality"),
            "water_residue": postpack_hbond_result.get("water_residue"),
            "strict_window_passed": _passes_hbond_ideal_window(postpack_hbond_result, params),
            "postpack_geometry_passed": postpack_pass,
            "postpack_ideal_passed": postpack_strict_ok,
            "output_pdb": out_path,
        })
        stats["total"] += 1

    # Write geometry CSV
    if geometry_csv_path:
        pd.DataFrame(geometry_rows).to_csv(geometry_csv_path, index=False)
        logger.info("Wrote H-bond geometry CSV: %s (%d rows)", geometry_csv_path, len(geometry_rows))

    return stats


# ═══════════════════════════════════════════════════════════════════
# MAIN FUNCTION
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Dock ligand conformers to pre-threaded mutant PYR1 (with water constraints)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('config', help='Configuration file (INI format)')
    parser.add_argument('array_index', nargs='?', type=int, help='SLURM array task index (optional)')

    args = parser.parse_args()

    # Load config
    config = ConfigParser()
    config.read(args.config)

    if 'mutant_docking' not in config:
        logger.error("Config file must have [mutant_docking] section")
        sys.exit(1)

    section = config['mutant_docking']
    def_section = config['DEFAULT']

    # Resolve array parameters
    array_index, array_count = _resolve_array_runtime(args.array_index, section)

    # Parse config parameters
    mutant_pdb = _cfg_str(section, 'MutantPDB', '')
    csv_file = _cfg_str(def_section, 'CSVFileName', '')
    path_to_conformers = _cfg_str(def_section, 'PathToConformers', '')
    output_dir = _cfg_str(section, 'OutputDir', './output')

    # Alignment parameters
    chain_letter = _cfg_str(def_section, 'ChainLetter', 'A')
    residue_number = _cfg_int(def_section, 'ResidueNumber', 1)
    lig_res_num = _cfg_int(def_section, 'LigandResidueNumber', 1)
    auto_align = _cfg_bool(def_section, 'AutoGenerateAlignment', False)

    # H-bond distance parameters (with smart defaults)
    hbond_ideal = _cfg_float(section, "HBondDistanceIdeal", 2.8)
    dist_buf = _cfg_float(section, "HBondDistanceIdealBuffer", 0.8)

    if "HBondDistanceMin" in section:
        hbond_min = _cfg_float(section, "HBondDistanceMin", hbond_ideal - dist_buf)
    else:
        hbond_min = max(0.0, hbond_ideal - dist_buf)

    if "HBondDistanceMax" in section:
        hbond_max = _cfg_float(section, "HBondDistanceMax", hbond_ideal + dist_buf)
    else:
        hbond_max = hbond_ideal + dist_buf

    # Quality minimum (optional)
    raw_quality_min = _cfg_clean(section.get("HBondQualityMin", ""))
    quality_min = float(raw_quality_min) if raw_quality_min else 0.0

    # Build run parameters dictionary
    run_params = {
        'lig_res_num': lig_res_num,
        'jump_num': _cfg_int(section, "JumpNum", 2),
        'rotation': _cfg_float(section, "Rotation", 25.0),
        'translation': _cfg_float(section, "Translation", 0.5),
        'max_pass_score': _cfg_float(section, "MaxScore", -300.0),
        'hbond_ideal': hbond_ideal,
        'hbond_min': hbond_min,
        'hbond_max': hbond_max,
        'donor_angle': _cfg_float(section, "HBondDonorAngleMin", 120.0),
        'acceptor_angle': _cfg_float(section, "HBondAcceptorAngleMin", 90.0),
        'quality_min': quality_min,
        'use_hbond_filter': _cfg_bool(section, "EnableHBondGeometryFilter", True),
        'use_closest_acceptor': _cfg_bool(section, "UseClosestLigandAcceptor", True),
        'hbond_constraint_weight': _cfg_float(section, "HBondConstraintWeight", 4.0),
        'hbond_constraint_sd': _cfg_float(section, "HBondConstraintSD", 0.25),
        'hbond_constraint_capture_max': _cfg_float(section, "HBondConstraintCaptureMax", 8.0),
        'dynamic_water_pull_constraint': _cfg_bool(section, "EnableDynamicWaterPullConstraint", True),
        'enforce_ideal_window': _cfg_bool(section, "EnforceHBondIdealWindow", False),
        'enforce_final_ideal_geometry': _cfg_bool(section, "EnforceFinalIdealGeometry", True),
        'distance_ideal_buffer': dist_buf,
        'donor_ideal': _cfg_float(section, "HBondDonorAngleIdeal", 180.0),
        'donor_ideal_buffer': _cfg_float(section, "HBondDonorAngleBuffer", 40.0),
        'acceptor_ideal': _cfg_float(section, "HBondAcceptorAngleIdeal", 180.0),
        'acceptor_ideal_buffer': _cfg_float(section, "HBondAcceptorAngleBuffer", 30.0),
        'require_acceptor_ideal_window': _cfg_bool(section, "RequireAcceptorIdealWindow", False),
        'postprocess_rename_water_sweep': _cfg_bool(section, "PostprocessRenameWaterSweep", True),
        'cluster_enabled': _cfg_bool(section, "EnablePoseClusteringInArrayTask", False),
        'cluster_rmsd_cutoff': _cfg_float(section, "ClusterRMSDCutoff", 0.75),
        'max_tries': _cfg_int(section, "MaxPerturbTries", 30),
        'debug_every': _cfg_int(section, "DebugEveryN", 10),
        'slow_min_threshold_sec': _cfg_float(section, "SlowMinimizationSeconds", 5.0)
    }

    # Validate inputs
    if not os.path.exists(mutant_pdb):
        logger.error(f"Mutant PDB not found: {mutant_pdb}")
        sys.exit(1)

    if not os.path.exists(csv_file):
        logger.error(f"CSV table not found: {csv_file}")
        sys.exit(1)

    # Initialize PyRosetta
    params_list = _cfg_str(def_section, 'ParamsList', '').split()
    if params_list:
        pyrosetta.init(f"-extra_res_fa {' '.join(params_list)} -mute all")
    else:
        pyrosetta.init("-mute all")

    # Load mutant structure
    mutant_pose = Pose()
    if params_list:
        res_set = pyrosetta.generate_nonstandard_residue_set(mutant_pose, params_list)
        pyrosetta.pose_from_file(mutant_pose, res_set, mutant_pdb)
    else:
        pyrosetta.pose_from_file(mutant_pose, mutant_pdb)

    logger.info(f"Loaded mutant pose: {mutant_pose.total_residue()} residues")

    # Resolve target residue for alignment atom definitions
    target_idx = mutant_pose.pdb_info().pdb2pose(chain_letter, residue_number)
    if target_idx <= 0 or target_idx > mutant_pose.total_residue():
        logger.error(
            f"Could not map target residue {chain_letter}{residue_number} into mutant pose. "
            f"pdb2pose={target_idx}, total_residue={mutant_pose.total_residue()}"
        )
        sys.exit(1)

    target_res = mutant_pose.residue(target_idx)
    logger.info(f"Target residue: {chain_letter}{residue_number} -> pose index {target_idx}")

    # Load conformer table
    df = prepare_dataframe(csv_file, path_to_conformers, target_res, lig_res_num, auto_align)

    # Slice for array task
    if array_count > 1:
        df = _slice_dataframe_for_array(df, array_index, array_count)
        logger.info(
            f"Array task {array_index}/{array_count}: assigned {len(df)} rows"
        )

    if df.empty:
        logger.info("No conformer rows assigned to this array task; exiting.")
        return

    # Output setup
    os.makedirs(output_dir, exist_ok=True)

    base_geometry_name = _cfg_str(section, "GeometryCSVName", "hbond_geometry_summary.csv")
    if array_count > 1:
        stem, ext = os.path.splitext(base_geometry_name)
        base_geometry_name = f"{stem}_array{array_index:04d}{ext or '.csv'}"
    geometry_csv_path = os.path.join(output_dir, base_geometry_name)

    output_name_prefix = _cfg_str(section, "OutputPrefix", "")
    if not output_name_prefix and array_count > 1:
        output_name_prefix = f"a{array_index:04d}_"

    # Log configuration
    logger.info("=" * 60)
    logger.info("MUTANT DOCKING WORKFLOW")
    logger.info("=" * 60)
    logger.info(f"Mutant PDB: {mutant_pdb}")
    logger.info(f"CSV table: {csv_file}")
    logger.info(f"Conformers: {path_to_conformers}")
    logger.info(f"Output dir: {output_dir}")
    logger.info(f"Array task: {array_index}/{array_count}")
    logger.info(f"Rows assigned: {len(df)}")
    logger.info("=" * 60)
    logger.info(
        "HBond filter: enabled=%s closest_acceptor=%s cst_wt=%.2f "
        "min=%.3f max=%.3f ideal=%.3f donor_min=%.1f acceptor_min=%.1f "
        "quality_min=%.3f max_tries=%d",
        run_params['use_hbond_filter'],
        run_params['use_closest_acceptor'],
        run_params['hbond_constraint_weight'],
        run_params['hbond_min'],
        run_params['hbond_max'],
        run_params['hbond_ideal'],
        run_params['donor_angle'],
        run_params['acceptor_angle'],
        run_params['quality_min'],
        run_params['max_tries'],
    )
    if run_params['enforce_final_ideal_geometry']:
        logger.info(
            "Final ideal geometry enforcement: dist=%.2f±%.2f donor>=%.1f%s",
            run_params['hbond_ideal'],
            run_params['distance_ideal_buffer'],
            run_params['donor_ideal'] - run_params['donor_ideal_buffer'],
            (
                f" acceptor>={run_params['acceptor_ideal'] - run_params['acceptor_ideal_buffer']:.1f}"
                if run_params['require_acceptor_ideal_window']
                else ""
            ),
        )
    logger.info("=" * 60)

    # Run docking
    start_time = time.time()

    stats = align_and_dock_conformers(
        mutant_pose,
        target_res,
        target_idx,
        path_to_conformers,
        df,
        run_params,
        output_dir=output_dir,
        cluster_enabled=run_params['cluster_enabled'],
        cluster_rmsd_cutoff=run_params['cluster_rmsd_cutoff'],
        geometry_csv_path=geometry_csv_path,
        output_name_prefix=output_name_prefix,
    )

    elapsed = time.time() - start_time

    # Post-process: rename waters to TP3
    if run_params['postprocess_rename_water_sweep']:
        count = 0
        for name in os.listdir(output_dir):
            if not name.lower().endswith(".pdb"):
                continue
            path = os.path.join(output_dir, name)
            try:
                dpu.rename_water_resnames_to_tp3(path)
                count += 1
            except Exception as exc:
                logger.warning(f"Failed TP3 rename for {path}: {exc}")
        logger.info(f"Renamed waters to TP3 in {count} PDB files")

    # Summary
    logger.info("=" * 60)
    logger.info(f"Docking complete in {elapsed:.1f}s")
    logger.info(f"Total conformers: {stats['total']}")
    logger.info(f"Passed score: {stats['passed_score']}")
    logger.info(f"Clustered (saved): {stats['clustered']}")
    logger.info("=" * 60)

    if not run_params['cluster_enabled'] and array_count > 1:
        logger.info(
            "Note: In-array clustering disabled. "
            "After all array tasks complete, run cluster_docked_post_array.py for global clustering."
        )


if __name__ == '__main__':
    main()
