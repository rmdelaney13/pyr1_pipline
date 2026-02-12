#!/usr/bin/env python
"""
Glycine-shaved conformer grading with SLURM array support.

This script processes ligand conformers with glycine-shaved pocket mutations,
optionally running as part of a SLURM job array. When run as an array:
  - Each task processes a slice of the conformer table
  - In-loop clustering is typically disabled (EnablePoseClusteringInArrayTask=False)
  - After all array tasks complete, run cluster_docked_post_array.py to
    aggregate outputs and perform global clustering

Usage:
  Single job:   python grade_conformers_glycine_shaved.py config.txt
  Array task:   python grade_conformers_glycine_shaved.py config.txt [array_index]
                (or set SLURM_ARRAY_TASK_ID environment variable)
  Post-cluster: python cluster_docked_post_array.py config.txt
"""
import argparse
import json
import os
import shutil
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

# --- Helper Functions ---

def _cfg_clean(raw):
    if raw is None:
        return ""
    # Support inline comments in config values (e.g., "0.5  # note")
    return str(raw).split("#", 1)[0].strip()

def _cfg_float(section, key, default):
    raw = section.get(key, None)
    cleaned = _cfg_clean(raw)
    if cleaned == "":
        cleaned = str(default)
    return float(cleaned)

def _cfg_int(section, key, default):
    raw = section.get(key, None)
    cleaned = _cfg_clean(raw)
    if cleaned == "":
        cleaned = str(default)
    return int(float(cleaned))

def _cfg_bool(section, key, default):
    raw = section.get(key, None)
    cleaned = _cfg_clean(raw).lower()
    if cleaned == "":
        return bool(default)
    return cleaned in {"1", "true", "yes", "on"}

def _resolve_array_runtime(cli_array_index, section):
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
    if array_count <= 1:
        return df.reset_index(drop=True)
    keep = (np.arange(len(df)) % array_count) == array_index
    return df.loc[keep].reset_index(drop=True)

def setup_packer_task():
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

def add_hbond_constraint(pose, lig_res_idx, lig_atom_name, target_res_idx, target_atom_name, dist_ideal=2.8, dist_sd=0.2):
    from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.id import AtomID

    try:
        if not pose.residue(lig_res_idx).has(lig_atom_name):
            return
        if not pose.residue(target_res_idx).has(target_atom_name):
            return
            
        lig_atom_id = AtomID(pose.residue(lig_res_idx).atom_index(lig_atom_name), lig_res_idx)
        target_atom_id = AtomID(pose.residue(target_res_idx).atom_index(target_atom_name), target_res_idx)
        
        func = HarmonicFunc(dist_ideal, dist_sd)
        cst = AtomPairConstraint(lig_atom_id, target_atom_id, func)
        pose.add_constraint(cst)
    except Exception:
        pass

def add_hbond_constraint_to_water(
    pose,
    lig_res_idx,
    lig_atom_name,
    dist_ideal=2.8,
    dist_sd=0.3,
    capture_max=8.0,
):
    """
    Add a direct pull between ligand acceptor and nearest water oxygen.
    This helps minimization move the ligand toward ideal H-bonding distance.
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
    from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.id import AtomID

    lig_res = pose.residue(lig_res_idx)
    water_indices = [i for i in range(1, pose.total_residue() + 1) 
                     if pose.residue(i).is_water() and 
                     pose.residue(i).nbr_atom_xyz().distance(lig_res.nbr_atom_xyz()) < 6.0]

    for wat_res_idx in water_indices:
        wat_res = pose.residue(wat_res_idx)
        potential_acceptors, potential_donors = [], []

        for i in range(1, pose.total_residue() + 1):
            if i == wat_res_idx: continue
            curr_res = pose.residue(i)
            if curr_res.nbr_atom_xyz().distance(wat_res.nbr_atom_xyz()) > (dist_cutoff + 5.0): continue

            for atom_i in range(1, curr_res.natoms() + 1):
                dist = curr_res.xyz(atom_i).distance(wat_res.xyz("O"))
                if dist > dist_cutoff: continue
                
                atom_name = curr_res.atom_name(atom_i).strip()
                is_lig = (i == lig_res_idx)
                if curr_res.atom_is_polar_hydrogen(atom_i):
                    potential_donors.append((dist - (0.5 if is_lig else 0.0), i, atom_name))
                elif curr_res.heavyatom_is_an_acceptor(atom_i):
                    potential_acceptors.append((dist - (0.5 if is_lig else 0.0), i, atom_name))

        potential_donors.sort()
        potential_acceptors.sort()

        if potential_donors:
            best = potential_donors[0]
            pose.add_constraint(AtomPairConstraint(AtomID(wat_res.atom_index("O"), wat_res_idx), 
                                                 AtomID(pose.residue(best[1]).atom_index(best[2]), best[1]), 
                                                 HarmonicFunc(2.0, 0.2)))
        if potential_acceptors:
            best = potential_acceptors[0]
            pose.add_constraint(AtomPairConstraint(AtomID(wat_res.atom_index("H1"), wat_res_idx), 
                                                 AtomID(pose.residue(best[1]).atom_index(best[2]), best[1]), 
                                                 HarmonicFunc(2.0, 0.2)))

def prepare_dataframe(csv_file, path_to_conformers, target_res, lig_res_num, auto_align):
    import alignment, pyrosetta
    from pyrosetta.rosetta.core.pose import Pose

    df = pd.read_csv(csv_file, index_col=False)
    df["Molecule File Stem"] = df["Molecule ID"].apply(lambda name: f"{name}/{name}")
    df["Conformer Range"] = df["Conformer Range"].apply(lambda v: tuple(str(v).split("_")[:2]) if "_" in str(v) else (0,0))

    if auto_align:
        for i, row in df.iterrows():
            if pd.notnull(row.get("Molecule Atoms")): continue
            mol_id = row["Molecule ID"]
            params = f"{path_to_conformers}/{mol_id}/{mol_id}.params"
            pdb = f"{path_to_conformers}/{mol_id}/{mol_id}_0001.pdb"
            if not os.path.exists(params): continue
            lig = Pose()
            pyrosetta.pose_from_file(lig, pyrosetta.generate_nonstandard_residue_set(lig, [params]), pdb)
            mol_ats, tar_ats = alignment.auto_align_residue_to_residue(lig, lig.residue(lig_res_num), target_res)
            df.at[i, "Molecule Atoms"], df.at[i, "Target Atoms"] = "-".join(mol_ats), "-".join(tar_ats)

    df["Molecule Atoms"] = df["Molecule Atoms"].apply(lambda x: tuple(str(x).split("-")) if pd.notnull(x) else None)
    df["Target Atoms"] = df["Target Atoms"].apply(lambda x: tuple(str(x).split("-")) if pd.notnull(x) else None)
    return df.dropna(subset=["Molecule Atoms", "Target Atoms"]).copy()

def _find_existing_cluster(coords, clusters, rmsd_cutoff):
    import docking_pipeline_utils as dpu
    for idx, cluster in enumerate(clusters):
        if dpu.ligand_rmsd(coords, cluster["coords"]) <= rmsd_cutoff: return idx
    return None

def _water_diagnostics(pose, lig_res_idx, acceptor_atom_name):
    """
    Return quick water diagnostics for debugging H-bond failures:
    (num_waters, nearest_water_O_distance_to_acceptor_or_None).
    """
    import numpy as np
    import docking_pipeline_utils as dpu

    lig = pose.residue(lig_res_idx)
    acc_idx = dpu.atom_index_by_name(lig, acceptor_atom_name) if acceptor_atom_name else None
    if acc_idx is None:
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

def _closest_ligand_acceptor_to_any_water(pose, lig_res_idx):
    """
    Find the ligand acceptor atom nearest to any water oxygen.
    Returns (atom_name_or_None, min_distance_or_None, water_residue_or_None).
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
        acc_xyz = np.array([float(lig.xyz(acc_idx).x), float(lig.xyz(acc_idx).y), float(lig.xyz(acc_idx).z)], dtype=float)
        for resid in range(1, pose.total_residue() + 1):
            water = pose.residue(resid)
            if not water.is_water():
                continue
            for o_idx in range(1, water.natoms() + 1):
                if water.atom_type(o_idx).element() != "O":
                    continue
                o_xyz = np.array([float(water.xyz(o_idx).x), float(water.xyz(o_idx).y), float(water.xyz(o_idx).z)], dtype=float)
                dist = float(np.linalg.norm(o_xyz - acc_xyz))
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    best_atom = lig.atom_name(acc_idx).strip()
                    best_water = resid

    return (best_atom, best_dist, best_water)

def _rename_water_resnames_in_output_dir(output_dir, dpu_module):
    """
    Post-run compatibility sweep for downstream tools (e.g., LigandMPNN):
    enforce TP3 water residue names across all written PDBs.
    """
    count = 0
    for name in os.listdir(output_dir):
        if not name.lower().endswith(".pdb"):
            continue
        path = os.path.join(output_dir, name)
        try:
            dpu_module.rename_water_resnames_to_tp3(path)
            count += 1
        except Exception as exc:
            logger.warning("Failed TP3 rename sweep for %s: %s", path, exc)
    return count

def _passes_hbond_ideal_window(hbond_result, params):
    """
    Optional stricter gate around ideal geometry.
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
    Always apply ideal-geometry checks (used for final post-pack filtering).
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

def align_and_save_conformers(
    pose,
    res,
    target_res_idx,
    path_to_conformers,
    df,
    params,
    output_dirs=None,
    rename_water_to_tp3=True,
    cluster_enabled=True,
    cluster_rmsd_cutoff=0.75,
    geometry_csv_path=None,
    output_name_prefix="",
):
    import pyrosetta
    import pyrosetta.rosetta.protocols.grafting as graft
    import pyrosetta.rosetta.protocols.rigid as rigid_moves
    import pyrosetta.rosetta.core.scoring as scoring
    import docking_pipeline_utils as dpu
    import conformer_prep

    packer = setup_packer_task()
    sf_all = pyrosetta.get_fa_scorefxn()
    sf_cst = pyrosetta.get_fa_scorefxn()
    sf_cst.set_weight(scoring.atom_pair_constraint, float(params.get("hbond_constraint_weight", 1.0)))

    min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mm = pyrosetta.MoveMap()
    mm.set_jump(True)
    min_mover.movemap(mm)
    min_mover.score_function(sf_cst)

    stats = {"total": 0, "passed_score": 0, "clustered": 0}
    clusters = []
    geometry_rows = []
    max_tries = int(params.get("max_tries", 30))
    debug_every = int(params.get("debug_every", 10))

    for conf_idx, conf in enumerate(
        conformer_prep.yield_ligand_poses(df, path_to_conformers, False, params['lig_res_num']),
        start=1
    ):
        if not conf: continue
        if conf_idx == 1 or conf_idx % debug_every == 0:
            logger.info(
                "Processing conformer %d (conf_num=%s). totals: seen=%d passed_score=%d clustered=%d",
                conf_idx, getattr(conf, "conf_num", "NA"), stats["total"], stats["passed_score"], stats["clustered"]
            )

        # Support both old and new conformer attribute names.
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

        if molecule_atoms is None or target_atoms is None or len(molecule_atoms) < 1 or len(target_atoms) < 1:
            logger.error("Conformer has empty alignment atom lists; skipping.")
            continue

        conf.align_to_target(res)
        grafted_pose = graft.insert_pose_into_pose(pose, conf.pose, len(pose.chain_sequence(1)))
        lig_idx = len(pose.chain_sequence(1)) + 1

        keep_candidate = False
        accepted_try_idx = None
        last_hbond_result = None
        strict_ok = False
        out_path = None
        score = None
        for try_idx in range(1, max_tries + 1):
            copy_pose = grafted_pose.clone()
            rigid_moves.RigidBodyPerturbMover(params['jump_num'], params['rotation'], params['translation']).apply(copy_pose)

            acceptor_name = molecule_atoms[0]
            neighbor_names = list(molecule_atoms[1:3])
            if params.get("use_closest_acceptor", True):
                closest_atom, _, _ = _closest_ligand_acceptor_to_any_water(copy_pose, lig_idx)
                if closest_atom is not None:
                    acceptor_name = closest_atom
                # No robust generic neighbor definition for a dynamic atom pick.
                neighbor_names = []

            if (not params.get("use_closest_acceptor", True)) and target_res_idx is not None and target_res_idx > 0:
                add_hbond_constraint(copy_pose, lig_idx, acceptor_name, target_res_idx, target_atoms[0], params['hbond_ideal'])
            if params.get("use_closest_acceptor", True) and params.get("dynamic_water_pull_constraint", True):
                add_hbond_constraint_to_water(
                    copy_pose,
                    lig_idx,
                    acceptor_name,
                    dist_ideal=params['hbond_ideal'],
                    dist_sd=params.get('hbond_constraint_sd', 0.3),
                    capture_max=params.get('hbond_constraint_capture_max', 8.0),
                )
            auto_setup_water_constraints(copy_pose, lig_idx)
            t_min0 = time.time()
            min_mover.apply(copy_pose)
            min_dt = time.time() - t_min0
            if min_dt > float(params.get("slow_min_threshold_sec", 5.0)):
                logger.info(
                    "Slow minimization: conformer %d (conf_num=%s), try %d/%d took %.2fs",
                    conf_idx, getattr(conf, "conf_num", "NA"), try_idx, max_tries, min_dt
                )

            hbond_result = dpu.evaluate_hbond_geometry(copy_pose, lig_idx, acceptor_name, neighbor_names,
                params['hbond_min'], params['hbond_max'], params['hbond_ideal'], params['donor_angle'], params['acceptor_angle'], params['quality_min'])
            last_hbond_result = hbond_result
            strict_ok = _passes_hbond_ideal_window(hbond_result, params)
            if ((hbond_result["passed"] and strict_ok) or (not params['use_hbond_filter'])):
                keep_candidate = True
                accepted_try_idx = try_idx
                break
            if try_idx == max_tries:
                num_waters, nearest_water_dist = _water_diagnostics(copy_pose, lig_idx, molecule_atoms[0])
                logger.info(
                    "Conformer %d (conf_num=%s) failed hbond filter after %d tries; best quality=%.3f dist=%s donor_angle=%s acceptor_angle=%s water_res=%s waters=%d nearest_water_O=%.3f strict_window=%s",
                    conf_idx,
                    getattr(conf, "conf_num", "NA"),
                    max_tries,
                    float(hbond_result.get("quality", 0.0)),
                    f"{hbond_result.get('distance', None):.3f}" if hbond_result.get("distance", None) is not None else "None",
                    f"{hbond_result.get('donor_angle', None):.1f}" if hbond_result.get("donor_angle", None) is not None else "None",
                    f"{hbond_result.get('acceptor_angle', None):.1f}" if hbond_result.get("acceptor_angle", None) is not None else "None",
                    hbond_result.get("water_residue", None),
                    num_waters,
                    nearest_water_dist if nearest_water_dist is not None else -1.0,
                    strict_ok,
                )

        if not keep_candidate:
            geometry_rows.append({
                "conf_idx": conf_idx,
                "conf_num": getattr(conf, "conf_num", "NA"),
                "accepted": False,
                "accepted_try": accepted_try_idx,
                "passed_score": False,
                "saved_cluster": False,
                "score": None,
                "distance": None if last_hbond_result is None else last_hbond_result.get("distance", None),
                "donor_angle": None if last_hbond_result is None else last_hbond_result.get("donor_angle", None),
                "acceptor_angle": None if last_hbond_result is None else last_hbond_result.get("acceptor_angle", None),
                "quality": None if last_hbond_result is None else last_hbond_result.get("quality", None),
                "water_residue": None if last_hbond_result is None else last_hbond_result.get("water_residue", None),
                "strict_window_passed": strict_ok,
                "output_pdb": None,
            })
            stats["total"] += 1; continue

        packer.apply(copy_pose)
        score = sf_all(copy_pose)
        postpack_hbond_result = dpu.evaluate_hbond_geometry(
            copy_pose,
            lig_idx,
            acceptor_name,
            neighbor_names,
            params['hbond_min'],
            params['hbond_max'],
            params['hbond_ideal'],
            params['donor_angle'],
            params['acceptor_angle'],
            params['quality_min'],
        )
        saved_cluster = False
        postpack_pass = bool(postpack_hbond_result.get("passed", False))
        postpack_strict_ok = (
            _passes_hbond_ideal_window_explicit(postpack_hbond_result, params)
            if params.get("enforce_final_ideal_geometry", True) else True
        )
        if (not params['use_hbond_filter']) or (postpack_pass and postpack_strict_ok):
            if score < params['max_pass_score']:
                stats["passed_score"] += 1
                coords = dpu.ligand_heavy_atom_coords(copy_pose, lig_idx)
                c_idx = _find_existing_cluster(coords, clusters, cluster_rmsd_cutoff) if cluster_enabled else None
                if c_idx is None:
                    out_dir = output_dirs["pass_score_repacked_dirs"][output_dirs["current_pass"]] if output_dirs else "output"
                    if not os.path.exists(out_dir): os.makedirs(out_dir, exist_ok=True)
                    out = os.path.join(out_dir, f"{output_name_prefix}rep_{stats['total']}.pdb")
                    dpu.dump_pose_pdb(copy_pose, out, rename_water=rename_water_to_tp3)
                    clusters.append({"coords": coords, "score": score, "path": out})
                    stats["clustered"] += 1
                    out_path = out
                    saved_cluster = True
        else:
            logger.info(
                "Post-pack geometry reject: conformer %d (conf_num=%s) dist=%s donor=%s acceptor=%s passed=%s strict=%s",
                conf_idx,
                getattr(conf, "conf_num", "NA"),
                f"{postpack_hbond_result.get('distance', None):.3f}" if postpack_hbond_result.get("distance", None) is not None else "None",
                f"{postpack_hbond_result.get('donor_angle', None):.1f}" if postpack_hbond_result.get("donor_angle", None) is not None else "None",
                f"{postpack_hbond_result.get('acceptor_angle', None):.1f}" if postpack_hbond_result.get("acceptor_angle", None) is not None else "None",
                postpack_pass,
                postpack_strict_ok,
            )
        geometry_rows.append({
            "conf_idx": conf_idx,
            "conf_num": getattr(conf, "conf_num", "NA"),
            "accepted": True,
            "accepted_try": accepted_try_idx,
            "passed_score": bool(score is not None and score < params['max_pass_score']),
            "saved_cluster": saved_cluster,
            "score": score,
            "distance": postpack_hbond_result.get("distance", None),
            "donor_angle": postpack_hbond_result.get("donor_angle", None),
            "acceptor_angle": postpack_hbond_result.get("acceptor_angle", None),
            "quality": postpack_hbond_result.get("quality", None),
            "water_residue": postpack_hbond_result.get("water_residue", None),
            "strict_window_passed": _passes_hbond_ideal_window(postpack_hbond_result, params),
            "postpack_geometry_passed": postpack_pass,
            "postpack_ideal_passed": postpack_strict_ok,
            "output_pdb": out_path,
        })
        stats["total"] += 1
    if geometry_csv_path:
        pd.DataFrame(geometry_rows).to_csv(geometry_csv_path, index=False)
        logger.info("Wrote H-bond geometry CSV: %s (%d rows)", geometry_csv_path, len(geometry_rows))
    return stats

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("config_file", nargs="?", default="config.txt")
    parser.add_argument("array_index", nargs="?", default=None, type=int)
    args = parser.parse_args(argv)

    # --- Setup Python Paths ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Add script dir and parent/legacy dirs to path
    for p in [script_dir, os.path.join(script_dir, ".."), os.path.join(script_dir, "..", "legacy")]:
        if p not in sys.path: sys.path.insert(0, p)

    config = ConfigParser()
    config.read(args.config_file)
    def_cfg, spec_cfg = config["DEFAULT"], config["grade_conformers"]
    array_index, array_count = _resolve_array_runtime(args.array_index, spec_cfg)

    import docking_pipeline_utils as dpu, alignment, conformer_prep, collision_check, pyrosetta
    from pyrosetta.rosetta.core.pose import Pose
    pyrosetta.init("-mute all")

    logger.info("Setting up PDBs...")
    params_list = def_cfg.get("ParamsList", "").split()
    pre_pose = Pose()
    post_pose = Pose()
    
    # Load with non-standard residues if params exist
    if params_list:
        pre_res_set = pyrosetta.generate_nonstandard_residue_set(pre_pose, params_list)
        pyrosetta.pose_from_file(pre_pose, pre_res_set, def_cfg["PrePDBFileName"])
        post_res_set = pyrosetta.generate_nonstandard_residue_set(post_pose, params_list)
        pyrosetta.pose_from_file(post_pose, post_res_set, def_cfg["PostPDBFileName"])
    else:
        pyrosetta.pose_from_file(pre_pose, def_cfg["PrePDBFileName"])
        pyrosetta.pose_from_file(post_pose, def_cfg["PostPDBFileName"])

    # Shave
    for i in [int(x) for x in spec_cfg.get("GlycineShavePositions", "").split()]:
        pyrosetta.toolbox.mutate_residue(post_pose, i, "G")

    # Allow a compact distance config:
    # if HBondDistanceMin/Max are omitted, derive them from ideal +/- buffer.
    hbond_ideal = _cfg_float(spec_cfg, "HBondDistanceIdeal", 2.8)
    dist_buf = _cfg_float(spec_cfg, "HBondDistanceIdealBuffer", 0.5)
    if "HBondDistanceMin" in spec_cfg:
        hbond_min = _cfg_float(spec_cfg, "HBondDistanceMin", hbond_ideal - dist_buf)
    else:
        hbond_min = max(0.0, hbond_ideal - dist_buf)
    if "HBondDistanceMax" in spec_cfg:
        hbond_max = _cfg_float(spec_cfg, "HBondDistanceMax", hbond_ideal + dist_buf)
    else:
        hbond_max = hbond_ideal + dist_buf

    raw_quality_min = _cfg_clean(spec_cfg.get("HBondQualityMin", ""))
    quality_min = float(raw_quality_min) if raw_quality_min else 0.0

    run_params = {
        'bin_width': _cfg_float(spec_cfg, "BinWidth", 1.0), 'vdw_modifier': _cfg_float(spec_cfg, "VDW_Modifier", 0.7),
        'lig_res_num': _cfg_int(def_cfg, "LigandResidueNumber", 1), 'jump_num': _cfg_int(spec_cfg, "JumpNum", 2),
        'rotation': _cfg_float(spec_cfg, "Rotation", 25.0), 'translation': _cfg_float(spec_cfg, "Translation", 0.5),
        'max_pass_score': _cfg_float(spec_cfg, "MaxScore", -300.0),
        'hbond_ideal': hbond_ideal,
        'hbond_min': hbond_min,
        'hbond_max': hbond_max,
        'donor_angle': _cfg_float(spec_cfg, "HBondDonorAngleMin", 100.0),
        'acceptor_angle': _cfg_float(spec_cfg, "HBondAcceptorAngleMin", 90.0),
        'quality_min': quality_min,
        'use_hbond_filter': _cfg_bool(spec_cfg, "EnableHBondGeometryFilter", True),
        'use_closest_acceptor': _cfg_bool(spec_cfg, "UseClosestLigandAcceptor", True),
        'hbond_constraint_weight': _cfg_float(spec_cfg, "HBondConstraintWeight", 1.0),
        'hbond_constraint_sd': _cfg_float(spec_cfg, "HBondConstraintSD", 0.3),
        'hbond_constraint_capture_max': _cfg_float(spec_cfg, "HBondConstraintCaptureMax", 8.0),
        'dynamic_water_pull_constraint': _cfg_bool(spec_cfg, "EnableDynamicWaterPullConstraint", True),
        'enforce_ideal_window': _cfg_bool(spec_cfg, "EnforceHBondIdealWindow", False),
        'enforce_final_ideal_geometry': _cfg_bool(spec_cfg, "EnforceFinalIdealGeometry", True),
        'distance_ideal_buffer': dist_buf,
        'donor_ideal': _cfg_float(spec_cfg, "HBondDonorAngleIdeal", 180.0),
        'donor_ideal_buffer': _cfg_float(spec_cfg, "HBondDonorAngleBuffer", 25.0),
        'acceptor_ideal': _cfg_float(spec_cfg, "HBondAcceptorAngleIdeal", 180.0),
        'acceptor_ideal_buffer': _cfg_float(spec_cfg, "HBondAcceptorAngleBuffer", 30.0),
        'require_acceptor_ideal_window': _cfg_bool(spec_cfg, "RequireAcceptorIdealWindow", False),
        'postprocess_rename_water_sweep': _cfg_bool(spec_cfg, "PostprocessRenameWaterSweep", True),
        'cluster_enabled': _cfg_bool(spec_cfg, "EnablePoseClusteringInArrayTask", False),
        'cluster_rmsd_cutoff': _cfg_float(spec_cfg, "ClusterRMSDCutoff", 0.75),
        'max_tries': _cfg_int(spec_cfg, "MaxPerturbTries", 30),
        'debug_every': _cfg_int(spec_cfg, "DebugEveryN", 10),
        'slow_min_threshold_sec': _cfg_float(spec_cfg, "SlowMinimizationSeconds", 5.0)
    }
    logger.info(
        "HBond filter settings: enabled=%s closest_acceptor=%s cst_wt=%.2f min=%.3f max=%.3f ideal=%.3f donor_min=%.1f acceptor_min=%.1f quality_min=%.3f ideal_window=%s max_tries=%d",
        run_params['use_hbond_filter'],
        run_params['use_closest_acceptor'],
        run_params['hbond_constraint_weight'],
        run_params['hbond_min'],
        run_params['hbond_max'],
        run_params['hbond_ideal'],
        run_params['donor_angle'],
        run_params['acceptor_angle'],
        run_params['quality_min'],
        run_params['enforce_ideal_window'],
        run_params['max_tries'],
    )
    if not raw_quality_min:
        logger.info("HBondQualityMin not set; quality gate disabled (using 0.0).")
    if run_params['enforce_ideal_window']:
        logger.info(
            "HBond ideal window: dist=%.2f±%.2f donor>=%.1f%s",
            run_params['hbond_ideal'],
            run_params['distance_ideal_buffer'],
            run_params['donor_ideal'] - run_params['donor_ideal_buffer'],
            (
                f" acceptor>={run_params['acceptor_ideal'] - run_params['acceptor_ideal_buffer']:.1f}"
                if run_params['require_acceptor_ideal_window'] else ""
            ),
        )
    logger.info(
        "Final post-pack ideal geometry enforcement: %s",
        run_params['enforce_final_ideal_geometry'],
    )

    chain_letter = def_cfg["ChainLetter"]
    residue_number = _cfg_int(def_cfg, "ResidueNumber", 1)

    # Output directory: explicit grade_conformers.OutputDir wins,
    # otherwise default to SCRATCH_ROOT/docked if SCRATCH_ROOT is defined.
    output_dir = spec_cfg.get("OutputDir", "").strip()
    if not output_dir:
        scratch_root = def_cfg.get("SCRATCH_ROOT", "").strip()
        output_dir = os.path.join(scratch_root, "docked") if scratch_root else "output"
    os.makedirs(output_dir, exist_ok=True)
    output_dirs = {
        "current_pass": 0,
        "pass_score_repacked_dirs": [output_dir],
    }
    logger.info("Writing docked outputs to: %s", output_dir)
    base_geometry_name = (
        _cfg_clean(spec_cfg.get("GeometryCSVName", "hbond_geometry_summary.csv"))
        or "hbond_geometry_summary.csv"
    )
    if array_count > 1:
        stem, ext = os.path.splitext(base_geometry_name)
        base_geometry_name = f"{stem}_array{array_index:04d}{ext or '.csv'}"
    geometry_csv_path = os.path.join(output_dir, base_geometry_name)
    output_name_prefix = _cfg_clean(spec_cfg.get("OutputPrefix", ""))
    if not output_name_prefix and array_count > 1:
        output_name_prefix = f"a{array_index:04d}_"

    # Resolve target on pre_pose for alignment atom definitions.
    pre_target_idx = pre_pose.pdb_info().pdb2pose(chain_letter, residue_number)
    post_target_idx = post_pose.pdb_info().pdb2pose(chain_letter, residue_number)
    if pre_target_idx <= 0 or pre_target_idx > pre_pose.total_residue():
        raise ValueError(
            f"Could not map target residue {chain_letter}{residue_number} into pre pose. "
            f"pre_pose pdb2pose={pre_target_idx}, post_pose pdb2pose={post_target_idx}, "
            f"pre_pose.total_residue={pre_pose.total_residue()}"
        )
    target_res = pre_pose.residue(pre_target_idx)

    # Optional: only apply direct ligand-target atom constraints if that residue exists in post_pose.
    target_idx = post_target_idx if post_target_idx > 0 else None

    df = prepare_dataframe(
        def_cfg["CSVFileName"],
        def_cfg["PathToConformers"],
        target_res,
        run_params['lig_res_num'],
        _cfg_bool(def_cfg, "AutoGenerateAlignment", False),
    )
    df = _slice_dataframe_for_array(df, array_index, array_count)
    logger.info(
        "Array runtime: index=%d count=%d rows_assigned=%d",
        array_index,
        array_count,
        len(df),
    )
    if not run_params['cluster_enabled'] and array_count > 1:
        logger.info(
            "In-array clustering disabled (EnablePoseClusteringInArrayTask=False). "
            "After all array tasks complete, run: python cluster_docked_post_array.py %s",
            args.config_file
        )
    if df.empty:
        logger.info("No conformer rows assigned to this array task; exiting.")
        return
    
    stats = align_and_save_conformers(
        post_pose,
        target_res,
        target_idx,
        def_cfg["PathToConformers"],
        df,
        run_params,
        output_dirs=output_dirs,
        cluster_enabled=run_params['cluster_enabled'],
        cluster_rmsd_cutoff=run_params['cluster_rmsd_cutoff'],
        geometry_csv_path=geometry_csv_path,
        output_name_prefix=output_name_prefix,
    )

    if run_params['postprocess_rename_water_sweep']:
        renamed = _rename_water_resnames_in_output_dir(output_dir, dpu)
        logger.info("Postprocess TP3 rename sweep completed for %d PDB files in %s", renamed, output_dir)

    logger.info(f"Done. Clustered {stats['clustered']} results.")

if __name__ == "__main__":
    main(sys.argv[1:])
