#!/usr/bin/env python
"""
Master orchestrator for PYR1 ML dataset generation.

CORRECT WORKFLOW (Direct-to-Mutant Docking):
1. Generate ligand conformers (RDKit ETKDG + MMFF)
2. Create alignment table (identify H-bond acceptor atoms for docking)
3. Thread mutations onto WT PYR1 → mutant.pdb
4. Backbone-constrained relax of mutant → mutant_relaxed.pdb
5. Dock conformers to MUTANT pocket (NO glycine shaving!)
6. Cluster poses → extract statistics (convergence, clashes)
7. Relax best pose (skip if severe clash)
8. Run AF3 binary + ternary predictions / Aggregate all features into ML table

This script handles:
- Workflow orchestration with dependency management
- SLURM job submission (optional)
- Caching (skip completed pairs)
- Progress tracking and logging

Author: Claude Code (Whitehead Lab PYR1 Pipeline)
Date: 2026-02-16
"""

import argparse
import json
import os
import sys
import time
import hashlib
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# Project root: two levels up from this script (ml_modelling/scripts/ -> project root)
_SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = _SCRIPT_DIR.parent.parent


# ═══════════════════════════════════════════════════════════════════
# CACHE KEY GENERATION
# ═══════════════════════════════════════════════════════════════════

def get_pair_cache_key(ligand_smiles: str, variant_signature: str) -> str:
    """Generate unique cache key for (ligand, variant) pair."""
    combined = f"{ligand_smiles}_{variant_signature}"
    return hashlib.md5(combined.encode()).hexdigest()


def is_stage_complete(pair_cache: Path, stage: str) -> bool:
    """
    Check if a pipeline stage has completed successfully.

    Args:
        pair_cache: Path to pair cache directory
        stage: Stage name ('conformers', 'threading', 'docking', 'relax', 'af3_binary', 'af3_ternary')

    Returns:
        True if stage completed, False otherwise
    """
    metadata_path = pair_cache / 'metadata.json'

    if not metadata_path.exists():
        return False

    try:
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        stages = metadata.get('stages', {})
        if stage not in stages:
            return False

        return stages[stage].get('status') == 'complete'

    except Exception:
        return False


def mark_stage_complete(pair_cache: Path, stage: str, output_dir: str = None) -> None:
    """Mark a stage as complete in metadata.json."""
    metadata_path = pair_cache / 'metadata.json'

    # Load existing metadata
    if metadata_path.exists():
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
    else:
        metadata = {'stages': {}}

    # Update stage status
    metadata['stages'][stage] = {
        'status': 'complete',
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
        'output_dir': str(output_dir) if output_dir else None
    }

    # Save metadata
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)


def update_stage_metadata(pair_cache: Path, stage: str, extra_data: Dict) -> None:
    """Merge extra key-value pairs into a stage's entry in metadata.json."""
    metadata_path = pair_cache / 'metadata.json'

    if not metadata_path.exists():
        return

    with open(metadata_path, 'r') as f:
        metadata = json.load(f)

    if stage not in metadata.get('stages', {}):
        metadata.setdefault('stages', {})[stage] = {}

    metadata['stages'][stage].update(extra_data)

    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)


def extract_docking_stats(docking_dir: Path) -> Optional[Dict]:
    """
    Parse hbond_geometry_summary CSV(s) to extract docking statistics.

    For SLURM runs with multiple array tasks, concatenates all
    hbond_geometry_summary_array*.csv files.

    Returns dict with:
        - total_docking_attempts: row count
        - passed_score_count: rows with passed_score=True
        - saved_cluster_count: rows with saved_cluster=True
        - pass_rate: passed/total (KEY ML feature)
        - best_score: min score among passed rows
        - mean_quality: mean H-bond quality among passed rows
        - best_cluster_pdb: path to best-scoring saved PDB
    """
    # Find all geometry CSV files (single file or array shards)
    csv_files = sorted(docking_dir.glob('hbond_geometry_summary_array*.csv'))
    single_csv = docking_dir / 'hbond_geometry_summary.csv'
    if not csv_files and single_csv.exists():
        csv_files = [single_csv]

    if not csv_files:
        logger.warning(f"  No geometry CSV found in {docking_dir}")
        return None

    # Concatenate all CSV files
    dfs = []
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            dfs.append(df)
        except Exception as e:
            logger.warning(f"  Failed to parse {csv_file}: {e}")
    if not dfs:
        return None

    all_rows = pd.concat(dfs, ignore_index=True)

    total = len(all_rows)
    if total == 0:
        return None

    # Parse boolean columns (may be stored as strings)
    for col in ['passed_score', 'saved_cluster']:
        if col in all_rows.columns:
            all_rows[col] = all_rows[col].astype(str).str.strip().str.lower().isin(['true', '1'])

    passed = all_rows[all_rows.get('passed_score', pd.Series(dtype=bool)) == True]
    saved = all_rows[all_rows.get('saved_cluster', pd.Series(dtype=bool)) == True]

    passed_count = len(passed)
    saved_count = len(saved)
    pass_rate = passed_count / total if total > 0 else 0.0

    # Count actual unique cluster rep PDB files (each array task clusters independently,
    # so saved_cluster_count is the sum across tasks, NOT the global cluster count)
    num_unique_reps = 0
    if saved_count > 0 and 'output_pdb' in saved.columns:
        unique_pdbs = saved['output_pdb'].dropna().unique()
        num_unique_reps = len(unique_pdbs)

    # Best score among passed rows
    best_score = float(passed['score'].min()) if passed_count > 0 and 'score' in passed.columns else np.nan
    mean_quality = float(passed['quality'].mean()) if passed_count > 0 and 'quality' in passed.columns else np.nan

    # Best cluster PDB (lowest score among saved clusters)
    best_cluster_pdb = None
    if saved_count > 0 and 'score' in saved.columns and 'output_pdb' in saved.columns:
        best_row = saved.loc[saved['score'].idxmin()]
        pdb_path = best_row.get('output_pdb', '')
        if pdb_path and str(pdb_path) != 'nan':
            best_cluster_pdb = str(pdb_path)

    stats = {
        'total_docking_attempts': int(total),
        'passed_score_count': int(passed_count),
        'saved_cluster_count': int(saved_count),
        'num_cluster_reps': int(num_unique_reps),
        'pass_rate': round(pass_rate, 4),
        'best_score': round(best_score, 3) if not np.isnan(best_score) else None,
        'mean_quality': round(mean_quality, 4) if not np.isnan(mean_quality) else None,
        'best_cluster_pdb': best_cluster_pdb,
    }

    logger.info(f"  Docking stats: {passed_count}/{total} passed ({pass_rate:.1%}), "
                f"{num_unique_reps} cluster reps, best score={best_score:.2f}")

    return stats


# ═══════════════════════════════════════════════════════════════════
# RELAX STAGE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

def find_top_docked_pdbs(pair_cache: Path, max_n: int = 20) -> List[Path]:
    """
    Find top N docked PDBs for relaxation, preferring global cluster representatives.

    Priority:
      1. Global clustering reps (cluster_*.pdb from run_clustering), sorted by best_score
      2. Geometry CSV saved_cluster=True rows, sorted by score
      3. Glob fallback (*rep_*.pdb, then *.pdb)
    """
    docking_dir = pair_cache / 'docking'

    # --- Source 1: Global cluster representatives (from run_clustering) ---
    stats_json = docking_dir / 'clustering_stats.json'
    if stats_json.exists():
        try:
            with open(stats_json, 'r') as f:
                stats = json.load(f)
            clusters = stats.get('clusters', [])
            # Sort by best_score (ascending = lowest energy first)
            clusters_sorted = sorted(clusters, key=lambda c: c.get('best_score', float('inf')))

            pdbs = []
            for c in clusters_sorted[:max_n]:
                pdb_path = docking_dir / f"cluster_{c['cluster_id']}.pdb"
                if pdb_path.exists():
                    pdbs.append(pdb_path)
            if pdbs:
                logger.info(f"  Found {len(pdbs)} cluster reps from global clustering "
                            f"({len(clusters)} total clusters)")
                return pdbs
        except Exception as e:
            logger.warning(f"  Failed to read clustering_stats.json: {e}")

    # --- Source 2: Geometry CSV saved_cluster=True rows ---
    csv_files = sorted(docking_dir.glob('hbond_geometry_summary_array*.csv'))
    single_csv = docking_dir / 'hbond_geometry_summary.csv'
    if not csv_files and single_csv.exists():
        csv_files = [single_csv]

    if csv_files:
        dfs = []
        for csv_file in csv_files:
            try:
                dfs.append(pd.read_csv(csv_file))
            except Exception:
                continue

        if dfs:
            all_rows = pd.concat(dfs, ignore_index=True)
            for col in ['saved_cluster']:
                if col in all_rows.columns:
                    all_rows[col] = all_rows[col].astype(str).str.strip().str.lower().isin(['true', '1'])

            saved = all_rows[all_rows.get('saved_cluster', pd.Series(dtype=bool)) == True]

            if len(saved) > 0 and 'score' in saved.columns and 'output_pdb' in saved.columns:
                saved_sorted = saved.sort_values('score', ascending=True).head(max_n)
                pdbs = []
                for _, row in saved_sorted.iterrows():
                    pdb_path = Path(str(row['output_pdb']))
                    if pdb_path.exists():
                        pdbs.append(pdb_path)
                if pdbs:
                    logger.info(f"  Found {len(pdbs)} top docked PDBs from geometry CSV")
                    return pdbs

    # --- Source 3: Glob fallback ---
    rep_pdbs = sorted(docking_dir.glob('*rep_*.pdb'))
    if not rep_pdbs:
        rep_pdbs = sorted(docking_dir.glob('*.pdb'))

    result = rep_pdbs[:max_n]
    if result:
        logger.info(f"  Found {len(result)} docked PDBs (glob fallback)")
    return result


def find_ligand_params(pair_cache: Path) -> Optional[Path]:
    """Find the ligand .params file in conformers_params directory."""
    params_dir = pair_cache / 'conformers_params'
    params_files = list(params_dir.glob('*/*.params'))
    if params_files:
        return params_files[0]
    # Also check directly
    params_files = list(params_dir.glob('*.params'))
    if params_files:
        return params_files[0]
    return None


def run_relax_single(
    input_pdb: Path,
    output_pdb: Path,
    ligand_params: Path,
    xml_path: str,
    ligand_chain: str = 'B',
    water_chain: str = 'D',
    timeout: int = 600
) -> bool:
    """
    Run relax on a single docked PDB.

    Args:
        input_pdb: Docked structure PDB
        output_pdb: Output relaxed PDB
        ligand_params: Ligand .params file
        xml_path: Interface scoring XML path
        ligand_chain: Ligand chain ID
        water_chain: Water chain ID
        timeout: Max seconds per structure

    Returns:
        True if successful
    """
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    relax_script = str(PROJECT_ROOT / 'design' / 'rosetta' / 'relax_general_universal.py')

    cmd = [
        'python',
        relax_script,
        str(input_pdb),
        str(output_pdb),
        str(ligand_params),
        '--xml_path', xml_path,
        '--ligand_chain', ligand_chain,
        '--water_chain', water_chain,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
        return True
    except subprocess.TimeoutExpired:
        logger.warning(f"  Relax timed out ({timeout}s): {input_pdb.name}")
        return False
    except subprocess.CalledProcessError as e:
        logger.warning(f"  Relax failed for {input_pdb.name}: {e.stderr[:200] if e.stderr else 'unknown error'}")
        return False


def run_relax_slurm(
    pdbs: List[Path],
    relax_dir: Path,
    ligand_params: Path,
    xml_path: str,
    ligand_chain: str = 'B',
    water_chain: str = 'D'
) -> Optional[str]:
    """
    Submit relax jobs to SLURM via manifest file.

    Writes a TSV manifest and submits an array job.

    Returns:
        SLURM job ID if successful, None otherwise
    """
    relax_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = relax_dir / 'relax_manifest.tsv'

    # Write manifest: one row per PDB
    with open(manifest_path, 'w') as f:
        for i, pdb in enumerate(pdbs):
            output_pdb = relax_dir / f'relaxed_{i}.pdb'
            f.write(f"{pdb}\t{output_pdb}\t{ligand_params}\t{xml_path}\t{ligand_chain}\t{water_chain}\n")

    n_tasks = len(pdbs)
    slurm_script = str(PROJECT_ROOT / 'ml_modelling' / 'scripts' / 'submit_relax_ml.sh')

    cmd = [
        'sbatch',
        f'--array=1-{n_tasks}',
        slurm_script,
        str(manifest_path),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        job_id = result.stdout.strip().split()[-1]
        logger.info(f"  Submitted relax array job: {job_id} ({n_tasks} tasks)")
        return job_id
    except subprocess.CalledProcessError as e:
        logger.error(f"  SLURM relax submission failed: {e.stderr}")
        return None


def parse_relax_scores(score_file: Path) -> Optional[Dict]:
    """
    Parse a relax score file (colon-delimited key: value format).

    Returns dict with all score keys, numeric values converted to float.
    """
    if not score_file.exists():
        return None

    scores = {}
    try:
        with open(score_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('SCORES:') or not line:
                    continue
                if ':' in line:
                    key, _, value = line.partition(':')
                    key = key.strip()
                    value = value.strip()
                    # Try numeric conversion
                    try:
                        scores[key] = float(value)
                    except ValueError:
                        scores[key] = value  # Keep as string ("yes"/"no"/"N/A")
        return scores if scores else None
    except Exception as e:
        logger.warning(f"  Failed to parse score file {score_file}: {e}")
        return None


def aggregate_relax_scores(score_dicts: List[Dict]) -> Dict:
    """
    Compute distributional features across all relaxed structures.

    For numeric metrics: best (min), mean, std, median
    For string metrics ("yes"/"no"): fraction_yes
    """
    if not score_dicts:
        return {}

    # Collect all keys
    all_keys = set()
    for d in score_dicts:
        all_keys.update(d.keys())

    aggregated = {}
    individual_scores = score_dicts  # Store raw scores for metadata

    for key in sorted(all_keys):
        values = [d[key] for d in score_dicts if key in d]
        if not values:
            continue

        # Check if numeric
        numeric_vals = []
        string_vals = []
        for v in values:
            if isinstance(v, (int, float)) and not np.isnan(v):
                numeric_vals.append(float(v))
            elif isinstance(v, str) and v.lower() in ('yes', 'no', '1', '0'):
                string_vals.append(v)
            elif isinstance(v, str):
                try:
                    numeric_vals.append(float(v))
                except ValueError:
                    string_vals.append(v)

        if numeric_vals:
            arr = np.array(numeric_vals)
            aggregated[f'{key}_best'] = round(float(np.min(arr)), 4)
            aggregated[f'{key}_mean'] = round(float(np.mean(arr)), 4)
            aggregated[f'{key}_std'] = round(float(np.std(arr)), 4)
            aggregated[f'{key}_median'] = round(float(np.median(arr)), 4)

        if string_vals:
            yes_count = sum(1 for v in string_vals if v.lower() in ('yes', '1'))
            aggregated[f'{key}_fraction_yes'] = round(yes_count / len(string_vals), 4) if string_vals else 0.0

    aggregated['n_structures_relaxed'] = len(score_dicts)

    return {
        'aggregated': aggregated,
        'individual_scores': individual_scores,
    }


# ═══════════════════════════════════════════════════════════════════
# PIPELINE STAGE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

def run_conformer_generation(
    ligand_name: str,
    ligand_smiles: str,
    output_dir: Path,
    num_conformers: int = 10
) -> bool:
    """
    Generate ligand conformers using ligand_conformers module.

    Args:
        ligand_name: Ligand name (for logging)
        ligand_smiles: SMILES string
        output_dir: Output directory
        num_conformers: Number of final conformers to generate

    Returns:
        True if successful, False otherwise
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  Generating {num_conformers} conformers for {ligand_name}...")

    cmd = [
        'python', '-m', 'ligand_conformers',
        '--input', ligand_smiles,
        '--input-type', 'smiles',
        '--outdir', str(output_dir),
        '--num-confs', str(num_conformers * 3),  # Generate 3× and select best
        '--k-final', str(num_conformers),
        '--cluster-rmsd', '1.0'
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  ✓ Conformers generated: {output_dir / 'conformers_final.sdf'}")
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Conformer generation failed: {e.stderr}")
        return False


def run_alignment_table_generation(
    conformers_sdf: Path,
    output_csv: Path,
    output_conformers_dir: Path,
    target_atom_triplets: str = "O2-C11-C9;O2-C9-C11"
) -> bool:
    """
    Generate alignment table CSV identifying H-bond acceptor atoms.

    Args:
        conformers_sdf: Path to conformers SDF file
        output_csv: Path for output alignment CSV
        output_conformers_dir: Directory for params/PDB outputs
        target_atom_triplets: Semicolon-separated target triplets on template ligand

    Returns:
        True if successful, False otherwise
    """
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_conformers_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  Generating alignment table for docking...")

    # Create temporary config file for create_table.py
    config_path = output_csv.parent / 'create_table_config.txt'
    with open(config_path, 'w') as f:
        f.write(f"""[DEFAULT]
PathToConformers = {output_conformers_dir}
CSVFileName = {output_csv}

[create_table]
MoleculeSDFs = {conformers_sdf}
CSVFileName = {output_csv}
PathToConformers = {output_conformers_dir}
UseMoleculeID = False
NoName = True
DynamicAcceptorAlignment = True
IncludeReverseNeighborOrder = False
MaxDynamicAlignments = 20
AcceptorMode = generic
TargetAtomTriplets = {target_atom_triplets}
DynamicAlignmentDebug = False
""")

    # Run create_table.py
    cmd = [
        'python',
        str(PROJECT_ROOT / 'docking' / 'scripts' / 'create_table.py'),
        str(config_path)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  ✓ Alignment table created: {output_csv}")
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Alignment table generation failed: {e.stderr}")
        return False


def run_mutation_threading(
    template_pdb: str,
    variant_signature: str,
    output_pdb: Path,
    chain: str = 'A'
) -> bool:
    """
    Thread mutations onto WT PYR1 template.

    Args:
        template_pdb: Path to WT PYR1 template
        variant_signature: Mutation signature (e.g., "59K;120A;160G") or empty for wildtype
        output_pdb: Path for output mutant PDB
        chain: Chain ID to mutate

    Returns:
        True if successful, False otherwise
    """
    import shutil
    import pandas as pd

    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    # Handle wildtype (no mutations)
    if pd.isna(variant_signature) or not variant_signature or variant_signature.strip() == "":
        logger.info(f"  Wildtype (no mutations) - copying template")
        try:
            shutil.copy2(template_pdb, output_pdb)
            logger.info(f"  ✓ Wildtype structure created: {output_pdb}")
            return True
        except Exception as e:
            logger.error(f"  ✗ Failed to copy template: {e}")
            return False

    logger.info(f"  Threading mutations: {variant_signature}")

    cmd = [
        'python',
        str(PROJECT_ROOT / 'scripts' / 'thread_variant_to_pdb.py'),
        '--template', template_pdb,
        '--signature', str(variant_signature),
        '--output', str(output_pdb),
        '--chain', chain,
        '--no-deletion',  # 3QN1 PDB preserves original WT numbering (gap at 69-70)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  ✓ Mutant structure created: {output_pdb}")
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Threading failed: {e.stderr}")
        return False


def run_constrained_relax(
    input_pdb: Path,
    output_pdb: Path,
    variant_signature: str,
    params: Optional[Path] = None,
    timeout: int = 300
) -> bool:
    """
    Repack sidechains in a 10 A shell around mutated residues.

    Uses PackRotamersMover (not FastRelax) for fast discrete rotamer
    optimization. Completes in seconds.

    Args:
        input_pdb: Threaded mutant PDB (protein-only, from threading stage)
        output_pdb: Output repacked PDB path
        variant_signature: Mutation signature (e.g., "59K;120A;160G")
        params: Optional .params file for extra_res_fa
        timeout: Max seconds (default 300 = 5 min)

    Returns:
        True if successful, False otherwise
    """
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    relax_script = str(PROJECT_ROOT / 'ml_modelling' / 'scripts' / 'constrained_relax.py')

    cmd = [
        'python',
        relax_script,
        str(input_pdb),
        str(output_pdb),
        '--mutations', variant_signature,
    ]

    if params:
        cmd.extend(['--params', str(params)])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
        return True
    except subprocess.TimeoutExpired:
        logger.warning(f"  Constrained relax timed out ({timeout}s): {input_pdb.name}")
        return False
    except subprocess.CalledProcessError as e:
        logger.warning(f"  Constrained relax failed: {e.stderr[:200] if e.stderr else 'unknown error'}")
        return False


def run_constrained_relax_slurm(
    input_pdb: Path,
    output_pdb: Path,
    variant_signature: str,
    params: Optional[Path] = None
) -> Optional[str]:
    """
    Submit sidechain repack to SLURM.

    Args:
        input_pdb: Threaded mutant PDB
        output_pdb: Output repacked PDB path
        variant_signature: Mutation signature (e.g., "59K;120A;160G")
        params: Optional .params file

    Returns:
        SLURM job ID if successful, None otherwise
    """
    slurm_script = str(PROJECT_ROOT / 'ml_modelling' / 'scripts' / 'submit_constrained_relax.sh')

    cmd = [
        'sbatch',
        '--output', str(output_pdb.parent / 'constrained_relax_%j.log'),
        '--error', str(output_pdb.parent / 'constrained_relax_%j.err'),
        slurm_script,
        str(input_pdb),
        str(output_pdb),
        variant_signature,
    ]

    if params:
        cmd.append(str(params))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        job_id = result.stdout.strip().split()[-1]
        logger.info(f"  Submitted constrained relax job: {job_id}")
        return job_id
    except subprocess.CalledProcessError as e:
        logger.error(f"  SLURM constrained relax submission failed: {e.stderr}")
        return None


def run_docking(
    mutant_pdb: Path,
    conformers_sdf: Path,
    alignment_csv: Path,
    conformers_dir: Path,
    output_dir: Path,
    reference_pdb: str,
    docking_repeats: int = 50,
    use_slurm: bool = False,
    array_tasks: int = 10,
    reference_chain: str = 'X',
    reference_residue: int = 1
) -> Optional[str]:
    """
    Run docking to mutant pocket.

    Args:
        mutant_pdb: Pre-threaded mutant structure (NO ligand)
        conformers_sdf: Ligand conformers SDF
        alignment_csv: Alignment table CSV
        conformers_dir: Conformers params/PDB directory
        output_dir: Output directory
        reference_pdb: Template PDB WITH ligand for alignment (e.g., 3QN1_H2O.pdb)
        docking_repeats: Number of docking attempts per conformer
        use_slurm: Submit as SLURM array job
        array_tasks: Number of SLURM array tasks
        reference_chain: Chain of template ligand in reference PDB (default: 'X')
        reference_residue: Residue number of template ligand (default: 1)

    Returns:
        Job ID if SLURM, 'local' if local, None if failed
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  Docking to mutant pocket ({docking_repeats} repeats)...")

    # Create config file
    # CRITICAL: ChainLetter/ResidueNumber point to the template LIGAND in reference_pdb
    # This ligand provides the target atom coordinates (O2, C11, C9) for SVD alignment
    # The mutant_pdb is used as the docking target (has no ligand)
    config_path = output_dir / 'docking_config.txt'

    # Enable in-loop clustering for local runs (skip external clustering script)
    enable_clustering = not use_slurm

    # Need A8T.params for loading the reference PDB with template ligand
    params_file = str(PROJECT_ROOT / "docking" / "ligand_alignment" / "files_for_PYR1_docking" / "A8T.params")

    with open(config_path, 'w') as f:
        f.write(f"""[DEFAULT]
CSVFileName = {alignment_csv}
PathToConformers = {conformers_dir}
ChainLetter = {reference_chain}
ResidueNumber = {reference_residue}
LigandResidueNumber = 1
AutoGenerateAlignment = False
PrePDBFileName = {reference_pdb}
ParamsList = {params_file}

[mutant_docking]
MutantPDB = {mutant_pdb}
LigandSDF = {conformers_sdf}
OutputDir = {output_dir}
DockingRepeats = {docking_repeats}
ArrayTaskCount = {array_tasks}
EnablePoseClusteringInArrayTask = {str(enable_clustering)}
ClusterRMSDCutoff = 2.0
EnablePocketProximityFilter = True
PocketMaxDistance = 5.0
MaxScore = -200
UseRelativeScoring = False
""")

    if use_slurm:
        # Estimate walltime based on effective repeats per array task.
        # Work-unit distribution: each task handles ~(docking_repeats / array_tasks)
        # repeats per conformer. Each docking attempt ~3-5 sec, plus conformer
        # loading and PyRosetta overhead. Generous to avoid timeouts.
        effective_repeats = max(1, -(-docking_repeats // array_tasks))  # ceil division
        if effective_repeats <= 5:
            walltime = '02:00:00'
        elif effective_repeats <= 20:
            walltime = '04:00:00'
        elif effective_repeats <= 50:
            walltime = '06:00:00'
        else:
            walltime = '12:00:00'

        # Submit SLURM array job
        cmd = [
            'sbatch',
            '--partition', 'amilan',
            '--qos', 'normal',
            '--job-name', f"dock_{mutant_pdb.stem}",
            '--output', str(output_dir / 'docking_%a.log'),
            '--array', f'0-{array_tasks-1}',
            '--time', walltime,
            '--cpus-per-task', '4',
            '--mem', '8G',
            str(PROJECT_ROOT / 'docking' / 'scripts' / 'submit_docking_mutant.sh'),
            str(config_path)
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            job_id = result.stdout.strip().split()[-1]
            logger.info(f"  ✓ Submitted docking job: {job_id}")
            return job_id

        except subprocess.CalledProcessError as e:
            logger.error(f"  ✗ SLURM submission failed: {e.stderr}")
            return None

    else:
        # Run locally (single task, array index 0)
        cmd = [
            'python',
            str(PROJECT_ROOT / 'docking' / 'scripts' / 'grade_conformers_mutant_docking.py'),
            str(config_path),
            '0'  # Array index 0
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"  ✓ Docking complete (local)")
            return 'local'

        except subprocess.CalledProcessError as e:
            logger.error(f"  ✗ Docking failed: {e.stderr}")
            return None


def run_clustering(
    docking_output_dir: Path,
    rmsd_cutoff: float = 2.0
) -> bool:
    """
    Global re-clustering of docked poses with statistics for ML features.

    Reads geometry CSVs from the docking directory, filters to rows with
    valid PDB files (cluster reps from in-loop clustering), and runs global
    clustering to produce clustering_stats.json with convergence metrics.

    Args:
        docking_output_dir: Directory with docking results (geometry CSVs + PDBs)
        rmsd_cutoff: RMSD clustering cutoff (Å)

    Returns:
        True if successful, False otherwise
    """
    import tempfile

    # Read geometry CSVs
    csv_files = sorted(docking_output_dir.glob('hbond_geometry_summary*.csv'))
    if not csv_files:
        logger.warning("  No geometry CSVs found for clustering")
        return False

    dfs = []
    for csv_file in csv_files:
        try:
            dfs.append(pd.read_csv(csv_file))
        except Exception as e:
            logger.warning(f"  Failed to read {csv_file}: {e}")
    if not dfs:
        return False

    all_rows = pd.concat(dfs, ignore_index=True)

    # Filter to rows with existing PDB files (these are the in-loop cluster reps)
    valid = all_rows.dropna(subset=['output_pdb'])
    valid = valid[valid['output_pdb'].apply(lambda p: Path(str(p)).exists())]

    if len(valid) == 0:
        logger.warning("  No valid PDB files found for global re-clustering")
        return False

    logger.info(f"  Global re-clustering {len(valid)} poses (RMSD cutoff {rmsd_cutoff} Å)...")

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write filtered CSV for clustering script
            filtered_csv = Path(tmpdir) / 'poses_0.csv'
            valid.to_csv(filtered_csv, index=False)

            cmd = [
                'python',
                str(PROJECT_ROOT / 'docking' / 'scripts' / 'cluster_docked_with_stats.py'),
                '--input-dir', tmpdir,
                '--output-dir', str(docking_output_dir),
                '--rmsd-cutoff', str(rmsd_cutoff),
                '--pattern-prefix', 'poses_',
                '--save-all-clusters',
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Verify output
        stats_json = docking_output_dir / 'clustering_stats.json'
        if stats_json.exists():
            with open(stats_json, 'r') as f:
                stats = json.load(f)
            global_stats = stats.get('global_stats', {})
            logger.info(f"  ✓ Clustering: {global_stats.get('num_clusters', '?')} clusters, "
                        f"convergence {global_stats.get('convergence_ratio', 0):.1%}")

            if global_stats.get('clash_count', 0) > 0:
                logger.warning(f"    ⚠ {global_stats['clash_count']} clusters with clashes!")
            return True
        else:
            logger.warning("  clustering_stats.json not created")
            return False

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Clustering failed: {e.stderr}")
        return False


def should_skip_relax(cluster_output_dir: Path, pair_cache: Path = None) -> bool:
    """
    Check if relax should be skipped (severe clash detected).

    Checks clustering_stats.json first, then falls back to docking stats
    in metadata.json (for local runs where clustering_stats.json doesn't exist).

    Args:
        cluster_output_dir: Clustering output directory
        pair_cache: Pair cache directory (for metadata.json fallback)

    Returns:
        True if should skip (best_score > 0 indicates clash), False otherwise
    """
    # Try clustering_stats.json first
    stats_json = cluster_output_dir / 'clustering_stats.json'
    if stats_json.exists():
        try:
            with open(stats_json, 'r') as f:
                stats = json.load(f)
                best_score = stats.get('global_stats', {}).get('best_overall_score', 0)
                if best_score > 0:
                    logger.info(f"  Skipping relax: clash detected (score = {best_score:.2f})")
                    return True
                return False
        except Exception:
            pass

    # Fallback: check docking stats in metadata.json
    if pair_cache is not None:
        metadata_path = pair_cache / 'metadata.json'
        if metadata_path.exists():
            try:
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                best_score = metadata.get('stages', {}).get('docking', {}).get('best_score')
                if best_score is not None and best_score > 0:
                    logger.info(f"  Skipping relax: clash detected from docking stats (score = {best_score:.2f})")
                    return True
            except Exception:
                pass

    return False


# ═══════════════════════════════════════════════════════════════════
# MAIN ORCHESTRATION
# ═══════════════════════════════════════════════════════════════════

def process_single_pair(
    pair: Dict,
    cache_dir: Path,
    template_pdb: str,
    reference_pdb: str,
    use_slurm: bool = False,
    docking_repeats: int = 50,
    docking_arrays: int = 10,
    reference_chain: str = 'X',
    reference_residue: int = 1,
    af3_args: Optional[Dict] = None
) -> Dict:
    """
    Process a single (ligand, variant) pair through the full pipeline.

    Args:
        pair: Dict with pair_id, ligand_name, ligand_smiles, variant_signature, etc.
        cache_dir: Base cache directory
        template_pdb: WT PYR1 template PDB path (NO ligand, for threading)
        reference_pdb: Reference PDB WITH ligand for alignment (e.g., 3QN1_H2O.pdb)
        use_slurm: Use SLURM for parallelization
        docking_repeats: Number of docking repeats
        docking_arrays: Number of SLURM array tasks for docking
        reference_chain: Chain of template ligand in reference PDB (default: 'X')
        reference_residue: Residue number of template ligand (default: 1)
        af3_args: Dict with AF3 configuration (wt_sequence, templates, skip_af3, etc.)

    Returns:
        Dict with status and any errors
    """
    pair_id = pair['pair_id']
    pair_cache = cache_dir / pair_id

    # Inject AF3 args into pair for Stage 8a to access
    if af3_args:
        pair['_af3_args'] = af3_args

    logger.info(f"\n{'='*60}")
    logger.info(f"Processing pair: {pair_id}")
    logger.info(f"  Ligand: {pair['ligand_name']}")
    logger.info(f"  Variant: {pair['variant_name']} ({pair['variant_signature']})")
    logger.info(f"{'='*60}")

    # Create metadata file
    pair_cache.mkdir(parents=True, exist_ok=True)
    metadata_path = pair_cache / 'metadata.json'

    if not metadata_path.exists():
        with open(metadata_path, 'w') as f:
            json.dump({
                'pair_id': pair_id,
                'ligand_name': pair['ligand_name'],
                'ligand_smiles': pair['ligand_smiles'],
                'variant_name': pair['variant_name'],
                'variant_signature': pair['variant_signature'],
                'label': pair.get('label', np.nan),
                'label_tier': pair.get('label_tier', ''),
                'stages': {}
            }, f, indent=2)

    # ──────────────────────────────────────────────────────────────
    # STAGE 1: Conformer Generation
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'conformers'):
        logger.info("[1/8] Conformer Generation")
        conformer_dir = pair_cache / 'conformers'

        success = run_conformer_generation(
            ligand_name=pair['ligand_name'],
            ligand_smiles=pair['ligand_smiles'],
            output_dir=conformer_dir,
            num_conformers=10
        )

        if success:
            mark_stage_complete(pair_cache, 'conformers', str(conformer_dir))
        else:
            return {'status': 'FAILED', 'stage': 'conformers'}

    else:
        logger.info("[1/8] Conformer Generation: ✓ CACHED")

    conformers_sdf = pair_cache / 'conformers' / 'conformers_final.sdf'

    # ──────────────────────────────────────────────────────────────
    # STAGE 2: Alignment Table Generation
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'alignment_table'):
        logger.info("[2/8] Alignment Table Generation")
        alignment_csv = pair_cache / 'alignment_table.csv'
        conformers_params_dir = pair_cache / 'conformers_params'

        success = run_alignment_table_generation(
            conformers_sdf=conformers_sdf,
            output_csv=alignment_csv,
            output_conformers_dir=conformers_params_dir
        )

        if success:
            mark_stage_complete(pair_cache, 'alignment_table', str(alignment_csv))
        else:
            return {'status': 'FAILED', 'stage': 'alignment_table'}

    else:
        logger.info("[2/8] Alignment Table Generation: ✓ CACHED")

    alignment_csv = pair_cache / 'alignment_table.csv'
    conformers_params_dir = pair_cache / 'conformers_params'

    # ──────────────────────────────────────────────────────────────
    # STAGE 3: Mutation Threading
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'threading'):
        logger.info("[3/8] Mutation Threading")
        mutant_pdb = pair_cache / 'mutant.pdb'

        success = run_mutation_threading(
            template_pdb=template_pdb,
            variant_signature=pair['variant_signature'],
            output_pdb=mutant_pdb,
            chain='A'
        )

        if success:
            mark_stage_complete(pair_cache, 'threading', str(mutant_pdb))
        else:
            return {'status': 'FAILED', 'stage': 'threading'}

    else:
        logger.info("[3/8] Mutation Threading: ✓ CACHED")

    mutant_pdb = pair_cache / 'mutant.pdb'

    # ──────────────────────────────────────────────────────────────
    # STAGE 4: Backbone-Constrained Relax of Threaded Mutant
    # ──────────────────────────────────────────────────────────────
    mutant_relaxed_pdb = pair_cache / 'mutant_relaxed.pdb'
    a8t_params = str(PROJECT_ROOT / 'docking' / 'ligand_alignment' / 'files_for_PYR1_docking' / 'A8T.params')

    variant_sig = pair['variant_signature']
    is_wildtype = pd.isna(variant_sig) or not str(variant_sig).strip()

    if not is_stage_complete(pair_cache, 'constrained_relax'):
        if is_wildtype:
            # Wildtype: no mutations to repack — just copy mutant.pdb
            import shutil
            logger.info("[4/8] Sidechain Repack: SKIPPED (wildtype)")
            shutil.copy2(str(mutant_pdb), str(mutant_relaxed_pdb))
            mark_stage_complete(pair_cache, 'constrained_relax', str(mutant_relaxed_pdb))
        elif mutant_relaxed_pdb.exists():
            # SLURM re-entry: output file exists from a previous SLURM job
            logger.info("[4/8] Sidechain Repack: ✓ DONE (SLURM re-entry)")
            mark_stage_complete(pair_cache, 'constrained_relax', str(mutant_relaxed_pdb))
        else:
            logger.info("[4/8] Sidechain Repack")

            if use_slurm:
                job_id = run_constrained_relax_slurm(
                    input_pdb=mutant_pdb,
                    output_pdb=mutant_relaxed_pdb,
                    variant_signature=str(variant_sig),
                    params=Path(a8t_params),
                )
                if job_id:
                    logger.info(f"  ⏳ SLURM job submitted: {job_id}")
                    logger.info(f"  → Re-run orchestrator after job finishes")
                    return {'status': 'SLURM_SUBMITTED', 'job_id': job_id, 'stage': 'constrained_relax'}
                else:
                    return {'status': 'FAILED', 'stage': 'constrained_relax'}
            else:
                success = run_constrained_relax(
                    input_pdb=mutant_pdb,
                    output_pdb=mutant_relaxed_pdb,
                    variant_signature=str(variant_sig),
                    params=Path(a8t_params),
                )
                if success:
                    mark_stage_complete(pair_cache, 'constrained_relax', str(mutant_relaxed_pdb))
                else:
                    return {'status': 'FAILED', 'stage': 'constrained_relax'}
    else:
        logger.info("[4/8] Sidechain Repack: ✓ CACHED")

    # Use relaxed mutant for all downstream stages
    mutant_pdb = mutant_relaxed_pdb

    # ──────────────────────────────────────────────────────────────
    # STAGE 5: Docking to Mutant Pocket
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'docking'):
        docking_dir = pair_cache / 'docking'

        # SLURM re-entry: detect output from a previous SLURM run
        # (geometry CSVs or PDB files already exist → mark complete)
        if docking_dir.exists() and use_slurm:
            existing_csvs = list(docking_dir.glob('hbond_geometry_summary*.csv'))
            existing_pdbs = list(docking_dir.glob('*rep_*.pdb'))
            if existing_csvs or existing_pdbs:
                logger.info("[5/8] Docking to Mutant: ✓ DONE (SLURM re-entry, found %d CSVs, %d PDBs)",
                            len(existing_csvs), len(existing_pdbs))
                mark_stage_complete(pair_cache, 'docking', str(docking_dir))
            else:
                logger.info("[5/8] Docking to Mutant")

        if not is_stage_complete(pair_cache, 'docking'):
            logger.info("[5/8] Docking to Mutant")

            job_id = run_docking(
                mutant_pdb=mutant_pdb,
                conformers_sdf=conformers_sdf,
                alignment_csv=alignment_csv,
                conformers_dir=conformers_params_dir,
                output_dir=docking_dir,
                reference_pdb=reference_pdb,
                docking_repeats=docking_repeats,
                use_slurm=use_slurm,
                array_tasks=docking_arrays,
                reference_chain=reference_chain,
                reference_residue=reference_residue
            )

            if job_id:
                if use_slurm:
                    logger.info(f"  ⏳ SLURM job submitted: {job_id}")
                    logger.info(f"  → Re-run orchestrator after jobs finish to continue")
                    return {'status': 'SLURM_SUBMITTED', 'job_id': job_id, 'stage': 'docking'}
                else:
                    mark_stage_complete(pair_cache, 'docking', str(docking_dir))
            else:
                return {'status': 'FAILED', 'stage': 'docking'}

    else:
        logger.info("[5/8] Docking to Mutant: ✓ CACHED")

    docking_dir = pair_cache / 'docking'

    # ── Extract docking stats (for both fresh and cached runs) ──
    docking_stats_in_meta = False
    metadata_path_check = pair_cache / 'metadata.json'
    if metadata_path_check.exists():
        with open(metadata_path_check, 'r') as f:
            _meta = json.load(f)
        docking_stats_in_meta = 'pass_rate' in _meta.get('stages', {}).get('docking', {})

    if not docking_stats_in_meta:
        docking_stats = extract_docking_stats(docking_dir)
        if docking_stats:
            update_stage_metadata(pair_cache, 'docking', docking_stats)
        else:
            logger.warning("  Could not extract docking stats (no geometry CSV found)")

    # ──────────────────────────────────────────────────────────────
    # STAGE 6: Clustering with Statistics
    # ──────────────────────────────────────────────────────────────
    # Check if cluster representatives already exist (from in-loop clustering
    # during docking, or from a previous cluster_docked_post_array.py run)
    if is_stage_complete(pair_cache, 'clustering'):
        logger.info("[6/8] Clustering & Statistics: ✓ CACHED")
    else:
        # Check if geometry CSVs contain saved_cluster=True rows (in-loop clustering)
        has_cluster_reps = bool(find_top_docked_pdbs(pair_cache, max_n=1))
        if has_cluster_reps:
            logger.info("[6/8] Clustering & Statistics: in-loop clustering detected")
            mark_stage_complete(pair_cache, 'clustering', str(docking_dir))
        elif use_slurm:
            logger.info("[6/8] Clustering & Statistics: PENDING")
            logger.info("  → No cluster reps found. Run cluster_docked_post_array.py first")
            return {'status': 'NEEDS_CLUSTERING', 'stage': 'clustering'}
        else:
            logger.info("[6/8] Clustering & Statistics: ✓ DONE (in-loop)")
            mark_stage_complete(pair_cache, 'clustering', str(docking_dir))

    # Generate clustering statistics if missing (for ML feature extraction)
    stats_json = docking_dir / 'clustering_stats.json'
    if not stats_json.exists():
        logger.info("  → Generating clustering statistics...")
        run_clustering(docking_dir, rmsd_cutoff=2.0)

    cluster_dir = pair_cache / 'docking'  # Use docking dir since clustering is in-loop

    # ──────────────────────────────────────────────────────────────
    # STAGE 7: Relax (skip if severe clash)
    # ──────────────────────────────────────────────────────────────
    if should_skip_relax(cluster_dir, pair_cache=pair_cache):
        logger.info("[7/8] Rosetta Relax: SKIPPED (clash detected)")
        mark_stage_complete(pair_cache, 'relax', None)
        update_stage_metadata(pair_cache, 'relax', {'skipped_reason': 'clash'})

    elif not is_stage_complete(pair_cache, 'relax'):
        logger.info("[7/8] Rosetta Relax")
        relax_dir = pair_cache / 'relax'
        relax_dir.mkdir(parents=True, exist_ok=True)
        xml_path = str(PROJECT_ROOT / 'docking' / 'ligand_alignment' / 'scripts' / 'interface_scoring.xml')

        # Find top docked PDBs
        top_pdbs = find_top_docked_pdbs(pair_cache, max_n=20)
        if not top_pdbs:
            logger.warning("  No docked PDBs found - skipping relax")
            mark_stage_complete(pair_cache, 'relax', None)
            update_stage_metadata(pair_cache, 'relax', {'skipped_reason': 'no_pdbs'})
        else:
            # Find ligand params
            ligand_params = find_ligand_params(pair_cache)
            if not ligand_params:
                logger.error("  No ligand .params file found - cannot relax")
                return {'status': 'FAILED', 'stage': 'relax', 'error': 'no_params'}

            if use_slurm:
                # ── SLURM mode: check for existing scores (re-entry) or submit ──
                existing_scores = sorted(relax_dir.glob('relaxed_*_score.sc'))
                if existing_scores:
                    # Re-entry: SLURM jobs completed, aggregate scores
                    logger.info(f"  Found {len(existing_scores)} relax score files (SLURM re-entry)")
                    score_dicts = []
                    for sf in existing_scores:
                        parsed = parse_relax_scores(sf)
                        if parsed:
                            score_dicts.append(parsed)

                    if score_dicts:
                        agg = aggregate_relax_scores(score_dicts)
                        # mark_stage_complete first (it replaces the stage dict),
                        # then update_stage_metadata to merge scores into it
                        mark_stage_complete(pair_cache, 'relax', str(relax_dir))
                        update_stage_metadata(pair_cache, 'relax', agg.get('aggregated', {}))
                        update_stage_metadata(pair_cache, 'relax', {
                            'individual_scores': agg.get('individual_scores', [])
                        })
                        logger.info(f"  ✓ Relax aggregation complete: {len(score_dicts)} structures")
                    else:
                        logger.warning("  No valid relax scores parsed from SLURM outputs")
                else:
                    # Submit SLURM array job
                    job_id = run_relax_slurm(
                        pdbs=top_pdbs,
                        relax_dir=relax_dir,
                        ligand_params=ligand_params,
                        xml_path=xml_path,
                    )
                    if job_id:
                        return {'status': 'SLURM_SUBMITTED', 'job_id': job_id, 'stage': 'relax'}
                    else:
                        return {'status': 'FAILED', 'stage': 'relax', 'error': 'slurm_submit'}
            else:
                # ── Local mode: relax each PDB sequentially ──
                logger.info(f"  Relaxing {len(top_pdbs)} structures locally...")
                score_dicts = []
                for i, pdb in enumerate(top_pdbs):
                    output_pdb = relax_dir / f'relaxed_{i}.pdb'
                    logger.info(f"  [{i+1}/{len(top_pdbs)}] Relaxing {pdb.name}...")
                    success = run_relax_single(
                        input_pdb=pdb,
                        output_pdb=output_pdb,
                        ligand_params=ligand_params,
                        xml_path=xml_path,
                    )
                    if success:
                        score_file = relax_dir / f'relaxed_{i}_score.sc'
                        parsed = parse_relax_scores(score_file)
                        if parsed:
                            score_dicts.append(parsed)

                if score_dicts:
                    agg = aggregate_relax_scores(score_dicts)
                    # mark_stage_complete first (it replaces the stage dict),
                    # then update_stage_metadata to merge scores into it
                    mark_stage_complete(pair_cache, 'relax', str(relax_dir))
                    update_stage_metadata(pair_cache, 'relax', agg.get('aggregated', {}))
                    update_stage_metadata(pair_cache, 'relax', {
                        'individual_scores': agg.get('individual_scores', [])
                    })
                    logger.info(f"  ✓ Relax complete: {len(score_dicts)}/{len(top_pdbs)} structures scored")
                else:
                    logger.warning("  No relax scores obtained")
                    mark_stage_complete(pair_cache, 'relax', str(relax_dir))
                    update_stage_metadata(pair_cache, 'relax', {'skipped_reason': 'all_failed'})

    else:
        logger.info("[7/8] Rosetta Relax: ✓ CACHED")

    # ──────────────────────────────────────────────────────────────
    # STAGE 8a: AF3 JSON Preparation (per-pair)
    # ──────────────────────────────────────────────────────────────
    af3_args = pair.get('_af3_args', {})
    skip_af3 = af3_args.get('skip_af3', False)

    if skip_af3:
        logger.info("[8/8] AF3 Preparation: SKIPPED (--skip-af3)")
    elif is_stage_complete(pair_cache, 'af3_binary_prep') and is_stage_complete(pair_cache, 'af3_ternary_prep'):
        logger.info("[8/8] AF3 Preparation: ✓ CACHED")
    else:
        logger.info("[8/8] AF3 JSON Preparation")
        from prepare_af3_ml import thread_mutations_to_sequence, create_af3_json

        wt_sequence = af3_args.get('wt_sequence', '')
        binary_template = af3_args.get('af3_binary_template', '')
        ternary_template = af3_args.get('af3_ternary_template', '')

        if not wt_sequence or not binary_template or not ternary_template:
            logger.warning("  AF3 templates/WT sequence not configured; skipping AF3 prep")
        else:
            variant_sig_str = str(pair.get('variant_signature', ''))
            ligand_smiles = pair.get('ligand_smiles', '')

            mutant_seq = thread_mutations_to_sequence(wt_sequence, variant_sig_str)

            for mode, template_path in [('binary', binary_template), ('ternary', ternary_template)]:
                stage_name = f'af3_{mode}_prep'
                if is_stage_complete(pair_cache, stage_name):
                    logger.info(f"  AF3 {mode} prep: ✓ CACHED")
                    continue

                output_json = pair_cache / f'af3_{mode}' / f'{pair_id}_{mode}.json'
                success = create_af3_json(
                    template_path=template_path,
                    mutant_sequence=mutant_seq,
                    ligand_smiles=ligand_smiles,
                    pair_id=pair_id,
                    output_path=str(output_json),
                    mode=mode,
                )
                if success:
                    mark_stage_complete(pair_cache, stage_name, str(output_json))
                    logger.info(f"  ✓ AF3 {mode} JSON: {output_json.name}")
                else:
                    logger.error(f"  ✗ Failed to create AF3 {mode} JSON")

    logger.info(f"✓ Pair {pair_id} processed successfully\n")

    return {'status': 'SUCCESS'}



# ═══════════════════════════════════════════════════════════════════
# AF3 POST-LOOP FUNCTIONS (Stage 8b)
# ═══════════════════════════════════════════════════════════════════

def collect_pending_af3_jsons(cache_dir: Path) -> Dict[str, List[Path]]:
    """
    Scan all pair caches for prepared AF3 JSONs that don't yet have output.

    Returns:
        {'binary': [path1, ...], 'ternary': [path1, ...]}
    """
    pending = {'binary': [], 'ternary': []}

    for pair_dir in sorted(cache_dir.iterdir()):
        if not pair_dir.is_dir() or pair_dir.name.startswith('af3_staging'):
            continue

        for mode in ('binary', 'ternary'):
            af3_dir = pair_dir / f'af3_{mode}'
            json_files = list(af3_dir.glob(f'*_{mode}.json'))
            summary_file = af3_dir / 'summary.json'

            # Pending = has input JSON but no summary.json yet
            if json_files and not summary_file.exists():
                pending[mode].extend(json_files)

    return pending


def analyze_af3_outputs(cache_dir: Path, af3_staging_dir: Path, af3_args: Dict) -> int:
    """
    Re-entry: scan AF3 output for completed predictions, extract metrics,
    compute ligand RMSD, write summary.json per pair.

    Returns:
        Number of pairs analyzed
    """
    from prepare_af3_ml import extract_af3_metrics, compute_all_ligand_rmsds_to_rosetta, \
        find_all_relaxed_pdbs, find_best_relaxed_pdb, write_summary_json, \
        compute_binary_ternary_ligand_rmsd

    analyzed = 0

    for pair_dir in sorted(cache_dir.iterdir()):
        if not pair_dir.is_dir() or pair_dir.name.startswith('af3_staging'):
            continue
        pair_id = pair_dir.name

        # Collect per-mode results for this pair: mode -> (metrics, ligand_rmsds, cif_path)
        mode_results = {}

        for mode in ('binary', 'ternary'):
            summary_file = pair_dir / f'af3_{mode}' / 'summary.json'
            if summary_file.exists():
                continue  # Already analyzed

            stage_name = f'af3_{mode}'
            if is_stage_complete(pair_dir, stage_name):
                continue

            # Check if AF3 output exists
            af3_output_dir = af3_staging_dir / f'{mode}_output'
            name = f"{pair_id}_{mode}"

            # Try both flat and nested output layouts
            summary_conf = af3_output_dir / f"{name}_summary_confidences.json"
            if not summary_conf.exists():
                summary_conf = af3_output_dir / name / f"{name}_summary_confidences.json"
            if not summary_conf.exists():
                continue  # Not yet completed

            logger.info(f"  Analyzing AF3 {mode} for {pair_id}...")
            metrics = extract_af3_metrics(
                af3_output_dir=str(af3_output_dir),
                pair_id=pair_id,
                mode=mode,
            )

            if metrics is None:
                logger.warning(f"  Failed to extract metrics for {pair_id} {mode}")
                continue

            # Locate CIF file
            cif_path = af3_output_dir / f"{name}_model.cif"
            if not cif_path.exists():
                cif_path = af3_output_dir / name / f"{name}_model.cif"

            # Compute ligand RMSD (min across all + best-dG structure)
            ligand_rmsds = {'min': None, 'best_dG': None}
            relaxed_pdbs = find_all_relaxed_pdbs(str(pair_dir))
            best_dg_pdb = find_best_relaxed_pdb(str(pair_dir))
            if relaxed_pdbs and cif_path.exists():
                ligand_rmsds = compute_all_ligand_rmsds_to_rosetta(
                    af3_cif_path=str(cif_path),
                    relaxed_pdbs=relaxed_pdbs,
                    best_dg_pdb=best_dg_pdb,
                )

            mode_results[mode] = (metrics, ligand_rmsds, cif_path if cif_path.exists() else None)

        if not mode_results:
            continue

        # Compute binary-to-ternary ligand RMSD if both modes available
        bt_rmsd = None
        if 'binary' in mode_results and 'ternary' in mode_results:
            b_cif = mode_results['binary'][2]
            t_cif = mode_results['ternary'][2]
            if b_cif and t_cif:
                bt_rmsd = compute_binary_ternary_ligand_rmsd(
                    binary_cif_path=str(b_cif),
                    ternary_cif_path=str(t_cif),
                )
                if bt_rmsd is not None:
                    logger.info(f"  {pair_id} binary-to-ternary ligand RMSD: {bt_rmsd:.3f} A")

        # Write summary.json for each newly-analyzed mode
        for mode, (metrics, ligand_rmsds, _) in mode_results.items():
            write_summary_json(str(pair_dir / f'af3_{mode}'), metrics, ligand_rmsds, bt_rmsd)
            mark_stage_complete(pair_dir, f'af3_{mode}', str(pair_dir / f'af3_{mode}'))
            analyzed += 1

            logger.info(f"  ✓ {pair_id} {mode}: ipTM={metrics.get('ipTM')}, "
                         f"RMSD_min={ligand_rmsds.get('min')}, "
                         f"RMSD_bestdG={ligand_rmsds.get('best_dG')}")

    return analyzed


def batch_and_submit_af3(
    cache_dir: Path,
    af3_staging_dir: Path,
    af3_args: Dict
) -> Dict[str, Optional[str]]:
    """
    Batch remaining AF3 JSONs into subdirectories and submit SLURM GPU array jobs.

    Returns:
        {'binary': job_id_or_None, 'ternary': job_id_or_None}
    """
    import shutil

    pending = collect_pending_af3_jsons(cache_dir)
    job_ids = {}

    for mode in ('binary', 'ternary'):
        json_files = pending[mode]
        if not json_files:
            logger.info(f"  No pending AF3 {mode} JSONs to submit")
            job_ids[mode] = None
            continue

        logger.info(f"  Batching {len(json_files)} {mode} AF3 JSONs...")

        # Create staging directory
        staging_dir = af3_staging_dir / mode
        staging_dir.mkdir(parents=True, exist_ok=True)

        # Copy JSONs to staging (flat)
        for jf in json_files:
            shutil.copy2(str(jf), str(staging_dir / jf.name))

        # Batch into subdirectories
        batch_size = af3_args.get('af3_batch_size', 60)
        all_jsons = sorted(staging_dir.glob('*.json'))
        num_batches = max(1, -(-len(all_jsons) // batch_size))  # ceil division

        for batch_idx in range(num_batches):
            batch_name = f"batch_{batch_idx + 1:02d}"
            batch_dir = staging_dir / batch_name
            batch_dir.mkdir(parents=True, exist_ok=True)

            start = batch_idx * batch_size
            end = min(start + batch_size, len(all_jsons))
            for jf in all_jsons[start:end]:
                shutil.move(str(jf), str(batch_dir / jf.name))

        logger.info(f"  Created {num_batches} batches of up to {batch_size} JSONs")

        # Create output directory
        output_dir = af3_staging_dir / f'{mode}_output'
        output_dir.mkdir(parents=True, exist_ok=True)
        log_dir = output_dir / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        # Generate SLURM script
        partition = af3_args.get('af3_partition', 'aa100')
        time_limit = af3_args.get('af3_time_limit', '01:00:00')
        model_params = af3_args.get('af3_model_params_dir', '/projects/ryde3462/software/af3')
        job_name = f"af3_ml_{mode}"

        script_content = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time={time_limit}
#SBATCH --partition={partition}
#SBATCH --qos=normal
#SBATCH --gres=gpu:1
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/af3_{mode}_%A_%a.out
#SBATCH --ntasks=10
#SBATCH --account=ucb472_asc2
#SBATCH --array=1-{num_batches}

module purge
module load alphafold/3.0.0

BASE_INPUT_DIR="{staging_dir}"
BASE_OUTPUT_DIR="{output_dir}"
export AF3_MODEL_PARAMETERS_DIR="{model_params}"

BATCH_ID=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
BATCH_FOLDER="batch_${{BATCH_ID}}"

export INPUT_DIR="${{BASE_INPUT_DIR}}/${{BATCH_FOLDER}}"
export OUTPUT_DIR="${{BASE_OUTPUT_DIR}}"

echo "Job Array ID: $SLURM_ARRAY_TASK_ID"
echo "Processing Batch: $BATCH_FOLDER"
echo "Input Dir:  $INPUT_DIR"
echo "Output Dir: $OUTPUT_DIR"

mkdir -p $OUTPUT_DIR

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist!"
    exit 1
fi

run_alphafold \\
  --input_dir=$INPUT_DIR \\
  --output_dir=$OUTPUT_DIR \\
  --model_dir=$AF3_MODEL_PARAMETERS_DIR \\
  --run_data_pipeline=false
"""

        script_path = af3_staging_dir / f'submit_af3_{mode}.sh'
        with open(script_path, 'w', newline='\n') as f:
            f.write(script_content)

        logger.info(f"  Generated SLURM script: {script_path}")

        # Submit
        try:
            result = subprocess.run(
                ['sbatch', str(script_path)],
                capture_output=True, text=True, check=True
            )
            job_id = result.stdout.strip().split()[-1]
            logger.info(f"  ✓ AF3 {mode} submitted: job {job_id} ({num_batches} array tasks)")
            job_ids[mode] = job_id
        except subprocess.CalledProcessError as e:
            logger.error(f"  ✗ AF3 {mode} submission failed: {e.stderr}")
            job_ids[mode] = None
        except FileNotFoundError:
            logger.error("  ✗ sbatch not found (not on SLURM cluster?)")
            job_ids[mode] = None

    return job_ids


def main():
    # Import WT sequence default from prepare_af3_ml
    from prepare_af3_ml import WT_PYR1_SEQUENCE

    parser = argparse.ArgumentParser(
        description='Orchestrate ML dataset generation pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--pairs-csv', required=True, help='CSV with pair_id, ligand_smiles, variant_signature, label')
    parser.add_argument('--cache-dir', required=True, help='Cache directory for outputs')
    parser.add_argument('--template-pdb', required=True, help='WT PYR1 template PDB (NO ligand, for threading)')
    parser.add_argument('--reference-pdb', required=True, help='Reference PDB WITH template ligand for alignment (e.g., 3QN1_H2O.pdb)')
    parser.add_argument('--docking-repeats', type=int, default=50, help='Docking repeats per conformer')
    parser.add_argument('--docking-arrays', type=int, default=10, help='Number of SLURM array tasks for docking (default: 10)')
    parser.add_argument('--use-slurm', action='store_true', help='Submit jobs to SLURM')
    parser.add_argument('--max-pairs', type=int, help='Limit number of pairs (for testing)')
    parser.add_argument('--reference-chain', default='X', help='Chain of template ligand in reference PDB (default: X)')
    parser.add_argument('--reference-residue', type=int, default=1, help='Residue number of template ligand (default: 1)')

    # AF3 arguments
    af3_group = parser.add_argument_group('AF3 options')
    af3_group.add_argument('--wt-sequence', default=WT_PYR1_SEQUENCE,
                           help='WT PYR1 protein sequence (181 residues)')
    af3_group.add_argument('--af3-binary-template',
                           default=str(PROJECT_ROOT / 'design' / 'templates' / 'pyr1_binary_template.json'),
                           help='Path to AF3 binary template JSON')
    af3_group.add_argument('--af3-ternary-template',
                           default=str(PROJECT_ROOT / 'design' / 'templates' / 'pyr1_ternary_template.json'),
                           help='Path to AF3 ternary template JSON')
    af3_group.add_argument('--af3-batch-size', type=int, default=60,
                           help='Number of AF3 JSONs per SLURM batch (default: 60)')
    af3_group.add_argument('--af3-partition', default='aa100',
                           help='GPU partition for AF3 jobs (default: aa100)')
    af3_group.add_argument('--af3-time-limit', default='01:00:00',
                           help='Walltime per AF3 array task (default: 01:00:00)')
    af3_group.add_argument('--af3-model-params-dir', default='/projects/ryde3462/software/af3',
                           help='Path to AF3 model parameters')
    af3_group.add_argument('--skip-af3', action='store_true',
                           help='Skip AF3 preparation and submission')

    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Build AF3 args dict to pass through
    af3_args = {
        'wt_sequence': args.wt_sequence,
        'af3_binary_template': args.af3_binary_template,
        'af3_ternary_template': args.af3_ternary_template,
        'af3_batch_size': args.af3_batch_size,
        'af3_partition': args.af3_partition,
        'af3_time_limit': args.af3_time_limit,
        'af3_model_params_dir': args.af3_model_params_dir,
        'skip_af3': args.skip_af3,
    }

    # Load pairs
    logger.info(f"Loading pairs from {args.pairs_csv}...")
    pairs_df = pd.read_csv(args.pairs_csv)

    if args.max_pairs:
        pairs_df = pairs_df.head(args.max_pairs)

    logger.info(f"Processing {len(pairs_df)} pairs")

    # Process each pair (Stages 1-8a)
    results = []

    for idx, row in pairs_df.iterrows():
        result = process_single_pair(
            pair=row.to_dict(),
            cache_dir=cache_dir,
            template_pdb=args.template_pdb,
            reference_pdb=args.reference_pdb,
            use_slurm=args.use_slurm,
            docking_repeats=args.docking_repeats,
            docking_arrays=args.docking_arrays,
            reference_chain=args.reference_chain,
            reference_residue=args.reference_residue,
            af3_args=af3_args,
        )

        results.append({
            'pair_id': row['pair_id'],
            **result
        })

    # ══════════════════════════════════════════════════════════════
    # STAGE 8b: AF3 Batch Submit + Analyze (post-loop)
    # ══════════════════════════════════════════════════════════════
    if not args.skip_af3 and args.use_slurm:
        af3_staging_dir = cache_dir / 'af3_staging'
        af3_staging_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"\n{'='*60}")
        logger.info("STAGE 8b: AF3 Batch Submit & Analyze")
        logger.info(f"{'='*60}")

        # Re-entry: analyze any completed AF3 output first
        n_analyzed = analyze_af3_outputs(cache_dir, af3_staging_dir, af3_args)
        if n_analyzed > 0:
            logger.info(f"  Analyzed {n_analyzed} AF3 predictions from previous runs")

        # Collect and submit remaining
        pending = collect_pending_af3_jsons(cache_dir)
        total_pending = sum(len(v) for v in pending.values())

        if total_pending > 0:
            logger.info(f"  {total_pending} AF3 predictions pending submission")
            job_ids = batch_and_submit_af3(cache_dir, af3_staging_dir, af3_args)

            for mode, jid in job_ids.items():
                if jid:
                    logger.info(f"  AF3 {mode} job: {jid}")
            logger.info("  → Re-run orchestrator after AF3 jobs complete to analyze output")
        else:
            logger.info("  No pending AF3 submissions (all complete or not yet prepared)")

    # Save results summary
    results_df = pd.DataFrame(results)
    results_csv = cache_dir / 'processing_summary.csv'
    results_df.to_csv(results_csv, index=False)

    logger.info(f"\n{'='*60}")
    logger.info("PROCESSING SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"Total pairs: {len(results_df)}")
    logger.info(f"Successful: {(results_df['status'] == 'SUCCESS').sum()}")
    logger.info(f"Failed: {(results_df['status'] == 'FAILED').sum()}")
    logger.info(f"Summary saved: {results_csv}")
    logger.info(f"{'='*60}\n")


if __name__ == '__main__':
    main()
