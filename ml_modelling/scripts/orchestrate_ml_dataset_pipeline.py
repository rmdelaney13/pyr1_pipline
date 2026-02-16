#!/usr/bin/env python
"""
Master orchestrator for PYR1 ML dataset generation.

CORRECT WORKFLOW (Direct-to-Mutant Docking):
1. Generate ligand conformers (RDKit ETKDG + MMFF)
2. Create alignment table (identify H-bond acceptor atoms for docking)
3. Thread mutations onto WT PYR1 → mutant.pdb
4. Dock conformers to MUTANT pocket (NO glycine shaving!)
5. Cluster poses → extract statistics (convergence, clashes)
6. Relax best pose (skip if severe clash)
7. Run AF3 binary + ternary predictions
8. Aggregate all features into ML table

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
        'docking/scripts/create_table.py',
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
        'scripts/thread_variant_to_pdb.py',
        '--template', template_pdb,
        '--signature', str(variant_signature),
        '--output', str(output_pdb),
        '--chain', chain
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  ✓ Mutant structure created: {output_pdb}")
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Threading failed: {e.stderr}")
        return False


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
    params_file = "docking/ligand_alignment/files_for_PYR1_docking/A8T.params"

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
ClusterRMSDCutoff = 0.75
""")

    if use_slurm:
        # Estimate walltime based on docking repeats
        # ~2-5 seconds per docking attempt, so:
        # 10 repeats: 30 min, 50 repeats: 2 hours, 100 repeats: 4 hours
        if docking_repeats <= 10:
            walltime = '00:30:00'
        elif docking_repeats <= 50:
            walltime = '02:00:00'
        elif docking_repeats <= 100:
            walltime = '04:00:00'
        else:
            walltime = '08:00:00'

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
            'docking/scripts/submit_docking_mutant.sh',
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
            'docking/scripts/grade_conformers_mutant_docking.py',
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
    cluster_output_dir: Path,
    rmsd_cutoff: float = 2.0
) -> bool:
    """
    Cluster docked poses and extract statistics.

    Args:
        docking_output_dir: Directory with docking results
        cluster_output_dir: Output directory for clustering
        rmsd_cutoff: RMSD clustering cutoff (Å)

    Returns:
        True if successful, False otherwise
    """
    cluster_output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  Clustering docked poses (RMSD cutoff {rmsd_cutoff} Å)...")

    cmd = [
        'python',
        'docking/scripts/cluster_docked_with_stats.py',
        '--input-dir', str(docking_output_dir),
        '--output-dir', str(cluster_output_dir),
        '--rmsd-cutoff', str(rmsd_cutoff),
        '--stats-csv', 'clustering_stats.csv'
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"  ✓ Clustering complete")

        # Check cluster statistics
        stats_json = cluster_output_dir / 'clustering_stats.json'
        if stats_json.exists():
            with open(stats_json, 'r') as f:
                stats = json.load(f)
                global_stats = stats.get('global_stats', {})
                logger.info(f"    Convergence: {global_stats.get('convergence_ratio', 0):.2%}")
                logger.info(f"    Best score: {global_stats.get('best_overall_score', np.nan):.2f}")

                # Warn on clashes
                if global_stats.get('clash_count', 0) > 0:
                    logger.warning(f"    ⚠ {global_stats['clash_count']} clusters with clashes!")

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"  ✗ Clustering failed: {e.stderr}")
        return False


def should_skip_relax(cluster_output_dir: Path) -> bool:
    """
    Check if relax should be skipped (severe clash detected).

    Args:
        cluster_output_dir: Clustering output directory

    Returns:
        True if should skip, False otherwise
    """
    stats_json = cluster_output_dir / 'clustering_stats.json'

    if not stats_json.exists():
        return False

    try:
        with open(stats_json, 'r') as f:
            stats = json.load(f)
            best_score = stats.get('global_stats', {}).get('best_overall_score', 0)

            # Skip if severe clash (score > 5)
            if best_score > 5:
                logger.info(f"  Skipping relax: severe clash (score = {best_score:.2f})")
                return True

            return False

    except Exception:
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
    reference_chain: str = 'X',
    reference_residue: int = 1
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
        reference_chain: Chain of template ligand in reference PDB (default: 'X')
        reference_residue: Residue number of template ligand (default: 1)

    Returns:
        Dict with status and any errors
    """
    pair_id = pair['pair_id']
    pair_cache = cache_dir / pair_id

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
        logger.info("[1/7] Conformer Generation")
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
        logger.info("[1/7] Conformer Generation: ✓ CACHED")

    conformers_sdf = pair_cache / 'conformers' / 'conformers_final.sdf'

    # ──────────────────────────────────────────────────────────────
    # STAGE 2: Alignment Table Generation
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'alignment_table'):
        logger.info("[2/7] Alignment Table Generation")
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
        logger.info("[2/7] Alignment Table Generation: ✓ CACHED")

    alignment_csv = pair_cache / 'alignment_table.csv'
    conformers_params_dir = pair_cache / 'conformers_params'

    # ──────────────────────────────────────────────────────────────
    # STAGE 3: Mutation Threading
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'threading'):
        logger.info("[3/7] Mutation Threading")
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
        logger.info("[3/7] Mutation Threading: ✓ CACHED")

    mutant_pdb = pair_cache / 'mutant.pdb'

    # ──────────────────────────────────────────────────────────────
    # STAGE 4: Docking to Mutant Pocket
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'docking'):
        logger.info("[4/7] Docking to Mutant")
        docking_dir = pair_cache / 'docking'

        job_id = run_docking(
            mutant_pdb=mutant_pdb,
            conformers_sdf=conformers_sdf,
            alignment_csv=alignment_csv,
            conformers_dir=conformers_params_dir,
            output_dir=docking_dir,
            reference_pdb=reference_pdb,
            docking_repeats=docking_repeats,
            use_slurm=use_slurm,
            reference_chain=reference_chain,
            reference_residue=reference_residue
        )

        if job_id:
            if use_slurm:
                # For SLURM, job is submitted but not finished
                # Don't mark as complete yet - skip to next pair
                logger.info(f"  ⏳ SLURM job submitted: {job_id}")
                logger.info(f"  → Clustering will run after jobs finish")
                return {'status': 'SLURM_SUBMITTED', 'job_id': job_id, 'stage': 'docking'}
            else:
                mark_stage_complete(pair_cache, 'docking', str(docking_dir))
        else:
            return {'status': 'FAILED', 'stage': 'docking'}

    else:
        logger.info("[4/7] Docking to Mutant: ✓ CACHED")

    docking_dir = pair_cache / 'docking'

    # ──────────────────────────────────────────────────────────────
    # STAGE 5: Clustering with Statistics
    # ──────────────────────────────────────────────────────────────
    # Skip external clustering for local runs (mutant docking does in-loop clustering)
    # For SLURM runs, clustering happens after all array tasks complete
    if use_slurm:
        logger.info("[5/7] Clustering & Statistics: SKIPPED (SLURM - run after jobs complete)")
        logger.info("  → After SLURM jobs finish, run: cluster_docked_post_array.py")
        # Don't mark as complete - user needs to run clustering manually
    elif not is_stage_complete(pair_cache, 'clustering'):
        logger.info("[5/7] Clustering & Statistics: ✓ DONE (in-loop)")
        # For local runs, clustering was done during docking (EnablePoseClusteringInArrayTask=True)
        # Mark as complete automatically
        mark_stage_complete(pair_cache, 'clustering', str(docking_dir))
    else:
        logger.info("[5/7] Clustering & Statistics: ✓ CACHED")

    cluster_dir = pair_cache / 'docking'  # Use docking dir since clustering is in-loop

    # ──────────────────────────────────────────────────────────────
    # STAGE 6: Relax (skip if severe clash)
    # ──────────────────────────────────────────────────────────────
    if should_skip_relax(cluster_dir):
        logger.info("[6/7] Rosetta Relax: SKIPPED (severe clash)")
        mark_stage_complete(pair_cache, 'relax', None)

    elif not is_stage_complete(pair_cache, 'relax'):
        logger.info("[6/7] Rosetta Relax")
        # TODO: Implement relax call
        logger.warning("  (Relax not yet implemented in orchestrator)")

    else:
        logger.info("[6/7] Rosetta Relax: ✓ CACHED")

    # ──────────────────────────────────────────────────────────────
    # STAGE 7: AF3 Predictions
    # ──────────────────────────────────────────────────────────────
    logger.info("[7/7] AF3 Predictions")
    logger.warning("  (AF3 not yet implemented in orchestrator)")

    logger.info(f"✓ Pair {pair_id} processed successfully\n")

    return {'status': 'SUCCESS'}


def main():
    parser = argparse.ArgumentParser(
        description='Orchestrate ML dataset generation pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--pairs-csv', required=True, help='CSV with pair_id, ligand_smiles, variant_signature, label')
    parser.add_argument('--cache-dir', required=True, help='Cache directory for outputs')
    parser.add_argument('--template-pdb', required=True, help='WT PYR1 template PDB (NO ligand, for threading)')
    parser.add_argument('--reference-pdb', required=True, help='Reference PDB WITH template ligand for alignment (e.g., 3QN1_H2O.pdb)')
    parser.add_argument('--docking-repeats', type=int, default=50, help='Docking repeats per conformer')
    parser.add_argument('--use-slurm', action='store_true', help='Submit jobs to SLURM')
    parser.add_argument('--max-pairs', type=int, help='Limit number of pairs (for testing)')
    parser.add_argument('--reference-chain', default='X', help='Chain of template ligand in reference PDB (default: X)')
    parser.add_argument('--reference-residue', type=int, default=1, help='Residue number of template ligand (default: 1)')

    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Load pairs
    logger.info(f"Loading pairs from {args.pairs_csv}...")
    pairs_df = pd.read_csv(args.pairs_csv)

    if args.max_pairs:
        pairs_df = pairs_df.head(args.max_pairs)

    logger.info(f"Processing {len(pairs_df)} pairs")

    # Process each pair
    results = []

    for idx, row in pairs_df.iterrows():
        result = process_single_pair(
            pair=row.to_dict(),
            cache_dir=cache_dir,
            template_pdb=args.template_pdb,
            reference_pdb=args.reference_pdb,
            use_slurm=args.use_slurm,
            docking_repeats=args.docking_repeats,
            reference_chain=args.reference_chain,
            reference_residue=args.reference_residue
        )

        results.append({
            'pair_id': row['pair_id'],
            **result
        })

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
