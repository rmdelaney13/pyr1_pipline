#!/usr/bin/env python
"""
Master orchestrator for PYR1 ML dataset generation.

CORRECT WORKFLOW (Direct-to-Mutant Docking):
1. Generate ligand conformers (RDKit ETKDG + MMFF)
2. Thread mutations onto WT PYR1 → mutant.pdb
3. Dock conformers to MUTANT pocket (NO glycine shaving!)
4. Cluster poses → extract statistics (convergence, clashes)
5. Relax best pose (skip if severe clash)
6. Run AF3 binary + ternary predictions
7. Aggregate all features into ML table

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
        'pyr1_pipeline/scripts/thread_variant_to_pdb.py',
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
    output_dir: Path,
    docking_repeats: int = 50,
    use_slurm: bool = False,
    array_tasks: int = 10
) -> Optional[str]:
    """
    Run docking to mutant pocket.

    Args:
        mutant_pdb: Pre-threaded mutant structure
        conformers_sdf: Ligand conformers SDF
        output_dir: Output directory
        docking_repeats: Number of docking attempts per conformer
        use_slurm: Submit as SLURM array job
        array_tasks: Number of SLURM array tasks

    Returns:
        Job ID if SLURM, 'local' if local, None if failed
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  Docking to mutant pocket ({docking_repeats} repeats)...")

    # Create config file
    config_path = output_dir / 'docking_config.txt'
    with open(config_path, 'w') as f:
        f.write(f"""[mutant_docking]
MutantPDB = {mutant_pdb}
LigandSDF = {conformers_sdf}
OutputDir = {output_dir}
DockingRepeats = {docking_repeats}
ArrayTaskCount = {array_tasks}
""")

    if use_slurm:
        # Submit SLURM array job
        cmd = [
            'sbatch',
            '--job-name', f"dock_{mutant_pdb.stem}",
            '--output', str(output_dir / 'docking_%a.log'),
            '--array', f'0-{array_tasks-1}',
            '--time', '04:00:00',
            '--cpus-per-task', '4',
            '--mem', '8G',
            'pyr1_pipeline/docking/scripts/submit_docking_mutant.sh',  # TODO: Create this wrapper
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
            'pyr1_pipeline/docking/scripts/grade_conformers_mutant_docking.py',
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
        'pyr1_pipeline/docking/scripts/cluster_docked_with_stats.py',
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
    use_slurm: bool = False,
    docking_repeats: int = 50
) -> Dict:
    """
    Process a single (ligand, variant) pair through the full pipeline.

    Args:
        pair: Dict with pair_id, ligand_name, ligand_smiles, variant_signature, etc.
        cache_dir: Base cache directory
        template_pdb: WT PYR1 template PDB path
        use_slurm: Use SLURM for parallelization
        docking_repeats: Number of docking repeats

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
        logger.info("[1/6] Conformer Generation")
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
        logger.info("[1/6] Conformer Generation: ✓ CACHED")

    conformers_sdf = pair_cache / 'conformers' / 'conformers_final.sdf'

    # ──────────────────────────────────────────────────────────────
    # STAGE 2: Mutation Threading
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'threading'):
        logger.info("[2/6] Mutation Threading")
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
        logger.info("[2/6] Mutation Threading: ✓ CACHED")

    mutant_pdb = pair_cache / 'mutant.pdb'

    # ──────────────────────────────────────────────────────────────
    # STAGE 3: Docking to Mutant Pocket
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'docking'):
        logger.info("[3/6] Docking to Mutant")
        docking_dir = pair_cache / 'docking'

        job_id = run_docking(
            mutant_pdb=mutant_pdb,
            conformers_sdf=conformers_sdf,
            output_dir=docking_dir,
            docking_repeats=docking_repeats,
            use_slurm=use_slurm
        )

        if job_id:
            mark_stage_complete(pair_cache, 'docking', str(docking_dir))
        else:
            return {'status': 'FAILED', 'stage': 'docking'}

    else:
        logger.info("[3/6] Docking to Mutant: ✓ CACHED")

    docking_dir = pair_cache / 'docking'

    # ──────────────────────────────────────────────────────────────
    # STAGE 4: Clustering with Statistics
    # ──────────────────────────────────────────────────────────────
    if not is_stage_complete(pair_cache, 'clustering'):
        logger.info("[4/6] Clustering & Statistics")
        cluster_dir = pair_cache / 'clustered'

        success = run_clustering(
            docking_output_dir=docking_dir,
            cluster_output_dir=cluster_dir,
            rmsd_cutoff=2.0
        )

        if success:
            mark_stage_complete(pair_cache, 'clustering', str(cluster_dir))
        else:
            return {'status': 'FAILED', 'stage': 'clustering'}

    else:
        logger.info("[4/6] Clustering & Statistics: ✓ CACHED")

    cluster_dir = pair_cache / 'clustered'

    # ──────────────────────────────────────────────────────────────
    # STAGE 5: Relax (skip if severe clash)
    # ──────────────────────────────────────────────────────────────
    if should_skip_relax(cluster_dir):
        logger.info("[5/6] Rosetta Relax: SKIPPED (severe clash)")
        mark_stage_complete(pair_cache, 'relax', None)

    elif not is_stage_complete(pair_cache, 'relax'):
        logger.info("[5/6] Rosetta Relax")
        # TODO: Implement relax call
        logger.warning("  (Relax not yet implemented in orchestrator)")

    else:
        logger.info("[5/6] Rosetta Relax: ✓ CACHED")

    # ──────────────────────────────────────────────────────────────
    # STAGE 6: AF3 Predictions
    # ──────────────────────────────────────────────────────────────
    logger.info("[6/6] AF3 Predictions")
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
    parser.add_argument('--template-pdb', required=True, help='WT PYR1 template PDB')
    parser.add_argument('--docking-repeats', type=int, default=50, help='Docking repeats per conformer')
    parser.add_argument('--use-slurm', action='store_true', help='Submit jobs to SLURM')
    parser.add_argument('--max-pairs', type=int, help='Limit number of pairs (for testing)')

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
            use_slurm=args.use_slurm,
            docking_repeats=args.docking_repeats
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
