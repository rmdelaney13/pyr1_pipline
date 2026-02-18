#!/usr/bin/env python
"""
Aggregate all Rosetta + AF3 outputs into ML feature table.

This script extracts features from:
  - Conformer generation (diversity, energy)
  - Docking (cluster statistics, best score, convergence)
  - Relax (interface energy, buried unsats, H-bonds)
  - AF3 binary (ipTM, pLDDT, interface PAE)
  - AF3 ternary (ipTM, pLDDT, water H-bonds)

Output: features_table.csv with ~45 columns per (ligand, variant) pair

Author: Claude Code (Whitehead Lab PYR1 Pipeline)
Date: 2026-02-16
"""

import argparse
import json
import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════
# FEATURE EXTRACTION FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

def extract_conformer_features(pair_cache: Path) -> Dict:
    """
    Extract ligand conformer generation features.

    Expected files:
        pair_cache/conformers/conformer_report.csv
        pair_cache/conformers/metadata.json

    Returns dict with:
        - conformer_count: Number of final conformers
        - conformer_min_energy: Lowest MMFF energy
        - conformer_max_rmsd: Maximum pairwise RMSD (diversity)
    """
    conformer_dir = pair_cache / 'conformers'
    report_csv = conformer_dir / 'conformer_report.csv'
    metadata_json = conformer_dir / 'metadata.json'

    if not report_csv.exists():
        return {
            'conformer_status': 'MISSING',
            'conformer_count': np.nan,
            'conformer_min_energy': np.nan,
            'conformer_max_rmsd': np.nan,
        }

    try:
        df = pd.read_csv(report_csv)

        return {
            'conformer_status': 'COMPLETE',
            'conformer_count': len(df),
            'conformer_min_energy': df['mmff_energy'].min() if 'mmff_energy' in df.columns else np.nan,
            'conformer_max_rmsd': df['rmsd_to_centroid'].max() if 'rmsd_to_centroid' in df.columns else np.nan,
        }
    except Exception as e:
        logger.warning(f"Failed to parse conformer report: {e}")
        return {
            'conformer_status': 'ERROR',
            'conformer_count': np.nan,
            'conformer_min_energy': np.nan,
            'conformer_max_rmsd': np.nan,
        }


def extract_docking_features(pair_cache: Path) -> Dict:
    """
    Extract docking cluster statistics (ML-critical features!).

    Primary source: pair_cache/docking/clustering_stats.json (SLURM post-clustering)
    Fallback source: pair_cache/metadata.json stages.docking (local runs with in-loop clustering)

    Returns dict with:
        - docking_best_score: Best Rosetta score across all poses
        - docking_cluster_1_size: Size of largest cluster
        - docking_convergence_ratio: Fraction in largest cluster (KEY METRIC)
        - docking_num_clusters: Total number of clusters
        - docking_score_range: Energy landscape breadth
        - docking_clash_flag: 1 if best_score > 0 (clash)
        - docking_cluster_1_rmsd: Intra-cluster RMSD (tightness)
        - docking_pass_rate: Fraction of docking attempts that passed score threshold
        - docking_total_attempts: Total docking attempts
        - docking_mean_quality: Mean H-bond quality among passing poses
    """
    nan_result = {
        'docking_status': 'MISSING',
        'docking_best_score': np.nan,
        'docking_cluster_1_size': np.nan,
        'docking_convergence_ratio': np.nan,
        'docking_num_clusters': np.nan,
        'docking_score_range': np.nan,
        'docking_clash_flag': np.nan,
        'docking_cluster_1_rmsd': np.nan,
        'docking_pass_rate': np.nan,
        'docking_total_attempts': np.nan,
        'docking_mean_quality': np.nan,
    }

    docking_dir = pair_cache / 'docking'
    stats_json = docking_dir / 'clustering_stats.json'

    features = dict(nan_result)

    # Source 1: clustering_stats.json (SLURM post-clustering)
    if stats_json.exists():
        try:
            with open(stats_json, 'r') as f:
                data = json.load(f)

            global_stats = data.get('global_stats', {})
            clusters = data.get('clusters', [])

            best_score = global_stats.get('best_overall_score', np.nan)
            features.update({
                'docking_status': 'COMPLETE',
                'docking_best_score': best_score,
                'docking_convergence_ratio': global_stats.get('convergence_ratio', np.nan),
                'docking_num_clusters': global_stats.get('num_clusters', np.nan),
                'docking_score_range': global_stats.get('score_range', np.nan),
                'docking_clash_flag': 1 if isinstance(best_score, (int, float)) and best_score > 0 else 0,
            })

            if clusters:
                top_cluster = clusters[0]
                features['docking_cluster_1_size'] = top_cluster.get('size', np.nan)
                features['docking_cluster_1_rmsd'] = top_cluster.get('intra_cluster_rmsd', np.nan)

        except Exception as e:
            logger.warning(f"Failed to parse clustering_stats.json: {e}")

    # Source 2: metadata.json stages.docking (from extract_docking_stats in orchestrator)
    metadata_path = pair_cache / 'metadata.json'
    if metadata_path.exists():
        try:
            with open(metadata_path, 'r') as f:
                metadata = json.load(f)
            dock_meta = metadata.get('stages', {}).get('docking', {})

            if 'pass_rate' in dock_meta:
                features['docking_pass_rate'] = dock_meta.get('pass_rate', np.nan)
                features['docking_total_attempts'] = dock_meta.get('total_docking_attempts', np.nan)
                features['docking_mean_quality'] = dock_meta.get('mean_quality', np.nan)

                # Fill in best_score from metadata if not from clustering_stats
                if features['docking_status'] == 'MISSING':
                    features['docking_status'] = 'COMPLETE'
                    best = dock_meta.get('best_score')
                    if best is not None:
                        features['docking_best_score'] = best
                        features['docking_clash_flag'] = 1 if best > 0 else 0
                    features['docking_num_clusters'] = dock_meta.get('num_cluster_reps',
                                                                      dock_meta.get('saved_cluster_count', np.nan))

        except Exception as e:
            logger.warning(f"Failed to read docking metadata: {e}")

    return features


def extract_rosetta_relax_features(pair_cache: Path) -> Dict:
    """
    Extract distributional Rosetta relax features from metadata.json.

    The orchestrator aggregates scores from multiple relaxed structures
    (up to 20 cluster representatives) and stores distributional stats
    (best, mean, std, median) in metadata.json stages.relax.

    Returns dict with distributional features:
        - rosetta_dG_sep_best/mean/std
        - rosetta_buried_unsats_best/mean
        - rosetta_shape_comp_best/mean
        - rosetta_cms_best/mean (ContactMolecularSurface - works for ligands)
        - rosetta_dsasa_int_mean
        - rosetta_hbonds_to_ligand_mean
        - rosetta_charge_satisfied_frac
        - rosetta_nres_int_mean, rosetta_hbonds_int_mean
        - rosetta_n_structures_relaxed
    """
    nan_features = {
        'rosetta_status': 'MISSING',
        'rosetta_dG_sep_best': np.nan,
        'rosetta_dG_sep_mean': np.nan,
        'rosetta_dG_sep_std': np.nan,
        'rosetta_buried_unsats_best': np.nan,
        'rosetta_buried_unsats_mean': np.nan,
        'rosetta_shape_comp_best': np.nan,
        'rosetta_shape_comp_mean': np.nan,
        'rosetta_cms_best': np.nan,
        'rosetta_cms_mean': np.nan,
        'rosetta_dsasa_int_mean': np.nan,
        'rosetta_hbonds_to_ligand_mean': np.nan,
        'rosetta_charge_satisfied_frac': np.nan,
        'rosetta_nres_int_mean': np.nan,
        'rosetta_hbonds_int_mean': np.nan,
        'rosetta_n_structures_relaxed': np.nan,
    }

    metadata_path = pair_cache / 'metadata.json'
    if not metadata_path.exists():
        return nan_features

    try:
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        relax_meta = metadata.get('stages', {}).get('relax', {})

        # Check if relax was skipped
        if relax_meta.get('skipped_reason'):
            return {**nan_features, 'rosetta_status': 'SKIPPED'}

        if relax_meta.get('status') != 'complete':
            return nan_features

        # Map from aggregated key names to output feature names
        # The orchestrator stores keys like "dG_sep_best", "dG_sep_mean", etc.
        key_map = {
            'rosetta_dG_sep_best': 'dG_sep_best',
            'rosetta_dG_sep_mean': 'dG_sep_mean',
            'rosetta_dG_sep_std': 'dG_sep_std',
            'rosetta_buried_unsats_best': 'buried_unsatisfied_polars_best',
            'rosetta_buried_unsats_mean': 'buried_unsatisfied_polars_mean',
            'rosetta_shape_comp_best': 'shape_complementarity_best',
            'rosetta_shape_comp_mean': 'shape_complementarity_mean',
            'rosetta_cms_best': 'contact_molecular_surface_best',
            'rosetta_cms_mean': 'contact_molecular_surface_mean',
            'rosetta_dsasa_int_mean': 'dsasa_int_mean',
            'rosetta_hbonds_to_ligand_mean': 'total_hbonds_to_ligand_mean',
            'rosetta_charge_satisfied_frac': 'charge_satisfied_mean',
            'rosetta_nres_int_mean': 'nres_int_mean',
            'rosetta_hbonds_int_mean': 'hbonds_int_mean',
            'rosetta_n_structures_relaxed': 'n_structures_relaxed',
        }

        features = {'rosetta_status': 'COMPLETE'}
        for feature_name, meta_key in key_map.items():
            val = relax_meta.get(meta_key, np.nan)
            if val is None:
                val = np.nan
            features[feature_name] = val

        return features

    except Exception as e:
        logger.warning(f"Failed to parse relax metadata: {e}")
        return {**nan_features, 'rosetta_status': 'ERROR'}


def extract_af3_features(pair_cache: Path, mode: str = 'binary') -> Dict:
    """
    Extract AlphaFold3 prediction features.

    Expected files:
        pair_cache/af3_binary/summary.json
        pair_cache/af3_ternary/summary.json

    Args:
        pair_cache: Path to pair cache directory
        mode: 'binary' or 'ternary'

    Returns dict with:
        - af3_{mode}_ipTM: Interface predicted TM-score
        - af3_{mode}_pLDDT_protein: Mean protein pLDDT
        - af3_{mode}_pLDDT_ligand: Mean ligand pLDDT
        - af3_{mode}_interface_PAE: Mean interface PAE
        - af3_{mode}_ligand_RMSD_min: Min ligand RMSD across all Rosetta structures
        - af3_{mode}_ligand_RMSD_bestdG: Ligand RMSD to best-energy Rosetta structure
        - af3_{mode}_ligand_RMSD_bt: Binary-to-ternary ligand RMSD
    """
    af3_dir = pair_cache / f'af3_{mode}'
    summary_json = af3_dir / 'summary.json'

    prefix = f'af3_{mode}_'

    nan_result = {
        f'{prefix}ipTM': np.nan,
        f'{prefix}pLDDT_protein': np.nan,
        f'{prefix}pLDDT_ligand': np.nan,
        f'{prefix}interface_PAE': np.nan,
        f'{prefix}ligand_RMSD_min': np.nan,
        f'{prefix}ligand_RMSD_bestdG': np.nan,
        f'{prefix}ligand_RMSD_bt': np.nan,
    }

    if not summary_json.exists():
        return {f'{prefix}status': 'MISSING', **nan_result}

    try:
        with open(summary_json, 'r') as f:
            data = json.load(f)

        # Support both old single-key and new dual-key formats
        rmsd_min = data.get('ligand_RMSD_to_template_min',
                            data.get('ligand_RMSD_to_template', np.nan))
        rmsd_bestdG = data.get('ligand_RMSD_to_template_bestdG', np.nan)

        return {
            f'{prefix}status': 'COMPLETE',
            f'{prefix}ipTM': data.get('ipTM', np.nan),
            f'{prefix}pLDDT_protein': data.get('mean_pLDDT_protein', np.nan),
            f'{prefix}pLDDT_ligand': data.get('mean_pLDDT_ligand', np.nan),
            f'{prefix}interface_PAE': data.get('mean_interface_PAE', np.nan),
            f'{prefix}ligand_RMSD_min': rmsd_min,
            f'{prefix}ligand_RMSD_bestdG': rmsd_bestdG,
            f'{prefix}ligand_RMSD_bt': data.get('ligand_RMSD_binary_vs_ternary', np.nan),
        }

    except Exception as e:
        logger.warning(f"Failed to parse AF3 {mode} summary: {e}")
        return {f'{prefix}status': 'ERROR', **nan_result}


def extract_all_features(pair_cache: Path, pair_metadata: Dict) -> Dict:
    """
    Extract all features for a single (ligand, variant) pair.

    Args:
        pair_cache: Path to pair cache directory
        pair_metadata: Dict with pair_id, ligand_name, variant_name, etc.

    Returns:
        Dict with ~45 features
    """
    features = {}

    # Add input metadata
    features.update({
        'pair_id': pair_metadata.get('pair_id', ''),
        'ligand_name': pair_metadata.get('ligand_name', ''),
        'ligand_smiles': pair_metadata.get('ligand_smiles', ''),
        'variant_name': pair_metadata.get('variant_name', ''),
        'variant_signature': pair_metadata.get('variant_signature', ''),
        'label': pair_metadata.get('label', np.nan),
        'label_tier': pair_metadata.get('label_tier', ''),
        'label_source': pair_metadata.get('label_source', ''),
        'label_confidence': pair_metadata.get('label_confidence', ''),
        'affinity_EC50_uM': pair_metadata.get('affinity_EC50_uM', np.nan),
    })

    # Extract conformer features
    features.update(extract_conformer_features(pair_cache))

    # Extract docking cluster statistics (CRITICAL for ML!)
    features.update(extract_docking_features(pair_cache))

    # Extract Rosetta relax features
    features.update(extract_rosetta_relax_features(pair_cache))

    # Extract AF3 binary features
    features.update(extract_af3_features(pair_cache, mode='binary'))

    # Extract AF3 ternary features
    features.update(extract_af3_features(pair_cache, mode='ternary'))

    return features


# ═══════════════════════════════════════════════════════════════════
# MAIN AGGREGATION
# ═══════════════════════════════════════════════════════════════════

def aggregate_features(
    cache_dir: Path,
    pairs_csv: str,
    output_csv: str
) -> pd.DataFrame:
    """
    Aggregate features for all pairs in dataset.

    Args:
        cache_dir: Base cache directory containing pair subdirectories
        pairs_csv: CSV with pair metadata (pair_id, ligand_name, variant_name, label, etc.)
        output_csv: Output features CSV path

    Returns:
        DataFrame with features for all pairs
    """
    logger.info(f"Loading pairs metadata from {pairs_csv}...")
    pairs_df = pd.read_csv(pairs_csv)
    logger.info(f"Loaded {len(pairs_df)} pairs")

    all_features = []

    for idx, row in pairs_df.iterrows():
        pair_id = row['pair_id']
        pair_cache = cache_dir / pair_id

        logger.info(f"[{idx+1}/{len(pairs_df)}] Extracting features for {pair_id}...")

        # Extract features
        features = extract_all_features(pair_cache, row.to_dict())
        all_features.append(features)

    # Convert to DataFrame
    features_df = pd.DataFrame(all_features)

    # Save output
    features_df.to_csv(output_csv, index=False)
    logger.info(f"✓ Saved {len(features_df)} feature rows to {output_csv}")

    # Print summary statistics
    logger.info("\n" + "=" * 60)
    logger.info("FEATURE EXTRACTION SUMMARY")
    logger.info("=" * 60)

    for stage in ['conformer', 'docking', 'rosetta', 'af3_binary', 'af3_ternary']:
        status_col = f'{stage}_status'
        if status_col in features_df.columns:
            counts = features_df[status_col].value_counts()
            logger.info(f"{stage.upper()}: {counts.to_dict()}")

    # Completeness analysis
    logger.info("\nCOMPLETENESS:")
    for col in features_df.columns:
        if col.endswith('_status') or col in ['pair_id', 'ligand_name', 'variant_name']:
            continue
        missing_pct = (features_df[col].isna().sum() / len(features_df)) * 100
        if missing_pct > 0:
            logger.info(f"  {col}: {100-missing_pct:.1f}% complete")

    logger.info("=" * 60)

    return features_df


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate ML features from pipeline outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--cache-dir', required=True, help='Base cache directory')
    parser.add_argument('--pairs-csv', required=True, help='Pairs metadata CSV (with pair_id, label, etc.)')
    parser.add_argument('--output', required=True, help='Output features CSV')

    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)

    if not cache_dir.exists():
        logger.error(f"Cache directory not found: {cache_dir}")
        sys.exit(1)

    if not os.path.exists(args.pairs_csv):
        logger.error(f"Pairs CSV not found: {args.pairs_csv}")
        sys.exit(1)

    # Aggregate features
    features_df = aggregate_features(
        cache_dir=cache_dir,
        pairs_csv=args.pairs_csv,
        output_csv=args.output
    )

    logger.info(f"\n✓ Feature aggregation complete! Output: {args.output}")


if __name__ == '__main__':
    main()
