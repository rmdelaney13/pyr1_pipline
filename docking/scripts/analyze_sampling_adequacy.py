#!/usr/bin/env python
"""
Analyze whether you have sufficient docking sampling.

This script helps you determine if you've run enough conformers by analyzing:
1. Cluster saturation (are new clusters still being discovered?)
2. Score distribution and convergence
3. Sampling efficiency metrics

Usage:
    python analyze_sampling_adequacy.py config.txt
    python analyze_sampling_adequacy.py config.txt --geometry-csvs "hbond_*.csv"
"""

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict
from configparser import ConfigParser

try:
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    print("Warning: matplotlib/numpy not available. Will skip plotting.")


def analyze_geometry_csvs(csv_pattern, output_dir):
    """Analyze H-bond geometry CSV files to assess sampling."""
    csv_files = sorted(glob.glob(csv_pattern))

    if not csv_files:
        print(f"No CSV files found matching: {csv_pattern}")
        return None

    print(f"\nFound {len(csv_files)} geometry CSV files")

    all_data = []
    total_conformers = 0
    passed_geometry = 0
    passed_score = 0
    saved_clusters = 0

    for csv_file in csv_files:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                all_data.append(row)
                total_conformers += 1
                if row.get('accepted') == 'True':
                    passed_geometry += 1
                if row.get('passed_score') == 'True':
                    passed_score += 1
                if row.get('saved_cluster') == 'True':
                    saved_clusters += 1

    print(f"\n{'='*60}")
    print("SAMPLING SUMMARY")
    print(f"{'='*60}")
    print(f"Total conformers processed:     {total_conformers:6d}")
    print(f"Passed geometry filter:         {passed_geometry:6d}  ({100*passed_geometry/max(1,total_conformers):.1f}%)")
    print(f"Passed score cutoff:            {passed_score:6d}  ({100*passed_score/max(1,total_conformers):.1f}%)")
    print(f"Unique clusters saved:          {saved_clusters:6d}  ({100*saved_clusters/max(1,total_conformers):.1f}%)")
    print(f"{'='*60}")

    # Assess sampling adequacy
    print(f"\n{'='*60}")
    print("SAMPLING ADEQUACY ASSESSMENT")
    print(f"{'='*60}")

    geometry_rate = passed_geometry / max(1, total_conformers)
    score_rate = passed_score / max(1, total_conformers)
    cluster_rate = saved_clusters / max(1, total_conformers)

    if geometry_rate < 0.01:
        print("⚠️  WARNING: Very low geometry pass rate (<1%)")
        print("   → Consider relaxing H-bond filter settings")
        print("   → Check MaxPerturbTries, HBondDistanceIdealBuffer, HBondDonorAngleMin")
    elif geometry_rate < 0.05:
        print("⚠️  Low geometry pass rate (<5%) - filters are quite strict")
    else:
        print(f"✓  Geometry pass rate is reasonable ({geometry_rate*100:.1f}%)")

    if score_rate < 0.001:
        print("⚠️  WARNING: Very low score pass rate (<0.1%)")
        print("   → Consider relaxing MaxScore cutoff")
        print("   → Check if binding site is too constrained")
    elif score_rate < 0.01:
        print("⚠️  Low score pass rate (<1%) - very selective scoring")
    else:
        print(f"✓  Score pass rate seems reasonable ({score_rate*100:.1f}%)")

    if saved_clusters < 10:
        print(f"⚠️  WARNING: Only {saved_clusters} unique clusters found")
        print("   → Insufficient sampling diversity")
        print("   → Consider running more conformers or relaxing filters")
    elif saved_clusters < 50:
        print(f"⚠️  Moderate cluster count ({saved_clusters}) - may want more sampling")
    else:
        print(f"✓  Good cluster diversity ({saved_clusters} unique clusters)")

    # Cluster saturation analysis
    if passed_score > 0:
        redundancy = passed_score / max(1, saved_clusters)
        print(f"\nCluster redundancy: {redundancy:.1f} poses per cluster")
        if redundancy < 1.5:
            print("   → Low redundancy: May need MORE conformers for saturation")
        elif redundancy < 5:
            print("   → Moderate redundancy: Sampling is reasonable")
        else:
            print("   → High redundancy: Likely approaching saturation")

    print(f"{'='*60}")

    # Score distribution analysis
    scores = []
    for row in all_data:
        if row.get('passed_score') == 'True' and row.get('score'):
            try:
                scores.append(float(row['score']))
            except ValueError:
                pass

    if scores and HAS_PLOTTING:
        scores = np.array(scores)
        print(f"\nSCORE STATISTICS (n={len(scores)})")
        print(f"  Best (min):     {np.min(scores):.1f}")
        print(f"  Median:         {np.median(scores):.1f}")
        print(f"  Mean:           {np.mean(scores):.1f}")
        print(f"  Worst (max):    {np.max(scores):.1f}")
        print(f"  Std Dev:        {np.std(scores):.1f}")

        # Plot score distribution
        plt.figure(figsize=(10, 6))
        plt.hist(scores, bins=50, edgecolor='black', alpha=0.7)
        plt.xlabel('Rosetta Score')
        plt.ylabel('Count')
        plt.title(f'Score Distribution (n={len(scores)} passing poses)')
        plt.axvline(np.median(scores), color='red', linestyle='--', label=f'Median: {np.median(scores):.1f}')
        plt.legend()
        plt.grid(True, alpha=0.3)

        plot_file = os.path.join(output_dir, 'score_distribution.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"\n✓ Score distribution plot saved: {plot_file}")
        plt.close()

    return {
        'total_conformers': total_conformers,
        'passed_geometry': passed_geometry,
        'passed_score': passed_score,
        'saved_clusters': saved_clusters,
        'scores': scores if scores else []
    }


def recommend_rmsd_cutoff(ligand_heavy_atoms=None):
    """Provide RMSD cutoff recommendations based on ligand size."""
    print(f"\n{'='*60}")
    print("RMSD CUTOFF RECOMMENDATIONS")
    print(f"{'='*60}")

    if ligand_heavy_atoms is not None:
        print(f"Your ligand: ~{ligand_heavy_atoms} heavy atoms")
        if ligand_heavy_atoms < 10:
            rec_rmsd = 0.5
            print(f"  Recommended RMSD cutoff: {rec_rmsd} Å")
            print(f"  Rationale: Small rigid ligands need tight clustering")
        elif ligand_heavy_atoms < 20:
            rec_rmsd = 0.75
            print(f"  Recommended RMSD cutoff: {rec_rmsd} Å")
            print(f"  Rationale: Medium-sized ligands (standard)")
        elif ligand_heavy_atoms < 30:
            rec_rmsd = 1.0
            print(f"  Recommended RMSD cutoff: {rec_rmsd} Å")
            print(f"  Rationale: Larger ligands have more conformational space")
        else:
            rec_rmsd = 1.5
            print(f"  Recommended RMSD cutoff: {rec_rmsd} Å")
            print(f"  Rationale: Very large/flexible ligands")
    else:
        print("General guidelines (by ligand heavy atom count):")
        print("  < 10 atoms:   0.5 Å  (small, rigid)")
        print("  10-20 atoms:  0.75 Å (standard, medium flexibility)")
        print("  20-30 atoms:  1.0 Å  (larger, more flexible)")
        print("  > 30 atoms:   1.5 Å  (very large/flexible)")

    print("\nAdjustment criteria:")
    print("  • If you get too many clusters (>200): INCREASE cutoff")
    print("  • If you get too few clusters (<10):   DECREASE cutoff")
    print("  • If poses in same cluster look different: DECREASE cutoff")
    print("  • If clearly different poses clustered together: INCREASE cutoff")
    print(f"{'='*60}")


def main(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument('config_file', help='Path to config file')
    parser.add_argument(
        '--geometry-csvs',
        default=None,
        help='Glob pattern for geometry CSV files (default: auto-detect from config OutputDir)'
    )
    parser.add_argument(
        '--ligand-heavy-atoms',
        type=int,
        default=None,
        help='Number of heavy atoms in your ligand (for RMSD cutoff recommendations)'
    )

    args = parser.parse_args(argv)

    # Load config
    config = ConfigParser()
    with open(args.config_file, 'r', encoding='utf-8-sig') as f:
        config.read_file(f)

    # Determine output directory
    try:
        import docking_pipeline_utils as dpu
        output_dir = dpu.cfg_get(config, "grade_conformers", "OutputDir", "output")
    except ImportError:
        spec = config.get("grade_conformers", {})
        output_dir = spec.get("OutputDir", config.get("DEFAULT", {}).get("SCRATCH_ROOT", "output"))
        if "docked" not in output_dir:
            output_dir = os.path.join(output_dir, "docked")

    output_dir = os.path.abspath(output_dir)

    # Determine geometry CSV pattern
    if args.geometry_csvs:
        csv_pattern = args.geometry_csvs
    else:
        csv_pattern = os.path.join(output_dir, "hbond_geometry_summary*.csv")

    print(f"Configuration file: {args.config_file}")
    print(f"Output directory:   {output_dir}")
    print(f"CSV pattern:        {csv_pattern}")

    # Analyze geometry CSVs
    results = analyze_geometry_csvs(csv_pattern, output_dir)

    # RMSD recommendations
    recommend_rmsd_cutoff(args.ligand_heavy_atoms)

    # Final recommendations
    if results:
        print(f"\n{'='*60}")
        print("RECOMMENDATIONS")
        print(f"{'='*60}")

        if results['saved_clusters'] < 20:
            print("\n🔴 INSUFFICIENT SAMPLING")
            print("   Recommendations:")
            print("   1. Increase ArrayTaskCount to run more conformers")
            print("   2. OR relax filtering parameters (MaxScore, H-bond geometry)")
            print("   3. Ensure your input SDF has diverse conformers")
            print("   4. Consider adjusting RMSD cutoff if poses are over-clustered")

        elif results['saved_clusters'] > 200:
            print("\n🟡 POSSIBLE OVER-CLUSTERING")
            print("   Recommendations:")
            print("   1. INCREASE ClusterRMSDCutoff (try 1.0 or 1.5 Å)")
            print("   2. Visually inspect if clusters contain similar poses")
            print("   3. You may have very good sampling diversity!")

        else:
            print("\n🟢 SAMPLING APPEARS ADEQUATE")
            print("   Recommendations:")
            print("   1. Visually inspect top clusters in PyMOL/Chimera")
            print("   2. Check if binding modes cover expected interactions")
            print("   3. If top clusters are very similar, may be converged")
            print("   4. If you need more diversity, run additional conformers")

        redundancy = results['passed_score'] / max(1, results['saved_clusters'])
        if redundancy > 10:
            print("\n   Note: High redundancy suggests you're approaching saturation")
            print("         (many conformers clustering to same binding modes)")

        print(f"{'='*60}")


if __name__ == '__main__':
    main(sys.argv[1:])
