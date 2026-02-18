#!/usr/bin/env python3
"""
Collect best Rosetta relaxed PDBs and AF3 binary CIFs into a single folder.

Usage:
    python collect_models.py --cache-dir /scratch/alpine/ryde3462/ml_dataset_test \
                             --output-dir ./collected_models
"""
from __future__ import annotations

import argparse
import shutil
import json
from pathlib import Path


def find_best_relaxed_pdb(pair_cache: Path) -> Path | None:
    relax_dir = pair_cache / 'relax'
    if not relax_dir.exists():
        return None
    best_score, best_pdb = float('inf'), None
    for sf in relax_dir.glob('relaxed_*_score.sc'):
        try:
            scores = {}
            for line in sf.read_text().splitlines():
                line = line.strip()
                if line.startswith('SCORES:') or not line or ':' not in line:
                    continue
                k, _, v = line.partition(':')
                try:
                    scores[k.strip()] = float(v.strip())
                except ValueError:
                    pass
            score = scores.get('dG_sep', scores.get('post_relax_total_score', float('inf')))
            if score < best_score:
                pdb = relax_dir / (sf.stem.replace('_score', '') + '.pdb')
                if pdb.exists():
                    best_score, best_pdb = score, pdb
        except Exception:
            continue
    return best_pdb


def find_af3_binary_cif(pair_cache: Path, pair_id: str) -> Path | None:
    name = f"{pair_id}_binary"
    # Check pair-local af3_binary dir
    for pattern in [
        pair_cache / 'af3_binary' / f'{name}_model.cif',
        pair_cache / 'af3_binary' / name / f'{name}_model.cif',
    ]:
        if pattern.exists():
            return pattern
    # Check staging dir (sibling of pair dirs)
    staging = pair_cache.parent / 'af3_staging' / 'binary_output'
    if staging.exists():
        for pattern in [
            staging / f'{name}_model.cif',
            staging / name / f'{name}_model.cif',
        ]:
            if pattern.exists():
                return pattern
    return None


def main():
    parser = argparse.ArgumentParser(description='Collect Rosetta + AF3 models')
    parser.add_argument('--cache-dir', required=True)
    parser.add_argument('--output-dir', default='./collected_models')
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for pair_dir in sorted(cache_dir.iterdir()):
        if not pair_dir.is_dir() or pair_dir.name.startswith('af3_staging'):
            continue
        pair_id = pair_dir.name

        pdb = find_best_relaxed_pdb(pair_dir)
        cif = find_af3_binary_cif(pair_dir, pair_id)

        if pdb:
            shutil.copy2(pdb, out / f'{pair_id}_rosetta.pdb')
            print(f"  {pair_id}: rosetta PDB copied (dG_sep best)")
        else:
            print(f"  {pair_id}: NO rosetta PDB found")

        if cif:
            shutil.copy2(cif, out / f'{pair_id}_af3_binary.cif')
            print(f"  {pair_id}: AF3 binary CIF copied")
        else:
            print(f"  {pair_id}: NO AF3 binary CIF found")

    print(f"\nCollected to: {out}")


if __name__ == '__main__':
    main()
