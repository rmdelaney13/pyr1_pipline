#!/usr/bin/env python3
"""
take_top_plddt_from_top_binary.py

Given:
  - a metrics CSV produced by sort_binary_models.py (already filtered by your criteria)
  - a directory containing the corresponding selected model files (top_binary)

This script selects the top N models by ligand pLDDT and copies/hardlinks/symlinks
the corresponding CIF files into a new output directory.

Expected filenames in --models_dir:
  <target><model_suffix>
Default: <target>_model.cif

Example:
  python take_top_plddt_from_top_binary.py \
    --csv        /scratch/alpine/ryde3462/DCA_design_flipped/af3/binary_outputs/top_binary_metrics.csv \
    --models_dir /scratch/alpine/ryde3462/DCA_design_flipped/af3/binary_outputs/top_binary \
    --out_dir    /scratch/alpine/ryde3462/DCA_design_flipped/af3/binary_outputs/top20_plddt \
    --out_csv    /scratch/alpine/ryde3462/DCA_design_flipped/af3/binary_outputs/top20_plddt/top20_plddt_metrics.csv \
    --n 20 \
    --id_col target \
    --plddt_col ligand_plddt \
    --copy_mode hardlink
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from pathlib import Path
from typing import Optional, List, Dict


def safe_float(v: Optional[str]) -> Optional[float]:
    if v is None:
        return None
    s = str(v).strip()
    if s == "" or s.lower() in {"nan", "none", "null", "na"}:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def link_or_copy(src: Path, dst: Path, mode: str) -> None:
    if mode == "copy":
        shutil.copy2(src, dst)
    elif mode == "hardlink":
        os.link(src, dst)
    elif mode == "symlink":
        dst.symlink_to(src)
    else:
        raise ValueError(f"Unknown copy mode: {mode}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Take top-N pLDDT models from a top_binary directory + metrics CSV.")
    ap.add_argument("--csv", required=True, help="Filtered metrics CSV (e.g., top_binary_metrics.csv).")
    ap.add_argument("--models_dir", required=True, help="Directory containing selected model files (e.g., top_binary).")
    ap.add_argument("--out_dir", required=True, help="Directory to write the top-N model files.")
    ap.add_argument("--out_csv", help="Optional output CSV for the top-N subset.")
    ap.add_argument("--n", type=int, default=20, help="Number of top models to keep (default 20).")

    ap.add_argument("--id_col", default="target", help="ID column name (default target).")
    ap.add_argument("--plddt_col", default="ligand_plddt", help="pLDDT column name (default ligand_plddt).")

    ap.add_argument("--model_suffix", default="_model.cif", help="Model filename suffix (default _model.cif).")
    ap.add_argument("--copy_mode", choices=["copy", "hardlink", "symlink"], default="copy", help="How to move files.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing files in out_dir.")
    ap.add_argument("--dedupe_by_id", action="store_true",
                    help="If set, keep only the best row per ID before taking top N (useful if duplicates exist).")
    ap.add_argument("--dry_run", action="store_true", help="Print actions without writing files.")

    args = ap.parse_args()

    csv_path = Path(args.csv).resolve()
    models_dir = Path(args.models_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_csv = Path(args.out_csv).resolve() if args.out_csv else None

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    if not models_dir.exists():
        raise FileNotFoundError(f"models_dir not found: {models_dir}")

    # Read rows
    rows: List[Dict[str, str]] = []
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV appears to have no header.")
        fieldnames = reader.fieldnames

        required = {args.id_col, args.plddt_col}
        missing = required - set(fieldnames)
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}")

        for r in reader:
            seq_id = (r.get(args.id_col) or "").strip()
            plddt = safe_float(r.get(args.plddt_col))
            if not seq_id or plddt is None:
                continue
            r["_plddt_float"] = str(plddt)  # stash for debugging if needed
            rows.append(r)

    if not rows:
        print("No usable rows found (missing IDs or pLDDT). Nothing to do.")
        return

    # Optionally dedupe by ID (keep best pLDDT row per ID)
    if args.dedupe_by_id:
        best: Dict[str, Dict[str, str]] = {}
        for r in rows:
            seq_id = r[args.id_col].strip()
            p = float(r["_plddt_float"])
            if (seq_id not in best) or (p > float(best[seq_id]["_plddt_float"])):
                best[seq_id] = r
        rows = list(best.values())

    # Sort by pLDDT desc
    rows.sort(key=lambda r: float(r["_plddt_float"]), reverse=True)

    # Take top N that actually have files
    selected: List[Dict[str, str]] = []
    missing_files = 0

    for r in rows:
        if len(selected) >= args.n:
            break
        seq_id = r[args.id_col].strip()
        src = models_dir / f"{seq_id}{args.model_suffix}"
        if not src.exists():
            missing_files += 1
            continue
        selected.append(r)

    print("=== TOP pLDDT SELECTION ===")
    print(f"Rows in CSV (usable):      {len(rows)}")
    print(f"Requested top N:           {args.n}")
    print(f"Selected with files found: {len(selected)}")
    print(f"Missing model files:       {missing_files}")

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    # Copy/link models
    for r in selected:
        seq_id = r[args.id_col].strip()
        src = models_dir / f"{seq_id}{args.model_suffix}"
        dst = out_dir / src.name

        if dst.exists():
            if args.overwrite:
                if not args.dry_run:
                    dst.unlink()
            else:
                continue

        if args.dry_run:
            print(f"[DRY] {args.copy_mode}: {src} -> {dst}")
            continue

        try:
            link_or_copy(src, dst, args.copy_mode)
        except FileExistsError:
            pass

    # Write out subset CSV
    if out_csv and selected:
        if not args.dry_run:
            out_csv.parent.mkdir(parents=True, exist_ok=True)
            with out_csv.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=[fn for fn in fieldnames if fn is not None])
                w.writeheader()
                for r in selected:
                    rr = {k: v for k, v in r.items() if k in fieldnames}
                    w.writerow(rr)
        print(f"Top-N CSV written to: {out_csv}")

    print(f"Top-N models written to: {out_dir}")


if __name__ == "__main__":
    main()

