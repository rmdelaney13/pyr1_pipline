#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import shutil
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))


def safe_float(x) -> Optional[float]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if s == "" or s.lower() in {"nan", "none", "null", "na"}:
            return None
        return float(s)
    except Exception:
        return None


def normalize_ternary_target(t: str) -> str:
    t = t.strip()
    return t[:-len("_model")] if t.endswith("_model") else t


def minmax(vals: List[Optional[float]]) -> Tuple[Optional[float], Optional[float]]:
    v = [x for x in vals if x is not None]
    if not v:
        return None, None
    return min(v), max(v)


def norm01(x: Optional[float], vmin: Optional[float], vmax: Optional[float], invert: bool = False) -> Optional[float]:
    if x is None or vmin is None or vmax is None:
        return None
    if vmax == vmin:
        n = 0.5
    else:
        n = (x - vmin) / (vmax - vmin)
        n = max(0.0, min(1.0, n))
    return 1.0 - n if invert else n


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def stage_file(src: Path, dst: Path, mode: str) -> None:
    if dst.exists():
        return
    if mode == "symlink":
        dst.symlink_to(src)
    elif mode == "hardlink":
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)
    else:
        shutil.copy2(src, dst)


def find_cif_in_dir(model_dir: Path) -> Optional[Path]:
    cifs = sorted(model_dir.glob("*_model.cif"))
    if cifs:
        return cifs[0]
    cifs = sorted(model_dir.glob("*.cif"))
    if cifs:
        return cifs[0]
    return None


def find_binary_cif(binary_dir: Path, design_id: str) -> Optional[Path]:
    """
    Supports:
      1) directory-per-target: binary_dir/<design_id>/*.cif
      2) flat directory: binary_dir/<design_id>_model.cif
      3) recursive search
    """
    flat_cif = binary_dir / f"{design_id}_model.cif"
    if flat_cif.exists():
        return flat_cif

    candidate_dir = binary_dir / design_id
    if candidate_dir.exists() and candidate_dir.is_dir():
        cif = find_cif_in_dir(candidate_dir)
        if cif:
            return cif

    patterns = [
        f"{design_id}*_model.cif",
        f"{design_id}*.cif",
    ]
    for pat in patterns:
        hits = sorted(binary_dir.rglob(pat))
        if hits:
            hits.sort(key=lambda p: (len(str(p)), str(p)))
            return hits[0]

    return None


def get_dist(row: Dict[str, str]) -> Optional[float]:
    keys = ["min_dist_to_ligand_O_aligned", "O3_O1_dist_aligned", "o3_o1_dist_aligned", "Dist"]
    for k in keys:
        if k in row and row[k]:
            val = safe_float(row[k])
            if val is not None:
                return val
    return None


def get_rmsd(row: Dict[str, str], rmsd_key: Optional[str] = None) -> Optional[float]:
    """
    If --rmsd_key is provided, use it (and only it).
    Otherwise try common column names.
    """
    if rmsd_key:
        return safe_float(row.get(rmsd_key))

    keys = [
        "ligand_rmsd",
        "ligand_RMSD",
        "rmsd",
        "RMSD",
        "rmsd_aligned",
        "RMSD_aligned",
        "aligned_rmsd",
        "Aligned_RMSD",
        "rmsd_to_ref",
        "RMSD_to_ref",
    ]
    for k in keys:
        if k in row and row[k]:
            val = safe_float(row[k])
            if val is not None:
                return val
    return None


def main():
    ap = argparse.ArgumentParser(
        description="Pair binary+ternary metrics, rank (RMSD-heavy), and stage both CIFs into one directory."
    )
    ap.add_argument("--ternary_csv", required=True, type=Path)
    ap.add_argument("--binary_csv", required=True, type=Path)

    ap.add_argument("--ternary_dir", required=True, type=Path)
    ap.add_argument("--binary_dir", required=True, type=Path)

    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--out_csv", default="paired_ranked_metrics.csv")
    ap.add_argument("--stage_dirname", default="top_cifs")
    ap.add_argument("--top_n", type=int, default=10)
    ap.add_argument("--copy_mode", choices=["copy", "hardlink", "symlink"], default="hardlink")

    # RMSD handling
    ap.add_argument(
        "--rmsd_key",
        default=None,
        help="Exact RMSD column name to use (recommended once you know it). If omitted, tries common names."
    )
    ap.add_argument(
        "--require_rmsd",
        action="store_true",
        help="If set, drop rows that lack ternary RMSD (instead of giving them score 0)."
    )

    # weights (default: mostly RMSD)
    ap.add_argument("--w_rmsd", type=float, default=0.80)
    ap.add_argument("--w_dist", type=float, default=0.10)
    ap.add_argument("--w_plddt", type=float, default=0.07)
    ap.add_argument("--w_iptm", type=float, default=0.03)

    args = ap.parse_args()

    ensure_dir(args.out_dir)
    stage_dir = args.out_dir / args.stage_dirname
    ensure_dir(stage_dir)

    ternary_rows = read_csv(args.ternary_csv)
    binary_rows = read_csv(args.binary_csv)

    # Index binary metrics by target
    bin_index: Dict[str, Dict[str, str]] = {}
    for r in binary_rows:
        tid = (r.get("target") or "").strip()
        if tid:
            bin_index[tid] = r

    paired: List[Dict[str, Any]] = []
    dropped_no_rmsd = 0

    for tr in ternary_rows:
        t_target = (tr.get("target") or "").strip()
        if not t_target:
            continue

        design_id = normalize_ternary_target(t_target)

        # Ternary metrics
        t_plddt = safe_float(tr.get("ligand_plddt"))
        t_iptm = safe_float(tr.get("ligand_iptm"))
        t_dist = get_dist(tr)
        t_rmsd = get_rmsd(tr, rmsd_key=args.rmsd_key)

        if args.require_rmsd and t_rmsd is None:
            dropped_no_rmsd += 1
            continue

        # Binary metrics (carried along; not required)
        br = bin_index.get(design_id)
        if br:
            b_plddt = safe_float(br.get("ligand_plddt"))
            b_iptm = safe_float(br.get("ligand_iptm") or br.get("chainB_iptm"))
            b_dist = get_dist(br)
            b_rmsd = get_rmsd(br, rmsd_key=args.rmsd_key)
        else:
            b_plddt, b_iptm, b_dist, b_rmsd = None, None, None, None

        paired.append({
            "ternary_target": t_target,
            "design_id": design_id,
            "ternary_ligand_plddt": t_plddt,
            "ternary_ligand_iptm": t_iptm,
            "ternary_dist": t_dist,
            "ternary_rmsd": t_rmsd,
            "binary_ligand_plddt": b_plddt,
            "binary_ligand_iptm": b_iptm,
            "binary_dist": b_dist,
            "binary_rmsd": b_rmsd,
            "has_binary_match": bool(br),
        })

    # Normalization ranges (ternary-based)
    trmsd_min, trmsd_max = minmax([p["ternary_rmsd"] for p in paired])
    tdist_min, tdist_max = minmax([p["ternary_dist"] for p in paired])
    tpld_min, tpld_max = minmax([p["ternary_ligand_plddt"] for p in paired])
    tipm_min, tipm_max = minmax([p["ternary_ligand_iptm"] for p in paired])

    def composite(p: Dict[str, Any]) -> float:
        # RMSD is primary. If missing and not required, just score 0.
        if p["ternary_rmsd"] is None:
            return 0.0

        s = 0.0
        s += args.w_rmsd * (norm01(p["ternary_rmsd"], trmsd_min, trmsd_max, invert=True) or 0.0)

        # Secondary terms (optional if present)
        if p["ternary_dist"] is not None:
            s += args.w_dist * (norm01(p["ternary_dist"], tdist_min, tdist_max, invert=True) or 0.0)
        if p["ternary_ligand_plddt"] is not None:
            s += args.w_plddt * (norm01(p["ternary_ligand_plddt"], tpld_min, tpld_max) or 0.0)
        if p["ternary_ligand_iptm"] is not None:
            s += args.w_iptm * (norm01(p["ternary_ligand_iptm"], tipm_min, tipm_max) or 0.0)

        return s

    for p in paired:
        p["composite_score"] = composite(p)

    # Sort: composite desc, then RMSD asc, then dist asc, then pLDDT desc
    paired.sort(key=lambda p: (
        -(p["composite_score"] if p["composite_score"] is not None else -1e9),
        (p["ternary_rmsd"] if p["ternary_rmsd"] is not None else 1e9),
        (p["ternary_dist"] if p["ternary_dist"] is not None else 1e9),
        -(p["ternary_ligand_plddt"] if p["ternary_ligand_plddt"] is not None else -1e9),
    ))

    # Write aggregated CSV
    out_csv = args.out_dir / args.out_csv
    fieldnames = [
        "rank",
        "composite_score",
        "design_id",
        "ternary_target",
        "ternary_rmsd",
        "ternary_ligand_plddt",
        "ternary_ligand_iptm",
        "ternary_dist",
        "binary_rmsd",
        "binary_ligand_plddt",
        "binary_ligand_iptm",
        "binary_dist",
        "has_binary_match",
        "ternary_cif",
        "binary_cif",
    ]

    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()

        staged_t = 0
        staged_b = 0

        for rank, p in enumerate(paired, start=1):
            if rank > args.top_n:
                ternary_cif_path = ""
                binary_cif_path = ""
            else:
                # Find ternary CIF
                t_sub = args.ternary_dir / p["ternary_target"]
                if t_sub.exists() and t_sub.is_dir():
                    t_cif = find_cif_in_dir(t_sub)
                else:
                    t_cif = args.ternary_dir / f"{p['ternary_target']}_model.cif"
                    if not t_cif.exists():
                        t_cif = None

                # Find binary CIF
                b_cif = find_binary_cif(args.binary_dir, p["design_id"])

                ternary_cif_path = str(t_cif) if t_cif else ""
                binary_cif_path = str(b_cif) if b_cif else ""

                # Stage files
                rank_tag = f"{rank:04d}"

                if t_cif and t_cif.exists():
                    dst = stage_dir / f"{rank_tag}_ternary_{p['design_id']}.cif"
                    stage_file(t_cif, dst, args.copy_mode)
                    staged_t += 1

                if b_cif and b_cif.exists():
                    dst = stage_dir / f"{rank_tag}_binary_{p['design_id']}.cif"
                    stage_file(b_cif, dst, args.copy_mode)
                    staged_b += 1

            row_out = {
                "rank": rank,
                "composite_score": f"{p['composite_score']:.4f}",
                "design_id": p["design_id"],
                "ternary_target": p["ternary_target"],
                "ternary_rmsd": p["ternary_rmsd"],
                "ternary_ligand_plddt": p["ternary_ligand_plddt"],
                "ternary_ligand_iptm": p["ternary_ligand_iptm"],
                "ternary_dist": p["ternary_dist"],
                "binary_rmsd": p["binary_rmsd"],
                "binary_ligand_plddt": p["binary_ligand_plddt"],
                "binary_ligand_iptm": p["binary_ligand_iptm"],
                "binary_dist": p["binary_dist"],
                "has_binary_match": p["has_binary_match"],
                "ternary_cif": ternary_cif_path,
                "binary_cif": binary_cif_path,
            }
            w.writerow(row_out)

    print(f"[OK] Ranked CSV: {out_csv}")
    print(f"[OK] Staged directory: {stage_dir} (mode={args.copy_mode})")
    print(f"[INFO] Top-N: {args.top_n}")
    print(f"[INFO] Staged ternary CIFs: {staged_t}")
    print(f"[INFO] Staged binary  CIFs: {staged_b}")
    if args.require_rmsd:
        print(f"[INFO] Dropped (no ternary RMSD): {dropped_no_rmsd}")
    else:
        print("[INFO] Designs with missing ternary RMSD score 0. Use --require_rmsd to drop them.")


if __name__ == "__main__":
    main()

