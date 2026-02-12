#!/usr/bin/env python3
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
        if s == "" or s.lower() in {"nan", "none", "null"}:
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
        os.link(src, dst)
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
    candidate_dir = binary_dir / design_id
    if candidate_dir.exists() and candidate_dir.is_dir():
        cif = find_cif_in_dir(candidate_dir)
        if cif:
            return cif
    patterns = [f"{design_id}*_model.cif", f"{design_id}*.cif"]
    for pat in patterns:
        hits = sorted(binary_dir.rglob(pat))
        if hits:
            hits.sort(key=lambda p: (len(str(p)), str(p)))
            return hits[0]
    return None

def aa3_to_aa1(resname: str) -> str:
    resname = (resname or "").strip().upper()
    mapping = {
        "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I",
        "LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S",
        "THR":"T","VAL":"V","TRP":"W","TYR":"Y",
        "MSE":"M","SEC":"U","PYL":"O",
    }
    return mapping.get(resname, "X")

def extract_masked_sequence_from_cif(cif_path: Path, chain_id: str, positions: List[int]) -> Optional[str]:
    """
    Robust masked extraction:
      - polymer residues only
      - map seqid.num -> AA1
      - masked sequence keeps order/length of positions; missing => 'X'
    """
    try:
        import gemmi
        st = gemmi.read_structure(str(cif_path))
        if len(st) == 0:
            return None
        model = st[0]

        chain = None
        for ch in model:
            if ch.name == chain_id:
                chain = ch
                break
        if chain is None:
            return None

        # polymer residues only (key fix)
        poly = chain.get_polymer()
        m: Dict[int, str] = {}
        for res in poly:
            m[res.seqid.num] = aa3_to_aa1(res.name)

        masked = "".join([m.get(p, "X") for p in positions])
        return masked if masked else None
    except Exception:
        return None

def seq_identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    if len(a) != len(b):
        n = min(len(a), len(b))
        matches = sum(1 for i in range(n) if a[i] == b[i])
        return matches / float(n) if n else 0.0
    matches = sum(1 for i in range(len(a)) if a[i] == b[i])
    return matches / float(len(a)) if a else 0.0

def greedy_cluster(items: List[Dict[str, Any]], id_threshold: float) -> List[List[Dict[str, Any]]]:
    clusters: List[List[Dict[str, Any]]] = []
    for it in items:
        placed = False
        for cl in clusters:
            rep = cl[0]
            if seq_identity(it["cluster_seq"], rep["cluster_seq"]) >= id_threshold:
                cl.append(it)
                placed = True
                break
        if not placed:
            clusters.append([it])
    return clusters

def main():
    ap = argparse.ArgumentParser(description="Pair binary+ternary metrics; diversify by clustering masked sequence.")
    ap.add_argument("--ternary_csv", required=True, type=Path)
    ap.add_argument("--binary_csv", required=True, type=Path)
    ap.add_argument("--ternary_dir", required=True, type=Path)
    ap.add_argument("--binary_dir", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)

    ap.add_argument("--out_csv", default="paired_ranked_metrics.csv")
    ap.add_argument("--diverse_csv", default="paired_ranked_diverse.csv")
    ap.add_argument("--diverse_dirname", default="top_cifs_diverse")

    ap.add_argument("--copy_mode", choices=["copy", "hardlink", "symlink"], default="hardlink")

    ap.add_argument("--w_dist", type=float, default=0.45)
    ap.add_argument("--w_plddt", type=float, default=0.35)
    ap.add_argument("--w_iptm", type=float, default=0.20)

    ap.add_argument("--preselect_n_plddt", type=int, default=100)
    ap.add_argument("--final_n_diverse", type=int, default=10)
    ap.add_argument("--cluster_identity", type=float, default=0.90)

    ap.add_argument("--seq_chain", default="A")
    ap.add_argument("--cluster_positions", default="59,81,83,92,94,108,110,117,120,122,141,159,160,163,164,167",
                    help="Comma-separated residue numbers used for clustering (masked).")
    ap.add_argument("--debug_seq", type=int, default=0,
                    help="Print first N extracted masked sequences for sanity check.")

    args = ap.parse_args()

    cluster_positions = [int(x.strip()) for x in args.cluster_positions.split(",") if x.strip()]

    ensure_dir(args.out_dir)
    diverse_dir = args.out_dir / args.diverse_dirname
    ensure_dir(diverse_dir)

    ternary_rows = read_csv(args.ternary_csv)
    binary_rows  = read_csv(args.binary_csv)

    bin_index: Dict[str, Dict[str, str]] = {}
    for r in binary_rows:
        tid = (r.get("target") or "").strip()
        if tid:
            bin_index[tid] = r

    paired: List[Dict[str, Any]] = []
    for tr in ternary_rows:
        t_target = (tr.get("target") or "").strip()
        if not t_target:
            continue

        design_id = normalize_ternary_target(t_target)

        t_plddt = safe_float(tr.get("ligand_plddt"))
        t_iptm  = safe_float(tr.get("ligand_iptm"))
        t_dist  = safe_float(tr.get("O3_O1_dist_aligned") or tr.get("o3_o1_dist_aligned"))

        br = bin_index.get(design_id)
        b_iptm = safe_float(br.get("chainB_iptm")) if br else None
        b_dist = safe_float(br.get("o3_o1_dist_aligned")) if br else None

        paired.append({
            "ternary_target": t_target,
            "design_id": design_id,
            "ternary_ligand_plddt": t_plddt,
            "ternary_ligand_iptm": t_iptm,
            "ternary_O3_O1_dist_aligned": t_dist,
            "binary_chainB_iptm": b_iptm,
            "binary_o3_o1_dist_aligned": b_dist,
            "has_binary_match": bool(br),
        })

    # composite score
    tdist_min, tdist_max = minmax([p["ternary_O3_O1_dist_aligned"] for p in paired])
    tpld_min, tpld_max   = minmax([p["ternary_ligand_plddt"] for p in paired])
    tipm_min, tipm_max   = minmax([p["ternary_ligand_iptm"] for p in paired])

    def composite(p: Dict[str, Any]) -> float:
        s = 0.0
        s += args.w_dist  * (norm01(p["ternary_O3_O1_dist_aligned"], tdist_min, tdist_max, invert=True) or 0.0)
        s += args.w_plddt * (norm01(p["ternary_ligand_plddt"], tpld_min, tpld_max) or 0.0)
        s += args.w_iptm  * (norm01(p["ternary_ligand_iptm"],  tipm_min, tipm_max) or 0.0)
        return s

    for p in paired:
        p["composite_score"] = composite(p)

    # write full ranked CSV (composite sort)
    paired_sorted = sorted(paired, key=lambda p: (
        -(p["composite_score"] if p["composite_score"] is not None else -1e9),
        (p["ternary_O3_O1_dist_aligned"] if p["ternary_O3_O1_dist_aligned"] is not None else 1e9),
        -(p["ternary_ligand_plddt"] if p["ternary_ligand_plddt"] is not None else -1e9),
        -(p["ternary_ligand_iptm"] if p["ternary_ligand_iptm"] is not None else -1e9),
        p["ternary_target"],
    ))

    out_csv = args.out_dir / args.out_csv
    base_fields = [
        "rank","composite_score","design_id","ternary_target",
        "ternary_ligand_plddt","ternary_ligand_iptm","ternary_O3_O1_dist_aligned",
        "binary_chainB_iptm","binary_o3_o1_dist_aligned","has_binary_match"
    ]
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=base_fields)
        w.writeheader()
        for rank, p in enumerate(paired_sorted, start=1):
            row = dict(p)
            row["rank"] = rank
            w.writerow({k: row.get(k, "") for k in base_fields})

    # preselect by pLDDT
    pre = sorted(paired, key=lambda p: (
        -(p["ternary_ligand_plddt"] if p["ternary_ligand_plddt"] is not None else -1e9),
        -(p["composite_score"] if p["composite_score"] is not None else -1e9),
        p["ternary_target"],
    ))[:args.preselect_n_plddt]

    enriched: List[Dict[str, Any]] = []
    missing_t = missing_seq = missing_b = 0

    for p in pre:
        t_dir = args.ternary_dir / p["ternary_target"]
        t_cif = find_cif_in_dir(t_dir) if (t_dir.exists() and t_dir.is_dir()) else None
        b_cif = find_binary_cif(args.binary_dir, p["design_id"])

        if not t_cif:
            missing_t += 1
            continue

        cluster_seq = extract_masked_sequence_from_cif(t_cif, args.seq_chain, cluster_positions)
        if not cluster_seq:
            missing_seq += 1
            continue

        q = dict(p)
        q["ternary_cif"] = str(t_cif)
        q["binary_cif"] = str(b_cif) if b_cif else ""
        q["cluster_seq"] = cluster_seq
        enriched.append(q)

        if not b_cif:
            missing_b += 1

    if args.debug_seq > 0:
        print("---- DEBUG masked sequences (cluster_seq) ----")
        for i, q in enumerate(enriched[:args.debug_seq], start=1):
            print(i, q["design_id"], q["cluster_seq"])
        print("--------------------------------------------")

    clusters = greedy_cluster(enriched, args.cluster_identity)
    reps = [cl[0] for cl in clusters]
    reps_sorted = sorted(reps, key=lambda p: (
        -(p["ternary_ligand_plddt"] if p["ternary_ligand_plddt"] is not None else -1e9),
        -(p["composite_score"] if p["composite_score"] is not None else -1e9),
        p["ternary_target"],
    ))
    chosen = reps_sorted[:args.final_n_diverse]

    # cluster metadata
    rep_meta: Dict[str, Tuple[int, int]] = {}
    for i, cl in enumerate(clusters, start=1):
        rep_meta[cl[0]["design_id"]] = (i, len(cl))

    diverse_fields = base_fields + ["cluster_id","cluster_size","ternary_cif","binary_cif","cluster_seq"]
    diverse_csv = args.out_dir / args.diverse_csv
    with diverse_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=diverse_fields)
        w.writeheader()

        for idx, p in enumerate(chosen, start=1):
            cid, csize = rep_meta.get(p["design_id"], (-1, 0))
            tag = f"{idx:04d}_cl{cid:03d}"

            t_cif = Path(p["ternary_cif"])
            b_cif = Path(p["binary_cif"]) if p.get("binary_cif") else None

            stage_file(t_cif, diverse_dir / f"{tag}_ternary_{p['design_id']}.cif", args.copy_mode)
            if b_cif and b_cif.exists():
                stage_file(b_cif, diverse_dir / f"{tag}_binary_{p['design_id']}.cif", args.copy_mode)

            row = dict(p)
            row["rank"] = ""
            row["cluster_id"] = cid
            row["cluster_size"] = csize
            w.writerow({k: row.get(k, "") for k in diverse_fields})

    print(f"[OK] Ranked CSV (all): {out_csv}")
    print(f"[OK] Diverse CSV: {diverse_csv}")
    print(f"[OK] Staged diverse CIFs into: {diverse_dir} (mode={args.copy_mode})")
    print(f"[INFO] Preselected by pLDDT: {len(pre)}")
    print(f"[INFO] Enriched (has ternary cif + seq): {len(enriched)}")
    print(f"[INFO] Missing ternary cif during enrich: {missing_t}")
    print(f"[INFO] Missing chain sequence during enrich: {missing_seq}")
    print(f"[INFO] Missing binary cif (still kept): {missing_b}")
    print(f"[INFO] Clusters formed: {len(clusters)} at identity >= {args.cluster_identity}")
    print(f"[INFO] Final diverse selected: {len(chosen)}")

if __name__ == "__main__":
    main()

