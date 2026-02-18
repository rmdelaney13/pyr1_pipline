#!/usr/bin/env python3
"""
Prepare, submit, and analyze AlphaFold3 predictions for the ML dataset pipeline.

This script handles all AF3 logic for the ML pipeline:
  - Generating AF3 input JSONs from variant signatures + SMILES
  - Extracting metrics from AF3 output (ipTM, pLDDT, interface PAE)
  - Computing ligand RMSD to Rosetta-relaxed reference structures

Subcommands:
    prepare   Generate AF3 input JSONs for a single pair
    analyze   Extract metrics from AF3 output for a single pair

Author: Claude Code (Whitehead Lab PYR1 Pipeline)
Date: 2026-02-17
"""

import argparse
import copy
import json
import logging
import re
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# WT PYR1 sequence (181 residues) — 3QN1 stabilized construct.
# Includes 6 backbone-stabilizing mutations (29S, 43E, 80D, 90N, 118T, 134R)
# and the WT pocket residues (59K, 159F, 160A) as resolved in 3QN1.
WT_PYR1_SEQUENCE = (
    "MASELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKHFIKSCSV"
    "EQNFEMRVGCTRDVIVISGLPANTSTERLDILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENRI"
    "WTVVLESYVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMARN"
)
assert len(WT_PYR1_SEQUENCE) == 181, f"WT sequence length {len(WT_PYR1_SEQUENCE)} != 181"


# ═══════════════════════════════════════════════════════════════════
# VARIANT THREADING (sequence-level, no PyRosetta needed)
# ═══════════════════════════════════════════════════════════════════

def parse_variant_signature(signature: str) -> Dict[int, str]:
    """
    Parse variant signature into {position: target_amino_acid} dict.

    Supports formats:
        "59K;120A;160G"     → {59: 'K', 120: 'A', 160: 'G'}
        "K59Q;Y120A;A160G"  → {59: 'Q', 120: 'A', 160: 'G'}
        "K59Q_Y120A_A160G"  → {59: 'Q', 120: 'A', 160: 'G'}
        ""                  → {}  (wildtype)
    """
    if not signature or str(signature).strip() in ('', 'nan', 'None'):
        return {}

    mutations = {}
    normalized = str(signature).replace('_', ';').replace(' ', ';')
    normalized = re.sub(r';+', ';', normalized)

    for mut in normalized.split(';'):
        mut = mut.strip()
        if not mut:
            continue
        match = re.match(r'^([A-Z])?(\d+)([A-Z])$', mut)
        if match:
            _, pos, target_aa = match.groups()
            mutations[int(pos)] = target_aa
        else:
            logger.warning(f"Could not parse mutation: '{mut}' in '{signature}'")

    return mutations


def thread_mutations_to_sequence(
    wt_sequence: str,
    variant_signature: str
) -> str:
    """
    Thread mutations onto WT sequence by simple character substitution.

    Args:
        wt_sequence: WT PYR1 amino acid sequence (181 chars)
        variant_signature: Mutation signature (e.g., "59K;120A;160G") or empty

    Returns:
        Mutant sequence string
    """
    mutations = parse_variant_signature(variant_signature)
    if not mutations:
        return wt_sequence

    seq_list = list(wt_sequence)
    for pos, target_aa in mutations.items():
        idx = pos - 1  # Convert 1-indexed position to 0-indexed
        if 0 <= idx < len(seq_list):
            seq_list[idx] = target_aa
        else:
            logger.warning(f"Position {pos} out of range for sequence length {len(seq_list)}")

    return ''.join(seq_list)


# ═══════════════════════════════════════════════════════════════════
# AF3 JSON GENERATION
# ═══════════════════════════════════════════════════════════════════

def create_af3_json(
    template_path: str,
    mutant_sequence: str,
    ligand_smiles: str,
    pair_id: str,
    output_path: str,
    mode: str = 'binary'
) -> bool:
    """
    Create an AF3 input JSON from a template, replacing protein A sequence
    and ligand B SMILES.

    Args:
        template_path: Path to AF3 JSON template (binary or ternary)
        mutant_sequence: Mutant protein A sequence
        ligand_smiles: SMILES string for ligand B
        pair_id: Unique pair identifier (used in JSON "name" field)
        output_path: Path to write the output JSON
        mode: 'binary' or 'ternary' (for name suffix)

    Returns:
        True if successful, False otherwise
    """
    try:
        with open(template_path, 'r') as f:
            template = json.load(f)
    except Exception as e:
        logger.error(f"Failed to load template {template_path}: {e}")
        return False

    new_json = copy.deepcopy(template)

    # Find and update protein A sequence
    protein_a_found = False
    for entry in new_json.get('sequences', []):
        if 'protein' in entry and entry['protein'].get('id') == 'A':
            entry['protein']['sequence'] = mutant_sequence
            protein_a_found = True
            break

    if not protein_a_found:
        logger.error(f"No protein entry with id='A' found in template {template_path}")
        return False

    # Find and update ligand B SMILES
    ligand_b_found = False
    for entry in new_json.get('sequences', []):
        if 'ligand' in entry and entry['ligand'].get('id') == 'B':
            entry['ligand']['smiles'] = ligand_smiles
            ligand_b_found = True
            break

    if not ligand_b_found:
        logger.warning(f"No ligand entry with id='B' in template; SMILES not updated")

    # Set name
    new_json['name'] = f"{pair_id}_{mode}"

    # Write output
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(new_json, f, indent=2)

    return True


# ═══════════════════════════════════════════════════════════════════
# AF3 METRICS EXTRACTION
# ═══════════════════════════════════════════════════════════════════

def _get_chain_order(confidences: dict) -> List[str]:
    """Get ordered list of unique chain IDs from confidences data."""
    for key in ('token_chain_ids', 'atom_chain_ids'):
        if key in confidences and isinstance(confidences[key], list):
            seen = set()
            order = []
            for c in confidences[key]:
                c = str(c)
                if c not in seen:
                    order.append(c)
                    seen.add(c)
            return order
    return ['A', 'B']


def extract_af3_metrics(
    af3_output_dir: str,
    pair_id: str,
    mode: str,
    protein_chain: str = 'A',
    ligand_chain: str = 'B'
) -> Optional[Dict]:
    """
    Extract AF3 prediction metrics from output files.

    Looks for:
        {pair_id}_{mode}_summary_confidences.json
        {pair_id}_{mode}_confidences.json

    Returns dict with:
        ipTM, mean_pLDDT_protein, mean_pLDDT_ligand,
        mean_interface_PAE
    """
    output_dir = Path(af3_output_dir)
    name = f"{pair_id}_{mode}"

    summary_path = output_dir / f"{name}_summary_confidences.json"
    conf_path = output_dir / f"{name}_confidences.json"

    # AF3 may nest output in a subdirectory named after the prediction
    if not summary_path.exists():
        summary_path = output_dir / name / f"{name}_summary_confidences.json"
        conf_path = output_dir / name / f"{name}_confidences.json"

    if not summary_path.exists():
        logger.warning(f"Summary confidences not found for {name} in {af3_output_dir}")
        return None

    if not conf_path.exists():
        logger.warning(f"Full confidences not found for {name} in {af3_output_dir}")
        return None

    try:
        with open(summary_path, 'r') as f:
            summary = json.load(f)
        with open(conf_path, 'r') as f:
            confidences = json.load(f)
    except Exception as e:
        logger.error(f"Failed to parse AF3 output for {name}: {e}")
        return None

    metrics = {}

    # ── ipTM ──
    metrics['ipTM'] = summary.get('iptm', None)

    # ── Per-chain pLDDT from atom-level data ──
    atom_plddts = np.array(confidences.get('atom_plddts', []), dtype=float)
    atom_chains = np.array(confidences.get('atom_chain_ids', []), dtype=str)

    if len(atom_plddts) > 0 and len(atom_chains) == len(atom_plddts):
        protein_mask = atom_chains == protein_chain
        ligand_mask = atom_chains == ligand_chain

        metrics['mean_pLDDT_protein'] = (
            round(float(atom_plddts[protein_mask].mean()), 2)
            if protein_mask.any() else None
        )
        metrics['mean_pLDDT_ligand'] = (
            round(float(atom_plddts[ligand_mask].mean()), 2)
            if ligand_mask.any() else None
        )
    else:
        metrics['mean_pLDDT_protein'] = None
        metrics['mean_pLDDT_ligand'] = None

    # ── Interface PAE (cross-chain blocks from PAE matrix) ──
    pae_matrix = confidences.get('pae', None)
    token_chains = confidences.get('token_chain_ids', None)

    if pae_matrix is not None and token_chains is not None:
        pae = np.array(pae_matrix, dtype=float)
        tok_chains = np.array(token_chains, dtype=str)

        protein_tok_mask = tok_chains == protein_chain
        ligand_tok_mask = tok_chains == ligand_chain

        if protein_tok_mask.any() and ligand_tok_mask.any():
            # Cross-chain PAE: protein→ligand and ligand→protein blocks
            pae_prot_to_lig = pae[np.ix_(protein_tok_mask, ligand_tok_mask)]
            pae_lig_to_prot = pae[np.ix_(ligand_tok_mask, protein_tok_mask)]
            interface_pae = np.concatenate([
                pae_prot_to_lig.ravel(),
                pae_lig_to_prot.ravel()
            ])
            metrics['mean_interface_PAE'] = round(float(interface_pae.mean()), 2)
        else:
            metrics['mean_interface_PAE'] = None
    else:
        metrics['mean_interface_PAE'] = None

    return metrics


# ═══════════════════════════════════════════════════════════════════
# LIGAND RMSD CALCULATION
# ═══════════════════════════════════════════════════════════════════

def _resolve_element(atom) -> str:
    """Get element symbol from an atom, with fallback to atom name."""
    elem = atom.element.strip().upper() if atom.element else ''
    if not elem:
        # PyRosetta PDBs often lack the element column; infer from atom name
        name = atom.get_name().strip()
        if name:
            elem = name[0].upper()
    return elem


def _find_rosetta_ligand_chain(structure, protein_chain: str = 'A') -> Optional[Tuple[str, str]]:
    """
    Auto-detect ligand chain and residue name in a Rosetta PDB.

    PyRosetta may place the ligand on the same chain as the protein (e.g. chain A).
    This function searches ALL chains (including the protein chain) for non-amino-acid,
    non-water residues with the most heavy atoms.

    Returns:
        (chain_id, residue_name) tuple, or None if not found.
        The residue_name is used to filter atoms when ligand shares a chain with protein.
    """
    from Bio.PDB.Polypeptide import is_aa

    # Collect all non-protein, non-water residues across ALL chains
    candidates = []  # (chain_id, resname, heavy_count)
    for model in structure:
        for chain in model:
            for res in chain:
                resname = res.get_resname().strip()
                # Skip water
                if resname in ('HOH', 'WAT', 'TP3', 'TIP', 'TIP3'):
                    continue
                # Skip standard amino acids
                if is_aa(res, standard=True):
                    continue
                # Count heavy atoms
                heavy_count = sum(
                    1 for atom in res
                    if _resolve_element(atom) not in ('H', '')
                )
                if heavy_count > 0:
                    candidates.append((chain.id, resname, heavy_count))
        break  # Only first model

    if not candidates:
        # Log structure contents for debugging
        chain_info = []
        for model in structure:
            for chain in model:
                res_names = [r.get_resname().strip() for r in chain]
                chain_info.append(f"  Chain {chain.id}: {len(res_names)} residues: "
                                  f"{', '.join(res_names[:8])}{'...' if len(res_names) > 8 else ''}")
            break
        logger.warning(f"No ligand residues found. Structure contents:\n" + "\n".join(chain_info))
        return None

    # Pick the candidate with the most heavy atoms (the ligand, not ions/small molecules)
    best = max(candidates, key=lambda c: c[2])
    chain_id, resname, heavy_count = best
    logger.info(f"  Auto-detected Rosetta ligand: chain={chain_id}, resname={resname}, "
                f"heavy_atoms={heavy_count}")
    return (chain_id, resname)


def compute_ligand_rmsd_to_rosetta(
    af3_cif_path: str,
    rosetta_pdb_path: str,
    protein_chain: str = 'A',
    af3_ligand_chain: str = 'B',
    rosetta_ligand_chain: Optional[str] = None
) -> Optional[float]:
    """
    Compute ligand heavy-atom RMSD between AF3 prediction and Rosetta-relaxed
    structure, after CA-based superposition of the protein.

    Uses element-based matching: groups heavy atoms by element, then within
    each element group uses coordinate-proximity matching (Hungarian algorithm)
    to handle different atom orderings.

    Args:
        af3_cif_path: Path to AF3 output CIF file
        rosetta_pdb_path: Path to Rosetta-relaxed PDB file
        protein_chain: Protein chain ID (for CA alignment)
        af3_ligand_chain: Ligand chain ID in the AF3 CIF (default: 'B')
        rosetta_ligand_chain: Ligand chain ID in the Rosetta PDB. If None,
            auto-detects by finding the non-protein, non-water chain with
            the most heavy atoms.

    Returns:
        Ligand RMSD in Angstroms, or None on failure
    """
    try:
        from Bio.PDB import MMCIFParser, PDBParser, Superimposer
    except ImportError:
        logger.error("Biopython required for RMSD calculation (pip install biopython)")
        return None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Parse structures
        cif_parser = MMCIFParser(QUIET=True)
        pdb_parser = PDBParser(QUIET=True)

        try:
            af3_struct = cif_parser.get_structure("af3", af3_cif_path)
            ros_struct = pdb_parser.get_structure("ros", rosetta_pdb_path)
        except Exception as e:
            logger.error(f"Failed to parse structures: {e}")
            return None

    def _get_ca_atoms(structure, chain_id):
        """Extract CA atoms sorted by residue ID."""
        ca_atoms = []
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        for atom in res:
                            if atom.get_name() == 'CA':
                                ca_atoms.append(atom)
            break  # Only first model
        return ca_atoms

    def _get_ligand_heavy_atoms(structure, chain_id, ligand_resname=None):
        """Extract ligand heavy atoms (exclude H), robust to missing element fields.

        If ligand_resname is provided, only atoms from residues matching that
        name are returned (needed when ligand shares a chain with protein).
        """
        atoms = []
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        if ligand_resname and res.get_resname().strip() != ligand_resname:
                            continue
                        for atom in res:
                            elem = _resolve_element(atom)
                            if elem != 'H':
                                atoms.append(atom)
            break
        return atoms

    # Auto-detect Rosetta ligand chain if not specified
    rosetta_ligand_resname = None
    if rosetta_ligand_chain is None:
        detection = _find_rosetta_ligand_chain(ros_struct, protein_chain)
        if detection is None:
            logger.warning("Could not auto-detect ligand chain in Rosetta PDB")
            return None
        rosetta_ligand_chain, rosetta_ligand_resname = detection
        logger.info(f"  Rosetta ligand: chain={rosetta_ligand_chain}, resname={rosetta_ligand_resname}")

    # Get CA atoms for alignment
    af3_ca = _get_ca_atoms(af3_struct, protein_chain)
    ros_ca = _get_ca_atoms(ros_struct, protein_chain)

    if not af3_ca or not ros_ca:
        logger.warning(f"No CA atoms found for chain {protein_chain}")
        return None

    # Use min count for alignment (structures may differ slightly in length)
    n_ca = min(len(af3_ca), len(ros_ca))
    if n_ca < 10:
        logger.warning(f"Too few CA atoms for alignment: {n_ca}")
        return None

    af3_ca = af3_ca[:n_ca]
    ros_ca = ros_ca[:n_ca]

    # Superimpose AF3 onto Rosetta using CA atoms
    sup = Superimposer()
    try:
        sup.set_atoms(ros_ca, af3_ca)  # ros is fixed, af3 is mobile
        sup.apply(list(af3_struct.get_atoms()))  # Transform all AF3 atoms
    except Exception as e:
        logger.error(f"Superposition failed: {e}")
        return None

    # Get ligand heavy atoms (different chains for AF3 vs Rosetta)
    af3_lig = _get_ligand_heavy_atoms(af3_struct, af3_ligand_chain)
    ros_lig = _get_ligand_heavy_atoms(ros_struct, rosetta_ligand_chain, rosetta_ligand_resname)

    if not af3_lig or not ros_lig:
        logger.warning(f"No ligand heavy atoms: AF3 chain {af3_ligand_chain}={len(af3_lig) if af3_lig else 0}, "
                        f"Rosetta chain {rosetta_ligand_chain}={len(ros_lig) if ros_lig else 0}")
        return None

    # Check heavy atom count compatibility
    if abs(len(af3_lig) - len(ros_lig)) / max(len(af3_lig), len(ros_lig)) > 0.2:
        logger.warning(
            f"Ligand heavy atom count mismatch >20%: "
            f"AF3({af3_ligand_chain})={len(af3_lig)}, Rosetta({rosetta_ligand_chain})={len(ros_lig)}"
        )
        return None

    # Element-based matching with Hungarian algorithm
    # Group atoms by element
    af3_by_elem = {}
    for atom in af3_lig:
        elem = _resolve_element(atom)
        af3_by_elem.setdefault(elem, []).append(atom)

    ros_by_elem = {}
    for atom in ros_lig:
        elem = _resolve_element(atom)
        ros_by_elem.setdefault(elem, []).append(atom)

    # Match atoms within each element group
    all_sq_dists = []
    for elem in set(af3_by_elem) & set(ros_by_elem):
        af3_atoms = af3_by_elem[elem]
        ros_atoms = ros_by_elem[elem]

        if len(af3_atoms) == 1 and len(ros_atoms) == 1:
            # Trivial case
            d = af3_atoms[0].get_vector() - ros_atoms[0].get_vector()
            all_sq_dists.append(d.norm()**2)
        else:
            # Build pairwise distance matrix, use Hungarian algorithm
            try:
                from scipy.optimize import linear_sum_assignment
            except ImportError:
                # Fallback: sort by distance to centroid
                af3_coords = np.array([a.get_coord() for a in af3_atoms])
                ros_coords = np.array([a.get_coord() for a in ros_atoms])
                af3_centroid = af3_coords.mean(axis=0)
                ros_centroid = ros_coords.mean(axis=0)
                af3_dists = np.linalg.norm(af3_coords - af3_centroid, axis=1)
                ros_dists = np.linalg.norm(ros_coords - ros_centroid, axis=1)
                af3_order = np.argsort(af3_dists)
                ros_order = np.argsort(ros_dists)
                n = min(len(af3_order), len(ros_order))
                for i in range(n):
                    d = af3_atoms[af3_order[i]].get_coord() - ros_atoms[ros_order[i]].get_coord()
                    all_sq_dists.append(np.sum(d**2))
                continue

            af3_coords = np.array([a.get_coord() for a in af3_atoms])
            ros_coords = np.array([a.get_coord() for a in ros_atoms])

            # Cost matrix: pairwise squared distances
            n_af3, n_ros = len(af3_coords), len(ros_coords)
            cost = np.zeros((n_af3, n_ros))
            for i in range(n_af3):
                for j in range(n_ros):
                    diff = af3_coords[i] - ros_coords[j]
                    cost[i, j] = np.sum(diff**2)

            row_ind, col_ind = linear_sum_assignment(cost)
            for r, c in zip(row_ind, col_ind):
                all_sq_dists.append(cost[r, c])

    if not all_sq_dists:
        logger.warning("No matched atom pairs for RMSD calculation")
        return None

    rmsd = float(np.sqrt(np.mean(all_sq_dists)))
    return round(rmsd, 3)


# ═══════════════════════════════════════════════════════════════════
# BINARY vs TERNARY LIGAND RMSD
# ═══════════════════════════════════════════════════════════════════

def compute_binary_ternary_ligand_rmsd(
    binary_cif_path: str,
    ternary_cif_path: str,
    protein_chain: str = 'A',
    ligand_chain: str = 'B',
) -> Optional[float]:
    """
    Compute ligand RMSD between AF3 binary and ternary predictions.

    Aligns binary protein (chain A) onto ternary protein (chain A)
    using CA atoms, then measures heavy-atom RMSD of the ligand (chain B).
    This measures how consistently AF3 places the ligand with vs without
    the water molecule.

    Args:
        binary_cif_path: Path to AF3 binary output CIF file
        ternary_cif_path: Path to AF3 ternary output CIF file
        protein_chain: Protein chain ID (default: 'A')
        ligand_chain: Ligand chain ID (default: 'B')

    Returns:
        Ligand RMSD in Angstroms, or None on failure
    """
    try:
        from Bio.PDB import MMCIFParser, Superimposer
    except ImportError:
        logger.error("Biopython required for RMSD calculation (pip install biopython)")
        return None

    parser = MMCIFParser(QUIET=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            t_struct = parser.get_structure("ternary", ternary_cif_path)
            b_struct = parser.get_structure("binary", binary_cif_path)
        except Exception as e:
            logger.error(f"Failed to parse CIF structures: {e}")
            return None

    def _get_ca(struct, chain_id):
        ca = []
        for model in struct:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        for atom in res:
                            if atom.get_name() == 'CA':
                                ca.append(atom)
            break
        return ca

    def _get_ligand_heavy(struct, chain_id):
        atoms = []
        for model in struct:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        for atom in res:
                            elem = _resolve_element(atom)
                            if elem and elem != 'H':
                                atoms.append(atom)
            break
        atoms.sort(key=lambda a: a.get_name())
        return atoms

    t_ca = _get_ca(t_struct, protein_chain)
    b_ca = _get_ca(b_struct, protein_chain)

    if not t_ca or not b_ca:
        logger.warning("No CA atoms for binary-ternary protein alignment")
        return None

    if len(t_ca) != len(b_ca):
        logger.warning(f"CA count mismatch: binary={len(b_ca)}, ternary={len(t_ca)}")
        return None

    t_lig = _get_ligand_heavy(t_struct, ligand_chain)
    b_lig = _get_ligand_heavy(b_struct, ligand_chain)

    if not t_lig or not b_lig:
        logger.warning(f"No ligand atoms for binary-ternary RMSD: "
                       f"binary={len(b_lig) if b_lig else 0}, "
                       f"ternary={len(t_lig) if t_lig else 0}")
        return None

    if len(t_lig) != len(b_lig):
        logger.warning(f"Ligand atom count mismatch: binary={len(b_lig)}, ternary={len(t_lig)}")
        return None

    # Align binary onto ternary using protein CA
    sup = Superimposer()
    try:
        sup.set_atoms(t_ca, b_ca)  # ternary is fixed, binary is mobile
        sup.apply(list(b_struct.get_atoms()))
    except Exception as e:
        logger.error(f"Binary-ternary superposition failed: {e}")
        return None

    # Compute ligand RMSD
    diffs = []
    for a_t, a_b in zip(t_lig, b_lig):
        d = a_t.get_coord() - a_b.get_coord()
        diffs.append(np.sum(d * d))

    return round(float(np.sqrt(np.mean(diffs))), 3)


# ═══════════════════════════════════════════════════════════════════
# RELAXED PDB FINDERS
# ═══════════════════════════════════════════════════════════════════

def find_all_relaxed_pdbs(pair_cache: str) -> List[Path]:
    """
    Find all relaxed PDBs in the relax directory.

    Returns:
        List of Paths to relaxed PDBs (may be empty)
    """
    relax_dir = Path(pair_cache) / 'relax'
    if not relax_dir.exists():
        return []
    return sorted(p for p in relax_dir.glob('relaxed_*.pdb')
                  if '_score' not in p.stem)


def find_best_relaxed_pdb(pair_cache: str) -> Optional[Path]:
    """
    Find the best relaxed PDB (lowest dG_sep or total score) from score files.

    Args:
        pair_cache: Path to pair cache directory

    Returns:
        Path to best relaxed PDB, or None if not found
    """
    relax_dir = Path(pair_cache) / 'relax'
    if not relax_dir.exists():
        return None

    score_files = sorted(relax_dir.glob('relaxed_*_score.sc'))
    if not score_files:
        return None

    best_score = float('inf')
    best_pdb = None

    for sf in score_files:
        try:
            scores = {}
            with open(sf, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('SCORES:') or not line:
                        continue
                    if ':' in line:
                        key, _, value = line.partition(':')
                        try:
                            scores[key.strip()] = float(value.strip())
                        except ValueError:
                            pass

            # Prefer dG_sep, fallback to post_relax_total_score
            score = scores.get('dG_sep', scores.get('post_relax_total_score', float('inf')))
            if score < best_score:
                best_score = score
                # Derive PDB path from score file name
                pdb_name = sf.stem.replace('_score', '') + '.pdb'
                pdb_path = relax_dir / pdb_name
                if pdb_path.exists():
                    best_pdb = pdb_path

        except Exception:
            continue

    return best_pdb


def compute_min_ligand_rmsd_to_rosetta(
    af3_cif_path: str,
    relaxed_pdbs: List[Path],
    protein_chain: str = 'A',
    af3_ligand_chain: str = 'B',
) -> Optional[float]:
    """
    Compute ligand RMSD from AF3 prediction to each relaxed Rosetta PDB,
    return the minimum. This finds the Rosetta pose that best agrees with AF3.

    Returns:
        Minimum ligand RMSD in Angstroms, or None if all fail
    """
    rmsds = []
    for pdb_path in relaxed_pdbs:
        rmsd = compute_ligand_rmsd_to_rosetta(
            af3_cif_path=af3_cif_path,
            rosetta_pdb_path=str(pdb_path),
            protein_chain=protein_chain,
            af3_ligand_chain=af3_ligand_chain,
        )
        if rmsd is not None:
            rmsds.append(rmsd)

    if not rmsds:
        return None

    min_rmsd = min(rmsds)
    logger.info(f"  Ligand RMSD to Rosetta: min={min_rmsd:.3f} A "
                f"(across {len(rmsds)}/{len(relaxed_pdbs)} structures, "
                f"mean={sum(rmsds)/len(rmsds):.3f})")
    return min_rmsd


# ═══════════════════════════════════════════════════════════════════
# SUMMARY JSON WRITER
# ═══════════════════════════════════════════════════════════════════

def write_summary_json(
    output_dir: str,
    metrics: Dict,
    ligand_rmsd: Optional[float],
    binary_ternary_rmsd: Optional[float] = None,
) -> None:
    """
    Write summary.json with keys matching aggregate_ml_features.py expectations.

    Output format:
        {"ipTM": 0.85, "mean_pLDDT_protein": 78.5, "mean_pLDDT_ligand": 62.3,
         "mean_interface_PAE": 4.2, "ligand_RMSD_to_template": 1.8,
         "ligand_RMSD_binary_vs_ternary": 0.45}
    """
    summary = {
        'ipTM': metrics.get('ipTM'),
        'mean_pLDDT_protein': metrics.get('mean_pLDDT_protein'),
        'mean_pLDDT_ligand': metrics.get('mean_pLDDT_ligand'),
        'mean_interface_PAE': metrics.get('mean_interface_PAE'),
        'ligand_RMSD_to_template': ligand_rmsd,
        'ligand_RMSD_binary_vs_ternary': binary_ternary_rmsd,
    }

    output_path = Path(output_dir) / 'summary.json'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"  Wrote {output_path}")


# ═══════════════════════════════════════════════════════════════════
# CLI SUBCOMMANDS
# ═══════════════════════════════════════════════════════════════════

def cmd_prepare(args):
    """Generate AF3 input JSONs for a pair."""
    mutant_seq = thread_mutations_to_sequence(args.wt_sequence, args.variant_signature)

    logger.info(f"Variant: {args.variant_signature or 'WT'}")
    logger.info(f"SMILES: {args.ligand_smiles}")

    results = {}
    for mode in ('binary', 'ternary'):
        template = args.binary_template if mode == 'binary' else args.ternary_template
        output_path = Path(args.output_dir) / f'af3_{mode}' / f'{args.pair_id}_{mode}.json'

        success = create_af3_json(
            template_path=template,
            mutant_sequence=mutant_seq,
            ligand_smiles=args.ligand_smiles,
            pair_id=args.pair_id,
            output_path=str(output_path),
            mode=mode,
        )
        results[mode] = success
        if success:
            logger.info(f"  Created {mode} JSON: {output_path}")
        else:
            logger.error(f"  Failed to create {mode} JSON")

    return all(results.values())


def _find_af3_cif(af3_dir: str, pair_id: str, mode: str) -> Optional[Path]:
    """Locate the AF3 model CIF file for a given pair/mode."""
    name = f"{pair_id}_{mode}"
    cif_path = Path(af3_dir) / f"{name}_model.cif"
    if not cif_path.exists():
        cif_path = Path(af3_dir) / name / f"{name}_model.cif"
    return cif_path if cif_path.exists() else None


def cmd_analyze(args):
    """Extract AF3 metrics and compute ligand RMSD for a pair."""
    logger.info(f"Analyzing AF3 output for {args.pair_id}")

    # First pass: extract metrics and template RMSD for each mode
    mode_results = {}  # mode -> (metrics, ligand_rmsd, cif_path)

    for mode in ('binary', 'ternary'):
        af3_dir = args.af3_output_dir
        logger.info(f"\n  Processing {mode}...")

        metrics = extract_af3_metrics(
            af3_output_dir=af3_dir,
            pair_id=args.pair_id,
            mode=mode,
            protein_chain=args.protein_chain,
            ligand_chain=args.ligand_chain,
        )

        if metrics is None:
            logger.warning(f"  No AF3 metrics found for {mode}")
            continue

        logger.info(f"  ipTM={metrics.get('ipTM')}, "
                     f"pLDDT_prot={metrics.get('mean_pLDDT_protein')}, "
                     f"pLDDT_lig={metrics.get('mean_pLDDT_ligand')}, "
                     f"PAE={metrics.get('mean_interface_PAE')}")

        # Compute min ligand RMSD across all relaxed Rosetta structures
        ligand_rmsd = None
        cif_path = _find_af3_cif(af3_dir, args.pair_id, mode)

        if args.pair_cache:
            relaxed_pdbs = find_all_relaxed_pdbs(args.pair_cache)
            if relaxed_pdbs and cif_path:
                ligand_rmsd = compute_min_ligand_rmsd_to_rosetta(
                    af3_cif_path=str(cif_path),
                    relaxed_pdbs=relaxed_pdbs,
                    protein_chain=args.protein_chain,
                    af3_ligand_chain=args.ligand_chain,
                )
            elif not relaxed_pdbs:
                logger.warning(f"  No relaxed PDBs found for RMSD calculation")
            elif not cif_path:
                logger.warning(f"  AF3 model CIF not found for {mode}")

        mode_results[mode] = (metrics, ligand_rmsd, cif_path)

    # Compute binary-to-ternary ligand RMSD (consistency across water conditions)
    bt_rmsd = None
    if 'binary' in mode_results and 'ternary' in mode_results:
        b_cif = mode_results['binary'][2]
        t_cif = mode_results['ternary'][2]
        if b_cif and t_cif:
            bt_rmsd = compute_binary_ternary_ligand_rmsd(
                binary_cif_path=str(b_cif),
                ternary_cif_path=str(t_cif),
                protein_chain=args.protein_chain,
                ligand_chain=args.ligand_chain,
            )
            if bt_rmsd is not None:
                logger.info(f"\n  Binary-to-ternary ligand RMSD: {bt_rmsd:.3f} A")
            else:
                logger.warning("\n  Could not compute binary-to-ternary ligand RMSD")

    # Write summary.json for each mode
    for mode, (metrics, ligand_rmsd, _) in mode_results.items():
        output_dir = Path(args.pair_cache) / f'af3_{mode}' if args.pair_cache else Path(args.af3_output_dir)
        write_summary_json(str(output_dir), metrics, ligand_rmsd, bt_rmsd)


def main():
    parser = argparse.ArgumentParser(
        description='AF3 preparation and analysis for ML dataset pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='command', help='Subcommand')

    # ── prepare ──
    p_prep = subparsers.add_parser('prepare', help='Generate AF3 input JSONs')
    p_prep.add_argument('--pair-id', required=True, help='Unique pair ID')
    p_prep.add_argument('--variant-signature', default='', help='Mutation signature')
    p_prep.add_argument('--ligand-smiles', required=True, help='Ligand SMILES')
    p_prep.add_argument('--binary-template', required=True, help='Path to binary AF3 template JSON')
    p_prep.add_argument('--ternary-template', required=True, help='Path to ternary AF3 template JSON')
    p_prep.add_argument('--output-dir', required=True, help='Pair cache directory')
    p_prep.add_argument('--wt-sequence', default=WT_PYR1_SEQUENCE, help='WT protein sequence')

    # ── analyze ──
    p_ana = subparsers.add_parser('analyze', help='Extract AF3 metrics')
    p_ana.add_argument('--pair-id', required=True, help='Unique pair ID')
    p_ana.add_argument('--af3-output-dir', required=True, help='AF3 inference output directory')
    p_ana.add_argument('--pair-cache', default='', help='Pair cache directory (for RMSD + summary output)')
    p_ana.add_argument('--protein-chain', default='A', help='Protein chain ID')
    p_ana.add_argument('--ligand-chain', default='B', help='Ligand chain ID')

    args = parser.parse_args()

    if args.command == 'prepare':
        success = cmd_prepare(args)
        sys.exit(0 if success else 1)
    elif args.command == 'analyze':
        cmd_analyze(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
