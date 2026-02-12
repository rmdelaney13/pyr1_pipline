#!/usr/bin/env python3

import os
import argparse
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

# Force matplotlib to not use Xserver backend (prevents display errors on clusters)
import matplotlib
matplotlib.use('Agg')

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y"
}


def extract_sequence_from_pdb(pdb_path, positions, chain_id=None):
    """
    Extract amino acids at specified residue numbers.
    Chain ID is optional (None = accept any chain).
    """
    seq = {}

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            # Chain handling (robust)
            chain = line[21].strip()
            if chain_id is not None and chain != chain_id:
                continue

            # Residue number handling (robust)
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue

            if resnum not in positions:
                continue

            resname = line[17:20].strip()
            aa = THREE_TO_ONE.get(resname)
            if aa:
                # Only record once per residue
                seq[resnum] = aa

    return seq


def build_frequency_table(pdb_dir, positions, chain_id=None):
    # Use a dictionary that preserves insertion order of 'positions'
    counts = {pos: Counter() for pos in positions}

    pdbs = sorted([f for f in os.listdir(pdb_dir) if f.endswith(".pdb")])
    print(f"[INFO] Found {len(pdbs)} PDBs in {pdb_dir}")

    total_hits = 0

    for pdb in pdbs:
        seq = extract_sequence_from_pdb(
            os.path.join(pdb_dir, pdb),
            positions,
            chain_id
        )
        for pos, aa in seq.items():
            if pos in counts:
                counts[pos][aa] += 1
                total_hits += 1

    print(f"[DEBUG] Total residue observations counted: {total_hits}")

    if total_hits == 0:
        raise RuntimeError(
            "No residues were extracted. "
            "Check residue numbering and chain handling."
        )

    df = pd.DataFrame.from_dict(counts, orient="index").fillna(0)
    df = df.reindex(columns=AA_ORDER, fill_value=0)
    
    # Ensure the dataframe follows the order of 'positions' list exactly
    df = df.reindex(positions)
    df.index.name = "Position"

    return df


def plot_logo(freq_df, positions, title, out_png):
    # Convert to probabilities
    prob_df = freq_df.div(freq_df.sum(axis=1), axis=0)

    # --- FIX START ---
    # Reset index so rows are 0, 1, 2... instead of 59, 81...
    # This aligns the data with the categorical x-axis 0..N
    plot_df = prob_df.reset_index(drop=True)
    # --- FIX END ---

    fig, ax = plt.subplots(figsize=(0.6 * len(positions), 4))

    logomaker.Logo(
        plot_df,
        ax=ax,
        color_scheme="chemistry",
        stack_order="small_on_top"
    )

    # --- categorical axis (no spacing gaps) ---
    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels([str(p) for p in positions], rotation=45)
    ax.set_xlim(-0.5, len(positions) - 0.5)
    ax.set_ylim(0, 1.0)

    ax.set_ylabel("AA Probability")
    ax.set_xlabel("Residue Position")
    ax.set_title(title)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"[RESULT] Logo written to {out_png}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate categorical sequence logo from a directory of PDBs"
    )
    parser.add_argument("pdb_dir", help="Directory containing filtered PDBs")
    parser.add_argument(
        "--positions",
        nargs="+",
        type=int,
        required=True,
        help="Residue positions to include (space-separated)"
    )
    parser.add_argument(
        "--chain",
        default=None,
        help="Protein chain ID (default: accept any chain)"
    )
    parser.add_argument("--title", default="Sequence Logo")
    parser.add_argument("--out", default="sequence_logo.png")

    args = parser.parse_args()

    positions = list(args.positions)

    freq_df = build_frequency_table(
        pdb_dir=args.pdb_dir,
        positions=positions,
        chain_id=args.chain
    )

    plot_logo(
        freq_df=freq_df,
        positions=positions,
        title=args.title,
        out_png=args.out
    )


if __name__ == "__main__":
    main()
