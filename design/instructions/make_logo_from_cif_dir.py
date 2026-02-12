#!/usr/bin/env python3

import os
import argparse
import glob
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import gemmi

# Force matplotlib to not use Xserver backend
import matplotlib
matplotlib.use('Agg')

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

def extract_sequence_from_cif(cif_path, positions, target_chain=None):
    """
    Extract amino acids at specified residue numbers using Gemmi.
    """
    seq = {}
    
    try:
        # Parse CIF
        doc = gemmi.cif.read_file(cif_path)
        block = doc.sole_block()
        structure = gemmi.make_structure_from_block(block)
        
        # Use first model
        model = structure[0]
        
        for chain in model:
            # Check Chain ID if specified
            if target_chain and chain.name != target_chain:
                continue
                
            for res in chain:
                # Check residue number (seqid)
                res_num = res.seqid.num
                
                if res_num in positions:
                    # Get 1-letter code
                    # Gemmi handles non-standard residues gracefully, returning 'X' if unknown
                    aa = gemmi.find_tabulated_residue(res.name).one_letter_code
                    
                    if aa in AA_ORDER:
                        seq[res_num] = aa

    except Exception as e:
        print(f"[WARN] Failed to parse {cif_path}: {e}")

    return seq


def build_frequency_table(file_dir, positions, chain_id=None):
    # Use a dictionary that preserves insertion order of 'positions'
    counts = {pos: Counter() for pos in positions}

    # Find CIF files
    files = sorted(glob.glob(os.path.join(file_dir, "*.cif")))
    
    # Fallback to PDB if no CIFs found (just in case)
    if not files:
        files = sorted(glob.glob(os.path.join(file_dir, "*.pdb")))
        if files:
            print("[INFO] No CIF files found, falling back to PDBs...")
    
    print(f"[INFO] Found {len(files)} structure files in {file_dir}")

    total_hits = 0

    for fpath in files:
        # Helper to decide which parser to use based on extension
        if fpath.endswith(".cif"):
            seq = extract_sequence_from_cif(fpath, positions, chain_id)
        else:
            # Reuse CIF parser for PDB if gemmi can handle it (Gemmi reads PDBs too)
            # But specific PDB parser is safer if sticking to strict legacy formats
            # Here we assume CIF for AF3
            continue 

        for pos, aa in seq.items():
            if pos in counts:
                counts[pos][aa] += 1
                total_hits += 1

    print(f"[DEBUG] Total residue observations counted: {total_hits}")

    if total_hits == 0:
        raise RuntimeError(
            "No residues were extracted. "
            "Check residue numbering, chain ID, or file format."
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

    # Reset index so rows are 0, 1, 2... for plotting
    plot_df = prob_df.reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(0.6 * len(positions), 4))

    logomaker.Logo(
        plot_df,
        ax=ax,
        color_scheme="chemistry",
        stack_order="small_on_top"
    )

    # Categorical axis
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
        description="Generate categorical sequence logo from a directory of CIFs using Gemmi"
    )
    parser.add_argument("file_dir", help="Directory containing .cif files")
    parser.add_argument(
        "--positions",
        nargs="+",
        type=int,
        required=True,
        help="Residue positions to include (space-separated)"
    )
    parser.add_argument(
        "--chain",
        default="A",
        help="Protein chain ID (default: A)"
    )
    parser.add_argument("--title", default="Sequence Logo")
    parser.add_argument("--out", default="sequence_logo.png")

    args = parser.parse_args()

    positions = list(args.positions)

    freq_df = build_frequency_table(
        file_dir=args.file_dir,
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
