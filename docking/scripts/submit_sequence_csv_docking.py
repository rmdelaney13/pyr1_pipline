#!/usr/bin/env python
"""
Helper script to submit sequence-CSV docking jobs to SLURM.

Automatically determines the array size from the CSV and generates submission scripts.

Usage:
    # Interactive: prepare and show submission command
    python submit_sequence_csv_docking.py config.txt

    # Auto-submit to SLURM
    python submit_sequence_csv_docking.py config.txt --submit

    # Preview sequences without submitting
    python submit_sequence_csv_docking.py config.txt --preview

    # Run locally for testing (first 2 sequences)
    python submit_sequence_csv_docking.py config.txt --local --max-local 2
"""

import argparse
import os
import sys
import subprocess
from configparser import ConfigParser

import pandas as pd


def count_sequences_in_csv(csv_path):
    """Count number of sequences in CSV."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Sequence CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)

    required_cols = {"name", "signature"}
    if not required_cols.issubset(set(df.columns)):
        raise ValueError(
            f"CSV must contain columns: {required_cols}. Found: {df.columns.tolist()}"
        )

    print(f"\nFound {len(df)} sequences in {csv_path}:")
    for idx, row in df.iterrows():
        print(f"  [{idx}] {row['name']}: {row['signature']}")

    return len(df), df


def generate_slurm_script(
    config_path,
    num_sequences,
    job_name="seq_dock",
    time="24:00:00",
    partition="amilan",
    qos="normal",
    ntasks=1,
    mem="16G",
    conda_env="ligand_alignment",
    output_dir="slurm_logs",
):
    """Generate SLURM submission script."""

    script_dir = os.path.dirname(os.path.abspath(__file__))
    docking_script = os.path.join(script_dir, "grade_conformers_sequence_csv_docking_multiple_slurm.py")

    if not os.path.exists(docking_script):
        raise FileNotFoundError(f"Docking script not found: {docking_script}")

    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --array=0-{num_sequences - 1}
#SBATCH --time={time}
#SBATCH --partition={partition}
#SBATCH --qos={qos}
#SBATCH --ntasks={ntasks}
#SBATCH --mem={mem}
#SBATCH --output={output_dir}/{job_name}_%A_%a.out
#SBATCH --error={output_dir}/{job_name}_%A_%a.err

# Load conda environment
source ~/.bashrc
conda activate {conda_env}

# Print job info
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=========================================="

# Run docking for this sequence
python {docking_script} {config_path} $SLURM_ARRAY_TASK_ID

echo "=========================================="
echo "End time: $(date)"
echo "=========================================="
"""

    return slurm_script


def run_local_test(config_path, max_sequences=None):
    """Run docking locally for testing (non-SLURM)."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    docking_script = os.path.join(script_dir, "grade_conformers_sequence_csv_docking_multiple_slurm.py")

    config = ConfigParser()
    with open(config_path, "r", encoding="utf-8-sig") as f:
        config.read_file(f)

    seq_section = config.get("SEQUENCE_DOCKING", fallback={})
    csv_path = seq_section.get("SequenceCSV", None)

    if not csv_path:
        raise ValueError("SequenceCSV not specified in config [SEQUENCE_DOCKING] section")

    num_sequences, df = count_sequences_in_csv(csv_path)

    if max_sequences is not None:
        num_sequences = min(num_sequences, max_sequences)

    print(f"\n{'='*80}")
    print(f"Running LOCAL test for {num_sequences} sequence(s)")
    print(f"{'='*80}\n")

    for idx in range(num_sequences):
        seq_name = df.iloc[idx]["name"]
        print(f"\n--- Running sequence {idx}: {seq_name} ---")

        cmd = [sys.executable, docking_script, config_path, str(idx)]
        print(f"Command: {' '.join(cmd)}\n")

        result = subprocess.run(cmd)

        if result.returncode != 0:
            print(f"\nERROR: Sequence {idx} ({seq_name}) failed with return code {result.returncode}")
            sys.exit(1)

    print(f"\n{'='*80}")
    print(f"LOCAL TEST COMPLETE")
    print(f"{'='*80}")


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )
    parser.add_argument(
        "config_file",
        help="Path to config file with [SEQUENCE_DOCKING] section"
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Automatically submit to SLURM (default: just generate script and show command)"
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Preview sequences in CSV without generating submission script"
    )
    parser.add_argument(
        "--local",
        action="store_true",
        help="Run locally for testing (no SLURM)"
    )
    parser.add_argument(
        "--max-local",
        type=int,
        default=None,
        help="Maximum number of sequences to run locally (default: all)"
    )
    parser.add_argument(
        "--job-name",
        default="seq_dock",
        help="SLURM job name (default: seq_dock)"
    )
    parser.add_argument(
        "--time",
        default="24:00:00",
        help="SLURM time limit (default: 24:00:00)"
    )
    parser.add_argument(
        "--partition",
        default="amilan",
        help="SLURM partition (default: amilan)"
    )
    parser.add_argument(
        "--qos",
        default="normal",
        help="SLURM QoS (default: normal)"
    )
    parser.add_argument(
        "--mem",
        default="16G",
        help="Memory per task (default: 16G)"
    )
    parser.add_argument(
        "--conda-env",
        default="ligand_alignment",
        help="Conda environment to activate (default: ligand_alignment)"
    )
    parser.add_argument(
        "--output-dir",
        default="slurm_logs",
        help="Directory for SLURM log files (default: slurm_logs)"
    )

    args = parser.parse_args()

    # Load config
    config = ConfigParser()
    with open(args.config_file, "r", encoding="utf-8-sig") as f:
        config.read_file(f)

    # Get CSV path from config
    seq_section = config["SEQUENCE_DOCKING"] if "SEQUENCE_DOCKING" in config else {}
    csv_path = seq_section.get("SequenceCSV", None)

    if not csv_path:
        print("ERROR: SequenceCSV must be specified in [SEQUENCE_DOCKING] section of config")
        sys.exit(1)

    # Count sequences
    num_sequences, df = count_sequences_in_csv(csv_path)

    if args.preview:
        print("\nPreview mode: no submission script generated.")
        return

    # Local testing mode
    if args.local:
        run_local_test(args.config_file, args.max_local)
        return

    # Generate SLURM script
    slurm_script = generate_slurm_script(
        config_path=os.path.abspath(args.config_file),
        num_sequences=num_sequences,
        job_name=args.job_name,
        time=args.time,
        partition=args.partition,
        qos=args.qos,
        mem=args.mem,
        conda_env=args.conda_env,
        output_dir=args.output_dir,
    )

    # Create output directory for logs
    os.makedirs(args.output_dir, exist_ok=True)

    # Write SLURM script
    script_name = f"submit_{args.job_name}.sh"
    with open(script_name, "w") as f:
        f.write(slurm_script)

    print(f"\n{'='*80}")
    print(f"SLURM submission script generated: {script_name}")
    print(f"{'='*80}")
    print(f"Array size: 0-{num_sequences - 1} ({num_sequences} sequences)")
    print(f"Job name: {args.job_name}")
    print(f"Time limit: {args.time}")
    print(f"Partition: {args.partition}")
    print(f"Memory: {args.mem}")
    print(f"Log directory: {args.output_dir}/")
    print(f"{'='*80}\n")

    if args.submit:
        print("Submitting to SLURM...")
        result = subprocess.run(["sbatch", script_name], capture_output=True, text=True)
        print(result.stdout)
        if result.returncode != 0:
            print(f"ERROR: {result.stderr}")
            sys.exit(1)
        print("Submission complete!")
    else:
        print("To submit this job, run:")
        print(f"  sbatch {script_name}")
        print("\nOr re-run with --submit flag:")
        print(f"  python {sys.argv[0]} {args.config_file} --submit")


if __name__ == "__main__":
    main()
