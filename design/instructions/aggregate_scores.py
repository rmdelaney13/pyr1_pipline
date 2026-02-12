#!/usr/bin/env python3

import os
import csv
import argparse

def aggregate_scores(input_dir, output_csv):
    rows = []
    all_keys = set()

    # Verify input directory exists
    if not os.path.exists(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        return

    print(f"Scanning {input_dir} for .sc files...")

    for fname in os.listdir(input_dir):
        # CHANGED: Look for .sc instead of _score.sc to match relax_threaded.py output
        if fname.endswith(".sc"):
            try:
                with open(os.path.join(input_dir, fname)) as f:
                    data = {}
                    for line in f:
                        line = line.strip()
                        if line.startswith("SCORES:") or not line:
                            continue
                        if ":" in line:
                            key, value = line.split(":", 1)
                            data[key.strip()] = value.strip()

                    # Always include the filename
                    data["filename"] = fname

                    # Compute the new "all_polar_and_charge" column:
                    #  - O1_polar_contact must be "yes"
                    #  - O2_polar_contact must be "yes"
                    #  - (Optional) charge_satisfied must be "1"
                    if (
                        data.get("O1_polar_contact", "").lower() == "yes"
                        and data.get("O2_polar_contact", "").lower() == "yes"
                        # Uncomment the line below if you want to enforce charge satisfaction too
                        # and data.get("charge_satisfied", "") == "1"
                    ):
                        data["all_polar_and_charge"] = "yes"
                    else:
                        data["all_polar_and_charge"] = "no"

                    rows.append(data)
                    all_keys.update(data.keys())
            except Exception as e:
                print(f"Warning: Could not read {fname}: {e}")

    if not rows:
        print("No .sc files found.")
        return

    # Ensure our new column is included in the header
    fieldnames = ["filename"] + sorted(k for k in all_keys if k not in {"filename"})

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Aggregated {len(rows)} score files into {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Aggregate Rosetta-style .sc files into a CSV."
    )
    parser.add_argument(
        "input_dir", help="Directory containing .sc files"
    )
    parser.add_argument(
        "--output",
        default="aggregated_scores.csv",
        help="Output CSV file name (default: aggregated_scores.csv)",
    )
    args = parser.parse_args()

    aggregate_scores(args.input_dir, args.output)
