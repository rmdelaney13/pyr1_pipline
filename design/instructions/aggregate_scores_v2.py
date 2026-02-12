import os
import argparse
import pandas as pd
import glob

def parse_kv_score_file(filepath):
    """
    Parses a score file where each line is 'Key: Value'.
    """
    data = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines or the header "SCORES:"
                if not line or line == "SCORES:":
                    continue
                
                # Split by the first colon only
                if ":" in line:
                    key, val = line.split(":", 1)
                    data[key.strip()] = val.strip()
        return data
    except Exception as e:
        print(f"Warning: Failed to parse {os.path.basename(filepath)}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Aggregate Key-Value style score files.")
    parser.add_argument("input_dir", help="Directory containing .sc files")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    args = parser.parse_args()

    # Find files
    sc_files = glob.glob(os.path.join(args.input_dir, "*.sc"))
    print(f"Found {len(sc_files)} score files in {args.input_dir}")
    
    if not sc_files:
        print("Error: No .sc files found.")
        return

    all_rows = []
    for f in sc_files:
        row_dict = parse_kv_score_file(f)
        if row_dict:
            # Force the filename to match the file on disk
            row_dict['filename'] = os.path.basename(f)
            all_rows.append(row_dict)
            
    if all_rows:
        final_df = pd.DataFrame(all_rows)
        
        # Reorder columns: filename first, then alphabetical
        cols = ['filename'] + sorted([c for c in final_df.columns if c != 'filename'])
        final_df = final_df[cols]
        
        final_df.to_csv(args.output, index=False)
        print(f"Success! Aggregated {len(final_df)} rows to {args.output}")
    else:
        print("Failed: No valid data extracted.")

if __name__ == "__main__":
    main()
