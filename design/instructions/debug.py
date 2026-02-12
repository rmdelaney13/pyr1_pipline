import pandas as pd
import os

# --- UPDATE THESE PATHS IF NEEDED ---
csv_path = "/scratch/alpine/ryde3462/CA_design_normal/relax_threaded_2/relax_2.csv"
relax_dir = "/scratch/alpine/ryde3462/CA_design_normal/relax_threaded_2"
# ------------------------------------

if not os.path.exists(csv_path):
    print(f"Error: CSV not found at {csv_path}")
else:
    df = pd.read_csv(csv_path)
    print(f"Loaded CSV with {len(df)} rows.")
    
    # Check the first 5 entries
    print("\n--- DEBUGGING FILENAMES ---")
    for i, row in df.head(5).iterrows():
        csv_filename = str(row['filename'])
        
        # Mimic the filtering script logic
        base = os.path.basename(csv_filename).strip()
        stem, _ext = os.path.splitext(base)
        if stem.endswith("_score"): stem = stem[:-len("_score")]
        if stem.endswith("_relaxed_score"): stem = stem[:-len("_relaxed_score")]
        
        expected_pdb = f"{stem}.pdb"
        full_path = os.path.join(relax_dir, expected_pdb)
        found = os.path.exists(full_path)
        
        print(f"\n1. CSV says:      {csv_filename}")
        print(f"2. Script looks for: {expected_pdb}")
        print(f"3. Found on disk?    {'YES' if found else 'NO'}")
