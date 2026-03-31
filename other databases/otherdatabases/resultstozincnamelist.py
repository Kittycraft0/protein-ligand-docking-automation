# 3/29/2026

import os
import glob
import time
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Extract ZINC IDs from AutoDock Vina pdbqt files.")
parser.add_argument("input_dir", help="The folder containing your .pdbqt files")

# Added nargs='?' to make it optional, and defined the default filename
parser.add_argument("output_file", nargs='?', default="zinc_ids_extracted.txt", 
                    help="The output text file name (defaults to 'zinc_ids_extracted.txt')")

args = parser.parse_args()

# Map the arguments to your existing variables
input_dir = args.input_dir
master_file = args.output_file

# 1. Find all files recursively
search_pattern = os.path.join(input_dir, "**", "*.pdbqt")
pdbqt_files = glob.glob(search_pattern, recursive=True)
total_files = len(pdbqt_files)

if total_files == 0:
    print(f"No .pdbqt files found in '{input_dir}' or its subdirectories.")
else:
    print(f"Found {total_files} files. Extracting ZINC IDs...\n")
    
    unique_zinc_ids = set()
    start_time = time.time()

    for i, pdbqt_file in enumerate(pdbqt_files):
        try:
            # Open and read the file line by line
            with open(pdbqt_file, 'r') as f:
                for line in f:
                    # Look for the target line
                    if line.startswith("REMARK  Name ="):
                        # Split by the equals sign and strip away whitespace
                        # Turns "REMARK  Name = ZINC000005530195" into "ZINC000005530195"
                        zinc_id = line.split("=")[1].strip()
                        unique_zinc_ids.add(zinc_id)
                        
                        # Break out of the loop early since we found the ID
                        break 
        except Exception:
            # Silently pass if a file is corrupted or unreadable
            pass

        # --- PROGRESS BAR WITH ETA ---
        current = i + 1
        percent = (current / total_files) * 100
        
        elapsed_time = time.time() - start_time
        time_per_file = elapsed_time / current
        remaining_files = total_files - current
        eta_seconds = time_per_file * remaining_files
        
        eta_formatted = time.strftime('%H:%M:%S', time.gmtime(eta_seconds))
        
        filled_length = int(50 * current // total_files)
        bar = '█' * filled_length + '-' * (50 - filled_length)
        
        print(f'\rProcessing: |{bar}| {percent:.1f}% ({current}/{total_files}) | ETA: {eta_formatted}', end='', flush=True)

    # 2. Write the unique IDs to the final master text file
    with open(master_file, 'w') as out_f:
        for z_id in unique_zinc_ids:
            out_f.write(f"{z_id}\n")

    print(f"\n\nSuccess! Extracted {len(unique_zinc_ids)} unique ZINC IDs into '{master_file}'.")