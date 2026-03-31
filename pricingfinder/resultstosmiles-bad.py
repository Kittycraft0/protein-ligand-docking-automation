# 3/29/2026

import os
import glob
import subprocess
import time

# Define your directories and file names
input_dir = "results-12-12-2025-run3trial4"
output_dir = "results-12-12-2025-run3trial4-smiles"
master_file = "molport_bulk_search.txt"

# 1. Create the results-smiles folder if it doesn't already exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Created directory: {output_dir}")

# 2 & 3. Process files recursively and compile the master .txt document
with open(master_file, 'w') as master_f:
    
    # The ** wildcard and recursive=True tell glob to search all subfolders
    search_pattern = os.path.join(input_dir, "**", "*.pdbqt")
    pdbqt_files = glob.glob(search_pattern, recursive=True)
    
    total_files = len(pdbqt_files)

    if total_files == 0:
        print(f"No .pdbqt files found in '{input_dir}' or its subdirectories.")
    else:
        print(f"Found {total_files} files. Starting conversion...\n")
        
        start_time = time.time() # <-- ADD THIS HERE

        for i, pdbqt_file in enumerate(pdbqt_files):
            
            # Get the base filename (e.g., 'ligand_1' from 'ligand_1.pdbqt')
            base_name = os.path.splitext(os.path.basename(pdbqt_file))[0]
            
            # Append the unique index to prevent identically named files from overwriting each other
            smi_file = os.path.join(output_dir, f"{base_name}_{i}.smi")
            
            # Construct the OpenBabel command
            command = ["obabel", "-i", "pdbqt", pdbqt_file, "-o", "smi", "-O", smi_file]
            
            try:
                # Run the command quietly
                subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                
                # Open the newly generated .smi file to grab the text
                with open(smi_file, 'r') as sf:
                    line = sf.readline().strip()
                    if line:
                        # Isolate just the SMILES string and write it to the master text file
                        smiles_only = line.split()[0]
                        master_f.write(f"{smiles_only}\n")
                        
            except subprocess.CalledProcessError:
                print(f"[-] OpenBabel failed to convert: {pdbqt_file}. Skipping.")
            except FileNotFoundError:
                print("\n[!] Error: 'obabel' command not found.")
                print("Make sure OpenBabel is installed and added to your system's PATH.")
                break 
            
            current = i + 1
            percent = (current / total_files) * 100
            
            # ETA Math
            elapsed_time = time.time() - start_time
            time_per_file = elapsed_time / current
            remaining_files = total_files - current
            eta_seconds = time_per_file * remaining_files
            
            # Format ETA as HH:MM:SS
            eta_formatted = time.strftime('%H:%M:%S', time.gmtime(eta_seconds))
            
            # Visual bar
            filled_length = int(50 * current // total_files)
            bar = '█' * filled_length + '-' * (50 - filled_length)
            
            # Print with ETA appended
            print(f'\rProcessing: |{bar}| {percent:.1f}% ({current}/{total_files}) | ETA: {eta_formatted}', end='', flush=True)

        # Added a few newlines here so the final success message doesn't overwrite the progress bar
        print(f"\n\nSuccess! All SMILES compiled into '{master_file}'.")