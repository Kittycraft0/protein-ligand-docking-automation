# 5/7/2026
import os
import subprocess
import argparse
from pathlib import Path
import concurrent.futures
import threading
from tqdm import tqdm

# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert PDBQT files and compile them into a single master SMILES list.")
parser.add_argument(
    "--input", 
    type=str, 
    default="./dock/zinc_id_ligands", 
    help="Path to the directory containing your PDBQT files."
)
parser.add_argument(
    "--output", 
    type=str, 
    default="./dock/initial_zinc_smiles_full_list.txt", 
    help="The master output text file for your SMILES strings."
)
args = parser.parse_args()

ZINC_LIGAND_DIR = Path(args.input)
OUTPUT_FILE = Path(args.output)

# Verify input directory exists before starting
if not ZINC_LIGAND_DIR.exists():
    print(f"[!] Error: The input directory '{ZINC_LIGAND_DIR}' does not exist.")
    exit(1)

print(f"Targeting input directory: {ZINC_LIGAND_DIR}")
print(f"Outputting to master file: {OUTPUT_FILE}")

# -----------------------------------------
# BLUESCREEN PROTECTION (Memory Load)
# -----------------------------------------
processed_ids = set()
if OUTPUT_FILE.exists():
    print("Found existing master file. Loading previously processed IDs...")
    with open(OUTPUT_FILE, 'r') as f:
        for line in f:
            parts = line.strip().split()
            # The format is expected to be: SMILES \t ZINC_ID
            if len(parts) >= 2:
                processed_ids.add(parts[1])
    print(f"Loaded {len(processed_ids)} already processed molecules. These will be skipped.")

# Get the list of all PDBQT files
print("Getting list of PDBQT files... (This will take a few seconds)")
pdbqt_files = list(ZINC_LIGAND_DIR.glob("*.pdbqt"))

if not pdbqt_files:
    print(f"[!] No .pdbqt files found in '{ZINC_LIGAND_DIR}'.")
    exit(1)

print(f"Found {len(pdbqt_files)} total files.")

# Create a lock to prevent threads from writing over each other in the master file
write_lock = threading.Lock()

def convert_to_master_list(pdbqt_path):
    zinc_id = pdbqt_path.stem 
    
    # Skip if we already did this one before a crash
    if zinc_id in processed_ids:
        return True
        
    try:
        # Construct the Open Babel command, but DO NOT save to a file
        cmd = [
            "obabel", 
            "-i", "pdbqt", str(pdbqt_path), 
            "-o", "smi"
        ]
        
        # Capture the output directly into Python's memory
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Extract just the SMILES string (ignoring the file path obabel appends)
        smiles_string = result.stdout.split()[0] if result.stdout else None
        
        if smiles_string:
            # Wait for the lock, write to the file, then release the lock
            with write_lock:
                # Open in append mode ('a') so it adds to the bottom of the file
                with open(OUTPUT_FILE, "a") as f:
                    f.write(f"{smiles_string}\t{zinc_id}\n")
            return True
        else:
            return False
            
    except subprocess.CalledProcessError:
        return False
    except FileNotFoundError:
        return "OBABEL_MISSING"

print("Converting and compiling to master list...")

# Using 32 workers to process in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
    results = list(tqdm(executor.map(convert_to_master_list, pdbqt_files), total=len(pdbqt_files)))

# Check for specific obabel installation errors
if "OBABEL_MISSING" in results:
    print("\n[!] Error: 'obabel' command not found.")
    print("Make sure Open Babel is installed and added to your system's PATH.")
else:
    failures = results.count(False)
    if failures > 0:
        print(f"\nFinished with {failures} errors (files skipped due to formatting issues).")
    else:
        print("\nFinished successfully with 0 errors!")