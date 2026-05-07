# 5/6/2026
import os
import subprocess
import argparse
from pathlib import Path
import concurrent.futures
from tqdm import tqdm

# Set up argument parsing for flexible directories
parser = argparse.ArgumentParser(description="Convert a folder of PDBQT files to SMILES using Open Babel.")
parser.add_argument(
    "--input", 
    type=str, 
    default="./dock/zinc_id_ligands", 
    help="Path to the directory containing your PDBQT files."
)
parser.add_argument(
    "--output", 
    type=str, 
    default="./dock/zinc_id_smiles", 
    help="Path to the directory where SMILES files will be saved."
)
args = parser.parse_args()

# Directory paths based on arguments
ZINC_LIGAND_DIR = Path(args.input)
SMILES_DIR = Path(args.output)

# Verify input directory exists before starting
if not ZINC_LIGAND_DIR.exists():
    print(f"[!] Error: The input directory '{ZINC_LIGAND_DIR}' does not exist.")
    exit(1)

# Create the new directory for SMILES files
print(f"Targeting input directory: {ZINC_LIGAND_DIR}")
print(f"Targeting output directory: {SMILES_DIR}")
SMILES_DIR.mkdir(parents=True, exist_ok=True)

# Get the list of all PDBQT files
print("Getting list of PDBQT files... (This will take a few seconds)")
pdbqt_files = list(ZINC_LIGAND_DIR.glob("*.pdbqt"))

if not pdbqt_files:
    print(f"[!] No .pdbqt files found in '{ZINC_LIGAND_DIR}'.")
    exit(1)

print(f"Found {len(pdbqt_files)} files to convert.")

def convert_to_smiles(pdbqt_path):
    # .stem gets the filename without the extension (e.g., "ZINC12345")
    zinc_id = pdbqt_path.stem 
    output_path = SMILES_DIR / f"{zinc_id}.smi"
    
    # -----------------------------------------
    # BLUESCREEN PROTECTION (Crash Recovery)
    # -----------------------------------------
    if output_path.exists():
        return True
        
    try:
        # Construct the Open Babel command
        cmd = [
            "obabel", 
            "-i", "pdbqt", str(pdbqt_path), 
            "-o", "smi", 
            "-O", str(output_path)
        ]
        
        # Run the command in the background quietly
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return True
        
    except subprocess.CalledProcessError:
        # Catch errors if Open Babel fails on a specific, weirdly-formatted molecule
        return False
    except FileNotFoundError:
        # This triggers if obabel isn't installed or in your system PATH
        return "OBABEL_MISSING"

print("Converting PDBQT to SMILES using Open Babel...")

# Using 32 workers to process in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
    results = list(tqdm(executor.map(convert_to_smiles, pdbqt_files), total=len(pdbqt_files)))

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