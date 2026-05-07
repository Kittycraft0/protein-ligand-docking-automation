# 5/6/2026
# create folder dock/zincidligands if it doesn't already exist
from pathlib import Path
# Directory paths
BASE_DIR = Path("./dock")
PROTEIN_DIR = BASE_DIR / "proteins"
LIGAND_DIR = BASE_DIR / "ligands"
CONFIG_DIR = BASE_DIR / "config"
CACHE_DIR = BASE_DIR / "cache"
RESULTS_DIR = BASE_DIR / "results"
COMPARISON_LIGAND_DIR = BASE_DIR / "comparison_ligands"

ZINC_LIGAND_DIR = BASE_DIR / "zinc_id_ligands"

# Cache files
LIGAND_NAMES = CACHE_DIR / "ligandNames.txt"
PROTEIN_NAMES = CACHE_DIR / "proteinNames.txt"
PROGRESS_CACHE = CACHE_DIR / "progress_cache.txt"
COMPARISON_LIGAND_NAMES = CACHE_DIR / "comparisonLigandNames.txt"

# Ensure necessary directories exist
for directory in [ZINC_LIGAND_DIR]:
    print(f"Creating directory {directory}")
    directory.mkdir(parents=True, exist_ok=True)

import sys
import shutil
# add --clear-everything option to empty dock/zincligands with progress bar before running
if "--clear-everything" in sys.argv:
    # TODO: add progress bar to rmtree
    print("Clearing all ligands in \"dock/zinc_id_ligands\"...")
    shutil.rmtree(ZINC_LIGAND_DIR, ignore_errors=True)
    ZINC_LIGAND_DIR.mkdir(parents=True, exist_ok=True)

#print(list(LIGAND_DIR.glob("*.pdbqt")))
# find ligands in cache, so:
print("Getting ligand directories... (5-10 seconds)")
ligand_dirs=list(CACHE_DIR.rglob("*.pdbqt"))
print(f"Received {len(ligand_dirs)} ligand directories!")
#print(list(CACHE_DIR.rglob("*.pdbqt"))) #do NOT print, it will take a solid 30 seconds

import os

## ai:
#def extract_pdbqt_name(filepath):
#    """
#    Extracts the molecule name from a PDBQT file.
#    
#    Args:
#        filepath (str): The full path to the file.
#        
#    Returns:
#        str: The extracted molecule name (e.g., 'ZINC000065090419').
#        
#    Raises:
#        FileNotFoundError: If the file path is invalid.
#        ValueError: If the 'Name =' string is missing from the file.
#    """
#    # Check if the file actually exists before trying to open it
#    if not os.path.exists(filepath):
#        raise FileNotFoundError(f"Error: The file at '{filepath}' was not found.")
#
#    # Open and read the file line by line
#    with open(filepath, 'r', encoding='utf-8') as file:
#        for line in file:
#            # Look for the specific REMARK line containing the name
#            if line.startswith("REMARK") and "Name =" in line:
#                # Split at 'Name =' and take everything to the right, stripping whitespace
#                molecule_name = line.split("Name =")[1].strip()
#                return molecule_name
#                
#    # If the loop finishes without returning, the name wasn't in the file
#    raise ValueError(f"Error: Molecule name not found in '{filepath}'.")
#
## mine again:
#from tqdm import tqdm
## for every ligand in cache (recursive search through folders ig idek) (with progress bar):
#for ligand_dir in tqdm(ligand_dirs):
#    # get zinc id
#    filename=extract_pdbqt_name(ligand_dir)
#    # extract zinc id from file
#    # copy to dock/zincligands/[zinc id].pdbqt
#    pass

import os
import shutil
import concurrent.futures
from tqdm import tqdm

def extract_pdbqt_name(filepath):
    # Check if the file actually exists before trying to open it
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Error: The file at '{filepath}' was not found.")

    # Open and read the file line by line
    with open(filepath, 'r', encoding='utf-8') as file:
        for line in file:
            # Look for the specific REMARK line containing the name
            if line.startswith("REMARK") and "Name =" in line:
                # Split at 'Name =' and take everything to the right, stripping whitespace
                return line.split("Name =")[1].strip()
                
    # If the loop finishes without returning, the name wasn't in the file
    raise ValueError(f"Error: Molecule name not found in '{filepath}'.")

def process_ligand(ligand_path):
    try:
        # Extract the ID
        zinc_id = extract_pdbqt_name(ligand_path)
        target_path = ZINC_LIGAND_DIR / f"{zinc_id}.pdbqt"
        
        # Skip if it already exists
        if target_path.exists():
            return True
            
        # Try to hardlink for instant "copying" (zero disk space used)
        try:
            os.link(ligand_path, target_path)
        except OSError:
            # Fallback to standard copy if hardlinking fails (e.g., cross-drive transfers)
            shutil.copy(ligand_path, target_path)
            
        return True
    
    except Exception:
        # Catch extraction errors or missing names
        return False

# Execute the tasks in parallel
print("Processing and linking ligands...")

# 32 max_workers is a solid sweet spot for heavy I/O operations
with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
    # Wrap the executor map in tqdm. Because we know the length of ligand_dirs, 
    # tqdm will give a highly accurate progress bar and ETA.
    results = list(tqdm(executor.map(process_ligand, ligand_dirs), total=len(ligand_dirs)))

# Optional: Check for failures
failures = results.count(False)
if failures > 0:
    print(f"Finished with {failures} errors (ligands skipped).")
else:
    print("Finished successfully with 0 errors!")