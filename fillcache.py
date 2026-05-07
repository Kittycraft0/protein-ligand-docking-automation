# 5/1/2026
# writing it myself from connecting things from other scripts because the ai kinda sucks
# lmao i was copying it over from dock.py and then i realized 
# "but wait, it's all compartmentalized! i can just import a module and run it!" 
# and so that's what i'm doing here
from pathlib import Path
# Directory paths
BASE_DIR = Path("./dock")
LIGAND_DIR = BASE_DIR / "ligands"
CACHE_DIR = BASE_DIR / "cache"
LIGAND_NAMES = CACHE_DIR / "ligandNames.txt"
CONFIG_DIR = BASE_DIR / "config"
RESULTS_DIR = BASE_DIR / "results"
PROTEIN_DIR = BASE_DIR / "proteins"
COMPARISON_LIGAND_DIR = BASE_DIR / "comparison_ligands"

from dock import LigandManager
from dock import DisplayManager
from dock import CacheManager

cache_manager = CacheManager(CACHE_DIR, RESULTS_DIR, LIGAND_DIR, PROTEIN_DIR, COMPARISON_LIGAND_DIR)
ligand_manager = LigandManager(LIGAND_DIR, CACHE_DIR)
display_manager=DisplayManager()


import sys
# Handle command-line arguments for clearing cache
#cache_manager = CacheManager(CACHE_DIR, RESULTS_DIR, LIGAND_DIR, PROTEIN_DIR, COMPARISON_LIGAND_DIR)
if "--clear-cache" in sys.argv:
    cache_manager.initialize_cache("clear")
elif "--clear-everything" in sys.argv:
    cache_manager.initialize_cache("clear-everything")
else:
    cache_manager.initialize_cache()




# Load ligands, proteins, and comparison ligands
ligands = [line.strip() for line in open(LIGAND_NAMES)]

total_ligand_size = sum(Path(f).stat().st_size for f in ligands) if ligands else 0
processed_size = 0
import time
extraction_start_time = time.time()

for ligand_file in ligands:
    current_file_path = Path(ligand_file)
    current_file_size = current_file_path.stat().st_size
    # We create a dictionary of arguments to pass to the display function
    display_args = {
        "title": "Step 1: Preparing Ligand Models",
        "processed_size": processed_size,
        "total_size": total_ligand_size,
        "item_name": current_file_path.name,
        "start_time": extraction_start_time,
        "current_file_size": current_file_size,
    }
    # The extract_models function will now call the display function internally
    # using the 'display_callback' argument we added to it.
    ligand_manager.extract_models(
        ligand_file,
        display_callback=display_manager.display_extraction_progress,
        display_args=display_args
    )
    # YOU NEED THIS LINE HERE:
    processed_size += current_file_size

# Final 100% display after the loop is done
display_manager.display_extraction_progress(
    title="Step 1: Preparing Ligand Models",
    processed_size=total_ligand_size,
    total_size=total_ligand_size,
    item_name="All models prepared.",
    start_time=extraction_start_time,
    current_file_size=0,
    current_file_progress=1.0 # Show 100% for the "current file" bar
)