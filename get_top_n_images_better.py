# 12/12/2025
# Riley Mohr

import sys
import os
import shutil
import glob
import subprocess
import re
import atexit
import argparse
from pathlib import Path

# --- Global variable for temp directory cleanup ---
TEMP_DIR = Path("./temp_ligand_processing")

def cleanup():
    """Remove the temporary directory if it exists."""
    if TEMP_DIR.exists():
        try:
            shutil.rmtree(TEMP_DIR, ignore_errors=True)
        except Exception:
            pass

atexit.register(cleanup)


def die(message):
    """Prints a detailed error message and terminates the script."""
    print(f"\nFATAL ERROR: {message}", file=sys.stderr)
    print("Program execution has been terminated.", file=sys.stderr)
    sys.exit(1)


def run_command(command, ignore_errors=False):
    """Runs an external command."""
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, shell=True)
        return result
    except subprocess.CalledProcessError as e:
        if not ignore_errors:
            print(f"Warning: Command failed: {e.cmd}")
        return None


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process docking results, copy files, and generate 2D images.")
    parser.add_argument('-d', '--dir', type=Path, required=True, help="The input path to the base results directory.")
    parser.add_argument('-o', '--out', type=Path, default=Path("."), help="The output directory.")
    parser.add_argument('-n', '--num', type=int, default=100, help="Number of top ligands to process.")
    parser.add_argument('-s', '--scores', type=Path, default=None, help="Path to the scores file.")
    parser.add_argument('-v', '--vina', type=Path, default=Path("C:/Program Files (x86)/The Scripps Research Institute/Vina/vina_split.exe"), help="Path to vina_split.exe")
    return parser.parse_args()


def main():
    global TEMP_DIR
    args = parse_arguments()

    base_results_dir = args.dir
    output_base_dir = args.out
    num_ligands = args.num
    vina_split_exe = args.vina

    if args.scores:
        scores_file = args.scores
    else:
        scores_file = base_results_dir / "scores/best_ligands_overall.txt"

    top_n_dir = output_base_dir / f"top{num_ligands}"
    top_n_images_dir = output_base_dir / f"top{num_ligands}images"
    TEMP_DIR = output_base_dir / "temp_ligand_processing"

    print(f"Starting script to process top {num_ligands} ligands.")
    print(f"Input Directory:  {base_results_dir}")
    print(f"Output Directory: {output_base_dir.resolve()}")
    
    zfill_width = len(str(num_ligands))

    # --- Read Scores File ---
    if not scores_file.is_file():
        die(f"Scores file not found at: {scores_file}")

    top_results = []
    with open(scores_file, 'r') as f:
        lines = [line for line in f if not line.strip().startswith('#')]
        
        if len(lines) < num_ligands:
            num_ligands = len(lines)
            
        for i, line in enumerate(lines[:num_ligands]):
            parts = line.strip().split()
            if len(parts) < 2: 
                continue
            try:
                score = float(parts[0])
                name = parts[1]
                rank = i + 1
                top_results.append({'rank': rank, 'score': score, 'name': name})
            except ValueError:
                continue

    # --- Create Directories ---
    output_base_dir.mkdir(parents=True, exist_ok=True)
    top_n_dir.mkdir(exist_ok=True)
    top_n_images_dir.mkdir(exist_ok=True)
    TEMP_DIR.mkdir(exist_ok=True)

    # --- Part 1: Copying Ligand Files ---
    print("\n--- Part 1: Copying Ligand Files ---")
    
    search_root = base_results_dir / "docked_ligands"
    if not search_root.exists():
        search_root = base_results_dir

    for ligand in top_results:
        rank_str = str(ligand['rank']).zfill(zfill_width)
        name = ligand['name']
        score = ligand['score']
        
        # Format score consistently for folder name (6 decimal places)
        score_str = f"{score:.6f}"
        
        files_to_copy = []

        match_full = re.match(r'(.+?)\.xaa_(\d+)_model_(\d+)', name)
        match_short = re.match(r'(.+?)\.xaa_(\d+)', name)

        if match_full:
            cutname1 = match_full.group(1)
            model_num = match_full.group(3)
            cutname2 = f"{cutname1}_model{model_num}"
            cutname3 = f"{cutname1}.xaa_{match_full.group(2)}"
            
            p1 = search_root / f"docked_{cutname1}/docked_{cutname2}/{name}*"
            p2 = search_root / f"docked_{cutname3}/*"
            files_to_copy = glob.glob(str(p1)) + glob.glob(str(p2))

        if not files_to_copy:
             # Fallback: Recursive search for filename
             search_name = name
             found_files = list(search_root.rglob(f"{search_name}*"))
             files_to_copy = [str(f) for f in found_files if f.is_file()]

        dest_folder = top_n_dir / f"{rank_str} {score_str} {name}"
        dest_folder.mkdir(exist_ok=True)

        if not files_to_copy:
            print(f"Warning: Could not find ANY files for '{name}' (Rank {ligand['rank']}).")
            continue

        for f_path_str in files_to_copy:
            f_path = Path(f_path_str)
            try:
                shutil.copy(f_path, dest_folder)
            except Exception:
                pass

    # --- Part 2: Generating 2D Images ---
    print("\n--- Part 2: Generating 2D Images ---")
    if not vina_split_exe.is_file():
        die(f"Vina split executable not found at '{vina_split_exe}'.")
    if not shutil.which("obabel"):
        die("'obabel' command not found.")

    for ligand in top_results:
        rank_str = str(ligand['rank']).zfill(zfill_width)
        name = ligand['name']
        score = ligand['score']
        
        # Use consistent score formatting here too
        score_str = f"{score:.6f}"

        # Note: We must reconstruct the folder name exactly as we made it in Part 1
        source_folder = top_n_dir / f"{rank_str} {score_str} {name}"
        
        if not source_folder.is_dir():
             # If folder doesn't exist, maybe it failed copy or we used the float version in Part 1? 
             # (This ensures we look for the exact same string)
             continue

        pdbqt_files = list(source_folder.glob("*.pdbqt"))
        if not pdbqt_files:
            print(f"Skipping Rank {rank_str}: No .pdbqt files found in '{source_folder.name}'.")
            continue
        
        source_pdbqt = pdbqt_files[0]
        temp_pdbqt = TEMP_DIR / source_pdbqt.name
        shutil.copy(source_pdbqt, temp_pdbqt)

        # 1. Try to split (ignore errors if it's already a single model)
        split_command = f'"{vina_split_exe}" --input "{temp_pdbqt}"'
        run_command(split_command, ignore_errors=True)

        # 2. Determine target file
        split_files = [f for f in TEMP_DIR.glob("*.pdbqt") if f.name != temp_pdbqt.name]
        
        pdbqt_target = None
        if split_files:
            pdbqt_target = split_files[0]
        else:
            pdbqt_target = temp_pdbqt

        # 3. Standardize Image Name
        # We enforce specific naming: NAME_ID (e.g., CEABMN_0655)
        
        match_full = re.match(r'(.+?)\.xaa_(\d+)_model_\d+', name)
        match_short = re.match(r'(.+?)\.xaa_(\d+)', name)

        image_name_base = ""
        
        if match_full:
            # Matches: NAME.xaa_1234_model_1 -> NAME_1234
            image_name_base = f"{match_full.group(1)}_{match_full.group(2)}"
        elif match_short:
            # Matches: NAME.xaa_1234 -> NAME_1234
            image_name_base = f"{match_short.group(1)}_{match_short.group(2)}"
        else:
            # Fallback: Just use the raw name, but replace invalid chars just in case
            image_name_base = name.replace(".xaa_", "_").replace(".", "_")

        image_name_full = f"{rank_str} {score_str} {image_name_base}.png"
        image_output_path = top_n_images_dir / image_name_full

        # 4. Run Obabel
        print(f"Rank {rank_str}: Generating image '{image_output_path.name}'...")
        obabel_command = f'obabel "{pdbqt_target}" -o png -O "{image_output_path}"'
        run_command(obabel_command)

        # Clear temp folder
        for item in TEMP_DIR.iterdir():
            if item.is_file(): item.unlink()
            elif item.is_dir(): shutil.rmtree(item)

    print("\n--- Script Finished Successfully ---")


if __name__ == "__main__":
    main()




"""
How to use the new script
You must now provide the folder path when running the script. You can also override the other settings if you want.

1. Basic Usage (Required): This uses the default scores file location and default number of ligands (100).

Bash

python script.py -d "results-7-20-2025-run2trial3/results"
2. Changing the Number of Ligands: To process the top 50 instead of 100:

Bash

python script.py -d "results-7-20-2025-run2trial3/results" -n 50
3. Specifying a different Scores File: If your scores file is in a weird place:

Bash

python script.py -d "results-folder" -s "C:/somewhere/else/my_scores.txt"
4. Specifying a different Vina Executable: If you run this on a different computer where Vina is installed somewhere else:

Bash

python script.py -d "results-folder" -v "C:/Apps/Vina/vina_split.exe"
5. Need Help? You can now run this to see all available options:

Bash

python script.py --help
"""