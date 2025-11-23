import sys
import os
import shutil
import glob
import subprocess
import re
from pathlib import Path
import atexit

# --- Configuration ---
# You can modify these paths if needed.
VINA_SPLIT_EXE = Path("C:/Program Files (x86)/The Scripps Research Institute/Vina/vina_split.exe")
BASE_RESULTS_DIR = Path("results-7-20-2025-run2trial3/results")
SCORES_FILE = BASE_RESULTS_DIR / "scores/best_ligands_overall.txt"

# --- Global variable for temp directory cleanup ---
TEMP_DIR = Path("./temp_ligand_processing")

def cleanup():
    """Remove the temporary directory if it exists."""
    if TEMP_DIR.exists():
        print(f"\nCleaning up and deleting temporary directory: {TEMP_DIR}")
        shutil.rmtree(TEMP_DIR, ignore_errors=True)

# Register the cleanup function to be called on script exit
atexit.register(cleanup)


def die(message):
    """Prints a detailed error message and terminates the script."""
    print(f"\nFATAL ERROR: {message}", file=sys.stderr)
    print("Program execution has been terminated.", file=sys.stderr)
    sys.exit(1)


def run_command(command):
    """Runs an external command and handles errors."""
    try:
        # Using shell=True for convenience with complex commands and paths
        result = subprocess.run(command, check=True, capture_output=True, text=True, shell=True)
        return result
    except FileNotFoundError:
        die(f"Command not found. Please ensure '{command.split()[0]}' is installed and in your PATH.")
    except subprocess.CalledProcessError as e:
        error_message = f"An external command failed:\n"
        error_message += f"  Command: {e.cmd}\n"
        error_message += f"  Exit Code: {e.returncode}\n"
        error_message += f"  Stdout: {e.stdout.strip()}\n"
        error_message += f"  Stderr: {e.stderr.strip()}"
        die(error_message)


def main():
    """Main function to orchestrate the ligand processing workflow."""
    # --- 1. Initialization and Argument Parsing ---
    try:
        num_ligands = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    except ValueError:
        die("Invalid argument. The number of ligands must be an integer.")

    print(f"Starting script to process top {num_ligands} ligands.")
    zfill_width = len(str(num_ligands))

    # Define output directories based on N
    top_n_dir = Path(f"top{num_ligands}")
    top_n_images_dir = Path(f"top{num_ligands}images")
    
    # --- 2. Read and Parse Scores File ---
    if not SCORES_FILE.is_file():
        die(f"Scores file not found at the expected location: {SCORES_FILE}")

    print(f"Reading scores from: {SCORES_FILE}")
    top_results = []
    with open(SCORES_FILE, 'r') as f:
        lines = [line for line in f if not line.strip().startswith('#')]
        
        if len(lines) < num_ligands:
            die(f"Scores file contains only {len(lines)} data lines, but {num_ligands} were requested.")
            
        for i, line in enumerate(lines[:num_ligands]):
            parts = line.strip().split()
            if len(parts) != 2:
                print(f"Warning: Skipping malformed line {i+1}: '{line.strip()}'")
                continue
            
            try:
                score = float(parts[0])
                name = parts[1]
                rank = i + 1
                top_results.append({'rank': rank, 'score': score, 'name': name})
            except ValueError:
                print(f"Warning: Skipping line {i+1} due to invalid score: '{line.strip()}'")

    print(f"Successfully parsed {len(top_results)} ligand records.")

    # --- 3. Create Destination Directories ---
    top_n_dir.mkdir(exist_ok=True)
    top_n_images_dir.mkdir(exist_ok=True)
    TEMP_DIR.mkdir(exist_ok=True)
    print(f"Output folders '{top_n_dir}' and '{top_n_images_dir}' are ready.")
    print(f"Temporary folder '{TEMP_DIR}' created.")


    # --- 4. Copy Ligand Files ---
    print("\n--- Part 1: Copying Ligand Files ---")
    for ligand in top_results:
        rank_str = str(ligand['rank']).zfill(zfill_width)
        name = ligand['name']
        score = ligand['score']

        # Use regex to parse ligand name components
        match = re.match(r'(.+?)\.xaa_(\d+)_model_(\d+)', name)
        if not match:
            die(f"Ligand name '{name}' (Rank {ligand['rank']}) does not match the expected pattern 'NAME.xaa_####_model_#'.")
        
        cutname1 = match.group(1)
        cutname2 = f"{cutname1}_model{match.group(3)}"
        cutname3 = f"{cutname1}.xaa_{match.group(2)}"

        # Define destination folder
        dest_folder = top_n_dir / f"{rank_str} {score} {name}"
        dest_folder.mkdir(exist_ok=True)

        # Define source paths
        source_pattern_1 = BASE_RESULTS_DIR / f"docked_ligands/docked_{cutname1}/docked_{cutname2}/{name}*"
        source_pattern_2 = BASE_RESULTS_DIR / f"docked_ligands/docked_{cutname3}/*"
        
        files_to_copy = glob.glob(str(source_pattern_1)) + glob.glob(str(source_pattern_2))
        
        if not files_to_copy:
            print(f"Warning: No source files found for ligand '{name}' (Rank {ligand['rank']}) using patterns:\n  1. {source_pattern_1}\n  2. {source_pattern_2}")

        # Copy files
        files_copied_count = 0
        for f_path_str in files_to_copy:
            f_path = Path(f_path_str)
            shutil.copy(f_path, dest_folder)
            files_copied_count += 1
        
        print(f"Rank {rank_str}: Copied {files_copied_count} files for '{name}' to '{dest_folder.name}'.")

    # --- 5. Generate 2D Images ---
    print("\n--- Part 2: Generating 2D Images ---")
    if not VINA_SPLIT_EXE.is_file():
        die(f"Vina split executable not found at '{VINA_SPLIT_EXE}'. Please check the path.")
    if not shutil.which("obabel"):
        die("'obabel' command not found. Please ensure Open Babel is installed and its command is in your system's PATH.")

    for ligand in top_results:
        rank_str = str(ligand['rank']).zfill(zfill_width)
        name = ligand['name']
        score = ligand['score']

        # Find the source folder from Part 1
        source_folder = top_n_dir / f"{rank_str} {score} {name}"
        if not source_folder.is_dir():
             die(f"The ligand folder for Rank {ligand['rank']} ('{source_folder}') does not exist. Cannot generate image.")

        # Find a .pdbqt file in the folder
        pdbqt_files = list(source_folder.glob("*.pdbqt"))
        if not pdbqt_files:
            die(f"No .pdbqt files found in folder '{source_folder}' for ligand Rank {ligand['rank']}.")
        
        # Take the first one as requested
        source_pdbqt = pdbqt_files[0]
        
        # Copy to temp folder for processing
        temp_pdbqt = TEMP_DIR / source_pdbqt.name
        shutil.copy(source_pdbqt, temp_pdbqt)

        # Run vina_split.exe
        print(f"Rank {rank_str}: Splitting {source_pdbqt.name} with vina_split...")
        split_command = f'"{VINA_SPLIT_EXE}" --input "{temp_pdbqt}"'
        run_command(split_command)

        # Find one of the newly generated single-conformation files
        split_files = [f for f in TEMP_DIR.glob("*.pdbqt") if f.name != temp_pdbqt.name]
        if not split_files:
            die(f"vina_split did not produce any output files in '{TEMP_DIR}' for Rank {ligand['rank']}.")
        
        pdbqt_single_file = split_files[0]

        # Generate the final image name
        match = re.match(r'(.+?)\.xaa_(\d+)_model_\d+', name)
        # We already validated the pattern, so we can assume it matches here
        image_name_base = f"{match.group(1)}{match.group(2)}"
        image_name_full = f"{rank_str} {score} {image_name_base}.png"
        image_output_path = top_n_images_dir / image_name_full

        # Run obabel command to generate the image
        print(f"Rank {rank_str}: Generating image '{image_output_path.name}' with obabel...")
        obabel_command = f'obabel "{pdbqt_single_file}" -o png -O "{image_output_path}"'
        run_command(obabel_command)

        # Clear the temp folder for the next iteration
        for item in TEMP_DIR.iterdir():
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)

    print("\n--- Script Finished Successfully ---")


if __name__ == "__main__":
    main()