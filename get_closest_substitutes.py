# 11/22/2025
# Riley Mohr



# Given that two iterations have now passed
# Reading off of 
# "D:\Users\iwbmo\Desktop\work\dockingprogram\protein-ligand-docking-automation\dock\results\ranked_best_ligands.txt"
# or ".\dock\results\ranked_best_ligands.txt"
# Finding the top so many that are within so many units (inputted number) from ABA
# The list is ranked by distance from the score of the target molecule to the closeness score by RMS I think
# Essentially just like

# make a tempt folder for unpacking the model files. make sure it does not exist before running, 
# and remove itafter program execution. 
# if a folder with the same name does exist, then add a _n after it, where n is a number.
# make n be an integer that goes up by 1 for each time it fails to find a folder without the same name.
# so then the folder is like temp_get_closest_substitutes_py_[timestamp]_[n] or something.
# while the number at the beginning of the line is smaller than or equal to the inputted number:
#   get the name on a line
#   find the pdbqt file associated with the name
#       should be located in... complex naming scheme below
#   that file is a bunch of the models and needs to be extracted, so that's why a temp folder for those was made


# e.g. line has this:
# 0.00000000 BCABMM.xaa_085_model_1
# then the score is 0.00000000
# and the name is BCABMM.xaa_085_model_1
# name="BCABMM.xaa_085_model_1"
# we may get a substring of the text to get "BCABMM"
# add "docked_" to the beginning to get "docked_BCAMBB", store that in another variable, call it folder1.
# get the substring of the name to get "_model_1", store that in a string variable called modelnumber.
# folder2 = folder1 + modelnumber (string addition)
# ligandfolder = "./dock/results/docked_ligands/folder1/folder2"
# search the folder for the file that starts with name and ends with ".pdbqt"
# there might be multiple, this would signify multiple proteins; 
# the data is technically different in spacial orientation, but orientation doesn't matter--just the model itself. 
# so just pick one. could be the first, who cares, so long as it ends in ".pdbqt" and isn't a .log file.

# once that file has been found, it has a long name. copy it to the temp folder for model extration, 
# renaming the file to for example BCABMM.xaa_085 in the process of copying it.
# then run the vina extraction line thing, probably can be foudn in dock.py, to extract the model in the folder.
# i forget what the model files are called post exptraion, probably liek _model_n where n is an integer 
# to signify which extracted model the file is. take the first one and move it to the inputted output folder,
# renaming the file to for example BCABMM.xaa_085 in the process of copying it.

# and ofc make a progress bar for it.



"""import os
import shutil
import random
import string
import subprocess
from tqdm import tqdm
import sys


def processFiles(ranked_best_ligands_text,output_folder_path):
    # starting
    print("Getting closest substitutes file started")

    #

    # check if path exists
    apath="./dock"
    print(os.path.exists(apath))
    
    print(f"argument 1: {sys.argv[1]}")


def main():
    # 1. Check for command-line argument for the input file
    if len(sys.argv) != 3 and len(sys.argv) != 2:
        # Provide usage instructions if the argument is missing
        print("Usage: python get_closest_substitutes.py " \
        "<output folder path and name> " \
        "<input path to ranked_best_ligands.txt, defaults to \".\\dock\\results\\ranked_best_ligands.txt\">")
        sys.exit(1) # Exit the script indicating an error
    
    if len(sys.argv) == 3:
        ranked_best_ligands_path = sys.argv[2]
    elif len(sys.argv)==2:
        ranked_best_ligands_path="./dock/results/ranked_best_ligands.txt"
        # make sure user didn't input this for the output directory
        if sys.argv[1]=="./dock/results/ranked_best_ligands.txt" or sys.argv[1]==".\\dock\\results\\ranked_best_ligands.txt":
            print("I don't think you have your inputs right; you're supposed to put your output directory there," \
            "not your input directory")
            sys.exit(1)

    output_dir_path=sys.argv[2]

    #print(f"sys.argv: {sys.argv}")

    # 2. Read the content from the specified input file
    try:
        with open(ranked_best_ligands_path, 'r') as f:
            ranked_best_ligands_text = f.read()
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: The file '{ranked_best_ligands_path}' was not found.")
        sys.exit(1)
    
    processFiles(ranked_best_ligands_text,output_dir_path)


if __name__ == "__main__":
    main()"""



import os
import shutil
import time
import subprocess
import argparse
from tqdm import tqdm

# --- CONFIGURATION ---
# Path to the results file
RESULTS_FILE = r".\dock\results\ranked_best_ligands.txt"
# Base directory where the docking folders are located
DOCK_BASE_DIR = r".\dock\results\docked_ligands"
# The command to run vina_split.
VINA_SPLIT_CMD = "vina_split" 

# Global debug flag
DEBUG_MODE = False
# Global counter for log categories
LOG_COUNTS = {}

def log(msg, category=None, limit=50):
    """
    Helper to print only if debug mode is on.
    Suppresses output if a specific category has been logged 'limit' times.
    """
    if DEBUG_MODE:
        if category:
            count = LOG_COUNTS.get(category, 0)
            LOG_COUNTS[category] = count + 1
            if count == limit:
                print(f"[DEBUG] ... (Limit reached for '{category}' logs. Suppressing further output.)")
                return
            elif count > limit:
                return
        print(f"[DEBUG] {msg}")

def create_temp_folder():
    timestamp = str(int(time.time()))
    n = 0
    while True:
        folder_name = f"temp_extract_{timestamp}_{n}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
            log(f"Created temp folder: {os.path.abspath(folder_name)}")
            return folder_name
        n += 1

def get_pdb_code(name_string):
    return name_string.split('.')[0]

def get_model_suffix(name_string):
    split_marker = "_model_"
    if split_marker in name_string:
        index = name_string.rfind(split_marker)
        # Fix: Change "_model_1" to "_model1" to match actual folder structure
        return name_string[index:].replace("_model_", "_model")
    return ""

def get_clean_name(name_string):
    return name_string.split('_model_')[0]

def main():
    global DEBUG_MODE

    # 1. Set up Argument Parsing
    parser = argparse.ArgumentParser(description="Extract top docked ligands based on score distance.")
    parser.add_argument("cutoff", type=float, help="The maximum score distance")
    parser.add_argument("output_dir", type=str, help="The folder where extracted models will be saved")
    parser.add_argument("--debug", action="store_true", help="Enable detailed logging for troubleshooting")

    args = parser.parse_args()
    
    cutoff_score = args.cutoff
    output_dir = args.output_dir
    DEBUG_MODE = args.debug

    if DEBUG_MODE:
        print("\n*** DEBUG MODE ENABLED ***")
        print(f"Current Working Directory: {os.getcwd()}")
        print(f"Looking for results file at: {os.path.abspath(RESULTS_FILE)}")
        print(f"Looking for dock folders at: {os.path.abspath(DOCK_BASE_DIR)}")
        print("------------------------------------------------\n")
        print("Note: Repeating debug logs will be suppressed after 50 occurrences.")

    # Create output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        log(f"Created output directory: {os.path.abspath(output_dir)}")

    # 2. Setup Temp Directory
    temp_dir = create_temp_folder()

    # 3. Read and Filter
    if not os.path.exists(RESULTS_FILE):
        print(f"CRITICAL ERROR: Could not find results file: {os.path.abspath(RESULTS_FILE)}")
        return

    tasks = []
    
    log("Reading results file line by line...")
    with open(RESULTS_FILE, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            
            try:
                score = float(parts[0])
                name = parts[1]
            except ValueError:
                continue

            if score <= cutoff_score:
                tasks.append((score, name))
            else:
                # log(f"Skipped {name} (Score {score} > {cutoff_score})")
                pass

    print(f"Found {len(tasks)} ligands to extract.")

    # 4. Process Files
    # If debug is on, we disable the tqdm progress bar so it doesn't mess up the print logs
    iterator = tasks if DEBUG_MODE else tqdm(tasks, desc="Extracting", unit="file")

    for score, name in iterator:
        if DEBUG_MODE:
            log(f"\n--- Processing: {name} ---", category="header")

        # Construct Path logic
        pdb_code = get_pdb_code(name)
        folder1 = f"docked_{pdb_code}"
        model_suffix = get_model_suffix(name)
        folder2 = f"{folder1}{model_suffix}"
        
        target_folder = os.path.join(DOCK_BASE_DIR, folder1, folder2)
        
        # Combine path details into one log to count it as one "event"
        path_msg = (f"Path construction details:\n"
                    f"  1. PDB Code: '{pdb_code}'\n"
                    f"  2. Folder 1: '{folder1}'\n"
                    f"  3. Folder 2: '{folder2}'\n"
                    f"  4. Full Search Path: '{os.path.abspath(target_folder)}'")
        log(path_msg, category="path_details")

        # Find the .pdbqt file
        found_file_path = None
        if os.path.exists(target_folder):
            files_in_dir = os.listdir(target_folder)
            # log(f"  Files found in folder: {files_in_dir}") # Uncomment if you suspect file naming issues
            
            for file in files_in_dir:
                if file.startswith(name) and file.endswith(".pdbqt") and not file.endswith(".log"):
                    found_file_path = os.path.join(target_folder, file)
                    log(f"  MATCH FOUND: {file}", category="match_found")
                    break
        else:
            log(f"  !! FOLDER NOT FOUND !! The directory does not exist.", category="folder_error")
        
        if not found_file_path:
            log(f"  SKIPPING: Could not locate source .pdbqt file for {name}", category="skip_missing")
            continue

        # Copy to Temp & Extract
        clean_name = get_clean_name(name)
        temp_pdbqt_name = f"{clean_name}.pdbqt"
        temp_file_path = os.path.join(temp_dir, temp_pdbqt_name)

        shutil.copy(found_file_path, temp_file_path)
        log(f"  Copied to temp: {temp_file_path}", category="copy_action")

        # Run Vina Split
        # Updated: Always capture output to prevent console spam, even in debug mode.
        log("  Running vina_split...", category="vina_run")
        try:
            result = subprocess.run([VINA_SPLIT_CMD, "--input", temp_pdbqt_name], 
                           cwd=temp_dir, 
                           stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE,
                           text=True)
            
            if result.returncode != 0 and DEBUG_MODE:
                print(f"  [ERROR] vina_split failed for {name}:")
                print(result.stderr)

        except FileNotFoundError:
            print("CRITICAL ERROR: 'vina_split' command not found. Is it installed or in your PATH?")
            break

        # Move the first model to Output
        extracted_model = None
        files_in_temp = os.listdir(temp_dir)
        
        # log(f"  Files in temp after split: {files_in_temp}")

        for f in files_in_temp:
            if f.startswith(clean_name) and "ligand" in f and f.endswith(".pdbqt"):
                extracted_model = os.path.join(temp_dir, f)
                break
        
        if extracted_model:
            destination = os.path.join(output_dir, temp_pdbqt_name)
            if os.path.exists(destination):
                log(f"  Overwriting existing file: {destination}", category="overwrite")
            shutil.move(extracted_model, destination)
            log(f"  SUCCESS: Moved extracted model to {destination}", category="success")
            
            # Clean temp folder
            for f in os.listdir(temp_dir):
                os.remove(os.path.join(temp_dir, f))
        else:
            log(f"  FAILURE: vina_split ran, but no output file matching '{clean_name}...ligand...' was found.", category="vina_fail")

    # 5. Cleanup
    try:
        shutil.rmtree(temp_dir)
        log("Temp directory removed.")
    except Exception as e:
        print(f"Warning: Could not fully remove temp dir: {e}")

    print("Done!")

if __name__ == "__main__":
    main()












"""
Bash example:
python extract_best.py 0.5 "./final_models"

0.5 is your cutoff score (the "so many units from ABA").
"./final_models" is where the files will go.
"""