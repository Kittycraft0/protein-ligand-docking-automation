import os
import shutil
import subprocess
import math
import time
import signal
from pathlib import Path
from time import sleep
from collections import deque

# Directory paths
BASE_DIR = Path("./dock")
PROTEIN_DIR = BASE_DIR / "proteins"
LIGAND_DIR = BASE_DIR / "ligands"
CONFIG_DIR = BASE_DIR / "config"
CACHE_DIR = BASE_DIR / "cache"
RESULTS_DIR = BASE_DIR / "results"
COMPARISON_LIGAND_DIR = BASE_DIR / "comparison_ligands"

# Cache files
LIGAND_NAMES = CACHE_DIR / "ligandNames.txt"
PROTEIN_NAMES = CACHE_DIR / "proteinNames.txt"
PROGRESS_CACHE = CACHE_DIR / "progress_cache.txt"
COMPARISON_LIGAND_NAMES = CACHE_DIR / "comparisonLigandNames.txt"

# Ensure necessary directories exist
for directory in [PROTEIN_DIR, LIGAND_DIR, CONFIG_DIR, CACHE_DIR, RESULTS_DIR, COMPARISON_LIGAND_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Initialize DEBUG variable
DEBUG = "--debug" in os.sys.argv

# Global variable to track if display_progress is rendering
IS_RENDERING = False

# Global variable to track the start time of the docking process
START_TIME = time.time()

# Global variable to track task durations
TASK_DURATIONS = deque(maxlen=10)  # Store durations for the last 10 tasks

# Function to wait until IS_RENDERING is false
def wait_for_render_stop():
    """
    Wait until the IS_RENDERING flag is set to False.
    """
    global IS_RENDERING
    while IS_RENDERING:
        time.sleep(0.01)  # Wait briefly for rendering to complete

# Function to initialize or clear cache
def initialize_cache(mode=None):
    if DEBUG:
        print("Debug mode enabled. Skipping cache and results clearing.")
        return

    if mode == "clear":
        backup_dir = CACHE_DIR / "cache_backup"
        backup_dir.mkdir(parents=True, exist_ok=True)
        for item in CACHE_DIR.iterdir():
            # Skip the cache_backup directory itself
            if item == backup_dir:
                continue

            destination = backup_dir / item.name
            if destination.exists():
                # Rename the file to avoid overwriting
                counter = 1
                new_name = f"{item.stem}_copy{counter}{item.suffix}"
                new_destination = backup_dir / new_name
                while new_destination.exists():
                    counter += 1
                    new_name = f"{item.stem}_copy{counter}{item.suffix}"
                    new_destination = backup_dir / new_name
                destination = new_destination
            shutil.move(str(item), str(destination))
        print("Cache cleared and backed up.")
    elif mode == "clear-everything":
        shutil.rmtree(CACHE_DIR, ignore_errors=True)
        shutil.rmtree(RESULTS_DIR, ignore_errors=True)
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        print("Cache and results cleared.")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Generate ligand and protein name lists
    if not list(LIGAND_DIR.glob("*.pdbqt")):
        print(f"Error: No ligand files found in {LIGAND_DIR}")
        exit(1)
    with open(LIGAND_NAMES, "w") as f:
        f.writelines([str(file) + "\n" for file in LIGAND_DIR.glob("*.pdbqt")])

    if not list(PROTEIN_DIR.glob("*.pdbqt")):
        print(f"Error: No protein files found in {PROTEIN_DIR}")
        exit(1)
    with open(PROTEIN_NAMES, "w") as f:
        f.writelines([str(file) + "\n" for file in PROTEIN_DIR.glob("*.pdbqt")])

    # Generate comparison ligand name list
    if not list(COMPARISON_LIGAND_DIR.glob("*.pdbqt")):
        print(f"Error: No comparison ligand files found in {COMPARISON_LIGAND_DIR}")
        exit(1)
    with open(COMPARISON_LIGAND_NAMES, "w") as f:
        f.writelines([str(file) + "\n" for file in COMPARISON_LIGAND_DIR.glob("*.pdbqt")])

    # Initialize progress cache if it doesn't exist
    if not PROGRESS_CACHE.exists():
        with open(PROGRESS_CACHE, "w") as f:
            f.write("COMPARISON_LIGAND_INDEX=0\nCOMPARISON_PROTEIN_INDEX=0\nLIGAND_INDEX=0\nMODEL_INDEX=0\nPROTEIN_INDEX=0\n")

# Function to read progress from cache
def read_progress():
    if PROGRESS_CACHE.exists():
        progress = {}
        with open(PROGRESS_CACHE) as f:
            for line in f:
                key, value = line.strip().split("=")
                progress[key] = int(value)
        return progress
    else:
        print("Error: Progress cache file not found.")
        exit(1)

# Function to update progress in cache
def update_progress(progress):
    with open(PROGRESS_CACHE, "w") as f:
        for key, value in progress.items():
            f.write(f"{key}={value}\n")

# Function to calculate docking box parameters
def calculate_docking_box(protein_file):
    center_x, center_y, center_z = 0, 0, 0
    count = 0
    with open(protein_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                center_x += float(parts[6])
                center_y += float(parts[7])
                center_z += float(parts[8])
                count += 1
    if count > 0:
        center_x /= count
        center_y /= count
        center_z /= count
    size_x, size_y, size_z = 20, 20, 20
    return center_x, center_y, center_z, size_x, size_y, size_z

# Function to extract models from ligand files
def extract_models(ligand_file):
    ligand_name = Path(ligand_file).stem  # Extract the base name without extension
    output_dir = CACHE_DIR / f"models_{ligand_name}"  # Output directory for models
    output_dir.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist

    # Check if the models are already extracted
    if list(output_dir.glob(f"{ligand_name}_model*.pdbqt")):
        print(f"Models for {ligand_name} already extracted.")
        return

    # Check if the ligand file contains multiple models
    with open(ligand_file) as f:
        if "MODEL" in f.read():
            # Use vina_split to split the ligand file into multiple PDBQT files
            result = subprocess.run(
                ["vina_split", "--input", ligand_file, "--ligand", str(output_dir / f"{ligand_name}_model")],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                print(f"Error: Failed to split ligand file {ligand_file}")
                exit(1)

            # Check if the output directory contains any split files
            if not list(output_dir.glob(f"{ligand_name}_model*.pdbqt")):
                print(f"Error: No models generated for ligand file {ligand_file}")
                exit(1)

            # Rename the split files to follow the expected naming convention
            for model_file in output_dir.glob(f"{ligand_name}_model*.pdbqt"):
                model_index = model_file.stem.split("_")[-1]  # Extract the model index
                model_file.rename(output_dir / f"{ligand_name}_model_{model_index}.pdbqt")
        else:
            # Copy the single model file to the output directory
            shutil.copy(ligand_file, output_dir / f"{ligand_name}_model_1.pdbqt")

# Function to draw a progress bar
def draw_progress_bar(row, col, length, current, total, waits=False):
    """
    Draw a progress bar in the terminal.
    :param row: Cursor row position for the progress bar.
    :param col: Cursor column position for the progress bar.
    :param length: Length of the progress bar.
    :param current: Current progress value.
    :param total: Total progress value.
    :param waits: Optional flag to wait after each character (default: False).
    """
    # Check if the current value exceeds the total value
    if current > total:
        current = total

    # Check if the current value is negative
    if current < 0:
        current = 0

    # Check if the progress bar length is valid
    if length <= 0:
        print("Error: Progress bar length must be greater than 0.")
        return

    # Check if the current and total values are valid
    if current < 0 or total <= 0:
        print("Error: Progress bar current and total values must be non-negative, and total must be greater than 0.")
        return

    # Calculate progress as a fraction (0 to 1)
    progress = current / total

    # Calculate the number of filled and empty segments
    filled = int(progress * length)
    empty = length - filled

    # Draw the progress bar
    bar = "[" + "\033[42m \033[0m" * filled + "\033[41m \033[0m" * empty + "]"
    percent = f"{progress * 100:.2f}%"

    # Move the cursor to the specified position and print the progress bar
    print(f"\033[{row};{col}H{bar} {percent}", end="", flush=True)

# Function to calculate RMS relative to ABA number
def calculate_rms_relative_to_aba(aba_score, score):
    """
    Calculate the RMS relative to the ABA score.
    :param aba_score: The ABA score (float).
    :param score: The docking score (float).
    :return: RMS relative to ABA (float).
    """
    rms_relative_to_aba = math.sqrt((score - aba_score) ** 2)
    return rms_relative_to_aba

# Calculate RMS relative to comparison ligands
def calculate_rms_relative_to_comparison(score, comparison_scores):
    """
    Calculate the RMS relative to comparison ligands.
    :param score: The docking score for the ligand.
    :param comparison_scores: List of comparison ligand scores.
    :return: RMS relative to comparison ligands.
    """
    if not comparison_scores:
        return float("inf")  # Return infinity if no comparison scores are available
    squared_differences = [(float(score) - float(c_score)) ** 2 for c_score in comparison_scores]
    rms = math.sqrt(sum(squared_differences) / len(squared_differences))
    return rms

# Function to monitor the progress of the current docking process with a text-based progress bar
def monitor_docking_progress(log_file):
    """
    Monitor the docking progress by reading the log file and updating a text-based progress bar.
    :param log_file: Path to the log file being monitored.
    """
    console_width = shutil.get_terminal_size((80, 20)).columns  # Get terminal width
    progress_bar_length = console_width - 30  # Adjust progress bar length
    progress = 0

    # Wait for the log file to be created
    while not Path(log_file).exists():
        time.sleep(0.2)  # Check every 0.2 seconds

    dock_progress_row=15 # Row for the progress bar
    
    # Monitor the log file for progress updates
    while progress < 100:
        if Path(log_file).exists():
            # Count the number of '*' characters in the log file to determine progress
            with open(log_file, "r") as f:
                progress = f.read().count("*") * 2  # Each '*' represents 2% progress

            # Ensure progress does not exceed 100%
            progress = min(progress, 100)

            # Update the progress bar
            move_cursor(dock_progress_row, 0)  # Move cursor to the progress bar row
            print(f"\rCurrent docking progress: {progress}%   ", end="")
            draw_progress_bar(dock_progress_row+1, 0, progress_bar_length, progress, 100)  # Fixed call

        time.sleep(0.2)  # Update display every 0.2 seconds

    # Ensure the progress bar is fully green at the end
    move_cursor(dock_progress_row, 0)  # Move cursor to the progress bar row
    print("\nCurrent docking progress: 100%")
    draw_progress_bar(dock_progress_row+1, 0, progress_bar_length, 100, 100)  # Fixed call

# Function to perform docking asynchronously
def perform_docking(ligand_file, protein_file, model_index):
    """
    Perform docking asynchronously using AutoDock Vina.
    :param ligand_file: Path to the ligand file.
    :param protein_file: Path to the protein file.
    :param model_index: Index of the ligand model.
    """
    ligand_name = Path(ligand_file).stem  # Extract the base name without extension
    protein_name = Path(protein_file).stem  # Extract the base name without extension
    output_file = RESULTS_DIR / f"temp/{ligand_name}_model{model_index}_vs_{protein_name}.pdbqt"
    log_file = RESULTS_DIR / f"temp/{ligand_name}_model{model_index}_vs_{protein_name}.log"

    # Ensure the temp directory exists
    (RESULTS_DIR / "temp").mkdir(parents=True, exist_ok=True)

    # Calculate docking box parameters
    center_x, center_y, center_z, size_x, size_y, size_z = calculate_docking_box(protein_file)

    # Debugging information
    if DEBUG:
        print("Running AutoDock Vina with the following parameters:")
        print(f"Receptor: {protein_file}")
        print(f"Ligand: {ligand_file}")
        print(f"Center: {center_x}, {center_y}, {center_z}")
        print(f"Size: {size_x}, {size_y}, {size_z}")
        print(f"Output: {output_file}")
        print(f"Log: {log_file}")

    # Run AutoDock Vina asynchronously
    vina_command = [
        "vina",
        "--receptor", str(protein_file),
        "--ligand", str(ligand_file),
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--out", str(output_file)
    ]

    with open(log_file, "w") as log:
        vina_process = subprocess.Popen(vina_command, stdout=log, stderr=log)

    # Monitor the docking progress
    monitor_docking_progress(log_file)

    # Wait for the docking process to complete
    vina_process.wait()

    # Check if the log file was created
    if not log_file.exists():
        print(f"Error: Log file not found: {log_file}")
        exit(1)

    # Display the contents of the log file for debugging
    if DEBUG:
        print(f"Contents of log file {log_file}:")
        with open(log_file, "r") as log:
            print(log.read())

    # Extract score from log file and store it
    score = None
    with open(log_file, "r") as log:
        for line in log:
            if line.startswith("   1"):  # Look for the first result line
                score = line.split()[1]  # Extract the score
                break

    if score and score.replace(".", "", 1).isdigit():  # Check if the score is a valid number
        if DEBUG:
            print(f"Extracted score: {score} for ligand: {ligand_name} with protein: {protein_name}")
        scores_file = RESULTS_DIR / f"scores_{protein_name}.txt"
        with open(scores_file, "a") as f:
            f.write(f"{score} {ligand_name}\n")
        if DEBUG:
            print(f"Score stored in {scores_file}")
    else:
        print(f"Failed to extract a valid score from log file: {log_file}")

# Function to display progress
def display_progress(current_task, total_tasks, ligand_index, total_ligands, model_index, total_models, protein_index, total_proteins, ligand_name, protein_name):
    """
    Display the progress of the docking process with an estimated time to completion (ETC).
:param current_task: Current task number.
    :param total_tasks: Total number of tasks.
    :param ligand_index: Current ligand index.
    :param total_ligands: Total number of ligands.
    :param model_index: Current model index.
    :param total_models: Total number of models.
    :param protein_index: Current protein index.
    :param total_proteins: Total number of proteins.
    :param ligand_name: Name of the current ligand.
    :param protein_name: Name of the current protein.
    """
    console_width = shutil.get_terminal_size((80, 20)).columns
    progress_bar_length = console_width - 30

    # Calculate elapsed time for the last 10 tasks
    if TASK_DURATIONS:
        avg_task_time = sum(TASK_DURATIONS) / len(TASK_DURATIONS)
    else:
        avg_task_time = 0

    # Calculate remaining tasks and ETC
    remaining_tasks = total_tasks - current_task
    remaining_time = avg_task_time * remaining_tasks if avg_task_time > 0 else 0

    # Format remaining time as HH:MM:SS
    hours, rem = divmod(int(remaining_time), 3600)
    minutes, seconds = divmod(rem, 60)
    etc = f"{hours:02}:{minutes:02}:{seconds:02}"

    print("\033[H\033[J", end="")  # Clear the terminal
    print("========================================")
    print("         Protein-Ligand Docking         ")
    print("========================================")
    print()
    print(f"Docking:")
    print(f"Ligand file {ligand_index + 1}/{total_ligands}")
    draw_progress_bar(7, 0, progress_bar_length, ligand_index + 1, total_ligands, False)

    move_cursor(9, 0)
    print(f"Ligand model: model {model_index + 1}/{total_models}")
    draw_progress_bar(10, 0, progress_bar_length, model_index + 1, total_models, False)

    move_cursor(12, 0)
    print(f"Ligand progress: protein {protein_index + 1}/{total_proteins}")
    draw_progress_bar(13, 0, progress_bar_length, protein_index + 1, total_proteins, False)

    # Current docking progress bar
    move_cursor(15, 0)
    print("Current docking progress:")
    draw_progress_bar(16, 0, progress_bar_length, 0, 100, False)

    move_cursor(19, 0)
    print(f"Ligand name: {ligand_name}")
    print(f"Protein name: {protein_name}")

    percent_through_dock = (current_task * 100) / total_tasks
    move_cursor(21, 0)
    print(f"Total progress: {current_task}/{total_tasks}")
    draw_progress_bar(22, 0, progress_bar_length, current_task, total_tasks, False)

    # Display estimated time to completion
    move_cursor(24, 0)
    print(f"Estimated time to completion: {etc}")

    move_cursor(26, 0)
    print("Do CTRL+C to exit")

# Function to display comparison progress
def display_comparison_progress(current_task, total_tasks, comparison_ligand_index, total_comparison_ligands, protein_index, total_proteins, comparison_ligand_name, protein_name):
    """
    Display the progress of docking for comparison ligands.
    :param current_task: Current task number.
    :param total_tasks: Total number of tasks.
    :param comparison_ligand_index: Current comparison ligand index.
    :param total_comparison_ligands: Total number of comparison ligands.
    :param protein_index: Current protein index.
    :param total_proteins: Total number of proteins.
    :param comparison_ligand_name: Name of the current comparison ligand.
    :param protein_name: Name of the current protein.
    """
    console_width = shutil.get_terminal_size((80, 20)).columns
    progress_bar_length = console_width - 30

    print("\033[H\033[J", end="")  # Clear the terminal
    print("========================================")
    print("   Initial Calculation of Comparison    ")
    print("========================================")
    print()
    print(f"Docking comparisons:")
    print(f"Total progress: comparison ligand {comparison_ligand_index + 1}/{total_comparison_ligands}")
    draw_progress_bar(7, 0, progress_bar_length, comparison_ligand_index + 1, total_comparison_ligands, False)
    move_cursor(9, 0)
    print(f"Comparison ligand progress: protein {protein_index + 1}/{total_proteins}")
    draw_progress_bar(10, 0, progress_bar_length, protein_index + 1, total_proteins, False)
    move_cursor(19, 0)
    print(f"Comparison ligand name: {comparison_ligand_name}")
    print(f"Protein name: {protein_name}")

    percent_through_dock = (current_task * 100) / total_tasks
    move_cursor(21, 0)
    print(f"% through dock: {percent_through_dock:.2f}%")
    draw_progress_bar(22, 0, progress_bar_length, current_task, total_tasks, False)

    move_cursor(24, 0)
    print("Do CTRL+C to exit")
    print()

# Function to add empty lines
def add_empty_lines():
    """
    Add empty lines to create space for display updates.
    """
    console_height = shutil.get_terminal_size((80, 20)).lines
    for _ in range(console_height):
        print()

# Function to clear the display area
def clear_display_area():
    """
    Clear the display area in the terminal.
    """
    console_height = shutil.get_terminal_size((80, 20)).lines
    console_width = shutil.get_terminal_size((80, 20)).columns
    print("\033[H", end="")  # Move cursor to the top left corner
    for _ in range(console_height):
        print(" " * console_width)
    print("\033[H", end="")  # Move cursor back to the top left corner

# Function to move the cursor to a specific position
def move_cursor(row, col):
    """
    Move the cursor to a specific row and column in the terminal.
    :param row: The row number (1-based).
    :param col: The column number (1-based).
    """
    print(f"\033[{row};{col}H", end="", flush=True)

# Function to handle termination signals
def terminate_script(signal_received, frame):
    """
    Handle termination signals and save progress before exiting.
    """
    print("\nTermination signal received. Saving progress...")
    progress = {
        "COMPARISON_LIGAND_INDEX": COMPARISON_LIGAND_INDEX,
        "COMPARISON_PROTEIN_INDEX": COMPARISON_PROTEIN_INDEX,
        "LIGAND_INDEX": LIGAND_INDEX,
        "MODEL_INDEX": MODEL_INDEX,
        "PROTEIN_INDEX": PROTEIN_INDEX
    }
    update_progress(progress)
    exit(0)

# Move all .log and .pdbqt files to the dedicated subfolder
def move_temp_files():
    """
    Move all .log and .pdbqt files from the temp directory to the results directory.
    If a file with the same name already exists in the destination, rename the file being moved.
    """
    temp_dir = RESULTS_DIR / "temp"
    for file in temp_dir.glob("*"):
        destination = RESULTS_DIR / file.name
        if destination.exists():
            # Rename the file to avoid overwriting
            counter = 1
            new_name = f"{file.stem}_copy{counter}{file.suffix}"
            new_destination = RESULTS_DIR / new_name
            while new_destination.exists():
                counter += 1
                new_name = f"{file.stem}_copy{counter}{file.suffix}"
                new_destination = RESULTS_DIR / new_name
            destination = new_destination
        shutil.move(str(file), str(destination))
    print("Moved all temporary files to the results directory.")

# Calculate and display top dockers for each protein
def calculate_top_dockers(proteins):
    """
    Sort and save the top dockers for each protein.
    :param proteins: List of protein file paths.
    """
    for protein_file in proteins:
        protein_name = Path(protein_file).stem
        scores_file = RESULTS_DIR / f"scores_{protein_name}.txt"
        top_dockers_file = RESULTS_DIR / f"top_dockers_{protein_name}.txt"

        if scores_file.exists():
            with open(scores_file, "r") as f:
                sorted_scores = sorted(f.readlines(), key=lambda x: float(x.split()[0]))
            with open(top_dockers_file, "w") as f:
                f.writelines(sorted_scores)
            print(f"Top dockers for {protein_name} saved to {top_dockers_file}.")
        else:
            print(f"No scores file found for protein: {protein_name}. Skipping sorting.")

# Calculate the best ligands using RMS relative to comparison ligands
def calculate_best_ligands(ligands, proteins, comparison_ligands):
    """
    Calculate the best ligands using RMS relative to comparison ligands.
    :param ligands: List of ligand file paths.
    :param proteins: List of protein file paths.
    :param comparison_ligands: List of comparison ligand file paths.
    """
    best_ligands_file = RESULTS_DIR / "best_ligands.txt"
    with open(best_ligands_file, "w") as best_ligands:
        for ligand_file in ligands:
            ligand_name = Path(ligand_file).stem
            total_score = 0
            count = 0

            for protein_file in proteins:
                protein_name = Path(protein_file).stem
                scores_file = RESULTS_DIR / f"scores_{protein_name}.txt"

                if not scores_file.exists():
                    continue

                # Get the score for the ligand
                score = None
                with open(scores_file, "r") as f:
                    for line in f:
                        if ligand_name in line:
                            score = float(line.split()[0])
                            break

                if score is not None:
                    # Get comparison scores
                    comparison_scores = []
                    for comparison_ligand_file in comparison_ligands:
                        comparison_ligand_name = Path(comparison_ligand_file).stem
                        with open(scores_file, "r") as f:
                            for line in f:
                                if comparison_ligand_name in line:
                                    comparison_scores.append(float(line.split()[0]))
                                    break

                    # Calculate RMS relative to comparison ligands
                    rms_relative_to_comparison = calculate_rms_relative_to_comparison(score, comparison_scores)
                    if rms_relative_to_comparison > 0:
                        total_score += 1 / (rms_relative_to_comparison ** 2)
                        count += 1

            if count > 0:
                final_score = math.sqrt(count / total_score)
                best_ligands.write(f"{final_score:.8f} {ligand_name}\n")
                print(f"Final score for {ligand_name}: {final_score:.8f}")

# Rank and display the best ligands
def rank_and_display_best_ligands():
    """
    Rank the best ligands and display the top 5.
    """
    best_ligands_file = RESULTS_DIR / "best_ligands.txt"
    ranked_best_ligands_file = RESULTS_DIR / "ranked_best_ligands.txt"

    if best_ligands_file.exists():
        with open(best_ligands_file, "r") as f:
            sorted_ligands = sorted(f.readlines(), key=lambda x: float(x.split()[0]))
        with open(ranked_best_ligands_file, "w") as f:
            f.writelines(sorted_ligands)
        print(f"Best ligands ranked and saved to {ranked_best_ligands_file}.")

        # Display the top 5 ligands
        print("Top 5 ligands:")
        for line in sorted_ligands[:5]:
            print(line.strip())
    else:
        print("No best ligands file found. Skipping ranking.")

# Register the termination signal handler
signal.signal(signal.SIGINT, terminate_script)
signal.signal(signal.SIGTERM, terminate_script)

# Main script execution
if __name__ == "__main__":
    import sys

    # Handle command-line arguments for clearing cache
    if "--clear-cache" in sys.argv:
        initialize_cache("clear")
    elif "--clear-everything" in sys.argv:
        initialize_cache("clear-everything")
    else:
        initialize_cache()

    # Read progress from cache
    progress = read_progress()
    COMPARISON_LIGAND_INDEX = progress.get("COMPARISON_LIGAND_INDEX", 0)
    COMPARISON_PROTEIN_INDEX = progress.get("COMPARISON_PROTEIN_INDEX", 0)
    LIGAND_INDEX = progress.get("LIGAND_INDEX", 0)
    MODEL_INDEX = progress.get("MODEL_INDEX", 0)
    PROTEIN_INDEX = progress.get("PROTEIN_INDEX", 0)

    # Load ligands, proteins, and comparison ligands
    ligands = [line.strip() for line in open(LIGAND_NAMES)]
    proteins = [line.strip() for line in open(PROTEIN_NAMES)]
    comparison_ligands = [line.strip() for line in open(COMPARISON_LIGAND_NAMES)]

    total_ligands = len(ligands)
    total_proteins = len(proteins)
    total_comparison_ligands = len(comparison_ligands)

    # Create a dedicated subfolder for output files
    (RESULTS_DIR / "temp").mkdir(parents=True, exist_ok=True)

    # Perform docking for comparison ligands
    while COMPARISON_LIGAND_INDEX < total_comparison_ligands:
        comparison_ligand_file = comparison_ligands[COMPARISON_LIGAND_INDEX]
        comparison_ligand_name = Path(comparison_ligand_file).stem

        while COMPARISON_PROTEIN_INDEX < total_proteins:
            protein_file = proteins[COMPARISON_PROTEIN_INDEX]
            protein_name = Path(protein_file).stem

            # Check if docking for this comparison ligand and protein has already been completed
            scores_file = RESULTS_DIR / f"scores_{protein_name}.txt"
            if scores_file.exists():
                with open(scores_file, "r") as f:
                    if any(comparison_ligand_name in line for line in f):
                        print(f"Skipping docking for {comparison_ligand_name} with {protein_name} (already completed).")
                        COMPARISON_PROTEIN_INDEX += 1
                        update_progress({
                            "COMPARISON_LIGAND_INDEX": COMPARISON_LIGAND_INDEX,
                            "COMPARISON_PROTEIN_INDEX": COMPARISON_PROTEIN_INDEX,
                            "LIGAND_INDEX": LIGAND_INDEX,
                            "MODEL_INDEX": MODEL_INDEX,
                            "PROTEIN_INDEX": PROTEIN_INDEX
                        })
                        continue

            # Display progress for comparison ligands
            display_comparison_progress(
                current_task=COMPARISON_LIGAND_INDEX * total_proteins + COMPARISON_PROTEIN_INDEX + 1,
                total_tasks=total_comparison_ligands * total_proteins,
                comparison_ligand_index=COMPARISON_LIGAND_INDEX,
                total_comparison_ligands=total_comparison_ligands,
                protein_index=COMPARISON_PROTEIN_INDEX,
                total_proteins=total_proteins,
                comparison_ligand_name=comparison_ligand_name,
                protein_name=protein_name
            )

            # Perform docking
            perform_docking(comparison_ligand_file, protein_file, 0)

            # Update progress after each docking task
            COMPARISON_PROTEIN_INDEX += 1
            update_progress({
                "COMPARISON_LIGAND_INDEX": COMPARISON_LIGAND_INDEX,
                "COMPARISON_PROTEIN_INDEX": COMPARISON_PROTEIN_INDEX,
                "LIGAND_INDEX": LIGAND_INDEX,
                "MODEL_INDEX": MODEL_INDEX,
                "PROTEIN_INDEX": PROTEIN_INDEX
            })

        # Reset protein index and move to the next comparison ligand
        COMPARISON_PROTEIN_INDEX = 0
        COMPARISON_LIGAND_INDEX += 1
        update_progress({
            "COMPARISON_LIGAND_INDEX": COMPARISON_LIGAND_INDEX,
            "COMPARISON_PROTEIN_INDEX": COMPARISON_PROTEIN_INDEX,
            "LIGAND_INDEX": LIGAND_INDEX,
            "MODEL_INDEX": MODEL_INDEX,
            "PROTEIN_INDEX": PROTEIN_INDEX
        })

    # Clear the display area after comparison ligand docking
    clear_display_area()

    # start time for docking
    start_time = time.time()

    # Perform docking for ligands
    while LIGAND_INDEX < total_ligands:
        ligand_file = ligands[LIGAND_INDEX]
        extract_models(ligand_file)
        models = list((CACHE_DIR / f"models_{Path(ligand_file).stem}").glob("*.pdbqt"))
        total_models = len(models)

        # Calculate total docking tasks
        total_tasks = total_ligands * total_models * total_proteins
        current_task = (
            LIGAND_INDEX * total_models * total_proteins +
            MODEL_INDEX * total_proteins +
            PROTEIN_INDEX + 1
        )

        while MODEL_INDEX < total_models:
            model_file = models[MODEL_INDEX]

            while PROTEIN_INDEX < total_proteins:
                protein_file = proteins[PROTEIN_INDEX]

                # Update ligand and protein names for progress display
                ligand_name = Path(model_file).stem
                protein_name = Path(protein_file).stem

                # Display progress before starting docking
                display_progress(
                    current_task=current_task,
                    total_tasks=total_tasks,
                    ligand_index=LIGAND_INDEX,
                    total_ligands=total_ligands,
                    model_index=MODEL_INDEX,
                    total_models=total_models,
                    protein_index=PROTEIN_INDEX,
                    total_proteins=total_proteins,
                    ligand_name=ligand_name,
                    protein_name=protein_name
                )

                # Perform docking
                perform_docking(model_file, protein_file, MODEL_INDEX)


                end_time = time.time()
                # Record task duration
                TASK_DURATIONS.append(end_time - start_time)
                start_time = end_time
                # Update progress
                current_task += 1
                PROTEIN_INDEX += 1
                update_progress({
                    "COMPARISON_LIGAND_INDEX": COMPARISON_LIGAND_INDEX,
                    "COMPARISON_PROTEIN_INDEX": COMPARISON_PROTEIN_INDEX,
                    "LIGAND_INDEX": LIGAND_INDEX,
                    "MODEL_INDEX": MODEL_INDEX,
                    "PROTEIN_INDEX": PROTEIN_INDEX
                })

            PROTEIN_INDEX = 0
            MODEL_INDEX += 1
            update_progress({
                "COMPARISON_LIGAND_INDEX": COMPARISON_LIGAND_INDEX,
                "COMPARISON_PROTEIN_INDEX": COMPARISON_PROTEIN_INDEX,
                "LIGAND_INDEX": LIGAND_INDEX,
                "MODEL_INDEX": MODEL_INDEX,
                "PROTEIN_INDEX": PROTEIN_INDEX
            })

        MODEL_INDEX = 0
        LIGAND_INDEX += 1
        update_progress({
            "COMPARISON_LIGAND_INDEX": COMPARISON_LIGAND_INDEX,
            "COMPARISON_PROTEIN_INDEX": COMPARISON_PROTEIN_INDEX,
            "LIGAND_INDEX": LIGAND_INDEX,
            "MODEL_INDEX": MODEL_INDEX,
            "PROTEIN_INDEX": PROTEIN_INDEX
        })

    # Move temporary files to the results directory
    move_temp_files()

    # Calculate and save top dockers for each protein
    calculate_top_dockers(proteins)

    # Calculate and save the best ligands
    calculate_best_ligands(ligands, proteins, comparison_ligands)

    # Rank and display the best ligands
    rank_and_display_best_ligands()