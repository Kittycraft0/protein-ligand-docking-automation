import argparse
from pathlib import Path
import sys
import shutil
import re
import os
import time

def create_file_index(source_dir: Path):
    """
    Performs a one-time scan of the source directory to create a fast lookup
    index of all PDBQT files.
    
    Returns:
        A dictionary mapping filename (str) to its full Path object.
    """
    print("Pre-scanning results directory to build a file index (this may take a moment)...")
    start_time = time.time()
    
    # The key is the filename, the value is the full path object
    file_index = {p.name: p for p in source_dir.glob("**/*.pdbqt")}
    
    end_time = time.time()
    print(f"Indexing complete. Found {len(file_index)} files in {end_time - start_time:.2f} seconds.\n")
    
    # 11/6/2025
    print(f"First five files of file_index:")
    print(list(file_index.items())[:5])
    
    return file_index

def collect_docking_hits(results_file: Path, source_results_dir: Path, output_dir: Path, file_index: dict, top_n: int = None, score_range: list = None):
    """
    Reads a ranked results file and copies the corresponding docked PDBQT files
    using the pre-built file index for high performance.
    """
    # --- 1. Validation and Setup ---
    if not results_file.is_file():
        print(f"Error: Results file not found at '{results_file}'")
        sys.exit(1)

    # --- 2. Extract Protein Name ---
    protein_name_match = re.search(r'top_dockers_(.+?)(?:\.txt|$)', results_file.name)
    if not protein_name_match:
        # Fallback for RMS files or other formats
        protein_name_match = re.search(r'_(PYR1_3K3K_DOCK|PYL2_3KDI_DOCK|1stp_DOCK)', results_file.name)
    
    if not protein_name_match:
        protein_name = input(f"Could not automatically determine protein name from '{results_file.name}'. Please enter the target protein name (e.g., PYR1_3K3K_DOCK): ")
        if not protein_name:
             sys.exit("Operation cancelled.")
    else:
        protein_name = protein_name_match.group(1)
    
    print(f"Detected Protein Target: {protein_name}")

    # --- 3. Create Output Directory and Read Results ---
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: '{output_dir}'")

    with open(results_file, 'r') as f:
        lines = f.readlines()

    # --- 4. Filter Lines ---
    lines_to_process = []
    if score_range:
        min_score, max_score = score_range
        print(f"\nFiltering for ligands with scores between {min_score} and {max_score}...")
        for line in lines:
            if line.startswith("#") or not line.strip(): continue
            try:
                score = float(line.split()[0])
                if min_score <= score <= max_score:
                    lines_to_process.append(line)
            except (ValueError, IndexError): continue
    else: # Default to top_n
        print(f"\nSelecting the top {top_n} ligands from '{results_file.name}'...")
        valid_lines = [line for line in lines if not line.startswith("#") and line.strip()]
        lines_to_process = valid_lines[:top_n]

    print(f"Found {len(lines_to_process)} docked poses matching criteria.")
    if not lines_to_process: return

    # --- 5. Find and Copy Files using the Index ---
    copied_count = 0
    failed_to_find = []
    print("\nStarting file collection using the index...")

    for i, line in enumerate(lines_to_process):
        sys.stdout.write(f"\rProcessing ligand {i+1}/{len(lines_to_process)}...")
        sys.stdout.flush()

        parts = line.split()
        if len(parts) < 2: continue
        model_name = parts[1]
        
        # --- NEW, EFFICIENT LOOKUP LOGIC ---
        # Construct the expected filename and look for it in our pre-built index.
        expected_filename = f"{model_name}_vs_{protein_name}.pdbqt"
        source_file_path = file_index.get(expected_filename)

        # Fallback for names with extra text (e.g., _model0)
        if not source_file_path:
            for fname, fpath in file_index.items():
                if fname.startswith(model_name) and fname.endswith(f"_vs_{protein_name}.pdbqt"):
                    source_file_path = fpath
                    break
        # --- END OF NEW LOGIC ---

        if source_file_path:
            try:
                shutil.copy(source_file_path, output_dir)
                copied_count += 1
            except Exception as e:
                print(f"\nERROR: Could not copy {source_file_path}. Reason: {e}")
        else:
            print("Failed to find a molecule:")
            print(f"Name: {model_name}")
            print(f"File index: {file_index}")
            print(f"Expected filename: {expected_filename}")
            print(f"Source file path: {source_file_path}")
            print(f"Output directory: {output_dir}")
            print("Remember, you put in the WHOLE suffix, i.e. KCNH2_7CNOB_DOCK, not just what you want it to be. " \
            "The program is looking for files with the name you type in, so make sure it matches.")
            failed_to_find.append(model_name)
            break
            
    # --- 6. Final Report ---
    print(f"\n\nProcess complete. Successfully copied {copied_count} of {len(lines_to_process)} files to '{output_dir}'.")
    if failed_to_find:
        print(f"\nWARNING: Could not find result files for the following {len(failed_to_find)} models:")
        # 11/6/2025 added iterator
        iterator=0
        for name in failed_to_find:
            print(f"  - {name}")
            iterator+=1
            if iterator > 10:
                print(f"  - and {len(failed_to_find)-10} others")
                break

def main():
    parser = argparse.ArgumentParser(
        description="Collects the final docked PDBQT files for top-scoring ligands from a results file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'results_file', 
        type=str, 
        help="Path to the ranked results file (e.g., 'dock/results/scores/top_dockers_PYR1_3K3K_DOCK.txt')."
    )
    parser.add_argument(
        'output_dir', 
        type=str, 
        help="Path for the new directory where the top ligand files will be copied."
    )
    parser.add_argument(
        '--source_dir', 
        type=str, 
        default="dock/results/docked_ligands",
        help="The source directory of your final docked ligand results (default: 'dock/results/docked_ligands')."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-n', '--number', 
        type=int, 
        help="Collect the top N ligands from the results file."
    )
    group.add_argument(
        '--score_range',
        nargs=2,
        type=float,
        metavar=('MIN_SCORE', 'MAX_SCORE'),
        help="Collects all ligands with scores between MIN_SCORE and MAX_SCORE (e.g., --score_range -10.0 -8.5)."
    )
    
    args = parser.parse_args()
    
    # Create the file index once
    file_index = create_file_index(Path(args.source_dir))
    
    collect_docking_hits(
        results_file=Path(args.results_file),
        source_results_dir=Path(args.source_dir),
        output_dir=Path(args.output_dir),
        file_index=file_index,
        top_n=args.number,
        score_range=args.score_range
    )

if __name__ == "__main__":
    main()
