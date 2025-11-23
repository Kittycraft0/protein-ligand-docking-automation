import argparse
from pathlib import Path
import sys
import shutil
from functools import reduce
import time

# --- Helper Functions (reused from other scripts for file copying) ---

def create_file_index(source_dir: Path):
    """
    Performs a one-time scan of the source directory to create a fast lookup
    index of all PDBQT files.

    Returns:
        A dictionary mapping filename (str) to its full Path object.
    """
    if not source_dir.is_dir():
        print(f"Error: Source directory for docked files not found at '{source_dir}'")
        sys.exit(1)
        
    print(f"Scanning source directory '{source_dir}' to build a file index (this may take a moment)...")
    start_time = time.time()
    
    # The key is the filename, the value is the full path object
    file_index = {p.name: p for p in source_dir.glob("**/*.pdbqt")}
    
    end_time = time.time()
    print(f"Indexing complete. Found {len(file_index)} files in {end_time - start_time:.2f} seconds.\n")
    return file_index

def copy_intersecting_files(ligands_to_copy: set, source_dir: Path, output_dir: Path, file_index: dict):
    """
    Copies the PDBQT files for the intersecting ligands from the source directory
    to the output directory.

    Args:
        ligands_to_copy: A set of ligand names to be copied.
        source_dir: The directory where all docked PDBQT files are stored.
        output_dir: The destination directory for the copied files.
        file_index: A pre-built index of files in the source_dir.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    copied_count = 0
    total_to_find = len(ligands_to_copy)
    print(f"Searching for PDBQT files for {total_to_find} intersecting ligands...")

    # Since one ligand name can correspond to multiple files (docked against different proteins),
    # we iterate through the file index to find all matches.
    for i, ligand_name in enumerate(sorted(list(ligands_to_copy))):
        sys.stdout.write(f"\rProcessing ligand {i+1}/{total_to_find}: {ligand_name}")
        sys.stdout.flush()
        
        found_match_for_ligand = False
        for indexed_filename, full_path in file_index.items():
            # Check if the file in the index belongs to the current intersecting ligand
            if indexed_filename.startswith(ligand_name):
                try:
                    shutil.copy(full_path, output_dir)
                    copied_count += 1
                    found_match_for_ligand = True
                except Exception as e:
                    print(f"\nERROR: Could not copy {full_path}. Reason: {e}")
    
    print(f"\n\nProcess complete. Copied {copied_count} PDBQT files to '{output_dir}'.")


# --- Main Script Logic ---

class AddFilterAction(argparse.Action):
    """Custom argparse action to handle per-file filtering arguments."""
    def __call__(self, parser, namespace, values, option_string=None):
        # --- THIS IS THE CORRECTED PART ---
        # Check if the destination list exists, and if not, create it.
        # This prevents the 'NoneType' object has no attribute 'append' error.
        if getattr(namespace, self.dest, None) is None:
            setattr(namespace, self.dest, [])
        
        # Get the list we will append to.
        filters_list = getattr(namespace, self.dest)
        # --- END OF CORRECTION ---

        if len(values) < 2: # A file and a filter type are the minimum
            raise argparse.ArgumentError(self, "Each --input must be followed by a file and a filter specification.")
        
        file_path, filter_type = values[0], values[1]
        filter_spec = {'path': Path(file_path)}
        
        if filter_type == 'n':
            if len(values) != 3: parser.error(f"Filter 'n' for '{file_path}' requires exactly one number.")
            try:
                filter_spec.update({'type': 'number', 'value': int(values[2])})
            except ValueError: parser.error(f"Invalid integer '{values[2]}' for 'n' filter.")
        
        elif filter_type == 's':
            if len(values) != 4: parser.error(f"Filter 's' for '{file_path}' requires two numbers (min max).")
            try:
                min_s, max_s = float(values[2]), float(values[3])
                if min_s > max_s: parser.error(f"min_score ({min_s}) cannot be > max_score ({max_s}).")
                filter_spec.update({'type': 'score_range', 'value': (min_s, max_s)})
            except ValueError: parser.error(f"Invalid float for 's' filter.")
        
        else:
            parser.error(f"Unknown filter type '{filter_type}'. Use 'n' for top number or 's' for score range.")
            
        filters_list.append(filter_spec)

def find_intersecting_ligands(filter_specs: list[dict]):
    """Finds the set of ligands common to all filtered sets."""
    all_ligand_sets = []
    print(f"Finding intersection from {len(filter_specs)} files with individual filters...\n")
    for i, spec in enumerate(filter_specs):
        file_path, filter_type, filter_value = spec['path'], spec['type'], spec['value']
        if not file_path.is_file():
            print(f"Warning: File not found at '{file_path}'. Skipping.")
            continue
        print(f"Processing file {i+1}/{len(filter_specs)}: {file_path.name}")
        with open(file_path, 'r') as f:
            lines = f.readlines()
        current_file_ligands = set()
        if filter_type == 'number':
            print(f"  -> Applying filter: Top {filter_value} ligands.")
            names = [l.split()[1] for l in lines if not l.startswith("#") and l.strip()]
            current_file_ligands = set(names[:filter_value])
        elif filter_type == 'score_range':
            min_score, max_score = filter_value
            print(f"  -> Applying filter: Score range [{min_score}, {max_score}].")
            for line in lines:
                if line.startswith("#") or not line.strip(): continue
                try:
                    parts = line.split(); score, name = float(parts[0]), parts[1]
                    if min_score <= score <= max_score: current_file_ligands.add(name)
                except (ValueError, IndexError): continue
        print(f"  -> Found {len(current_file_ligands)} ligands matching criteria.")
        all_ligand_sets.append(current_file_ligands)
    if not all_ligand_sets: return set()
    return reduce(lambda set1, set2: set1.intersection(set2), all_ligand_sets)

def main():
    parser = argparse.ArgumentParser(
        description="Finds intersecting ligands from multiple score files and either lists them or copies their PDBQT files.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example Usage:
  # Output a TEXT FILE with names of intersecting ligands
  python %(prog)s -o "dock/results/intersection.txt" --output-mode file --input file1.txt n 10000 --input file2.txt n 5000

  # Copy intersecting ligand PDBQT FILES to a new directory
  python %(prog)s -o "dock/results/intersecting_hits" --output-mode directory --source-dir "dock/results/docked_ligands" --input file1.txt s -12.0 -10.5 --input file2.txt s -11.0 -9.5
"""
    )
    parser.add_argument(
        '-o', '--output',
        dest='output_path',
        type=str,
        required=True,
        help="Path for the output. A file path for 'file' mode, a directory path for 'directory' mode."
    )
    parser.add_argument(
        '--output-mode',
        choices=['file', 'directory'],
        default='file',
        help="Choose output type: 'file' (default) creates a text list of ligand names. 'directory' copies PDBQT files."
    )
    parser.add_argument(
        '--source-dir',
        type=str,
        default="dock/results/docked_ligands",
        help="[Required for 'directory' mode] The source directory of all docked PDBQT files."
    )
    parser.add_argument(
        '--input',
        dest='filters',
        action=AddFilterAction,
        nargs='+',
        metavar='FILE FILTER VAL(S)',
        required=True,
        help="Define an input file and its filter. Can be used multiple times.\n"
             "  Filter types:\n"
             "    n N:         Select top N ligands.\n"
             "    s MIN MAX:   Select ligands with scores between MIN and MAX."
    )
    
    args = parser.parse_args()
    
    # --- Main Logic ---
    common_ligands = find_intersecting_ligands(args.filters)

    if not common_ligands:
        print("\nNo common ligands found matching the specified criteria for all files.")
        return

    print(f"\nFound {len(common_ligands)} intersecting ligands.")
    output_path = Path(args.output_path)

    # --- Handle Output Based on Mode ---
    if args.output_mode == 'file':
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(f"# Intersection found from {len(args.filters)} files with the following filters:\n")
            for spec in args.filters:
                if spec['type'] == 'number':
                    f.write(f"# - {spec['path'].name}: Top {spec['value']} ligands\n")
                else:
                    f.write(f"# - {spec['path'].name}: Scores between {spec['value'][0]} and {spec['value'][1]}\n")
            f.write("#\n")
            for ligand in sorted(list(common_ligands)):
                f.write(f"{ligand}\n")
        print(f"List of intersecting ligands saved to: '{output_path}'")

    elif args.output_mode == 'directory':
        source_path = Path(args.source_dir)
        file_index = create_file_index(source_path)
        copy_intersecting_files(common_ligands, source_path, output_path, file_index)

if __name__ == "__main__":
    main()
