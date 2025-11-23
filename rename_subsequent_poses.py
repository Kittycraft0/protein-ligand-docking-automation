import argparse
from pathlib import Path
import sys
import os
import shutil

def rename_and_copy_files(source_dir: Path, output_dir: Path):
    """
    Recursively finds all .pdbqt files in a source directory. It renames
    files by removing protein and model-specific suffixes from the name.

    The renaming logic transforms a name like:
    'LIGAND.xaa_06_model_1_vs_PROTEIN_DOCK.pdbqt' into 'LIGAND.xaa_06.pdbqt'.
    
    It also prevents creating duplicate files if a simplified name already
    exists in the output directory.

    Args:
        source_dir: The directory to search for PDBQT files.
        output_dir: The directory where the renamed files will be saved.
    """
    # --- 1. Setup ---
    if not source_dir.is_dir():
        print(f"Error: Source directory not found at '{source_dir}'")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory created/ensured at: '{output_dir}'")

    # Use .rglob to find files recursively in all subdirectories
    source_files = list(source_dir.rglob("*.pdbqt"))
    total_files = len(source_files)
    if not source_files:
        print(f"No .pdbqt files found in '{source_dir}' or its subdirectories.")
        return

    print(f"\nFound {total_files} .pdbqt files. Starting processing...")
    
    # --- 2. Main Loop ---
    copied_count = 0
    renamed_count = 0
    duplicates_skipped_count = 0
    error_skipped_count = 0
    
    for i, file_path in enumerate(source_files):
        # Progress Indicator
        progress_percent = (i + 1) / total_files * 100
        sys.stdout.write(f"\rProcessing file {i + 1}/{total_files} ({progress_percent:.1f}%)")
        sys.stdout.flush()

        original_name = file_path.name

        try:
            # --- 3. Renaming Logic ---
            # Start with the filename without the .pdbqt extension
            # Note: removesuffix requires Python 3.9+
            base_name = original_name.removesuffix('.pdbqt')

            # First, remove the protein part if it exists
            if '_vs_' in base_name:
                base_name = base_name.split('_vs_')[0]

            # Next, remove the model part if it exists
            if '_model_' in base_name:
                base_name = base_name.split('_model_')[0]
            
            new_name = f"{base_name}.pdbqt"
            
            # --- 4. Copy the file, handling duplicates ---
            destination_file = output_dir / new_name
            
            # If the simplified filename already exists, skip it to avoid duplicates.
            if destination_file.exists():
                duplicates_skipped_count += 1
                continue # Skip to the next file in the loop

            shutil.copy2(file_path, destination_file)
            
            copied_count += 1
            if new_name != original_name:
                renamed_count += 1

        except Exception as e:
            print(f"\nERROR: Could not process file '{original_name}'. Reason: {e}")
            error_skipped_count += 1

    # --- 5. Final Summary ---
    print(f"\n\nProcess complete.")
    print(f"Scanned {total_files} source files.")
    print(f"- Copied {copied_count} files to the output directory:")
    print(f"  - Renamed: {renamed_count} files.")
    print(f"  - Copied without changes: {copied_count - renamed_count} files.")
    if duplicates_skipped_count > 0:
        print(f"- Skipped {duplicates_skipped_count} duplicate files (target name already existed).")
    if error_skipped_count > 0:
        print(f"- Skipped {error_skipped_count} files due to errors.")

def main():
    """
    Sets up argument parsing and runs the main file processing function.
    """
    parser = argparse.ArgumentParser(
        description="Copies and renames PDBQT files from a complex docking result format to a simple ligand name.\n"
                    "Example: 'LIGAND_MODEL_vs_PROTEIN_DOCK.pdbqt' -> 'LIGAND.pdbqt'",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'source_directory', 
        type=str, 
        help="The source directory containing the PDBQT files to be renamed."
    )
    parser.add_argument(
        'output_directory', 
        type=str, 
        help="The directory where the clean, renamed PDBQT files will be saved."
    )
    args = parser.parse_args()

    # Convert string paths to Path objects for easier handling
    source_path = Path(args.source_directory)
    output_path = Path(args.output_directory)

    rename_and_copy_files(source_path, output_path)

if __name__ == "__main__":
    main()
