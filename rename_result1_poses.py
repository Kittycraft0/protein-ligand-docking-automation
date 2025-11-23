# 7/14/2025

import argparse
from pathlib import Path
import sys
import os
import shutil

def rename_and_copy_files(source_dir: Path, output_dir: Path):
    """
    Recursively finds all .pdbqt files in a source directory. It renames
    files matching a specific pattern and copies all files to an output directory.

    The renaming logic transforms a name like:
    'LIGAND.xaa_model_001_..._DOCK.pdbqt' into 'LIGAND.xaa_001.pdbqt'.
    
    Files that do not match the pattern are copied with their original names.

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
    processed_count = 0
    renamed_count = 0
    skipped_count = 0
    for i, file_path in enumerate(source_files):
        # Progress Indicator
        progress_percent = (i + 1) / total_files * 100
        sys.stdout.write(f"\rProcessing file {i + 1}/{total_files} ({progress_percent:.1f}%)")
        sys.stdout.flush()

        original_name = file_path.name
        new_name = original_name  # Default to the original name
        was_renamed = False

        try:
            # --- 3. Renaming Logic ---
            # Check if the specific renaming pattern '_model_' exists in the filename.
            if '_model_' in original_name:
                parts = original_name.split('_model_')
                # Ensure the split was successful and makes sense
                if len(parts) >= 2:
                    base_name = parts[0]
                    model_number_part = parts[1]
                    
                    # The model number should be the first part after the split
                    model_number = model_number_part.split('_')[0]

                    # Construct the new, simplified filename
                    new_name = f"{base_name}_{model_number}.pdbqt"
                    was_renamed = True
            
            # If the pattern was not found, new_name remains the same as original_name.

            # --- 4. Copy the file (with either the new or original name) ---
            source_file = file_path
            destination_file = output_dir / new_name
            
            shutil.copy2(source_file, destination_file)
            
            processed_count += 1
            if was_renamed:
                renamed_count += 1

        except Exception as e:
            print(f"\nERROR: Could not process file '{original_name}'. Reason: {e}")
            skipped_count += 1

    # --- 5. Final Summary ---
    print(f"\n\nProcess complete.")
    print(f"Successfully copied {processed_count} of {total_files} files.")
    print(f"- Renamed: {renamed_count} files.")
    print(f"- Copied without changes: {processed_count - renamed_count} files.")
    if skipped_count > 0:
        print(f"Skipped {skipped_count} files due to errors.")

def main():
    """
    Sets up argument parsing and runs the main file processing function.
    """
    parser = argparse.ArgumentParser(
        description="Copies and renames PDBQT files from a complex format to a simple one.\n"
                    "Example: 'LIGAND.xaa_model_001_..._DOCK.pdbqt' -> 'LIGAND.xaa_001.pdbqt'",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'source_directory', 
        type=str, 
        help="The source directory containing the original multi-model PDBQT files."
    )
    parser.add_argument(
        'output_directory', 
        type=str, 
        help="The directory where the renamed, single-pose PDBQT files will be saved."
    )
    args = parser.parse_args()

    # Convert string paths to Path objects for easier handling
    source_path = Path(args.source_directory)
    output_path = Path(args.output_directory)

    rename_and_copy_files(source_path, output_path)

if __name__ == "__main__":
    main()
