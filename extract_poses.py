import argparse
from pathlib import Path
import sys
import os

def extract_first_model(source_dir: Path, output_dir: Path):
    """
    Scans a directory of PDBQT files, extracts the first MODEL from each,
    and saves it as a new single-pose PDBQT file. If a file has no MODEL
    records, it is copied as-is.

    Args:
        source_dir: The directory containing multi-model PDBQT files.
        output_dir: The directory where the new single-pose files will be saved.
    """
    # --- 1. Setup ---
    if not source_dir.is_dir():
        print(f"Error: Source directory not found at '{source_dir}'")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory created/ensured at: '{output_dir}'")

    source_files = list(source_dir.glob("*.pdbqt"))
    total_files = len(source_files)
    if not source_files:
        print(f"No .pdbqt files found in '{source_dir}'.")
        return

    print(f"\nProcessing {total_files} files from '{source_dir.name}'...")
    
    # --- 2. Main Loop ---
    for i, file_path in enumerate(source_files):
        # Progress Indicator
        progress_percent = (i + 1) / total_files * 100
        sys.stdout.write(f"\rProcessing file {i + 1}/{total_files} ({progress_percent:.1f}%) - {file_path.name}")
        sys.stdout.flush()

        try:
            with open(file_path, 'r') as f:
                content = f.read()

            output_content = []
            
            # --- 3. Extraction Logic ---
            if "MODEL" in content:
                in_first_model = False
                # Split the file into lines to process it
                for line in content.splitlines(True): # Keep newlines
                    if line.strip().startswith("MODEL 1"):
                        in_first_model = True
                    
                    if in_first_model:
                        output_content.append(line)
                    
                    if in_first_model and line.strip() == "ENDMDL":
                        break # Stop after finding the end of the first model
            else:
                # If no MODEL keyword, it's already a single-pose file
                output_content = content

            # --- 4. Save the new file ---
            if output_content:
                new_file_path = output_dir / file_path.name
                with open(new_file_path, 'w') as f_out:
                    # If it was a list of lines, join them
                    if isinstance(output_content, list):
                        f_out.writelines(output_content)
                    else: # Otherwise it was a single string
                        f_out.write(output_content)
            
        except Exception as e:
            print(f"\nERROR: Could not process file '{file_path.name}'. Reason: {e}")

    print(f"\n\nProcess complete. Extracted first model from {total_files} files into '{output_dir}'.")


def main():
    parser = argparse.ArgumentParser(
        description="Extracts the first MODEL (best pose) from multi-model PDBQT files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'source_directory', 
        type=str, 
        help="The input directory containing your multi-model top hit PDBQT files."
    )
    parser.add_argument(
        'output_directory', 
        type=str, 
        help="The new directory where the clean, single-pose PDBQT files will be saved."
    )
    args = parser.parse_args()

    extract_first_model(Path(args.source_directory), Path(args.output_directory))

if __name__ == "__main__":
    main()

