# 6/7/2025
import argparse
from pathlib import Path
import sys
import shutil
import os

VALID_AUTODOCK_TYPES = {
    'A', 'C', 'NA', 'N', 'OA', 'O', 'SA', 'S', 
    'HD', 'H', 'P', 'F', 'Cl', 'Br', 'I', 
    'Mg', 'Mn', 'Zn', 'Ca', 'Fe'
}

SUGGESTED_REPLACEMENTS = {
    'Si': 'C'
}

def check_and_fix_pdbqt_files(directory: Path, perform_fix: bool = False):
    """
    Scans a directory of PDBQT files for invalid atom types and optionally fixes them,
    backing up original files before modification using a robust replacement method.
    """
    print(f"Recursively scanning directory: {directory}\n")
    found_bad_files = False

    script_location = Path(os.path.dirname(os.path.realpath(__file__)))
    backup_dir = script_location / "pdbqt_backups"
    if perform_fix:
        backup_dir.mkdir(exist_ok=True)
        print(f"Original files will be backed up to: {backup_dir}\n")

    pdbqt_files = list(directory.glob("**/*.pdbqt"))
    total_files = len(pdbqt_files)
    if not pdbqt_files:
        print(f"No .pdbqt files found in {directory} or its subdirectories.")
        return

    for i, file_path in enumerate(pdbqt_files):
        if (i + 1) % 10 == 0 or i == total_files - 1:
            progress_percent = (i + 1) / total_files * 100
            sys.stdout.write(f"\rScanning file {i + 1}/{total_files} ({progress_percent:.1f}%)")
            sys.stdout.flush()

        is_bad_file = False
        lines_to_write = []
        
        try:
            with open(file_path, 'r') as f:
                original_lines = f.readlines()

            for line_num, line in enumerate(original_lines):
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    parts = line.split()
                    if len(parts) < 2: continue
                    
                    atom_type = parts[-1]
                    
                    if atom_type not in VALID_AUTODOCK_TYPES:
                        if not is_bad_file:
                            sys.stdout.write("\n")
                        is_bad_file = True
                        found_bad_files = True
                        
                        relative_path = file_path.relative_to(directory)
                        print(f"--- Invalid Atom Found in: {relative_path} ---")
                        print(f"  Line {line_num+1}: Found type '{atom_type}'")
                        
                        if perform_fix and atom_type in SUGGESTED_REPLACEMENTS:
                            replacement = SUGGESTED_REPLACEMENTS[atom_type]
                            print(f"  Action: Replacing '{atom_type}' with '{replacement}'.")
                            
                            # --- NEW, SAFER REPLACEMENT LOGIC ---
                            # This preserves original spacing by only replacing the last word.
                            # It splits the line only once, from the right, at the atom_type.
                            line_start = line.rsplit(atom_type, 1)[0]
                            new_line = line_start + replacement + "\n"
                            lines_to_write.append(new_line)
                            # ------------------------------------

                        else:
                            lines_to_write.append(line)
                            if atom_type in SUGGESTED_REPLACEMENTS:
                                print(f"  Suggestion: Re-run with the --fix flag to replace with '{SUGGESTED_REPLACEMENTS[atom_type]}'.")
                            else:
                                print("  Suggestion: No automatic replacement available. You may need to edit this file manually.")
                        print("-" * (25 + len(str(relative_path))))
                        
                    else:
                        lines_to_write.append(line)
                else:
                    lines_to_write.append(line)

            if is_bad_file and perform_fix:
                backup_file_path = backup_dir / file_path.name
                if not backup_file_path.exists():
                    shutil.copy(file_path, backup_file_path)
                    print(f"  INFO: Original file backed up to '{backup_file_path.name}'.")
                
                with open(file_path, 'w') as f:
                    f.writelines(lines_to_write)
                print(f"SUCCESS: Corrected and overwrote '{file_path.name}'.\n")

        except Exception as e:
            sys.stdout.write("\n")
            print(f"Could not process file {file_path.name}. Error: {e}")

    sys.stdout.write("\n")
    if not found_bad_files:
        print("Scan complete. All PDBQT files appear to have valid atom types!")
    else:
        print("Scan complete. Found invalid files listed above.")


def main():
    parser = argparse.ArgumentParser(
        description="Scans and validates AutoDock Vina PDBQT files for correct atom types.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'directory', 
        type=str, 
        help="The directory containing the .pdbqt files to check (e.g., 'dock/cache')."
    )
    parser.add_argument(
        '--fix', 
        action='store_true', 
        help="If set, automatically replaces known bad atom types and backs up the originals."
    )
    args = parser.parse_args()

    target_dir = Path(args.directory)
    if not target_dir.is_dir():
        print(f"Error: Directory not found at '{target_dir}'")
        sys.exit(1)
        
    check_and_fix_pdbqt_files(target_dir, args.fix)

if __name__ == "__main__":
    main()


# to run: python check_cache.py dock/cache