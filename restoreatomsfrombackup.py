import argparse
from pathlib import Path
import sys
import shutil
import os

def restore_pdbqt_backups(backup_dir: Path, cache_dir: Path):
    """
    Restores PDBQT files from a backup directory to their original
    locations within the cache directory structure.
    """
    if not backup_dir.is_dir():
        print(f"Error: Backup directory not found at '{backup_dir}'")
        return

    if not cache_dir.is_dir():
        print(f"Error: Cache directory not found at '{cache_dir}'. Cannot restore.")
        return

    print(f"Restoring files from: {backup_dir}")
    print(f"Target cache directory: {cache_dir}\n")

    backup_files = list(backup_dir.glob("*.pdbqt"))
    total_files = len(backup_files)
    if not backup_files:
        print("No backup files (.pdbqt) found to restore.")
        return

    restored_count = 0
    for i, backup_file_path in enumerate(backup_files):
        # Print progress indicator
        progress_percent = (i + 1) / total_files * 100
        sys.stdout.write(f"\rProcessing file {i + 1}/{total_files} ({progress_percent:.1f}%)")
        sys.stdout.flush()

        # --- Logic to determine original location ---
        # Example filename: ACEBRN.xaa(1)_model_0463.pdbqt
        # We need to extract the parent name: "ACEBRN.xaa(1)"
        
        filename = backup_file_path.name
        parent_ligand_name = ""

        if "_model_" in filename:
            parent_ligand_name = filename.rsplit("_model_", 1)[0]
        else:
            # Handle cases where the ligand had no models
            parent_ligand_name = backup_file_path.stem

        if not parent_ligand_name:
            print(f"\nCould not determine parent ligand for '{filename}'. Skipping.")
            continue

        # Reconstruct the original destination directory path
        destination_dir = cache_dir / f"models_{parent_ligand_name}"
        
        if not destination_dir.is_dir():
            # This can happen if the original cache folder was deleted.
            # We'll recreate it to be safe.
            destination_dir.mkdir(parents=True, exist_ok=True)

        destination_file_path = destination_dir / filename
        
        # Copy the file back, overwriting the destination
        try:
            shutil.copy(backup_file_path, destination_file_path)
            restored_count += 1
        except Exception as e:
            print(f"\nFailed to restore '{filename}' to '{destination_file_path}'. Error: {e}")

    print(f"\n\nRestore complete. {restored_count}/{total_files} files were successfully restored.")


def main():
    # Get the directory where this script is located
    script_dir = Path(os.path.dirname(os.path.realpath(__file__)))

    parser = argparse.ArgumentParser(
        description="Restores PDBQT files from the backup folder to their original locations in the cache.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--backup_dir', 
        type=str, 
        default=str(script_dir / "pdbqt_backups"),
        help="Path to the backup directory (defaults to 'pdbqt_backups' in the script's folder)."
    )
    parser.add_argument(
        '--cache_dir', 
        type=str, 
        default=str(script_dir / "dock" / "cache"),
        help="Path to the main cache directory (defaults to 'dock/cache' in the script's folder)."
    )
    args = parser.parse_args()
    
    restore_pdbqt_backups(Path(args.backup_dir), Path(args.cache_dir))


if __name__ == "__main__":
    main()
