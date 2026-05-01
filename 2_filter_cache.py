# 5/1/2026
# ai
import os
import csv
import shutil
from pathlib import Path

def purge_unwanted_cache_files(cache_dir, approved_csv):
    """Deletes any .pdbqt file in the cache that isn't listed in the approved CSV."""
    cache_path = Path(cache_dir)
    
    if not Path(approved_csv).exists():
        print(f"[!] Error: Approved list '{approved_csv}' not found.")
        return

    # 1. Load the approved files into a Set for fast lookup
    approved_files = set()
    with open(approved_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # The CSV saved the relative path (e.g., 'models_ACBBMM/ACBBMM_model_1.pdbqt')
            approved_files.add(row["File_Path"])
            
    print(f"[+] Loaded {len(approved_files)} approved molecules from your list.")

    # 2. Iterate through the cache and delete unapproved files
    all_models = list(cache_path.glob("models_*/*.pdbqt"))
    deleted_count = 0
    kept_count = 0
    
    print("Purging expensive/unwanted models from cache...")
    
    for model_file in all_models:
        relative_path = model_file.relative_to(cache_path).as_posix()
        
        if relative_path not in approved_files:
            try:
                os.remove(model_file)
                deleted_count += 1
            except Exception as e:
                print(f"[-] Could not delete {model_file.name}: {e}")
        else:
            kept_count += 1

    # 3. Clean up empty directories
    for folder in cache_path.glob("models_*"):
        if folder.is_dir() and not any(folder.iterdir()):
            shutil.rmtree(folder)

    print(f"\n[+] Cleanup Complete!")
    print(f"    Kept: {kept_count} files")
    print(f"    Deleted: {deleted_count} files")
    print("\nYour cache is now filtered by price. You can safely run dock.py!")

if __name__ == "__main__":
    CACHE_DIRECTORY = "dock/cache"
    APPROVED_CSV = "cache_ranked_by_price.csv" # The file you manually trimmed
    
    purge_unwanted_cache_files(CACHE_DIRECTORY, APPROVED_CSV)