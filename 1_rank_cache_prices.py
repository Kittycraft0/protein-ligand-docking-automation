# 5/1/2026
# ai
import os
import subprocess
import pandas as pd
import csv
from pathlib import Path

def load_pricing(pricing_file):
    """Loads SMILES to Price mapping from the master pricing spreadsheet."""
    price_map = {}
    if Path(pricing_file).exists():
        df = pd.read_csv(pricing_file)
        df_valid = df[df['Price_Per_Gram'] != 'N/A']
        for _, row in df_valid.iterrows():
            # Ensure we don't trip over empty SMILES
            smiles = str(row.get('Ligand_SMILES', '')).strip()
            if smiles:
                price_map[smiles] = float(row['Price_Per_Gram'])
    else:
        print(f"[!] Pricing file '{pricing_file}' not found. Cannot map prices.")
    return price_map

def get_smiles_via_obabel(pdbqt_file):
    """Uses OpenBabel to convert a PDBQT file to a SMILES string."""
    try:
        result = subprocess.run(
            ['obabel', '-i', 'pdbqt', str(pdbqt_file), '-o', 'smi'],
            capture_output=True, text=True, check=True
        )
        # obabel outputs "SMILES filename" or just "SMILES". We just want the SMILES part.
        output = result.stdout.strip().split()[0] 
        return output
    except Exception:
        return None

def scan_cache_and_rank(cache_dir, output_csv, price_map):
    cache_path = Path(cache_dir)
    progress_file = Path("price_scan_progress.txt")
    
    # Bluescreen resistance: Load already processed files
    completed_files = set()
    if progress_file.exists():
        with open(progress_file, "r") as f:
            completed_files = set(line.strip() for line in f)

    # Find all extracted models in the cache
    all_models = list(cache_path.glob("models_*/*.pdbqt"))
    total_models = len(all_models)
    print(f"Found {total_models} models in the cache to analyze.")

    file_exists = Path(output_csv).exists()
    
    with open(output_csv, "a", newline='') as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:
            writer.writerow(["File_Path", "SMILES", "Price_Per_Gram"])
            
        for i, model_file in enumerate(all_models):
            relative_path = model_file.relative_to(cache_path).as_posix()
            
            if relative_path in completed_files:
                continue
                
            print(f"[{i+1}/{total_models}] Processing {model_file.name}...")
            
            smiles = get_smiles_via_obabel(model_file)
            price = "N/A"
            
            if smiles:
                # Find the exact or closest match in our price dictionary
                price = price_map.get(smiles, "N/A")
                
            writer.writerow([relative_path, smiles, price])
            csvfile.flush() # Ensure it writes immediately in case of crash
            
            # Record progress
            with open(progress_file, "a") as pf:
                pf.write(f"{relative_path}\n")

    print(f"\n[+] Done! Output saved to '{output_csv}'.")

if __name__ == "__main__":
    CACHE_DIRECTORY = "dock/cache"
    PRICING_FILE = "MASTER_PRICING_SORTED.csv"
    OUTPUT_CSV = "cache_ranked_by_price.csv"
    
    prices = load_pricing(PRICING_FILE)
    if prices:
        scan_cache_and_rank(CACHE_DIRECTORY, OUTPUT_CSV, prices)