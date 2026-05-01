# 5/1/2026
# ai
import os
import re
import csv
import pandas as pd
from pathlib import Path

def load_pricing(pricing_file):
    """Loads the pricing data into a fast dictionary lookup."""
    price_map = {}
    if Path(pricing_file).exists():
        df = pd.read_csv(pricing_file)
        # Filter out N/A and create mapping: ZINC_ID -> Price
        df_valid = df[df['Price_Per_Gram'] != 'N/A']
        for _, row in df_valid.iterrows():
            price_map[str(row['Original_ID']).strip()] = float(row['Price_Per_Gram'])
        print(f"[+] Loaded {len(price_map)} priced molecules from {pricing_file}")
    else:
        print(f"[!] Pricing file '{pricing_file}' not found. Will output 'N/A' for prices.")
    return price_map

def scan_bulk_files(source_dir, output_csv, price_map):
    source_path = Path(source_dir)
    progress_file = source_path / "scan_progress.txt"
    
    # Bluescreen resistance: load already scanned files
    completed_files = set()
    if progress_file.exists():
        with open(progress_file, "r") as f:
            completed_files = set(line.strip() for line in f)

    all_bulk_files = list(source_path.glob("*.pdbqt"))
    
    # Setup CSV Writer
    file_exists = Path(output_csv).exists()
    with open(output_csv, "a", newline='') as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:
            writer.writerow(["Bulk_File", "Model_Number", "Molecule_ID", "Price_Per_Gram"])
            
        print(f"Scanning {len(all_bulk_files)} bulk files...")
        
        for i, bulk_file in enumerate(all_bulk_files):
            if bulk_file.name in completed_files:
                continue
                
            print(f"Scanning [{i+1}/{len(all_bulk_files)}]: {bulk_file.name}")
            
            current_model = "Unknown"
            current_zinc = "Unknown"
            
            with open(bulk_file, "r") as f:
                for line in f:
                    # Find Model Number
                    if line.startswith("MODEL"):
                        current_model = line.split()[1].strip()
                        current_zinc = "Unknown" # Reset for new model
                    
                    # Find ZINC ID (Usually in a REMARK line)
                    elif line.startswith("REMARK") and "ZINC" in line:
                        match = re.search(r'(ZINC[a-zA-Z0-9\-]+)', line)
                        if match:
                            current_zinc = match.group(1)
                    
                    # End of Model: Time to record it
                    elif line.startswith("ENDMDL"):
                        price = price_map.get(current_zinc, "N/A")
                        writer.writerow([bulk_file.name, current_model, current_zinc, price])
            
            # Update progress file
            with open(progress_file, "a") as pf:
                pf.write(f"{bulk_file.name}\n")

    print(f"\n[+] Scan complete. Ranked list saved to '{output_csv}'.")
    print("-> Open it in Excel, sort by price, delete the ones you don't want, and save it for the next script.")

if __name__ == "__main__":
    # CONFIGURATION
    SOURCE_DIRECTORY = "dock/unzipped" # Where your bulk .pdbqt files are
    PRICING_FILE = "MASTER_PRICING_SORTED.csv"
    OUTPUT_CSV = "pre_docking_price_ranks.csv"
    
    prices = load_pricing(PRICING_FILE)
    scan_bulk_files(SOURCE_DIRECTORY, OUTPUT_CSV, prices)