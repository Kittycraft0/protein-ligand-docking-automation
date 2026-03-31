# 3/29/2026


import pandas as pd
import glob
import re
import os

def process_pricing_files():
    # Load the SMILES to ZINC bridge first
    smiles_to_zinc = {}
    try:
        with open("final_smiles.txt", "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    smiles_to_zinc[parts[0]] = parts[1]
    except FileNotFoundError:
        print("[!] 'final_smiles.txt' not found. Cannot map IDs.")

    # Grab ALL possible spreadsheet formats, not just CSVs!
    all_files = []
    for ext in ["*.csv", "*.tsv", "*.xls", "*.xlsx"]:
        all_files.extend(glob.glob(ext))
        
    all_data = []

    for file in all_files:
        if any(keyword in file for keyword in ["Duplicate", "Unavailable", "Not Selected", "MASTER", "OMNI"]):
            continue
            
        try:
            # --- INTELLIGENT FILE LOADER ---
            if file.endswith('.xls') or file.endswith('.xlsx'):
                df = pd.read_excel(file)
            elif file.endswith('.tsv'):
                df = pd.read_csv(file, sep='\t')
            else:
                try:
                    df = pd.read_csv(file)
                    if len(df.columns) < 2:
                        df = pd.read_csv(file, sep=';')
                except Exception:
                    df = pd.read_csv(file, sep=';')
                    
            # --- 1. PROCESS MOLPORT FORMAT ---
            if 'Molport ID' in df.columns:
                print(f"[+] Processing MolPort file: {file}")
                for idx, row in df.iterrows():
                    price = row.get('Price, USD', None)
                    qty = row.get('Quantity', None)
                    unit = str(row.get('Packing', '')).lower().strip()
                    smiles = row.get('Search Criteria', '')
                    
                    if pd.isna(price) or pd.isna(qty): continue
                    
                    try:
                        qty = float(qty)
                        if unit == 'mg': grams = qty / 1000.0
                        elif unit == 'kg': grams = qty * 1000.0
                        else: grams = qty
                            
                        ppg = float(price) / grams if grams > 0 else None
                        
                        all_data.append({
                            'Source': 'MolPort',
                            'Vendor_Name': row.get('Supplier Name', 'MolPort Vendor'),
                            'Ligand_SMILES': smiles,
                            'Original_ID': smiles_to_zinc.get(smiles, 'Unknown ID'),
                            'Vendor_ID': row.get('Molport ID', ''),
                            'Pack_Size': f"{qty} {unit}",
                            'Price_USD': float(price),
                            'Price_Per_Gram': ppg
                        })
                    except ValueError: continue
                    
            # --- 2. PROCESS CHEMSPACE FORMAT ---
            elif 'Chemspace ID' in df.columns:
                print(f"[+] Processing Chemspace file: {file}")
                for idx, row in df.iterrows():
                    price_col = [c for c in df.columns if 'Cheapest Pack Price' in c]
                    price = row[price_col[0]] if price_col else None
                    pack_str = str(row.get('Cheapest Pack', '')).lower().strip()
                    smiles = row.get('SMILES', '')
                    
                    if pd.isna(price): continue
                    
                    match = re.search(r'([\d\.]+)\s*(mg|g|kg)', pack_str)
                    if not match: continue
                    
                    try:
                        qty = float(match.group(1))
                        unit = match.group(2)
                        
                        if unit == 'mg': grams = qty / 1000.0
                        elif unit == 'kg': grams = qty * 1000.0
                        else: grams = qty
                            
                        ppg = float(price) / grams if grams > 0 else None
                        
                        all_data.append({
                            'Source': 'Chemspace',
                            'Vendor_Name': 'Chemspace Aggregator',
                            'Ligand_SMILES': smiles,
                            'Original_ID': smiles_to_zinc.get(smiles, 'Unknown ID'),
                            'Vendor_ID': row.get('Chemspace ID', ''),
                            'Pack_Size': pack_str,
                            'Price_USD': float(price),
                            'Price_Per_Gram': ppg
                        })
                    except ValueError: continue
                        
            # --- 3. PROCESS ZINC CATALOGS (TSV) ---
            elif 'ZINC_ID' in df.columns and 'Vendors' in df.columns:
                print(f"[+] Processing ZINC Catalogs file: {file}")
                for idx, row in df.iterrows():
                    all_data.append({
                        'Source': 'ZINC API',
                        'Vendor_Name': row.get('Vendors', 'No Vendor Listed'),
                        'Ligand_SMILES': '', 
                        'Original_ID': row.get('ZINC_ID', ''), # Use ZINC_ID directly
                        'Vendor_ID': row.get('ZINC_ID', ''),
                        'Pack_Size': 'Unknown',
                        'Price_USD': None,      
                        'Price_Per_Gram': None
                    })

            # --- 4. PROCESS EMOLECULES FORMAT ---
            elif any(col in df.columns for col in ['eMolecules ID', 'eMolecules Number']):
                print(f"[+] Processing eMolecules file: {file}")
                
                price_col = next((c for c in df.columns if 'Price' in c or 'Cost' in c), None)
                qty_col = next((c for c in df.columns if 'Qty' in c or 'Amount' in c or 'Quantity' in c), None)
                unit_col = next((c for c in df.columns if 'Unit' in c), None)
                supplier_col = next((c for c in df.columns if 'Supplier' in c or 'Vendor' in c), 'eMolecules Vendor')
                id_col = next((c for c in df.columns if 'eMolecules' in c), '')
                smiles_col = next((c for c in df.columns if 'SMILES' in c.upper()), '')
                
                for idx, row in df.iterrows():
                    price = row.get(price_col, None) if price_col else None
                    qty = row.get(qty_col, None) if qty_col else None
                    unit = str(row.get(unit_col, 'g')).lower().strip() if unit_col else 'g'
                    smiles = row.get(smiles_col, '')
                    
                    if pd.isna(price) or pd.isna(qty): continue
                    
                    try:
                        if isinstance(price, str):
                            price = float(re.sub(r'[^\d\.]', '', price))
                            
                        qty = float(qty)
                        if unit == 'mg': grams = qty / 1000.0
                        elif unit == 'kg': grams = qty * 1000.0
                        else: grams = qty
                            
                        ppg = price / grams if grams > 0 else None
                        
                        all_data.append({
                            'Source': 'eMolecules',
                            'Vendor_Name': row.get(supplier_col, 'eMolecules Vendor'),
                            'Ligand_SMILES': smiles,
                            'Original_ID': smiles_to_zinc.get(smiles, 'Unknown ID'),
                            'Vendor_ID': row.get(id_col, ''),
                            'Pack_Size': f"{qty} {unit}",
                            'Price_USD': float(price),
                            'Price_Per_Gram': ppg
                        })
                    except (ValueError, TypeError):
                        continue
                        
            else:
                print(f"[-] Unrecognized format (Skipping): {file}")

        except Exception as e:
            print(f"[!] Error reading {file}: {e}")

    # --- 5. COMBINE, CLEAN, AND SORT ---
    if all_data:
        print("\nMerging and calculating best prices...")
        master_df = pd.DataFrame(all_data)
        
        # Sort so ZINC blank prices are at the bottom
        master_df = master_df.sort_values(by='Price_Per_Gram', ascending=True, na_position='last')
        
        master_df['Price_USD'] = master_df['Price_USD'].apply(lambda x: round(x, 2) if pd.notnull(x) else 'N/A')
        master_df['Price_Per_Gram'] = master_df['Price_Per_Gram'].apply(lambda x: round(x, 2) if pd.notnull(x) else 'N/A')
        
        # Save output so the dashboard can read it!
        output_name = "MASTER_PRICING_SORTED.csv"
        master_df.to_csv(output_name, index=False)
        print(f"Success! {len(master_df)} options combined and saved to '{output_name}'.")
    else:
        print("No valid data could be extracted.")

if __name__ == "__main__":
    process_pricing_files()