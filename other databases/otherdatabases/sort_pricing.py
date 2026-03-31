# 3/29/2026
import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(description="Calculate price-per-gram and sort chemical pricing data.")
parser.add_argument("input_csv", help="The raw CSV exported from Chemspace/MolPort")
parser.add_argument("output_csv", nargs='?', default="sorted_best_prices.csv", help="Sorted output file")
args = parser.parse_args()

print(f"Loading data from {args.input_csv}...")
df = pd.read_csv(args.input_csv)

# Note: You may need to change these column names to exactly match the headers in your exported CSV
price_col = 'Price' 
pack_size_col = 'Pack Size'

def calculate_price_per_gram(row):
    try:
        price_str = str(row[price_col]).replace('$', '').replace(',', '').strip()
        price = float(price_str)
        
        size_str = str(row[pack_size_col]).lower()
        
        # Extract the number and the unit (mg, g, kg)
        match = re.search(r'([\d\.]+)\s*(mg|g|kg)', size_str)
        if not match:
            return None
            
        amount = float(match.group(1))
        unit = match.group(2)
        
        # Convert everything to grams
        if unit == 'mg':
            mass_in_grams = amount / 1000.0
        elif unit == 'kg':
            mass_in_grams = amount * 1000.0
        else:
            mass_in_grams = amount
            
        return price / mass_in_grams
    except Exception:
        return None

print("Standardizing units and calculating price-per-gram...")
df['Price_Per_Gram'] = df.apply(calculate_price_per_gram, axis=1)

# Drop rows where we couldn't calculate a price (e.g., out of stock, call for quote)
clean_df = df.dropna(subset=['Price_Per_Gram'])

# Sort from cheapest to most expensive
sorted_df = clean_df.sort_values(by='Price_Per_Gram', ascending=True)

# Save the sorted data
sorted_df.to_csv(args.output_csv, index=False)
print(f"Success! Sorted pricing data saved to {args.output_csv}")