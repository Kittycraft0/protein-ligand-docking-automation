# 3/29/2026

import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Format ZINC IDs to bypass the ZINC web uploader bug.")
parser.add_argument("input_file", nargs='?', default="zinc_ids_extracted.txt", help="Input text file with ZINC IDs")
parser.add_argument("output_file", nargs='?', default="zinc_ids_web_compatible.txt", help="Output text file for the web")
args = parser.parse_args()

try:
    with open(args.input_file, 'r') as f:
        # Read lines and strip whitespace
        zinc_ids = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(f"Error: Could not find '{args.input_file}'.")
    exit(1)

cleaned_ids = []
skipped_ids = []

for z_id in zinc_ids:
    # Strip out the "ZINC" letters
    stripped_str = z_id.upper().replace("ZINC", "")
    
    try:
        # Casting to int automatically drops the leading zeros
        clean_num = str(int(stripped_str))
        cleaned_ids.append(clean_num)
    except ValueError:
        # THE FIX: If it can't be turned into a number (like "ABA"), skip it entirely
        skipped_ids.append(z_id)
        continue

# Write the web-compatible IDs to the new file
with open(args.output_file, 'w') as f:
    for c_id in cleaned_ids:
        f.write(f"{c_id}\n")

print(f"Success! Formatted {len(cleaned_ids)} IDs. Saved to '{args.output_file}'.")
if skipped_ids:
    print(f"Skipped {len(skipped_ids)} non-numeric IDs: {', '.join(skipped_ids)}")