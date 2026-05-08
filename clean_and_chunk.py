# 5/8/2026
import re
from pathlib import Path
from tqdm import tqdm

# Input your Open Babel generated file
INPUT_FILE = Path("dock/upload_ready_smiles.txt")
CHUNK_SIZE = 10000

if not INPUT_FILE.exists():
    print(f"Error: Could not find {INPUT_FILE}")
    exit(1)

BASE_DIR = Path("./dock")
SMILES_CHUNKS_DIR = BASE_DIR / "smiles_chunks"

# Ensure necessary directories exist
for directory in [SMILES_CHUNKS_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

print("Scrubbing explicit brackets and chunking files...")

# This Regex specifically finds standard atoms trapped in brackets 
# without special charges or stereochem. It matches [C], [c], [N], [O], [Cl], etc.
bracket_pattern = re.compile(r'\[([CcnNoOsSFPClBrI]+)\]')

# Read all lines into memory (1.2 million lines is only ~60MB, so this is safe)
with open(INPUT_FILE, 'r') as f:
    lines = f.readlines()

current_chunk = 1
current_lines = []

for line in tqdm(lines):
    parts = line.split()
    if len(parts) >= 2:
        smiles = parts[0]
        zinc_id = parts[1]
        
        # Strip the fake brackets out, leaving normal atoms
        clean_smiles = bracket_pattern.sub(r'\1', smiles)
        
        current_lines.append(f"{clean_smiles}\t{zinc_id}\n")

    # If we hit 10,000 lines, save to a new file and reset
    if len(current_lines) == CHUNK_SIZE:
        # THE FIX: Use the Path object / filename
        out_name = SMILES_CHUNKS_DIR / f"cleaned_chunk_{current_chunk}.txt"
        with open(out_name, "w") as out:
            out.writelines(current_lines)
        current_chunk += 1
        current_lines = []

# Write any leftover lines at the very end
if current_lines:
    # THE FIX: Added the directory path here as well
    out_name = SMILES_CHUNKS_DIR / f"cleaned_chunk_{current_chunk}.txt"
    with open(out_name, "w") as out:
        out.writelines(current_lines)

print(f"\nSuccess! Cleaned the SMILES and generated {current_chunk} chunk files in {SMILES_CHUNKS_DIR}.")