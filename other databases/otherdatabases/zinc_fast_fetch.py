# 3/29/2026

import requests
import concurrent.futures
import time
import argparse
import sys

parser = argparse.ArgumentParser(description="Fast Multithreaded ZINC API Fetcher")
parser.add_argument("input_file", nargs='?', default="zinc_ids_extracted.txt", help="Input text file with ZINC IDs")
parser.add_argument("output_file", nargs='?', default="final_smiles.txt", help="Output text file for SMILES")
args = parser.parse_args()

try:
    with open(args.input_file, 'r') as f:
        zinc_ids = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(f"Error: Could not find '{args.input_file}'.")
    sys.exit(1)

total = len(zinc_ids)
print(f"Found {total} unique ZINC IDs. Firing up the API fetcher...\n")

start_time = time.time()
successful_smiles = []
failed_ids = []

# The magic key. Without this User-Agent, Cloudflare blocks Python instantly.
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
}

def get_smiles(z_id):
    url = f"https://zinc20.docking.org/substances/{z_id}.smi"
    
    # Try up to 3 times per molecule if the server drops us
    for attempt in range(3):
        try:
            resp = requests.get(url, headers=HEADERS, timeout=10)
            if resp.status_code == 200:
                return z_id, resp.text.split()[0], None
            
            # If we get an HTTP error (like 429 Too Many Requests), wait and retry
            time.sleep(1)
        except Exception as e:
            # If it's a hard connection reset, back off for 2 seconds and retry
            time.sleep(2)
            
    return z_id, None, "Connection forcibly closed after 3 retries"

# --- Replace your ThreadPoolExecutor line with this to throttle the concurrency ---
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    futures = {executor.submit(get_smiles, z_id): z_id for z_id in zinc_ids}
    
    for i, future in enumerate(concurrent.futures.as_completed(futures)):
        z_id, smiles, error = future.result()
        
        if smiles:
            # We are writing BOTH the smiles and the ID so Chemspace knows what is what
            successful_smiles.append(f"{smiles}\t{z_id}")
        else:
            failed_ids.append(f"{z_id} (Error: {error})")
        
        current = i + 1
        percent = (current / total) * 100
        
        elapsed = time.time() - start_time
        time_per_request = elapsed / current if current > 0 else 0
        eta_seconds = int(time_per_request * (total - current))
        eta_str = time.strftime('%H:%M:%S', time.gmtime(eta_seconds))
        
        bar_len = int(50 * current // total)
        bar = '█' * bar_len + '-' * (50 - bar_len)
        
        print(f'\rFetching: |{bar}| {percent:.1f}% | ETA: {eta_str}', end='', flush=True)

with open(args.output_file, 'w') as f:
    for s in successful_smiles:
        f.write(f"{s}\n")

print(f"\n\nDone! Saved {len(successful_smiles)} valid SMILES strings to '{args.output_file}'.")

# If it fails, print out exactly why so we can debug it
if failed_ids:
    print(f"\nFailed to fetch {len(failed_ids)} IDs. First 5 errors:")
    for err in failed_ids[:5]:
        print(f" - {err}")