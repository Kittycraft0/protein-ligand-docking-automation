# 3/29/2026
import requests
import concurrent.futures
import time
import argparse
import sys

parser = argparse.ArgumentParser(description="Fetch Purchasing Catalogs from ZINC API")
parser.add_argument("input_file", nargs='?', default="zinc_ids_extracted.txt", help="Input file with ZINC IDs")
parser.add_argument("output_file", nargs='?', default="zinc_catalogs.tsv", help="Output TSV file")
args = parser.parse_args()

try:
    with open(args.input_file, 'r') as f:
        zinc_ids = [line.strip() for line in f if line.strip() and line.startswith("ZINC")]
except FileNotFoundError:
    print(f"Error: Could not find '{args.input_file}'.")
    sys.exit(1)

total = len(zinc_ids)
print(f"Found {total} valid ZINC IDs. Fetching catalog data...\n")

start_time = time.time()
results = []
HEADERS = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'}

def get_catalogs(z_id):
    # Querying the specific catalogs text endpoint
    url = f"https://zinc15.docking.org/substances/{z_id}/catalogs.txt"
    
    for attempt in range(3):
        try:
            resp = requests.get(url, headers=HEADERS, timeout=10)
            if resp.status_code == 200:
                # The API returns multiple lines. We just extract the vendor names.
                lines = resp.text.strip().split('\n')
                vendors = list(set([line.split('\t')[0] for line in lines if line]))
                vendor_str = ", ".join(vendors) if vendors else "No Vendors Listed"
                return z_id, vendor_str
            time.sleep(1)
        except Exception:
            time.sleep(2)
            
    return z_id, "ERROR: Connection Failed"

# Keeping workers at 2 to stay under the Cloudflare radar
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
    futures = {executor.submit(get_catalogs, z_id): z_id for z_id in zinc_ids}
    
    for i, future in enumerate(concurrent.futures.as_completed(futures)):
        z_id, catalogs = future.result()
        results.append(f"{z_id}\t{catalogs}")
        
        current = i + 1
        percent = (current / total) * 100
        elapsed = time.time() - start_time
        time_per_request = elapsed / current if current > 0 else 0
        eta_seconds = int(time_per_request * (total - current))
        eta_str = time.strftime('%H:%M:%S', time.gmtime(eta_seconds))
        
        bar_len = int(50 * current // total)
        bar = '█' * bar_len + '-' * (50 - bar_len)
        print(f'\rFetching: |{bar}| {percent:.1f}% | ETA: {eta_str}', end='', flush=True)

# Save as a Tab-Separated Values (TSV) file so it opens cleanly in Excel
with open(args.output_file, 'w') as f:
    f.write("ZINC_ID\tVendors\n") # Header
    for res in results:
        f.write(f"{res}\n")

print(f"\n\nDone! Saved catalog data to '{args.output_file}'.")