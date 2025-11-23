import argparse
from pathlib import Path
import sys
import os
import subprocess
import re
import time

def get_rmsd_from_chimerax(chimerax_path: str, protein_file: Path, ref_ligand_file: Path, test_ligand_file: Path):
    """
    Uses ChimeraX in nogui mode to calculate the RMSD between the best poses
    of two docked ligand files, within the context of the protein.
    """
    # Create a temporary command script for ChimeraX
    script_content = f"""
# Open all three components: the protein, the reference ligand, and the test ligand.
# They will be loaded as models #1, #2, and #3 respectively.
open "{protein_file.as_posix()}"
open "{ref_ligand_file.as_posix()}"
open "{test_ligand_file.as_posix()}"

# All models are now in the same coordinate system. We can directly
# calculate the RMSD between the best pose of the reference ligand (#2.1)
# and the best pose of the test ligand (#3.1).
# The `rmsd` command is correct here because we want to measure the
# distance between the atoms in their current positions without a new fit.
rmsd #2.1 to #3.1

# Exit ChimeraX
exit
"""
    script_path = Path("temp_rmsd_script.cxc")
    with open(script_path, "w") as f:
        f.write(script_content)

    # Run ChimeraX from the command line
    try:
        result = subprocess.run(
            [chimerax_path, "--nogui", script_path.as_posix()],
            capture_output=True, text=True, timeout=60 # 60-second timeout
        )
        os.remove(script_path)

        # Parse the output to find the RMSD value
        # The line we want looks like: "RMSD between #2.1 and #3.1, 19 atoms = 1.234 angstroms"
        for line in result.stdout.splitlines():
            match = re.search(r"RMSD between.*?=\s*([0-9\.]+)\s*angstroms", line)
            if match:
                return float(match.group(1))
        
        # For debugging: if RMSD is not found, show the output from ChimeraX
        # print(f"\n--- ChimeraX output for {test_ligand_file.name} ---")
        # print(result.stdout)
        # print("--------------------------------------------------")
        return None

    except Exception as e:
        print(f"\nAn error occurred while running ChimeraX for {test_ligand_file.name}: {e}")
        if script_path.exists():
            os.remove(script_path)
        return None


def rank_poses_by_rmsd(protein_file: Path, ref_ligand_file: Path, test_ligand_files: list, output_file: Path):
    """
    Calculates the RMSD of a list of ligands compared to a single reference ligand.
    """
    # --- Find ChimeraX executable ---
    possible_paths = [
        "C:/Program Files/ChimeraX 1.9/bin/ChimeraX.exe",
        "C:/Program Files/ChimeraX/bin/ChimeraX.exe"
    ]
    chimerax_path = next((p for p in possible_paths if Path(p).exists()), None)
    if not chimerax_path:
        print("Error: Could not find ChimeraX.exe. Please edit 'possible_paths' in the script.")
        sys.exit(1)
    
    print(f"Using ChimeraX found at: {chimerax_path}")

    # --- Validation ---
    if not protein_file.is_file():
        print(f"Error: Protein file '{protein_file}' not found.")
        sys.exit(1)
    if not ref_ligand_file.is_file():
        print(f"Error: Reference ligand file '{ref_ligand_file}' not found.")
        sys.exit(1)
        
    print(f"Using reference pose: '{ref_ligand_file.name}'")
    print(f"Comparing against {len(test_ligand_files)} other poses.")

    # --- Main Calculation Loop ---
    rmsd_results = []
    total_files = len(test_ligand_files)
    start_time = time.time()

    for i, hit_file in enumerate(test_ligand_files):
        # Don't compare a file to itself
        if hit_file.name == ref_ligand_file.name:
            continue

        progress_percent = (i + 1) / total_files * 100
        sys.stdout.write(f"\rProcessing file {i + 1}/{total_files} ({progress_percent:.1f}%)")
        sys.stdout.flush()

        # Pass the protein file along with the ligands to the function
        rmsd_val = get_rmsd_from_chimerax(chimerax_path, protein_file, ref_ligand_file, hit_file)
        
        if rmsd_val is not None:
            parent_folder_name = hit_file.parent.name
            rmsd_results.append((rmsd_val, parent_folder_name, hit_file.name))

    # --- Sort and Save Results ---
    print(f"\n\nCalculation complete in {time.time() - start_time:.2f} seconds. Sorting results...")
    rmsd_results.sort()

    with open(output_file, 'w') as f:
        f.write(f"# Poses ranked by structural RMSD to reference: {ref_ligand_file.name}\n")
        f.write(f"# Protein context: {protein_file.name}\n")
        f.write("# Lower RMSD = More Similar Pose\n")
        f.write("# RMSD (A)   Ligand_Folder                                 Source_File\n")
        for rmsd, folder, name in rmsd_results:
            f.write(f"{rmsd:<10.4f} {folder:<45} {name}\n")

    print(f"\nProcess complete. Ranked list saved to '{output_file}'.")
    print("\n--- Top 5 Most Similar Poses ---")
    for rmsd, folder, name in rmsd_results[:5]:
        print(f"  {rmsd:<10.4f} {folder}")


def main():
    parser = argparse.ArgumentParser(
        description="Ranks docked ligand poses by their structural RMSD to a reference ligand pose.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('protein_file', help="The full path to the protein PDBQT file used for the docking run.")
    parser.add_argument('reference_file', help="The full path to the single, docked reference ligand PDBQT file.")
    parser.add_argument('output_file', help="Name of the output text file for the ranked list.")
    parser.add_argument('test_files', nargs='+', help="A list of all other docked ligand PDBQT files to compare.")
    
    args = parser.parse_args()

    rank_poses_by_rmsd(
        protein_file=Path(args.protein_file),
        ref_ligand_file=Path(args.reference_file),
        test_ligand_files=[Path(f) for f in args.test_files],
        output_file=Path(args.output_file)
    )

if __name__ == "__main__":
    main()
