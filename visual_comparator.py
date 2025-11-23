import argparse
from pathlib import Path
import sys
import os
import subprocess
import time

def show_visual_comparison(chimerax_path: str, protein_file: Path, ref_ligand_file: Path, test_ligand_file: Path):
    """
    Opens a ChimeraX GUI window to display a side-by-side visual comparison of two poses.
    """
    # This block of text is the exact set of ChimeraX commands we confirmed works.
    script_content = f"""
# Step 0: Start fresh
close session

# Step 1: Open the three necessary files
open "{protein_file.as_posix()}"
open "{ref_ligand_file.as_posix()}"
open "{test_ligand_file.as_posix()}"

# Step 2: Style the scene for clarity
# Hide extra poses from the test ligand file and style the main ones
hide #3
show #3.1
cartoon #1; color #1&C gray
style #2 ball; color #2 gold transparency 60; graphics silhouettes #2 color gold width 2
style #3.1 ball; color #3.1 cyan

# Step 3: Highlight the binding pocket
select zone #2,#3.1 5 #1&protein residues true
show sel atoms
style sel stick
color sel&C tan

# Step 4: Finalize the view
# First, select the two ligands to center the view on them
select #2 #3.1
# Now, view that selection with padding and no clipping fog
view sel clip false pad 0.3
# Finally, clear the selection for a clean view
~select

# Optional: Add labels for clarity
label #2.1 name "Reference (ABA)" offset 0,0,1 height 0.8 color gold
label #3.1 name "Test Hit" offset 0,0,-1 height 0.8 color cyan
"""
    script_path = Path("temp_visualizer_script.cxc")
    with open(script_path, "w") as f:
        f.write(script_content)

    try:
        # Run ChimeraX WITH the GUI. This will open a window for you to see.
        # We use Popen so the python script can continue and wait for user input.
        print("Displaying comparison in ChimeraX... Close the ChimeraX window to continue to the next ligand.")
        process = subprocess.Popen([chimerax_path, script_path.as_posix()])
        process.wait()  # This makes the script pause until you close the ChimeraX window.
        
    except Exception as e:
        print(f"\nAn error occurred while running ChimeraX for {test_ligand_file.name}: {e}")
    finally:
        # Always clean up the temporary script file
        if script_path.exists():
            os.remove(script_path)


def visual_ranker(protein_file: Path, ref_ligand_file: Path, test_ligand_files: list):
    """
    Steps through a list of test ligands, showing a visual comparison for each.
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
    if not protein_file.is_file() or not ref_ligand_file.is_file():
        print(f"Error: A required file was not found. Check paths for protein and reference ligand.")
        sys.exit(1)
        
    print(f"\nStarting visual comparison slide show.")
    print(f"Reference pose is: '{ref_ligand_file.name}'")

    # --- Main Loop ---
    total_files = len(test_ligand_files)
    for i, hit_file in enumerate(test_ligand_files):
        print("\n" + "="*50)
        print(f"Displaying hit {i+1}/{total_files}: {hit_file.name}")
        
        show_visual_comparison(chimerax_path, protein_file, ref_ligand_file, hit_file)
        
        if i < total_files - 1:
             # Wait for the user to press Enter before loading the next one
             input("Press Enter to continue to the next ligand...")
        else:
             print("All ligands viewed.")


def main():
    parser = argparse.ArgumentParser(
        description="Interactively view docked poses for visual comparison against a reference.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('protein_file', help="The full path to the protein PDBQT file.")
    parser.add_argument('reference_file', help="The full path to the docked reference ligand PDBQT file (e.g., the ABA pose).")
    parser.add_argument('test_files', nargs='+', help="A list of all other docked ligand PDBQT files to compare.")
    
    args = parser.parse_args()

    visual_ranker(
        protein_file=Path(args.protein_file),
        ref_ligand_file=Path(args.reference_file),
        test_ligand_files=[Path(f) for f in args.test_files]
    )

if __name__ == "__main__":
    main()

