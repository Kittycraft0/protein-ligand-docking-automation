# File: get_both_boxes.py
# This script calculates and draws BOTH the tight and buffered docking boxes
# around the current selection in one go.

from chimerax.core.commands import run
import numpy
import os

# --- START: User Configuration ---
BUFFER_SIZE = 15.0      # Total size to add to each dimension for the larger box
TIGHT_BOX_COLOR = "gold"
BUFFERED_BOX_COLOR = "cyan"
# --- END: User Configuration ---


def draw_wireframe_box(session, file_name, size, center, color):
    """A helper function to write and open a BILD file for a box."""
    script_dir = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(script_dir, file_name)

    half_size = size / 2.0
    min_c = center - half_size
    max_c = center + half_size

    v = [
        (min_c[0], min_c[1], min_c[2]), (max_c[0], min_c[1], min_c[2]),
        (min_c[0], max_c[1], min_c[2]), (max_c[0], max_c[1], min_c[2]),
        (min_c[0], min_c[1], max_c[2]), (max_c[0], min_c[1], max_c[2]),
        (min_c[0], max_c[1], max_c[2]), (max_c[0], max_c[1], max_c[2])
    ]
    edges = [
        (0,1), (0,2), (0,4), (1,3), (1,5), (2,3),
        (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)
    ]

    with open(file_path, "w") as f:
        f.write(f".color {color}\n")
        for edge in edges:
            start_v, end_v = v[edge[0]], v[edge[1]]
            f.write(f".move {start_v[0]} {start_v[1]} {start_v[2]}\n")
            f.write(f".draw {end_v[0]} {end_v[1]} {end_v[2]}\n")
            
    run(session, f"open \"{file_path}\"")
    print(f"SUCCESS: Drew '{file_name}'")


# --- Main script logic starts here ---
current_selection_object = run(session, "select sel")
selected_residues = current_selection_object.residues

if not selected_residues:
    print("\nERROR: No residues are selected.\n")
else:
    all_atoms_in_selection = [a for res in selected_residues for a in res.atoms]
    
    if not all_atoms_in_selection:
        print("\nERROR: The selection contains no atoms.\n")
    else:
        # --- Perform calculations (only need to do this once) ---
        coords = numpy.array([a.coord for a in all_atoms_in_selection])
        center = coords.mean(axis=0)
        span = coords.max(axis=0) - coords.min(axis=0)
        
        tight_box_size = span
        buffered_box_size = span + BUFFER_SIZE

        # --- Draw BOTH boxes using the helper function ---
        draw_wireframe_box(session, "tight_box.bild", tight_box_size, center, TIGHT_BOX_COLOR)
        draw_wireframe_box(session, "buffered_box.bild", buffered_box_size, center, BUFFERED_BOX_COLOR)

        # --- Print BOTH sets of parameters to the Log ---
        print("\n--- Tight Box Parameters (for visualization) ---")
        print(f"center_x = {center[0]:.3f}, center_y = {center[1]:.3f}, center_z = {center[2]:.3f}")
        print(f"size_x = {tight_box_size[0]:.3f}, size_y = {tight_box_size[1]:.3f}, size_z = {tight_box_size[2]:.3f}")
        print("--------------------------------------------------")

        print("\n--- Docking Box Parameters (for Vina conf.txt) ---")
        print(f"center_x = {center[0]:.3f}")
        print(f"center_y = {center[1]:.3f}")
        print(f"center_z = {center[2]:.3f}")
        print("")
        print(f"size_x = {buffered_box_size[0]:.3f}")
        print(f"size_y = {buffered_box_size[1]:.3f}")
        print(f"size_z = {buffered_box_size[2]:.3f}")
        print("--------------------------------------------------\n")