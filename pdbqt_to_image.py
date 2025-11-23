#!/usr/bin/env python3
import os
import subprocess
import argparse
import sys

def generate_output_name(pdbqt_filename):
    """
    Generates the output PNG filename based on the specified renaming rule.
    Example: CCEBRN.xaa_0725.pdbqt -> CCEBRNxaa0725.png
    
    Args:
        pdbqt_filename (str): The base name of the input .pdbqt file.
        
    Returns:
        str: The generated .png filename.
    """
    # Remove the .pdbqt extension to work with the base name
    base_name = pdbqt_filename.replace('.pdbqt', '')
    
    # Remove the "." and "_" characters from the base name
    renamed_base = base_name.replace('.', '').replace('_', '')
    
    # Return the new name with the .png extension
    return f"{renamed_base}.png"

def process_file(file_path, output_dir=None):
    """
    Runs the obabel command on a single .pdbqt file.
    
    Args:
        file_path (str): The full path to the .pdbqt file.
        output_dir (str, optional): The directory to save the output PNG. 
                                    If None, saves in the same directory as the input file.
    """
    # Ensure the file is a .pdbqt file before processing
    if not file_path.lower().endswith('.pdbqt'):
        print(f"Skipping non-pdbqt file: {os.path.basename(file_path)}")
        return

    # Get the directory and filename from the full path
    input_directory, filename = os.path.split(file_path)
    
    # Determine the final output directory
    # If an output_dir is provided, use it. Otherwise, use the input file's directory.
    final_output_dir = output_dir if output_dir else input_directory
    
    # Generate the desired output filename
    output_filename = generate_output_name(filename)
    output_path = os.path.join(final_output_dir, output_filename)
    
    # Construct the full obabel command
    command = [
        'obabel',
        file_path,
        '-o', 'png',
        '-O', output_path
    ]
    
    print(f"Processing: {filename} -> {output_path}")
    
    try:
        # Execute the command
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"Successfully created: {output_path}")
    except FileNotFoundError:
        print("\nError: 'obabel' command not found.", file=sys.stderr)
        print("Please ensure Open Babel is installed and accessible in your system's PATH.", file=sys.stderr)
        sys.exit(1) # Exit the script if obabel is not found
    except subprocess.CalledProcessError as e:
        # This error occurs if obabel returns a non-zero exit code (i.e., an error)
        print(f"Error processing {filename}:", file=sys.stderr)
        print(e.stderr, file=sys.stderr)

def main():
    """
    Main function to parse arguments and orchestrate file processing.
    """
    parser = argparse.ArgumentParser(
        description="Batch convert .pdbqt files to .png using Open Babel with specific renaming.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "input_path",
        help="Path to a single .pdbqt file or a directory containing .pdbqt files."
    )
    parser.add_argument(
        "-o", "--output_dir",
        help="Path to the output directory. If not specified, PNGs are saved in the input directory.\nThe directory will be created if it doesn't exist."
    )
    
    args = parser.parse_args()
    input_path = args.input_path
    output_dir = args.output_dir

    # Check if the provided input path exists
    if not os.path.exists(input_path):
        print(f"Error: The input path '{input_path}' does not exist.", file=sys.stderr)
        return

    # If an output directory is specified, check if it exists and create it if not.
    if output_dir and not os.path.isdir(output_dir):
        print(f"Output directory '{output_dir}' not found. Creating it.")
        os.makedirs(output_dir)

    # --- Process based on whether the path is a file or directory ---
    if os.path.isfile(input_path):
        process_file(input_path, output_dir)
    elif os.path.isdir(input_path):
        print(f"Scanning directory: {input_path}")
        # Walk through the directory and process each .pdbqt file found
        for filename in os.listdir(input_path):
            file_path = os.path.join(input_path, filename)
            if os.path.isfile(file_path):
                process_file(file_path, output_dir)
    else:
        print(f"Error: The path '{input_path}' is not a valid file or directory.", file=sys.stderr)

if __name__ == "__main__":
    main()
