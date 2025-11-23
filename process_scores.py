import datetime
import sys
import os # Import the 'os' module to handle file paths correctly

def process_and_sort_data(input_text):
    """
    Parses docking score data, groups scores by model name, sorts them,
    and returns a formatted string for output.

    Args:
        input_text: A string containing the raw data, with each line
                    formatted as "[score] [model_name]".

    Returns:
        A string containing the formatted and sorted data, ready to be
        written to a file.
    """
    # Dictionary to hold the data: { 'name': [score1, score2, ...] }
    model_scores = {}

    # 1. Iterate through each line of the input data
    for line in input_text.strip().split('\n'):
        # Skip empty or malformed lines, or lines that are just the date
        if not line.strip() or ' ' not in line or line.count('/') == 2:
            continue

        try:
            # Split the line at the first space to separate score and name
            parts = line.strip().split(' ', 1)
            score_str, name = parts
            score = float(score_str)

            # Add the score to the list for that model name.
            model_scores.setdefault(name, []).append(score)

        except (ValueError, IndexError) as e:
            print(f"Skipping malformed line: '{line}'. Error: {e}")
            continue

    # Before sorting the models, first sort the scores within each model's list
    # in ascending order. This is crucial for the secondary sort key to work correctly.
    for name in model_scores:
        model_scores[name].sort()

    # 2. Convert dictionary to a list of (name, scores) items for sorting
    #    The sort key is a tuple:
    #    - First, sort by the number of scores in descending order (most scores first).
    #    - Second, sort by the list of scores itself in ascending order.
    sorted_items = sorted(
        model_scores.items(),
        key=lambda item: (-len(item[1]), item[1])
    )

    # 3. Prepare the output string
    today_date = datetime.date.today().strftime("%m/%d/%Y")
    output_lines = [today_date]

    # --- MODIFICATION START ---
    SUMMARY_THRESHOLD = 30 # Set the threshold for summarizing scores

    # Format each sorted item for the output file
    for name, scores in sorted_items:
        num_scores = len(scores)

        # If the number of scores for a model is greater than the threshold, summarize it.
        if num_scores > SUMMARY_THRESHOLD:
            min_score = scores[0]
            max_score = scores[-1]
            # Format with count, min, and max score. Using an f-string for clarity.
            output_line = f"{num_scores} scores (range: {min_score:.2f} to {max_score:.2f}) {name}"
        else:
            # Otherwise, list all the scores as before.
            scores_str = [str(s) for s in scores]
            output_line = ' '.join(scores_str) + ' ' + name
        
        output_lines.append(output_line)
    # --- MODIFICATION END ---

    return '\n'.join(output_lines)

def main():
    """
    Main function to read data from an input file specified on the command line,
    process it, and write the sorted results to a new file.
    """
    # 1. Check for command-line argument for the input file
    if len(sys.argv) != 2:
        # Provide usage instructions if the argument is missing
        print("Usage: python process_scores.py <input_filename>")
        sys.exit(1) # Exit the script indicating an error

    input_filename = sys.argv[1]

    # 2. Read the content from the specified input file
    try:
        with open(input_filename, 'r') as f:
            input_text = f.read()
    except FileNotFoundError:
        print(f"Error: The file '{input_filename}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        sys.exit(1)

    # Process the data read from the file
    output_content = process_and_sort_data(input_text)
    
    # --- FIX ---
    # Separate the directory path and the filename from the input path
    dir_path = os.path.dirname(input_filename)
    base_name = os.path.basename(input_filename)
    
    # Create the new filename and join it with the original directory path
    output_filename = os.path.join(dir_path, f"sorted_{base_name}")
    # --- END FIX ---
    
    # 3. Create and write the results to the new file
    try:
        with open(output_filename, "w") as f:
            f.write(output_content)
    except Exception as e:
        print(f"An error occurred while writing to the output file: {e}")
        sys.exit(1)

    print(f"Successfully processed '{input_filename}' and created '{output_filename}' with the sorted data.")


# Run the main function when the script is executed
if __name__ == "__main__":
    main()