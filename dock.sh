# 3/18/2025
#!/bin/bash

# Directory paths
PROTEIN_DIR="./dock/proteins"
LIGAND_DIR="./dock/ligands"
CONFIG_DIR="./dock/config"
CACHE_DIR="./dock/cache"
RESULTS_DIR="./dock/results"
COMPARISON_LIGAND_DIR="./dock/comparison_ligands"

# Cache files
LIGAND_NAMES="$CACHE_DIR/ligandNames.txt"
PROTEIN_NAMES="$CACHE_DIR/proteinNames.txt"
PROGRESS_CACHE="$CACHE_DIR/progress_cache.txt"

# Ensure necessary directories exist
mkdir -p "$PROTEIN_DIR" "$LIGAND_DIR" "$CONFIG_DIR" "$CACHE_DIR" "$RESULTS_DIR" "$COMPARISON_LIGAND_DIR"

# Initialize DEBUG variable
DEBUG=0

# Check for --debug flag
for arg in "$@"; do
    if [[ "$arg" == "--debug" ]]; then
        DEBUG=1
    fi
done

# Function to initialize or clear cache
initialize_cache() {
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "Debug mode enabled. Skipping cache and results clearing."
        return
    fi

    if [[ "$1" == "clear" ]]; then
        mkdir -p "$CACHE_DIR/cache_backup"
        mv "$CACHE_DIR"/* "$CACHE_DIR/cache_backup/" 2>/dev/null
        echo "Cache cleared and backed up."
    elif [[ "$1" == "clear-everything" ]]; then
        rm -rf "$CACHE_DIR"
        rm -rf "$RESULTS_DIR"
        mkdir -p "$CACHE_DIR"
        mkdir -p "$RESULTS_DIR"
        echo "Cache and results cleared."
    fi

    mkdir -p "$CACHE_DIR"
    mkdir -p "$RESULTS_DIR"

    # Generate ligand and protein name lists
    if ! ls "$LIGAND_DIR"/*.pdbqt > "$LIGAND_NAMES" 2>/dev/null; then
        echo "Error: No ligand files found in $LIGAND_DIR"
        exit 1
    fi

    if ! ls "$PROTEIN_DIR"/*.pdbqt > "$PROTEIN_NAMES" 2>/dev/null; then
        echo "Error: No protein files found in $PROTEIN_DIR"
        exit 1
    fi

    if ! ls "$COMPARISON_LIGAND_DIR"/*.pdbqt > "$CACHE_DIR/comparisonLigandNames.txt" 2>/dev/null; then
        echo "Error: No comparison ligand files found in $COMPARISON_LIGAND_DIR"
        exit 1
    fi

    # Initialize progress cache if it doesn't exist
    if [[ ! -f "$PROGRESS_CACHE" ]]; then
        echo "LIGAND_INDEX=0" > "$PROGRESS_CACHE"
        echo "MODEL_INDEX=0" >> "$PROGRESS_CACHE"
        echo "PROTEIN_INDEX=0" >> "$PROGRESS_CACHE"
    fi
}

# Function to read progress from cache
read_progress() {
    if [[ -f "$PROGRESS_CACHE" ]]; then
        source "$PROGRESS_CACHE"
    else
        echo "Error: Progress cache file not found."
        exit 1
    fi
}

# Function to update progress in cache
update_progress() {
    echo "LIGAND_INDEX=$LIGAND_INDEX" > "$PROGRESS_CACHE"
    echo "MODEL_INDEX=$MODEL_INDEX" >> "$PROGRESS_CACHE"
    echo "PROTEIN_INDEX=$PROTEIN_INDEX" >> "$PROGRESS_CACHE"
}

# Function to calculate docking box parameters
calculate_docking_box() {
    local protein_file="$1"
    local center_x center_y center_z size_x size_y size_z

    # Calculate the geometric center and size of the docking box
    center_x=$(awk '/^ATOM/ {x+=$7; count++} END {print x/count}' "$protein_file")
    center_y=$(awk '/^ATOM/ {y+=$8; count++} END {print y/count}' "$protein_file")
    center_z=$(awk '/^ATOM/ {z+=$9; count++} END {print z/count}' "$protein_file")
    size_x=20
    size_y=20
    size_z=20

    echo "$center_x $center_y $center_z $size_x $size_y $size_z"
}

# Function to extract models from ligand files
extract_models() {
    local ligand_file="$1"
    ligand_name=$(basename "$ligand_file" .pdbqt)  # Remove 'local' to make it global
    local output_dir="$CACHE_DIR/models_$ligand_name"
    mkdir -p "$output_dir"

    # Check if the models are already extracted
    if ls "$output_dir/${ligand_name}_model"*.pdbqt 1> /dev/null 2>&1; then
        echo "Models for $ligand_name already extracted."
        return
    fi

    # Check if the ligand file contains multiple models
    if grep -q "MODEL" "$ligand_file"; then
        # Use vina_split to split the ligand file into multiple PDBQT files
        if ! vina_split --input "$ligand_file" --ligand "$output_dir/${ligand_name}_model"; then
            echo "Error: Failed to split ligand file $ligand_file"
            exit 1
        fi

        # Check if the output directory contains any split files
        if [ ! "$(ls -A $output_dir)" ]; then
            echo "Error: No models generated for ligand file $ligand_file"
            exit 1
        fi

        # Rename the split files to follow the expected naming convention
        for model_file in "$output_dir/${ligand_name}_model"*.pdbqt; do
            model_index=$(basename "$model_file" | sed -E 's/.*_model([0-9]+)\.pdbqt/\1/')
            mv "$model_file" "$output_dir/${ligand_name}_model_${model_index}.pdbqt"
        done
    else
        # Copy the single model file to the output directory
        cp "$ligand_file" "$output_dir/${ligand_name}_model_1.pdbqt"
    fi
}

# Function to calculate RMS relative to ABA number
calculate_rms_relative_to_aba() {
    local aba_score="$1"
    local score="$2"
    local rms_relative_to_aba=$(echo "scale=8; sqrt(($score - $aba_score)^2)" | bc -l)
    echo "$rms_relative_to_aba"
}

# Function to monitor the progress of the current docking process
monitor_docking_progress() {
    local log_file="$1"
    local progress=0

    while [[ ! -f "$log_file" ]]; do
        sleep 0.2  # Check every 0.2 seconds
    done

    while [[ $progress -lt 100 ]]; do
        if [[ -f "$log_file" ]]; then
            progress=$(grep -oP '(?<=\|)[\*]{1,}' "$log_file" | wc -c)
            progress=$((progress * 5))  # Each '*' represents 5% progress
            echo -ne "Current docking progress: $progress%\r"
        fi
        sleep 0.2  # Update display every 0.2 seconds
    done
    echo -ne "Current docking progress: 100%\n"
}

# Function to perform docking asynchronously
perform_docking() {
    local ligand_file="$1"
    local protein_file="$2"
    local model_index="$3"
    ligand_name=$(basename "$ligand_file" .pdbqt)  # Remove 'local' to make it global
    protein_name=$(basename "$protein_file" .pdbqt)  # Remove 'local' to make it global
    local output_file="$RESULTS_DIR/temp/${ligand_name}_model${model_index}_vs_${protein_name}.pdbqt"
    local log_file="$RESULTS_DIR/temp/${ligand_name}_model${model_index}_vs_${protein_name}.log"

    # Calculate docking box parameters
    read center_x center_y center_z size_x size_y size_z <<< $(calculate_docking_box "$protein_file")

    # Debugging information
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "Running AutoDock Vina with the following parameters:"
        echo "Receptor: $protein_file"
        echo "Ligand: $ligand_file"
        echo "Center: $center_x, $center_y, $center_z"
        echo "Size: $size_x, $size_y, $size_z"
        echo "Output: $output_file"
        echo "Log: $log_file"
    fi

    # Run AutoDock Vina asynchronously
    vina --receptor "$protein_file" \
         --ligand "$ligand_file" \
         --center_x "$center_x" \
         --center_y "$center_y" \
         --center_z "$center_z" \
         --size_x "$size_x" \
         --size_y "$size_y" \
         --size_z "$size_z" \
         --out "$output_file" &>> "$log_file" &

    # Get the PID of the background process
    local vina_pid=$!

    # Monitor the docking progress
    monitor_docking_progress "$log_file" &

    # Wait for the docking process to complete
    wait $vina_pid

    # Check if the log file was created
    if [[ ! -f "$log_file" ]]; then
        echo "Error: Log file not found: $log_file"
        exit 1
    fi

    # Display the contents of the log file for debugging
    if [[ "$DEBUG" -eq 1 ]]; then
        echo "Contents of log file $log_file:"
        cat "$log_file"
    fi

    # Extract score from log file and store it
    local score=$(grep -m 1 "^   1" "$log_file" | awk '{print $2}')
    if [[ -n "$score" && "$score" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
        if [[ "$DEBUG" -eq 1 ]]; then
            echo "Extracted score: $score for ligand: $ligand_name with protein: $protein_name"  # Debugging information
        fi
        echo "$score $ligand_name" >> "$RESULTS_DIR/scores_${protein_name}.txt"
        if [[ "$DEBUG" -eq 1 ]]; then
            echo "Score stored in $RESULTS_DIR/scores_${protein_name}.txt"  # Debugging information
        fi
    else
        echo "Failed to extract a valid score from log file: $log_file"  # Debugging information
    fi
}

# Function to display progress
display_progress() {
    local total_ligands total_proteins total_models
    total_ligands=$(wc -l < "$LIGAND_NAMES")
    total_proteins=$(wc -l < "$PROTEIN_NAMES")
    total_models=${#models[@]}

    # clear
    echo "========================================"
    echo "         Protein-Ligand Docking         "
    echo "========================================"
    echo
    echo "Docking:"
    echo "Total progress: ligand $((LIGAND_INDEX + 1))/$total_ligands"
    echo -n "["
    for ((i = 0; i < $(((LIGAND_INDEX + 1) * 20 / total_ligands)); i++)); do echo -n -e "\e[42m \e[0m"; done
    for ((i = $(((LIGAND_INDEX + 1) * 20 / total_ligands)); i < 20; i++)); do echo -n -e "\e[41m \e[0m"; done
    echo "]"
    echo

    echo "Ligand model: model $((MODEL_INDEX + 1))/$total_models"
    echo -n "["
    if (( total_models > 0 )); then
        for ((i = 0; i < $(((MODEL_INDEX + 1) * 20 / total_models)); i++)); do echo -n -e "\e[42m \e[0m"; done
        for ((i = $(((MODEL_INDEX + 1) * 20 / total_models)); i < 20; i++)); do echo -n -e "\e[41m \e[0m"; done
    else
        for ((i = 0; i < 20; i++)); do echo -n -e "\e[41m \e[0m"; done
    fi
    echo "]"
    echo

    echo "Ligand progress: protein $((PROTEIN_INDEX + 1))/$total_proteins"
    echo -n "["
    for ((i = 0; i < $(((PROTEIN_INDEX + 1) * 20 / total_proteins)); i++)); do echo -n -e "\e[42m \e[0m"; done
    for ((i = $(((PROTEIN_INDEX + 1) * 20 / total_proteins)); i < 20; i++)); do echo -n -e "\e[41m \e[0m"; done
    echo "]"
    echo

    echo "Ligand name: $ligand_name"
    echo "Protein name: $protein_name"
    local percent_through_dock=$(printf "%.8f" $(echo "scale=8; $current_task * 100 / $total_tasks" | bc))
    echo "% through dock: $percent_through_dock%"
    echo -n "["
    for ((i = 0; i < $((current_task * 40 / total_tasks)); i++)); do echo -n -e "\e[42m \e[0m"; done
    for ((i = $((current_task * 40 / total_tasks)); i < 40; i++)); do echo -n -e "\e[41m \e[0m"; done
    echo "]"
    echo
    echo "Do CTRL+C to exit"
    echo
}

# Function to display progress for comparison ligands
display_comparison_progress() {
    local total_comparison_ligands total_proteins
    total_comparison_ligands=$(wc -l < "$CACHE_DIR/comparisonLigandNames.txt")
    total_proteins=$(wc -l < "$PROTEIN_NAMES")

    echo "========================================"
    echo "   Initial Calculation of Comparison    "
    echo "========================================"
    echo
    echo "Docking:"
    echo "Total progress: comparison ligand $((COMPARISON_LIGAND_INDEX + 1))/$total_comparison_ligands"
    echo -n "["
    for ((i = 0; i < $(((COMPARISON_LIGAND_INDEX + 1) * 20 / total_comparison_ligands)); i++)); do echo -n -e "\e[42m \e[0m"; done
    for ((i = $(((COMPARISON_LIGAND_INDEX + 1) * 20 / total_comparison_ligands)); i < 20; i++)); do echo -n -e "\e[41m \e[0m"; done
    echo "]"
    echo

    echo "Comparison ligand progress: protein $((PROTEIN_INDEX + 1))/$total_proteins"
    echo -n "["
    for ((i = 0; i < $(((PROTEIN_INDEX + 1) * 20 / total_proteins)); i++)); do echo -n -e "\e[42m \e[0m"; done
    for ((i = $(((PROTEIN_INDEX + 1) * 20 / total_proteins)); i < 20; i++)); do echo -n -e "\e[41m \e[0m"; done
    echo "]"
    echo

    echo "Comparison ligand name: $comparison_ligand_name"
    echo "Protein name: $protein_name"
    local percent_through_dock=$(printf "%.8f" $(echo "scale=8; $current_task * 100 / $total_tasks" | bc))
    echo "% through dock: $percent_through_dock%"
    echo -n "["
    for ((i = 0; i < $((current_task * 40 / total_tasks)); i++)); do echo -n -e "\e[42m \e[0m"; done
    for ((i = $((current_task * 40 / total_tasks)); i < 40; i++)); do echo -n -e "\e[41m \e[0m"; done
    echo "]"
    echo
    echo "Do CTRL+C to exit"
    echo
}

# Function to handle termination signals
terminate_script() {
    echo "Termination signal received. Saving progress..."
    update_progress
    exit 0
}

# Trap SIGINT (Ctrl+C) and SIGTERM to execute terminate_script
trap terminate_script SIGINT SIGTERM

# Main script execution
if [[ "$1" == "--clear-cache" ]]; then
    initialize_cache "clear"
elif [[ "$1" == "--clear-everything" ]]; then
    initialize_cache "clear-everything"
else
    initialize_cache
fi

read_progress

ligands=($(cat "$LIGAND_NAMES"))
proteins=($(cat "$PROTEIN_NAMES"))
comparison_ligands=($(cat "$CACHE_DIR/comparisonLigandNames.txt"))

total_ligands=${#ligands[@]}
total_proteins=${#proteins[@]}
total_comparison_ligands=${#comparison_ligands[@]}

# Create a dedicated subfolder for output files
mkdir -p "$RESULTS_DIR/temp"

# Display initial progress for comparison ligands
COMPARISON_LIGAND_INDEX=0
for comparison_ligand_file in "${comparison_ligands[@]}"; do
    comparison_ligand_name=$(basename "$comparison_ligand_file" .pdbqt)
    for protein_file in "${proteins[@]}"; do
        protein_name=$(basename "$protein_file" .pdbqt)
        display_comparison_progress
        perform_docking "$comparison_ligand_file" "$protein_file" 0
        ((PROTEIN_INDEX++))
    done
    PROTEIN_INDEX=0
    ((COMPARISON_LIGAND_INDEX++))
done

# Perform docking for ligands
while (( LIGAND_INDEX < total_ligands )); do
    ligand_file="${ligands[$LIGAND_INDEX]}"
    extract_models "$ligand_file"
    models=("$CACHE_DIR/models_$(basename "$ligand_file" .pdbqt)"/*.pdbqt)
    total_models=${#models[@]}

    # Calculate total docking tasks
    total_tasks=$((total_ligands * total_models * total_proteins))
    current_task=$((LIGAND_INDEX * total_models * total_proteins + MODEL_INDEX * total_proteins + PROTEIN_INDEX + 1))

    while (( MODEL_INDEX < total_models )); do
        model_file="${models[$MODEL_INDEX]}"

        while (( PROTEIN_INDEX < total_proteins )); do
            protein_file="${proteins[$PROTEIN_INDEX]}"

            # Update ligand and protein names for progress display
            ligand_name=$(basename "$model_file" .pdbqt)
            protein_name=$(basename "$protein_file" .pdbqt)

            # Display progress before starting docking
            display_progress

            # Perform docking
            perform_docking "$model_file" "$protein_file" "$MODEL_INDEX"

            # Update progress
            ((current_task++))
            ((PROTEIN_INDEX++))
            update_progress
        done

        PROTEIN_INDEX=0
        ((MODEL_INDEX++))
        update_progress
    done

    MODEL_INDEX=0
    ((LIGAND_INDEX++))
    update_progress
done

# Move all .log and .pdbqt files to the dedicated subfolder
mv "$RESULTS_DIR/temp"/* "$RESULTS_DIR/"

# Calculate and display top dockers for each protein
for protein_file in "${proteins[@]}"; do
    protein_name=$(basename "$protein_file" .pdbqt)
    if [[ -f "$RESULTS_DIR/scores_${protein_name}.txt" ]]; then
        sort -n "$RESULTS_DIR/scores_${protein_name}.txt" > "$RESULTS_DIR/top_dockers_${protein_name}.txt"
    else
        echo "No scores file found for protein: $protein_name. Skipping sorting."
    fi
done

# Calculate the best ligands using RMS relative to comparison ligands
echo "Calculating best ligands using RMS relative to comparison ligands..."
for ligand_file in "${ligands[@]}"; do
    ligand_name=$(basename "$ligand_file" .pdbqt)
    total_score=0
    count=0
    for protein_file in "${proteins[@]}"; do
        protein_name=$(basename "$protein_file" .pdbqt)
        score=$(grep "$ligand_name" "$RESULTS_DIR/scores_${protein_name}.txt" | awk '{print $1}')
        if [[ -n "$score" && "$score" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
            echo "Score for $ligand_name with $protein_name: $score"  # Debugging information
            comparison_scores=()
            for comparison_ligand_file in "${comparison_ligands[@]}"; do
                comparison_ligand_name=$(basename "$comparison_ligand_file" .pdbqt)
                comparison_score=$(grep "$comparison_ligand_name" "$RESULTS_DIR/scores_${protein_name}.txt" | awk '{print $1}')
                if [[ -n "$comparison_score" && "$comparison_score" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
                    comparison_scores+=("$comparison_score")
                fi
            done
            rms_relative_to_comparison=$(calculate_rms_relative_to_comparison "$score" "${comparison_scores[@]}")
            echo "RMS relative to comparison for $ligand_name with $protein_name: $rms_relative_to_comparison"  # Debugging information
            total_score=$(echo "$total_score + 1/($rms_relative_to_comparison^2)" | bc -l 2>/dev/null)
            ((count++))
        fi
    done
    if (( count > 0 )); then
        final_score=$(echo "scale=8; sqrt($count / $total_score)" | bc -l 2>/dev/null)
        if [[ -n "$final_score" && "$final_score" != "0" ]]; then
            echo "Final score for $ligand_name: $final_score"  # Debugging information
            echo "$final_score $ligand_name" >> "$RESULTS_DIR/best_ligands.txt"
        fi
    fi
done

# Sort the best ligands file
if [[ -f "$RESULTS_DIR/best_ligands.txt" ]]; then
    sort -n "$RESULTS_DIR/best_ligands.txt" > "$RESULTS_DIR/ranked_best_ligands.txt"
    echo "Best ligands ranked and saved to $RESULTS_DIR/ranked_best_ligands.txt"
else
    echo "No best ligands file found. Skipping ranking."
fi

# Display the top 5 ligands and their overall scores
if [[ -f "$RESULTS_DIR/ranked_best_ligands.txt" ]]; then
    echo "Top 5 ligands:"
    head -n 5 "$RESULTS_DIR/ranked_best_ligands.txt"
else
    echo "No ranked best ligands file found. Cannot display top 5 ligands."
fi

echo "Docking process completed. Results saved in $RESULTS_DIR."
