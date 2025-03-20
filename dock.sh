# 3/18/2025
#!/bin/bash

# Directory paths
PROTEIN_DIR="./dock/proteins"
LIGAND_DIR="./dock/ligands"
CONFIG_DIR="./dock/config"
CACHE_DIR="./dock/cache"
RESULTS_DIR="./dock/results"

# Cache files
LIGAND_NAMES="$CACHE_DIR/ligandNames.txt"
PROTEIN_NAMES="$CACHE_DIR/proteinNames.txt"
PROGRESS_CACHE="$CACHE_DIR/progress_cache.txt"

# Ensure necessary directories exist
mkdir -p "$PROTEIN_DIR" "$LIGAND_DIR" "$CONFIG_DIR" "$CACHE_DIR" "$RESULTS_DIR"

# Function to initialize or clear cache
initialize_cache() {
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

# Function to perform docking
perform_docking() {
    local ligand_file="$1"
    local protein_file="$2"
    local model_index="$3"
    ligand_name=$(basename "$ligand_file" .pdbqt)  # Remove 'local' to make it global
    protein_name=$(basename "$protein_file" .pdbqt)  # Remove 'local' to make it global
    local output_file="$RESULTS_DIR/${ligand_name}_model${model_index}_vs_${protein_name}.pdbqt"
    local log_file="$RESULTS_DIR/${ligand_name}_model${model_index}_vs_${protein_name}.log"

    # Calculate docking box parameters
    read center_x center_y center_z size_x size_y size_z <<< $(calculate_docking_box "$protein_file")

    # Debugging information
    echo "Running AutoDock Vina with the following parameters:"
    echo "Receptor: $protein_file"
    echo "Ligand: $ligand_file"
    echo "Center: $center_x, $center_y, $center_z"
    echo "Size: $size_x, $size_y, $size_z"
    echo "Output: $output_file"
    echo "Log: $log_file"

    # Run AutoDock Vina
    if ! vina --receptor "$protein_file" \
              --ligand "$ligand_file" \
              --center_x "$center_x" \
              --center_y "$center_y" \
              --center_z "$center_z" \
              --size_x "$size_x" \
              --size_y "$size_y" \
              --size_z "$size_z" \
              --out "$output_file" &>> "$log_file"; then
        echo "Error: Docking failed for $ligand_file with $protein_file. Check log file: $log_file"
        exit 1
    fi
}

# Function to display progress
display_progress() {
    local total_ligands total_proteins total_models
    total_ligands=$(wc -l < "$LIGAND_NAMES")
    total_proteins=$(wc -l < "$PROTEIN_NAMES")
    total_models=${#models[@]}

    clear
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
    for ((i = 0; i < $(((MODEL_INDEX + 1) * 20 / total_models)); i++)); do echo -n -e "\e[42m \e[0m"; done
    for ((i = $(((MODEL_INDEX + 1) * 20 / total_models)); i < 20; i++)); do echo -n -e "\e[41m \e[0m"; done
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

total_ligands=${#ligands[@]}
total_proteins=${#proteins[@]}

# Display initial progress
display_progress

while (( LIGAND_INDEX < total_ligands )); do
    ligand_file="${ligands[$LIGAND_INDEX]}"
    extract_models "$ligand_file"
    models=("$CACHE_DIR/models_$(basename "$ligand_file" .pdbqt)"/*.pdbqt)
    total_models=${#models[@]}

    # Calculate total docking tasks
    total_tasks=$((total_ligands * total_models * total_proteins))
    current_task=$((LIGAND_INDEX * total_models * total_proteins + MODEL_INDEX * total_proteins + PROTEIN_INDEX + 1))

    # Display progress before starting docking
    display_progress

    while (( MODEL_INDEX < total_models )); do
        model_file="${models[$MODEL_INDEX]}"

        while (( PROTEIN_INDEX < total_proteins )); do
            protein_file="${proteins[$PROTEIN_INDEX]}"

            # Perform docking
            perform_docking "$model_file" "$protein_file" "$MODEL_INDEX"

            # Display progress
            display_progress

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

echo "Docking process completed."
