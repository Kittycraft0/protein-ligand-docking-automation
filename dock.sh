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

    # Initialize progress cache
    echo "LIGAND_INDEX=0" > "$PROGRESS_CACHE"
    echo "MODEL_INDEX=0" >> "$PROGRESS_CACHE"
    echo "PROTEIN_INDEX=0" >> "$PROGRESS_CACHE"
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
    local ligand_name=$(basename "$ligand_file" .pdbqt)
    local model_index=0
    local output_dir="$CACHE_DIR/models_$ligand_name"
    mkdir -p "$output_dir"

    awk -v output_dir="$output_dir" -v ligand_name="$ligand_name" '
    /^MODEL/ {
        model_index = $2
        output_file = sprintf("%s/%s_model_%d.pdbqt", output_dir, ligand_name, model_index)
        print > output_file
        next
    }
    /^ENDMDL/ {
        print >> output_file
        close(output_file)
        next
    }
    {
        print >> output_file
    }
    ' "$ligand_file"
}

# Function to perform docking
perform_docking() {
    local ligand_file="$1"
    local protein_file="$2"
    local model_index="$3"
    local ligand_name
    local protein_name
    ligand_name=$(basename "$ligand_file" .pdbqt)
    protein_name=$(basename "$protein_file" .pdbqt)
    local output_file="$RESULTS_DIR/${ligand_name}_model${model_index}_vs_${protein_name}.pdbqt"
    local log_file="$RESULTS_DIR/${ligand_name}_model${model_index}_vs_${protein_name}.log"

    # Calculate docking box parameters
    read center_x center_y center_z size_x size_y size_z <<< $(calculate_docking_box "$protein_file")

    # Run AutoDock Vina
    if ! vina --receptor "$protein_file" \
              --ligand "$ligand_file" \
              --center_x "$center_x" \
              --center_y "$center_y" \
              --center_z "$center_z" \
              --size_x "$size_x" \
              --size_y "$size_y" \
              --size_z "$size_z" \
              --out "$output_file" \
              --log "$log_file" &>> "$log_file"; then
        echo "Error: Docking failed for $ligand_file with $protein_file"
        exit 1
    fi
}

# Function to display progress
display_progress() {
    local total_ligands total_proteins
    total_ligands=$(wc -l < "$LIGAND_NAMES")
    total_proteins=$(wc -l < "$PROTEIN_NAMES")

    echo "Docking Progress:"
    echo "Total ligands: $total_ligands"
    echo "Total proteins: $total_proteins"
    echo "Current ligand index: $LIGAND_INDEX"
    echo "Current model index: $MODEL_INDEX"
    echo "Current protein index: $PROTEIN_INDEX"
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
else
    initialize_cache
fi

read_progress

ligands=($(cat "$LIGAND_NAMES"))
proteins=($(cat "$PROTEIN_NAMES"))

total_ligands=${#ligands[@]}
total_proteins=${#proteins[@]}

# Calculate total docking tasks
total_tasks=$((total_ligands * total_proteins))
current_task=1

while (( LIGAND_INDEX < total_ligands )); do
    ligand_file="${ligands[$LIGAND_INDEX]}"
    extract_models "$ligand_file"
    models=("$CACHE_DIR/models_$(basename "$ligand_file" .pdbqt)"/*.pdbqt)
    total_models=${#models[@]}

    while (( MODEL_INDEX < total_models )); do
        model_file="${models[$MODEL_INDEX]}"

        while (( PROTEIN_INDEX < total_proteins )); do
            protein_file="${proteins[$PROTEIN_INDEX]}"

            # Perform docking
            perform_docking "$model_file" "$protein_file" "$MODEL_INDEX"

            # Display progress
            echo "Task $current_task of $total_tasks completed."
            ((current_task++))

            # Update progress
            ((PROTEIN_INDEX++))
            update_progress
            display_progress
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
