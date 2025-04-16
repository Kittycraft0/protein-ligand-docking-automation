# Protein-Ligand Docking Automation

This project automates protein-ligand docking using [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina), organizing results and generating comprehensive score and ranking files for analysis. Progress is saved after each docking operation, allowing seamless resumption after interruptions.

---

## Features

- **Automated Docking:** Batch docking of multiple ligands to multiple proteins.
- **Progress Persistence:** Resume from last successful step after interruptions.
- **Dynamic Configuration:** Customize docking parameters via a config file.
- **Comprehensive Logging:** Detailed logs for each docking operation.
- **Structured Results:** Output files and folders organized for easy analysis.

---

## Prerequisites

- **AutoDock Vina:** [Download and install](https://github.com/ccsb-scripps/AutoDock-Vina).
- **Python 3.7+**
- **Windows or Unix-like environment**

---

## Setup

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/yourusername/protein-ligand-docking-automation.git
   cd protein-ligand-docking-automation
   ```

2. **Prepare Input Files:**
   - Place `.pdbqt` protein files in `dock/proteins/`
   - Place `.pdbqt` ligand files in `dock/ligands/`
   - Place `.pdbqt` comparison ligand files in `dock/comparison_ligands/`

3. **Configure Docking Parameters:**
   - Edit `dock/config/config.txt` to set Vina parameters (see below).

---

## Folder Structure

```
dock/
  config/
    config.txt           # Docking parameters
  proteins/              # Protein .pdbqt files
  ligands/               # Ligand .pdbqt files
  comparison_ligands/    # Comparison ligand .pdbqt files
  cache/                 # Internal cache files
  results/
    scores/
      scores_<protein_name>.txt                # Scores for each protein
      scores_<comparisonLigand_name>_RMS.txt   # RMS scores for each comparison ligand
      scores_<comparisonLigand_name>/
        scores_<comparisonLigand_name>_in_<protein_name>.txt
      scores_proteins/
        scores_<protein_name>/
          scores_<protein_name>_unsorted.txt
          scores_<protein_name>_sorted.txt
    docked_ligands/
      docked_<ligand_name>/
        <ligand_name>_scores.txt
        <ligand_name>_in_<protein_name>.log
        <ligand_name>_in_<protein_name>.pdbqt
      docked_<ligand_name>_modelXX/
        ...
  log/
    docking_log.txt      # Log of all docking runs
```

---

## Configuration

Edit `dock/config/config.txt` to set Vina parameters:

- `cpu`: Number of CPU threads to use (e.g., 8)
- `exhaustiveness`: Search exhaustiveness (higher = more accurate, slower)
- `energy_range`: Energy range for output modes
- `num_modes`: Number of binding modes to generate

Example:
```
cpu=8
exhaustiveness=16
energy_range=4
num_modes=9
```

---

## Output Files

- **scores_<protein_name>.txt**: Scores for all ligands docked to a protein.
- **scores_<comparisonLigand_name>_RMS.txt**: RMS of score differences for each ligand vs. a comparison ligand.
- **scores_<comparisonLigand_name>/scores_<comparisonLigand_name>_in_<protein_name>.txt**: Score differences for each ligand vs. a comparison ligand for a specific protein.
- **scores_proteins/scores_<protein_name>_unsorted.txt**: Unsorted scores for a protein.
- **scores_proteins/scores_<protein_name>_sorted.txt**: Sorted scores for a protein.
- **docked_ligands/docked_<ligand_name>/**: Contains all docking result files and a summary for each ligand.

---

## Log Files

- **log/docking_log.txt**: Chronological log of all docking runs, errors, and events.

---

## Additional Data & Recommendations

- **Docking parameters** for each run are stored in logs or as metadata.
- **Timestamps** for each docking operation are logged.
- **Failed docking attempts** are recorded in `docking_log.txt`.
- **Summary statistics** (mean, median, stddev) for each ligand/protein/comparison ligand can be generated from scores.
- **README** (this file) describes the structure and usage.
- **Consider adding:**  
  - Per-ligand/protein statistics files  
  - A summary file listing all ligands, proteins, and comparison ligands used  
  - A metadata file for each docking run with parameters and environment info

---

## Usage

1. Place your protein, ligand, and comparison ligand `.pdbqt` files in their respective folders.
2. Edit `config/config.txt` as needed.
3. Run the script:
   ```bash
   python dock.py
   ```
4. Results will appear in the `results/` folder.

---

## Contact

For questions or issues, contact [iwbmo05@gmail.com](mailto:iwbmo05@gmail.com) or iwbmo on Discord.