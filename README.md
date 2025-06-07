# Protein-Ligand Docking Automation Pipeline

A robust Python pipeline to automate high-throughput molecular docking using AutoDock Vina. This script manages file organization, executes docking runs, saves progress, and generates comprehensive analysis files to rank ligands based on their performance against a set of proteins and user-defined comparison ligands.

---

## Key Features

-   **High-Throughput Docking:** Automatically runs a matrix of docking jobs for multiple ligands against multiple proteins.
-   **Persistent & Resumable Workflow:** Saves progress after every single docking operation. If the script is stopped (even with `Ctrl+C`), it can be restarted and will resume exactly where it left off.
-   **Comparative Analysis Engine:** Goes beyond simple docking scores by calculating rankings based on score differences and Root Mean Square (RMS) values relative to known "comparison ligands."
-   **Automated File Management:** Intelligently organizes all inputs, outputs, logs, and results into a clean and predictable folder structure. Handles multi-model ligands by splitting them into individual files.
-   **Dynamic Configuration:** All Vina parameters (`cpu`, `exhaustiveness`, etc.) and protein-specific docking box coordinates are managed through simple text configuration files. The script will auto-generate templates for you.
-   **Live Progress Monitoring:** Detailed, multi-bar progress display in the terminal shows overall progress, per-file progress, and an estimated time to completion (ETC).

---

## Prerequisites

-   **Python 3.7+**: The script is written in Python.
-   **AutoDock Vina**: The core docking engine. Both `vina.exe` and `vina_split.exe` must be installed and accessible via your system's PATH.
    -   [Download and Install AutoDock Vina](https://vina.scripps.edu/downloads/)
    -   **Important:** Ensure the Vina installation directory is added to your system's `PATH` environment variable to avoid `FileNotFoundError`.

---

## Installation & Setup

1.  **Clone the Repository (or place the script):**
    ```bash
    git clone [https://github.com/yourusername/protein-ligand-docking-automation.git](https://github.com/yourusername/protein-ligand-docking-automation.git)
    cd protein-ligand-docking-automation
    ```

2.  **Prepare Input Files (Crucial Step):**
    This script expects that your protein and ligand files are already in the correct **`.pdbqt`** format. AutoDock Vina requires files with proper partial charges (like Gasteiger charges) and specific atom types.
    -   Use tools like **AutoDock Tools (ADT)** or **Open Babel** to convert your source files (e.g., `.pdb`) into properly prepared `.pdbqt` files before placing them in the folders below.

3.  **Set Up Folder Structure:**
    The script operates within a main `dock/` directory. Place your prepared `.pdbqt` files into the appropriate subdirectories:
    -   `dock/proteins/`
    -   `dock/ligands/`
    -   `dock/comparison_ligands/`

4.  **Initial Run & Configuration:**
    -   The **first time** you run the script, it will see that the configuration files are missing and will automatically create templates for you.
    -   It will then stop and instruct you to fill them out.
    -   **`dock/config/config.txt`**: Open this file and set your desired Vina parameters (cpu, exhaustiveness, etc.).
    -   **`dock/config/config_<protein_name>.txt`**: For each protein, open its corresponding config file and paste in the `center_x/y/z` and `size_x/y/z` values for the docking box.

---

## Folder Structure Overview

```
dock/
├── cache/
│   ├── models_<ligand_name>/ # Sub-folder for split ligand models
│   │   └── *.pdbqt
│   ├── ligandNames.txt
│   ├── proteinNames.txt
│   └── progress_cache.txt
│
├── config/
│   ├── config.txt                 # Main Vina parameters
│   └── config_<protein_name>.txt  # Per-protein box coordinates
│
├── comparison_ligands/
│   └── *.pdbqt                    # Your reference ligands
│
├── ligands/
│   └── *.pdbqt                    # The ligands you want to test
│
├── log/
│   └── docking_log.txt            # Master log of all script events
│
├── proteins/
│   └── *.pdbqt                    # Your receptor proteins
│
└── results/
    ├── docked_ligands/
    │   └── docked_<ligand_name>/  # Sub-folder for each ligand's results
    │       ├── <ligand>_in_<protein>.log
    │       ├── <ligand>_in_<protein>.pdbqt
    │       └── <ligand>_scores.txt
    │
    └── scores/
        ├── scores_<protein_name>.txt          # Master sorted scores for a protein
        ├── scores_<comparison>_RMS.txt      # Overall ranking vs. a comparison ligand
        └── scores_<comparison>/             # Per-protein rankings vs. a comparison ligand
            └── scores_<comp>_in_<prot>.txt
```

---

## Usage

Navigate to the script's directory in your terminal and run it using Python.

#### Basic Run
This will start or resume the docking process using the existing configuration.
```bash
python dock.py
```

#### Command-Line Arguments
Control the script's behavior with optional flags:

-   **`--clear-cache`**: (Safe Reset) Archives the contents of `dock/cache/` to a backup folder. This does not delete your results. Use this if you have updated your input ligand files.
    ```bash
    python dock.py --clear-cache
    ```
-   **`--clear-everything`**: (Destructive Reset) **Permanently deletes both the `dock/cache/` and `dock/results/` folders.** Use this with caution when starting a completely new project.
    ```bash
    python dock.py --clear-everything
    ```
-   **`--debug`**: Runs in debug mode, printing extra information about Vina commands and log file contents.
    ```bash
    python dock.py --debug
    ```

---

## Troubleshooting

-   **Error: `FileNotFoundError: [WinError 2] The system cannot find the file specified` (referring to `vina` or `vina_split`)**
    -   **Cause:** The folder where you installed AutoDock Vina is not in your system's PATH.
    -   **Solution:** Find your Vina installation directory (e.g., `C:\Program Files (x86)\The Scripps Research Institute\Vina`) and add this full path to your Windows/Linux/macOS `Path` environment variable. You must restart your terminal after making this change.

-   **Error: `Parse error on line X in file ...`**
    -   **Cause:** One of your input `.pdbqt` files (either a ligand or a protein) is corrupted or was not prepared correctly.
    -   **Solution:** Delete the faulty `.pdbqt` file and re-create it using a reliable tool like Open Babel or AutoDock Tools, ensuring proper charges and atom types are assigned.

---

## Contact

For questions or issues, feel free to reach out.
-   **Email:** [iwbmo05@gmail.com](mailto:iwbmo05@gmail.com)
-   **Discord:** iwbmo