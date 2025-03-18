# Protein-Ligand Docking Automation

This project automates the process of protein-ligand docking using [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina). It facilitates batch processing of multiple protein and ligand files, ensuring progress is saved after each docking operation and allowing the process to resume seamlessly after interruptions.

I wrote up what I wanted the program to do, but I mostly had AI write everything. Including this entire readme, except for this line.

## Features

- **Automated Docking:** Sequentially docks multiple ligands to multiple proteins without manual intervention.
- **Progress Persistence:** Saves progress after each docking operation, enabling resumption from the last successful step after interruptions.
- **Dynamic Configuration:** Allows customization of docking parameters through a configuration file.
- **Comprehensive Logging:** Provides detailed logs for each docking operation, facilitating easy monitoring and debugging.

## Prerequisites

- **AutoDock Vina:** Ensure that [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) is installed and accessible in your system's PATH.
- **Bash Shell:** The automation script is written in Bash; a Unix-like environment is required.

## Setup

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/protein-ligand-docking-automation.git
   cd protein-ligand-docking-automation
Prepare Protein Files:

Download protein structures and prepare them using ChimeraX.
Identify and log docking zones for each protein.
Convert protein files to .pdbqt format and place them in the proteins directory.
Prepare Ligand Files:

Download ligand libraries from the ZINC Database.
Extract and convert ligand files to .pdbqt format, placing them in the ligands directory.
Configure Docking Parameters:

Edit the config/vina_config.txt file to set desired docking parameters.
Usage
Run the Automation Script:

bash
Copy
Edit
./dock_automation.sh
Monitor Progress:

The script displays real-time progress, including:

Total progress across all ligands and proteins.
Progress of the current ligand being processed.
Progress of the current docking operation between a ligand and a protein.
Handling Interruptions:

In case of interruptions (e.g., system shutdowns), simply rerun the script. It will resume from the last saved progress point.

Results
Docking results are saved in the results directory, organized by protein and ligand names. Each result includes:

Docked Structures: Output files in .pdbqt format.
Logs: Detailed logs of each docking operation for review.
License
This project is licensed under the MIT License.

Acknowledgements
AutoDock Vina for the docking engine.
ChimeraX for protein preparation.
ZINC Database for ligand libraries.
pgsql
Copy
Edit


This `README.md` provides a comprehensive overview of the project, guiding users through setup, usage, and understanding the results. It ensures that users with the necessary prerequisites can effectively utilize the automation script for protein-ligand docking tasks.
::contentReference[oaicite:0]{index=0}
 
