import os
import shutil
import subprocess
import math
import time
import signal
from pathlib import Path
from time import sleep
from collections import deque
import statistics

# Directory paths
BASE_DIR = Path("./dock")
PROTEIN_DIR = BASE_DIR / "proteins"
LIGAND_DIR = BASE_DIR / "ligands"
CONFIG_DIR = BASE_DIR / "config"
CACHE_DIR = BASE_DIR / "cache"
RESULTS_DIR = BASE_DIR / "results"
COMPARISON_LIGAND_DIR = BASE_DIR / "comparison_ligands"

# Cache files
LIGAND_NAMES = CACHE_DIR / "ligandNames.txt"
PROTEIN_NAMES = CACHE_DIR / "proteinNames.txt"
PROGRESS_CACHE = CACHE_DIR / "progress_cache.txt"
COMPARISON_LIGAND_NAMES = CACHE_DIR / "comparisonLigandNames.txt"

# Ensure necessary directories exist
for directory in [PROTEIN_DIR, LIGAND_DIR, CONFIG_DIR, CACHE_DIR, RESULTS_DIR, COMPARISON_LIGAND_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Initialize DEBUG variable
DEBUG = "--debug" in os.sys.argv

if DEBUG:
    print("Debug mode enabled.")

# Global variable to track if display_progress is rendering
IS_RENDERING = False

# Global variable to track the start time of the docking process
START_TIME = time.time()

# Global variable to track task durations
TASK_DURATIONS = deque(maxlen=30)  # Store durations for the last 30 tasks

# Global variable to track failed docking attempts
FAILED_DOCKINGS = []

CONFIG_PATH = CONFIG_DIR / "config.txt"
if not CONFIG_PATH.exists():
    with open(CONFIG_PATH, "w") as f:
        self.results_dir = results_dir
        self.ligand_dir = ligand_dir
        self.protein_dir = protein_dir
        self.comparison_ligand_dir = comparison_ligand_dir

    def initialize_cache(self, mode=None):
        if mode == "clear":
            backup_dir = self.cache_dir / "cache_backup"
            backup_dir.mkdir(parents=True, exist_ok=True)
            for item in self.cache_dir.iterdir():
                if item == backup_dir:
                    continue
                destination = backup_dir / item.name
                if destination.exists():
                    counter = 1
                    new_name = f"{item.stem}_copy{counter}{item.suffix}"
                    new_destination = backup_dir / new_name
                    while new_destination.exists():
                        counter += 1
                        new_name = f"{item.stem}_copy{counter}{item.suffix}"
                        new_destination = backup_dir / new_name
                    destination = new_destination
                shutil.move(str(item), str(destination))
            print("Cache cleared and backed up.")
        elif mode == "clear-everything":
            shutil.rmtree(self.cache_dir, ignore_errors=True)
            shutil.rmtree(self.results_dir, ignore_errors=True)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self.results_dir.mkdir(parents=True, exist_ok=True)
            print("Cache and results cleared.")

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # Generate ligand and protein name lists
        ligand_names = self.cache_dir / "ligandNames.txt"
        protein_names = self.cache_dir / "proteinNames.txt"
        comparison_ligand_names = self.cache_dir / "comparisonLigandNames.txt"
        progress_cache = self.cache_dir / "progress_cache.txt"

        if not list(self.ligand_dir.glob("*.pdbqt")):
            print(f"Error: No ligand files found in {self.ligand_dir}")
            exit(1)
        with open(ligand_names, "w") as f:
            f.writelines([str(file) + "\n" for file in self.ligand_dir.glob("*.pdbqt")])

        if not list(self.protein_dir.glob("*.pdbqt")):
            print(f"Error: No protein files found in {self.protein_dir}")
            exit(1)
        with open(protein_names, "w") as f:
            f.writelines([str(file) + "\n" for file in self.protein_dir.glob("*.pdbqt")])

        if not list(self.comparison_ligand_dir.glob("*.pdbqt")):
            print(f"Error: No comparison ligand files found in {self.comparison_ligand_dir}")
            exit(1)
        with open(comparison_ligand_names, "w") as f:
            f.writelines([str(file) + "\n" for file in self.comparison_ligand_dir.glob("*.pdbqt")])

        if not progress_cache.exists():
            with open(progress_cache, "w") as f:
                f.write("COMPARISON_LIGAND_INDEX=0\nCOMPARISON_PROTEIN_INDEX=0\nLIGAND_INDEX=0\nMODEL_INDEX=0\nPROTEIN_INDEX=0\n")


class RMSCalculator:
    @staticmethod
    def calculate_rms_relative_to_aba(aba_score, score):
        return math.sqrt((score**2 - aba_score**2))

    @staticmethod
    def calculate_rms_relative_to_comparison(score, comparison_scores):
        if not comparison_scores:
            return float("inf")
        squared_differences = [(float(score) - float(c_score)) ** 2 for c_score in comparison_scores]
        return math.sqrt(sum(squared_differences) / len(squared_differences))

class DockingProgressMonitor:
    def __init__(self, display_manager):
        self.display_manager = display_manager

    def monitor(self, log_file):
        """
        Monitor the docking progress by reading the log file and updating a text-based progress bar.
        :param log_file: Path to the log file being monitored.
        """
        console_width = shutil.get_terminal_size((80, 20)).columns  # Get terminal width
        progress_bar_length = console_width - 30  # Adjust progress bar length
        progress = 0

        # Wait for the log file to be created
        while not Path(log_file).exists():
            time.sleep(0.2)  # Check every 0.2 seconds

        dock_progress_row = 15  # Row for the progress bar

        # Monitor the log file for progress updates
        while progress < 100:
            if Path(log_file).exists():
                # Count the number of '*' characters in the log file to determine progress
                with open(log_file, "r") as f:
                    progress = f.read().count("*") * 2  # Each '*' represents 2% progress

                # Ensure progress does not exceed 100%
                progress = min(progress, 100)

                # Update the progress bar
                self.display_manager.move_cursor(dock_progress_row, 0)
                print(f"\rCurrent docking progress: {progress}%   ", end="")
                self.display_manager.draw_progress_bar(dock_progress_row + 1, 0, progress_bar_length, progress, 100)

            time.sleep(0.2)  # Update display every 0.2 seconds

        # Ensure the progress bar is fully green at the end
        self.display_manager.move_cursor(dock_progress_row, 0)
        print("\nCurrent docking progress: 100%")
        self.display_manager.draw_progress_bar(dock_progress_row + 1, 0, progress_bar_length, 100, 100)

class DockingConfig:
    def __init__(self, config_path):
        self.config = self.read_config(config_path)
        self.validate_and_report()

    def read_config(self, config_path):
        config = {}
        if not Path(config_path).exists():
            print(f"Warning: Config file {config_path} not found. Using defaults.")
            return config
        with open(config_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "=" in line:
                    key, value = line.split("=", 1)
                    config[key.strip()] = value.strip()
        return config

    def validate_and_report(self):
        expected_keys = ["cpu", "exhaustiveness", "energy_range", "num_modes"]
        print("Configuration loaded:")
        for key in expected_keys:
            if key in self.config:
                print(f"  {key}: {self.config[key]}")
            else:
                print(f"  WARNING: {key} not set in config file!")

class ScoreManager:
    def __init__(self, results_dir):
        self.results_dir = Path(results_dir)
        self.scores_dir = self.results_dir / "scores"

    def write_score(self, protein_name, ligand_name, score):
        # Write or append to scores_<protein_name>.txt
        scores_file = self.scores_dir / f"scores_{protein_name}.txt"
        header = "# Score  Ligand\n"
        if not scores_file.exists() or os.path.getsize(scores_file) == 0:
            with open(scores_file, "w") as f:
                f.write(header)
        with open(scores_file, "a") as f:
            f.write(f"{score} {ligand_name}\n")
        self.update_sorted_scores(scores_file)
        self.write_stats(scores_file)
        self.write_scores_csv(scores_file)

    def update_sorted_scores(self, scores_file):
        # Sort and rewrite scores file
        with open(scores_file, "r") as f:
            lines = [line for line in f if not line.startswith("#")]
        sorted_lines = sorted(lines, key=lambda x: float(x.split()[0]))
        with open(scores_file, "w") as f:
            write_header(f, "Docking scores for all ligands against this protein.", columns=["Score", "Ligand"])
            for line in sorted_lines:
                f.write(line)

    def write_stats(self, scores_file):
        # Write statistics file
        with open(scores_file) as f:
            scores = [float(line.split()[0]) for line in f if not line.startswith("#")]
        if scores:
            mean = statistics.mean(scores)
            median = statistics.median(scores)
            stdev = statistics.stdev(scores) if len(scores) > 1 else 0
            stats_file = scores_file.parent / (scores_file.stem + "_stats.txt")
            with open(stats_file, "w") as f:
                f.write(f"Mean: {mean}\nMedian: {median}\nStdDev: {stdev}\n")

    def write_scores_csv(self, scores_file):
        csv_file = scores_file.with_suffix('.csv')
        with open(scores_file) as fin, open(csv_file, "w") as fout:
            write_header(fout, "CSV version of docking scores for this protein.", columns=["Score", "Ligand"])
            fout.write("Score,Ligand\n")
            for line in fin:
                if line.startswith("#"): continue
                fout.write(",".join(line.strip().split()) + "\n")

    def write_ligand_stats(self, ligand_file, proteins):
        # Write per-ligand stats
        ligand_name = Path(ligand_file).stem
        scores = []
        for protein_file in proteins:
            protein_name = Path(protein_file).stem
            scores_file = self.scores_dir / f"scores_{protein_name}.txt"
            if scores_file.exists():
                with open(scores_file) as f:
                    for line in f:
                        if line.startswith("#"): continue
                        if ligand_name in line:
                            scores.append(float(line.split()[0]))
                            break
        if scores:
            mean = statistics.mean(scores)
            median = statistics.median(scores)
            stdev = statistics.stdev(scores) if len(scores) > 1 else 0
            stats_file = self.scores_dir / "ligands" / f"{ligand_name}_stats.txt"
            stats_file.parent.mkdir(parents=True, exist_ok=True)
            with open(stats_file, "w") as f:
                f.write(f"Mean: {mean}\nMedian: {median}\nStdDev: {stdev}\n")

    def write_ligand_protein_scores(self, ligands, proteins):
        # Write per-ligand stats
        SCORES_DIR = self.scores_dir / "ligands"
        SCORES_DIR.mkdir(parents=True, exist_ok=True)
        for ligand_file in ligands:
            ligand_name = Path(ligand_file).stem
            out_file = SCORES_DIR / f"scores_{ligand_name}.txt"
            with open(out_file, "w") as f:
                f.write("# Score  Protein\n")
                for protein_file in proteins:
                    protein_name = Path(protein_file).stem
                    scores_file = self.scores_dir / f"scores_{protein_name}.txt"
                    score = None
                    if scores_file.exists():
                        with open(scores_file) as sf:
                            for line in sf:
                                if line.startswith("#"): continue
                                if ligand_name in line:
                                    score = line.split()[0]
                                    break
                    if score is not None:
                        f.write(f"{score} {protein_name}\n")

    def write_comparison_ligand_protein_scores(self, comparison_ligands, proteins):
        # Write per-comparison-ligand/protein scores
        SCORES_DIR = self.scores_dir / "comparison_ligands"
        SCORES_DIR.mkdir(parents=True, exist_ok=True)
        for comp_lig_file in comparison_ligands:
            comp_lig_name = Path(comp_lig_file).stem
            out_file = SCORES_DIR / f"scores_{comp_lig_name}.txt"
            with open(out_file, "w") as f:
                f.write("# Score  Protein\n")
                for protein_file in proteins:
                    protein_name = Path(protein_file).stem
                    scores_file = self.scores_dir / f"scores_{protein_name}.txt"
                    score = None
                    if scores_file.exists():
                        with open(scores_file) as sf:
                            for line in sf:
                                if line.startswith("#"): continue
                                if comp_lig_name in line:
                                    score = line.split()[0]
                                    break
                    if score is not None:
                        f.write(f"{score} {protein_name}\n")

    def write_summary_file(self, ligands, proteins, comparison_ligands):
        # Write summary file
        summary_file = self.results_dir / "summary.txt"
        with open(summary_file, "w") as f:
            f.write("Ligands:\n")
            for ligand in ligands:
                f.write(f"  {Path(ligand).stem}\n")
            f.write("\nProteins:\n")
            for protein in proteins:
                f.write(f"  {Path(protein).stem}\n")
            f.write("\nComparison Ligands:\n")
            for comp in comparison_ligands:
                f.write(f"  {Path(comp).stem}\n")

    def write_per_protein_comparison_ligand_rankings(self, ligands, proteins, comparison_ligands):
        # Write per-protein, per-comparison-ligand rankings
        SCORES_DIR = self.scores_dir
        for comp_lig_file in comparison_ligands:
            comp_lig_name = Path(comp_lig_file).stem
            for protein_file in proteins:
                protein_name = Path(protein_file).stem
                scores_file = SCORES_DIR / f"scores_{protein_name}.txt"
                comp_score = None
                if not scores_file.exists():
                    continue
                with open(scores_file) as f:
                    for line in f:
                        if line.startswith("#"): continue
                        if comp_lig_name in line:
                            comp_score = float(line.split()[0])
                            break
                if comp_score is None:
                    continue
                diffs = []
                with open(scores_file) as f:
                    for line in f:
                        if line.startswith("#"): continue
                        score, lig_name = line.split()
                        diff = float(score) - comp_score
                        diffs.append((abs(diff), diff, lig_name))
                diffs.sort()
                ranked_file = SCORES_DIR / f"scores_{comp_lig_name}_in_{protein_name}_ranked.txt"
                with open(ranked_file, "w") as f:
                    f.write(f"# Generated: {time.ctime()}\n")
                    f.write("# ScoreDiff  Ligand\n")
                    for _, diff, lig_name in diffs:
                        f.write(f"{diff:.4f} {lig_name}\n")

    def write_ligand_rms_to_comparison(self, ligands, proteins, comparison_ligands):
        # Write per-ligand RMS to comparison ligands
        SCORES_DIR = self.scores_dir
        OUT_DIR = SCORES_DIR / "ligands"
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        for ligand_file in ligands:
            ligand_name = Path(ligand_file).stem
            out_file = OUT_DIR / f"{ligand_name}_RMS_to_comparison_ligands.txt"
            with open(out_file, "w") as f:
                f.write(f"# Generated: {time.ctime()}\n")
                f.write("# ComparisonLigand  RMS\n")
                for comp_lig_file in comparison_ligands:
                    comp_lig_name = Path(comp_lig_file).stem
                    rms_file = SCORES_DIR / f"scores_{comp_lig_name}_RMS.txt"
                    rms_value = None
                    if rms_file.exists():
                        with open(rms_file) as rf:
                            for line in rf:
                                if line.startswith("#"): continue
                                if ligand_name in line:
                                    rms_value = float(line.split()[0])
                                    break
                    if rms_value is not None:
                        f.write(f"{comp_lig_name} {rms_value:.4f}\n")

    def write_failed_docking_summary(self, failed_dockings):
        # Write failed docking summary
        out_file = self.results_dir / "failed_docking_attempts.txt"
        with open(out_file, "w") as f:
            f.write(f"# Generated: {time.ctime()}\n")
            f.write("# Ligand  Protein\n")
            for ligand, protein in failed_dockings:
                f.write(f"{ligand} {protein}\n")

    @staticmethod
    def write_metadata(dest_dir, vina_command, env_info):
        with open(dest_dir / "metadata.txt", "w") as f:
            f.write("Vina command:\n")
            f.write(" ".join(vina_command) + "\n")
            f.write("Environment:\n")
            f.write(env_info + "\n")

    def calculate_best_ligands(self, ligands, proteins, comparison_ligands):
        best_ligands_file = self.results_dir / "best_ligands.txt"
        with open(best_ligands_file, "w") as best_ligands:
            for ligand_file in ligands:
                models = list((CACHE_DIR / f"models_{Path(ligand_file).stem}").glob("*.pdbqt"))
                for model_file in models:
                    model_name = Path(model_file).stem
                    diffs = []
                    for protein_file in proteins:
                        protein_name = Path(protein_file).stem
                        scores_file = self.scores_dir / f"scores_{protein_name}.txt"
                        model_score = self.get_score(scores_file, model_name)
                        if model_score is None:
                            continue
                        comp_scores = []
                        for comp_lig_file in comparison_ligands:
                            comp_lig_name = Path(comp_lig_file).stem
                            comp_score = self.get_score(scores_file, comp_lig_name)
                            if comp_score is not None:
                                comp_scores.append(comp_score)
                        for comp_score in comp_scores:
                            diffs.append(model_score - comp_score)
                    if diffs:
                        rms = math.sqrt(sum(d**2 for d in diffs) / len(diffs))
                        best_ligands.write(f"{rms:.8f} {model_name}\n")

    def rank_and_display_best_ligands(self):
        """
        Rank the best ligands and display the top 5.
        """
        best_ligands_file = self.results_dir / "best_ligands.txt"
        ranked_best_ligands_file = self.results_dir / "ranked_best_ligands.txt"

        if best_ligands_file.exists():
            with open(best_ligands_file, "r") as f:
                sorted_ligands = sorted(f.readlines(), key=lambda x: float(x.split()[0]))
            with open(ranked_best_ligands_file, "w") as f:
                f.writelines(sorted_ligands)
            print(f"Best ligands ranked and saved to {ranked_best_ligands_file}.")

            # Display the top 5 ligands
            print("Top 5 ligands:")
            for line in sorted_ligands[:5]:
                print(line.strip())  # This will now show model names, e.g., "0.27196004 AAEAMN_model01"
        else:
            print("No best ligands file found. Skipping ranking.")

    def generate_comparison_scores(self, ligands, proteins, comparison_ligands):
        SCORES_DIR = self.scores_dir
        for comp_lig_file in comparison_ligands:
            comp_lig_name = Path(comp_lig_file).stem
            rms_file = SCORES_DIR / f"scores_{comp_lig_name}_RMS.txt"
            rms_scores = []
            for protein_file in proteins:
                protein_name = Path(protein_file).stem
                comp_scores_file = SCORES_DIR / f"scores_{protein_name}.txt"
                comp_score = None
                with open(comp_scores_file) as f:
                    for line in f:
                        if line.startswith("#"): continue
                        if comp_lig_name in line:
                            comp_score = float(line.split()[0])
                            break
                if comp_score is None:
                    continue
                diffs = []
                with open(comp_scores_file) as f:
                    for line in f:
                        if line.startswith("#"): continue
                        score, lig_name = line.split()
                        diff = float(score) - comp_score
                        diffs.append((abs(diff), diff, lig_name))
                diffs.sort()
                per_protein_file = SCORES_DIR / f"scores_{comp_lig_name}" / f"scores_{comp_lig_name}_in_{protein_name}.txt"
                per_protein_file.parent.mkdir(parents=True, exist_ok=True)
                with open(per_protein_file, "w") as f:
                    f.write("# ScoreDiff  Ligand\n")
                    for _, diff, lig_name in diffs:
                        f.write(f"{diff:.4f} {lig_name}\n")
                for _, diff, lig_name in diffs:
                    rms_scores.append((lig_name, diff))
            ligand_rms = {}
            for lig_name, diff in rms_scores:
                ligand_rms.setdefault(lig_name, []).append(diff)
            rms_list = []
            for lig_name, diffs in ligand_rms.items():
                rms = math.sqrt(sum(d**2 for d in diffs) / len(diffs))
                rms_list.append((rms, lig_name))
            rms_list.sort()
            with open(rms_file, "w") as f:
                write_header(f, f"RMS values for all ligands relative to {comp_lig_name}.", columns=["RMS", "Ligand"])
                for rms, lig_name in rms_list:
                    f.write(f"{rms:.4f} {lig_name}\n")

    def generate_ligand_score_summaries(self, ligands, proteins, comparison_ligands):
        SCORES_DIR = self.scores_dir
        DOCKED_LIGANDS_DIR = self.results_dir / "docked_ligands"
        for ligand_file in ligands:
            ligand_name = Path(ligand_file).stem
            # Handle models
            model_dirs = list((DOCKED_LIGANDS_DIR / f"docked_{ligand_name}").glob("docked_*"))
            if model_dirs:
                for model_dir in model_dirs:
                    model_name = model_dir.name.replace("docked_", "")
                    summary_file = model_dir / f"{model_name}_scores.txt"
                    self.write_ligand_summary(summary_file, model_name, proteins, comparison_ligands, SCORES_DIR)
            else:
                ligand_dir = DOCKED_LIGANDS_DIR / f"docked_{ligand_name}"
                summary_file = ligand_dir / f"{ligand_name}_scores.txt"
                self.write_ligand_summary(summary_file, ligand_name, proteins, comparison_ligands, SCORES_DIR)

    def write_ligand_summary(self, summary_file, ligand_name, proteins, comparison_ligands, SCORES_DIR):
        summary_file.parent.mkdir(parents=True, exist_ok=True)  # Ensure directory exists
        with open(summary_file, "w") as f:
            for comp_lig_file in comparison_ligands:
                comp_lig_name = Path(comp_lig_file).stem
                # RMS and rank
                rms_file = SCORES_DIR / f"scores_{comp_lig_name}_RMS.txt"
                rms_rank = None
                rms_value = None
                if rms_file.exists():
                    with open(rms_file) as rf:
                        rms_lines = [line.strip() for line in rf if line.strip()]
                        for idx, line in enumerate(rms_lines):
                            if ligand_name in line:
                                rms_value = float(line.split()[0])
                                rms_rank = idx + 1
                                break
                f.write(f"{comp_lig_name}:\n")
                if rms_value is not None:
                    f.write(f"  RMS: {rms_value:.4f} - #{rms_rank}/{len(rms_lines)}\n")
                else:
                    f.write("  RMS: N/A\n")
                # Per-protein scores and ranks
                for protein_file in proteins:
                    protein_name = Path(protein_file).stem
                    per_protein_file = SCORES_DIR / f"scores_{comp_lig_name}" / f"scores_{comp_lig_name}_in_{protein_name}.txt"
                    if per_protein_file.exists():
                        with open(per_protein_file) as pf:
                            lines = [line.strip() for line in pf if line.strip()]
                            for idx, line in enumerate(lines):
                                if ligand_name in line:
                                    diff = line.split()[0]
                                    rank = idx + 1
                                    f.write(f"  {protein_name}: {diff} - #{rank}/{len(lines)}\n")
                                    break

    def calculate_top_dockers(self, proteins):
        SCORES_DIR = self.scores_dir / "scores_proteins"
        SCORES_DIR.mkdir(parents=True, exist_ok=True)
        for protein_file in proteins:
            protein_name = Path(protein_file).stem
            scores_file = self.scores_dir / f"scores_{protein_name}.txt"
            top_dockers_file = self.scores_dir / f"top_dockers_{protein_name}.txt"
            protein_dir = SCORES_DIR / f"scores_{protein_name}"
            protein_dir.mkdir(parents=True, exist_ok=True)
            unsorted_file = protein_dir / f"scores_{protein_name}_unsorted.txt"
            sorted_file = protein_dir / f"scores_{protein_name}_sorted.txt"

            if scores_file.exists():
                with open(scores_file, "r") as f:
                    lines = [line for line in f if not line.startswith("#")]
                    sorted_scores = sorted(lines, key=lambda x: float(x.split()[0]))
                # Write unsorted with header
                with open(unsorted_file, "w") as f:
                    f.write("# Score  Ligand\n")
                    f.writelines(lines)
                # Write sorted with header
                with open(sorted_file, "w") as f:
                    write_header(f, "Docking scores for all ligands against this protein.", columns=["Score", "Ligand"])
                    for line in sorted_scores:
                        f.write(line)
                # Also update top_dockers_file with header
                with open(top_dockers_file, "w") as f:
                    write_header(f, "Docking scores for all ligands against this protein.", columns=["Score", "Ligand"])
                    f.writelines(sorted_scores)
                self.write_stats(scores_file)
                self.write_scores_csv(scores_file)
                print(f"Top dockers for {protein_name} saved to {top_dockers_file}.")
            else:
                print(f"No scores file found for protein: {protein_name}. Skipping sorting.")

    def get_score(self, scores_file, ligand_name):
        if scores_file.exists():
            with open(scores_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    score, name = line.split()
                    if name == ligand_name:
                        return float(score)
        return None

# create a log folder in the dock folder
if not (BASE_DIR / "log").exists():
    (BASE_DIR / "log").mkdir(parents=True, exist_ok=True)
# create a log file in the log folder
LOG_FILE = BASE_DIR / "log" / "docking_log.txt"
# open the log file in append mode
with open(LOG_FILE, "a") as log:
    log.write("\n")
    log.write(f"Docking log - Start time: {time.ctime(START_TIME)}\n")
    log.write("========================================\n")

class LogManager:
    def __init__(self, log_file):
        self.log_file = log_file

    def log(self, message):
        with open(self.log_file, "a") as log:
            log.write(f"{time.ctime(time.time())}: {message}\n")
            log.flush()

class LigandManager:
    def __init__(self, ligand_dir, cache_dir):
        self.ligand_dir = ligand_dir
        self.cache_dir = cache_dir

    def extract_models(self, ligand_file):
        ligand_name = Path(ligand_file).stem
        output_dir = self.cache_dir / f"models_{ligand_name}"
        output_dir.mkdir(parents=True, exist_ok=True)
        if list(output_dir.glob("*.pdbqt")):
            print(f"Models for {ligand_name} already extracted.")
            return
        with open(ligand_file) as f:
            if "MODEL" in f.read():
                result = subprocess.run(
                    ["vina_split", "--input", ligand_file, "--ligand", str(output_dir / f"{ligand_name}_model")],
                    capture_output=True,
                    text=True
                )
                if result.returncode != 0:
                    print(f"Error: Failed to split ligand file {ligand_file}")
                    exit(1)
                if not list(output_dir.glob(f"{ligand_name}_model*.pdbqt")):
                    print(f"Error: No models generated for ligand file {ligand_file}")
                    exit(1)
                for model_file in output_dir.glob(f"{ligand_name}_model*.pdbqt"):
                    stem_parts = model_file.stem.split("_model")
                    if len(stem_parts) == 2 and stem_parts[1].isdigit():
                        model_index = stem_parts[1]
                    else:
                        model_index = model_file.stem.split("_")[-1]
                    new_name = f"{ligand_name}_model_{model_index}.pdbqt"
                    model_file.rename(output_dir / new_name)
            else:
                # Use base name for single-model ligands
                shutil.copy(ligand_file, output_dir / f"{ligand_name}.pdbqt")

    @staticmethod
    def clean_ligand_name(name):
        return name.split(".")[0].replace(" ", "_")

class ProteinManager:
    def __init__(self, protein_dir, config_dir):
        self.protein_dir = protein_dir
        self.config_dir = config_dir

    def get_box_from_config(self, protein_name):
        config_path = self.config_dir / f"config_{protein_name}.txt"
        if config_path.exists():
            params = {}
            with open(config_path) as f:
                for line in f:
                    if "=" in line:
                        k, v = line.strip().split("=", 1)
                        params[k.strip()] = float(v.strip())
            return (
                params.get("center_x"),
                params.get("center_y"),
                params.get("center_z"),
                params.get("size_x"),
                params.get("size_y"),
                params.get("size_z"),
            )
        return None

    def calculate_docking_box(self, protein_file):
        protein_name = Path(protein_file).stem
        box = self.get_box_from_config(protein_name)
        if box and all(v is not None for v in box):
            return box
        # fallback to automatic calculation
        center_x, center_y, center_z = 0, 0, 0
        count = 0
        with open(protein_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    parts = line.split()
                    center_x += float(parts[6])
                    center_y += float(parts[7])
                    center_z += float(parts[8])
                    count += 1
        if count > 0:
            center_x /= count
            center_y /= count
            center_z /= count
        size_x, size_y, size_z = 20, 20, 20
        return center_x, center_y, center_z, size_x, size_y, size_z

class ProgressManager:
    def __init__(self, cache_file, score_manager):
        self.cache_file = cache_file
        self.score_manager = score_manager

    def read_progress(self):
        if self.cache_file.exists():
            progress = {}
            with open(self.cache_file) as f:
                for line in f:
                    key, value = line.strip().split("=")
                    progress[key] = int(value)
            return progress
        else:
            print("Error: Progress cache file not found.")
            exit(1)

    def update_progress(self, progress):
        with open(self.cache_file, "w") as f:
            for key, value in progress.items():
                f.write(f"{key}={value}\n")

    def terminate_script(self, signal_received, frame):
        print("\nTermination signal received. Saving progress...")
        progress = {
            "COMPARISON_LIGAND_INDEX": globals().get("COMPARISON_LIGAND_INDEX", 0),
            "COMPARISON_PROTEIN_INDEX": globals().get("COMPARISON_PROTEIN_INDEX", 0),
            "LIGAND_INDEX": globals().get("LIGAND_INDEX", 0),
            "MODEL_INDEX": globals().get("MODEL_INDEX", 0),
            "PROTEIN_INDEX": globals().get("PROTEIN_INDEX", 0)
        }
        self.update_progress(progress)
        self.score_manager.write_failed_docking_summary(FAILED_DOCKINGS)
        exit(0)

    def update_from_globals(self):
        self.update_progress({
            "COMPARISON_LIGAND_INDEX": globals().get("COMPARISON_LIGAND_INDEX", 0),
            "COMPARISON_PROTEIN_INDEX": globals().get("COMPARISON_PROTEIN_INDEX", 0),
            "LIGAND_INDEX": globals().get("LIGAND_INDEX", 0),
            "MODEL_INDEX": globals().get("MODEL_INDEX", 0),
            "PROTEIN_INDEX": globals().get("PROTEIN_INDEX", 0)
        })

class DockingTask:
    def __init__(self, ligand_file, protein_file, model_index, config, results_dir, debug, score_manager):
        self.ligand_file = ligand_file
        self.protein_file = protein_file
        self.model_index = model_index
        self.config = config
        self.results_dir = results_dir
        self.debug = debug
        self.score_manager = score_manager

    def run(self):
        ligand_name = Path(self.ligand_file).stem
        protein_name = Path(self.protein_file).stem
        output_file = self.results_dir / f"temp/{ligand_name}_model{self.model_index}_vs_{protein_name}.pdbqt"
        log_file = self.results_dir / f"temp/{ligand_name}_model{self.model_index}_vs_{protein_name}.log"

        (self.results_dir / "temp").mkdir(parents=True, exist_ok=True)
        protein_manager = ProteinManager(PROTEIN_DIR, CONFIG_DIR)
        center_x, center_y, center_z, size_x, size_y, size_z = protein_manager.calculate_docking_box(self.protein_file)

        if self.debug:
            print("Running AutoDock Vina with the following parameters:")
            print(f"Receptor: {self.protein_file}")
            print(f"Ligand: {self.ligand_file}")
            print(f"Center: {center_x}, {center_y}, {center_z}")
            print(f"Size: {size_x}, {size_y}, {size_z}")
            print(f"Output: {output_file}")
            print(f"Log: {log_file}")

        vina_command = [
            "vina",
            "--receptor", str(self.protein_file),
            "--ligand", str(self.ligand_file),
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--out", str(output_file)
        ]

        if "cpu" in self.config:
            vina_command += ["--cpu", self.config["cpu"]]
        if "exhaustiveness" in self.config:
            vina_command += ["--exhaustiveness", self.config["exhaustiveness"]]
        if "energy_range" in self.config:
            vina_command += ["--energy_range", self.config["energy_range"]]
        if "num_modes" in self.config:
            vina_command += ["--num_modes", self.config["num_modes"]]

        with open(log_file, "w") as log:
            vina_process = subprocess.Popen(vina_command, stdout=log, stderr=log)

        docking_monitor = DockingProgressMonitor(display_manager)
        docking_monitor.monitor(log_file)
        vina_process.wait()

        if not log_file.exists():
            print(f"Error: Log file not found: {log_file}")
            FAILED_DOCKINGS.append((ligand_name, protein_name))
            return

        if self.debug:
            display_manager.move_cursor(shutil.get_terminal_size((80, 20)).lines - 1, 0)
            print(f"Contents of log file {log_file}:")
            with open(log_file, "r") as log:
                print(log.read())

        score = None
        with open(log_file, "r") as log:
            for line in log:
                if line.startswith("   1"):
                    score = line.split()[1]
                    break

        log_manager.log(f"Docking completed for {ligand_name} with {protein_name}. Score: {score}")

        if score and (score.replace(".", "", 1).replace("-", "", 1).isdigit() or score == "-"):
            log_manager.log(f"   Valid score extracted: {score} for ligand: {ligand_name} with protein: {protein_name}")
            if self.debug:
                print(f"Extracted score: {score} for ligand: {ligand_name} with protein: {protein_name}")
            SCORES_DIR = self.results_dir / "scores"
            if not SCORES_DIR.exists():
                SCORES_DIR.mkdir(parents=True, exist_ok=True)
            scores_file = SCORES_DIR / f"scores_{protein_name}.txt"
            header = "# Score  Ligand\n"
            if not scores_file.exists() or os.path.getsize(scores_file) == 0:
                with open(scores_file, "w") as f:
                    f.write(header)
            with open(scores_file, "a") as f:
                f.write(f"{score} {ligand_name}\n")
            self.score_manager.update_sorted_scores(scores_file)
            self.score_manager.write_stats(scores_file)
            self.score_manager.write_scores_csv(scores_file)
            if self.debug:
                print(f"Score stored in {scores_file}")
        else:
            log_manager.log(f"Failed to extract a valid score from log file: {log_file}")
            print(f"Failed to extract a valid score from log file: {log_file}")
            FAILED_DOCKINGS.append((ligand_name, protein_name))

        ligand_name = Path(self.ligand_file).stem.split(".")[0]
        model_part = ""
        if "_model_" in Path(self.ligand_file).stem:
            model_part = Path(self.ligand_file).stem.split("_model_")[-1]
            dest_dir = self.results_dir / "docked_ligands" / f"docked_{ligand_name}" / f"docked_{ligand_name}_model{model_part}"
        else:
            dest_dir = self.results_dir / "docked_ligands" / f"docked_{ligand_name}"
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.move(str(output_file), dest_dir / output_file.name)
        shutil.move(str(log_file), dest_dir / log_file.name)

        env_info = f"Python version: {sys.version}"
        self.score_manager.write_metadata(dest_dir, vina_command, env_info)

class ResultOrganizer:
    def __init__(self, results_dir):
        self.results_dir = results_dir

    def move_temp_files(self):
        temp_dir = self.results_dir / "temp"
        DOCKED_LIGANDS_DIR = self.results_dir / "docked_ligands"
        DOCKED_LIGANDS_DIR.mkdir(parents=True, exist_ok=True)
        for file in temp_dir.glob("*"):
            parts = file.stem.split("_vs_")[0].split("_model")
            ligand_base = parts[0]
            model_part = f"model{parts[1]}" if len(parts) > 1 else None
            ligand_base_clean = ligand_base.split(".")[0]
            if model_part:
                dest_dir = DOCKED_LIGANDS_DIR / f"docked_{ligand_base_clean}" / f"docked_{ligand_base_clean}_{model_part}"
            else:
                dest_dir = DOCKED_LIGANDS_DIR / f"docked_{ligand_base_clean}"
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(file), dest_dir / file.name)

class DisplayManager:
    def __init__(self):
        pass

    def display_progress(self, current_task, total_tasks, ligand_index, total_ligands, model_index, total_models, protein_index, total_proteins, ligand_name, protein_name):
        """
        Display the progress of the docking process with an estimated time to completion (ETC).
        :param current_task: Current task number.
        :param total_tasks: Total number of tasks.
        :param ligand_index: Current ligand index.
        :param total_ligands: Total number of ligands.
        :param model_index: Current model index.
        :param total_models: Total number of models.
        :param protein_index: Current protein index.
        :param total_proteins: Total number of proteins.
        :param ligand_name: Name of the current ligand.
        :param protein_name: Name of the current protein.
        """
        console_width = shutil.get_terminal_size((80, 20)).columns
        progress_bar_length = console_width - 30

        # Calculate elapsed time for the last 10 tasks
        if TASK_DURATIONS:
            avg_task_time = sum(TASK_DURATIONS) / len(TASK_DURATIONS)
        else:
            avg_task_time = 0

        remaining_tasks = total_tasks - current_task
        remaining_time = avg_task_time * remaining_tasks if avg_task_time > 0 else 0

        # Format remaining time as HH:MM:SS
        hours, rem = divmod(int(remaining_time), 3600)
        minutes, seconds = divmod(rem, 60)
        etc = f"{hours:02}:{minutes:02}:{seconds:02}"

        print("\033[H\033[J", end="")  # Clear the terminal
        print("========================================")
        print("         Protein-Ligand Docking         ")
        print("========================================")
        print()
        print(f"Docking:")
        print(f"Ligand file {ligand_index + 1}/{total_ligands}")
        self.draw_progress_bar(7, 0, progress_bar_length, ligand_index + 1, total_ligands, False)

        self.move_cursor(9, 0)
        print(f"Ligand model: model {model_index + 1}/{total_models}")
        self.draw_progress_bar(10, 0, progress_bar_length, model_index + 1, total_models, False)

        self.move_cursor(12, 0)
        print(f"Ligand progress: protein {protein_index + 1}/{total_proteins}")
        self.draw_progress_bar(13, 0, progress_bar_length, protein_index + 1, total_proteins, False)

        # Current docking progress bar
        self.move_cursor(15, 0)
        print("Current docking progress:")
        self.draw_progress_bar(16, 0, progress_bar_length, 0, 100, False)

        self.move_cursor(19, 0)
        print(f"Ligand name: {ligand_name}")
        print(f"Protein name: {protein_name}")

        percent_through_dock = (current_task * 100) / total_tasks
        self.move_cursor(21, 0)
        print(f"Total progress: {current_task}/{total_tasks}")
        self.draw_progress_bar(22, 0, progress_bar_length, current_task, total_tasks, False)

        # Display estimated time to completion
        self.move_cursor(24, 0)
        print(f"Estimated time to completion: {etc}")

        self.move_cursor(26, 0)
        print("Do CTRL+C to exit")

    def display_comparison_progress(self, current_task, total_tasks, comparison_ligand_index, total_comparison_ligands, protein_index, total_proteins, comparison_ligand_name, protein_name):
        """
        Display the progress of docking for comparison ligands.
        :param current_task: Current task number.
        :param total_tasks: Total number of tasks.
        :param comparison_ligand_index: Current comparison ligand index.
        :param total_comparison_ligands: Total number of comparison ligands.
        :param protein_index: Current protein index.
        :param total_proteins: Total number of proteins.
        :param comparison_ligand_name: Name of the current comparison ligand.
        :param protein_name: Name of the current protein.
        """
        console_width = shutil.get_terminal_size((80, 20)).columns
        progress_bar_length = console_width - 30

        print("\033[H\033[J", end="")  # Clear the terminal
        print("========================================")
        print("   Initial Calculation of Comparison    ")
        print("========================================")
        print()
        print(f"Docking comparisons:")
        print(f"Total progress: comparison ligand {comparison_ligand_index + 1}/{total_comparison_ligands}")
        self.draw_progress_bar(7, 0, progress_bar_length, comparison_ligand_index + 1, total_comparison_ligands, False)
        self.move_cursor(9, 0)
        print(f"Comparison ligand progress: protein {protein_index + 1}/{total_proteins}")
        self.draw_progress_bar(10, 0, progress_bar_length, protein_index + 1, total_proteins, False)
        self.move_cursor(19, 0)
        print(f"Comparison ligand name: {comparison_ligand_name}")
        print(f"Protein name: {protein_name}")

        percent_through_dock = (current_task * 100) / total_tasks
        self.move_cursor(21, 0)
        print(f"% through dock: {percent_through_dock:.2f}%")
        self.draw_progress_bar(22, 0, progress_bar_length, current_task, total_tasks, False)

        self.move_cursor(24, 0)
        print("Do CTRL+C to exit")
        print()

    def draw_progress_bar(self, row, col, length, current, total, waits=False):
        """
        Draw a progress bar in the terminal.
        :param row: Cursor row position for the progress bar.
        :param col: Cursor column position for the progress bar.
        :param length: Length of the progress bar.
        :param current: Current progress value.
        :param total: Total progress value.
        :param waits: Optional flag to wait after each character (default: False).
        """
        # Check if the current value exceeds the total value
        if current > total:
            current = total

        # Check if the current value is negative
        if current < 0:
            current = 0

        # Check if the progress bar length is valid
        if length <= 0:
            print("Error: Progress bar length must be greater than 0.")
            return

        # Check if the current and total values are valid
        if current < 0 or total <= 0:
            print("Error: Progress bar current and total values must be non-negative, and total must be greater than 0.")
            return

        # Calculate progress as a fraction (0 to 1)
        progress = current / total

        # Calculate the number of filled and empty segments
        filled = int(progress * length)
        empty = length - filled

        # Draw the progress bar
        bar = "[" + "\033[42m \033[0m" * filled + "\033[41m \033[0m" * empty + "]"
        percent = f"{progress * 100:.2f}%"

        # Move the cursor to the specified position and print the progress bar
        print(f"\033[{row};{col}H{bar} {percent}", end="", flush=True)

    def clear_display_area(self):
        """
        Clear the display area in the terminal.
        """
        console_height = shutil.get_terminal_size((80, 20)).lines
        console_width = shutil.get_terminal_size((80, 20)).columns
        print("\033[H", end="")  # Move cursor to the top left corner
        for _ in range(console_height):
            print(" " * console_width)
        print("\033[H", end="")  # Move cursor back to the top left corner

    def move_cursor(self, row, col):
        """
        Move the cursor to a specific row and column in the terminal.
        :param row: The row number (1-based).
        :param col: The column number (1-based).
        """
        print(f"\033[{row};{col}H", end="", flush=True)

    def add_empty_lines(self):
        """
        Add empty lines to create space for display updates.
        """
        console_height = shutil.get_terminal_size((80, 20)).lines
        for _ in range(console_height):
            print()

    @staticmethod
    def wait_for_render_stop():
        global IS_RENDERING
        while IS_RENDERING:
            time.sleep(0.01)

class CacheManager:
    def __init__(self, cache_dir, results_dir, ligand_dir, protein_dir, comparison_ligand_dir):
        self.cache_dir = cache_dir
        self.results_dir = results_dir
        self.ligand_dir = ligand_dir
        self.protein_dir = protein_dir
        self.comparison_ligand_dir = comparison_ligand_dir

    def initialize_cache(self, mode=None):
        if mode == "clear":
            backup_dir = self.cache_dir / "cache_backup"
            backup_dir.mkdir(parents=True, exist_ok=True)
            for item in self.cache_dir.iterdir():
                if item == backup_dir:
                    continue
                destination = backup_dir / item.name
                if destination.exists():
                    counter = 1
                    new_name = f"{item.stem}_copy{counter}{item.suffix}"
                    new_destination = backup_dir / new_name
                    while new_destination.exists():
                        counter += 1
                        new_name = f"{item.stem}_copy{counter}{item.suffix}"
                        new_destination = backup_dir / new_name
                    destination = new_destination
                shutil.move(str(item), str(destination))
            print("Cache cleared and backed up.")
        elif mode == "clear-everything":
            shutil.rmtree(self.cache_dir, ignore_errors=True)
            shutil.rmtree(self.results_dir, ignore_errors=True)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            self.results_dir.mkdir(parents=True, exist_ok=True)
            print("Cache and results cleared.")

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # Generate ligand and protein name lists
        ligand_names = self.cache_dir / "ligandNames.txt"
        protein_names = self.cache_dir / "proteinNames.txt"
        comparison_ligand_names = self.cache_dir / "comparisonLigandNames.txt"
        progress_cache = self.cache_dir / "progress_cache.txt"

        if not list(self.ligand_dir.glob("*.pdbqt")):
            print(f"Error: No ligand files found in {self.ligand_dir}")
            exit(1)
        with open(ligand_names, "w") as f:
            f.writelines([str(file) + "\n" for file in self.ligand_dir.glob("*.pdbqt")])

        if not list(self.protein_dir.glob("*.pdbqt")):
            print(f"Error: No protein files found in {self.protein_dir}")
            exit(1)
        with open(protein_names, "w") as f:
            f.writelines([str(file) + "\n" for file in self.protein_dir.glob("*.pdbqt")])

        if not list(self.comparison_ligand_dir.glob("*.pdbqt")):
            print(f"Error: No comparison ligand files found in {self.comparison_ligand_dir}")
            exit(1)
        with open(comparison_ligand_names, "w") as f:
            f.writelines([str(file) + "\n" for file in self.comparison_ligand_dir.glob("*.pdbqt")])

        if not progress_cache.exists():
            with open(progress_cache, "w") as f:
                f.write("COMPARISON_LIGAND_INDEX=0\nCOMPARISON_PROTEIN_INDEX=0\nLIGAND_INDEX=0\nMODEL_INDEX=0\nPROTEIN_INDEX=0\n")

def write_header(f, description, columns=None):
    f.write(f"# Generated: {time.ctime()}\n")
    f.write(f"# {description}\n")
    if columns:
        if f.name.endswith('.csv'):
            f.write("# " + ",".join(columns) + "\n")
        else:
            f.write("# " + "  ".join(columns) + "\n")

# Register the termination signal handler
progress_manager = ProgressManager(PROGRESS_CACHE, ScoreManager(RESULTS_DIR))
signal.signal(signal.SIGINT, progress_manager.terminate_script)
signal.signal(signal.SIGTERM, progress_manager.terminate_script)

# Initialize DisplayManager
display_manager = DisplayManager()

# Initialize LogManager
log_manager = LogManager(LOG_FILE)

# Instantiate LigandManager
ligand_manager = LigandManager(LIGAND_DIR, CACHE_DIR)

# Main script execution
if __name__ == "__main__":
    import sys

    # Handle command-line arguments for clearing cache
    cache_manager = CacheManager(CACHE_DIR, RESULTS_DIR, LIGAND_DIR, PROTEIN_DIR, COMPARISON_LIGAND_DIR)
    if "--clear-cache" in sys.argv:
        cache_manager.initialize_cache("clear")
    elif "--clear-everything" in sys.argv:
        cache_manager.initialize_cache("clear-everything")
    else:
        cache_manager.initialize_cache()

    # Read progress from cache
    progress = progress_manager.read_progress()
    COMPARISON_LIGAND_INDEX = progress.get("COMPARISON_LIGAND_INDEX", 0)
    COMPARISON_PROTEIN_INDEX = progress.get("COMPARISON_PROTEIN_INDEX", 0)
    LIGAND_INDEX = progress.get("LIGAND_INDEX", 0)
    MODEL_INDEX = progress.get("MODEL_INDEX", 0)
    PROTEIN_INDEX = progress.get("PROTEIN_INDEX", 0)

    # Load ligands, proteins, and comparison ligands
    ligands = [line.strip() for line in open(LIGAND_NAMES)]
    proteins = [line.strip() for line in open(PROTEIN_NAMES)]
    comparison_ligands = [line.strip() for line in open(COMPARISON_LIGAND_NAMES)]

    total_ligands = len(ligands)
    total_proteins = len(proteins)
    total_comparison_ligands = len(comparison_ligands)

    # After loading proteins
    missing_configs = []
    unfilled_configs = []

    for protein_file in proteins:
        protein_name = Path(protein_file).stem
        config_path = CONFIG_DIR / f"config_{protein_name}.txt"
        if not config_path.exists():
            with open(config_path, "w") as f:
                f.write(
                    "center_x=\n"
                    "center_y=\n"
                    "center_z=\n"
                    "size_x=\n"
                    "size_y=\n"
                    "size_z=\n"
                    "# Fill in the bounding box values above for this protein.\n"
                )
            missing_configs.append(str(config_path))
        else:
            # Check if all required fields are filled
            with open(config_path) as f:
                content = f.read()
                for key in ["center_x=", "center_y=", "center_z=", "size_x=", "size_y=", "size_z="]:
                    if f"{key}\n" in content or f"{key}\r\n" in content:  # still blank
                        unfilled_configs.append(str(config_path))
                        break

    if missing_configs:
        print("The following protein config files were created and must be filled in before running the program again:")
        for path in missing_configs:
            print("  ", path)
        exit(1)

    if unfilled_configs:
        print("The following protein config files are missing required bounding box values:")
        for path in unfilled_configs:
            print("  ", path)
        print("Please fill in all center_x, center_y, center_z, size_x, size_y, and size_z values before running the program again.")
        exit(1)

    # Create a dedicated subfolder for output files
    (RESULTS_DIR / "temp").mkdir(parents=True, exist_ok=True)

    # Extract models for all ligands first
    for ligand_file in ligands:
        ligand_manager.extract_models(ligand_file)

    # Now count models for all ligands
    all_ligands_models = []
    for ligand_file in ligands:
        models = list((CACHE_DIR / f"models_{Path(ligand_file).stem}").glob("*.pdbqt"))
        all_ligands_models.append(models)

    total_tasks = sum(len(models_for_ligand) * total_proteins for models_for_ligand in all_ligands_models)

    if total_tasks == 0:
        print("No docking tasks to perform. Check your ligand/model extraction.")
        exit(1)

    # Perform docking for comparison ligands
    while COMPARISON_LIGAND_INDEX < total_comparison_ligands:
        comparison_ligand_file = comparison_ligands[COMPARISON_LIGAND_INDEX]
        comparison_ligand_name = Path(comparison_ligand_file).stem

        while COMPARISON_PROTEIN_INDEX < total_proteins:
            protein_file = proteins[COMPARISON_PROTEIN_INDEX]
            protein_name = Path(protein_file).stem

            # Check if docking for this comparison ligand and protein has already been completed
            scores_file = RESULTS_DIR / f"scores_{protein_name}.txt"
            if scores_file.exists():
                with open(scores_file, "r") as f:
                    if any(comparison_ligand_name in line for line in f):
                        print(f"Skipping docking for {comparison_ligand_name} with {protein_name} (already completed).")
                        COMPARISON_PROTEIN_INDEX += 1
                        progress_manager.update_from_globals()
                        continue

            # Display progress for comparison ligands
            display_manager.display_comparison_progress(
                current_task=COMPARISON_LIGAND_INDEX * total_proteins + COMPARISON_PROTEIN_INDEX + 1,
                total_tasks=total_comparison_ligands * total_proteins,
                comparison_ligand_index=COMPARISON_LIGAND_INDEX,
                total_comparison_ligands=total_comparison_ligands,
                protein_index=COMPARISON_PROTEIN_INDEX,
                total_proteins=total_proteins,
                comparison_ligand_name=comparison_ligand_name,
                protein_name=protein_name
            )

            # Perform docking
            docking_task = DockingTask(comparison_ligand_file, protein_file, 0, {}, RESULTS_DIR, DEBUG, ScoreManager(RESULTS_DIR))
            docking_task.run()

            # Update progress after each docking task
            COMPARISON_PROTEIN_INDEX += 1
            progress_manager.update_from_globals()

        # Reset protein index and move to the next comparison ligand
        COMPARISON_PROTEIN_INDEX = 0
        COMPARISON_LIGAND_INDEX += 1
        progress_manager.update_from_globals()

    # Clear the display area after comparison ligand docking
    if(not DEBUG):
        display_manager.clear_display_area()

    # start time for docking
    start_time = time.time()

    # Perform docking for ligands
    while LIGAND_INDEX < total_ligands:
        ligand_file = ligands[LIGAND_INDEX]
        ligand_manager.extract_models(ligand_file)  # <-- FIXED
        models = list((CACHE_DIR / f"models_{Path(ligand_file).stem}").glob("*.pdbqt"))
        for model_file in models:
            model_name = Path(model_file).stem
        total_models = len(models)

        current_task = (
            LIGAND_INDEX * total_models * total_proteins +
            MODEL_INDEX * total_proteins +
            PROTEIN_INDEX + 1
        )

        while MODEL_INDEX < total_models:
            model_file = models[MODEL_INDEX]

            while PROTEIN_INDEX < total_proteins:
                protein_file = proteins[PROTEIN_INDEX]

                # Update ligand and protein names for progress display
                ligand_name = Path(model_file).stem
                protein_name = Path(protein_file).stem

                # Display progress before starting docking
                display_manager.display_progress(
                    current_task=current_task,
                    total_tasks=total_tasks,
                    ligand_index=LIGAND_INDEX,
                    total_ligands=total_ligands,
                    model_index=MODEL_INDEX,
                    total_models=total_models,
                    protein_index=PROTEIN_INDEX,
                    total_proteins=total_proteins,
                    ligand_name=ligand_name,
                    protein_name=protein_name
                )

                # Perform docking
                docking_task = DockingTask(model_file, protein_file, MODEL_INDEX, {}, RESULTS_DIR, DEBUG, ScoreManager(RESULTS_DIR))
                docking_task.run()

                end_time = time.time()
                # Record task duration
                TASK_DURATIONS.append(end_time - start_time)
                start_time = end_time
                # Update progress
                current_task += 1
                PROTEIN_INDEX += 1
                progress_manager.update_from_globals()

            PROTEIN_INDEX = 0
            MODEL_INDEX += 1
            progress_manager.update_from_globals()

        MODEL_INDEX = 0
        LIGAND_INDEX += 1
        progress_manager.update_from_globals()

    # Move temporary files to the results directory
    result_organizer = ResultOrganizer(RESULTS_DIR)
    result_organizer.move_temp_files()

    # Initialize ScoreManager
    score_manager = ScoreManager(RESULTS_DIR)

    # Calculate and save top dockers for each protein
    score_manager.calculate_top_dockers(proteins)

    # Generate comparison scores
    score_manager.generate_comparison_scores(ligands, proteins, comparison_ligands)

    # Calculate and save the best ligands
    score_manager.calculate_best_ligands(ligands, proteins, comparison_ligands)

    # Rank and display the best ligands
    score_manager.rank_and_display_best_ligands()

    # Generate per-ligand/model score summaries
    score_manager.generate_ligand_score_summaries(ligands, proteins, comparison_ligands)

    # Write ligand-protein scores
    score_manager.write_ligand_protein_scores(ligands, proteins)

    # Write comparison ligand-protein scores
    score_manager.write_comparison_ligand_protein_scores(comparison_ligands, proteins)

    # Write summary file
    score_manager.write_summary_file(ligands, proteins, comparison_ligands)

    # Write per-protein comparison ligand rankings
    score_manager.write_per_protein_comparison_ligand_rankings(ligands, proteins, comparison_ligands)

    # Write RMS to comparison ligands for each ligand
    score_manager.write_ligand_rms_to_comparison(ligands, proteins, comparison_ligands)

    for ligand_file in ligands:
        score_manager.write_ligand_stats(ligand_file, proteins)

    # Write failed docking summary
    score_manager.write_failed_docking_summary(FAILED_DOCKINGS)

    # Write CSVs for all scores_<protein_name>.txt files
    for protein_file in proteins:
        protein_name = Path(protein_file).stem
        scores_file = RESULTS_DIR / "scores" / f"scores_{protein_name}.txt"
        if scores_file.exists():
            score_manager.write_scores_csv(scores_file)

    # After all scoring and summaries
    ranked_dir = RESULTS_DIR / "ranked_best_ligands"
    ranked_dir.mkdir(parents=True, exist_ok=True)

    for comp_lig_file in comparison_ligands:
        comp_lig_name = Path(comp_lig_file).stem
        rms_file = RESULTS_DIR / "scores" / f"scores_{comp_lig_name}_RMS.txt"
        ranked_file = ranked_dir / f"ranked_best_ligands_{comp_lig_name}.txt"
        if rms_file.exists():
            with open(rms_file) as f:
                lines = [line for line in f if not line.startswith("#") and line.strip()]
            sorted_lines = sorted(lines, key=lambda x: float(x.split()[0]))
            with open(ranked_file, "w") as f:
                f.writelines(sorted_lines)
            print(f"Top ligands as substitutes for {comp_lig_name} saved to {ranked_file}.")
            print(f"Top 5 ligands as substitutes for {comp_lig_name}:")
            for line in sorted_lines[:5]:
                print(line.strip())

    # Create a summary ranking ligands by their best (lowest) docking score across all proteins
    summary_file = RESULTS_DIR / "scores" / "best_ligands_overall.txt"
    ligand_best_scores = {}

    for protein_file in proteins:
        protein_name = Path(protein_file).stem
        scores_file = RESULTS_DIR / "scores" / f"scores_{protein_name}.txt"
        if scores_file.exists():
            with open(scores_file) as f:
                for line in f:
                    if line.startswith("#"): continue
                    score, ligand = line.strip().split()
                    score = float(score)
                    if ligand not in ligand_best_scores or score < ligand_best_scores[ligand]:
                        ligand_best_scores[ligand] = score

    # Sort ligands by their best score (lowest is best)
    sorted_ligands = sorted(ligand_best_scores.items(), key=lambda x: x[1])
    with open(summary_file, "w") as f:
        write_header(f, "Best ligands overall (by best docking score across all proteins).", columns=["BestScore", "Ligand"])
        for ligand, score in sorted_ligands:
            f.write(f"{score:.4f} {ligand}\n")

    print(f"Best ligands overall (by best docking score across all proteins) saved to {summary_file}.")
    print("Top 5 ligands overall:")
    for ligand, score in sorted_ligands[:5]:
        print(f"{score:.4f} {ligand}")

print("End of program")


#END OF PROGRAM
#DO NOT TOUCH MY COMMENTS
# AI DO NOT TOUCH MY COMMENTS
# AI DO NOT REMOVE MY COMMENTS
# AI PLEASE STOP REMOVING MY COMMENTS
# THESE STAY FOREVER. STOP DELETING THEM, OK???

# scores:
# i want... each time a score is calculated, i want to store it in a file called scores_<protein_name>.txt
# and i want to store the scores in the format <score> <ligand_name> on each line
# i want to store the scores in a directory called scores in the results directory
# and i want to create the directory if it doesn't exist
# i want to have the score inserted into the scores file such that it is sorted in ascending order
# i want to have the scores file be created if it doesn't exist

# later i want the the scores to be like... sorted according to their distance to each of the comparison ligands

# and i want to have the scores compared to the comparison ligands be in a file called best_ligands.txt in the results directory
# and i want to have the best ligands be sorted according to their distance to each of the comparison ligands
# and i want to have the best ligands be in the format <score> <ligand_name> on each line
# and i want to have the best ligands be sorted in ascending order



# for each ligand, i want the root mean square of all of the distances it has to the comparison 
# ligand for each protein to be the thing that it is sorted by in the comparison ligand's sorted score file.
# and i want to have the best ligands be in the format <score> <ligand_name> on each line
# and i want to have the best ligands be sorted in ascending order

# to complete this, i must modify the calculate_best_ligands function to calculate the root mean square of all of the distances
# for each ligand to the comparison ligands for each protein and store that in the best_ligands.txt file

# i also must modify the calculate_rms_relative_to_comparison function to calculate the root mean square of all of the distances
# for each ligand to the comparison ligands for each protein and store that in the best_ligands.txt file



# file structure of program

# or like... the scores thing

# so the scores... so
# they're calculated

# and then once they're calculated
# like you have for each comparison ligand, for each protein, there is a score
# now you have calculated each ligand to each protein
# so now for each of those scores you get to each protein, you need to calculate the distance or 
# whatever to like that of each comparison ligand



# process:::::::::

# Comparison Ligands protien scores per ligand:
# for each comparison ligand:
#   for each protein:
#       get score for that comparison ligand to that protein
#       save the score in directory dock/results/scores/comparison_ligands in file scores_<comparisonLigand_name>.txt
#       save the score in the format <score> <protein_name> on each line

# Ligands protien scores per ligand:
# for each ligand: 
#   for each protein:
#       get score for that ligand to that protein (score 1)
#       for each comparison ligand with current protein:
#           CompLigand[i]Protein[j] store comparisonLigandScore-ligandScore
#           
#           
#           score for comparison ligand to current protein is already saved (score 2)
#           calculate the distance between the two scores (score 1 and score 2)
#           store the distance in a list for that comparison ligand
#           
#           store the scores in a list for that ligand
#           store this score in a file called scores_<ligand_name>.txt in the dock/results/scores/ligands directory
#           store the score in the format <score> <protein_name> on each line


# for each ligand:
#   fetch the ligand 
#   for each comparison ligand:
#       
#       calculate the root mean square of the distances for each protein
#       have that be the value that it is sorted by in the best_ligands.txt file













# ok no more autocomplete with AI, that's very distracting

# desired results:
# see how good a ligand is as a substitute for a comparison ligand
# that is known by the root mean square of ligandScore-comparisonLigandScore for all proteins
# so, calculate that with each ligand-comparisonLigand pair, and generate a scoreboard or whatever for each substitute ligand,
# to find which ligands best work as substitutes for each comparisonLigand
# also generate a unique scoreboard for each comparisonLigand-protein pair, 
# such that it is shown which ligands are best as substitutes for a comparisonLigand for each individual protein, 
# because that might be useful as well


# so how will i do that?
# since there will likely be less comparisonLigands than proteins, i will have folders for each comparisonLigand
# that each contain all individual scores for each protein and the composite score using the RMS of the other files' data.
# Each file will store data as <score-comparisonScore> <protein_name> on each line, sorted by magnitude of <score-comparisonScore>,
# such that lower magnitudes are ranked as better, higher up in the list.

# I have it ordered score-comparisonScore, so that if score is lower than comparisonScore, then the resulting value is negative, showing
# that score is so much score below that of comparisonScore, and vice versa.

# The file names for the scores file will have a consistent structure:
# scores_<comparisonLigand_name>_<protein_name>.txt


# For the actual docked ligand files, there will exist a docked_ligands folder in the results folder.
# Inside of that folder, there will exist a subfolder for each docked ligand (yes, it could be hundreds or thousands of folders)
# that is named docked_<ligand_name>
# that will contain the .log and .pdbqt files for all proteins docked with that dedicated ligand.
# The folder should also contain a scores file that details the specific ligand's final ranking in each of the scores files,
# for each protein-comparisonLigand pair's score file and in each comparisonLigand_RMS score file.

# If the ligand is actually one of several models of a greater ligand file, then store each model as a subfolder of a folder with the 
# name of the overall as follows:
#docked_ligands
#   <ligand_name>
#       <ligand_name1>
#           ...
#       <ligand_name2>
#           ...
#       <ligand_name3>
#           ...
#       <ligand_name4>
#           ...
#       ...

# ex:
#docked_ligands
#   AAEAMN
#       AAEAMN_model01
#           ...
#       AAEAMN_model02
#           ...
#       AAEAMN_model03
#           ...
#       AAEAMN_model04
#           ...
#       AAEAMN_model05
#           ...
#       ...

# Make a similar system where like it makes a folder of folders whenever a model is encountered 
# that would otherwise get a dedicated folder for all other systems here as well

# also, have the program remove the weird .xaa thing and any other like extraneous file 
# extensions maybe when naming folders and final files i think



# <ligand#_name>_scores.txt should contain the following:
#<comparisonLigand1_name>:
#   RMS: <RMSOfScoresWithcomparisonLigand1> - #<rankOutOfTotalRMSWithcomparisonLigand1>
#   <protein1_name>: <score-comparison1ScoreWithProtein1> - #<rankOfLigand#OutOfTotalScoresWithcomparisonLigand1AndProtein1>
#   <protein2_name>: <score-comparison1ScoreWithProtein2> - #<rankOfLigand#OutOfTotalScoresWithcomparisonLigand1AndProtein2>
#   <protein3_name>: <score-comparison1ScoreWithProtein3> - #<rankOfLigand#OutOfTotalScoresWithcomparisonLigand1AndProtein3>
#   <protein4_name>: <score-comparison1ScoreWithProtein4> - #<rankOfLigand#OutOfTotalScoresWithcomparisonLigand1AndProtein4>
#   <protein5_name>: <score-comparison1ScoreWithProtein5> - #<rankOfLigand#OutOfTotalScoresWithcomparisonLigand1AndProtein5>
#<comparisonLigand2_name>:
#   RMS: <RMSOfScoresWithcomparisonLigand2> - #<rankOutOfTotalRMSWithcomparisonLigand2>
#   <protein1_name>: <score-comparison2ScoreWithProtein1> - #<rankOfLigand#OutOfTotalScoresWithcomparisonLigand2AndProtein1>


# So, the file structure will be as follows:
# dock
#   results
#       scores
#           scores_<comparisonLigand1_name>_RMS.txt
#           scores_<comparisonLigand2_name>_RMS.txt
#           scores_<comparisonLigand1_name>
#               scores_<comparisonLigand1_name>_in_<protein1_name>.txt
#               scores_<comparisonLigand1_name>_in_<protein2_name>.txt
#               scores_<comparisonLigand1_name>_in_<protein3_name>.txt
#               scores_<comparisonLigand1_name>_in_<protein4_name>.txt
#               scores_<comparisonLigand1_name>_in_<protein5_name>.txt
#           scores_<comparisonLigand2_name>
#               scores_<comparisonLigand2_name>_in_<protein1_name>.txt
#               scores_<comparisonLigand2_name>_in_<protein2_name>.txt
#               scores_<comparisonLigand2_name>_in_<protein3_name>.txt
#               scores_<comparisonLigand2_name>_in_<protein4_name>.txt
#               scores_<comparisonLigand2_name>_in_<protein5_name>.txt
#           scores_proteins
#               scores_<protein1_name>
#                   scores_<protein1_name>_unsorted.txt
#                   scores_<protein1_name>_sorted.txt
#               scores_<protein2_name>
#                   scores_<protein2_name>_unsorted.txt
#                   scores_<protein2_name>_sorted.txt
#               scores_<protein3_name>
#                   scores_<protein3_name>_unsorted.txt
#                   scores_<protein3_name>_sorted.txt
#               scores_<protein4_name>
#                   scores_<protein4_name>_unsorted.txt
#                   scores_<protein4_name>_sorted.txt
#               scores_<protein5_name>
#                   scores_<protein5_name>_unsorted.txt
#                   scores_<protein5_name>_sorted.txt
#               
#       docked_ligands
#           docked_<ligand1_name>
#               <ligand1_name>_scores.txt
#               <ligand1_name>_in_<protein1_name>.log
#               <ligand1_name>_in_<protein1_name>.pdbqt
#               <ligand1_name>_in_<protein2_name>.log
#               <ligand1_name>_in_<protein2_name>.pdbqt
#               <ligand1_name>_in_<protein3_name>.log
#               <ligand1_name>_in_<protein3_name>.pdbqt
#               <ligand1_name>_in_<protein4_name>.log
#               <ligand1_name>_in_<protein4_name>.pdbqt
#               <ligand1_name>_in_<protein5_name>.log
#               <ligand1_name>_in_<protein5_name>.pdbqt
#           docked_<ligand2_name>
#               <ligand2_name>_scores.txt
#               <ligand2_name>_in_<protein1_name>.log
#               <ligand2_name>_in_<protein1_name>.pdbqt
#               <ligand2_name>_in_<protein2_name>.log
#               <ligand2_name>_in_<protein2_name>.pdbqt
#               <ligand2_name>_in_<protein3_name>.log
#               <ligand2_name>_in_<protein3_name>.pdbqt
#               <ligand2_name>_in_<protein4_name>.log
#               <ligand2_name>_in_<protein4_name>.pdbqt
#               <ligand2_name>_in_<protein5_name>.log
#               <ligand2_name>_in_<protein5_name>.pdbqt
#           ...

# ex:
#dock
#   results
#       scores
#           scores_ligand_abscisic_acid_RMS.txt
#           scores_ligand_biotin_RMS.txt
#           scores_ligand_abscisic_acid
#               scores_ligand_abscisic_acid_in_1stp
#               scores_ligand_abscisic_acid_in_3k3k
#           scores_ligand_biotin
#               scores_ligand_biotin_in_1stp
#               scores_ligand_biotin_in_3k3k
#       docked_ligands
#           docked_biotin
#               biotin_scores.txt
#               biotin_in_1stp.log
#               biotin_in_1stp.pdbqt
#               biotin_in_3k3k.log
#               biotin_in_3k3k.pdbqt
#           docked_AAEAMN
#               docked_AAEAMN_model01
#                   AAEAMN_model01_scores.txt
#                   AAEAMN_model01_in_1stp.log
#                   AAEAMN_model01_in_1stp.pdbqt
#                   AAEAMN_model01_in_3k3k.log
#                   AAEAMN_model01_in_3k3k.pdbqt
#               docked_AAEAMN_model02
#                   AAEAMN_model02_scores.txt
#                   AAEAMN_model02_in_1stp.log
#                   AAEAMN_model02_in_1stp.pdbqt
#                   AAEAMN_model02_in_3k3k.log
#                   AAEAMN_model02_in_3k3k.pdbqt
#               docked_AAEAMN_model03
#                   AAEAMN_model03_scores.txt
#                   AAEAMN_model03_in_1stp.log
#                   AAEAMN_model03_in_1stp.pdbqt
#                   AAEAMN_model03_in_3k3k.log
#                   AAEAMN_model03_in_3k3k.pdbqt
#               docked_AAEAMN_model04
#                   AAEAMN_model04_scores.txt
#                   AAEAMN_model04_in_1stp.log
#                   AAEAMN_model04_in_1stp.pdbqt
#                   AAEAMN_model04_in_3k3k.log
#                   AAEAMN_model04_in_3k3k.pdbqt
#               docked_AAEAMN_model05
#                   AAEAMN_model05_scores.txt
#                   AAEAMN_model05_in_1stp.log
#                   AAEAMN_model05_in_1stp.pdbqt
#                   AAEAMN_model05_in_3k3k.log
#                   AAEAMN_model05_in_3k3k.pdbqt
#               


# Upon program completion, the top 5 ligands from each scores_<comparisonLigand#_name>_RMS.txt file will be displayed to the screen.


