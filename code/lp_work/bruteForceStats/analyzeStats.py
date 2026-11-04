from pathlib import Path
import csv
from itertools import product
import matplotlib.pyplot as plt
import gzip
import random
import shutil

DATA_DIRECTORY = Path(f"../../rshash/DNA_datasets")
N_SAMPLES = 100000 
K = 31

GENOMES = ["Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs", 
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]
MINIMIZER_LENGTHS = []
for m1 in range(15,26,2):
    for m2 in range(m1,26,2):
        for m3 in range(m2,26,2):
            MINIMIZER_LENGTHS.append([m1,m2,m3])

THRESHOLDS = []
for t1 in range(5,8):
    for t2 in range(t1,8):
        for t3 in range(t2,8):
            THRESHOLDS.append([2**t1,2**t2,2**t3])

NO_LEVELS = len(MINIMIZER_LENGTHS[0])

def getParameterDirectory(minimizer_lengths: list[int], thresholds: list[int]) -> Path:
    no_levels = len(minimizer_lengths)
    m_part = "m="
    t_part = "t="
    for i, (m, t) in enumerate(zip(minimizer_lengths, thresholds), start=1):
        m_part += f"{m}"
        t_part += f"{t}"
        if(i != no_levels):
            m_part += f","
            t_part += f","
    return Path(f"{m_part}_{t_part}")


def plotOperationTimes(genome_stats: Path, genome_plots: Path):
    print(f"plotting operation times for: {genome_stats}")
    operation_names = [
        "R_contains",
        "S_select",
        "O_access",
        "T_compare",
        "H_access",
    ]
    values_by_operation = {name: {} for name in operation_names}

    for minimizer_lengths, thresholds in product(MINIMIZER_LENGTHS, THRESHOLDS):
        parameters_path = (genome_stats / getParameterDirectory(minimizer_lengths, thresholds))
        if not parameters_path.is_dir():
            continue

        stats_path = parameters_path / "times_per_operation.csv"
        with open(stats_path, newline="") as f:
            rows = csv.DictReader(f, skipinitialspace=True)

            for row in rows:
                level = int(row["level"])

                for name in operation_names:
                    values_by_operation[name].setdefault(level, []).append(
                        float(row[name])
                    )

    for name in operation_names:
        levels = sorted(values_by_operation[name])
        data = [values_by_operation[name][level] for level in levels]

        fig, ax = plt.subplots()
        ax.boxplot(data, tick_labels=levels)
        ax.set_title(f"{name} time distribution by level")
        ax.set_xlabel("Level")
        ax.set_ylabel("ns")
        fig.tight_layout()

        plt.savefig(str(genome_plots / f"{name} time distribution by level"), dpi=300)
        plt.close()
        

def main() -> None:
    stats_directory = Path(f"allStats")  #directory for stats of all genomes
    if(not stats_directory.is_dir()):
        return
    
    plot_path = Path("plots")
    plot_path.mkdir(exist_ok=True)
    
    for genome in GENOMES:
        genome_stats = stats_directory / f"{genome}"  #directory for the stats of a specific genome
        if(not genome_stats.is_dir()):               
            continue
        
        genome_plots = plot_path / f"{genome}"
        if(genome_plots.is_dir()):
            continue
        
        genome_plots.mkdir(exist_ok=True)
        plotOperationTimes(genome_stats, genome_plots)


if __name__ == '__main__':
    main()


