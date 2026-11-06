from pathlib import Path
import csv
from itertools import product
import matplotlib.pyplot as plt

import numpy as np
from scipy.stats import gaussian_kde

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


def boxPlotOperationTimes(genome_stats: Path, genome_plots: Path):
    print(f"plotting operation times box plots for: {genome_stats}")
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
        
def plotDistributionsOfOperationTimes(genome_stats: Path, genome_plots: Path):
    print(f"plotting operation time distributions for: {genome_stats}")
    genome_plots.mkdir(parents=True, exist_ok=True)

    operation_names = [
        "R_contains",
        "S_select",
        "O_access",
        "T_compare",
        "H_access",
    ]
    values_by_operation = {name: {} for name in operation_names}

    for minimizer_lengths, thresholds in product(MINIMIZER_LENGTHS, THRESHOLDS):
        parameters_path = (
            genome_stats / getParameterDirectory(minimizer_lengths, thresholds)
        )
        if not parameters_path.is_dir():
            continue

        stats_path = parameters_path / "times_per_operation.csv"
        if not stats_path.is_file():
            continue

        with open(stats_path, newline="") as f:
            rows = csv.DictReader(f, skipinitialspace=True)

            for row in rows:
                level = int(row["level"])

                for name in operation_names:
                    values_by_operation[name].setdefault(level, []).append(
                        float(row[name])
                    )

    for name in operation_names:
        fig, ax = plt.subplots()
        plotted_any = False

        for level in sorted(values_by_operation[name]):
            values = np.asarray(values_by_operation[name][level], dtype=float)

            # KDE needs at least two non-identical values.
            if len(values) < 2 or np.ptp(values) == 0:
                continue

            kde = gaussian_kde(values)
            x = np.linspace(values.min(), values.max(), 300)
            ax.plot(x, kde(x), label=f"Level {level}")
            plotted_any = True

        ax.set_title(f"{name} time distribution by level")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Density")
        if plotted_any:
            ax.legend(title="Level")

        fig.tight_layout()
        fig.savefig(
            genome_plots / f"{name}_time_distribution_by_level.png",
            dpi=300,
        )
        plt.close(fig)
        

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
        
        box_plots = plot_path / "box_plots" / f"{genome}"
        box_plots.mkdir(parents=True ,exist_ok=True)
        boxPlotOperationTimes(genome_stats, box_plots)
        
        distribution_plots = plot_path / "distribution_plots" / f"{genome}"
        distribution_plots.mkdir(parents=True, exist_ok=True)
        plotDistributionsOfOperationTimes(genome_stats, distribution_plots)


if __name__ == '__main__':
    main()


