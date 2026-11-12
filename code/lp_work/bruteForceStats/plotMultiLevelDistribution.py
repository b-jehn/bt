from pathlib import Path
import csv
from itertools import product
import matplotlib.pyplot as plt

import numpy as np
from scipy.stats import gaussian_kde

DATA_DIRECTORY = Path(f"../../rshash/DNA_datasets")
N_SAMPLES = 100000 
K = 31

GENOMES = ["Falco_tinnunculus.FalTin1.0.dna.toplevel",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Gadus_morhua.gadMor3.0.dna.toplevel",
            "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]

SEQUENCES = GENOMES

MINIMIZER_LENGTHS = [16,16,16]
THRESHOLDS = [32,32,32]

NO_LEVELS = 3

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
    return f"{m_part}_{t_part}"



def plot_multilevel_distributions(stats_path: Path, result_path: Path):
    print(f"plotting multilevel distribution at: {result_path}")
    distributions = {}

    for level in range(1, NO_LEVELS + 1):
        stats_file = stats_path / f"{level}kmer_minimizer_distribution.csv"

        frequencies = []
        perc_m = []
        perc_k = []
        with open(stats_file, "r", newline="") as f:
            reader = csv.reader(f)
            next(reader, None)  # Skip header
            for row in reader:
                frequencies.append(float(row[0]))
                perc_m.append(float(row[1]))
                perc_k.append(float(row[2]))
                
        cdf_m = [sum(perc_m[:i + 1]) for i in range(len(perc_m))]
        cdf_k = [sum(perc_k[:i + 1]) for i in range(len(perc_k))]
        distributions[level] = {
            "frequency": frequencies,
            "perc_M": perc_m,
            "perc_K": perc_k,
            "cdf_m": cdf_m,
            "cdf_k": cdf_k
        }
        
    with open((result_path / f"F_M_{stats_path.parts[-2]}.csv"), "w") as data_file:
        headers = ["frequency", "F_M_1", "F_M_2", "F_M_3"]
        writer = csv.writer(data_file)
        writer.writerow(headers)
        for f,(F_M_1, F_M_2, F_M_3) in enumerate(zip(distributions[1]["cdf_m"],distributions[2]["cdf_m"],distributions[3]["cdf_m"]), start=1):
            writer.writerow([f,F_M_1, F_M_2, F_M_3])
            

    """
    # Probability distributions
    fig, axes = plt.subplots(2, 2, figsize=(14, 5))

    for level, data in distributions.items():
        axes[0].plot(data["frequency"], data["perc_M"], marker=".", label=f"Level {level}")
        axes[1].plot(data["frequency"], data["perc_K"], marker=".", label=f"Level {level}")

    axes[0].set(title="Probability distributions (M)", xlabel="Value", ylabel="Probability")
    axes[1].set(title="Probability distributions (K)", xlabel="Value", ylabel="Probability")

    for ax in axes:
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(result_path / "multilevel_distributions.png", dpi=300)
    plt.close(fig)
"""
    # Cumulative probability distributions
    #fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig, axes = plt.subplots(1, 1)
    
    for level, data in distributions.items():
        axes.plot(data["frequency"], data["cdf_m"], linestyle='dotted', color='b', marker=f"${level}$", label=f"Level {level}")
        #axes[1].plot(data["frequency"], data["cdf_k"], marker=".", label=f"Level {level}")

    axes.set(title=f"{stats_path.parts[-2]}\nCumulative distributions (M) for {stats_path.parts[-1]}", xlabel="Frequency", ylabel="Cumulative Probability")
    #axes[1].set(title="Cumulative distributions (K)", xlabel="Value", ylabel="Cumulative probability")

    #for ax in axes:
    #    ax.legend()
    #    ax.grid(True, alpha=0.3)
    
    axes.legend()
    axes.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(result_path / "multilevel_cumulative_distributions.png", dpi=300)
    plt.close(fig)
    
    

def main() -> None:
    stats_directory = Path(f"allStats")  #directory for stats of all genomes
    if(not stats_directory.is_dir()):
        return
    
    plot_path = Path("plots")
    plot_path.mkdir(exist_ok=True)
    
    for genome in SEQUENCES:
        genome_stats_path = stats_directory / f"{genome}" / getParameterDirectory(MINIMIZER_LENGTHS, THRESHOLDS)  #directory for the stats of a specific genome
        if(not genome_stats_path.is_dir()):
            print(f"no stats at: {genome_stats_path}")
            continue
        
        multilevel_distributions_path = plot_path / "multilevel_distributions" / f"{genome}"
        multilevel_distributions_path.mkdir(parents=True ,exist_ok=True)
        plot_multilevel_distributions(genome_stats_path, multilevel_distributions_path)


if __name__ == '__main__':
    main()


