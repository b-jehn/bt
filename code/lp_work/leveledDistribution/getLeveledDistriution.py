import subprocess
import os

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

genome = "cod" #also has to be manually changed in build.cpp - as of now
levels = 3
minimizer_sets = [[14,14,14], [5, 18, 18], [18, 18, 18], [14, 16, 18], [5,5,5], [7,7,7], [9,9,9], [11,11,11], [13,13,13], [15,15,15], [17,17,17]]
#minimizer_sets = [[19,19,19]]
t = [64, 64, 64]

for m in minimizer_sets:
    ###############
    # build rshash and compute stats
    ###############
    subprocess.run(
        [
            "../../rshash/build/source/rshash",
            "build",
            "-i", "../../rshash/DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa",
            "-d", "../../rshash/DNA_datasets/Gadus_morhua.rshash",
            "-k", "31",
            "--m1", str(m[0]),
            "--m2", str(m[1]),
            "--m3", str(m[2]),
            "--t1", str(t[0]),
            "--t2", str(t[1]),
            "--t3", str(t[2]),
            "--level", str(levels)
        ],
        check=True,
    )

    ###############
    # plot p_M stats
    ###############
    # directory containing level1.csv, level2.csv, and level3.csv
    stats_directory = Path(f"./{genome}/leveledMinimizerDistribution_m={m[0]},{m[1]},{m[2]}_t={t[0]},{t[1]},{t[2]}")
    plt.figure(figsize=(8, 6))

    for level in range(1, 4):
        stats_path = stats_directory / f"level{level}.csv"

        stats = pd.read_csv(stats_path)

        # Use the first column as x and the second column as y
        x = stats.iloc[:, 0]
        y = stats.iloc[:, 1]


        plt.scatter(x, y, label=f"Level {level}\nm={m[level-1]}; t= {t[level-1]}", s=35)

    plt.xlabel("threshold")
    plt.ylabel("percentage of minimizers with frequency <= threshold")
    plt.title("Minimizer Frequencies per Level")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 100)
    plt.tight_layout()

    #save the plot
    plt.savefig(stats_directory / "minimizer_frequencies_per_level.png", dpi=300)
    #plt.show()

    ###############
    # plot p_K stats
    ###############
    # directory containing level1.csv, level2.csv, and level3.csv
    stats_directory = Path(f"./{genome}/leveledMinimizerDistribution_m={m[0]},{m[1]},{m[2]}_t={t[0]},{t[1]},{t[2]}")
    plt.figure(figsize=(8, 6))

    for level in range(1, 4):
        stats_path = stats_directory / f"level{level}.csv"

        stats = pd.read_csv(stats_path)

        # Use the first column as x and the second column as y
        x = stats.iloc[:, 0]
        y = stats.iloc[:, 2]


        plt.scatter(x, y, label=f"Level {level}\nm={m[level-1]}; t= {t[level-1]}", s=35)

    plt.xlabel("threshold")
    plt.ylabel("percentage of kmers with minimizers of frequency <= threshold")
    plt.title("Covered k_mers per Level")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 100)
    plt.tight_layout()

    #save the plot
    plt.savefig(stats_directory / "kmer_coverage_by_level.png", dpi=300)
    #plt.show()
    plt.close()