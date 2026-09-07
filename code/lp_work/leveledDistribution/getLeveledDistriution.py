import subprocess
import os

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

genome = "cod" #also has to be manually changed in build.cpp - as of now
levels = 3
m1, m2, m3 = 5, 18, 18, 
t1, t2, t3 = 255, 255, 255

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
        "--m1", str(m1),
        "--m2", str(m2),
        "--m3", str(m3),
        "--t1", str(t1),
        "--t2", str(t2),
        "--t3", str(t3),
        "--level", str(levels)
    ],
    check=True,
)

###############
# plot the stats
###############
# directory containing level1.csv, level2.csv, and level3.csv
stats_directory = Path(f"./{genome}/leveledMinimizerDistribution_m={m1},{m2},{m3}_t={t1},{t2},{t3}")
plt.figure(figsize=(8, 6))

for level in range(1, 4):
    stats_path = stats_directory / f"level{level}.csv"

    stats = pd.read_csv(stats_path)

    # Use the first column as x and the second column as y
    x = stats.iloc[:, 0]
    y = stats.iloc[:, 1]

    plt.scatter(x, y, label=f"Level {level}", s=35)

plt.xlabel("threshold")
plt.ylabel("percentage of kmers with minimizers of frequency <= threshold")
plt.title("Kmer Coverage by Threshold and Level")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

#save the plot
plt.savefig(stats_directory / "kmer_coverage_by_level.png", dpi=300)
plt.show()