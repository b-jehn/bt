import subprocess
import os

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

genome = "Gadus_morhua.gadMor3.0.dna.toplevel" #also has to be manually changed in build.cpp - as of now
levels = 3
#minimizer_sets = [[14,14,14], [5, 18, 18], , [14, 16, 18], [5,5,5], [7,7,7], [9,9,9], [11,11,11], [13,13,13], [15,15,15], [17,17,17],[18, 18, 18],[19,19,19]]
minimizerLengths_sets = [[15, m, m] for m in range(10,20)]# + [[m,m+2,m+4] for m in range(5,16)] + [[m,m+4,m+8] for m in range(5,12)]
threshold_sets = [[64, 64, 64]]
#t = [32, 32, 32]

for t in threshold_sets:
    for m in minimizerLengths_sets:
        ###############
        # build rshash and compute stats
        ###############
        stats_directory = Path(f"./{genome}/m={m[0]},{m[1]},{m[2]}_t={t[0]},{t[1]},{t[2]}")
        if(stats_directory.is_dir()):
            print(f"stats for m: {m} and t: {t} already exist")
            continue
        stats_directory.mkdir()
        
        print(f"\n\ncomputing leveled stats for m: {m} and t: {t}")
        logfile_path = stats_directory / f"build.log"
        with logfile_path.open("w") as logfile:
            subprocess.run(
                [
                    "../../rshash/build/source/rshash",
                    "build",
                    "-i", f"../../rshash/DNA_datasets/{genome}.fa",
                    "-d", f"../../rshash/DNA_datasets/{genome}.rshash",
                    "-k", "31",
                    "--m1", str(m[0]),
                    "--m2", str(m[1]),
                    "--m3", str(m[2]),
                    "--t1", str(t[0]),
                    "--t2", str(t[1]),
                    "--t3", str(t[2]),
                    "--level", str(levels)
                ],
                stdout=logfile,
                stderr=subprocess.STDOUT,
                check=True,
            )

        for type in ["distinct", "cumulated"]:
            ###############
            # plot p_M stats
            ###############
            plt.figure(figsize=(8, 6))

            for level in range(1, 4):
                stats_path = stats_directory / f"{type}Stats_level{level}.csv"

                stats = pd.read_csv(stats_path)

                # Use the first column as x and the second column as y
                x = stats.iloc[:, 0]
                y = stats.iloc[:, 1]


                plt.scatter(x, y, label=f"Level {level}\nm={m[level-1]}; t= {t[level-1]}", s=35)

            plt.xlabel(stats.columns[0])
            plt.ylabel(stats.columns[1])
            plt.title(f"{type} minimizer frequencies per Level")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.ylim(0, 100)
            plt.tight_layout()

            #save the plot
            plt.savefig(stats_directory / f"{type}_p_M_per_level.png", dpi=300)
            #plt.show()
            plt.close()

            ###############
            # plot p_K stats
            ###############
            plt.figure(figsize=(8, 6))

            for level in range(1, 4):
                stats_path = stats_directory / f"{type}Stats_level{level}.csv"

                stats = pd.read_csv(stats_path)

                # Use the first column as x and the second column as y
                x = stats.iloc[:, 0]
                y = stats.iloc[:, 2]


                plt.scatter(x, y, label=f"Level {level}\nm={m[level-1]}; t= {t[level-1]}", s=35)

            plt.xlabel(stats.columns[0])
            plt.ylabel(stats.columns[2])
            plt.title(f"{type} k_mers covered per Level")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.ylim(0, 100)
            plt.tight_layout()

            #save the plot
            plt.savefig(stats_directory / f"{type}_p_K_per_level.png", dpi=300)
            #plt.show()
            plt.close()