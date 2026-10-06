import os

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

genome = "Gadus_morhua.gadMor3.0.dna.toplevel" #also has to be manually changed in build.cpp - as of now
print(f"ploting distributions for")


minimizer_lengths = [0,0,0]

stats_directories = [Path, Path, Path]

x = [i for i in range(1,65)]
for minimizer_lengths[0] in range(15,20):
    stats_directories[0] = Path(f"./{genome}/m1={minimizer_lengths[0]}")
    for minimizer_lengths[1] in range(15,20):
        stats_directories[1] = stats_directories[0] / Path(f"m2={minimizer_lengths[1]}")
        for minimizer_lengths[2] in range(15,20):
            stats_directories[2] = stats_directories[1] / Path(f"m3={minimizer_lengths[2]}")
            for filename in os.listdir(stats_directories[2]):
                if filename.endswith('.png'):
                    os.remove(os.path.join(stats_directories[2], filename))
                    print(f"Deleted: {filename}")
            """
            print(f"ploting distributions for m=[{minimizer_lengths}]")
            for type in ["distinct", "cumulated"]:
                
                ###############
                # plot p_M stats
                ###############
                plt.figure(figsize=(8, 6))

                for level in range(3):
                    stats_path = stats_directories[level] / f"{type}Stats.csv"

                    stats = pd.read_csv(stats_path)

                    # Use (the first column as x) and the second column as y
                    y = stats.iloc[:, 1]

                    plt.scatter(x, y, label=f"Level {level+1}\nm={minimizer_lengths[level]}; t= {64}", s=35)

                plt.xlabel(stats.columns[0])
                plt.ylabel(stats.columns[1])
                plt.title(f"{type} minimizer frequencies per Level")
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.ylim(0, 100)
                plt.tight_layout()

                #save the plot
                plt.savefig(stats_directories[2] / f"{type}_p_M_per_level.png", dpi=300)
                #plt.show()
                plt.close()
                
                ###############
                # plot p_K stats
                ###############
                plt.figure(figsize=(8, 6))

                for level in range(3):
                    stats_path = stats_directories[level] / f"{type}Stats.csv"

                    stats = pd.read_csv(stats_path)

                    # Use (the first column as x) and the second column as y
                    y = stats.iloc[:, 2]

                    plt.scatter(x, y, label=f"Level {level+1}\nm={minimizer_lengths[level]}; t= {64}", s=35)

                plt.xlabel(stats.columns[0])
                plt.ylabel(stats.columns[2])
                plt.title(f"{type} k_mers covered per Level")
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.ylim(0, 100)
                plt.tight_layout()

                #save the plot
                plt.savefig(stats_directories[2] / f"{type}_p_K_per_level.png", dpi=300)
                #plt.show()
                plt.close()
                """
