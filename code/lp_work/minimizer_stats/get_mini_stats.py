import subprocess
import os

T = "255"
genome = "frequent_skmers_level1"

# checking if the directory demo_folder 
# exist or not.
if not os.path.exists(genome):
    
    # if the demo_folder directory is not present 
    # then create it.
    os.makedirs(genome)

for m in range(18,19):
    print(f"writing {genome}/stats_m" + str(m) + "_t" + T + ".csv")
    with open(f"{genome}/stats_m" + str(m) + "_t" + T + ".csv", "w") as outfile:
        subprocess.run(
            [
                "../../rshash/build/test/ministats",
                "-i", f"../../rshash/DNA_datasets/{genome}.fa",
                "-k", "31",
                "-m", str(m),
                "-t", T,
            ],
            stdout=outfile,
            check=True,
        )