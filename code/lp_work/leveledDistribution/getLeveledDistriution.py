import subprocess
import os

m = 15
T = "255"
genome = "cod"

# checking if the directory demo_folder 
# exist or not.
if not os.path.exists(genome):
    
    # if the demo_folder directory is not present 
    # then create it.
    os.makedirs(genome)

outfile_path = f"{genome}/stats_m" + str(m) + "_t" + T + ".csv"
print(f"writing {outfile_path}")
with open(outfile_path, "w") as outfile:
    subprocess.run(
        [
            "../../rshash/build/test/leveledStats",
            "-i", f"../../rshash/DNA_datasets/{genome}.fa",
            "-k", "31",
            "-m", str(m),
            "-t", T,
        ],
        stdout=outfile,
        check=True,
    )