import os

T_MAX = "255"
MINIMIZER_LENGTHS = [m for m in range(5,25)]
genome = "homo_sapiens"

# checking if the directory demo_folder 
# exist or not.
if not os.path.exists(f"{genome}_leveled_stats"):
    # if the demo_folder directory is not present 
    # then create it.
    os.makedirs(f"{genome}_leveled_stats")

for m in MINIMIZER_LENGTHS:
    print(f"writing {genome}_leveled_stats/stats_m{m}_{T_MAX}.csv")
    with open(f"{genome}_leveled_stats/stats_m{m}_{T_MAX}.csv", "w") as outfile:
        subprocess.run(
            [
                "../../rshash/build/test/ministats",
                "-i", f"../../rshash/DNA_datasets/{genome}.fa",
                "-k", "31",
                "-m", str(m),
                "-t", T_MAX,
            ],
            stdout=outfile,
            check=True,
        )