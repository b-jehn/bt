import subprocess
T = "255"

for m in range(26,32):
    with open("stats_m" + str(m) + "_t" + T + ".csv", "w") as outfile:
        subprocess.run(
            [
                "../../rshash/build/test/ministats",
                "-i", "../../rshash/DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa",
                "-k", "31",
                "-m", str(m),
                "-t", T,
            ],
            stdout=outfile,
            check=True,
        )