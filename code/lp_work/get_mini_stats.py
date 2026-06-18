import subprocess

for m in range(15,25):
    with open("stats_m" + str(m) + ".csv", "w") as outfile:
        subprocess.run(
            [
                "../rshash/build/test/ministats",
                "-i", "../rshash/DNA_datasets/Gadus_morhua.gadMor3.0.dna.toplevel.fa",
                "-k", "31",
                "-m", str(m),
                "-t", "64",
            ],
            stdout=outfile,
            check=True,
        )