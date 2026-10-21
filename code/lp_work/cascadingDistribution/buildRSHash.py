import subprocess

GENOME = "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs"
MINIMIZER_LENGTHS = [15]
THRESHOLDS = [128]

args = [
    "../../rshash/build/source/rshash",
    "build",
    "-i", f"../../rshash/DNA_datasets/{GENOME}.fa",
    "-d", f"../../rshash/DNA_datasets/{GENOME}.rshash",
    "-k", "31",
]

for i, (m, t) in enumerate(zip(MINIMIZER_LENGTHS, THRESHOLDS), start=1):
    args.extend([f"--m{i}", str(m), f"--t{i}", str(t)])

args.extend(["--level", str(len(MINIMIZER_LENGTHS))])

subprocess.run(
    args,
    stderr=subprocess.STDOUT,
    check=True,
)