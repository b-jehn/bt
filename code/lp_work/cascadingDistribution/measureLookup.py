import subprocess
from pathlib import Path
import gzip
import random

GENOME = "Gadus_morhua.gadMor3.0.dna.toplevel"
DATA_DIRECTORY = Path(f"../../rshash/DNA_datasets")
N_QUERIES = 100000 
K = 31

def get_random_Kmers(source_path: Path, output_path: Path, n: int) -> None:
    """write n K-base DNA reads from random positions in source to output_path"""
    # Read sequence lines, ignoring FASTA headers and whitespace.
    parts = []
    with source_path.open() as source:
        for line in source:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append("".join(line.split()))

    sequence = "".join(parts).upper()
    if len(sequence) < K:
        raise ValueError("Source sequence must contain at least K bases")

    with gzip.open(output_path, "wt", encoding="ascii") as out:
        for i in range(n):
            start = random.randrange(len(sequence) - K + 1)
            out.write(f"@read_{i+1}\n{sequence[start:start + K]}\n+\n{'I' * K}\n")


query_path = (DATA_DIRECTORY / "query.fastq.gz")
get_random_Kmers((DATA_DIRECTORY / f"{GENOME}.fa"), query_path, N_QUERIES)

logfile_path = Path("./lookup_log.txt")

with logfile_path.open("w") as logfile:
    subprocess.run(
        [
            "../../rshash/build/source/rshash",
            "lookup",
            "-d", str(DATA_DIRECTORY / f"{GENOME}.rshash"),
            "-q", str(query_path),
        ],
        stdout=logfile,
        check=True,
    )