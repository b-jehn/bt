import subprocess
from pathlib import Path
import gzip
import secrets
from time import perf_counter


def generate_random_31mers(output_path: Path, n: int):
    """Write n random 31-base DNA reads to path/query.fastq.gz."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    bases = "ACGT"

    with gzip.open(output_path, "wt", encoding="ascii") as out:
        for i in range(n):
            sequence = "".join(secrets.choice(bases) for _ in range(31))
            out.write(f"@read_{i + 1}\n{sequence}\n+\n{'I' * 31}\n")

GENOME = "Gadus_morhua.gadMor3.0.dna.toplevel"
DATA_DIRECTORY = Path(f"../../rshash/DNA_datasets")
N_QUERIES = 100000 

query_path = (DATA_DIRECTORY / "query.fastq.gz")
generate_random_31mers(query_path, N_QUERIES)

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