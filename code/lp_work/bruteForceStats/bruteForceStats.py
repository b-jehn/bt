import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import gzip
import random
import shutil

DATA_DIRECTORY = Path(f"../../rshash/DNA_datasets")
N_SAMPLES = 100000 
K = 31

GENOMES = ["Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]

MINIMIZER_LENGTHS = []
for m1 in range(15,26,2):
    for m2 in range(m1,26,2):
        for m3 in range(m2,26,2):
            MINIMIZER_LENGTHS.append([m1,m2,m3])

THRESHOLDS = []
for t1 in range(5,8):
    for t2 in range(t1,8):
        for t3 in range(t2,8):
            THRESHOLDS.append([2**t1,2**t2,2**t3])


def get_random_Kmers(source_path: Path, output_path: Path, no_samples: int) -> None:
    """write n K-base DNA reads from random positions in source to output_path"""
    # Read sequence lines, ignoring FASTA headers and whitespace.
    parts = []
    with source_path.open() as source:
        for line in source:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append("".join(line.split()))

    sequence = "".join(parts)
    if len(sequence) < K:
        raise ValueError("Source sequence must contain at least K bases")

    with gzip.open(output_path, "wt", encoding="ascii") as out:
        for i in range(no_samples):
            skmer = parts[random.randrange(0, len(parts)-1, 1)]
            start = random.randrange(len(skmer) - K + 1)
            out.write(f"@read_{i+1}\n{skmer[start:start + K]}\n+\n{'I' * K}\n")

def getStats(genome: str, genome_path: Path, query_path: Path, minimizer_lengths: list[int], thresholds: list[int]) -> None:
    print(f"computating stats for minimizer lengths: {minimizer_lengths} and thresholds: {thresholds}")
    no_levels = len(minimizer_lengths)
    m_part = "m="
    t_part = "t="
    for i, (m, t) in enumerate(zip(minimizer_lengths, thresholds), start=1):
        m_part += f"{m}"
        t_part += f"{t}"
        if(i != no_levels):
            m_part += f","
            t_part += f","
            
    stats_path =  genome_path / f"{m_part}_{t_part}"    #directory for a specific set of level parameter
    if(not stats_path.is_dir()):
        stats_path.mkdir()

    args = [
        "../../rshash/build/source/rshash",
        "build",
        "-i", str(DATA_DIRECTORY / f"{genome}.fa"),
        "-d", str(DATA_DIRECTORY / f"{genome}_{m_part}_{t_part}.rshash"),
        "-k", str(K),
        "--stats_directory", str(stats_path),
        "--no_samples", str(N_SAMPLES),
        "--level", str(no_levels)
    ]
    for i, (m, t) in enumerate(zip(minimizer_lengths, thresholds), start=1):
        args.extend([f"--m{i}", str(m), f"--t{i}", str(t)])

    logfile_path = stats_path / "build.log"
    with logfile_path.open("w") as logfile:
        subprocess.run(
            args,
            stdout=logfile,
            stderr=subprocess.STDOUT,
            check=True,
        ) 

    ############################
    #measure average lookup time
    ############################
    logfile_path = stats_path / "lookup.log"
    with logfile_path.open("w") as logfile:
        subprocess.run(
            [
                "../../rshash/build/source/rshash",
                "lookup",
                "-d", str(DATA_DIRECTORY / f"{genome}_{m_part}_{t_part}.rshash"),
                "-q", str(query_path),
            ],
            stdout=logfile,
            check=True,
        )
    (DATA_DIRECTORY / f"{genome}_{m_part}_{t_part}.rshash").unlink() #delete rshash structures which are no longer needed

def main(n_threads: int) -> None:
    stats_directory = Path(f"stats_directory")  #directory for stats of all genomes
    if(not stats_directory.is_dir()):
            stats_directory.mkdir()
            
    for genome in GENOMES:
        print(f"start computation for threads of: {genome}")
        genome_path = stats_directory / f"{genome}"  #directory for the stats of a specific genome
        if(not genome_path.is_dir()):               
                genome_path.mkdir()

        #create query file with N_SAMPLES many random kmers that appear in the genome
        #needed for computing the average lookup time
        query_path = (DATA_DIRECTORY / f"{genome}_query.fastq.gz")
        get_random_Kmers((DATA_DIRECTORY / f"{genome}.fa"), query_path, N_SAMPLES)  
        
        level_parameters = [
            (minimizer_length, threshold)
            for minimizer_length in MINIMIZER_LENGTHS
            for threshold in THRESHOLDS
        ]
        
        # Generative AI use Disclaimer:
        #   the following code for threading the computation was generated using GPT-6 Luna
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = [
                executor.submit(
                    getStats,
                    genome,
                    genome_path,
                    query_path,
                    minimizer_length,
                    threshold,
                )
                for minimizer_length, threshold in level_parameters
            ]

            # Re-raise any exception raised by a worker.
            for future in futures:
                future.result()
        
        query_path.unlink() #delete query file, as it is no longer needed

    #zip stats_directory
    shutil.make_archive("allStats", 'zip', root_dir=str(stats_directory)) 
    #remove unziped stats
    if stats_directory.is_dir():
        for genome_path in stats_directory.iterdir():
            if genome_path.is_dir():
                for level_parameters in genome_path.iterdir():
                    if level_parameters.is_dir():
                        for file_path in level_parameters.iterdir():
                            if file_path.is_file():
                                file_path.unlink()
                                
                        level_parameters.rmdir()  # Removes the directory only if it is now empty
                        
                genome_path.rmdir()  # Removes the directory only if it is now empty
                
        stats_directory.rmdir()  # Removes the directory only if it is now empty
    
if __name__ == '__main__':
    main(int(sys.argv[1]))

