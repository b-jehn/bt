import subprocess
import csv

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from time import perf_counter

from concurrent.futures import ThreadPoolExecutor

N_THREADS = 1  # Set the number of concurrent workers.

GENOME = "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs"
MAX_LEVEL = 3
POSSIBLE_MINIMIZER_LENGTHS = [[m for m in range(15,16)] for _ in range(MAX_LEVEL)]
POSSIBLE_THRESHOLDS = [[2**t for t in range(7,8)] for _ in range(MAX_LEVEL)]

genome_directory = Path(f"./{GENOME}")
if(not genome_directory.is_dir()):
    print(f"creating directory for '{GENOME}' stats")
    genome_directory.mkdir()

print(f"computing cascading distribution on GENOME: {GENOME}")

def get_p_K_distinct(level_path, level):
    kmers_path = level_path / f"t{level}={POSSIBLE_THRESHOLDS[level-1][0]}" / "n_kmers.txt"
    if(not kmers_path.exists()):
        return
    n_kmers = int(kmers_path.open().readline())
    (level_path / "n_kmers.txt").write_text(str(n_kmers))
    with (level_path / "p_K_distinct.csv").open("w") as f:
        f.write("f, percentage of kmers having minimizers of frequency == f on the level\n")
        previous_threshold = 0
        for t in POSSIBLE_THRESHOLDS[level-1]:
            stats_path = level_path / f"t{level}={t}" / "infrequent_kmers.csv"
            reader = csv.reader(stats_path.open())
            for freq in range(previous_threshold+1): next(reader)              #Skip the header and all values already computed for lower thresholds
            for row in reader:
                f.write(f"{int(row[1])/n_kmers}\n")    #write distinct percentages for every row left over
            previous_threshold = t
            stats_path.unlink(missing_ok=True)  #delete file as it is no longer needed
            (level_path / f"t{level}={t}" / "n_kmers.txt").unlink(missing_ok=True)
        
def call_subprocess(level_path: Path, minimizer_length, threshold,level,previous_threshold,lastLevel=False):
    logfile_path = level_path / "build.log"
    if(lastLevel or threshold == POSSIBLE_THRESHOLDS[level-1][0]):
        if(level==1):
            infile = f"../../rshash/DNA_datasets/{GENOME}.fa"
        else:
            infile = str(level_path.parents[0 if lastLevel else 1] / "freq_skmers.fasta")
    else:
        infile = str(level_path.parent / f"t{level}={previous_threshold}" / "freq_skmers.fasta")
    
    start_time = perf_counter()
    with logfile_path.open("w") as logfile:
        subprocess.run(
            [
                "../../rshash/build/source/rshash",
                "build",
                "-i", infile,
                "-d", f"../../rshash/DNA_datasets/{GENOME}.rshash",
                "-k", "31",
                "--m1", str(minimizer_length),
                "--t1", str(threshold),
                "--level", str(level),
                "--stats_directory", str(level_path),
            ],
            stdout=logfile,
            stderr=subprocess.STDOUT,
            check=True,
        )
        elapsed_seconds = perf_counter() - start_time
        time_file = level_path / "seconds_for_computation.txt"
        time_file.write_text(f"{elapsed_seconds:.6f}\n")

"""
for the last level, the stats for every threshold can be computed directly
as it is no longer necessary to create freq_kmers.fa for every threshold
because there is no computation on higher levels
"""
def handleMaxLevel(level_path: Path, level: int, minimizers: list):
    for m in minimizers:
        level_path = level_path / f"m{level}={m}"
        if(level_path.is_dir()):
            return
        
        level_path.mkdir()
        print(str(level_path.parts[-2*level:]) + " - computing stats on last level")
        call_subprocess(level_path, 
                m, 
                POSSIBLE_THRESHOLDS[level-1][-1],
                level,
                0,
                lastLevel=True)
        #get p_K_distinct.csv from infrequent_kmers.csv
        n_kmers = int((level_path / "n_kmers.txt").open().readline())
        with (level_path / "p_K_distinct.csv").open("w") as f:
            f.write("f, percentage of kmers having minimizers of frequency == f on the level\n")
            stats_path = level_path / "infrequent_kmers.csv"
            reader = csv.reader(stats_path.open())
            next(reader)    #Skip the header
            for row in reader:
                f.write(f"{int(row[1])/n_kmers}\n")    #write distinct percentages for every row left over
            #delete files that are no longer needed
            stats_path.unlink(missing_ok=True)  
            (level_path / "freq_skmers.fasta").unlink(missing_ok=True) #TODO: dont even compute this in the first case
        level_path = level_path.parent 
        
def stats_recursive(level_path, level, maxlevel, first_level_minimizer_length):
    #needed for threading
    minimizers = (
        [first_level_minimizer_length]
        if level == 1 and first_level_minimizer_length is not None
        else POSSIBLE_MINIMIZER_LENGTHS[level-1]
    )
    
    if(level == maxlevel): #directly compute p_K instead of computing freq_kmer.fa for every possible t
        handleMaxLevel(level_path, level, minimizers)
        return
    
    previous_threshold = 0
    for m in minimizers:
        level_path = level_path / f"m{level}={m}"
        if(not level_path.is_dir()): level_path.mkdir()
        for t in POSSIBLE_THRESHOLDS[level-1]:
            level_path = level_path / f"t{level}={t}"
            if(not (level_path / "seconds_for_computation.txt").exists()):
                level_path.mkdir()
                print(str(level_path.parts[-2*level:]) + " - computing stats")
                call_subprocess(level_path, 
                                m, 
                                t,
                                level,
                                previous_threshold)
            else: print(str(level_path.parts[-2*level:]) + " - stats already exist")
            stats_recursive(level_path,
                            level+1,
                            maxlevel,
                            first_level_minimizer_length)
            level_path = level_path.parent
            #delete freq_skmers.fasta from previous computation, which is no longer needed
            previous_output = (level_path
                                / f"t{level}={previous_threshold}"
                                / "freq_skmers.fasta")
            previous_output.unlink(missing_ok=True)
            previous_threshold = t
            
        get_p_K_distinct(level_path, level)
            
        #delete freq_skmers.fasta from previous computation, which is no longer needed
        previous_output = (level_path 
                            / f"t{level}={previous_threshold}"
                            / "freq_skmers.fasta")
        previous_output.unlink(missing_ok=True)
        level_path = level_path.parent

def run_first_level_minimizer_length(minimizer):
    print(f"running thread for m1={minimizer}")
    stats_recursive(
        level_path=genome_directory,
        level=1,
        maxlevel=MAX_LEVEL,
        first_level_minimizer_length=minimizer
    )

if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        list(executor.map(run_first_level_minimizer_length, POSSIBLE_MINIMIZER_LENGTHS[0]))