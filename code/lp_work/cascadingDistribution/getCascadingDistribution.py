import subprocess
import os

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

"""
Idea for stats, that take different thresholds into consideration: p_K on the first level is not dependent on t1
only the freq_kmers.fa file is affected and thus affects further levels
- 
try out the following approach: create a file freq_kmers_t=i.fa for every i in thresholds
then use that file for calculating stats on the next level
-
or think a little further: technically k-mers should only have one possible frequency when building with the same seed
so it should be possible to 
first build for the smallest threshold whilst
    - computing n_kmers
    - saving the stats and
    - store freq_kmers_t=t1.fa
then use freq_kmers_t=t1.fa as input for the computation using the next threshold
    - compute its stats
        - BUT use n_kmers from fist step for percentage calculation!! TODO!!
        - pK for f <= the previous threshould should be 0 as those k-mers should not be part of the input
        - they can thus be simply copied from the previous computation
    - store freq_kmers_t=t2.fa
repeat ...
and finally freq_kmers_t=n-1.fa as input for the computation using highest threshold
    - compute its stats
        - BUT use n_kmers from fist step for percentage calculation!! TODO!!
        - pK for f <= the previous threshould should be 0 as those k-mers should not be part of the input
        - they can thus be simply copied from the previous computation
    - store freq_kmers_t=tn.fa
    
    
further idea: cut of branches of the computation, if freq k_mers is to large
- this would mean that those parameters would never be considered optimal
- this would also mean that the information about cut of branches and thus uninteresting 
    parameters would have to be saved and used by MINLP.py so that those -non existent- values
    are never used when solving the MINLP
""" 

GENOME = "freq_skmers"
genome_directory = Path(f"./{GENOME}")
if(not genome_directory.is_dir()):
    print(f"creating directory for '{GENOME}' stats")
    genome_directory.mkdir()

print(f"computing cascading distribution on GENOME: {GENOME}")

possible_minimizers = [m for m in range(15,20)]
minimizer_lengths = [int, int, int]

possible_thresholds = [2**t for t in range(7)]
thresholds = [int,int,int]

stats_directories = [Path, Path, Path]

def create_stats_files(stats_directory):
    stats_path = stats_directory / "distinctStats.csv"
    with stats_path.open("w", encoding ="utf-8") as stats:
        stats.write("f, ")
        stats.write("p_M: percentage of minimizers with frequency == f on the level, ")
        stats.write("p_k: percentage of kmers having minimizers of frequency == f on the level\n")
        
    stats_path = stats_directory / "cumulatedStats.csv"
    with stats_path.open("w", encoding ="utf-8") as stats:
        stats.write("f, ")
        stats.write("p_M: percentage of minimizers with frequency <= f on the level, ")
        stats.write("p_k: percentage of kmers having minimizers of frequency <= f on the level\n")
    
        
def call_subprocess(level):
    if(level==0):
        m_string = f"m1={minimizer_lengths[0]}"
    elif(level==1):
        m_string = f"(m1={minimizer_lengths[0]}), m2={minimizer_lengths[1]}"
    else:
        m_string = f"(m1={minimizer_lengths[0]}, m2={minimizer_lengths[1]}), m3={minimizer_lengths[2]}"
    
    if(stats_directories[level].is_dir()):
        print(f"{m_string} - stats already exist")
        return
    
    print(f"{m_string} - computing stats")
    stats_directories[level].mkdir()
    
    logfile_path = stats_directories[level] / "build.log"

    if(thresholds[level] == possible_thresholds[0]):
        create_stats_files(stats_directories[level])
        if(level==0):
            infile = f"../../rshash/DNA_datasets/{GENOME}.fa"
        else:
            infile = str(stats_directories[level-1] / f"freq_skmers_t={possible_thresholds[-1]}.fasta")
    else:
        infile = str(stats_directories[level] / f"freq_skmers_t={thresholds[level]/2}.fasta")
    
    with logfile_path.open("w") as logfile:
        subprocess.run(
            [
                "../../rshash/build/source/rshash",
                "build",
                "-i", infile,
                "-d", f"../../rshash/DNA_datasets/{GENOME}.rshash",
                "-k", "31",
                "--m1", str(minimizer_lengths[level]),
                "--t1", str(thresholds[level]),
                "--level", str(level+1),
                "--stats_directory", str(stats_directories[level]),
            ],
            stdout=logfile,
            stderr=subprocess.STDOUT,
            check=True,
        )

for minimizer_lengths[0] in possible_minimizers:
    stats_directories[0] = Path(f"./{GENOME}/m1={minimizer_lengths[0]}")
    for thresholds[0] in possible_thresholds:
        call_subprocess(0)

    for minimizer_lengths[1] in possible_minimizers:
        stats_directories[1] = stats_directories[0] / Path(f"m2={minimizer_lengths[1]}")
        call_subprocess(1)
        
        for minimizer_lengths[2] in possible_minimizers:
            stats_directories[2] = stats_directories[1] / Path(f"m3={minimizer_lengths[2]}")
            call_subprocess(2)