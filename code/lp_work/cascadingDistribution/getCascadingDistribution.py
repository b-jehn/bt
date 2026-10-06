import subprocess
import os

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

genome = "Gadus_morhua.gadMor3.0.dna.toplevel" #also has to be manually changed in build.cpp - as of now
print(f"computing cascading distribution on genome: {genome}")

minimizer_lengths = [int, int, int]
stats_directories = [Path, Path, Path]

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
    if(level==0):
        infile = f"../../rshash/DNA_datasets/{genome}.fa"
    else:
        infile = str(stats_directories[level-1] / "freq_skmers.fasta")
    
    with logfile_path.open("w") as logfile:
        subprocess.run(
            [
                "../../rshash/build/source/rshash",
                "build",
                "-i", infile,
                "-d", f"../../rshash/DNA_datasets/{genome}.rshash",
                "-k", "31",
                "--m1", str(minimizer_lengths[level]),
                "--t1", str(64),
                "--level", str(level+1),
                "--stats_directory", str(stats_directories[level]),
            ],
            stdout=logfile,
            stderr=subprocess.STDOUT,
            check=True,
        )
    

for minimizer_lengths[0] in range(15,20):
    stats_directories[0] = Path(f"./{genome}/m1={minimizer_lengths[0]}")
    call_subprocess(0)

    for minimizer_lengths[1] in range(15,20):
        stats_directories[1] = stats_directories[0] / Path(f"m2={minimizer_lengths[1]}")
        call_subprocess(1)
        
        for minimizer_lengths[2] in range(15,20):
            stats_directories[2] = stats_directories[1] / Path(f"m3={minimizer_lengths[2]}")
            call_subprocess(2)