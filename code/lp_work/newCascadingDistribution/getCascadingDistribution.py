import subprocess
import os

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from time import perf_counter

from concurrent.futures import ThreadPoolExecutor

N_THREADS = 1  # Set the number of concurrent workers.

GENOME = "Gadus_morhua.gadMor3.0.dna.toplevel"
genome_directory = Path(f"./{GENOME}")
if(not genome_directory.is_dir()):
    print(f"creating directory for '{GENOME}' stats")
    genome_directory.mkdir()

print(f"computing cascading distribution on GENOME: {GENOME}")

possible_minimizers = [m for m in range(15,20)]
possible_thresholds = [2**t for t in range(7)]

def get_parameter_info(minimizer_lengths, thresholds, level):
    parts = []

    for l in range(1,level):
        parts.append(f"(m{l}={minimizer_lengths[l]})")
        parts.append(f"(t{l}={thresholds[l]})")
        
    parts.append(f"m{level}={minimizer_lengths[level]}")
    parts.append(f"t{level}={thresholds[level]}")

    return f"{', '.join(parts)}\n"
    
        
def call_subprocess(level_path, minimizer_lengths, thresholds,level, previous_threshold):
    logfile_path = level_path / "build.log"
    if(thresholds[level] == possible_thresholds[0]):
        if(level==1):
            infile = f"../../rshash/DNA_datasets/{GENOME}.fa"
        else:
            infile = str(level_path.parents[1] / "freq_skmers.fasta")
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
                "--m1", str(minimizer_lengths[level]),
                "--t1", str(thresholds[level]),
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


def stats_recursive(level_path, minimizer_lengths, thresholds, level, maxlevel, top_level_minimizer):
    if(level>maxlevel):
        return
    
    minimizers = (
        [top_level_minimizer]
        if level == 1 and top_level_minimizer is not None
        else possible_minimizers
    )
    
    previous_threshold = 0
    minimizer_lengths.append(0)
    thresholds.append(0)
    for minimizer_lengths[level] in minimizers:
        level_path = level_path / f"m{level}={minimizer_lengths[level]}"
        if(not level_path.is_dir()): level_path.mkdir()
        for thresholds[level] in possible_thresholds:
            level_path = level_path / f"t{level}={thresholds[level]}"
            if(not level_path.is_dir()):
                level_path.mkdir()
                print(get_parameter_info(minimizer_lengths, thresholds,level) + " - computing stats")
                call_subprocess(level_path, 
                                minimizer_lengths, 
                                thresholds,
                                level,
                                previous_threshold)
            #else: print(get_parameter_info(minimizer_lengths, thresholds,level) + " - stats already exist")
            stats_recursive(level_path, 
                            minimizer_lengths, 
                            thresholds,
                            level+1,
                            maxlevel,
                            top_level_minimizer)
            level_path = level_path.parent
            #delete freq_skmers.fasta from previous computation, which is no longer needed
            previous_output = (level_path
                                / f"t{level}={previous_threshold}"
                                / "freq_skmers.fasta")
            previous_output.unlink(missing_ok=True)
            previous_threshold = thresholds[level]

            
        #delete freq_skmers.fasta from previous computation, which is no longer needed
        previous_output = (level_path 
                            / f"t{level}={previous_threshold}"
                            / "freq_skmers.fasta")
        previous_output.unlink(missing_ok=True)
        level_path = level_path.parent

def run_top_level_minimizer(minimizer):
    print(f"running thread for m1={minimizer}")
    stats_recursive(
        level_path=genome_directory,
        minimizer_lengths=[0],
        thresholds=[0],
        level=1,
        maxlevel=3,
        top_level_minimizer=minimizer
    )

if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=N_THREADS) as executor:
        list(executor.map(run_top_level_minimizer, possible_minimizers))