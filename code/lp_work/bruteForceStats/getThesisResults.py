from pathlib import Path
import csv
from itertools import product
import matplotlib.pyplot as plt

import numpy as np
from scipy.stats import gaussian_kde

from optimization.cascadingStats import *

DATA_DIRECTORY = Path(f"../../rshash/DNA_datasets")
N_SAMPLES = 100000 
K = 31

GENOMES = ["Falco_tinnunculus.FalTin1.0.dna.toplevel",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Gadus_morhua.gadMor3.0.dna.toplevel",
            "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel",
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

NO_LEVELS = 3

def getParameterDirectory(minimizer_lengths: list[int], thresholds: list[int]) -> Path:
    no_levels = len(minimizer_lengths)
    m_part = "m="
    t_part = "t="
    for i, (m, t) in enumerate(zip(minimizer_lengths, thresholds), start=1):
        m_part += f"{m}"
        t_part += f"{t}"
        if(i != no_levels):
            m_part += f","
            t_part += f","
    return Path(f"{m_part}_{t_part}")


def getAllOperationTimes(stats_directory: Path, thesis_results_path: Path):
    print(f"getting information needed for latex boxplots")
    #allStats = {genome: {
    #    field_name: {
    #        values: [],
    #        median: int,
    #    } for field_name in fields(operationTimes_Stats)
    #} for genome in GENOMES         
    #}
    
    with open((thesis_results_path / "allOperationTimes.csv"), "w") as csvfile:
        headers = [f"Rcontains{l}" for l in range(NO_LEVELS)]
        headers += [f"Sselect{l}" for l in range(NO_LEVELS)]
        headers += [f"Oaccess{l}" for l in range(NO_LEVELS)]
        headers += [f"Tcompare{l}" for l in range(NO_LEVELS)]
        headers += ["Haccess", "LookupTime"]
        
        writer = csv.writer(csvfile)
        writer.writerow([field_name_to_tex[header] for header in headers])
        for genome in GENOMES:
            print(f"writing stats for: {genome}")
            genome_stats_path = stats_directory / f"{genome}"  #directory for the stats of a specific genome
            if(not genome_stats_path.is_dir()):               
                continue
            
            for minimizers, thresholds in zip(MINIMIZER_LENGTHS, THRESHOLDS):
                _, times = getObservedStats(genome_stats_path, minimizers,thresholds)
                values = []
                for stat_field in fields(times):
                    name = stat_field.name
                    attr = getattr(times, name)
                    if isinstance(attr, list):
                        for level in range(NO_LEVELS):
                            values.append(attr[level])
                    else:
                        values.append(attr)
                writer.writerow(values)
        

def main() -> None:
    thesis_results_path = Path("thesis_results")
    thesis_results_path.mkdir(parents=True, exist_ok=True)
    
    stats_directory = Path("allStats")  #directory for stats of all genomes
    if(not stats_directory.is_dir()):
        return
    getAllOperationTimes(stats_directory, thesis_results_path)

if __name__ == '__main__':
    main()


