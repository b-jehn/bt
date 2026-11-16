import numpy as np
import json

import sys
import os
import matplotlib.pyplot as plt
from itertools import product
from datetime import datetime

from pathlib import Path
from dataclasses import fields

from cascadingStats import *
from MINLP import MINLP

#temporary!! TODO: measure and pass as part of cascading stats
K = 31  #kmer length that is/was used to gain minimizer_stats
LOAD_FACTOR = 0.75  #of the hash table

GENOMES = ["Gadus_morhua.gadMor3.0.dna.toplevel",
            "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]

MINIMIZER_LENGTHS_CONSIDERED = [m for m in range(15,26,2)]
THRESHOLDS_CONSIDERED = [2**i for i in range(5,8)]  

NO_LEVELS = 3
NON_DECREASING_PARAMETERS = True

def optimize(log_sequence_length: float, genome_stats: Distribution_Stats,  logPath: Path, csvPath: Path, opt_goal="t_lookup" , max_lookupTime=np.inf, max_bpkTotal=np.inf):
            
    print(f"optimizing MINLP with\nopt_goal: {opt_goal}\nmax lookupTime: {max_lookupTime}\nmax bpk_total: {max_bpkTotal}")
    modeledSpace, modeledTimes, info, solution, elapsed_seconds, optimal_parameters = MINLP(
        optimization_goal=opt_goal,
        max_lookup_time=max_lookupTime,
        max_bpk_total=max_bpkTotal,
        log_sequence_length=log_sequence_length,
        no_levels=NO_LEVELS, 
        minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED, 
        thresholdsConsidered=THRESHOLDS_CONSIDERED, 
        distribution_stats=genome_stats,
        load_factor=LOAD_FACTOR,
        nonDecreasingParameters=NON_DECREASING_PARAMETERS)

    timestamp = str(datetime.now())
    info.insert(0, f"{timestamp}---------------------------------------------\n")
    #============================#
    #OUTPUT INFO
    #============================#
    #for line in info:
    #    print(line, end='')
    #============================#
    #WRITE SOLUTION TO LOG
    #============================#
    with open(logPath, "a") as logfile:    #TODO: change to mode a - append, once consistent
        for line in info:
            logfile.write(line)
    
    print(f"appending hashSetSizes csv file")
    with open(csvPath, "a") as f:
        writer = csv.writer(f)
        values = [opt_goal, max_bpkTotal, max_lookupTime]
        values += [int(round(optimal_parameters[f"m{l}"])) for l in range(1, NO_LEVELS + 1)]
        values += [int(round(optimal_parameters[f"t{l}"]))for l in range(1, NO_LEVELS + 1)]
        values += [modeledSpace.bpkTotal, modeledTimes.LookupTime]
        writer.writerow(values)
        
    return modeledSpace.bpkTotal, modeledTimes.LookupTime
            
def main():
    for sequence in GENOMES[:1]:
        sequence_path = f'../../../rshash/DNA_datasets/{sequence}.fa'
        if not os.path.isfile(sequence_path):
            sys.exit(f"{sequence_path} is not a file")
                
        log_sequence_length = np.log2(os.path.getsize(sequence_path))

        genome_stats = Distribution_Stats(n_kmers=0, n_minimizers=0, F_M=0, F_K=0, f_avg=0, next_level={(m,t): None for m,t in product(MINIMIZER_LENGTHS_CONSIDERED, THRESHOLDS_CONSIDERED)}) 
        fill_cascading_stats(root_stats=genome_stats, 
                            genome_stats_path=Path(f"../allStats/{sequence}"), 
                            level=1, 
                            maxLevel=NO_LEVELS,
                            currentMinimizerLengths=[MINIMIZER_LENGTHS_CONSIDERED[-1]] * NO_LEVELS,
                            currentThresholds=[THRESHOLDS_CONSIDERED[-1]] * NO_LEVELS,
                            minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED,
                            thresholdsConsidered=THRESHOLDS_CONSIDERED,
                            nonDecreasingParameters=True)
        
        result_path = Path(f"./optimizationResults")
        result_path.mkdir(parents=True, exist_ok=True)
        
        logPath = result_path / f"{sequence}/optimization.log"
        csvPath = result_path / f"optimizationResults_{sequence}.csv"
        
        # initialise files / clear previous contend
        #with open(logPath, 'w') as file:
        #    pass
        #with open(csvPath, "w") as f:
        #    headers = ["optGoal", "maxBPK", "maxLookupTime", "m1", "m2", "m3", "t1", "t2", "t3","modeledBPK","modeledLookupTime"]
        #    writer = csv.writer(f)
        #    writer.writerow(headers)
        
        #_, min_lookupTime = optimize(sequence, logPath, csvPath, opt_goal="t_lookup")
        
        #min_bpk, _ = optimize(sequence, logPath, csvPath, opt_goal="bpk_total")
        min_lookupTime = 300.0
        min_bpk = 15.0
        
        optimize(log_sequence_length, genome_stats, logPath, csvPath, opt_goal="t_lookup", max_bpkTotal=1.5*min_bpk)
        optimize(log_sequence_length, genome_stats, logPath, csvPath, opt_goal="bpk_total", max_lookupTime=1.5*min_lookupTime)
        
if __name__ == '__main__':
    main()