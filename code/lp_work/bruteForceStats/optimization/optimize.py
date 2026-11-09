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

GENOMES = ["Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]

SEQUENCE_NAME = "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"
OPTIMIZATION_GOAL = "t_lookup"   #"t_lookup" or "bpk_total"
MAX_LOOKUP_TIME = np.inf
MAX_BPK_TOTAL = np.inf

MINIMIZER_LENGTHS_CONSIDERED = [m for m in range(15,26,2)]
THRESHOLDS_CONSIDERED = [2**i for i in range(5,8)]  

NO_LEVELS = 3
NON_DECREASING_PARAMETERS = True

sequence_path = f'../../../rshash/DNA_datasets/{SEQUENCE_NAME}.fa'
sequence_length = 0
if os.path.isfile(sequence_path):
    sequence_length = os.path.getsize(sequence_path)
else:
    sys.exit(f"{sequence_path} is not a file")
log_sequence_length = np.log2(sequence_length)

data_path = Path(f"{SEQUENCE_NAME}/MINLP_logs")
data_path.mkdir(exist_ok=True, parents=True)

genome_stats = Distribution_Stats(n_kmers=0, n_minimizers=0, F_M=0, F_K=0, f_avg=0, next_level={(m,t): None for m,t in product(MINIMIZER_LENGTHS_CONSIDERED, THRESHOLDS_CONSIDERED)}) 
genome_stats_directory = Path(f"../allStats/{SEQUENCE_NAME}")
fill_cascading_stats(root_stats=genome_stats, 
                    genome_stats_path=genome_stats_directory, 
                    level=1, 
                    maxLevel=NO_LEVELS,
                    currentMinimizerLengths=[MINIMIZER_LENGTHS_CONSIDERED[-1]] * NO_LEVELS,
                    currentThresholds=[THRESHOLDS_CONSIDERED[-1]] * NO_LEVELS,
                    minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED,
                    thresholdsConsidered=THRESHOLDS_CONSIDERED,
                    nonDecreasingParameters=True)


        
print(f"optimizing MINLP")
modeledSpace, modeledTimes, info, solution, elapsed_seconds = MINLP(
    optimization_goal="t_lookup",
    max_lookup_time=MAX_LOOKUP_TIME,
    max_bpk_total=MAX_BPK_TOTAL,
    log_sequence_length=log_sequence_length,
    no_levels=NO_LEVELS, 
    minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED, 
    thresholdsConsidered=THRESHOLDS_CONSIDERED, 
    distribution_stats=genome_stats,
    load_factor=LOAD_FACTOR,
    nonDecreasingParameters=NON_DECREASING_PARAMETERS)

#============================#
#OUTPUT INFO
#============================#
timestamp = str(datetime.now())
info.insert(0, f"{timestamp}---------------------------------------------\n")
info.insert(1, f"solved MINLP for sequence: {SEQUENCE_NAME}\n")
for line in info:
    print(line, end='')
#============================#
#WRITE SOLUTION TO LOG
#============================#
logfile_path = data_path / f"optimization.log"
with open(logfile_path, "a") as logfile:    #TODO: change to mode a - append, once consistent
    for line in info:
        logfile.write(line)