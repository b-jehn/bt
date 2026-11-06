import math
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

SEQUENCE_NAME = "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs"
OPTIMIZATION_GOAL = "t_lookup"   #"t_lookup" or "bpk_total"
MAX_LOOKUP_TIME = math.inf
MAX_BPK_TOTAL = math.inf

MINIMIZER_LENGTHS_CONSIDERED = [m for m in range(15,26,2)]
THRESHOLDS_CONSIDERED = [2**i for i in range(5,8)]  

NO_LEVELS = 3
NON_DECREASING_PARAMETERS = True

#TODO: tmp, for testing
forcedMinimizerLength = []
forcedThresholds = []

logfile_path = f'MINLP_logs/{SEQUENCE_NAME}.txt'
sequence_path = f'../../../rshash/DNA_datasets/{SEQUENCE_NAME}.fa'
sequence_length = 0
if os.path.isfile(sequence_path):
    sequence_length = os.path.getsize(sequence_path)
else:
    sys.exit(f"{sequence_path} is not a file")
log_sequence_length = math.log(sequence_length, 2)


#TODO: model or measure!
t_access_R_empirical = 30
t_access_O_empirical = 10
t_access_H_empirical = 20
t_rank_R_empirical = 35
t_select_S_empirical = 35
t_compare_T_empirical = 10

data_path = Path(f"{SEQUENCE_NAME}/data")
data_path.mkdir(exist_ok=True, parents=True)
error_file_path = data_path / "squarred_errors.json"

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
        
setOfMinimizerLengths = [[15, 19, 19]]
setOfThresholds = [[32, 32, 128]]

for forcedMinimizerLengths, forcedThresholds in product(setOfMinimizerLengths, setOfThresholds):
    print(f"running MINLP with forced minimizer lengths: {forcedMinimizerLengths} and thresholds: {forcedThresholds}")
    modeledSpace, modeledTimes, info, solution, elapsed_seconds = MINLP(
        optimization_goal="t_lookup",
        max_lookup_time=math.inf,
        max_bpk_total=math.inf,
        log_sequence_length=log_sequence_length,
        no_levels=NO_LEVELS, 
        minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED, 
        thresholdsConsidered=THRESHOLDS_CONSIDERED, 
        distribution_stats=genome_stats,
        load_factor=LOAD_FACTOR,
        nonDecreasingParameters=NON_DECREASING_PARAMETERS,
        forceParameters=[forcedMinimizerLengths,forcedThresholds],
        hideOutput=False)

    measuredSpace, measuredTimes = getMeasuredStats(genome_stats_directory, 
                                                    forcedMinimizerLengths,
                                                    forcedThresholds)
    
    for line in info:
        print(line, end='')
    
    print(f"minimizer_lengths: {forcedMinimizerLengths}, thresholds. {forcedThresholds}")
    print("stats:\nmeasured space\mmodeled space\n")
    for f in fields(modeledSpace):
        name = f.name
        #print(f"{name}:\n{getattr(measuredSpace, name)}\n{getattr(modeledSpace, name)}\n")
        
    print("stats:\nmeasured times\nmodeled times:")
    for f in fields(modeledTimes):
        name = f.name
        print(f"{name}:\n{getattr(measuredTimes, name)}\n{getattr(modeledTimes, name)}\n")