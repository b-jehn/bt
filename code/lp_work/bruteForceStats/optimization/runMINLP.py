import math

import sys
import os
from itertools import product
from datetime import datetime

from pathlib import Path

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

#setOfMinimizerLengths = []
#for m1 in range(15,26,2):
#    for m2 in range(m1,26,2):
#        for m3 in range(m2,26,2):
#            setOfMinimizerLengths.append([m1,m2,m3])
#            
#setOfThresholds = []
#for t1 in range(5,8):
#    for t2 in range(t1,8):
#        for t3 in range(t2,8):
#            setOfThresholds.append([2**t1,2**t2,2**t3])
#
#
#for minimizerLenths, thresholds in product(setOfMinimizerLengths, setOfThresholds):
space, times, info, solution, elapsed_seconds = MINLP(
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
    forceParameters=[[15,17,19],[32,64,128]])

print(space)

#============================#
#OUTPUT INFO
#============================#
#info.insert(0, f"{datetime.now()}---------------------------------------------\n")
#info.insert(1, f"solved MINLP for sequence: {SEQUENCE_NAME}\n")
#for line in info:
#    print(line, end='')
#============================#
#WRITE SOLUTION TO LOG
#============================#
#with open(logfile_path, "w") as logfile:    #TODO: change to mode a - append, once consistent
#    for line in info:
#        logfile.write(line)