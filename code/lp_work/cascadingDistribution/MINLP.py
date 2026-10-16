from pyscipopt import Model, quicksum, quickprod
from datetime import datetime
from time import perf_counter
import math
import csv
import sys
import os

from itertools import product
from dataclasses import dataclass, field
from pathlib import Path

#temporary!! TODO: measure and pass as part of cascading stats
K = 31  #kmer length that is/was used to gain minimizer_stats
LOAD_FACTOR = 0.75  #of the hash table

SEQUENCE_NAME = "Gadus_morhua.gadMor3.0.dna.toplevel"
OPTIMIZATION_GOAL = "t_lookup"   #"t_lookup" or "bpk_total"
MAX_LOOKUP_TIME = math.inf
MAX_BPK_TOTAL = math.inf

LEVELS = 3
MINIMIZER_LENGTHS = [15] #[m for m in range(15,20)]
THRESHOLDS = [64,128] #[2**i for i in range(0,7)]

logfile_path = f'MINLP_logs/{SEQUENCE_NAME}.txt'
sequence_path = f'../../rshash/DNA_datasets/{SEQUENCE_NAME}.fa'
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

@dataclass
class Stats:
    n_kmers: int = field(default_factory=int)
    p_K: float = field(default_factory=float)
    f_avg: float = field(default_factory=float)
    next_level: dict[(int, int), "Stats | None"] = field(default_factory=dict)  #links to next stats item depending on (m,t)


def readStatsFromFiles(stats_directory):
    n_kmers = int((stats_directory / "n_kmers.txt").read_text().strip())

    reader = csv.reader((stats_directory / "p_K_distinct.csv").open())
    next(reader)  # Skip the header
    p_K_distinct = [0]+[float(row[0]) for row in reader if row]
    
    p_K = [0]
    f_avg = [0]
    weighted_frequency_sum = 0
    for f, p in enumerate(p_K_distinct[1:]):
        p_K.append(p_K[-1] + p) 
        weighted_frequency_sum += p * (f+1)
        f_avg.append(weighted_frequency_sum/p_K[-1])
        
    return n_kmers, p_K, f_avg

def fill_cascading_stats(root_stats, level_path, level, maxLevel):
    if(level>maxLevel):
        return 
    previous_threshold = 0
    for m in MINIMIZER_LENGTHS:
        level_path = level_path / f"m{level}={m}"
        n_kmers, p_K, f_avg = readStatsFromFiles(level_path)
        for t in THRESHOLDS:
            current_stats = Stats(n_kmers, p_K[t], f_avg[t], next_level={(m,t): None for m,t in product(MINIMIZER_LENGTHS, THRESHOLDS)})
            root_stats.next_level[(m,t)] = current_stats
            level_path = level_path / f"t{level}={t}"
            fill_cascading_stats(current_stats, level_path, level+1, maxLevel)
            level_path = level_path.parent
        level_path = level_path.parent

cascading_stats = Stats(n_kmers=0, p_K=0, f_avg=0, next_level={(m,t): None for m,t in product(MINIMIZER_LENGTHS, THRESHOLDS)}) 
stats_directory = Path(f"./{SEQUENCE_NAME}")
fill_cascading_stats(cascading_stats, stats_directory, 1, 3)

def getCascadingStats(minimizer_lengths,thresholds):
    stats = cascading_stats
    for m, t in zip(minimizer_lengths, thresholds):
        stats = stats.next_level[(m,t)]
    return stats.n_kmers, stats.p_K, stats.f_avg

#============================#
#Set up SCIP Model
#============================#
print("setting up scip model")
scip = Model("Level")  # scip name is optional

#============================#
#VARIABLES OR MEASUREMENTS??
#============================#
t_access_R = scip.addVar("t_access_R")
t_access_O = scip.addVar("t_access_O")
t_access_H = scip.addVar("t_access_H")
t_rank_R = scip.addVar("t_rank_R")
t_select_S = scip.addVar("t_select_S")
t_compare_T = scip.addVar("t_compare_T")

bpk_T = 2   #bits per k-mer for storing the text T
#============================#
#VARIABLES
#============================#
im_vars = {                                             #binary indicators for m on each level l
    (m, l): scip.addVar(name=f"im_{m}_{l}", vtype="B")
    for m in MINIMIZER_LENGTHS
    for l in range(LEVELS)
}
it_vars = {                                             #binary indicators for t on each level l
    (t, l): scip.addVar(name=f"it_{t}_{l}", vtype="B")
    for t in THRESHOLDS
    for l in range(LEVELS)
}

m_vars = [scip.addVar(f"m_{l}") for l in range(LEVELS)]         #minimizer length on each level
t_vars = [scip.addVar(f"t_{l}") for l in range(LEVELS)]         #threshold on each level

n_kmers_vars = [scip.addVar(f"n_kmers_{l}") for l in range(LEVELS)] #number of kmers on level l for each level l
p_K_vars = [scip.addVar(f"p_K_{l}") for l in range(LEVELS)]         #percentage of kmers that have minimizers of length m with a frequency <= t for each level
f_avg_vars = [scip.addVar(f"f_avg_{l}") for l in range(LEVELS)]     #average frequency of minimizers of length m with a frequency <= t for each level

t_level_vars = [scip.addVar(f"t_level_{l}") for l in range(LEVELS)] #time spent for finding k-mer with an infrequent minimizer on level l
t_lookup = scip.addVar("t_lookup")  #total lookup time

bpk_R_vars = [scip.addVar(f"bpk_R_{l}") for l in range(LEVELS)] #bits per kmer in R on each level
bpk_S_vars = [scip.addVar(f"bpk_S_{l}") for l in range(LEVELS)] #bits per kmer in S on each level
bpk_O_vars = [scip.addVar(f"bpk_O_{l}") for l in range(LEVELS)] #bits per kmer in O on each level
bpk_H = scip.addVar("bpk_H")            #bits per kmer in H
bpk_total = scip.addVar("bpk_total")    #total bits per k-mer

#============================#
#OBJECTIVE
#============================#
if(OPTIMIZATION_GOAL == "t_lookup"):
    scip.setObjective(t_lookup, sense="minimize")   #t_lookup objective
elif(OPTIMIZATION_GOAL == "bpk_total"):
    scip.setObjective(bpk_total, sense="minimize")   #space objective
else:
    sys.exit(f"OPTIMIZATION_GOAL should be 't_lookup' or 'bpk_total' but is: {OPTIMIZATION_GOAL}")
#============================#
#CONSTRAINTS
#============================#
#constraints regarding the structure operation times
scip.addCons(t_access_R - t_access_R_empirical == 0)
scip.addCons(t_access_O - t_access_O_empirical == 0)
scip.addCons(t_access_H - t_access_H_empirical == 0)
scip.addCons(t_rank_R - t_rank_R_empirical == 0)
scip.addCons(t_select_S - t_select_S_empirical == 0)
scip.addCons(t_compare_T - t_compare_T_empirical == 0)

for l in range(LEVELS):
    scip.addCons(quicksum(im_vars[(m,l)] for m in MINIMIZER_LENGTHS) == 1, name=f"im_mutually_exclusive_{l}")   #mutually exclusive indicators for m
    scip.addCons(quicksum(it_vars[(t,l)] for t in THRESHOLDS) == 1, name=f"it__mutually_exclusive_{l}")         #mutually exclusive indicators for t
    scip.addCons(m_vars[l] - quicksum(im_vars[(m,l)] * m for m in MINIMIZER_LENGTHS) == 0, name=f"get_m_{l}")   #determine m
    scip.addCons(t_vars[l] - quicksum(it_vars[(t,l)] * t for t in THRESHOLDS) == 0, name=f"get_t_{l}")          #determine t
    
    #determine p_K_l and f_avg_l
    n_kmers_terms = []
    p_K_sum_terms = []
    f_avg_sum_terms = []
    for minimizer_lengths in product(MINIMIZER_LENGTHS, repeat=l+1):
        m_indicator = 1
        for l, m in enumerate(minimizer_lengths):
            m_indicator *= im_vars[(m, l)]
        for thresholds in product(THRESHOLDS, repeat=l+1):
            t_indicator = 1
            for i, t in enumerate(thresholds):
                t_indicator *= it_vars[(t, l)]
            n_kmers, p_K_empiric, f_avg_empiric = getCascadingStats(minimizer_lengths, thresholds)
            n_kmers_terms.append(m_indicator * t_indicator * n_kmers)
            p_K_sum_terms.append(m_indicator * t_indicator * p_K_empiric)
            f_avg_sum_terms.append(m_indicator * t_indicator * f_avg_empiric)
    
    scip.addCons(n_kmers_vars[l] - quicksum(n_kmers_terms) == 0)
    scip.addCons(p_K_vars[l] - quicksum(p_K_sum_terms) == 0)
    scip.addCons(f_avg_vars[l] - quicksum(f_avg_sum_terms) == 0)
    
    scip.addCons(t_level_vars[l] - (t_rank_R + 2*t_select_S + (f_avg_vars[l]+1)/2*(t_access_O+t_compare_T)) == 0)

    scip.addCons(bpk_R_vars[l] == 1)
    #TODO: use real calculation for R
    #n_DIM = p_K_vars[l]*n_kmers_vars[l] #number of distinct infrequent minimizers on level
    #scip.addCons(bpk_R_vars[l] - 
    #            (1/n_kmers_vars[l] * 
    #            (n_DIM
    #                * math.ceil(math.log(4**m_vars[l]/n_DIM, 2))
    #                + 3*n_DIM)
    #            ) == 0)
    
    scip.addCons(bpk_S_vars[l] - 
                (quickprod(1-p_K_vars[i] for i in range(l))
                * p_K_vars[l] * 2/(K-m_vars[l]+2)
                ) == 0)
    
    #bits per k-mer in O
    scip.addCons(bpk_O_vars[l] - (bpk_S_vars[l] * log_sequence_length) == 0)


#recursivly build the t_lookup constraint level by level
remaining_time = t_access_H
for l in reversed(range(LEVELS)):
    remaining_time = t_access_R + p_K_vars[l] * t_level_vars[l] + (1 - p_K_vars[l]) * remaining_time

scip.addCons(t_lookup - remaining_time == 0, name="t_lookup")
scip.addCons(t_lookup <= MAX_LOOKUP_TIME)

#bits per k-mer in H
scip.addCons(bpk_H - (quickprod(1-p_K_vars[l] for l in range(LEVELS)) * (64+log_sequence_length)/LOAD_FACTOR) == 0)
scip.addCons(bpk_total - (bpk_T + quicksum(bpk_R_vars[l] + bpk_S_vars[l] + bpk_O_vars[l] for l in range(LEVELS)) + bpk_H) == 0)
scip.addCons(bpk_total <= MAX_BPK_TOTAL)

#============================#
#SOLVE
#============================#
print("start solving")
start_time = perf_counter()
scip.optimize()
elapsed_seconds = perf_counter() - start_time
sol = scip.getBestSol()

#============================#
#OUTPUT SOLUTION
#============================#
print(f"\nsolved MINLP for sequence: {SEQUENCE_NAME}")
print(f"with optimization goal: {OPTIMIZATION_GOAL}, MAX_LOOKUP_TIME = {MAX_LOOKUP_TIME} and MAX_BPK_TOTAL = {MAX_BPK_TOTAL}")

print(f"for {LEVELS} levels and k = {K}")
print(f"considered minimizer lengths: {MINIMIZER_LENGTHS}")
print(f"considered thresholds: {THRESHOLDS}")

print("\nfound solution:")
for l in range(LEVELS):
    print(f"m_{l}: {format(sol[m_vars[l]])}")
print()
for l in range(LEVELS):
    print(f"t_{l}: {format(sol[t_vars[l]])}") 
print() 
for l in range(LEVELS):  
    print(f"p_K_{l}: {format(sol[p_K_vars[l]])}")
print() 
for l in range(LEVELS):  
    print(f"f_avg_{l}: {format(sol[f_avg_vars[l]])}")
    
print("\n modeled lookup time:")
print(f"t_lookup: {format(sol[t_lookup])}")

print("\n modeled space consumption:")
print(f"bpk_total: {format(sol[bpk_total])}")
print("more specifically:")
print(f"bpk_T: {bpk_T}")
for l in range(LEVELS):
    print(f"bpk_R_{l}: {format(sol[bpk_R_vars[l]])}")
    print(f"bpk_S_{l}: {format(sol[bpk_S_vars[l]])}")
    print(f"bpk_O_{l}: {format(sol[bpk_O_vars[l]])}")
print(f"bpk_H: {format(sol[bpk_H])}")


#============================#
#WRITE SOLUTION TO LOG
#============================#
with open(logfile_path, "w") as logfile:
    logfile.write(f"{datetime.now()}---------------------------------------------\n")
    logfile.write(f"solved MINLP for sequence: {SEQUENCE_NAME}\n")
    logfile.write(f"with optimization goal: {OPTIMIZATION_GOAL}, MAX_LOOKUP_TIME = {MAX_LOOKUP_TIME} and MAX_BPK_TOTAL = {MAX_BPK_TOTAL}\n")
    logfile.write(f"for {LEVELS} levels and k = {K}\nconsidered minimizer lengths: {MINIMIZER_LENGTHS}\nconsidered thresholds: {THRESHOLDS}\n")

    logfile.write(f"\ntotal time to solve MINLP: {elapsed_seconds:.6f}\n")
    logfile.write("\nfound solution:\n")
    for l in range(LEVELS):
        logfile.write(f"m_{l}: {format(sol[m_vars[l]])}\n")
    logfile.write("\n")
    for l in range(LEVELS):
        logfile.write(f"t_{l}: {format(sol[t_vars[l]])}\n") 
    logfile.write("\n") 
    for l in range(LEVELS):  
        logfile.write(f"p_K_{l}: {format(sol[p_K_vars[l]])}\n")
    logfile.write("\n")
    for l in range(LEVELS):  
        logfile.write(f"f_avg_{l}: {format(sol[f_avg_vars[l]])}\n")
    logfile.write("\n")
    
    logfile.write("\n modeled lookup time:\n")
    logfile.write(f"t_lookup: {format(sol[t_lookup])}\n")
    
    logfile.write("\n modeled space consumption:\n")
    logfile.write(f"bpk_total: {format(sol[bpk_total])}\n")
    logfile.write("more specifically:\n")
    logfile.write(f"bpk_T: {bpk_T}\n")
    for l in range(LEVELS):
        logfile.write(f"bpk_R_{l}: {format(sol[bpk_S_vars[l]])}\n")
        logfile.write(f"bpk_S_{l}: {format(sol[bpk_S_vars[l]])}\n")
        logfile.write(f"bpk_O_{l}: {format(sol[bpk_O_vars[l]])}\n")

    logfile.write(f"bpk_H: {format(sol[bpk_H])}\n\n\n")



    



    