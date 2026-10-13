from pyscipopt import Model, quicksum, quickprod
from datetime import datetime
import math
import csv
import sys
import os

from itertools import product
from dataclasses import dataclass, field

K = 31  #kmer length that is/was used to gain minimizer_stats
T = 255 #max threshold - the .csv files in ./minimizer_stats/{SEQUENCE_NAME}/ have to contain values for all t <= T
LOAD_FACTOR = 0.75  #of the hash table

SEQUENCE_NAME = "Gadus_morhua.gadMor3.0.dna.toplevel"
OPTIMIZATION_GOAL = "t_lookup"   #"t_lookup" or "bpk_total"
MAX_BPK_TOTAL = 5               #needed, when OPTIMIZATION_GOAL == "t_lookup"
LEVELS = 3
MINIMIZER_LENGTHS = [m for m in range(15,20)]
THRESHOLDS = [2**i for i in range(0,7)]

if THRESHOLDS[-1] >= T:
    sys.exit(f"maximum threshold t: {THRESHOLDS[-1]} has to be smaller than T: {T}")

logfile_path = f'logs/log_{SEQUENCE_NAME}'
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


"""
P_K_empirical[i][j]: 
the empirically measured percentage of 
    minimizers with length == MINIMIZER_LENGTHS[i] 
    which have a frequency < j
""" 
@dataclass
class Stats:
    p_K: list[float] = field(default_factory=list)
    next_level: dict[int, "Stats | None"] = field(default_factory=dict)

cascading_stats = Stats(p_K=[], next_level={m1: None for m1 in MINIMIZER_LENGTHS[0]})

for m1 in MINIMIZER_LENGTHS: 
    curr_level1_stats = Stats(p_K=[], next_level={m1: None for m2 in MINIMIZER_LENGTHS[1]})
    cascading_stats.next_level[m1] = curr_level1_stats
    reader = csv.reader(open(f'./{SEQUENCE_NAME}/m1={m1}/cumulatedStats.csv'), delimiter=',')
    for t, row in enumerate(reader):
        if(t==0):
            curr_level1_stats.p_K.append(0) #percentage of appearing minimizers with f==0 is 0
        else:
            curr_level1_stats.p_K.append(float(row[2])/100)   #convert to percentage
    for m2 in MINIMIZER_LENGTHS: 
        curr_level2_stats = Stats(p_K=[], next_level={m3: None for m3 in MINIMIZER_LENGTHS[2]})
        curr_level1_stats.next_level[m2] = curr_level2_stats
        reader = csv.reader(open(f'./{SEQUENCE_NAME}/m1={m1}/m2={m2}/cumulatedStats.csv'), delimiter=',')
        for t, row in enumerate(reader):
            if(t==0):
                curr_level2_stats.p_K.append(0) #percentage of appearing minimizers with f==0 is 0
            else:
                curr_level2_stats.p_K.append(float(row[2])/100)   #convert to percentage
        for m3 in MINIMIZER_LENGTHS: 
            curr_level3_stats = Stats(p_K=[])
            curr_level2_stats.next_level[m3] = curr_level3_stats
            reader = csv.reader(open(f'./{SEQUENCE_NAME}/m1={m1}/m2={m2}/m3={m3}/cumulatedStats.csv'), delimiter=',')
            for t, row in enumerate(reader):
                if(t==0):
                    curr_level3_stats.p_K.append(0) #percentage of appearing minimizers with f==0 is 0
                else:
                    curr_level3_stats.p_K.append(float(row[2])/100)   #convert to percentage

def getCascadingStats(minimizer_lengths,threshold):
    stats = cascading_stats
    for m in minimizer_lengths:
        stats = stats.next_level[m]
    return stats.p_K[threshold]

"""
P_K_needed[i][j]: only relevant values of p_K_empirical
the empirically measured percentage of 
    minimizers with length == MINIMIZER_LENGTHS[i] 
    which have a frequency < THRESHOLDS[j]
""" 
P_K_needed = [[P_K_empirical[m][t] for t in THRESHOLDS] for m in range(len(MINIMIZER_LENGTHS))]
#print(f"P_K_needed:\n {P_K_needed}")

"""
f_avg[i][j]: 
the empirically measured average frequency of 
    minimizers with length == MINIMIZER_LENGTHS[i] 
    which have a frequency < THRESHOLDS[j]
""" 
f_avg_empirical = [[] for m in MINIMIZER_LENGTHS] 
for m in range(len(MINIMIZER_LENGTHS)):
    current = 0
    for t in range(1, THRESHOLDS[-1]+1):
        current = current + (P_K_empirical[m][t]-P_K_empirical[m][t-1])*t
        if t in THRESHOLDS:
            f_avg_empirical[m].append(current/P_K_empirical[m][t])
print(f"f_avg_empirical:\n {f_avg_empirical}")

#============================#
#Set up SCIP Model
#============================#
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

p_K_vars = [scip.addVar(f"p_K_{l}") for l in range(LEVELS)]     #percentage of kmers that have minimizers of length m with a frequency <= t for each level
f_avg_vars = [scip.addVar(f"f_avg_{l}") for l in range(LEVELS)] #average frequency of minimizers of length m with a frequency <= t for each level

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
    p_K_sum_terms = []
    f_avg_sum_terms = []
    for minimizer_lengths in product(MINIMIZER_LENGTHS, repeat=l):
        m_indicator = 1
        for l, m in enumerate(minimizer_lengths):
            m_indicator *= im_vars[(m, l)]

        for thresholds in product(THRESHOLDS, repeat=l):
            t_indicator = 1
            for i, t in enumerate(thresholds):
                t_indicator *= it_vars[(t, l)]

            p_K_empiric = getCascadingStats(
                l, minimizer_lengths, thresholds
            )
            f_avg_empiric = #TODO!
            p_K_sum_terms.append(m_indicator * t_indicator * p_K_empiric)
            f_avg_sum_terms.append(m_indicator * t_indicator * f_avg_empiric)
    
    scip.addCons(p_K_vars[l] - quicksum(p_K_sum_terms) == 0)
    scip.addCons(f_avg_vars[l] - quicksum(f_avg_sum_terms) == 0)
    
    scip.addCons(t_level_vars[l] - (t_rank_R + 2*t_select_S + (f_avg_vars[l]+1)/2*(t_access_O+t_compare_T)) == 0)

#TODO: From here
    #bits per k-mer in O
    scip.addCons(bpk_O_vars[l] - (quickprod((1-p_K_vars[i]) for i in range(l)) * p_K_vars[l] * 2/(K-m_vars[l]+2) * log_sequence_length) == 0)


    
    
scip.addCons(bpk_O - quicksum(bpk_O_vars[l] for l in range(LEVELS)) == 0)

#recursivly build the t_lookup constraint level by level
remaining_time = t_access_H
for l in reversed(range(LEVELS)):
    remaining_time = t_access_R + p_K_vars[l] * t_level_vars[l] + (1 - p_K_vars[l]) * remaining_time

scip.addCons(t_lookup - remaining_time == 0, name="t_lookup")

#bits per k-mer in H
scip.addCons(bpk_H - (quickprod(((1-p_K_vars[l])) for l in range(LEVELS)) * (64+log_sequence_length)/LOAD_FACTOR) == 0)
scip.addCons(bpk_total - (bpk_O + bpk_H) == 0)
if(OPTIMIZATION_GOAL == "t_lookup"):
    scip.addCons(bpk_total <= MAX_BPK_TOTAL)

#============================#
#SOLVE
#============================#
scip.optimize()
sol = scip.getBestSol()

#============================#
#OUTPUT SOLUTION
#============================#
'''
for l in range(LEVELS):
    print(f"level {l}:")
    for i, m in enumerate(MINIMIZER_LENGTHS):
        print(f"m_{m}_{l}: {format(sol[im_vars[l][i]])}")
    for i, t in enumerate(THRESHOLDS):
        print(f"t_{t}_{l}: {format(sol[it_vars[l][i]])}")
    print("")
'''
print(f"\nsolved MINLP for sequence: {SEQUENCE_NAME}")
if(OPTIMIZATION_GOAL == "t_lookup"):
    print(f"with optimization goal: {OPTIMIZATION_GOAL} and MAX_BPK_TOTAL = {MAX_BPK_TOTAL}")
else:
    print(f"with optimization goal: {OPTIMIZATION_GOAL}")
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
    print(f"bpk_O_{l}: {format(sol[bpk_O_vars[l]])}")
print(f"bpk_H: {format(sol[bpk_H])}")
print(f"t_lookup: {format(sol[t_lookup])}")
print(f"bpk_total: {format(sol[bpk_total])}")

#============================#
#WRITE SOLUTION TO LOG
#============================#
with open(logfile_path, "a") as logfile:
    date_time = datetime.now().strftime("%m/%d/%Y %H:%M:%S")
    with open(logfile_path, "a") as logfile:
        logfile.write(f"{date_time}----------------------------------------\nsolved MINLP for sequence: {SEQUENCE_NAME}\n")
    if(OPTIMIZATION_GOAL == "t_lookup"):
        logfile.write(f"with optimization goal: {OPTIMIZATION_GOAL} and MAX_BPK_TOTAL = {MAX_BPK_TOTAL}\n")
    else:
        logfile.write(f"with optimization goal: {OPTIMIZATION_GOAL}\n")
    logfile.write(f"for {LEVELS} levels and k = {K}\nconsidered minimizer lengths: {MINIMIZER_LENGTHS}\nconsidered thresholds: {THRESHOLDS}\n")

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
        logfile.write(f"bpk_O_{l}: {format(sol[bpk_O_vars[l]])}\n")
    logfile.write(f"bpk_H: {format(sol[bpk_H])}\n")
    logfile.write(f"t_lookup: {format(sol[t_lookup])}\n")
    logfile.write(f"bpk_total: {format(sol[bpk_total])}\n\n\n")
    