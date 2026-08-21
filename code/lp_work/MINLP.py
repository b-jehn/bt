from pyscipopt import Model, quicksum, quickprod
import math
import csv
import sys
import os

K = 31  #kmer length that is/was used to gain minimizer_stats
T = 255 #max threshold - the .csv files in ./minimizer_stats/{SEQUENCE_NAME}/ have to contain values for all t <= T
LOAD_FACTOR = 0.75  #of the hash table

SEQUENCE_NAME = "cod"
OPTIMIZATION_GOAL = "bpk_total"   #"runtime" or "bpk_total"
MAX_BPK_TOTAL = 5               #needed, when OPTIMIZATION_GOAL == "runtime"
LEVELS = 3
MINIMIZER_LENGTHS = [m for m in range(5,20)]
THRESHOLDS = [2**i for i in range(0,8)]

if THRESHOLDS[-1] >= T:
    sys.exit(f"maximum threshold t: {THRESHOLDS[-1]} has to be smaller than T: {T}")

sequence_path = f'../rshash/DNA_datasets/{SEQUENCE_NAME}.fa'
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
P_K_empirical = [[0] for m in MINIMIZER_LENGTHS] #percentage of appearing minimizers with f==0 is 0
for i, m in enumerate(MINIMIZER_LENGTHS):  
    reader = csv.reader(open(f'./minimizer_stats/{SEQUENCE_NAME}/stats_m{m}_t{T}.csv'), delimiter=',')
    for t, row in enumerate(reader):
        if(t>0 and t<=T):
            P_K_empirical[i].append(float(row[2])/100)   #convert to percentage
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

scip = Model("Level")  # scip name is optional

#============================#
#VARIABLES
#============================#
#optimization goals
runtime = scip.addVar("runtime")
bpk_total = scip.addVar("bpk_total")
#times for vector/hash table/text operations
t_access_R = scip.addVar("t_access_R")
t_access_O = scip.addVar("t_access_O")
t_access_H = scip.addVar("t_access_H")
t_rank_R = scip.addVar("t_rank_R")
t_select_S = scip.addVar("t_select_S")
t_compare_T = scip.addVar("t_compare_T")

p_K_vars = [scip.addVar(f"p_K_{l}") for l in range(LEVELS)]     #percentage of kmers that have minimizers of length m with a frequency <= t for each level

f_avg_vars = [scip.addVar(f"f_avg_{l}") for l in range(LEVELS)]      #average frequency of minimizers of length m with a frequency <= t for each level


sm_vars = [[scip.addVar(vtype='B', name=f"sm_{m}_{l}") for m in MINIMIZER_LENGTHS] for l in range(LEVELS)] #binary selectors for m on each level
st_vars = [[scip.addVar(vtype='B', name=f"st_{t}_{l}") for t in THRESHOLDS] for l in range(LEVELS)]        #binary selectors for t on each level

m_vars = [scip.addVar(f"m_{l}") for l in range(LEVELS)]         #minimizer length on each level
t_vars = [scip.addVar(f"t_{l}") for l in range(LEVELS)]         #threshold on each level

bpk_O_vars = [scip.addVar(f"bpk_O_{l}") for l in range(LEVELS)] #bits per kmer in O on each level
bpk_O = scip.addVar("bpk_O")
bpk_H = scip.addVar("bpk_H")                                   #bits per kmer in H
bpk_total = scip.addVar("bpk_total")
#============================#
#OBJECTIVE
#============================#
if(OPTIMIZATION_GOAL == "runtime"):
    scip.setObjective(runtime, sense="minimize")   #runtime objective
elif(OPTIMIZATION_GOAL == "bpk_total"):
    scip.setObjective(bpk_total, sense="minimize")   #runtime objective
else:
    sys.exit(f"OPTIMIZATION_GOAL should be 'runtime' or 'bpk_total' but is: {OPTIMIZATION_GOAL}")
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

#recursivly build the runtime constraint level by level
remaining_time = t_access_H
for l in reversed(range(LEVELS)):
    remaining_time = t_access_R + p_K_vars[l] * (t_rank_R + 2*t_select_S + f_avg_vars[l]/2*(t_access_O+t_compare_T)) + (1 - p_K_vars[l]) * remaining_time

scip.addCons(runtime - remaining_time == 0, name="runtime")

for l in range(LEVELS):
    #determine p_K_l
    scip.addCons(p_K_vars[l] - quicksum(sm * quicksum((st * P_K_needed[i][j]) for j, st in enumerate(st_vars[l])) for i, sm in enumerate(sm_vars[l])) == 0)
    #determine f_avg_l
    scip.addCons(f_avg_vars[l] - quicksum(sm * quicksum((st * f_avg_empirical[i][j]) for j, st in enumerate(st_vars[l])) for i, sm in enumerate(sm_vars[l])) == 0)
    #scip.addCons(f_avg_vars[l] - 2 == 0)
    #find m and t
    scip.addCons(quicksum(sm for sm in sm_vars[l]) == 1, name=f"sm_mutually_exclusive_{l}")     #mutually exclusive selectors for m
    scip.addCons(quicksum(st for st in st_vars[l]) == 1, name=f"st__mutually_exclusive_{l}")    #mutually exclusive selectors for t
    #get m and t
    scip.addCons(m_vars[l] - quicksum(sm * m for sm, m in zip(sm_vars[l], MINIMIZER_LENGTHS)) == 0, name=f"get_m_{l}")  #determin m
    scip.addCons(t_vars[l] - quicksum(st * t for st, t in zip(st_vars[l], THRESHOLDS)) == 0, name=f"get_t_{l}")         #determin t

    #bits per k-mer in O
    scip.addCons(bpk_O_vars[l] - (quickprod((1-p_K_vars[i]) for i in range(l)) * p_K_vars[l] * 2/(K-m_vars[l]+2) * log_sequence_length) == 0)

scip.addCons(bpk_O - quicksum(bpk_O_vars[l] for l in range(LEVELS)) == 0)

#bits per k-mer in H
scip.addCons(bpk_H - (quickprod(((1-p_K_vars[l])) for l in range(LEVELS)) * (64+log_sequence_length)) == 0)
scip.addCons(bpk_total - (bpk_O + bpk_H) == 0)
if(OPTIMIZATION_GOAL == "runtime"):
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
        print(f"m_{m}_{l}: {format(sol[sm_vars[l][i]])}")
    for i, t in enumerate(THRESHOLDS):
        print(f"t_{t}_{l}: {format(sol[st_vars[l][i]])}")
    print("")
'''
print(f"\nsolved MINLP for sequence: {SEQUENCE_NAME}")
if(OPTIMIZATION_GOAL == "runtime"):
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
print(f"runtime: {format(sol[runtime])}")
print(f"bpk_total: {format(sol[bpk_total])}")
    