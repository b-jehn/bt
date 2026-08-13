from pyscipopt import Model, quicksum, quickprod
import math
import csv
import sys

genome = "homo_sapiens"
K = 31  #kmer length
T = 255 #max threshold
TEXT_LENGTH = 681137419
log_text_length = math.log(TEXT_LENGTH, 2)
LOAD_FACTOR = 0.75

T_SL = 2    #time for lookup in S
T_H = 1     #time for lookup in hashtable

levels = 3
minimizer_lengths = range(5,20)
thresholds = [2**i for i in range(1,7)]

if thresholds[-1] >= T:
    sys.exit(f"maximum threshold t: {thresholds[-1]} has to be smaller than T: {T}")

"""
P_Kempirical[i][j]: 
the empirically measured percentage of 
    minimizers with length == minimizer_lengths[i] 
    which have a frequency < thresholds[j]
""" 
P_Kempirical = [[] for m in minimizer_lengths] 
for i, m in enumerate(minimizer_lengths):
    reader = csv.reader(open(f'./minimizer_stats/{genome}/stats_m{m}_t{T}.csv'), delimiter=',')
    for t, row in enumerate(reader):
        if t not in thresholds:
            continue
        P_Kempirical[i].append(float(row[2])/100)   #convert to percentage
        
print(P_Kempirical)

scip = Model("Level")  # scip name is optional

#variables
runtime = scip.addVar("runtime")
bpk_total = scip.addVar("bpk_total")
p_K_vars = [scip.addVar(f"p_K_{l}") for l in range(levels)]     #percentage of kmers that have minimizers of length m with a frequency < t for each level
m_vars = [scip.addVar(f"m_{l}") for l in range(levels)]         #minimizer length on each level
t_vars = [scip.addVar(f"t_{l}") for l in range(levels)]         #threshold on each level
bpk_O_vars = [scip.addVar(f"bpk_O_{l}") for l in range(levels)] #bits per kmer in O on each level
bpk_O = scip.addVar(f"bpk_O")
bpk_H = scip.addVar(f"bpk_H")                                   #bits per kmer in H


sm_vars = [[scip.addVar(vtype='B', name=f"sm_{m}_{l}") for m in minimizer_lengths] for l in range(levels)] #binary selectors for m
st_vars = [[scip.addVar(vtype='B', name=f"st_{t}") for t in thresholds] for l in range(levels)]            #binary selectors for t

#OBJECTIVE
scip.setObjective(bpk_total, sense="minimize")   #runtime objective
#CONSTRAINTS
#runtime constraint
#recursivly build the constraint which is layerd level by level
remaining_time = T_H
for p_K in reversed(p_K_vars[:levels]):
    remaining_time = p_K * T_SL + (1 - p_K) * remaining_time

scip.addCons(runtime - remaining_time == 0)

for l in range(levels):
    #determine p_K_l
    scip.addCons(p_K_vars[l] - quicksum(sm * quicksum((st * P_Kempirical[i][j]) for j, st in enumerate(st_vars[l])) for i, sm in enumerate(sm_vars[l])) == 0)
    #find m and t
    scip.addCons(quicksum(sm for sm in sm_vars[l]) == 1, name=f"sm_mutually_exclusive_{l}")     #mutually exclusive selectors for m
    scip.addCons(quicksum(st for st in st_vars[l]) == 1, name=f"st__mutually_exclusive_{l}")    #mutually exclusive selectors for t
    #get m and t
    scip.addCons(m_vars[l] - quicksum(sm * m for sm, m in zip(sm_vars[l], minimizer_lengths)) == 0, name=f"get_m_{l}")  #determin m
    scip.addCons(t_vars[l] - quicksum(st * t for st, t in zip(st_vars[l], thresholds)) == 0, name=f"get_t_{l}")         #determin t

    #bits per k-mer in O
    scip.addCons(bpk_O_vars[l] - (quickprod((1-p_K_vars[i]) * 2/(K-m_vars[i]+2) for i in range(l)) * p_K_vars[l] * log_text_length) == 0)

scip.addCons(bpk_O - quicksum(bpk_O_vars[l] for l in range(levels)) == 0)

#bits per k-mer in H
scip.addCons(bpk_H - (quickprod(((1-p_K_vars[l])) for l in range(levels)) * (64+log_text_length)) == 0)
scip.addCons(bpk_total - (bpk_O + bpk_H) == 0)
#scip.addCons(bpk_total <= 40)

#FIND solution
scip.optimize()
sol = scip.getBestSol()

for l in range(levels):
    print(f"level {l}:")
    for i, m in enumerate(minimizer_lengths):
        print(f"m_{m}_{l}: {format(sol[sm_vars[l][i]])}")
    for i, t in enumerate(thresholds):
        print(f"t_{t}_{l}: {format(sol[st_vars[l][i]])}")
    print("")
    
print("found solution:")
for l in range(levels):
    print(f"m_{l}: {format(sol[m_vars[l]])}")
print()
for l in range(levels):
    print(f"t_{l}: {format(sol[t_vars[l]])}") 
print() 
for l in range(levels):  
    print(f"p_K_{l}: {format(sol[p_K_vars[l]])}")
print()
for l in range(levels):
    print(f"bpk_O_{l}: {format(sol[bpk_O_vars[l]])}")
print(f"bpk_H: {format(sol[bpk_H])}")
print(f"runtime: {format(sol[runtime])}")
print(f"bpk_total: {format(sol[bpk_total])}")
    