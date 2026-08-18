from pyscipopt import Model, quicksum, quickprod
import math
import csv
import sys

genome = "homo_sapiens"
K = 31  #kmer length
T = 255 #max threshold
LOAD_FACTOR = 0.75  #of the hash table

#TODO: model or measure!
t_access_R_empirical = 30
t_access_O_empirical = 30
t_access_H_empirical = 20
t_rank_R_empirical = 35
t_select_S_empirical = 35
t_compare_T_empirical = 30


TEXT_LENGTH = 681137419 #TODO: derive from file
log_text_length = math.log(TEXT_LENGTH, 2)


levels = 3
minimizer_lengths = range(5,20)
thresholds = [2**i for i in range(1,7)]

if thresholds[-1] >= T:
    sys.exit(f"maximum threshold t: {thresholds[-1]} has to be smaller than T: {T}")

"""
P_K_empirical[i][j]: 
the empirically measured percentage of 
    minimizers with length == minimizer_lengths[i] 
    which have a frequency < j
""" 
P_K_empirical = [[0] for m in minimizer_lengths] #percentage of appearing minimizers with f==0 is 0
for i, m in enumerate(minimizer_lengths):  
    reader = csv.reader(open(f'./minimizer_stats/{genome}/stats_m{m}_t{T}.csv'), delimiter=',')
    for t, row in enumerate(reader):
        if(t>0 and t<=T):
            P_K_empirical[i].append(float(row[2])/100)   #convert to percentage
"""
P_K_needed[i][j]: only relevant values of p_K_empirical
the empirically measured percentage of 
    minimizers with length == minimizer_lengths[i] 
    which have a frequency < thresholds[j]
""" 
P_K_needed = [[P_K_empirical[m][t] for t in thresholds] for m in range(len(minimizer_lengths))]
print(P_K_needed)

"""
f_avg[i][j]: 
the empirically measured average frequency of 
    minimizers with length == minimizer_lengths[i] 
    which have a frequency < thresholds[j]
""" 
f_avg_empirical = [[] for m in minimizer_lengths] 
for m in range(len(minimizer_lengths)):
    current = 0
    for t in range(1, thresholds[-1]+1):
        current = current + (P_K_empirical[m][t]-P_K_empirical[m][t-1])*t
        if t in thresholds:
            f_avg_empirical[m].append(current)
print(f_avg_empirical)

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

p_K_vars = [scip.addVar(f"p_K_{l}") for l in range(levels)]     #percentage of kmers that have minimizers of length m with a frequency <= t for each level

f_avg_vars = [scip.addVar(f"f_avg_{l}") for l in range(levels)]      #average frequency of minimizers of length m with a frequency <= t for each level


sm_vars = [[scip.addVar(vtype='B', name=f"sm_{m}_{l}") for m in minimizer_lengths] for l in range(levels)] #binary selectors for m
st_vars = [[scip.addVar(vtype='B', name=f"st_{t}") for t in thresholds] for l in range(levels)]            #binary selectors for t

m_vars = [scip.addVar(f"m_{l}") for l in range(levels)]         #minimizer length on each level
t_vars = [scip.addVar(f"t_{l}") for l in range(levels)]         #threshold on each level

bpk_O_vars = [scip.addVar(f"bpk_O_{l}") for l in range(levels)] #bits per kmer in O on each level
bpk_O = scip.addVar(f"bpk_O")
bpk_H = scip.addVar(f"bpk_H")                                   #bits per kmer in H

#============================#
#OBJECTIVE
#============================#
scip.setObjective(runtime, sense="minimize")   #runtime objective

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
for l in reversed(range(levels)):
    remaining_time = t_access_R * p_K_vars[l] * (t_rank_R + 2*t_select_S + f_avg_vars[l]/2*(t_access_O+t_compare_T)) + (1 - p_K_vars[l]) * remaining_time

scip.addCons(runtime - remaining_time == 0)

for l in range(levels):
    #determine p_K_l
    scip.addCons(p_K_vars[l] - quicksum(sm * quicksum((st * P_K_needed[i][j]) for j, st in enumerate(st_vars[l])) for i, sm in enumerate(sm_vars[l])) == 0)
    #determine f_avg_l
    scip.addCons(f_avg_vars[l] - quicksum(sm * quicksum((st * f_avg_empirical[i][j]) for j, st in enumerate(st_vars[l])) for i, sm in enumerate(sm_vars[l])) == 0)
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
scip.addCons(bpk_total <= 8)

#============================#
#SOLVE
#============================#
scip.optimize()
sol = scip.getBestSol()

for l in range(levels):
    print(f"level {l}:")
    for i, m in enumerate(minimizer_lengths):
        print(f"m_{m}_{l}: {format(sol[sm_vars[l][i]])}")
    for i, t in enumerate(thresholds):
        print(f"t_{t}_{l}: {format(sol[st_vars[l][i]])}")
    print("")

#============================#
#OUTPUT SOLUTION
#============================#
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
    