from pyscipopt import Model, quicksum
import math
import csv
import sys

K = 31  #kmer length
T = 255 #max threshold
TEXT_LENGTH = 681137419
log_text_length = math.log(TEXT_LENGTH, 2)
LOAD_FACTOR = 0.75

T_SL = 1    #time for lookup in S
T_H = 2     #time for lookup in hashtable
LEVELS = 1

minimizer_lengths = range(5,25)
thresholds = [2**i for i in range(1,8)]

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
    reader = csv.reader(open(f'./minimizer_stats/stats_m{m}_t{T}.csv'), delimiter=',')
    for t, row in enumerate(reader):
        if t not in thresholds:
            continue
        P_Kempirical[i].append(float(row[2])/100)   #convert to percentage
        
print(P_Kempirical)

scip = Model("Level1")  # scip name is optional

#variables
p_K = scip.addVar("p_K")    #percentage of kmers that have minimizers of length m with a frequency < t
m = scip.addVar("m")        #minimizer length
t = scip.addVar("t")        #threshold
bpk_O = scip.addVar("bpk_O")  #bits per kmer in O
bpk_H = scip.addVar("bpk_H")  #bits per kmer in O

sm_vars = [scip.addVar(vtype='B', name=f"sm_{m}") for m in minimizer_lengths]   #binary selectors for m
st_vars = [scip.addVar(vtype='B', name=f"st_{t}") for t in thresholds]          #binary selectors for t

"""
sm = [[None for i in range(len(minimizer_lengths))] for i in range(LEVELS)]
st = [[None for i in range(len(thresholds))] for i in range(LEVELS)]
for l in range(LEVELS):
    for i, m in enumerate(minimizer_lengths):
        sm[l][i] = scip.addVar(vtype='B', name=f"sm_L{l}_m{m}")
        
    for i, t in enumerate(thresholds):
        st[l][i] = scip.addVar(vtype='B', name=f"st_L{l}_t{t}")
"""

#OBJECTIVE
scip.setObjective(T_SL * p_K + T_H * (1 - p_K), sense="minimize")   #runtime objective

#CONSTRAINTS
scip.addCons(p_K - quicksum(sm * quicksum((st * P_Kempirical[i][j]) for j, st in enumerate(st_vars)) for i, sm in enumerate(sm_vars)) == 0)
#find m and t
scip.addCons(quicksum(sm for sm in sm_vars) == 1, name="sm_mutually_exclusive")     #mutually exclusive selectors for m
scip.addCons(quicksum(st for st in st_vars) == 1, name="st__mutually_exclusive")    #mutually exclusive selectors for t
#get m and t
scip.addCons(m - quicksum(sm * m for sm, m in zip(sm_vars, minimizer_lengths)) == 0, name="get_m")  #determin m
scip.addCons(t - quicksum(st * t for st, t in zip(st_vars, thresholds)) == 0, name="get_t")         #determin t

#bits per k-mer in O
scip.addCons(bpk_O - (2/(K-m+2) * p_K * log_text_length) == 0)
scip.addCons(bpk_O <= 10)

#bits per k-mer in H
scip.addCons(bpk_H - ((LOAD_FACTOR + (1-LOAD_FACTOR)/LOAD_FACTOR) * 2/(K-m+2) * (1-p_K) * log_text_length) == 0)
scip.addCons(bpk_H <= 10)

#FIND solution
scip.optimize()
sol = scip.getBestSol()

print(f"m: {format(sol[m])}")
print(f"t: {format(sol[t])}\n")

print(f"p_K: {format(sol[p_K])}")
print(f"bpk_O: {format(sol[bpk_O])}")
print(f"bpk_H: {format(sol[bpk_H])}\n")

for i, m in enumerate(minimizer_lengths):
    print(f"m_{m}: {format(sol[sm_vars[i]])}")
for i, t in enumerate(thresholds):
    print(f"t_{t}: {format(sol[st_vars[i]])}")
    