from pyscipopt import Model, quicksum
import csv

T_SL = 1    #time for lookup in S
T_H = 2     #time for lookup in hashtable
LEVELS = 1

minimizer_lengths = range(15,25)
thresholds = [16, 32, 64]

"""
P_Kempirical[i][j]: 
the empirically measured percentage of 
    minimizers with length == minimizer_lengths[i] 
    which have a frequency < thresholds[j]
""" 
P_Kempirical = [[] for m in minimizer_lengths] 
for i, m in enumerate(minimizer_lengths):
    reader = csv.reader(open(f'stats_m{m}.csv'), delimiter=',')
    for t, row in enumerate(reader):
        if t not in thresholds:
            continue
        P_Kempirical[i].append(float(row[2]))

scip = Model("Level1")  # scip name is optional

#variables
p_K = scip.addVar("p_K")

sm_vars = [scip.addVar(vtype='B', name=f"sm_{m}") for m in minimizer_lengths]
st_vars = [scip.addVar(vtype='B', name=f"st_{t}") for t in thresholds]

"""
sm = [[None for i in range(len(minimizer_lengths))] for i in range(LEVELS)]
st = [[None for i in range(len(thresholds))] for i in range(LEVELS)]
for l in range(LEVELS):
    for i, m in enumerate(minimizer_lengths):
        sm[l][i] = scip.addVar(vtype='B', name=f"sm_L{l}_m{m}")
        
    for i, t in enumerate(thresholds):
        st[l][i] = scip.addVar(vtype='B', name=f"st_L{l}_t{t}")
"""

#objective
scip.setObjective(T_SL * p_K + T_H * (1 - p_K), sense="minimize")

#constraints
scip.addCons(p_K - quicksum(sm + quicksum((st * P_Kempirical[i][j]) for j, st in enumerate(st_vars)) for i, sm in enumerate(sm_vars)) == 1)
scip.addCons(quicksum(sm for sm in sm_vars) == 1, name="sm_sum_cons")
scip.addCons(quicksum(st for st in st_vars) == 1, name="st_sum_cons")

scip.optimize()
sol = scip.getBestSol()
for i, m in enumerate(minimizer_lengths):
    print(f"m_{m}: {format(sol[sm_vars[i]])}")
for i, t in enumerate(thresholds):
    print(f"t_{t}: {format(sol[st_vars[i]])}")