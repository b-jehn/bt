import sys
import numpy as np
from pyscipopt import Model, quicksum, quickprod, log
from time import perf_counter
from itertools import product

from cascadingStats import *

def MINLP(
    optimization_goal: str,
    max_lookup_time: float,
    max_bpk_total: float,
    log_sequence_length: float,
    no_levels: int, 
    minimizerLengthsConsidered: list[int], 
    thresholdsConsidered: list[int], 
    distribution_stats: Distribution_Stats,
    load_factor=0.75,
    nonDecreasingParameters=False,
    forceParameters=[],
    hideOutput=True,
    operationTimes=operationTimes_Stats(),
    inverseTimes=False):
    #============================#
    #Set up SCIP Model
    #============================#
    scip = Model("Level")  # scip name is optional
    scip.hideOutput(hideOutput)

    #============================#
    #VARIABLES OR MEASUREMENTS??
    #============================#
    t_contains_R_vars = [scip.addVar(f"t_contains_R_{l}") for l in range(no_levels)]
    t_select_S_vars = [scip.addVar(f"t_select_S_{l}") for l in range(no_levels)]
    t_access_O_vars = [scip.addVar(f"t_access_O_{l}") for l in range(no_levels)]
    t_compare_T_vars = [scip.addVar(f"t_compare_T_{l}") for l in range(no_levels)]
    t_access_H = scip.addVar("t_access_H")

    bpk_T = 2.0   #bits per k-mer for storing the text T
    #============================#
    #VARIABLES
    #============================#
    im_vars = {                                             #binary indicators for m on each level l
        (m, l): scip.addVar(vtype="B", name=f"im_{m}_{l}")
        for m in minimizerLengthsConsidered
        for l in range(no_levels)
    }
    it_vars = {                                             #binary indicators for t on each level l
        (t, l): scip.addVar(vtype="B", name=f"it_{t}_{l}")
        for t in thresholdsConsidered
        for l in range(no_levels)
    }

    m_vars = [scip.addVar(vtype="I", name=f"m_{l}", lb=1, ub=minimizerLengthsConsidered[-1]) for l in range(no_levels)]         #minimizer length on each level
    t_vars = [scip.addVar(vtype="I", name=f"t_{l}", lb=1, ub=thresholdsConsidered[-1]) for l in range(no_levels)]         #threshold on each level
    
    #TODO: remove these two? if its possible to model without them
    # technically, these should be integers; however, the optimization runs much quickers with the variables being continues
    n_kmers_vars = [scip.addVar(vtype="C", name=f"n_kmers_{l}", lb=1) for l in range(no_levels)] #number of kmers on level l for each level l
    n_minimizers_vars = [scip.addVar(vtype="C", name=f"n_minimizers_{l}", lb=1) for l in range(no_levels)] #number of kmers on level l for each level l

    F_M_vars = [scip.addVar(vtype="C", name=f"F_M_{l}", lb=0.0, ub=1) for l in range(no_levels)]         #percentage of minimizers of length m with a frequency <= t for each level
    F_K_vars = [scip.addVar(vtype="C", name=f"F_K_{l}", lb=0.0, ub=1) for l in range(no_levels)]         #percentage of kmers that have minimizers of length m with a frequency <= t for each level
    f_avg_vars = [scip.addVar(vtype="C", name=f"f_avg_{l}", lb=1) for l in range(no_levels)]     #average frequency of minimizers of length m with a frequency <= t for each level

    t_level_vars = [scip.addVar(vtype="C", name=f"t_level_{l}", lb=0.0) for l in range(no_levels)] #time spent for finding k-mer with an infrequent minimizer on level l
    t_lookup = scip.addVar(vtype="C", name="t_lookup", lb=0.0)  #total lookup time

    n_DIM_vars = [scip.addVar(vtype="C", name=f"n_DIM_{l}", lb=1) for l in range(no_levels)] #number of distinct infrequent minimizers on each level

    bpk_R_vars = [scip.addVar(vtype="C", name=f"bpk_R_{l}", lb=0.0) for l in range(no_levels)] #bits per kmer in R on each level
    bpk_S_vars = [scip.addVar(vtype="C", name=f"bpk_S_{l}", lb=0.0) for l in range(no_levels)] #bits per kmer in S on each level
    bpk_O_vars = [scip.addVar(vtype="C", name=f"bpk_O_{l}", lb=0.0) for l in range(no_levels)] #bits per kmer in O on each level
    bpk_H = scip.addVar(vtype="C", name="bpk_H", lb=0.0)            #bits per kmer in H
    bpk_total = scip.addVar(vtype="C", name="bpk_total", lb=0.0)    #total bits per k-mer

    #============================#
    #OBJECTIVE
    #============================#
    if(optimization_goal == "t_lookup"):
        scip.setObjective(t_lookup, sense="minimize")   #t_lookup objective
    elif(optimization_goal == "bpk_total"):
        scip.setObjective(bpk_total, sense="minimize")   #space objective
    else:
        sys.exit(f"optimization_goal should be 't_lookup' or 'bpk_total' but is: {optimization_goal}")
    #============================#
    #CONSTRAINTS
    #============================#
    for l in range(no_levels):
        if(forceParameters):
            forcedMinimizerLengths, forcedThresholds = forceParameters
            scip.addCons(m_vars[l] == forcedMinimizerLengths[l])
            scip.addCons(t_vars[l] == forcedThresholds[l])
            
            n_kmers, n_minizers, F_M_empiric, F_K_empiric, f_avg_empiric = getCascadingDistribution_Stats(distribution_stats, forcedMinimizerLengths[:l+1], forcedThresholds[:l+1])
            scip.addCons(n_kmers_vars[l] == n_kmers)
            scip.addCons(n_minimizers_vars[l] == n_minizers)
            scip.addCons(F_M_vars[l] == F_M_empiric)
            scip.addCons(F_K_vars[l] == F_K_empiric)
            scip.addCons(f_avg_vars[l] == f_avg_empiric)
            
        else:
            if(nonDecreasingParameters and l!=0):
                scip.addCons(m_vars[l] >= m_vars[l-1])
                scip.addCons(t_vars[l] >= t_vars[l-1])

            scip.addCons(quicksum(im_vars[(m,l)] for m in minimizerLengthsConsidered) == 1, name=f"im_mutually_exclusive_{l}")   #mutually exclusive indicators for m
            scip.addCons(quicksum(it_vars[(t,l)] for t in thresholdsConsidered) == 1, name=f"it__mutually_exclusive_{l}")         #mutually exclusive indicators for t
            scip.addCons(m_vars[l] == quicksum(im_vars[(m,l)] * m for m in minimizerLengthsConsidered), name=f"get_m_{l}")   #determine m
            scip.addCons(t_vars[l] == quicksum(it_vars[(t,l)] * t for t in thresholdsConsidered), name=f"get_t_{l}")         #determine t
            
            #determine F_K_l and f_avg_l
            n_kmers_terms = []
            n_minimizers_terms = []
            F_M_sum_terms = []
            F_K_sum_terms = []
            f_avg_sum_terms = []
        
            for minimizer_lengths in product(minimizerLengthsConsidered, repeat=(l+1)):
                if(nonDecreasingParameters and not all(i <= j for i, j in zip(minimizer_lengths, minimizer_lengths[1:]))):
                    continue
                m_indicator = 1
                for m_l, m in enumerate(minimizer_lengths):
                    m_indicator *= im_vars[(m, m_l)]
                for thresholds in product(thresholdsConsidered, repeat=(l+1)):
                    if(nonDecreasingParameters and not all(i <= j for i, j in zip(thresholds, thresholds[1:]))):
                        continue
                    t_indicator = 1
                    for t_l, t in enumerate(thresholds):
                        t_indicator *= it_vars[(t, t_l)]
                    n_kmers, n_minizers, F_M_empiric, F_K_empiric, f_avg_empiric = getCascadingDistribution_Stats(distribution_stats, minimizer_lengths, thresholds)
                    n_kmers_terms.append(m_indicator * t_indicator * n_kmers)
                    n_minimizers_terms.append(m_indicator * t_indicator * n_minizers)
                    F_M_sum_terms.append(m_indicator * t_indicator * F_M_empiric)
                    F_K_sum_terms.append(m_indicator * t_indicator * F_K_empiric)
                    f_avg_sum_terms.append(m_indicator * t_indicator * f_avg_empiric)
        
            scip.addCons(n_kmers_vars[l] == quicksum(n_kmers_terms))
            scip.addCons(n_minimizers_vars[l] == quicksum(n_minimizers_terms))
            scip.addCons(F_M_vars[l] == quicksum(F_M_sum_terms))
            scip.addCons(F_K_vars[l] == quicksum(F_K_sum_terms))
            scip.addCons(f_avg_vars[l] == quicksum(f_avg_sum_terms))

        #scip.addCons(n_DIM_vars[l] == F_K_vars[l]*n_kmers_vars[l]/f_avg_vars[l] * 2/(K-m_vars[l]+2))
        scip.addCons(n_DIM_vars[l] == F_M_vars[l]*n_minimizers_vars[l]/f_avg_vars[l])

        #bits per k-mer in R,S and O
        scip.addCons(bpk_R_vars[l] == n_DIM_vars[l]/n_kmers_vars[0] * ((2*m_vars[l] - (log(n_DIM_vars[l])/log(2.0))) + 3))
        scip.addCons(bpk_S_vars[l] == F_M_vars[l]*n_minimizers_vars[l]/n_kmers_vars[0])
        scip.addCons(bpk_O_vars[l] == bpk_S_vars[l] * log_sequence_length)
        
        #outdated:
        #scip.addCons(bpk_S_vars[l] == 
        #            quickprod(1-F_K_vars[i] for i in range(l))
        #                * F_K_vars[l] * 2/(K-m_vars[l]+2)
        #            )
        
        #constraints regarding the structure operation times        
        if(inverseTimes):
            continue
        elif(len(operationTimes.Rcontains) == no_levels):
            scip.addCons(t_contains_R_vars[l] == operationTimes.Rcontains[l])
            scip.addCons(t_select_S_vars[l] == operationTimes.Sselect[l])
            scip.addCons(t_access_O_vars[l] == operationTimes.Oaccess[l])
            scip.addCons(t_compare_T_vars[l] == operationTimes.Tcompare[l])
            scip.addCons(t_access_H == operationTimes.Haccess)
        else: #fixed, really inaccurate values
            scip.addCons(t_contains_R_vars[l] == 230.0)
            scip.addCons(t_select_S_vars[l] == 22.0)
            scip.addCons(t_access_O_vars[l] == 6.0)
            scip.addCons(t_compare_T_vars[l] == 18.0)
            scip.addCons(t_access_H == 12.0)
    
        scip.addCons(t_level_vars[l] == 2*t_select_S_vars[l] + ((f_avg_vars[l]+1)/2)*(t_access_O_vars[l]+t_compare_T_vars[l]))

        
        #modeled values
        #scip.addCons(t_contains_R_vars[l] == 2*m_vars[l] - (log(n_DIM_vars[l])/log(2.0)))
        #scip.addCons(t_select_S_vars[l] == log(f_avg_vars[l])/log(2.0))
        

    remaining_time = t_access_H
    for l in reversed(range(no_levels)):
        remaining_time = t_contains_R_vars[l] + F_K_vars[l] * t_level_vars[l] + (1 - F_K_vars[l]) * remaining_time
        
    t_lookup_constraint = scip.addCons(t_lookup == remaining_time, name="t_lookup")
    scip.addCons(t_lookup <= max_lookup_time)
    if(inverseTimes):
        scip.addCons(t_lookup == operationTimes.Lookup_time)

    #bits per k-mer in H
    scip.addCons(bpk_H == quickprod(1-F_K_vars[l] for l in range(no_levels)) * (64 + 8)/load_factor)
    scip.addCons(bpk_total == bpk_T + quicksum(bpk_R_vars[l] + bpk_S_vars[l] + bpk_O_vars[l] for l in range(no_levels)) + bpk_H)
    scip.addCons(bpk_total <= max_bpk_total)

    #============================#
    #SOLVE
    #============================#
    if(not hideOutput): 
        scip.printCons(t_lookup_constraint)
        print("starting solving")
    start_time = perf_counter()
    scip.optimize()
    elapsed_seconds = perf_counter() - start_time
    
    status = scip.getStatus()
    if status != "optimal":
        print("SCIP status:", status)
        scip.writeProblem("model.cip")
        return None
    
    sol = scip.getBestSol()
    #============================#
    #CREATE INFO
    #============================#
    space = spaceConsumption_Stats()
    times = operationTimes_Stats()
    
    space.bpkT = bpk_T
    space.bpkH = sol[bpk_H]
    space.bpkTotal = sol[bpk_total]
    times.Haccess = sol[t_access_H]
    times.LookupTime = sol[t_lookup]
    for l in range(no_levels):
        space.bpkR.append(sol[bpk_R_vars[l]])
        space.bpkS.append(sol[bpk_S_vars[l]])
        space.bpkO.append(sol[bpk_O_vars[l]])
        
        times.Rcontains.append(sol[t_contains_R_vars[l]])
        times.Sselect.append(sol[t_select_S_vars[l]])
        times.Oaccess.append(sol[t_access_O_vars[l]])
        times.Tcompare.append(sol[t_compare_T_vars[l]])
        
    info = []
    info.append(f"with optimization goal: {optimization_goal},\nMAX_LOOKUP_TIME = {max_lookup_time} and\nMAX_BPK_TOTAL = {max_bpk_total}\n\n")

    info.append(f"for {no_levels} levels and k = {31}\n")
    info.append(f"considered minimizer lengths: {minimizerLengthsConsidered}\n")
    info.append(f"considered thresholds: {thresholdsConsidered}\n\n")

    info.append(f"total time to solve MINLP: {elapsed_seconds:.6f} seconds\n")
    info.append("found solution:\n")

    info.append("LEVEL PARAMETERS:\n")
    for l in range(no_levels):
        info.append(f"m_{l+1}: {format(sol[m_vars[l]])}\n")
    info.append("\n")
    for l in range(no_levels):
        info.append(f"t_{l+1}: {format(sol[t_vars[l]])}\n") 
    info.append("\n\n") 

    info.append("LEVEL STATS:\n")
    for l in range(no_levels):  
        info.append(f"n_kmers_{l+1}: {format(sol[n_kmers_vars[l]])}\n")
    info.append("\n") 
    for l in range(no_levels):  
        info.append(f"n_minimizers_{l+1}: {format(sol[n_minimizers_vars[l]])}\n")
    info.append("\n")
    for l in range(no_levels):  
        info.append(f"n_DIM_{l+1}: {format(sol[n_DIM_vars[l]])}\n")
    info.append("\n") 
    for l in range(no_levels):  
        info.append(f"F_M_{l+1}: {format(sol[F_M_vars[l]])}\n")
    info.append("\n") 
    for l in range(no_levels):  
        info.append(f"F_K_{l+1}: {format(sol[F_K_vars[l]])}\n")
    info.append("\n") 
    for l in range(no_levels):  
        info.append(f"f_avg_{l+1}: {format(sol[f_avg_vars[l]])}\n")
    info.append("\n\n")

    info.append("LOOKUP TIME:\n")
    info.append(f"t_lookup: {format(sol[t_lookup])} (ns for querying an appearing k-mer)\n\n")

    info.append("SPACE CONSUMPTION:\n")
    info.append(f"bpk_T: {bpk_T}\n")
    info.append("\n")
    for l in range(no_levels):
        info.append(f"bpk_R_{l+1}: {format(sol[bpk_R_vars[l]])}\n")
    info.append("\n")
    for l in range(no_levels):
        info.append(f"bpk_S_{l+1}: {format(sol[bpk_S_vars[l]])}\n")
    info.append("\n")
    for l in range(no_levels):
        info.append(f"bpk_O_{l+1}: {format(sol[bpk_O_vars[l]])}\n")
    info.append("\n")
    info.append(f"bpk_H: {format(sol[bpk_H])}\n")
    info.append("\n")
    info.append(f"=> bpk_total: {format(sol[bpk_total])}\n\n")
    
    return space, times, info, sol, elapsed_seconds