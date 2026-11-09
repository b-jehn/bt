import numpy as np
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

GENOMES = ["Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]

SEQUENCE_NAME = GENOMES[2]
OPTIMIZATION_GOAL = "t_lookup"   #"t_lookup" or "bpk_total"
MAX_LOOKUP_TIME = np.inf
MAX_BPK_TOTAL = np.inf

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
log_sequence_length = np.log2(sequence_length)


#TODO: model or measure!
t_access_R_empirical = 30
t_access_O_empirical = 10
t_access_H_empirical = 20
t_rank_R_empirical = 35
t_select_S_empirical = 35
t_compare_T_empirical = 10

data_path = Path(f"{SEQUENCE_NAME}/data")
data_path.mkdir(exist_ok=True, parents=True)


modeled_AND_measured_stats_jsonfile_path = data_path / 'modeled_AND_measured_stats.json'
modeled_AND_measured_stats = []
if(modeled_AND_measured_stats_jsonfile_path.exists()):
    print("loading modeled_AND_measured_stats from previous computations json file")
    # Load modeled_AND_measured_stats from json file
    with open(modeled_AND_measured_stats_jsonfile_path, "r", encoding="utf-8") as file:
        modeled_AND_measured_stats = [
            tuple(pair) for pair in modeled_AND_measured_stats_from_jsonable(json.load(file))
    ]
else:
    print("computing modeled_AND_measured_stats")
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

    setOfMinimizerLengths = []
    for m1 in range(15,26,2):
        for m2 in range(m1,26,2):
            for m3 in range(m2,26,2):
                setOfMinimizerLengths.append([m1,m2,m3])
                
    setOfThresholds = []
    for t1 in range(5,8):
        for t2 in range(t1,8):
            for t3 in range(t2,8):
                setOfThresholds.append([2**t1,2**t2,2**t3])
                
    #setOfMinimizerLengths = [[17,17,17]]
    #setOfThresholds = [[128,128,128]]

    for forcedMinimizerLengths, forcedThresholds in product(setOfMinimizerLengths, setOfThresholds):
        measuredSpace, measuredTimes = getMeasuredStats(genome_stats_directory, 
                                                        forcedMinimizerLengths,
                                                        forcedThresholds)
        
        print(f"running MINLP with forced minimizer lengths: {forcedMinimizerLengths} and thresholds: {forcedThresholds}")
        modeledSpace, modeledTimes, info, solution, elapsed_seconds = MINLP(
            optimization_goal="t_lookup",
            max_lookup_time=MAX_LOOKUP_TIME,
            max_bpk_total=MAX_BPK_TOTAL,
            log_sequence_length=log_sequence_length,
            no_levels=NO_LEVELS, 
            minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED, 
            thresholdsConsidered=THRESHOLDS_CONSIDERED, 
            distribution_stats=genome_stats,
            load_factor=LOAD_FACTOR,
            nonDecreasingParameters=NON_DECREASING_PARAMETERS,
            forceParameters=[forcedMinimizerLengths,forcedThresholds],
            operationTimes=measuredTimes)
        
        modeled_AND_measured_stats.append((modeledSpace, measuredSpace))
        modeled_AND_measured_stats.append((modeledTimes, measuredTimes))
        
        
    # write modeled_AND_measured_stats to json file
    with open(str(modeled_AND_measured_stats_jsonfile_path), "w", encoding="utf-8") as file:
        json.dump(modeled_AND_measured_stats_to_jsonable(modeled_AND_measured_stats), file, indent=2)


errors_jsonfile_path = data_path / "errors.json"
errors = {}
if(errors_jsonfile_path.exists()):
    # Load errors from json file
    print("loading model errors from previous computations json file")
    with errors_jsonfile_path.open() as json_file:
        pairs = json.load(json_file)
        errors = {tuple(k): v for k, v in pairs}
else:
    ##############
    # squared error calculation and plotting for every stats item and level
    ##############
    ##############
    # Collect errors by stats class, field, and level.
    # Generative AI use Disclaimer:
    #   the following code for computing the squarred error for every pai of modeled and measured stats
    #   and every value modeled was generated using GPT-6 Luna
    ##############
    print("computing model errors")
    for modeled, measured in modeled_AND_measured_stats:
        stats_name = type(modeled).__name__

        for stat_field in fields(modeled):
            name = stat_field.name
            modeled_value = getattr(modeled, name)
            measured_value = getattr(measured, name)

            if isinstance(modeled_value, list):
                for level, (m, obs) in enumerate(zip(modeled_value, measured_value)):
                    key = (stats_name, name, level)
                    errors.setdefault(key, []).append(m - obs)
            else:
                key = (stats_name, name, None)
                errors.setdefault(key, []).append(modeled_value - measured_value)

    # write errors to json file
    with open(str( data_path / 'errors.json'), 'w') as f:
        json.dump([[list(k), v] for k, v in errors.items()], f)

genome_plots = Path(f"{SEQUENCE_NAME}/plots")
genome_plots.mkdir(exist_ok=True, parents=True)

##############
# Compare modeled vs measured lookup time
# Generative AI use Disclaimer:
#   the following code for plotting the measured lookup times vs the the modeled ones was generated using GPT-6 Luna
##############
# extract (modeledTimes.Lookup_time, measuredTimes.Lookup_time) pairs from modeled_AND_measured_stats
lookup_time_pairs = [
    (modeled.Lookup_time, measured.Lookup_time)
    for modeled, measured in modeled_AND_measured_stats
    if isinstance(modeled, operationTimes_Stats)
    and isinstance(measured, operationTimes_Stats)
]

title = "Modeled vs Measured Lookup_time"
# x-axis: measured, y-axis: modeled
measured_values = [measured for modeled, measured in lookup_time_pairs]
modeled_values = [modeled for modeled, measured in lookup_time_pairs]

plt.scatter(measured_values, modeled_values)
plt.xlabel("Measured Lookup_time")
plt.ylabel("Modeled Lookup_time")
plt.title(title)
plt.grid(True)
plt.tight_layout()
plt.savefig(str(genome_plots / title), dpi=300)
plt.close()

##############
# Make one box plot per field/level pair. 
# Generative AI use Disclaimer:
#   the following code for creating box plots of the squarred_errors was generated using GPT-6 Luna
##############
from collections import defaultdict

for error_type in ("absolute_error", "squared_error"):
    print(f"plotting {error_type}s")
    plots_path = genome_plots / error_type
    plots_path.mkdir(exist_ok=True, parents=True)
    
    grouped_errors = defaultdict(dict)
    for (stats_name, name, level), error in errors.items():
        if error_type == "absolute_error":
            grouped_errors[(stats_name, name)][level] = [np.abs(e) for e in error]
        else: 
            grouped_errors[(stats_name, name)][level] = [e**2 for e in error]
            
    for (stats_name, name), errors_by_level in grouped_errors.items():
        # Sort levels consistently; None appears first if present.
        levels = sorted(errors_by_level, key=lambda level: (level is not None, level))
        labels = [str(level) if level is not None else "None" for level in levels]
        data = [errors_by_level[level] for level in levels]

        title = f"{stats_name}: {name}"

        fig, ax = plt.subplots()
        ax.boxplot(data, tick_labels=labels, showmeans=True, meanline=True)
        #ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.set_xlabel("Level")
        ax.set_ylabel(error_type)
        ax.set_title(title)

        # Label each box's median.
        for x, values in enumerate(data, start=1):
            median = np.median(values)
            ax.annotate(
                str(round(median, 3)),
                xy=(x + 0.2, median),
                va="center",
                ha="left",
            )

        fig.tight_layout()
        plt.savefig(str(plots_path / title), dpi=300)
        plt.close()


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