import numpy as np
import json

import sys
import os
import matplotlib.pyplot as plt
from itertools import product
from datetime import datetime

from pathlib import Path
from dataclasses import fields
from collections import defaultdict

from cascadingStats import *
from MINLP import MINLP

#temporary!! TODO: measure and pass as part of cascading stats
K = 31  #kmer length that is/was used to gain minimizer_stats
LOAD_FACTOR = 0.75  #of the hash table

GENOMES = ["Falco_tinnunculus.FalTin1.0.dna.toplevel",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Gadus_morhua.gadMor3.0.dna.toplevel",
            "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]


OPTIMIZATION_GOAL = "t_lookup"   #"t_lookup" or "bpk_total"
MAX_LOOKUP_TIME = np.inf
MAX_BPK_TOTAL = np.inf

MINIMIZER_LENGTHS_CONSIDERED = [m for m in range(15,26,2)]
THRESHOLDS_CONSIDERED = [2**i for i in range(5,8)]  

NO_LEVELS = 3
NON_DECREASING_PARAMETERS = True

OVERWRITE_PLOTS = False

#TODO: tmp, for testing
forcedMinimizerLength = []
forcedThresholds = []

evaluation_results = Path(f"evaluation_results")
evaluation_results.mkdir(exist_ok=True, parents=True)

for sequence_name in GENOMES:
    sequence_path = f'../../../rshash/DNA_datasets/{sequence_name}.fa'
    sequence_length = 0
    if os.path.isfile(sequence_path):
        sequence_length = os.path.getsize(sequence_path)
    else:
        sys.exit(f"{sequence_path} is not a file")
    log_sequence_length = np.log2(sequence_length)

    sequence_path = evaluation_results / f"{sequence_name}"
    sequence_path.mkdir(exist_ok=True, parents=True)
    
    print(f"evaluating model for sequence: {sequence_name}")
    data_path = sequence_path / "data"
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
        genome_stats_directory = Path(f"../allStats/{sequence_name}")
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
                forceParameters=[forcedMinimizerLengths,forcedThresholds])
            
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
            stats_class = type(modeled).__name__

            for stat_field in fields(modeled):
                name = stat_field.name
                modeled_value = getattr(modeled, name)
                measured_value = getattr(measured, name)

                if isinstance(modeled_value, list):
                    for level, (m, obs) in enumerate(zip(modeled_value, measured_value)):
                        key = (stats_class, name, level)
                        errors.setdefault(key, []).append(m - obs)
                else:
                    key = (stats_class, name, None)
                    errors.setdefault(key, []).append(modeled_value - measured_value)

        # write errors to json file
        with open(str( data_path / 'errors.json'), 'w') as f:
            json.dump([[list(k), v] for k, v in errors.items()], f)

    ##############
    # code for plotting data starts here
    ##############
    genome_plots = sequence_path / "plots"
    genome_plots.mkdir(exist_ok=True, parents=True)
    ##############
    # Compare modeled vs measured lookup time
    # Generative AI use Disclaimer:
    #   the following code for plotting the measured vs the the modeled VALUES was generated using GPT-6 Luna
    ##############
    # extract stats pairs from modeled_AND_measured_stats
    stats_classes = (spaceConsumption_Stats, operationTimes_Stats)
    # Prepare storage for each class and field.
    plot_data = {
        stats_class: {
            stat_field.name: {"modeled": [], "measured": []}
            for stat_field in fields(stats_class)
        }
        for stats_class in stats_classes
    }

    # Iterate over modeled_AND_measured_stats only once.
    for modeled, measured in modeled_AND_measured_stats:
        for stats_class in stats_classes:
            if isinstance(modeled, stats_class) and isinstance(measured, stats_class):
                for stat_field in fields(stats_class):
                    field_name = stat_field.name
                    modeled_value = getattr(modeled, field_name)
                    measured_value = getattr(measured, field_name)
                    data = plot_data[stats_class][field_name]
                    
                    #TODO: REMOVE, only temporary
                    #if(field_name == "Lookup_time" and (modeled_value < 300 or modeled_value > 500)):
                    #    continue

                    if isinstance(modeled_value, list) and isinstance(measured_value, list):
                        # Pair corresponding list elements.
                        for modeled_item, measured_item in zip(modeled_value, measured_value):
                            data["modeled"].append(modeled_item)
                            data["measured"].append(measured_item)
                    else:
                        data["modeled"].append(modeled_value)
                        data["measured"].append(measured_value)

    print("computing correlation coefficients - ATTENTION - plot_data NO LONGER DIFFERENTIATES BETWEEN LEVELS!!") 
    correlation_stats_path = sequence_path / "correlation_coefficients.csv"
    fieldnames = []
    correlation_coefficients = []
    for stats_class in stats_classes:
        for field_name, data in plot_data[stats_class].items():
            modeled_values = data["modeled"]
            measured_values = data["measured"]
            
            fieldnames.append(field_name)
            correlation_coefficients.append(np.corrcoef(modeled_values, measured_values)[0, 1])
            
    with open(correlation_stats_path, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(fieldnames)
        writer.writerow(correlation_coefficients)

    #create plots
    plots_path = genome_plots / "Modeled vs Measured"
    if(not OVERWRITE_PLOTS and plots_path.is_dir()):
        print(f"Modeled vs Measured plots already exists at {plots_path}")
    else:
        plots_path.mkdir(exist_ok=True, parents=True)
        print(f"plotting Modeled vs Measured values plots")
        
        # Create one plot for each field.
        for stats_class in stats_classes:
            for field_name, data in plot_data[stats_class].items():
                modeled_values = data["modeled"]
                measured_values = data["measured"]

                if not modeled_values:
                    continue

                title = f"{field_name} for {sequence_name.replace('.','_')}"

                #plot red diagonal where x=y
                plt.axline((0,0),slope=1,linestyle='-',color='red',label="modeled=measured")
                plt.legend()
                
                plt.scatter(measured_values, modeled_values)
                plt.xlabel(f"Measured {field_name}")
                plt.ylabel(f"Modeled {field_name}")
                plt.title(title)
                plt.grid(True)
                plt.tight_layout()
                plt.savefig(plots_path / f"{title}.png", dpi=300)
                plt.close()
    
    ##############
    # Make one box plot per field/level pair. 
    # Generative AI use Disclaimer:
    #   the following code for creating box plots of the squarred_errors was generated using GPT-6 Luna
    ##############
    for error_type in ("absolute_error", "squared_error"):
        plots_path = genome_plots / error_type
        if(not OVERWRITE_PLOTS and plots_path.is_dir()):
            print(f"{error_type} plots already exists at {plots_path}")
        else:
            plots_path.mkdir(exist_ok=True, parents=True)
            print(f"plotting {error_type} plots")
            
            grouped_errors = defaultdict(dict)
            for (stats_class, name, level), error in errors.items():
                if error_type == "absolute_error":
                    grouped_errors[(stats_class, name)][level] = [np.abs(e) for e in error]
                else: 
                    grouped_errors[(stats_class, name)][level] = [e**2 for e in error]
                    
            for (stats_class, name), errors_by_level in grouped_errors.items():
                # Sort levels consistently; None appears first if present.
                levels = sorted(errors_by_level, key=lambda level: (level is not None, level))
                labels = [str(level) if level is not None else "None" for level in levels]
                data = [errors_by_level[level] for level in levels]

                title = f"{name} for {sequence_name.replace('.','_')}"

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
                
    print("\n")


    #============================#
    #OUTPUT INFO
    #============================#
    #info.insert(0, f"{datetime.now()}---------------------------------------------\n")
    #info.insert(1, f"solved MINLP for sequence: {sequence_name}\n")
    #for line in info:
    #    print(line, end='')
    #============================#
    #WRITE SOLUTION TO LOG
    #============================#
    #logfile_path = f'MINLP_logs/{sequence_name}.txt'
    #with open(logfile_path, "w") as logfile:    #TODO: change to mode a - append, once consistent
    #    for line in info:
    #        logfile.write(line)