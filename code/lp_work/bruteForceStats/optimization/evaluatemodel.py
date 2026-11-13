import numpy as np
import json

import sys
import os
import matplotlib.pyplot as plt
from itertools import product
from datetime import datetime

from pathlib import Path
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

#GENOMES = [GENOMES[0]]


OPTIMIZATION_GOAL = "t_lookup"   #"t_lookup" or "bpk_total"
MAX_LOOKUP_TIME = np.inf
MAX_BPK_TOTAL = np.inf

MINIMIZER_LENGTHS_CONSIDERED = [m for m in range(15,26,2)]
THRESHOLDS_CONSIDERED = [2**i for i in range(5,8)]  

NO_LEVELS = 3
NON_DECREASING_PARAMETERS = True

OVERWRITE_PLOTS = False

setOfMinimizerLengths = []
for m1 in MINIMIZER_LENGTHS_CONSIDERED:
    for m2 in MINIMIZER_LENGTHS_CONSIDERED:
        if NON_DECREASING_PARAMETERS and m2 < m1:
            continue
        for m3 in MINIMIZER_LENGTHS_CONSIDERED:
            if NON_DECREASING_PARAMETERS and m3 < m2:
                continue
            setOfMinimizerLengths.append([m1,m2,m3])
            
setOfThresholds = []
for t1 in THRESHOLDS_CONSIDERED:
    for t2 in THRESHOLDS_CONSIDERED:
        if NON_DECREASING_PARAMETERS and t2 < t1:
            continue
        for t3 in THRESHOLDS_CONSIDERED:
            if NON_DECREASING_PARAMETERS and t3 < t2:
                continue
            setOfThresholds.append([t1,t2,t3])


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

    predicted_AND_observed_stats_jsonfile_path = data_path / 'predicted_AND_observed_stats.json'
    predicted_AND_observed_stats = []
    if(predicted_AND_observed_stats_jsonfile_path.exists()):
        print("loading predicted_AND_observed_stats from previous computations json file")
        # Load predicted_AND_observed_stats from json file
        with open(predicted_AND_observed_stats_jsonfile_path, "r", encoding="utf-8") as file:
            predicted_AND_observed_stats = [
                tuple(pair) for pair in predicted_AND_observed_stats_from_jsonable(json.load(file))
        ]
    else:
        print("computing predicted_AND_observed_stats")
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
                    
        #setOfMinimizerLengths = [[17,17,17]]
        #setOfThresholds = [[128,128,128]]

        for forcedMinimizerLengths, forcedThresholds in product(setOfMinimizerLengths, setOfThresholds):
            observedSpace, observedTimes = getObservedStats(genome_stats_directory, 
                                                            forcedMinimizerLengths,
                                                            forcedThresholds)
            
            print(f"running MINLP with forced minimizer lengths: {forcedMinimizerLengths} and thresholds: {forcedThresholds}")
            predictedSpace, predictedTimes, info, solution, elapsed_seconds = MINLP(
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
            
            predicted_AND_observed_stats.append((predictedSpace, observedSpace))
            predicted_AND_observed_stats.append((predictedTimes, observedTimes))
            
            
        # write predicted_AND_observed_stats to json file
        with open(str(predicted_AND_observed_stats_jsonfile_path), "w", encoding="utf-8") as file:
            json.dump(predicted_AND_observed_stats_to_jsonable(predicted_AND_observed_stats), file, indent=2)


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
        #   the following code for computing the squarred error for every pai of predicted and observed stats
        #   and every value predicted was generated using GPT-6 Luna
        ##############
        print("computing model errors")
        for predicted, observed in predicted_AND_observed_stats:
            stats_class = type(predicted).__name__

            for stat_field in fields(predicted):
                name = stat_field.name
                predicted_value = getattr(predicted, name)
                observed_value = getattr(observed, name)

                if isinstance(predicted_value, list):
                    for level, (pred, obs) in enumerate(zip(predicted_value, observed_value)):
                        key = (stats_class, name, level)
                        errors.setdefault(key, []).append(pred - obs)
                else:
                    key = (stats_class, name, None)
                    errors.setdefault(key, []).append(predicted_value - observed_value)

        # write errors to json file
        with open(str( data_path / 'errors.json'), 'w') as f:
            json.dump([[list(k), v] for k, v in errors.items()], f)
    
    ####################
    #sorted_lookup_times_path = sequence_path / "sorted_lookup_times.txt"
    #if(not sorted_lookup_times_path.is_file()):
    #    #sort parameter pairs by lookup time 
    #    print(f"writing sorted lookup times to: {sorted_lookup_times_path}")
    #    param_lookuptimes = {}
    #    for forcedMinimizerLengths, forcedThresholds in product(setOfMinimizerLengths, setOfThresholds):
    #        observedSpace, observedTimes = getObservedStats(Path(f"../allStats/{sequence_name}"), 
    #                                                        forcedMinimizerLengths,
    #                                                        forcedThresholds)
    #        params_string = str(getParameterDirectory(forcedMinimizerLengths, forcedThresholds))
    #        param_lookuptimes[params_string] = observedTimes.LookupTime
    #        
    #    with open(sorted_lookup_times_path, "w") as f:
    #        for params, lookuptime in sorted(param_lookuptimes.items(), key=lambda item: item[1]):
    #            f.write(f"{params}: {lookuptime}\n")
    #        f.close()
    #    
    
    #############
    # Stuff for Thesis gets stored here
    #############
    thesis_results_path = Path("../thesis_results")
    thesis_results_path.mkdir(parents=True, exist_ok=True)

    ##############
    # code for plotting data starts here
    ##############
    genome_plots = sequence_path / "plots"
    genome_plots.mkdir(exist_ok=True, parents=True)
    ##############
    # Compare predicted vs observed lookup time
    # Generative AI use Disclaimer:
    #   the following code for plotting the observed vs the the predicted VALUES was generated using GPT-6 Luna
    ##############
    # extract stats pairs from predicted_AND_observed_stats
    stats_classes = (spaceConsumption_Stats, operationTimes_Stats)
    
    # Prepare storage for each class and field.
    plot_data = {
        stats_class: {
            stat_field.name: [{"predicted": [], "observed": [], "diff": []}] * NO_LEVELS
            if get_origin(get_type_hints(stats_class)[stat_field.name]) is list
            else {"predicted": [], "observed": [], "diff": []}
            for stat_field in fields(stats_class)
        }
        for stats_class in stats_classes
    }

    # Iterate over predicted_AND_observed_stats only once.
    for predicted, observed in predicted_AND_observed_stats:
        for stats_class in stats_classes:
            if isinstance(predicted, stats_class) and isinstance(observed, stats_class):
                for stat_field in fields(stats_class):
                    field_name = stat_field.name
                    predicted_value = getattr(predicted, field_name)
                    observed_value = getattr(observed, field_name)
                    data = plot_data[stats_class][field_name]
                    
                    #TODO: REMOVE, only temporary
                    #if(field_name == "Lookup_time" and (predicted_value < 300 or predicted_value > 500)):
                    #    continue

                    if isinstance(predicted_value, list) and isinstance(observed_value, list):
                        # Pair corresponding list elements.
                        for level, (predicted_item, observed_item) in enumerate(zip(predicted_value, observed_value)):
                            data[level]["predicted"].append(predicted_item)
                            data[level]["observed"].append(observed_item)
                            data[level]["diff"].append(observed_item-predicted_item)
                    else:
                        data["predicted"].append(predicted_value)
                        data["observed"].append(observed_value)
                        data["diff"].append(observed_value-predicted_value)
                        
    correlation_stats_path = sequence_path / f"correlation_stats_{sequence_name}.csv"

    fieldnames = []
    correlation_coefficients = []
    raw_error_medians = []
    absolute_error_medians = []
    squared_error_medians = []
    for stats_class in stats_classes:
        for field_name, data in plot_data[stats_class].items():
            if stats_class == operationTimes_Stats and not field_name == "LookupTime":
                continue
            if get_origin(get_type_hints(stats_class)[field_name]) is list:
                for level in range(NO_LEVELS):
                    predicted_values = data[level]["predicted"]
                    observed_values = data[level]["observed"]
                    differences = data[level]["diff"]
                    fieldnames.append(f"{field_name}{level}")
                    correlation_coefficients.append(np.corrcoef(predicted_values, observed_values)[0, 1])
                    raw_error_medians.append(np.median(differences))
                    absolute_error_medians.append(np.median(np.absolute(differences)))
                    squared_error_medians.append(np.median(np.power(differences, 2)))
            else:
                predicted_values = data["predicted"]
                observed_values = data["observed"]
                differences = data["diff"]
                fieldnames.append(f"{field_name}")
                correlation_coefficients.append(np.corrcoef(predicted_values, observed_values)[0, 1])
                raw_error_medians.append(np.median(differences))
                absolute_error_medians.append(np.median(np.absolute(differences)))
                squared_error_medians.append(np.median(np.power(differences, 2)))
            
    with open(correlation_stats_path, 'w') as csvfile:
        corr_coeff_tex = "$r$"
        rawErrMedian_tex = f"$\\mathbb{{E}}(Pred-Obs)$"
        absErrMedian_tex = f"$\\mathbb{{E}}(|Pred-Obs|)$"
        sqrdErrMedian_tex = f"$\\mathbb{{E}}(Pred-Obs)^2$"
        headers = ["metric", corr_coeff_tex, rawErrMedian_tex, absErrMedian_tex, sqrdErrMedian_tex]
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows([field_name_to_tex[name], *(f"{x:.2e}" for x in (coeff, raw_e_m, abs_e_m, sqrd_e_m))] for name,coeff,raw_e_m,abs_e_m,sqrd_e_m in zip(fieldnames, correlation_coefficients, raw_error_medians, absolute_error_medians, squared_error_medians))
     
    """
    //transpose stats   
    stats_evaluating_space_consumption_path = sequence_path / f"stats_evaluating_space_consumption{sequence_name}.csv"
    with open(correlation_stats_path, newline="") as f:
        source_rows = list(csv.reader(f))
    with open(stats_evaluating_space_consumption_path, "w", newline="") as f:
        csv.writer(f).writerows(list(zip(*source_rows)))
    """
    
    #create plots
    plots_path = genome_plots / "Predicted vs Observed"
    if(not OVERWRITE_PLOTS and plots_path.is_dir()):
        print(f"Predicted vs Observed plots already exists at {plots_path}")
    else:
        plots_path.mkdir(exist_ok=True, parents=True)
        print(f"plotting Predicted vs Observed values plots")
        
        # Create one plot for each field.
        for stats_class in stats_classes:
            for field_name, data in plot_data[stats_class].items():
                if stats_class == operationTimes_Stats and not field_name == "LookupTime":
                    continue
                if get_origin(get_type_hints(stats_class)[field_name]) is list:
                    for level in range(NO_LEVELS):
                        predicted_values = data[level]["predicted"]
                        observed_values = data[level]["observed"]
                        
                        if not predicted_values:
                            continue
                        
                        title = f"{field_name}_{level} for {sequence_name.replace('.','_')}"

                        ax=plt.gca()
                        ax.set_aspect('equal')
                        #plot red diagonal where x=y
                        plt.axline((0,0),slope=1,linestyle='-',color='red',label="predicted=observed")
                        plt.legend()
                        
                        plt.scatter(predicted_values, observed_values)
                        plt.xlabel(f"Predicted {field_name}_{level}")
                        plt.ylabel(f"Observed {field_name}_{level}")
                        plt.title(title)
                        plt.grid(True)
                        plt.tight_layout()
                        plt.savefig(plots_path / f"{title}.png", dpi=300)
                        plt.close()

                else:
                    predicted_values = data["predicted"]
                    observed_values = data["observed"]
                    
                    if not predicted_values:
                        continue

                    title = f"{field_name} for {sequence_name.replace('.','_')}"
                    
                    ax=plt.gca()
                    ax.set_aspect('equal')
                    #plot red diagonal where x=y
                    plt.axline((0,0),slope=1,linestyle='-',color='red',label="predicted=observed")
                    plt.legend()
                    
                    plt.scatter(predicted_values, observed_values)
                    plt.xlabel(f"Predicted {field_name}")
                    plt.ylabel(f"Observed {field_name}")
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