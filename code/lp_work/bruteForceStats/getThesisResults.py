from pathlib import Path
import csv, sys, os
from itertools import product
import matplotlib.pyplot as plt
import itertools
from operator import itemgetter 

import numpy as np
from scipy.stats import gaussian_kde

from optimization.cascadingStats import *
#from optimization.MINLP import *


DATA_DIRECTORY = Path(f"../../rshash/DNA_datasets")
N_SAMPLES = 100000 
K = 31

GENOMES = ["Falco_tinnunculus.FalTin1.0.dna.toplevel",
            "Falco_tinnunculus.FalTin1.0.dna.toplevel.fa.unitigs",
            "Gadus_morhua.gadMor3.0.dna.toplevel",
            "Gadus_morhua.gadMor3.0.dna.toplevel.fa.unitigs",
            "Homo_sapiens.GRCh38.dna.toplevel",
            "Homo_sapiens.GRCh38.dna.toplevel.fa.unitigs"]

GENOME_NAMES = ["Kestrel (Raw)",
                "Kestrel (Unitigs)",
                "Cod (Raw)",
                "Cod (Unitigs)",
                "Human (Raw)",
                "Human (Unitigs)",]


MINIMIZER_LENGTHS_CONSIDERED = [m for m in range(15,26,2)]
THRESHOLDS_CONSIDERED = [2**i for i in range(5,8)]  

NON_DECREASING_PARAMETERS = True

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

NO_LEVELS = 3

def getParameterDirectory(minimizer_lengths: list[int], thresholds: list[int]) -> Path:
    no_levels = len(minimizer_lengths)
    m_part = "m="
    t_part = "t="
    for i, (m, t) in enumerate(zip(minimizer_lengths, thresholds), start=1):
        m_part += f"{m}"
        t_part += f"{t}"
        if(i != no_levels):
            m_part += f","
            t_part += f","
    return Path(f"{m_part}_{t_part}")

def writeBoxPlotTex(plotLabels: list, plotData: list, result_path: Path, colours=False) -> None:
    with open(result_path, "w") as f:
        #write latex code for a figure which contains a box plot for each plotLable that plots the corresponding data
        f.write(r"""\begin{tikzpicture}
            \begin{axis}[
                boxplot/draw direction=y,
                xtick={1,2,3,4,5},
                xticklabels={""" + ",".join("{" + label + "}" for label in plotLabels) + r"""},
                x tick label style={align=center},
                ylabel={Measured Time (ns)},
            ]
            """)

        for position, data in enumerate(plotData, start=1):
            values = np.asarray(data, dtype=float)
            q1, median, q3 = np.percentile(values, [25, 50, 75])
            iqr = q3 - q1
            lower_whisker = values[values >= q1 - 1.5 * iqr].min()
            upper_whisker = values[values <= q3 + 1.5 * iqr].max()

            f.write(
                rf"\addplot+"
                rf" [{'' if colours else 'draw=black, '}boxplot prepared={{"
                rf"draw position={position}, "
                rf"lower whisker={lower_whisker}, lower quartile={q1}, "
                rf"median={median}, upper quartile={q3}, "
                rf"upper whisker={upper_whisker}"
                rf"}}] coordinates {{}};" + "\n"
            )

        f.write(r"""\end{axis}
            \end{tikzpicture}
            """)

def getAllOperationTimes(stats_directory: Path, thesis_results_path: Path):
    print(f"getting information needed for latex boxplots")
    allTimes = {
        stats_field.name : ([[] for _ in range(NO_LEVELS)] 
                    if get_origin(get_type_hints(operationTimes_Stats)[stats_field.name]) is list
                    else [])
        for stats_field in fields(operationTimes_Stats)
    }
    
    for genome in GENOMES:
        print(f"getting stats for: {genome}")
        genome_stats_path = stats_directory / f"{genome}"  #directory for the stats of a specific genome
        if(not genome_stats_path.is_dir()):               
            continue
        
        for minimizers, thresholds in product(setOfMinimizerLengths, setOfThresholds):
            _, times = getObservedStats(genome_stats_path, minimizers,thresholds)
            for stats_field in fields(times):
                name = stats_field.name
                attr = getattr(times, name)
                if isinstance(attr, list):
                    for level in range(NO_LEVELS):
                        allTimes[name][level].append(attr[level])
                else:
                    allTimes[name].append(attr)
                    
    plotLabels = [
        f"$\\timeForOperationOnStructure{{\\contains}}{{R}}$",
        f"$\\timeForOperationOnStructure{{\\select}}{{S}}$",
        f"$\\timeForOperationOnStructure{{\\access}}{{O}}$",
        f"$\\timeForOperationOnStructure{{\\compare}}{{T}}$",
        f"$\\timeForOperationOnStructure{{check}}{{H}}$"
    ]

    #plot containing one boxplot for each operation time across all levels
    plotData = []
    for stats_field in fields(operationTimes_Stats):
        if stats_field.name == "LookupTime":
            continue
        if get_origin(get_type_hints(operationTimes_Stats)[stats_field.name]) is list:
            plotData.append(list(itertools.chain.from_iterable(allTimes[stats_field.name])))
        else: plotData.append(allTimes[stats_field.name])
    
    writeBoxPlotTex(plotLabels, plotData, (thesis_results_path / "boxPlotsCombined.tex"))
    
    #plot for Rcontains on level one, two and three
    plotLabels = [
        f"$\\timeForOperationOnStructure{{\\contains}}{{R^1}}$",
        f"$\\timeForOperationOnStructure{{\\contains}}{{R^2}}$",
        f"$\\timeForOperationOnStructure{{\\contains}}{{R^3}}$"
    ]

    plotData = []
    for stats_field in fields(operationTimes_Stats):
        if not stats_field.name == "Rcontains":
            continue
        for level in range(NO_LEVELS):
            plotData.append(allTimes[stats_field.name][level])
            
    writeBoxPlotTex(plotLabels, plotData, (thesis_results_path / "boxPlotsRcontains.tex"), colours=True)

    # Generative "AI" use Disclaimer:
    #   the following function was generated using GPT-6 Luna
def writeModelErrorOverview(result_path: Path, variables: list, columns=1, withMeasuredOperationTimes=False):        
    #headers = [
    #    "Dataset",
    #    r"$\NAMEcorrCoef$",
    #    r"$\NAMEmeanErrRaw$",
    #    r"$\NAMEmeanErrAbs$",
    #    r"$\NAMEmeanErrSqrd$",
    #]
    headers = [
        "Dataset",
        f"$\\frac{{cov(Pred,Obs)}}{{\\sigma_{{Pred}} \\sigma_{{Obs}}}}$",
        f"$\\mathbb{{E}}(Pred-Obs)$",
        f"$\\mathbb{{E}}(|Pred-Obs|)$",
        f"$\\mathbb{{E}}(Pred-Obs)^2$",
    ]

    # Load each genome's CSV rows once.
    stats_prefix=""
    if withMeasuredOperationTimes:
        stats_prefix="WITH_MEASURED_OPERATION_TIMES_"
        
    rows_by_genome = {}
    for genome in GENOMES:
        csv_path = Path(
            f"./optimization/{stats_prefix}evaluation_results/correlation_stats_{genome}.csv"
        )
        print(f"getting stats from: {csv_path}")
        with csv_path.open(newline="", encoding="utf-8") as f:
            rows_by_genome[genome] = {
                row["variable"]: row
                for row in csv.DictReader(f)
            }

    tables = []
    tables.append("\n".join([
        r"\begin{tabular}{lrrrr}",
        r"\toprule",
    ]))
    for variable in variables:
        lines=[]
        if withMeasuredOperationTimes:
            lines = [
                " & ".join(headers) + r" \\",
                r"\midrule",
            ]
        else:
            lines = [
                rf"\textbf{{{variable}}}" + " &" * (len(headers) - 1) + r" \\",
                r"\midrule",
                " & ".join(headers) + r" \\",
                r"\midrule",
            ]

        for dataset, genome in zip(GENOME_NAMES, GENOMES):
            rows = rows_by_genome[genome]
            if variable not in rows:
                raise ValueError(
                    f"Variable {variable!r} not found for {genome}"
                )

            row = rows[variable]
            r_value = row[r"$r$"]
            if r_value.lower() == "nan":
                r_value = "--"

            mean_error = row[r"$\mathbb{E}(Pred-Obs)$"]
            mean_abs_error = row[r"$\mathbb{E}(|Pred-Obs|)$"]
            mean_squared_error = row[r"$\mathbb{E}(Pred-Obs)^2$"]

            lines.append(
                f"{dataset} & {r_value} & {mean_error} & {mean_abs_error} & {mean_squared_error} \\\\"
            )
            
        if not variable == variables[-1]:
            lines.append(r"\midrule")
            lines.append(r"\midrule")
        
        tables.append("\n".join(lines))

    tables.append("\n".join([
        r"\bottomrule",
        r"\end{tabular}",
    ]))

    with (result_path).open("w", encoding="utf-8") as f:
        f.write("\n".join(tables))
        
def writeSortedLookupTimes(stats_directory: Path, result_path: Path, start_slice_genomes=0, end_slice_genomes=[len(GENOMES)]):
    sortedLookupTimes = {}
    
    for genome in GENOMES[start_slice_genomes:end_slice_genomes]:
        print(f"getting lookup times for: {genome}")
        genome_stats_path = stats_directory / f"{genome}"  #directory for the stats of a specific genome
        if(not genome_stats_path.is_dir()):               
            continue
        sequence_path = f'../../rshash/DNA_datasets/{genome}.fa'
        if not os.path.isfile(sequence_path):
            sys.exit(f"{sequence_path} is not a file")
        log_sequence_length = np.log2(os.path.getsize(sequence_path))
        
        genome_stats = Distribution_Stats(n_kmers=0, n_minimizers=0, F_M=0, F_K=0, f_avg=0, next_level={(m,t): None for m,t in product(MINIMIZER_LENGTHS_CONSIDERED, THRESHOLDS_CONSIDERED)})  
        fill_cascading_stats(root_stats=genome_stats, 
                            genome_stats_path=genome_stats_path, 
                            level=1, 
                            maxLevel=NO_LEVELS,
                            currentMinimizerLengths=[MINIMIZER_LENGTHS_CONSIDERED[-1]] * NO_LEVELS,
                            currentThresholds=[THRESHOLDS_CONSIDERED[-1]] * NO_LEVELS,
                            minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED,
                            thresholdsConsidered=THRESHOLDS_CONSIDERED,
                            nonDecreasingParameters=True)
                    
        genomeLookupTimes = []
        for forcedMinimizerLengths, forcedThresholds in product(setOfMinimizerLengths, setOfThresholds):
            observedSpace, observedTimes = getObservedStats(genome_stats_path, 
                                                forcedMinimizerLengths,
                                                forcedThresholds)
            
            print(f"running MINLP with forced minimizer lengths: {forcedMinimizerLengths} and thresholds: {forcedThresholds}")
            predictedSpace, predictedTimes, info, solution, elapsed_seconds, optimal_parameters = MINLP(
                optimization_goal="t_lookup",
                max_lookup_time=np.inf,
                max_bpk_total=np.inf,
                log_sequence_length=log_sequence_length,
                no_levels=NO_LEVELS, 
                minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED, 
                thresholdsConsidered=THRESHOLDS_CONSIDERED, 
                distribution_stats=genome_stats,
                load_factor=0.75,
                nonDecreasingParameters=NON_DECREASING_PARAMETERS,
                forceParameters=[forcedMinimizerLengths,forcedThresholds])

            genomeLookupTimes.append([forcedMinimizerLengths, forcedThresholds, predictedTimes.LookupTime, observedTimes.LookupTime])

        sortedLookupTimes[genome] = sorted(genomeLookupTimes, key=itemgetter(2), reverse=True)

        print(f"writing sortedLookupTimes_{genome}.csv")
        with open((result_path / f"sortedPredictedLookupTimes_{genome}.csv"), "w") as f:
            headers = ["m1","m2","m3","t1","t2","t3","predictedLookupTime", "observedLookupTime"]
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows([*minimizer_lengths, *thresholds, predLookupTime, obsLookupTime] for minimizer_lengths,thresholds,predLookupTime,obsLookupTime in sortedLookupTimes[genome])
        
def writeHashSetSizes(stats_directory: Path, result_path: Path, start_slice_genomes=0, end_slice_genomes=[len(GENOMES)]):
    for genome in GENOMES[start_slice_genomes:end_slice_genomes]:
        print(f"getting bpk_H for: {genome}")
        genome_stats_path = stats_directory / f"{genome}"  #directory for the stats of a specific genome
        if(not genome_stats_path.is_dir()):               
            continue
        sequence_path = f'../../rshash/DNA_datasets/{genome}.fa'
        if not os.path.isfile(sequence_path):
            sys.exit(f"{sequence_path} is not a file")
        log_sequence_length = np.log2(os.path.getsize(sequence_path))
        
        genome_stats = Distribution_Stats(n_kmers=0, n_minimizers=0, F_M=0, F_K=0, f_avg=0, next_level={(m,t): None for m,t in product(MINIMIZER_LENGTHS_CONSIDERED, THRESHOLDS_CONSIDERED)})  
        fill_cascading_stats(root_stats=genome_stats, 
                            genome_stats_path=genome_stats_path, 
                            level=1, 
                            maxLevel=NO_LEVELS,
                            currentMinimizerLengths=[MINIMIZER_LENGTHS_CONSIDERED[-1]] * NO_LEVELS,
                            currentThresholds=[THRESHOLDS_CONSIDERED[-1]] * NO_LEVELS,
                            minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED,
                            thresholdsConsidered=THRESHOLDS_CONSIDERED,
                            nonDecreasingParameters=True)
                    
        hashSetSizes = []
        for forcedMinimizerLengths, forcedThresholds in product(setOfMinimizerLengths, setOfThresholds):
            observedSpace, observedTimes = getObservedStats(genome_stats_path, 
                                                forcedMinimizerLengths,
                                                forcedThresholds)
            
            print(f"running MINLP with forced minimizer lengths: {forcedMinimizerLengths} and thresholds: {forcedThresholds}")
            predictedSpace, predictedTimes, info, solution, elapsed_seconds, optimal_parameters = MINLP(
                optimization_goal="t_lookup",
                max_lookup_time=np.inf,
                max_bpk_total=np.inf,
                log_sequence_length=log_sequence_length,
                no_levels=NO_LEVELS, 
                minimizerLengthsConsidered=MINIMIZER_LENGTHS_CONSIDERED, 
                thresholdsConsidered=THRESHOLDS_CONSIDERED, 
                distribution_stats=genome_stats,
                load_factor=0.75,
                nonDecreasingParameters=NON_DECREASING_PARAMETERS,
                forceParameters=[forcedMinimizerLengths,forcedThresholds])

            hashSetSizes.append([forcedMinimizerLengths[0], predictedSpace.bpkH, observedSpace.bpkH])

        print(f"writing hashSetSizes_{genome}.csv")
        with open((result_path / f"hashSetSizes_{genome}.csv"), "w") as f:
            headers = ["m1", "PredBpkH", "ObsBpkH"]
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(row for row in hashSetSizes)
            
def convert_csv_files(thesis_result_paths):
    LATEX_HEADERS = [
        r"$\optGoal{}$",
        r"$\minimizerLength^1$",
        r"$\minimizerLength^2$",
        r"$\minimizerLength^3$",
        r"$\threshold^1$",
        r"$\threshold^2$",
        r"$\threshold^3$",
        r"$bpk_{total}$",
        r"$\timeForOperationOnStructure{lookup}{}$",
    ]
    for genome, genome_name in zip(GENOMES, GENOME_NAMES):
        source_path = Path(f"./optimization/optimizationResults/optimizationResults_{genome}.csv")
        result_path = thesis_result_paths / f"optimizationResults_{genome}.csv"
        with open(source_path, newline="", encoding="utf-8") as source_file:
            rows = csv.reader(source_file)
            next(rows, None)  # Discard the source header

            output_rows = []
            for i, row in enumerate(rows, start=1):
                if len(row) < 11:
                    raise ValueError(f"Unexpected row in {source_path}: {row}")

                optGoal = i
                if i == 3:
                    optGoal = "fast"
                if i == 4:
                    optGoal = "compact"
                # Keep the minimizer lengths, thresholds, and predicted values.
                output_rows.append([
                    rf"$\optGoal{{{optGoal}}}$",
                    *row[3:9],
                    f"{float(row[9]):.4f}",
                    f"{float(row[10]):.4f}",
                ])
                
        Path(result_path).parent.mkdir(parents=True, exist_ok=True)
        with open(result_path, "w", newline="", encoding="utf-8") as result_file:
            writer = csv.writer(result_file)
            writer.writerow(LATEX_HEADERS)
            writer.writerows(output_rows)
            
import csv

def make_parameter_table(path_kestrel, path_cod, path_human, result_path):
    """Return LaTeX for a table of the fast and compact parameters in three CSVs."""

    genomes = [
        ("Kestrel", path_kestrel),
        ("Cod", path_cod),
        ("Human", path_human),
    ]
    m_columns = [
        r"$\minimizerLength^1$",
        r"$\minimizerLength^2$",
        r"$\minimizerLength^3$",
    ]
    t_columns = [
        r"$\threshold^1$",
        r"$\threshold^2$",
        r"$\threshold^3$",
    ]

    # Read the rows for each optimization goal from each CSV.
    data = {}
    for genome, path in genomes:
        with open(path, newline="", encoding="utf-8") as csvfile:
            print("file opened")
            rows = csv.DictReader(csvfile)
            data[genome] = {
                row[r"$\optGoal{}$"]: row
                for row in rows
            }

    lines = [
        r"\begin{table}[h]",
        r"    \centering",
        r"    \begin{tabular}{lc|lll}",
        r"        \toprule",
        r"        $\optGoal{}$ & & Kestrel & Cod & Human \\",
        r"        \midrule",
    ]

    for goal in ("compact", "fast"):
        for parameter, columns in (("m", m_columns), ("t", t_columns)):
            cells = []
            for genome, _ in genomes:
                row = data[genome][rf"$\optGoal{{{goal}}}$"]
                cells.append(",".join(row[column] for column in columns))

            version_cell = rf"$\optGoal{{{goal}}}$" if parameter == "m" else ""
            lines.append(
                f"        {version_cell} & {parameter} & "
                + " & ".join(cells)
                + r" \\"
            )

        if goal == "compact":
            lines.append(r"        \midrule")

    lines.extend([
        r"        \bottomrule",
        r"    \end{tabular}",
        r"    \caption{Parameter configurations for the fast or compact version.}",
        r"    \label{tab:MINLPOptimalParameters}",
        r"\end{table}",
    ])

    with (result_path).open("w", encoding="utf-8") as f:
        f.write("\n".join(lines))
        

def main() -> None:
    thesis_results_path = Path("thesis_results")
    thesis_results_path.mkdir(parents=True, exist_ok=True)
    
    stats_directory = Path("allStats")  #directory for stats of all genomes
    if(not stats_directory.is_dir()):
        return
    #getAllOperationTimes(stats_directory, thesis_results_path)
    
    #variables = [
    #    r"$\timeForOperationOnStructure{lookup}{}$",
    #    r"$bpk_{total}$",
    #]
    #writeModelErrorOverview((thesis_results_path / "modelErrorOverview.tex"), variables, columns=1)
    #
    #variables = [
    #    r"$\timeForOperationOnStructure{lookup}{}$",
    #]
    #writeModelErrorOverview((thesis_results_path / "WITH_MEASURED_OPERATION_TIMES_lookupTimeOverview.tex"), variables, columns=1, withMeasuredOperationTimes=True)
    #
    #variables = [
    #    "$bpk_T$",
    #    "$bpk_H$",
    #]
    #writeModelErrorOverview((thesis_results_path / "bpk_T_H_Overview.tex"), variables, columns=1)
    
    #writeSortedLookupTimes(stats_directory, thesis_results_path, start_slice_genomes=4)
    
    #writeHashSetSizes(stats_directory, thesis_results_path, start_slice_genomes=2, end_slice_genomes=4)
    
    #convert_csv_files(thesis_results_path)
    make_parameter_table(path_kestrel=Path(f"thesis_results/optimizationResults_{GENOMES[0]}.csv"), 
                        path_cod=Path(f"thesis_results/optimizationResults_{GENOMES[2]}.csv"), 
                        path_human=Path(f"thesis_results/optimizationResults_{GENOMES[4]}.csv"),
                        result_path=Path("thesis_results/modeledParameterChoices.tex"))

if __name__ == '__main__':
    main()


