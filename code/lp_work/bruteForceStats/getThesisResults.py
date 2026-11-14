from pathlib import Path
import csv
from itertools import product
import matplotlib.pyplot as plt
import itertools

import numpy as np
from scipy.stats import gaussian_kde

from optimization.cascadingStats import *

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

    plotData = []
    for stats_field in fields(operationTimes_Stats):
        if stats_field.name == "LookupTime":
            continue
        if get_origin(get_type_hints(operationTimes_Stats)[stats_field.name]) is list:
            plotData.append(list(itertools.chain.from_iterable(allTimes[stats_field.name])))
        else: plotData.append(allTimes[stats_field.name])
    
    writeBoxPlotTex(plotLabels, plotData, (thesis_results_path / "combinedBoxPlots.tex"))
        
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
            
    writeBoxPlotTex(plotLabels, plotData, (thesis_results_path / "RcontainsBoxPlot.tex"), colours=True)

    # Generative "AI" use Disclaimer:
    #   the following function was generated using GPT-6 Luna
def writeModelErrorOverview(result_path: Path, columns=1):
    variables = [
        r"$\timeForOperationOnStructure{lookup}{}$",
        r"$bpk_{total}$",
    ]

    headers = [
        "Dataset",
        r"$\NAMEcorrCoef$",
        r"$\NAMEmeanErrRaw$",
        r"$\NAMEmeanErrAbs$",
        r"$\NAMEmeanErrSqrd$",
    ]

    # Load each genome's CSV rows once.
    rows_by_genome = {}
    for genome in GENOMES:
        csv_path = Path(
            f"./optimization/evaluation_results/correlation_stats_{genome}.csv"
        )
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

    #separator = "\n\\hfill\n" if columns > 1 else "\n\\par\\medskip\n"
    separator = "\n"
    
    ##only necessary, if full figure is wanted
    #body = separator.join(tables)
    #tex_file_content = "\n".join([
    #    r"\begin{figure}[htbp]",
    #    r"\centering",
    #    body,
    #    r"\end{figure}",
    #])

    with (result_path / "modelErrorOverview.tex").open("w", encoding="utf-8") as f:
        f.write(separator.join(tables))
def main() -> None:
    thesis_results_path = Path("thesis_results")
    thesis_results_path.mkdir(parents=True, exist_ok=True)
    
    stats_directory = Path("allStats")  #directory for stats of all genomes
    if(not stats_directory.is_dir()):
        return
    #getAllOperationTimes(stats_directory, thesis_results_path)
    
    writeModelErrorOverview(thesis_results_path, columns=1)

if __name__ == '__main__':
    main()


