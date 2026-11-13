import csv
from dataclasses import dataclass, field, is_dataclass, fields
from pathlib import Path
from itertools import product
from typing import get_origin, get_type_hints
import re

@dataclass
class spaceConsumption_Stats:
    bpkT: float = field(default_factory=float)
    bpkR: list[float] = field(default_factory=list[float])
    bpkS: list[float] = field(default_factory=list[float])
    bpkO: list[float] = field(default_factory=list[float])
    bpkH: float = field(default_factory=float)
    bpkTotal: float = field(default_factory=float)

@dataclass
class operationTimes_Stats:
    Rcontains: list[float] = field(default_factory=list[float])
    Sselect: list[float] = field(default_factory=list[float])
    Oaccess: list[float] = field(default_factory=list[float])
    Tcompare: list[float] = field(default_factory=list[float])
    Haccess: float = field(default_factory=float)
    LookupTime: float = field(default_factory=float)

field_name_to_tex = {
    "bpkT":"$bpk_T$",
    "bpkH":"$bpk_H$",
    "bpkTotal":f"$bpk_{{total}}$",
    "Haccess":f"$\\timeForOperationOnStructure{{\\check}}{{H}}$",
    "LookupTime":f"$\\timeForOperationOnStructure{{lookup}}{{}}$"
}
for l in range(3):
    field_name_to_tex[f"bpkR{l}"] = f"$bpk_{{R^{l+1}}}$"
    field_name_to_tex[f"bpkS{l}"] = f"$bpk_{{S^{l+1}}}$"
    field_name_to_tex[f"bpkO{l}"] = f"$bpk_{{O^{l+1}}}$"
    field_name_to_tex[f"Rcontains{l}"] = f"$\\timeForOperationOnStructure{{\\contains}}{{R^{l+1}}}$"
    field_name_to_tex[f"Sselect{l}"] = f"$\\timeForOperationOnStructure{{\\select}}{{S^{l+1}}}$"
    field_name_to_tex[f"Oaccess{l}"] = f"$\\timeForOperationOnStructure{{\\access}}{{O^{l+1}}}$"
    field_name_to_tex[f"Tcompare{l}"] = f"$\\timeForOperationOnStructure{{\\compare}}{{T^{l+1}}}$"

@dataclass
class Distribution_Stats:
    n_kmers: int = field(default_factory=int)
    n_minimizers: int = field(default_factory=int)
    F_M: float = field(default_factory=float)
    F_K: float = field(default_factory=float)
    f_avg: float = field(default_factory=float)
    next_level: dict[(int, int), "Distribution_Stats | None"] = field(default_factory=dict)  #links to next stats item depending on (m,t)

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

def readDistribution_StatsFromFiles(stats_path: Path, level: int) -> int | int | list | list | list:                
    n_kmers = int((stats_path / f"{level}n_kmers.txt").read_text().strip())
    n_minimizers = int((stats_path / f"{level}n_minimisers.txt").read_text().strip())

    reader = csv.reader((stats_path / f"{level}kmer_minimizer_distribution.csv").open())
    next(reader)  # Skip the header
    
    think_of_better_name = 0.0
    
    F_M = [0.0]     #initialise cumulative minimizer frequency distribution
    F_K = [0.0]     #initialise cumulative k-mer (minimizer) frequency distribution
    f_avg = [0.0]   #initialise average frequency of k-mers
    for f, row in enumerate(reader, start=1):
        think_of_better_name += float(row[1])
        F_M.append(F_M[-1] + float(row[1])*f)
        F_K.append(F_K[-1] + float(row[2]))
        f_avg.append(F_M[-1]/think_of_better_name)
        
    return n_kmers, n_minimizers, F_K, F_M, f_avg

def fill_cascading_stats(
    root_stats: Distribution_Stats, 
    genome_stats_path: Path, 
    level: int, 
    maxLevel: int, 
    currentMinimizerLengths: list[int],
    currentThresholds: list[int],
    minimizerLengthsConsidered: list[int], 
    thresholdsConsidered: list[int],
    nonDecreasingParameters=False) -> None:
    
    if(level>maxLevel):
        return 
    for m in minimizerLengthsConsidered:
        if(nonDecreasingParameters and level!=1 and m<currentMinimizerLengths[level-2]):
            continue
        currentMinimizerLengths[level-1] = m
        stats_path =  genome_stats_path / getParameterDirectory(currentMinimizerLengths, currentThresholds)    #directory for a specific set of level parameter
        n_kmers, n_minimizers, F_K, F_M, f_avg = readDistribution_StatsFromFiles(stats_path, level)
        for t  in thresholdsConsidered:
            if(nonDecreasingParameters and level!=1 and t<currentThresholds[level-2]):
                        continue
            currentThresholds[level-1] = t
            current_stats = Distribution_Stats(n_kmers, 
                                n_minimizers, 
                                F_M[t], 
                                F_K[t], 
                                f_avg[t], 
                                next_level={None if level == maxLevel 
                                            else None if (nonDecreasingParameters and level!=1 and (m2<m or t2<t))
                                                else (m2,t2): None for m2,t2 in product(minimizerLengthsConsidered, thresholdsConsidered)})
            
            root_stats.next_level[(m,t)] = current_stats
            fill_cascading_stats(current_stats, 
                                genome_stats_path, 
                                level+1, 
                                maxLevel,
                                currentMinimizerLengths,
                                currentThresholds,
                                minimizerLengthsConsidered,
                                thresholdsConsidered,
                                nonDecreasingParameters)

        
def getCascadingDistribution_Stats(cascading_stats: Distribution_Stats, minimizer_lengths: list[int], thresholds: list[int]) -> Distribution_Stats:
    stats = cascading_stats
    for m, t in zip(minimizer_lengths, thresholds):
        stats = stats.next_level[(m,t)]
    return stats.n_kmers, stats.n_minimizers, stats.F_M, stats.F_K, stats.f_avg

def column_values(rows: list[dict[str, str]], column: str) -> list[float]:
    return [float(row[column]) for row in rows]

def getObservedStats(genome_stats_path: Path,
                    minimizer_lengths: list[int],
                    thresholds: list[int]) -> tuple[spaceConsumption_Stats, operationTimes_Stats]:
    stats_path =  genome_stats_path / getParameterDirectory(minimizer_lengths, thresholds) / "space_consumption.csv"
    with stats_path.open() as f:
        space_rows = list(csv.DictReader(f, skipinitialspace=True))
        space = spaceConsumption_Stats(
            bpkT=float(space_rows[0]["bpk_text"]),
            bpkR=column_values(space_rows, "bpk_R"),
            bpkS=column_values(space_rows, "bpk_S"),
            bpkO=column_values(space_rows, "bpk_O"),
            bpkH=float(space_rows[0]["bpk_H"]),
            bpkTotal=float(space_rows[0]["bpk_total"])
        )
    
    stats_path =  stats_path.parent / "times_per_operation.csv"
    with stats_path.open() as f:
        time_rows = list(csv.DictReader(f, skipinitialspace=True))
        times = operationTimes_Stats(
            Rcontains=column_values(time_rows, "R_contains"),
            Sselect=column_values(time_rows, "S_select"),
            Oaccess=column_values(time_rows, "O_access"),
            Tcompare=column_values(time_rows, "T_compare"),
            Haccess=float(time_rows[0]["H_access"]),
        )
    
    lookup_log =  stats_path.parent / "lookup.log"
    match=re.search(r"^time_per_kmer = ((?:\d+(?:\.\d*)?|\.\d+))", lookup_log.read_text(), re.MULTILINE)  #regular expression matches with the text after "time_per_kmer = "
    times.LookupTime = float(match.group(1))   #match group retrieves the text matching with the regular expression
    
    return space, times

##############
# helpers to write the array predicted_AND_observed_stats to a json file / recover it from a json file
# Generative AI use Disclaimer:
#   the following code was generated using GPT-6 Luna
##############
STATS_TYPES = {
    "spaceConsumption_Stats": spaceConsumption_Stats,
    "operationTimes_Stats": operationTimes_Stats,
}

def predicted_AND_observed_stats_to_jsonable(value):
    if is_dataclass(value):
        return {
            "__type__": type(value).__name__,
            **{f.name: predicted_AND_observed_stats_to_jsonable(getattr(value, f.name)) for f in fields(value)}
        }
    if isinstance(value, (list, tuple)):
        return [predicted_AND_observed_stats_to_jsonable(item) for item in value]
    return value

def predicted_AND_observed_stats_from_jsonable(value):
    if isinstance(value, list):
        return [predicted_AND_observed_stats_from_jsonable(item) for item in value]

    if isinstance(value, dict):
        decoded = {key: predicted_AND_observed_stats_from_jsonable(item) for key, item in value.items()}
        type_name = decoded.pop("__type__", None)

        if type_name is not None:
            return STATS_TYPES[type_name](**decoded)
        return decoded

    return value
