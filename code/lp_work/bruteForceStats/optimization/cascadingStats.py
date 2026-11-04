import csv
from dataclasses import dataclass, field
from pathlib import Path
from itertools import product

@dataclass
class spaceConsumption_Stats:
    bpk_T: float = field(default_factory=float)
    bpk_R: list[float] = field(default_factory=list[float])
    bpk_S: list[float] = field(default_factory=list[float])
    bpk_O: list[float] = field(default_factory=list[float])
    bpk_H: float = field(default_factory=float)
    bpk_total: float = field(default_factory=float)

@dataclass
class operationTimes_Stats:
    R_contains: list[float] = field(default_factory=list[float])
    S_select: list[float] = field(default_factory=list[float])
    O_access: list[float] = field(default_factory=list[float])
    T_compare: list[float] = field(default_factory=list[float])
    H_access: float = field(default_factory=float)

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