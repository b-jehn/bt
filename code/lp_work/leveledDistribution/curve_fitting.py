from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expit
from scipy.optimize import curve_fit

genome = "Gadus_morhua.gadMor3.0.dna.toplevel"
LEVELS = 3
threshold_sets = [[64,64,64]]
minimizerLengths_sets = [[m,m,m] for m in range(10,20)]

def logistic1(x, a, b, c):
    return a/(1+b*np.exp(c*x))  

#numerically stable equivalent to logistic1
def logistic2(x, a, b, c):
    return a * expit(-(c * x + np.log(b)))

parameter_dict = {
  0: 'a',
  1: 'b',
  2: 'c'
}

"""
t = threshold_sets[0]
m = minimizerLengths_sets[8]
stats_directory = Path(f"./{genome}/m={m[0]},{m[1]},{m[2]}_t={t[0]},{t[1]},{t[2]}")
stats_path = stats_directory / f"distinctStats_level{1}.csv"

stats = pd.read_csv(stats_path)
# Use the first column as x and the second column as y
x = np.array(stats.iloc[:, 0], dtype=np.float64) 
y = np.array(stats.iloc[:, 1], dtype=np.float64) 

parameters, covariance = curve_fit(logistic_func, x, y, bounds=(-1,1))

print(parameters)
"""

all_parameters = [[] for l in range(LEVELS)] 

for t in threshold_sets:
    #curve fitting for all specified minimizer length sets
    for m_index, m in enumerate(minimizerLengths_sets):
        stats_directory = Path(f"./{genome}/m={m[0]},{m[1]},{m[2]}_t={t[0]},{t[1]},{t[2]}")
        
        print(f"optimizing for m={m}")
        for level in range(1, LEVELS+1):
            stats_path = stats_directory / f"distinctStats_level{level}.csv"

            stats = pd.read_csv(stats_path)

            # Use the first column as x and the second column as y
            x = np.array(stats.iloc[:, 0], dtype=np.float64) 
            y = np.array(stats.iloc[:, 1], dtype=np.float64) 
            
            parameters, covariance = curve_fit(logistic2, x, y, maxfev=10000)
            
            all_parameters[level-1].append(parameters)
            print(f"level{level}: a={all_parameters[level-1][m_index][0]}, b={all_parameters[level-1][m_index][1]}, c={all_parameters[level-1][m_index][2]}")
            
    print("fitted data for all levels to logistic function: a/(1+b*np.exp(c*x))\n with parameters:")
    for i, m in enumerate(minimizerLengths_sets):
        print(f"m={m}:")
        for level in range(LEVELS):
            print(f"level{level+1}: a={all_parameters[level][i][0]}, b={all_parameters[level][i][1]}, c={all_parameters[level][i][2]}")
            

    #plot parameters across levels
    plt.figure(figsize=(8, 6))
    x = [m[0] for m in minimizerLengths_sets]
    for level in range(LEVELS):
        for parameter in range(3):
            # Use the first column as x and the second column as y
            y = [all_parameters[level][m][parameter] for m in range(len(minimizerLengths_sets))]

            plt.plot(x, y, "o-", label=f"{parameter_dict[parameter]} on level {level+1}, \n")

    plt.xlabel("minimizer lengths on all levels")
    plt.ylabel("parameter value")
    plt.title("logistic regression parameter values")
    plt.legend()
    plt.grid(True, alpha=0.3)
    #plt.ylim(0, 100)
    plt.tight_layout()

    #save the plot
    #plt.savefig("logistic regression", dpi=300)
    plt.show()
    plt.close()


