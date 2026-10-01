from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expit
from scipy.optimize import curve_fit, least_squares, minimize

genome = "Gadus_morhua.gadMor3.0.dna.toplevel"
LEVELS = 3
threshold_sets = [[64,64,64]]
minimizerLengths_sets = [[m,m,m] for m in range(10,20)]

"""
best solution so far but the curve is monotonous
"""
def loglogistic4P(x, params):
    """
    Four-parameter log-logistic curve.

    Parameters
    ----------
    x : array-like
        Positive predictor values.
    params : array-like
        [lower, upper, log_alpha, log_beta]

    Returns
    -------
    y : ndarray
        Predicted response values.
    """
    x = np.asarray(x, dtype=float)

    lower, upper, log_alpha, log_beta = params
    alpha = np.exp(log_alpha)
    beta = np.exp(log_beta)

    z = beta * (np.log(x) - np.log(alpha))

    # expit(z) is numerically stable even for very large |z|
    return lower + (upper - lower) * expit(z)

def fit_loglogistic4P(x, y):
    """
    Fit a four-parameter log-logistic curve by robust nonlinear least squares.

    Returns
    -------
    result : scipy.optimize.OptimizeResult
        Fitted result. The parameters are:
        result.x = [lower, upper, log_alpha, log_beta]
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be one-dimensional arrays.")

    if x.size != y.size:
        raise ValueError("x and y must have the same length.")

    if x.size < 4:
        raise ValueError("At least four observations are needed.")

    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("x and y must contain only finite values.")

    if np.any(x <= 0):
        raise ValueError("All x values must be strictly positive.")

    # Initial estimates
    y_min = np.min(y)
    y_max = np.max(y)
    y_range = max(y_max - y_min, np.finfo(float).eps)

    lower0 = y_min
    upper0 = y_max

    # Estimate alpha from the x-value near the middle response
    midpoint = lower0 + 0.5 * (upper0 - lower0)
    alpha0 = x[np.argmin(np.abs(y - midpoint))]

    # Reasonable initial slope
    beta0 = 1.0

    p0 = np.array([
        lower0,
        upper0,
        np.log(alpha0),
        np.log(beta0),
    ])

    def residuals(params):
        prediction = loglogistic4P(x, params)
        return prediction - y

    # log_alpha and log_beta are unbounded, while the asymptotes
    # are constrained to remain within a broad range around the data.
    y_scale = max(np.ptp(y), 1.0)

    lower_bounds = np.array([
        y_min - 10.0 * y_scale,
        y_min - 10.0 * y_scale,
        -np.inf,       # log(alpha)
        -20.0,         # beta >= exp(-20)
    ])

    upper_bounds = np.array([
        y_max + 10.0 * y_scale,
        y_max + 10.0 * y_scale,
        np.inf,        # log(alpha)
        20.0,          # beta <= exp(20)
    ])

    result = least_squares(
        residuals,
        x0=p0,
        bounds=(lower_bounds, upper_bounds),
        loss="soft_l1",       # robust against outliers
        f_scale=y_scale,
        max_nfev=20_000,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )

    return result.x


def plotDataAgainstModel(x_data, y_data, parameters, minimizer_lengths, thresholds, level):
    #plt.figure(figsize=(8, 6))
    #y_learned = [logistic(x, parameters[0], parameters[1], parameters[2]) for x in x_data]
    y_learned = [loglogistic4P(x, parameters) for x in x_data]
    #y_learned = [loglogistic_pdf(x, parameters) for x in x_data]
    
    plt.plot(x_data, y_data, "xr", label="measured")
    plt.plot(x_data, y_learned, "o-b", label="learned")

    plt.xlabel("frequency")
    plt.ylabel("percentage")
    plt.title(f"level={level}, m={minimizer_lengths}, t={thresholds}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    #plt.ylim(0, 100)
    plt.tight_layout()

    #save the plot
    plt.savefig(f"./loglogistic4PRegression/level{level}/m={minimizer_lengths}, t={thresholds}", dpi=300)
    #plt.show()
    plt.close()

parameter_dict = {
  0: "lower_asmpt",
  1: "upper_asmpt",
  2: "log_alpha",
  3: "log_beta"
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
            x = np.asarray(stats.iloc[:, 0], dtype=np.float64) 
            y = np.asarray(stats.iloc[:, 1], dtype=np.float64) 
            
            parameters = fit_loglogistic4P(x,y)
            all_parameters[level-1].append(parameters)
            plotDataAgainstModel(x,y,parameters, m,t, level)
            print(f"level{level}: a={all_parameters[level-1][m_index][0]}, b={all_parameters[level-1][m_index][1]}")
            residuals = loglogistic4P(x, parameters) - y
            print(f"residuals={np.sum(residuals**2)}")

    """ 
    print("fitted data for all levels to logistic function: a * expit(-(c * x + np.log(b))) = a/(1+b*np.exp(c*x))\n with parameters:")
    for i, m in enumerate(minimizerLengths_sets):
        print(f"m={m}:")
        for level in range(LEVELS):
            print(f"level{level+1}: a={all_parameters[level][i][0]}, b={all_parameters[level][i][1]}, c={all_parameters[level][i][2]}")
    """            

    #plot parameters across levels
    mWindowStart = 4    #default = 0, set to x to ignore the first x minimizer lengths
    plt.figure(figsize=(16, 12))
    x = [m[0] for m in minimizerLengths_sets[mWindowStart:]]
    
    #plot every parameter a for all levels
    for parameter in range(4):
        plt.subplot(2, 2, parameter+1)
        for level in range(LEVELS):
                # Use the first column as x and the second column as y
                y = [all_parameters[level][m+mWindowStart][parameter] for m in range(len(minimizerLengths_sets[mWindowStart:]))]

                plt.plot(x, y, "o-", label=f"{parameter_dict[parameter]} on level {level+1}, \n")

        plt.xlabel("minimizer lengths on all levels")
        plt.ylabel("parameter value")
        plt.title(f"4Ploglogistic regression parameter:{parameter_dict[parameter]} values")
        plt.legend()
        plt.grid(True, alpha=0.3)
        #plt.ylim(0, 100)
        plt.tight_layout()

    #save the plot
    plt.savefig("./loglogistic4PRegression/parameterValues", dpi=300)
    #plt.show()
    plt.close()