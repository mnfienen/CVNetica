import numpy as np
import scipy.stats.mstats as scpstat

class netica_binning:
    def __init__(self, x, n):
        x = x[np.isnan(x) == False]  # only retain values in x that are not NaN
        self.x = x  #vector of values to bin
        self.n = n  # number of bins

    def bin_thresholds(self):
        self.probs = np.linspace(0, 1.0, self.n+1)
        self.binlevels = scpstat.mquantiles(self.x, self.probs)

