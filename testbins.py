import neticaBinTools as nBT
import numpy as np
import pythonNetica as pyn

import CV_tools as CVT
import numpy as np
import pickle, gzip
import sys


indat = np.genfromtxt('glacial4433.cas', dtype=None, names=True)


cdat = pyn.pynetica()

cdat.probpars = CVT.input_parameters('example.xml')

# Initialize a pynetica instance/env using password in a text file
cdat.pyt.start_environment(cdat.probpars.pwdfile)

# read in the data from a base cas file
cdat.read_cas_file(cdat.probpars.baseCAS)

cdat.UpdateNeticaBinThresholds()

bins = nBT.netica_binning(indat['log_TRANGM'], 4)

bins.bin_thresholds()

