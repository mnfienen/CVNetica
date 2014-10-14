import pythonNetica as pyn

import CV_tools as CVT
import numpy as np
import pickle
import gzip
import sys
import shutil

'''
CV_driver.py

a cross-validation driver for Netica
a m!ke@usgs joint
'''


############
# CONFIGURATION FILE NAME

parfile = sys.argv[1]

############
# initialize
cdat = pyn.pynetica()

# read in the problem parameters
cdat.probpars = CVT.input_parameters(parfile)

# Initialize a pynetica instance/env using password in a text file
cdat.pyt.start_environment(cdat.probpars.pwdfile)
cdat.pyt.LimitMemoryUsage(5.0e16)  # --> crank up the memory available

# read in the data from a base cas file
cdat.read_cas_file(cdat.probpars.baseCAS)

# read in the data from a prediction cas file
cdat.read_pred_cas_file(cdat.probpars.valCAS)


# perform rebinning if requested
if cdat.probpars.rebin_flag:

    # note that in CV case we first copy originalNET over to baseNET --- here we start with baseNET
    # sets equiprobable bins for each node as requested
    cdat.UpdateNeticaBinThresholds(True)

# set up the experience node indexing
cdat.NodeParentIndexing(cdat.probpars.baseNET, cdat.probpars.baseCAS)

# associate the casefile with the net
print '*'*5 + 'Learning base casefile in base net' + '*'*5 + '\n\n'

cdat.pyt.rebuild_net(cdat.probpars.baseNET,
                         cdat.probpars.baseCAS,
                         cdat.probpars.voodooPar,
                         cdat.probpars.baseNET,
                         cdat.probpars.EMflag)
                         
# run the predictions using the base net -->
cdat.basepred, cdat.NETNODES = cdat.predictBayes(cdat.probpars.baseNET, cdat.N, cdat.casdata)

# run the predictions using the base net -->
cdat.VALpred, cdat.VALNETNODES = cdat.predictBayes(cdat.probpars.baseNET, cdat.N_val, cdat.val_casdata)


# Done with Netica so shut it down
cdat.pyt.CloseNetica()


# first need to sanitize away any ctypes/Netica pointers
cdat.pyt.sanitize()
# now dump into a pickle file
outfilename = parfile[:-4] + '_cdat.pklz'
print 'Dumping cdat to pickle file --> {0:s}'.format(outfilename)
ofp = gzip.open(outfilename, 'wb')
pickle.dump(cdat, ofp)
ofp.close()

