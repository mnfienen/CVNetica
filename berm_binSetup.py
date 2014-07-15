import neticaBinTools as nBT
import numpy as np
import pythonNetica as pyn

import CV_tools as CVT
import numpy as np
import pickle, gzip
import sys
import os
import shutil
print os.getcwd()

# set the base CAS and NETA file names
indat = np.genfromtxt('berm.cas', dtype=None, names=True)
basenet = 'berm.neta'
parfileroot = 'berm_cvDriverParam_'

cdat = pyn.pynetica()

# set the bin arrangements --> note that there must be already
# created XML files for these, including the bin XML files
# loop through parameter control xmls files like this...
sets = ['k2']

# now loop over each
for cset in sets:
    # read in the parameters from the XML file
    cdat.probpars = CVT.input_parameters('{0:s}{1:s}.xml'.format(parfileroot,cset))

    # copy over the base neta file to the one named in the XML file for this bin set
    # shutil.copyfile(basenet, cdat.probpars.baseNET)

    # Initialize a pynetica instance/env using password in a text file
    cdat.pyt.start_environment(cdat.probpars.pwdfile)

    # read in the data from a base cas file
    cdat.read_cas_file(cdat.probpars.baseCAS)
    
    # sets equiprobable bins for each node
    cdat.UpdateNeticaBinThresholds(basenet)

