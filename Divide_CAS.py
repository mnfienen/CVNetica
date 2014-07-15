import CV_tools as CVT
import pythonNeticaUtils as pyn
import numpy as np
import sys

# instantiate a pynetica object
cdat = pyn.pynetica()

# read in a CASfile name from the command line
basecas = sys.argv[1]

#read the CAS file
cdat.read_cas_file(basecas)

cdat.numfolds = int(sys.argv[2])

cdat.allfolds = CVT.all_folds()
cdat.allfolds.k_fold_maker(cdat.N, cdat.numfolds)

i = 0
header = open(basecas, 'r').readline()
useinds = np.empty(0)
for curruseinds in cdat.allfolds.leftout:
    i += 1
    useinds = (np.hstack((useinds, curruseinds)).astype(int))
    ofp = open(basecas[:-4] + '_set{0:d}_.cas'.format(i), 'w')
    ofp.write(header)
    np.savetxt(ofp, cdat.casdata[useinds])
    ofp.close()