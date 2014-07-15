import pythonNeticaUtils as pyn
import CV_tools as CVT
import numpy as np
import pickle, gzip
import sys


############
# CONFIGURATION FILE NAME
try:
    parfile = sys.argv[1]
except:
    parfile = 'example.xml'
############
# initialize
cdat = pyn.pynetica()

# read in the problem parameters
cdat.probpars = CVT.input_parameters(parfile)

# read in cas file
cdat.read_cas_file(cdat.probpars.baseCAS)

# Initialize a pynetica instance/env using password in a text file
cdat.start_environment(cdat.probpars.pwdfile)



junks = cdat.NodeParentIndexing(cdat.probpars.baseNET,cdat.probpars.baseCAS)



'''
# let's operate on a single example here for now
testnode = cdat.GetNodeNamed('islandwidth',cnet)

# get a list of the parents of the node
parents = cdat.GetNodeParents(testnode)

# for each parent, get its set of possible states

# now find how many possible node states there are
cartprod = cdat.SizeCartesianProduct(parents)

print cartprod
'''