__author__ = 'nplant'
import gzip
import pickle
import sys
from scipy import io

cpickle = sys.argv[1]
ifp = gzip.open(cpickle, 'rb')
cdat = pickle.load(ifp)
ifp.close()

'''
this is for figuring out how to save to matlab
io.savemat('dump.mat', dict(cdat=cdat))
'''


i=1
