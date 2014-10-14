import gzip
import pickle
import sys
import matplotlib.pyplot as plt
import numpy as np

cpickle = sys.argv[1]
ifp = gzip.open(cpickle, 'rb')
cdat = pickle.load(ifp)
ifp.close()


def one2one(meas,mod,title):
    plt.figure()
    plt.plot(meas,mod,'o')
    minmeas = np.min(meas)
    maxmeas = np.max(meas)
    minmod = np.min(mod)
    maxmod = np.max(mod)

    plt.plot([np.min((minmeas,minmod)),np.max((maxmeas,maxmod))],[np.min((minmeas,minmod)),np.max((maxmeas,maxmod))], 'r-')
    plt.xlabel('Measured')
    plt.ylabel('Forecast')
    plt.title(title)
    plt.savefig(title + '.pdf')

one2one(cdat.casdata['SW_SRC'],cdat.basepred['SW_SRC'].stats.median, 'SW_SRC_cal')
one2one(cdat.val_casdata['SW_SRC'],cdat.VALpred['SW_SRC'].stats.median, 'SW_SRC_val')

i=1
