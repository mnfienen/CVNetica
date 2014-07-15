import numpy as np
import matplotlib.pyplot as plt
import pickle, gzip, sys
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
#--modify default matplotlib settings
mpl.rcParams['font.sans-serif']          = 'Univers 57 Condensed'
mpl.rcParams['font.serif']               = 'Times'
mpl.rcParams['font.cursive']             = 'Zapf Chancery'
mpl.rcParams['font.fantasy']             = 'Comic Sans MS'
mpl.rcParams['font.monospace']           = 'Courier New'
mpl.rcParams['mathtext.default']         = 'regular'
mpl.rcParams['pdf.compression']          = 0
mpl.rcParams['pdf.fonttype']             = 42
#--figure text sizes
mpl.rcParams['legend.fontsize']  = 18
mpl.rcParams['axes.labelsize']   = 18
mpl.rcParams['xtick.labelsize']  = 18
mpl.rcParams['ytick.labelsize']  = 18
xmlfilebase = sys.argv[1]

# open up the pickle file
ifp = gzip.open(xmlfilebase +'_cdat.pklz','rb')
cdat = pickle.load(ifp)
ifp.close()

# loop over the response variables and plot experience vs. stddev
for cn in cdat.probpars.scenario.response:
    outfig = plt.figure()
    outfig.add_subplot(111)
    
    errs = cdat.basepred[cn].stats.std
    expers = cdat.BaseNeticaTests.experience[cn].case_experience
    
    plt.plot(errs,expers,'.')
    plt.xlabel('Stddev of prediction')
    plt.ylabel('Experience')
    plt.title('Experience for output node -> %s' %(cn))
    
    plt.savefig('%s_experience.pdf' %(cn))