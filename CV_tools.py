import xml.etree.ElementTree as ET
import numpy as np

###################
# Tools for reading XML input file
###################

class netica_scenario:
    def __init__(self):
        self.name = None
        self.nodesIn = None
        self.response = None
def tf2flag(intxt):
    # converts text written in XML file to True or False flag
    if intxt.lower() == 'true':
        return True
    else:
        return False

class input_parameters:
    # class and methods to read and parse XML input file
    def __init__(self, infile):
        self.infile = infile
        try:
            inpardat = ET.parse(self.infile)
        except:
            raise(FileOpenFail(self.infile)) 
        inpars = inpardat.getroot()
        self.baseNET = inpars.findall('.//control_data/baseNET')[0].text
        self.baseCAS = inpars.findall('.//control_data/baseCAS')[0].text
        self.pwdfile = inpars.findall('.//control_data/pwdfile')[0].text
        # see if rebinning is required and, if so, read in the relevant information
        try:
            self.rebin_flag = tf2flag(inpars.findall('.//control_data/rebin_flag')[0].text)
            self.originalNET = inpars.findall('.//control_data/originalNET')[0].text
        except IndexError:
            self.rebin_flag = False  # this is the default
        if self.rebin_flag:
            self.binsetup = dict()
            allbins = inpardat.findall('.//newbins/node')
            for cbin in allbins:
                tmp = cbin.attrib
                self.binsetup[cbin.text] = int(tmp['numbins'])

        self.CVflag = tf2flag(inpars.findall('.//kfold_data/CVflag')[0].text)
        self.numfolds = int(inpars.findall('.//kfold_data/numfolds')[0].text)
        self.scenario = netica_scenario()
        self.scenario.name = inpars.findall('.//scenario/name')[0].text
        self.scenario.nodesIn = []
        for cv in inpars.findall('.//scenario/input'):
            self.scenario.nodesIn.append(cv.text)
        self.scenario.response = []
        for cr in inpars.findall('.//scenario/response'):
            self.scenario.response.append(cr.text)        
        self.CASheader = list(self.scenario.nodesIn)
        self.CASheader.extend(self.scenario.response)    
        self.EMflag = tf2flag(inpars.findall('.//learnCPTdata/useEM')[0].text)
        self.report_sens = False
        self.report_sens = tf2flag(inpars.findall('.//sensitivity/report_sens')[0].text)
        self.voodooPar = float(inpars.findall('.//learnCPTdata/voodooPar')[0].text)



###################
# Tools for k-fold setup
###################
class all_folds:
    # a class containing leftout and retained indices for cross validation

    def __init__(self):
        self.leftout = list()
        self.retained = list()
        self.casfiles = list() #filenames for calibration data (retained indices only)
        self.caldata = list()  # calibration data (same as written to the calibration case file)
        self.valdata = list()  # validation data (the data left out -- will be used to calc predictions)
        # calibration and validation output from making predictions
        self.calpred = list()  
        self.valpred = list()
        self.calNODES = list()  
        self.valNODES = list()
        self.valN = list()
        self.calN = list()
        self.numfolds = None

    def k_fold_maker(self,n,k):
        # k_fold index maker
        # a m!ke@usgs joint
        # mnfienen@usgs.gov
        # k_fold_maker(n,k,allfolds)
        # input:
        #   n is the length of the sequence of indices
        #   k is the number of folds to split it into
        #   allfolds is an all_folds class
        # returns an all_folds with each member having k elements
        # allfolds.leftout[i] is the left out indices of fold i
        # allfolds.retained[i] is the kept indices of fold i
        currinds = np.arange(n)
        inds_per_fold = np.int(np.floor(n/k))
        dingleberry = np.remainder(n, k)
        for i in np.arange(k-1):
            allinds = currinds.copy()    
            np.random.shuffle(currinds)
            self.leftout.append(currinds[0:inds_per_fold].copy())
            self.retained.append(np.setdiff1d(np.arange(n), self.leftout[i]))
            currinds = currinds[inds_per_fold:]
        self.leftout.append(currinds)
        self.retained.append(np.setdiff1d(np.arange(n), self.leftout[-1]))
        self.numfolds = k

#################
# Error classes
#################

# -- cannot open an input file
class FileOpenFail(Exception):
    def __init__(self,filename):
        self.fn = filename
    def __str__(self):
        return('\n\nCould not open %s.' %(self.fn))    