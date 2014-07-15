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
        
class input_parameters:
    # class and methods to read and parse XML input file
    def __init__(self,infile):
        self.infile = infile
        try:
            inpardat = ET.parse(self.infile)
        except:
            raise(FileOpenFail(self.infile)) 
        inpars = inpardat.getroot()
        self.baseNET = inpars.findall('.//input_files/baseNET')[0].text
        self.baseCAS = inpars.findall('.//input_files/baseCAS')[0].text
        self.pwdfile = inpars.findall('.//input_files/pwdfile')[0].text
        self.numfolds = int(inpars.findall('.//kfold_data/numfolds')[0].text)
        self.scenario = netica_scenario()
        self.scenario.name =  inpars.findall('.//scenario/name')[0].text
        self.scenario.nodesIn = []
        for cv in inpars.findall('.//scenario/input'):
            self.scenario.nodesIn.append(cv.text)
 		self.scenario.response = []
        for cv in inpars.findall('.//scenario/response'):
            self.scenario.response.append(cv.text)
 