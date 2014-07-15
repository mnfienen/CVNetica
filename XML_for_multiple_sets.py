import re
import sys
import xml.etree.ElementTree as ET
import numpy as np

infile = sys.argv[1]  # xml file
numfolds = int(sys.argv[2])  # numfolds

inpardat = ET.parse(infile)
inpars = inpardat.getroot()
baseCAS = inpars.findall('.//input_files/baseCAS')[0].text
scenarioname = inpars.findall('.//scenario/name')[0].text

indat = open(infile, 'r').readlines()

for i in np.arange(numfolds):
    ofp = open(infile[:-4] + "_set{0:d}_.xml".format(i + 1), 'w')
    for line in indat:
        if '<baseCAS>' in line:
            ofp.write(re.sub(baseCAS, baseCAS[:-4] + "_set{0:d}_.cas".format(i + 1), line.strip()) + '\n')
        elif '<name>' in line:
            ofp.write(re.sub(scenarioname, scenarioname + "_set{0:d}".format(i + 1), line.strip()) + '\n')
        else:
            ofp.write(line.strip() + '\n')
    ofp.close()