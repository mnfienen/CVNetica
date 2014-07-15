import os

numfolds = 8
baseroot = 'SEEDA1A1B7D7E7'
basecases = ['sqrt_dist']
allfiles = os.listdir(os.getcwd())

for ccas in basecases:
    runline1 = 'python Divide_CAS.py {0:s}{1:s}.cas {2:d}'.format(baseroot, ccas, numfolds)
    print runline1
    os.system(runline1)
    runline2 = 'python XML_for_multiple_sets.py glass{0:s}.xml {1:d}'.format(ccas, numfolds)
    print runline2
    os.system(runline2)


allfiles = os.listdir(os.getcwd())
for cfile in allfiles:
    if cfile[-4:] == '.xml':
        runstring = u"python CV_driver.py {0:s}".format(cfile)
        print runstring
        os.system(runstring)

allfiles = os.listdir(os.getcwd())
for cfile in allfiles:
    if cfile[-4:] == '.dat' and ('CAL' in cfile or 'VAL' in cfile):
        runstring = "python postproc_summarize_kfold.py {0:s}".format(cfile)
        print runstring
        os.system(runstring)
