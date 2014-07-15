import pickle, gzip
ifp = gzip.open('example_cdat.pklz','rb')
cdat = pickle.load(ifp)
ifp.close()
