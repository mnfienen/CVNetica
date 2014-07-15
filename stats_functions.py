import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp1d
from numpy.linalg import lstsq as nplstsq

def makeInputPdf(pdfRanges, pdfParam,pdfType='norm',cont_discrete='continuous'):
    '''
    pdfRanges --> gives the range for each bin in the node
    pdfParam --> a Nx2 vector of mean and std
    pdfType --> indicates the distribution assumption (must be 'norm' for now)
    cont_discrete --> indicates of node is 'continuous' or 'discrete'
                                (MUST BE CONTINUOUS FOR NOW!)
    
    returns PDF
    '''
    [N,m] = pdfParam.shape
    r = len(pdfRanges)
    if cont_discrete == 'continuous':
        PDF = np.zeros((N,r-1))
        for i in np.arange(N):
            cdf = norm.cdf(pdfRanges,pdfParam[i,0],pdfParam[i,1])
            cdf[-1]=1.0
            PDF[i,:] = np.diff(cdf)
    return PDF

def getPy(p,pdf,ranges):
    '''
    return the value of y with probability p from the PDF supplied
    this is the percentile (p)
    '''
    M = len(ranges)
    N,Nbins = pdf.shape
    f = np.cumsum(pdf, 1)
    f = np.hstack((np.zeros((N, 1)), f))
    # now make special cases for discrete or continuous
    if Nbins == M:
        # D I S C R E T E ((( NOT TESTED!!!!!)))
        if p<0.5: # indicates lower bound
            f = sum(f<p,1)
        else: # upper bound
            f = f>p
            f = M-sum(f,1)+1
            # catch any overrun
            f[f>M] = M
        py = ranges[f-1]
    else:
        # C O N T I N U O U S
        py = np.nan * np.ones((N,1))
        for i in np.arange(N):
            funique,uniqueid = np.unique(f[i,:],return_index=True)
            # note strange syntax in the next line because interp1d returns
            # a function that then performs the interpolation
            py[i] = interp1d(funique,ranges[uniqueid])(p)
    return py

def getMeanStdMostProb(pdf,ranges,continuous,blank):
    '''
    return the Mean, StdDev, and MostProbably value from a pdf
    '''
    [Nlocs,Npdf] = pdf.shape
    id1 = np.arange(Npdf)
    if continuous:
        id2  = np.arange(Npdf)+1
    else:
        id2 = np.arange(Npdf)
    # MOST PROBABLE FIRST
    # preallocate output
    mostProb = np.nan * np.ones((Nlocs,1))
    for i in np.arange(Nlocs):
        cid = np.where(pdf[i,:]==np.max(pdf[i,:]))
        if len(cid)>1:
            cid = cid[-1]
        else:
            cid = cid[0]
        mostProb[i] = 0.5*(ranges[id1[cid]]+ranges[id2[cid]])[0]
    # NEXT, DO THE MEAN and STDDEV
    
    retMean = np.atleast_2d(np.dot(pdf,(ranges[id1]+ranges[id2])/2.0)).T*blank
    retStd = np.zeros((Nlocs,1))
    # set up the sub-bin variance
    binWidthVariance = 0.5*np.atleast_2d(np.dot(pdf,(ranges[id2]-ranges[id1])/np.sqrt(3.0))).T
    for i in np.arange(Nlocs):
        retStd[i] = np.dot(pdf[i,:],(((ranges[id1]+ranges[id2])/2.0)-retMean[i])**2.0)
        retStd[i] = np.sqrt(binWidthVariance[i] + retStd[i])
        
    return retMean, retStd, mostProb,
    
def LSQR_skill(x, z, w=None):
    '''
    function to calculate skill of a model. Performs (optionally) weighted
    regression to get the coefficients.
    input:
    x --> vector of expected values (ML or mean)
    z --> observations
    w --> optional vector of std_dev weights to go on the diagonal
      N.B. ==> these weights will be inverted as they are assumed to be in the
                form of standard deviation
    output:
    sk  --> the model skill
    '''
    # first, adjust for NaNs
    keepinds = np.where(np.isnan(z) == False)[0]
    x = x[keepinds]
    z = z[keepinds]
    n = len(x)
    X = np.hstack((np.ones((n, 1)), x))
    
    # handle the weights - assuming uncorrelated errors for now
    if w:
        w = w[keepinds]
        Q = np.diag(1.0/w)
    else:
        pass
        # form weighted covariance if using weights
        # Q = np.diag(np.ones(n))

    # form the normal equations including weights
    # reserved --> XtX = np.dot(np.dot(X.T, Q), X)
    
    # solve for b
    b = nplstsq(X, z)[0]
    
    # calculate the variance of the data N.B. the mean is assumed already subtracted off from z
    obsresid = z
    msz = np.dot(obsresid.T, obsresid)/n
    # calculate the variance of the residuals
    modresid = np.dot(X, b) - z
    msr = np.dot(modresid.T, modresid)/n
    # calculate the skill
    sk = 1-(msr/msz)
    #print 'sk  --> %f' %(sk)
    
    '''
    #special case if mean already off:
    msz2 = np.dot(z.T,z)/n
    msr2 = msz2 - np.dot(np.dot(b.T,XtX),b)
    sk2 = 1-(msr2/msz2)
    print 'sk2 --> %f' %(sk2)
    '''
    return sk