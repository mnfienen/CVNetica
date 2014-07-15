import numpy as np
import os
import re
import sys
import ctypes as ct
import pythonNeticaConstants as pnC
import pythonNeticaTools as pyT
import cthelper as cth
import neticaBinTools as nBT
import stats_functions as statfuns

class parent_inds:
    def __init__(self):
        self.parent_names = list()
        self.parent_indices = list()


class netica_test:
    def __init__(self):
        self.logloss = None
        self.errrate = None
        self.quadloss = None
        self.confusion_matrix = None
        self.experience = None

class pred_stats:
    def __init__(self):
        self.alpha = None
        self.palpha = None
        self.mean = None
        self.mostProb = None
        self.std = None
        self.median = None
        self.p025 = None
        self.p05 = None
        self.p25 = None
        self.p75 = None
        self.p95 = None
        self.p975 = None
        self.palphaPlus = None
        self.pAlphaMinus = None
        self.meanabserrM = None
        self.meanabserrML = None
        self.rmseM = None
        self.meaneM = None
        self.rmseML = None
        self.meaneML = None
        self.parent_states = None


class predictions:
    def __init__(self):
        self.z = None
        self.pdf = None
        self.pdfIn = None
        self.ranges = None
        self.rangesplt = None
        self.priorPDF = None
        self.probModelPrior = None
        self.probModelUpdate = None
        self.dataPDF = None
        self.ofp = None
        # statistics go here
        self.stats = None

class pynetica:
    def __init__(self):
        self.casdata = None
        self.pyt = pyT.pyneticaTools()
        self.basepred = None
        self.parent_inds = None
        self.NeticaTests = dict()
        self.NeticaTests['CAL'] = list()
        self.NeticaTests['VAL'] = list()
        self.probpars = None


    #############################################
    # Major validation and prediction functions #
    #############################################


    def NodeParentIndexing(self, netName, casfile):
        '''
        Find all the configurations of states in the parent nodes for each response node
        This is used only for
        '''
        # open the net stored in netName
        cnet = self.pyt.OpenNeticaNet(netName)
        #get the nodes and their number
        allnodes = self.pyt.GetNetNodes(cnet)
        numnodes = self.pyt.LengthNodeList(allnodes)

        #parent indices dictionary for the results
        parent_indices = dict()
        # now focus in on the response nodes only
        respnodes = self.probpars.scenario.response
        for cr in respnodes:
            parent_indices[cr] = parent_inds()
            crespnode = self.pyt.GetNodeNamed(cr, cnet)
            # get the parent nodes and their names
            cparents = self.pyt.GetNodeParents(crespnode)
            numparents = self.pyt.LengthNodeList(cparents)
            for cp in np.arange(numparents):
                tmpnode = self.pyt.NthNode(cparents, cp)
                parent_indices[cr].parent_names.append(
                    cth.c_char_p2str(self.pyt.GetNodeName(tmpnode)))

        # open a streamer to the CAS file we will read over
        cas_streamer = self.pyt.NewFileStreamer(casfile)

        # loop over the cases
        for ccas in np.arange(self.N):
            if ccas == 0:
                case_posn = pnC.netica_const.FIRST_CASE
            else:
                case_posn = pnC.netica_const.NEXT_CASE
            # first set the findings according to what's in the case file
            case_posn_out = self.pyt.ReadNetFindings2(case_posn, cas_streamer, allnodes)
            # now, for each parent, in order, read the states
            for cr in respnodes:
                tmpinds = list()
                for cp in parent_indices[cr].parent_names:
                    cnode = self.pyt.GetNodeNamed(cp,cnet)
                    tmpinds.append(self.pyt.GetNodeFinding(cnode))
                parent_indices[cr].parent_indices.append(tmpinds)
        for cr in respnodes:
            print 'making into an array --> %s' %(cr)
            parent_indices[cr].parent_indices = np.array(
                    parent_indices[cr].parent_indices, dtype=int)
        # clean up the temporary streamer and net
        self.pyt.DeleteNet(cnet)
        self.pyt.DeleteStream(cas_streamer)
        self.parent_inds = parent_indices

    def UpdateNeticaBinThresholds(self):
        '''
        Function that reads in new numbers of bins and sets each node to
        have equiprobable bins in that number
        '''
        # first open the net
        print "*"*5 + "Setting up a rebinned net {0:s} copying nodes from {1:s}".format(self.probpars.baseNET,
                                                                                self.probpars.originalNET) + "*"*5
        cnet = self.pyt.OpenNeticaNet(self.probpars.originalNET)
        for cbin in self.probpars.binsetup:
            cnodebins = nBT.netica_binning(self.casdata[cbin], self.probpars.binsetup[cbin])
            cnodebins.bin_thresholds()
            cnode = self.pyt.GetNodeNamed(cbin, cnet)
            if self.probpars.binsetup[cbin] != 0:
                # only do this if requested bins not zero. Else, use same bins as input net
                print "Setting node {0:s} to have {1:d} bins".format(cbin, self.probpars.binsetup[cbin])
                self.pyt.SetNodeLevels(cnode, cnodebins.binlevels)
            else:
                print "NOT Setting node {0:s} to have {1:d} bins".format(cbin, self.probpars.binsetup[cbin])
        outfile_streamer = self.pyt.NewFileStreamer(self.probpars.baseNET)
        self.pyt.CompileNet(cnet)
        print "Writing new bin configurations for net to: {0:s}".format(self.probpars.baseNET)
        self.pyt.WriteNet(cnet, outfile_streamer)
        self.pyt.DeleteStream(outfile_streamer)
        self.pyt.DeleteNet(cnet)

    def PredictBayesNeticaCV(self, cfold, cnetname, calval):
        '''
        function using Netica built-in testing functionality to evaluate Net
        '''
        ctestresults = netica_test()

        # open up the current net
        cnet = self.pyt.OpenNeticaNet(cnetname)
        #retract all the findings
        self.pyt.RetractNetFindings(cnet)

        # first create a caseset with the current leftout indices casefile
        if cfold > -10:
            if calval.upper() == 'CAL':
                ccasefile = '{0:s}_fold_{1:d}.cas'.format(self.probpars.scenario.name, cfold)
            elif calval.upper() == 'VAL':
                ccasefile = '{0:s}_fold_{1:d}_leftout.cas'.format(self.probpars.scenario.name, cfold)
            else:
                pass
        # unless this is the base case -->
        else:
            ccasefile = self.probpars.baseCAS
        currcases = self.pyt.NewCaseset('cval{0:d}'.format(np.abs(cfold)))
        ccaseStreamer = self.pyt.NewFileStreamer(ccasefile)
        self.pyt.AddFileToCaseset(currcases, ccaseStreamer, 100.0)

        # create a set of prediction nodes
        numprednodes = len(self.probpars.scenario.response)
        cnodelist = self.pyt.NewNodeList2(numprednodes, cnet)
        for i, cn in enumerate(self.probpars.scenario.response):
            cnode = self.pyt.GetNodeNamed(cn, cnet)
            self.pyt.SetNthNode(cnodelist, i, cnode)
        # create a tester object
        ctester = self.pyt.NewNetTester(cnodelist, cnodelist)
        self.pyt.DeleteNodeList(cnodelist)

        # test the network using the left-out cases
        # first retract all the findings and compile the net
        self.pyt.TestWithCaseset(ctester, currcases)
        self.pyt.DeleteCaseset(currcases)
        #
        # now get the results
        #
        ctestresults.logloss = dict()
        ctestresults.errrate = dict()
        ctestresults.quadloss = dict()
        ctestresults.confusion_matrix = dict()
        ctestresults.experience = dict()
        for cn in self.probpars.scenario.response:
            cnode = self.pyt.GetNodeNamed(cn, cnet)
            # get log loss
            ctestresults.logloss[cn] = self.pyt.GetTestLogLoss(ctester, cnode)
            print "LogLoss for {0:s} --> {1:f}".format(cn, ctestresults.logloss[cn])
            # get error rate
            ctestresults.errrate[cn] = self.pyt.GetTestErrorRate(ctester, cnode)
            print "ErrorRate for {0:s} --> {1:f}".format(cn, ctestresults.errrate[cn])
            # get quadratic loss
            ctestresults.quadloss[cn] = self.pyt.GetTestQuadraticLoss(ctester, cnode)
            print "QuadLoss for {0:s} --> {1:f}".format(cn, ctestresults.quadloss[cn])
            # write confusion matrix --- only for the base case
            if cfold < 0:
                print "Calculating confusion matrix for {0:s}".format(cn)
                ctestresults.confusion_matrix[cn] = self.pyt.ConfusionMatrix(ctester, cnode)
                # also calculate the experience for the node
                print "Calculating Experience for the base Net, node --> {0:s}".format(cn)
                ctestresults.experience[cn] = self.pyt.ExperienceAnalysis(cn, cnet)

        self.pyt.DeleteNetTester(ctester)
        self.pyt.DeleteNet(cnet)

        # write to the proper dictionary
        if cfold > -10:
            if calval.upper() == 'CAL':
                self.NeticaTests['CAL'].append(ctestresults)
            elif calval.upper() == 'VAL':
                self.NeticaTests['VAL'].append(ctestresults)
            else:
                pass
        else:
            self.BaseNeticaTests = ctestresults


    def ExperiencePostProc(self):
        print "Post-Processing Experience data, matching with predicted nodes and cases"
        for cn in self.probpars.scenario.response:
            for ccas in np.arange(self.N):
                testinds = self.parent_inds[cn].parent_indices[ccas, :]
                tmp = testinds-self.BaseNeticaTests.experience[cn].parent_states
                tmp = np.sum(np.abs(tmp), axis=1)
                cind = np.where(tmp == 0)
                self.BaseNeticaTests.experience[cn].case_experience.append(
                    self.BaseNeticaTests.experience[cn].node_experience[cind[0]])



    def predictBayes(self, netName, N, casdata):
        '''
        netName --> name of the built net to make predictions on
        '''
        # first read in the information about a Net's nodes
        cNETNODES = self.pyt.ReadNodeInfo(netName)
        '''
        Initialize output 
        '''
        # initialize dictionary of predictions objects
        cpred = dict()

        print "Making predictions for net named --> {0:s}".format(netName)
        cnet = self.pyt.OpenNeticaNet(netName)
        #retract all the findings
        self.pyt.RetractNetFindings(cnet)
        for CN in cNETNODES:
            CNODES = cNETNODES[CN]
            Cname = CNODES.name
            if Cname in self.probpars.scenario.response:
                cpred[Cname] = predictions()
                cpred[Cname].stats = pred_stats()
                Nbins = CNODES.Nbeliefs
                cpred[Cname].pdf = np.zeros((N, Nbins))
                cpred[Cname].ranges = np.array(CNODES.levels)
                # get plottable ranges
                if Nbins < len(CNODES.levels):
                    # continuous, so plot bin centers
                    CNODES.continuous = True
                    cpred[Cname].continuous = True
                    cpred[Cname].rangesplt = (cpred[Cname].ranges[1:] -
                                              0.5*np.diff(cpred[Cname].ranges))
                else:
                    #discrete so just use the bin values
                    cpred[Cname].rangesplt = cpred[Cname].ranges.copy()

                cpred[Cname].priorPDF = CNODES.beliefs

        allnodes = self.pyt.GetNetNodes(cnet)
        numnodes = self.pyt.LengthNodeList(allnodes)
        #
        # Now loop over each input and get the Netica predictions
        #
        for i in np.arange(N):
            sys.stdout.write('predicting value {0} of {1}\r'.format(i,N))
            sys.stdout.flush()
            # first have to enter the values for each node
            # retract all the findings again
            self.pyt.RetractNetFindings(cnet)
            for cn in np.arange(numnodes):
                cnode = self.pyt.NthNode(allnodes, cn)
                cnodename = cth.c_char_p2str(self.pyt.GetNodeName(cnode))
                # set the current node values
                if cnodename in self.probpars.scenario.nodesIn:
                    self.pyt.EnterNodeValue(cnode, casdata[cnodename][i])
            for cn in np.arange(numnodes):
            # obtain the updated beliefs from ALL nodes including input and output
                cnode = self.pyt.NthNode(allnodes, cn)
                cnodename = cth.c_char_p2str(self.pyt.GetNodeName(cnode))
                if cnodename in self.probpars.scenario.response:
                    # get the current belief values
                    cpred[cnodename].pdf[i, :] = cth.c_float_p2float(
                        self.pyt.GetNodeBeliefs(cnode),
                        self.pyt.GetNodeNumberStates(cnode))
        #
        # Do some postprocessing for just the output nodes
        #
        currstds = np.ones((N, 1))*1.0e-16
        for i in self.probpars.scenario.response:
            print 'postprocessing output node --> {0:s}'.format(i)
            # record whether the node is continuous or discrete
            if cpred[i].continuous:
                curr_continuous = 'continuous'
            else:
                curr_continuous = 'discrete'
            pdfRanges = cpred[i].ranges
            cpred[i].z = np.atleast_2d(casdata[i]).T
            pdfParam = np.hstack((cpred[i].z, currstds))
            pdfData = statfuns.makeInputPdf(pdfRanges, pdfParam, 'norm', curr_continuous)

            cpred[i].probModelUpdate = np.nansum(pdfData * cpred[i].pdf, 1)
            cpred[i].probModelPrior = np.nansum(pdfData * np.tile(cpred[i].priorPDF,
                                                (N, 1)), 1)
            cpred[i].logLikelihoodRatio = (np.log10(cpred[i].probModelUpdate + np.spacing(1)) -
                                           np.log10(cpred[i].probModelPrior + np.spacing(1)))
            cpred[i].dataPDF = pdfData.copy()
            # note --> np.spacing(1) is like eps in MATLAB
            # get the PDF stats here

            print 'getting stats'
            cpred = self.PDF2Stats(i, cpred, alpha=0.1)

        self.pyt.DeleteNet(cnet)
        return cpred, cNETNODES

    def PDF2Stats(self, nodename, cpred, alpha=None):
        '''
        extract statistics from the PDF informed by a Bayesian Net

        most information is contained in self which is a pynetica object
        however, the nodename indicates which node to calculate stats for
        '''

        # normalize the PDF in case it doesn't sum to unity        
        pdf = np.atleast_2d(cpred[nodename].pdf).copy()
        pdf /= np.tile(np.atleast_2d(np.sum(pdf, 1)).T, (1, pdf.shape[1]))

        # Start computing the statistics
        [Nlocs, Npdf] = pdf.shape
        blank = 0.0 + ~np.isnan(pdf[:, 0])
        blank[blank == 0] = np.nan
        blank = np.atleast_2d(blank).T

        # handle the specific case of a user-specified percentile range
        if alpha:
            self.alpha = alpha
            dalpha = (1.0 - alpha)/2.0
            # first return the percentile requested in bAlpha
            cpred[nodename].stats.palpha = blank*statfuns.getPy(
                alpha, pdf,
                cpred[nodename].ranges)

            # now get the tails from requested bAlpha
            # upper tail
            cpred[nodename].stats.palphaPlus = blank*statfuns.getPy(
                1.0-dalpha, pdf,
                cpred[nodename].ranges)

            # lower tail
            cpred[nodename].stats.palphaMinus = blank*statfuns.getPy(
                dalpha, pdf,
                cpred[nodename].ranges)
        # now handle the p75,p95, and p975 cases
        # 75th percentiles
        cpred[nodename].stats.p25 = blank*statfuns.getPy(0.25, pdf,
                                                         cpred[nodename].ranges)
        cpred[nodename].stats.p75 = blank*statfuns.getPy(0.75, pdf,
                                                         cpred[nodename].ranges)
        # 95th percentiles
        cpred[nodename].stats.p05 = blank*statfuns.getPy(0.05, pdf,
                                                         cpred[nodename].ranges)
        cpred[nodename].stats.p95 = blank*statfuns.getPy(0.95, pdf,
                                                         cpred[nodename].ranges)
        # 97.5th percentiles
        cpred[nodename].stats.p025 = blank*statfuns.getPy(0.025, pdf,
                                                          cpred[nodename].ranges)
        cpred[nodename].stats.p975 = blank*statfuns.getPy(0.975, pdf,
                                                          cpred[nodename].ranges)
        # MEDIAN  
        cpred[nodename].stats.median = blank*statfuns.getPy(0.5, pdf,
                                                            cpred[nodename].ranges)

        # now get the mean, ML, and std values
        (cpred[nodename].stats.mean,
         cpred[nodename].stats.std,
         cpred[nodename].stats.mostProb) = statfuns.getMeanStdMostProb(pdf,
                                                                       cpred[nodename].ranges,
                                                                       cpred[nodename].continuous,
                                                                       blank)

        cpred[nodename].stats.skMean = statfuns.LSQR_skill(
            cpred[nodename].stats.mean,
            cpred[nodename].z-np.nanmean(cpred[nodename].z))

        cpred[nodename].stats.skML = statfuns.LSQR_skill(
            cpred[nodename].stats.mostProb,
            cpred[nodename].z-np.nanmean(cpred[nodename].z))
        Mresid = (cpred[nodename].stats.mean -
                  cpred[nodename].z)
        cpred[nodename].stats.rmseM = (
            np.sqrt(np.nanmean(Mresid**2)))
        cpred[nodename].stats.meaneM = np.nanmean(Mresid)
        MLresid = (cpred[nodename].stats.mostProb -
                   cpred[nodename].z)
        cpred[nodename].stats.rmseML = (
            np.sqrt(np.nanmean(MLresid**2)))
        cpred[nodename].stats.meaneML = np.nanmean(MLresid)
        cpred[nodename].stats.meanabserrM = np.nanmean(np.abs(Mresid))
        cpred[nodename].stats.meanabserrML = np.nanmean(np.abs(MLresid))


        return cpred

    def PredictBayesPostProc(self, cpred, outname, casname, cNeticaTestStats):
        ofp = open(outname, 'w')
        ofp.write('Validation statistics for net --> {0:s} and casefile --> {1:s}\n'.format(outname, casname))
        ofp.write('%14s '*12
                  %('Response', 'skillMean', 'rmseMean', 'meanErrMean', 'meanAbsErrMean',
                    'skillML', 'rmseML', 'meanErrML', 'meanAbsErrML', 'LogLoss', 'ErrorRate', 'QuadraticLoss')
                  + '\n')
        for i in self.probpars.scenario.response:
            print 'writing output for --> {0:s}'.format(i)
            ofp.write('%14s %14.4f %14.6e %14.6e %14.6e %14.4f %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n'
                      %(i, cpred[i].stats.skMean,
                        cpred[i].stats.rmseM,
                        cpred[i].stats.meaneM,
                        cpred[i].stats.meanabserrM,
                        cpred[i].stats.skML,
                        cpred[i].stats.rmseML,
                        cpred[i].stats.meaneML,
                        cpred[i].stats.meanabserrML,
                        cNeticaTestStats.logloss[i],
                        cNeticaTestStats.errrate[i],
                        cNeticaTestStats.quadloss[i]))
        ofp.close()
        outfileConfusion = re.sub('base_stats', 'Confusion', outname)
        ofpC = open(outfileConfusion, 'w')
        ofpC.write('Confusion matrices for net --> %s and casefile --> %s\n'
                  %(outname, casname))
        for j in self.probpars.scenario.response:
            ofpC.write('*'*16 + '\nConfusion matrix for %s\n' %(j) + '*'*16 + '\n')
            numstates = len(self.NETNODES[j].levels)-1
            ofpC.write('%24s' %(''))
            for i in np.arange(numstates):
                ofpC.write('%24s' %('%8.4e--%8.4e'%(self.NETNODES[j].levels[i],
                                        self.NETNODES[j].levels[i+1])))
            ofpC.write('\n')
            for i in np.arange(numstates):
                ofpC.write('%24s' %('%8.4e--%8.4e'%(self.NETNODES[j].levels[i],
                                        self.NETNODES[j].levels[i+1])))
                for k in cNeticaTestStats.confusion_matrix[j][i,:]:
                    ofpC.write('%24d' %(k))
                ofpC.write('\n')
            ofpC.write('\n' * 2)
        ofpC.close()




    def PredictBayesPostProcCV(self, cpred, numfolds, ofp, calval, cNeticaTestStats):
        for cfold in np.arange(numfolds):
            for j in self.probpars.scenario.response:
                print 'writing %s cross-validation output for --> %s' %(calval, j)
                ofp.write('%14d %14s %14.4f %14.6e %14.6e %14.6e %14.4f %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n'
                      %(cfold, j, cpred[cfold][j].stats.skMean,
                        cpred[cfold][j].stats.rmseM,
                        cpred[cfold][j].stats.meaneM,
                        cpred[cfold][j].stats.meanabserrM,
                        cpred[cfold][j].stats.skML,
                        cpred[cfold][j].stats.rmseML,
                        cpred[cfold][j].stats.meaneML,
                        cpred[cfold][j].stats.meanabserrML,
                        cNeticaTestStats[cfold].logloss[j],
                        cNeticaTestStats[cfold].errrate[j],
                        cNeticaTestStats[cfold].quadloss[j]))

    def SummarizePostProcCV(self):
        '''
        Method to consolidate metrics accross all folds in a cross-validation into a single file
        This is done after the fully detailed files are already written
        '''
        for infile in [self.probpars.Cal_outfile, self.probpars.Val_outfile]:
            header = open(infile, 'r').readlines()[0:3]
            stats = ['min', 'max', 'mean', 'median', 'std']

            response_headers = ['skillMean', 'rmseMean', 'meanErrMean', 'meanAbsErrMean',
                                'skillML', 'rmseML', 'meanErrML', 'meanAbsErrML',
                                'LogLoss', 'ErrorRate', 'QuadraticLoss']

            indat = np.genfromtxt(infile, skiprows=3, dtype=None, names=True)
            unique_responses = np.unique(indat['Response'])
            numfolds = np.max(indat['Current_Fold'])+1
            outdat = dict()  # dictionary of responses

            for cres in unique_responses:
                outdat[cres] = dict()
                outdat[cres]['min'] = dict()
                outdat[cres]['max'] = dict()
                outdat[cres]['mean'] = dict()
                outdat[cres]['median'] = dict()
                outdat[cres]['std'] = dict()
                currinds = np.where(indat['Response'] == cres)[0]
                for cstat in response_headers:
                    outdat[cres]['min'][cstat] = np.min(indat[cstat][currinds])
                    outdat[cres]['max'][cstat] = np.max(indat[cstat][currinds])
                    outdat[cres]['mean'][cstat] = np.mean(indat[cstat][currinds])
                    outdat[cres]['median'][cstat] = np.median(indat[cstat][currinds])
                    outdat[cres]['std'][cstat] = np.std(indat[cstat][currinds])
            ofp = open(infile[:-4] + '_SUMMARY.dat', 'w')
            ofp.write('SUMMARY STATISTICS-->\n')
            for line in header:
                ofp.write(line)
            ofp.write('{0:>16s}{1:>16s}'.format('Stat', 'Response'))
            for chead in response_headers:
                ofp.write('{0:>16s}'.format(chead))
            ofp.write('\n')
            for currstat in stats:
                for cresp in unique_responses:
                    ofp.write('{0:>16s}'.format(currstat))
                    ofp.write('{0:>16s}'.format(cresp))
                    for cval in response_headers:
                        if 'skill' in cval:
                            ofp.write('{0:16.5f}'.format(outdat[cresp][currstat][cval]))
                        else:
                            ofp.write('{0:16.5e}'.format(outdat[cresp][currstat][cval]))
                    ofp.write('\n')
                ofp.write('\n')
            ofp.close()

    def SensitivityAnalysis(self):
        '''
        Peforms sensitivity analysis on each response node assuming all 
        input nodes are active (as defined in self.probpars.scenario)
        
        Reports results to a text file.
        '''
        print '\n' * 3 + '*' * 10 + '\n' + 'Performing Sensitity Analysis\n' + '*'*10
        # meke a streamer to the Net file
        net_streamer = self.pyt.NewFileStreamer(self.probpars.baseNET)
        # read in the net using the streamer        
        cnet = self.pyt.ReadNet(net_streamer)
        # remove the input net streamer
        self.pyt.DeleteStream(net_streamer)
        self.pyt.CompileNet(cnet)
        self.sensitivityvar = dict()
        self.sensitivityEntropy = dict()
        self.sensitivityEntropyNorm = dict()
        self.percentvarreduction = dict()
        allnodes = list()
        allnodes.extend(self.probpars.scenario.nodesIn)
        allnodes.extend(self.probpars.scenario.response)
        for cres in self.probpars.scenario.response:
            print "Calculating sensitivity to node --> %s" %(cres)
            # calculate the sensitivity for each response variable using all nodes  as Vnodes
            Qnode = self.pyt.GetNodeNamed(cres,cnet)
            Vnodes = self.pyt.GetNetNodes(cnet)
            self.sensitivityvar[cres] = dict()
            self.sensitivityEntropy[cres] = dict()
            self.sensitivityEntropyNorm[cres] = dict()
            self.percentvarreduction[cres] = dict()
            sensvar = self.pyt.NewSensvToFinding(Qnode,Vnodes,ct.c_int(pnC.netica_const.VARIANCE_OF_REAL_SENSV))
            sensmutual = self.pyt.NewSensvToFinding(Qnode,Vnodes,ct.c_int(pnC.netica_const.ENTROPY_SENSV))
            for cn in allnodes:
                Vnode = self.pyt.GetNodeNamed(cn,cnet)
                self.sensitivityvar[cres][cn] = self.pyt.GetVarianceOfReal(sensvar,Vnode)
                self.sensitivityEntropy[cres][cn] = self.pyt.GetMutualInfo(sensmutual,Vnode)
                # percent variance reduction is the variance reduction of a node divided by variance reduction of self
            for cn in allnodes:
                self.percentvarreduction[cres][cn] = self.sensitivityvar[cres][cn]/self.sensitivityvar[cres][cres]
                self.sensitivityEntropyNorm[cres][cn] = self.sensitivityEntropy[cres][cn]/self.sensitivityEntropy[cres][cres]
            print "Deleting sensitivity to --> %s" %(cres)
            self.pyt.DeleteSensvToFinding(sensvar)
            self.pyt.DeleteSensvToFinding(sensmutual)
        self.pyt.DeleteNet(cnet)

        # #### WRITE OUTPUT #### #
        ofp = open(self.probpars.scenario.name + 'Sensitivity.dat','w')
        ofp.write('Sensitivity analysis for scenario --> %s\n' %(self.probpars.scenario.name))
        ofp.write('Base Case Net: %s\nBase Case Casfile: %s\n' %(self.probpars.baseNET,self.probpars.baseCAS))
        # write out the raw variance reduction values        
        ofp.write('#'*10 + '   Raw Variance Reduction Values   ' + '#'*10 + '\n')
        ofp.write('{0:<14s}'.format('Response_node '))
        for cn in allnodes:
            ofp.write('%-14s' %(cn))
        ofp.write('\n')
        for cres in self.sensitivityvar:
            ofp.write('%-14s' %(cres))
            for cn in allnodes:
                ofp.write('%-14.5f' %(self.sensitivityvar[cres][cn]))
            ofp.write('\n')
        # write out the percent variance reduction values                    
        ofp.write('#'*10 + '   Percent Variance Reduction Values   ' + '#'*10 + '\n')
        ofp.write('%-14s' %('Response_node '))
        for cn in allnodes:
            ofp.write('%-14s' %(cn))
        ofp.write('\n')
        for cres in self.percentvarreduction:
            ofp.write('%-14s' %(cres))
            for cn in allnodes:
                ofp.write('%-14.5f' %(self.percentvarreduction[cres][cn]*100.0))
            ofp.write('\n')
        # write out the mutual information (Entropy) values
        ofp.write('#'*10 + '   Mutual Information (Entropy)   ' + '#'*10 + '\n')
        ofp.write('%-14s' %('Response_node '))
        for cn in allnodes:
            ofp.write('%-14s' %(cn))
        ofp.write('\n')
        for cres in self.sensitivityEntropy:
            ofp.write('%-14s' %(cres))
            for cn in allnodes:
                ofp.write('%-14.5f' %(self.sensitivityEntropy[cres][cn]))
            ofp.write('\n')

        # write out the normalized mutual information (Entropy) values
        ofp.write('#'*10 + '   Mutual Information (Entropy) Normalized  ' + '#'*10 + '\n')
        ofp.write('%-14s' %('Response_node '))
        for cn in allnodes:
            ofp.write('%-14s' %(cn))
        ofp.write('\n')
        for cres in self.sensitivityEntropyNorm:
            ofp.write('%-14s' %(cres))
            for cn in allnodes:
                ofp.write('%-14.5f' %(self.sensitivityEntropyNorm[cres][cn]))
            ofp.write('\n')
        ofp.close()



    def read_cas_file(self,casfilename):
        '''
        function to read in a casfile into a pynetica object.
        '''
        # first read in and strip out all comments and write out to a scratch file
        tmpdat = open(casfilename, 'r').readlines()
        ofp = open('###tmp###', 'w')
        for line in tmpdat:
            #line = re.sub('\?','*',line)
            if '//' not in line:
                ofp.write(line)
            elif line.strip().split()[0].strip() == '//':
                pass
            elif '//' in line:
                line = re.sub('//.*', '', line)
                if len(line.strip()) > 0:
                    ofp.write(line)
        ofp.close()
        self.casdata = np.genfromtxt('###tmp###', names=True,
                                     dtype=None, missing_values='*,?')
        os.remove('###tmp###')
        self.N = len(self.casdata)

    # cross validation driver
    def cross_val_setup(self):
        # open a file pointer to the stats output file for all the folds
        self.probpars.Val_outfile = '{0:s}_kfold_stats_VAL_{1:d}_folds.dat'.format(
            self.probpars.scenario.name, self.probpars.numfolds)
        kfoldOFP_Val = open(self.probpars.Val_outfile, 'w')
        kfoldOFP_Val.write(
            'Validation statistics for cross validation.\n' +
            'Base net --> {0:s} and casefile --> {1:s}\n'.format(self.probpars.baseNET,self.probpars.baseCAS) +
            'Current scenario is: {0:s}\n'.format(self.probpars.scenario.name))
        kfoldOFP_Val.write('%14s '*13
                  %('Current_Fold','Response','skillMean','rmseMean','meanErrMean','meanAbsErrMean',
                    'skillML','rmseML','meanErrML','meanAbsErrML','LogLoss','ErrorRate','QuadraticLoss')
                  + '\n')
        self.probpars.Cal_outfile = '%s_kfold_stats_CAL_%d_folds.dat' %(
            self.probpars.scenario.name,self.probpars.numfolds)
        kfoldOFP_Cal = open(self.probpars.Cal_outfile, 'w')
        kfoldOFP_Cal.write('Calibration statistics for cross validation.\nBase net --> %s and casefile --> %s\n'
                  %(self.probpars.baseNET, self.probpars.baseCAS) +
                  'Current scenario is: %s\n' %(self.probpars.scenario.name))
        kfoldOFP_Cal.write('%14s '*13
                %('Current_Fold','Response','skillMean','rmseMean','meanErrMean','meanAbsErrMean',
                  'skillML','rmseML','meanErrML','meanAbsErrML','LogLoss','ErrorRate','QuadraticLoss')
                + '\n')

        for cfold in np.arange(self.probpars.numfolds):
            self.allfolds.calNODES.append(None)
            self.allfolds.valNODES.append(None)
            self.allfolds.calpred.append(None)
            self.allfolds.valpred.append(None)

            cname = '{0:s}_fold_{1:d}.cas'.format(self.probpars.scenario.name, cfold)
            self.allfolds.casfiles.append(cname)
            retinds = np.array(self.allfolds.retained[cfold], dtype=int)
            # outdat only includes the columns that are in CASheader
            outdat = np.atleast_2d(self.casdata[self.probpars.CASheader[0]][retinds]).T
            # caldata and valdata both include all columns for simplicity
            self.allfolds.caldata.append(self.casdata[retinds])
            leftoutinds = np.array(self.allfolds.leftout[cfold], dtype=int)
            # we will also make a CAS file for the leftout data for using Netica's testing codes
            outdatLeftOut = np.atleast_2d(self.casdata[self.probpars.CASheader[0]][leftoutinds]).T
            self.allfolds.valdata.append(self.casdata[leftoutinds])
            self.allfolds.valN.append(len(leftoutinds))
            self.allfolds.calN.append(len(retinds))
            # concatenate together the columns of data that will make up the CAS files
            for i, chead in enumerate(self.probpars.CASheader):
                if i > 0:
                    outdat = np.hstack((outdat, np.atleast_2d(self.casdata[chead][retinds]).T))
                    outdatLeftOut = np.hstack((outdatLeftOut, np.atleast_2d(self.casdata[chead][leftoutinds]).T))
            # write out the retained casefile            
            ofp = open(cname, 'w')
            for cnode in self.probpars.CASheader:
                ofp.write('{0:s} '.format(cnode))
            ofp.write('\n')
            for line in outdat:
                for cv in line:
                    if np.isnan(cv):
                        ofp.write(' {0:20s} '.format('*'))
                    else:
                        ofp.write(' {0:20.8e} '.format(cv))
                ofp.write('\n')
            ofp.close()

            # write out the leftout casefile for later use with the Netica validation testing functions
            ofpLeftOut = open(cname[:-4] + '_leftout.cas', 'w')
            for cnode in self.probpars.CASheader:
                ofpLeftOut.write('{0:s} '.format(cnode))
            ofpLeftOut.write('\n')
            np.savetxt(ofpLeftOut, outdatLeftOut)
            ofpLeftOut.close()
            ofpLeftOut.close()
        return kfoldOFP_Val, kfoldOFP_Cal



