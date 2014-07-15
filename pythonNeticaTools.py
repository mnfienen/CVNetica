import numpy as np
import os
import ctypes as ct
import platform
import pythonNeticaConstants as pnC
import cthelper as cth


class statestruct:
    def __init__(self):
        self.obj = None
        self.name = None
        self.numeric = None

class experience:
    def __init__(self):
        self.parent_names = list()
        self.parent_states = None
        self.node_experience = list()
        self.case_experience = list()

class nodestruct:
    def __init__(self):
        self.name = None
        self.title = None
        self.beliefs = None
        self.Nbeliefs = None
        self.Likelihood = None
        self.continuous = False
        self.state = []

class pyneticaTools:
    def __init__(self):
        self.n = None
        self.mesg = ct.create_string_buffer('\000' * 1024)
        self.env = None
        
    def sanitize(self):
        print 'Sanitizing pynetica object to remove pointers'
        # code to strip out all ctypes information from SELF to 
        # allow for pickling
        self.n = None
        self.mesg = None
        self.env = None
        
    def start_environment(self, licfile):
        # read in the license file information
        self.licensefile = licfile
        if os.path.exists(self.licensefile):
            self.license = open(self.licensefile, 'r').readlines()[0].strip().split()[0]
        else:
            print ("Warning: License File [{0:s}] not found.\n".format(self.licensefile) +
                   "Opening Netica without licence, which will limit size of nets that can be used.\n" +
                   "Window may become unresponsive.")
            self.license = None         
        self.NewNeticaEnviron()
            
    #############################################
    # Major validation and prediction functions #
    #############################################
    def rebuild_net(self, NetName, newCaseFile, voodooPar, outfilename, EMflag=False):
        '''
         rebuild_net(NetName,newCaseFilename,voodooPar,outfilename)
         a m!ke@usgs joint <mnfienen@usgs.gov>
         function to build the CPT tables for a new CAS file on an existing NET
         (be existing, meaning that the nodes, edges, and bins are dialed)
         INPUT:
               NetName --> a filename, including '.neta' extension
               newCaseFilename --> new case file including '.cas' extension
               voodooPar --> the voodoo tuning parameter for building CPTs
               outfilename --> netica file for newly build net (including '.neta')
               EMflag --> if True, use EM to learn from casefile, else (default)
                         incorporate the CPT table directly
         '''   
        # create a Netica environment
        print 'Rebuilding net: {0:s} using Casefile: {1:s}'.format(NetName, newCaseFile)
        # make a streamer to the Net file
        net_streamer = self.NewFileStreamer(NetName)
        # read in the net using the streamer        
        cnet = self.ReadNet(net_streamer)
        # remove the input net streamer
        self.DeleteStream(net_streamer)  
        self.CompileNet(cnet)      
        #get the nodes and their number
        allnodes = self.GetNetNodes(cnet)
        numnodes = self.LengthNodeList(allnodes)

        # loop over the nodes deleting CPT
        for cn in np.arange(numnodes):
            cnode = self.NthNode(allnodes,cn)
            self.DeleteNodeTables(cnode)
        # make a streamer to the new cas file
        new_cas_streamer = self.NewFileStreamer(newCaseFile)

        if EMflag:
            print 'Learning new CPTs using EM algorithm'
            # to use EM learning, must first make a learner and set a couple options
            newlearner = self.NewLearner(pnC.learn_method_bn_const.EM_LEARNING)
            self.SetLearnerMaxTol(newlearner, 1.0e-6)
            self.SetLearnerMaxIters(newlearner, 1000)
            # now must associate the casefile with a caseset (weighted by unity)
            newcaseset = self.NewCaseset('currcases')
            self.AddFileToCaseset(newcaseset, new_cas_streamer, 1.0)
            self.LearnCPTs(newlearner, allnodes, newcaseset, voodooPar)
            self.DeleteCaseset(newcaseset)
            self.DeleteLearner(newlearner)

        else:
            print 'Learning new CPTs using ReviseCPTsByCaseFile'
            self.ReviseCPTsByCaseFile(new_cas_streamer, allnodes, voodooPar)
        outfile_streamer = self.NewFileStreamer(outfilename)
        self.CompileNet(cnet)


        outfile_streamer = self.NewFileStreamer(outfilename)
        print 'Writing new net to: %s' %(outfilename)
        self.WriteNet(cnet,outfile_streamer)
        self.DeleteStream(outfile_streamer)  
        self.DeleteNet(cnet)

    def OpenNeticaNet(self,netName):
        '''
        Open a net identified by netName.
        Returns a pointer to the opened net after it is compiled
        '''
        # meke a streamer to the Net file
        cname = netName
        if '.neta' not in netName:
            cname += '.neta'
        net_streamer = self.NewFileStreamer(cname)
        # read in the net using the streamer        
        cnet = self.ReadNet(net_streamer)
        # remove the input net streamer
        self.DeleteStream(net_streamer)  
        self.CompileNet(cnet)      
        return cnet

    def ReadNodeInfo(self, netName):
        '''
        Read in all information on beliefs, states, and likelihoods for all 
        nodes in the net called netName
        '''
        # open the net stored in netName
        cnet = self.OpenNeticaNet(netName)
        #get the nodes and their number
        allnodes = self.GetNetNodes(cnet)
        numnodes = self.LengthNodeList(allnodes)
        print 'Reading Node information from net --> {0:s}'.format(netName)
        cNETNODES = dict()
        # loop over the nodes
        for cn in np.arange(numnodes):
            cnode = self.NthNode(allnodes, cn)
            cnodename = cth.c_char_p2str(self.GetNodeName(cnode))
            cNETNODES[cnodename] = nodestruct()
            cNETNODES[cnodename].name = cth.c_char_p2str(self.GetNodeName(cnode))
            cNETNODES[cnodename].title = cth.c_char_p2str(self.GetNodeTitle(cnode))
            print '   Parsing node --> %s' %(cNETNODES[cnodename].title)
            cNETNODES[cnodename].Nbeliefs = self.GetNodeNumberStates(cnode)
            cNETNODES[cnodename].beliefs = cth.c_float_p2float(
                self.GetNodeBeliefs(cnode),
                cNETNODES[cnodename].Nbeliefs)
            cNETNODES[cnodename].likelihood = cth.c_float_p2float(
                self.GetNodeLikelihood(cnode),
                cNETNODES[cnodename].Nbeliefs)
            cNETNODES[cnodename].levels = cth.c_double_p2float(
                self.GetNodeLevels(cnode),
                cNETNODES[cnodename].Nbeliefs + 1)

            # loop over the states in each node
            for cs in range(cNETNODES[cnodename].Nbeliefs):
                cNETNODES[cnodename].state.append(statestruct())
                cNETNODES[cnodename].state[-1].name = cth.c_char_p2str(
                    self.GetNodeStateName(cnode,cs))   
       
        self.DeleteNet(cnet)
        return cNETNODES
    


    def ConfusionMatrix(self,ctester,cnode):
        '''
        Makes a confusion matrix for a particular node specified by name in cnode
        within the tester environment laid out in ctester
        '''
        numstates = self.GetNodeNumberStates(cnode)
        
        confusion_matrix = np.zeros((numstates,numstates))
        for a in np.arange(numstates):
            for p in np.arange(numstates):
                confusion_matrix[a,p] = self.GetTestConfusion(ctester,cnode,p,a)
        return confusion_matrix

    def ExperienceAnalysis(self,cn,cnet):
        '''
        calculate the experience for the node named in cn
        '''
        cnex = experience()
        # get a list of the parents of the node
        testnode = self.GetNodeNamed(cn,cnet)
        #start a list for the cartesian sum of node states
        allstates = list()
        cparents = self.GetNodeParents(testnode)    
        numnodes = self.LengthNodeList(cparents)
        for cp in np.arange(numnodes):
            # append the name to the list of returned names
            cnode = self.NthNode(cparents,cp)
            cnex.parent_names.append(cth.c_char_p2str(self.GetNodeName(cnode)))
            # find the number of states for each parent
            allstates.append(np.arange(self.GetNodeNumberStates(
                self.NthNode(cparents,cp))))
        if numnodes > 1:
            cnex.parent_states = self.cartesian(allstates)
        else:
            cnex.parent_states = allstates
        for cs in cnex.parent_states:
            cnex.node_experience.append(self.GetnodeExperience(
                testnode,cs.ctypes.data_as(ct.POINTER(ct.c_int))))
        cnex.node_experience = np.array(cnex.node_experience)
        # change the null pointers (meaning 
        cnex.node_experience[cnex.node_experience<1]=0.0
        
        return cnex

    ###################################
    # Key helper functions for Netica #   
    ###################################

    def NewNeticaEnviron(self):
        '''
        create a new Netica environment based on operating system
        '''
        # first access the .dll or .so
        try:
            if 'window' in platform.system().lower():
                self.n = ct.windll.Netica
            else:
                self.n = ct.cdll.LoadLibrary("./libnetica.so")
        except:
            raise(dllFail(platform.system()))
        #next try to establish an environment for Netica
        self.env = ct.c_void_p(self.n.NewNeticaEnviron_ns(self.license, None, None))
        # try to intialize Netica
        res = self.n.InitNetica2_bn(self.env, ct.byref(self.mesg))
        if res >= 0:
            print '\n'*2 + '#' * 40 + '\nOpening Netica:'
            print self.mesg.value
        else:
            raise(NeticaInitFail(res.value))    
        print 'Netica is open\n' + '#'*40 + '\n' * 2
        
    def CloseNetica(self):
        res = self.n.CloseNetica_bn(self.env, ct.byref(self.mesg))    
        if res >= 0:
            print "Closing Netica:"
            print self.mesg.value
        else:
            raise(NeticaCloseFail(res.value))    
        self.n = None
        
    def GetError(self, severity = pnC.errseverity_ns_const.ERROR_ERR, after = None):
        res = self.n.GetError_ns(self.env, severity, after)
        if res: return ct.c_void_p(res)
        else:   return None

    def ErrorMessage(self, error):
        return self.n.ErrorMessage_ns(error)

    # general error-checking function    
    def chkerr(self,err_severity = pnC.errseverity_ns_const.ERROR_ERR):
        if self.GetError(err_severity):
            exceptionMsg = ("\npythonNeticaUtils: \nError " + 
                            str(ct.cast(ct.c_void_p(self.ErrorMessage(self.GetError(err_severity))), ct.c_char_p).value))
            raise NeticaException(exceptionMsg)

    ################################################################
    # Small definitions and little functions in alphabetical order #  
    ################################################################   
    def AddFileToCaseset(self,caseset,streamer,degree):
        self.n.AddFileToCaseset_cs(caseset,streamer,ct.c_double(degree),None)
        self.chkerr()

       
    def cartesian(self,arrays,out=None):   
        '''
        function to calculate the Cartesian sum of multiple arrays.
        This is used to provide the permutations (odometer style) of all
        the possible parent states when calculating experience.
        See: http://stackoverflow.com/questions/1208118/
        using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
        '''

        arrays = [np.asarray(x) for x in arrays]
        dtype = arrays[0].dtype
    
        n = np.prod([x.size for x in arrays])
        if out is None:
            out = np.zeros([n, len(arrays)], dtype=dtype)
    
        m = n / arrays[0].size
        out[:,0] = np.repeat(arrays[0], m)
        if arrays[1:]:
            # recursive?
            self.cartesian(arrays[1:], out=out[0:m,1:])
            for j in xrange(1, arrays[0].size):
                out[j*m:(j+1)*m,1:] = out[0:m,1:]
        return out

    def CompileNet(self, net):
        self.n.CompileNet_bn(net)
        self.chkerr()

    def CopyNet(self,oldnet, newnetname,options):
        newnet = self.n.CopyNet_bn(oldnet,newnetname,self.env,options)
        self.chkerr()
        return newnet

    def CopyNodes(self,oldnodes,newnet,options):
        newnodes = self.n.CopyNodes_bn(oldnodes,newnet,options)
        self.chkerr()
        return newnodes 

    def DeleteCaseset(self,caseset):
        self.n.DeleteCaseset_cs(caseset)
        self.chkerr()

    def DeleteLearner(self,newlearner):
        self.n.DeleteLearner_bn(newlearner)
        self.chkerr()

    def DeleteNet(self,cnet):
        self.n.DeleteNet_bn(cnet)
        self.chkerr()
    
    def DeleteNetTester(self,ctester):
        self.n.DeleteNetTester_bn(ctester)
        self.chkerr()

    def DeleteNodeTables(self,node):
        self.n.DeleteNodeTables_bn(node)
        self.chkerr()

    def DeleteNodeList(self,cnodes):
        self.n.DeleteNodeList_bn(cnodes)
        self.chkerr()

    def DeleteStream(self,cstream):
        self.n.DeleteStream_ns(cstream)
        self.chkerr()

    def DeleteSensvToFinding(self,sens):
        self.n.DeleteSensvToFinding_bn(sens)
        self.chkerr()

    def EnterFinding(self,cnode,cval):
        self.n.EnterFinding_bn(cnode,ct.c_double(cval))
        self.chkerr()

    def EnterNodeValue(self,cnode,cval):
        self.n.EnterNodeValue_bn(cnode,ct.c_double(cval))
        self.chkerr()
        
    def GetMutualInfo(self,sensentrop,Vnode):
        tmpNeticaFun = self.n.GetMutualInfo_bn
        tmpNeticaFun.restype=ct.c_double
        retvar = self.n.GetMutualInfo_bn(sensentrop,Vnode)
        self.chkerr()
        return retvar        

    def GetNetNodes(self,cnet):
        allnodes = self.n.GetNetNodes2_bn(cnet,None)
        self.chkerr()
        return allnodes

    def GetNodeBeliefs(self,cnode):
        beliefs = self.n.GetNodeBeliefs_bn(cnode)
        self.chkerr()
        return beliefs

    def GetNodeExpectedValue(self,cnode):
        std_dev = ct.c_double()
        # make a temporary function variable to be able to set the
        # return value
        tmpNeticaFun = self.n.GetNodeExpectedValue_bn
        tmpNeticaFun.restype=ct.c_double
        expected_val = tmpNeticaFun(cnode,ct.byref(std_dev),
                                    None,None)
        self.chkerr()
        return expected_val, std_dev.value

    def GetnodeExperience(self,cnode,parent_states):
        tmpNeticaFun = self.n.GetNodeExperience_bn
        tmpNeticaFun.restype=ct.c_double
        experience = tmpNeticaFun(cnode,parent_states)
        self.chkerr()
        return experience
    
    def GetNodeFinding(self,cnode):
        cf = self.n.GetNodeFinding_bn(cnode)
        self.chkerr()
        return cf
    
    def GetNodeLevels(self,cnode):
        nodelevels = self.n.GetNodeLevels_bn(cnode)
        self.chkerr()
        return nodelevels

    def GetNodeLikelihood(self,cnode):
        nodelikelihood = self.n.GetNodeLikelihood_bn(cnode)
        self.chkerr()
        return nodelikelihood

    def GetNodeName(self,cnode):
        cname = self.n.GetNodeName_bn(cnode)
        self.chkerr()
        return cname

    def GetNodeNamed(self,nodename,cnet):
        retnode = self.n.GetNodeNamed_bn(nodename,cnet)
        self.chkerr()
        return(retnode)
    
    def GetNodeNumberStates(self,cnode):
        numstates = self.n.GetNodeNumberStates_bn(cnode)
        self.chkerr()
        return numstates

    def GetNodeParents(self,cnode):
        parents = self.n.GetNodeParents_bn(cnode)
        self.chkerr()
        return parents
    
    def GetNodeStateName(self,cnode,cstate):
        stname = self.n.GetNodeStateName_bn(cnode,ct.c_int(cstate))
        self.chkerr()
        return stname

    def GetNodeTitle(self,cnode):
        ctitle = self.n.GetNodeTitle_bn(cnode)
        self.chkerr()
        return ctitle

    def GetTestLogLoss(self,ctester,cnode):
        tmpNeticaFun = self.n.GetTestLogLoss_bn
        tmpNeticaFun.restype=ct.c_double
        logloss = self.n.GetTestLogLoss_bn(ctester,cnode)
        self.chkerr()
        return logloss
    
    def GetTestConfusion(self,ctester,cnode,predState,actualState):
        tmpNeticaFun = self.n.GetTestConfusion_bn
        tmpNeticaFun.restype=ct.c_double        
        confusion = tmpNeticaFun(ctester,cnode,ct.c_int(predState),ct.c_int(actualState))
        self.chkerr()
        return confusion
    
    def GetTestErrorRate(self,ctester,cnode):
        tmpNeticaFun = self.n.GetTestErrorRate_bn
        tmpNeticaFun.restype=ct.c_double
        errrate = tmpNeticaFun(ctester,cnode)
        self.chkerr()
        return errrate
    
    def GetTestQuadraticLoss(self,ctester,cnode):
        tmpNeticaFun = self.n.GetTestQuadraticLoss_bn
        tmpNeticaFun.restype=ct.c_double
        quadloss = self.n.GetTestQuadraticLoss_bn(ctester,cnode)
        self.chkerr()
        return quadloss
        
    def GetVarianceOfReal(self,sensv,Vnode):
        tmpNeticaFun = self.n.GetVarianceOfReal_bn
        tmpNeticaFun.restype=ct.c_double
        retvar = self.n.GetVarianceOfReal_bn(sensv,Vnode)
        self.chkerr()
        return retvar
        
    def LearnCPTs(self,learner,nodes,caseset,voodooPar):
        self.n.LearnCPTs_bn(learner,nodes,caseset,ct.c_double(voodooPar))
        self.chkerr()

    def LengthNodeList(self, nodelist):
        res = self.n.LengthNodeList_bn(nodelist)
        self.chkerr()
        return res    

    def LimitMemoryUsage(self, memlimit):
        self.n.LimitMemoryUsage_ns(ct.c_double(memlimit),self.env)
        print 'set memory limit to ---> %f bytes' %memlimit
        self.chkerr()
        
    def NewCaseset(self,name):
        newcaseset = self.n.NewCaseset_cs(name,self.env)
        self.chkerr()
        return newcaseset

    def NewFileStreamer(self,infile):
        streamer =  self.n.NewFileStream_ns (infile, self.env,None)
        self.chkerr()
        return streamer

    def NewLearner(self,method):
        newlearner = self.n.NewLearner_bn(method,None,self.env)
        self.chkerr()
        return newlearner

    def NewNet(self, netname):
        newnet = self.n.NewNet_bn(netname,self.env)
        self.chkerr()
        return newnet
    
    def NewNetTester(self,test_nodes,unobs_nodes):
        tester = self.n.NewNetTester_bn(test_nodes,unobs_nodes,ct.c_int(-1))
        self.chkerr()
        return tester

    def NewNodeList2(self,length,cnet):
        nodelist = self.n.NewNodeList2_bn(ct.c_int(length),cnet)
        self.chkerr()
        return nodelist
    
    def NewSensvToFinding(self,Qnode,Vnodes,what_find):
        sensv = self.n.NewSensvToFinding_bn(Qnode,Vnodes,what_find)
        self.chkerr()
        return sensv

    def NthNode(self,nodelist,index_n):
        cnode = self.n.NthNode_bn(nodelist,ct.c_int(index_n))
        self.chkerr()
        return cnode
        
    def ReadNet(self,streamer):
        cnet = self.n.ReadNet_bn(streamer,ct.c_int(pnC.netica_const.NO_WINDOW))
        # check for errors
        self.chkerr()
        # reset the findings
        self.n.RetractNetFindings_bn(cnet)
        self.chkerr()
        return cnet

    def ReadNetFindings2(self,case_posn,filestream,allnodes):
        #tmpneticafun = self.n.ReadNetFindings2_bn
        #tmpneticafun.restype = ct.c_long
        case_position = ct.c_long(case_posn)
        self.n.ReadNetFindings2_bn(ct.byref(case_position),
                            filestream,
                            ct.c_ubyte(0),
                            allnodes,None,None)
        self.chkerr()      
        return case_position.value
                                   
                                   
    def RetractNetFindings(self,cnet):
        self.n.RetractNetFindings_bn(cnet)
        self.chkerr()

    def ReviseCPTsByCaseFile(self,casStreamer,cnodes,voodooPar):
        self.n.ReviseCPTsByCaseFile_bn(casStreamer,cnodes,ct.c_int(0),
                                       ct.c_double(voodooPar))
        self.chkerr()

    def SetLearnerMaxIters(self,learner,maxiters):
        self.n.SetLearnerMaxIters_bn(learner,ct.c_int(maxiters))
        self.chkerr()    

    def SetLearnerMaxTol(self,learner,tol):
        self.n.SetLearnerMaxTol_bn(learner,ct.c_double(tol))
        self.chkerr()         
        
    def SetNetAutoUpdate(self,cnet,belief_value):
        self.n.SetNetAutoUpdate_bn(cnet,belief_value)
        self.chkerr()

    def SetNthNode(self, nodelist, position, cnode):
        self.n.SetNthNode_bn(nodelist, ct.c_int(position), cnode)
        self.chkerr()

    def SetNodeLevels(self, cnode, clevels):
        self.n.SetNodeLevels_bn(cnode, ct.c_int(len(clevels)-1),
                                                clevels.ctypes.data_as(ct.POINTER(ct.c_double)))
        self.chkerr()

    def TestWithCaseset(self, test, cases):
        self.n.TestWithCaseset_bn(test, cases)
        self.chkerr()
        
    def WriteNet(self, cnet, filename_streamer):
        self.n.WriteNet_bn(cnet, filename_streamer)
        self.chkerr()

#################
# Error Classes #
#################
# -- can't open external file
class dllFail(Exception):
    def __init__(self,cplat):
        self.cplat = cplat
    def __str__(self):
        if "windows" in self.cplat.lower():
            return("\n\nCannot open Netica.dll.\nBe sure it's in the path")
        else:
            return("\n\nCannot open libnetica.so.\nBe sure it's in the path")
# -- can't initialize Netica
class NeticaInitFail(Exception):
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return("\n\nCannot initialize Netica. Netica message is:\n%s\n" 
               %(self.msg))
# -- can't close Netica
class NeticaCloseFail(Exception):
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return("\n\nCannot properly close Netica. Netica message is:\n%s\n" 
               %(self.msg))
# -- General Netica Exception
class NeticaException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg