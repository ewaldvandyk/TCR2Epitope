import inspect
import os
import itertools
import numpy.matlib as matlib
import numpy as np
import ctypes as ct
import warnings

class Nmers:
    def __init(self):
        self.strList = list()
        self.posList = list()

class Gmm_model:
    def __init__(self, dim=None, numMix=None):
        if not isinstance(dim, (int, long)) or not isinstance(numMix, (int, long)):
            raise IOError("Both dim and numMix must be specified")
        self.dim = dim
        self.numMix = numMix
        self.mixPriors = None
        self.gaussCoeffs = None
        self.mus = None
        self.sig_inv = None
   
    def set_mixPriors(self, mixPriors=None):
        if not type(mixPriors).__module__ == "numpy.matrixlib.defmatrix":
            raise TypeError("mixPriors not of type numpy.matrixlib.defmatrix")
        (priorFold,numMix) = mixPriors.shape 
        if numMix != self.numMix:
            raise IndexError("Length of mixPriors doesn't match the number of mixtures") 
        if priorFold != 1:
            raise IndexError("mixPrior cannot contain multiple rows")
        self.mixPriors = mixPriors
            
    def set_gaussCoeffs(self, gaussCoeffs=None):
        if not type(gaussCoeffs).__module__ == "numpy.matrixlib.defmatrix":
            raise TypeError("gaussCoeffs not of type numpy.matrixlib.defmatrix")
        (gaussFold,numMix) = gaussCoeffs.shape
        if numMix != self.numMix:
            raise IndexError("Length of gaussCoeffs doesn't match the number of mixtures") 
        if gaussFold != 1:
            raise IndexError("gaussCoeffs cannot contain multiple rows")
        self.gaussCoeffs = gaussCoeffs
    
    def set_mus(self, mus=None):
        if not type(mus).__module__ == 'numpy.matrixlib.defmatrix':
            raise TypeError("Mu matrix not of type numpy.matrixlib.defmatrix")
        (nMix, dim) = mus.shape
        if dim != self.dim:
            raise IndexError("Mu vector does not have the proper dimension")
        if nMix != self.numMix:
            raise IndexError("Mu matrix does not have the proper number of mixtures")
        self.mus = mus
        
    def set_sig_inv(self, sig_inv=None):
        if not isinstance(sig_inv, list):
            raise TypeError("Gaussian sig_inv matrices not passed in list format")
        numMix = len(sig_inv)
        if numMix != self.numMix:
            raise IndexError("sig_inv does not have the proper number of mixtures")
        for mI in range(numMix):
            if not type(sig_inv[mI]).__module__ == "numpy.matrixlib.defmatrix":
                raise TypeError("sig_inv mixture %d not in numpy.matrixlib.defmatrix format"%(mI))
            (dimR, dimC) = sig_inv[mI].shape
            if dimR != self.dim or dimC != self.dim:
                raise IndexError("sig_inv mixture %d doesn't have the proper dimensions"%(mI))
        self.sig_inv = sig_inv
            
    def vec_2_posterior(self, vec=None):
        if self.mixPriors == None or self.gaussCoeffs == None or \
           self.mus == None or self.sig_inv == None:
            raise IOError("Model parameters not loaded properly")
        if not type(vec).__module__ == "numpy.matrixlib.defmatrix":
            raise TypeError("Input vector not of type numpy.matrixlib.defmatrix")
        
        numVec, vecLen = vec.shape
        if numVec != 1:
            raise IndexError("Only one vector allowed at a time")
        if vecLen != self.dim:
            raise IndexError("nMer vector dimension doesn't agree with gmm dimension")
        gaussP  = matlib.zeros(shape=[1, self.numMix])
        for bI in range(self.numMix):
            x = vec - self.mus[bI,:]
            x_sqr = x*self.sig_inv[bI]*x.transpose()
            pwr = -0.5*x_sqr
            gaussP[0,bI] = self.gaussCoeffs[0,bI] * np.exp(pwr)
        
        postNorm = 1/(gaussP*self.mixPriors.transpose())
        
        postVec = postNorm[0,0]*matlib.multiply(gaussP, self.mixPriors) 
#         print postVec
        return postVec

class Gmm_model_c:
    gmm_obj_file_name = "gmm.so"
    def __init__(self, dim=None, numMix=None, c_object_file=None):
        if not isinstance(dim, (int, long)) or not isinstance(numMix, (int, long)):
            raise IOError("Both dim and numMix must be specified")
        self.dim = dim
        self.numMix = numMix
        
        self.mixPriorFlag       = False
        self.gaussCoeffsFlag    = False
        self.musFlag            = False
        self.sigInvFlag         = False
            
        if c_object_file == None or not os.path.isfile(c_object_file):
            classFileName = inspect.getfile(self.__class__)
            classDir = os.path.dirname(classFileName)
            so_file = os.path.join(classDir, Gmm_model_c.gmm_obj_file_name)
            self.c_lib=ct.cdll.LoadLibrary(so_file)
        else:
            self.c_lib=ct.cdll.LoadLibrary(c_object_file)
        self.c_lib.init(self.numMix,self.dim)
        
    def __del__(self):
        self.c_lib.clean()
        
    def set_mixPriors(self, mixPriors=None):
        if not type(mixPriors).__module__ == "numpy.matrixlib.defmatrix":
            raise TypeError("mixPriors not of type numpy.matrixlib.defmatrix")
        (priorFold,numMix) = mixPriors.shape 
        if numMix != self.numMix:
            raise IndexError("Length of mixPriors doesn't match the number of mixtures") 
        if priorFold != 1:
            raise IndexError("mixPrior cannot contain multiple rows")
        cMixPriors = np.ctypeslib.as_ctypes(mixPriors.flatten())
        self.c_lib.set_mixPriors(cMixPriors)
        self.mixPriorFlag = True
        
    def set_gaussCoeffs(self, gaussCoeffs=None):
        if not type(gaussCoeffs).__module__ == "numpy.matrixlib.defmatrix":
            raise TypeError("gaussCoeffs not of type numpy.matrixlib.defmatrix")
        (gaussFold,numMix) = gaussCoeffs.shape
        if numMix != self.numMix:
            raise IndexError("Length of gaussCoeffs doesn't match the number of mixtures") 
        if gaussFold != 1:
            raise IndexError("gaussCoeffs cannot contain multiple rows")
        cGaussCoeffs = np.ctypeslib.as_ctypes(gaussCoeffs.flatten())
        self.c_lib.set_gaussCoeffs(cGaussCoeffs)
        self.gaussCoeffsFlag = True
    
    def set_mus(self, mus=None):
        if not type(mus).__module__ == 'numpy.matrixlib.defmatrix':
            raise TypeError("Mu matrix not of type numpy.matrixlib.defmatrix")
        (nMix, dim) = mus.shape
        if dim != self.dim:
            raise IndexError("Mu matrix does not have the proper dimension")
        if nMix != self.numMix:
            raise IndexError("Mu matrix does not have the proper number of mixtures")
        cMus = np.ctypeslib.as_ctypes(mus.flatten())
        self.c_lib.set_mus(cMus)
        self.musFlag = True

    def set_sig_inv(self, sig_inv=None):
        if not isinstance(sig_inv, list):
            raise TypeError("Gaussian sig_inv matrices not passed in list format")
        numMix = len(sig_inv)
        if numMix != self.numMix:
            raise IndexError("sig_inv does not have the proper number of mixtures")
        for mI in range(numMix):
            if not type(sig_inv[mI]).__module__ == "numpy.matrixlib.defmatrix":
                raise TypeError("sig_inv mixture %d not in numpy.matrixlib.defmatrix format"%(mI))
            (dimR, dimC) = sig_inv[mI].shape
            if dimR != self.dim or dimC != self.dim:
                raise IndexError("sig_inv mixture %d doesn't have the proper dimensions"%(mI))
        
        new_sigInv = np.zeros(shape = self.dim*self.dim*self.numMix)
        sigStart = 0
        sigEnd = self.dim*self.dim      
        for mat in sig_inv:
            new_sigInv[sigStart:sigEnd] = mat.flatten()
            sigStart    += self.dim*self.dim
            sigEnd      += self.dim*self.dim
        
        cSig_inv = np.ctypeslib.as_ctypes(new_sigInv)
        self.c_lib.set_sig_inv(cSig_inv)
#         print sig_inv[120]
        self.sigInvFlag = True
        
    def vec_2_posterior(self, vec=None):
        if self.mixPriorFlag == False or \
        self.gaussCoeffsFlag == False or \
        self.musFlag  == False or \
        self.sigInvFlag == False:
            raise IOError("Model parameters not loaded properly")

        # Convert to np.array because as_ctypes 
        # has a bug and does not process flattened 
        # matrix well
        vec2 = np.array([vec[0,i] for i in range(self.dim)]) 
        cVec     = np.ctypeslib.as_ctypes(vec2) 
        postVec = np.zeros(self.numMix)
        cPostVec  = np.ctypeslib.as_ctypes(postVec)
        self.c_lib.c_vec_2_posterior(cVec,cPostVec)
        postMat = matlib.matrix(postVec)
        return postMat
        
class AaVecSpace:
    nMerLenRange = [1, 6]
    nMerSpanRange = [1,8]
    cdrSeqLenRange = [3, 1000]
    
    def __init__(self, symVecFile=None, nMerLen=3, nMerSpan=4):
        
        classFileName = inspect.getfile(self.__class__)
        classDir = os.path.dirname(classFileName)
        self._dataDir = os.path.abspath(os.path.join(classDir, "..", "..", "data"))
        if symVecFile == None:
            symVecFile = os.path.join(self._dataDir, 'pmbecAAvecs.tsv')
        
        if not os.path.isfile(symVecFile):
            raise IOError("Symbol vector file file not found")
        self._read_aaVecs(symVecFile)
        self.set_nMerLen(nMerLen)
        self.set_nMerSpan(nMerSpan)
        
    
    def _read_aaVecs(self, vecFile):
        self._aaVecs = dict();
        vecDim = 0;
        with open(vecFile, 'rU') as f:
            for line in f:
                elems = line.strip().split("\t")
                self._aaVecs[elems[0]] = map(float, elems[1:])
                currVecDim = len(self._aaVecs[elems[0]])
                if vecDim != 0 and vecDim != currVecDim:
                    raise IOError("Dimension of vectors not consistent within symbol vector file") 
                vecDim = currVecDim 
        self._vecDim = vecDim
        try:
            self._nMerLen
        except AttributeError:
            self._nMerVecDim = None
        else:
            self._nMerVecDim = (self._vecDim+1)*self._nMerLen
        
    def is_proper_cdr_seq(self, cdrSeq):
        if not isinstance(cdrSeq, basestring):
            return False
        cdrSeqLen = len(cdrSeq)
        if cdrSeqLen < AaVecSpace.cdrSeqLenRange[0] or cdrSeqLen > AaVecSpace.cdrSeqLenRange[1]:
            return False
        if not set(cdrSeq) <= set(self._aaVecs.keys()):
            return False
        return True
    
    def is_proper_nMer_params(self, cdrSeqLen=None):
        if not isinstance(cdrSeqLen, (int, long)):
            return False
        if self._nMerLen > self._nMerSpan:
            return False
        if cdrSeqLen < self._nMerLen:
            return False
        return True
        
    def set_nMerLen(self, nMerLen):
        if not isinstance(nMerLen, (int, long)):
            raise IOError("nMerLen must be an integer")
        if nMerLen < AaVecSpace.nMerLenRange[0] or nMerLen > AaVecSpace.nMerLenRange[1]:
            raise IOError("nMerLen must be between %d and %d"%(AaVecSpace.nMerLenRange[0], AaVecSpace.nMerLenRange[1]))
        self._nMerLen = nMerLen
        try:
            self._vecDim
        except AttributeError:
            self._nMerVecDim = None
        else:
            self._nMerVecDim = (self._vecDim+1)*self._nMerLen
                  
    def set_nMerSpan(self, nMerSpan):
        if not isinstance(nMerSpan, (int, long)):
            raise IOError("nMerSpan must be an integer")
        if nMerSpan < AaVecSpace.nMerSpanRange[0] or nMerSpan > AaVecSpace.nMerSpanRange[1]:
            raise IOError("nMerSpan must be between %d and %d"%(AaVecSpace.nMerSpanRange[0], AaVecSpace.nMerSpanRange[1]))        
        self._nMerSpan = nMerSpan
                
    def get_nmer_aaVec(self, nMer=None):
        if nMer == None:
            raise IOError("nMer string not specified")
        n = len(nMer)
        nMerVec = [0] * (n*self._vecDim)
        for aaI, aa in enumerate(nMer):
            aaStart = aaI*self._vecDim
            aaEnd = aaStart+self._vecDim
            nMerVec[aaStart:aaEnd] = self._aaVecs[aa]
        return nMerVec
    
    def get_nmer_aaPos(self, nMerI=None, cdrLen=None):
        if nMerI == None:
            raise IOError("nMer position indices not specified")
        if cdrLen == None:
            raise IOError("Cdr length not specified")
        
        return [float(i)/(cdrLen-1) for i in nMerI]
    
    def get_nmer_vec(self, nMer=None, nMerI=None, cdrLen = None):
        nmer_aaVec = self.get_nmer_aaVec(nMer)
        nmer_aaPos = self.get_nmer_aaPos(nMerI, cdrLen)
#         nmer_vec = nmer_aaVec
        nmer_vec = nmer_aaVec + nmer_aaPos
        return nmer_vec
        
    def get_cdr_nmers(self, cdrSeq=None):
        nmers = Nmers()
        if not self.is_proper_cdr_seq(cdrSeq):
            raise IOError("CDR sequence not specified correctly.")
        cdrSeqLen = len(cdrSeq)
        if not self.is_proper_nMer_params(cdrSeqLen):
            raise IOError("Invalid combination for nMerLen, nMerSpan and CDR length")
        
        nMerFirstStart = 0;
        nMerFirstEnd = cdrSeqLen - self._nMerLen + 1
        
        nMerStrList = list()
        nMerPosList = list()
        for nMerStart in range(nMerFirstStart, nMerFirstEnd):
            combStart = nMerStart+1
            combEnd = min(nMerStart+self._nMerSpan, cdrSeqLen)
            comb = itertools.combinations(range(combStart, combEnd), self._nMerLen-1)
            for currComb in comb:
                nMerPos = [nMerStart] + list(currComb)
                nMerPosList.append(nMerPos)
                nMerStrList.append("".join([cdrSeq[i] for i in nMerPos]))
                
        nmers.strList = nMerStrList
        nmers.posList = nMerPosList
        
        return nmers

        
    def get_cdr_nmer_matrix(self, cdrSeq=None):
        
        nMers = self.get_cdr_nmers(cdrSeq)
        
        cdrSeqLen = len(cdrSeq)
        nMerNum = len(nMers.strList)

        nMerVecMat = matlib.zeros((nMerNum, self._nMerVecDim))
        for nMerI in range(nMerNum):
            nMerVecMat[nMerI,:] = self.get_nmer_vec(nMers.strList[nMerI], nMers.posList[nMerI], cdrSeqLen)
        
        return nMerVecMat


class Gmm:
    
    def __init__(self, modelFile = None):
        classFileName = inspect.getfile(self.__class__)
        classDir = os.path.dirname(classFileName) 
        
         
        self._dataDir = os.path.abspath(os.path.join(classDir, "..", "..", "data"))
        self._c_objFile = os.path.join(classDir, Gmm_model_c.gmm_obj_file_name)
        
        if modelFile == None:
            modelFile = os.path.join(self._dataDir, 'gmm_nMer03_nSpan04_nMix_600.tsv')
        if not os.path.isfile(modelFile):
            raise IOError("GMM model file file not found") 
        if not os.path.isfile(self._c_objFile):
            warnings.warn("C-code not compiled properly for GMM. Compiling a .so file will significantly speed up computations.", UserWarning)
            self.set_cFlag(False)
        else:
            self.set_cFlag(True)
        
        self._read_model(modelFile)
#         print self._model.mus
    
    def _read_model(self, modelFile):
        with open(modelFile, 'rU') as f:
            headerLine = f.next()
            fields = headerLine.strip().split("\t")
            if len(fields) != 4:
                raise IOError("Header of model file not correct")
            self.nMerLen  = int(fields[0])
            self.nMerSpan = int(fields[1])
            self.numMix = int(fields[2])
            self.dim = int(fields[3])
            
            priors      = matlib.zeros(shape=(1,self.numMix))
            gaussCoeffs = matlib.zeros(shape=(1,self.numMix))
            mus    = matlib.zeros(shape=(self.numMix, self.dim))
            sig_inv = []
            for _ in range(self.numMix):
                sig_inv.append(matlib.zeros(shape=(self.dim,self.dim)))
            
            if self.get_cFlag() == True:
                self._model = Gmm_model_c(dim=self.dim, numMix=self.numMix, c_object_file=self._c_objFile)
            else:
                self._model = Gmm_model(dim=self.dim, numMix=self.numMix)
            
            nextRow = f.next().strip().split("\t")
            priors[:] = [float(i) for i in nextRow]
            self._model.set_mixPriors(priors)
            
            nextRow = f.next().strip().split("\t")
            gaussCoeffs[:] = [float(i) for i in nextRow]
            self._model.set_gaussCoeffs(gaussCoeffs)
            
            for bI in range(self.numMix):
                nextRow = f.next().strip().split("\t")
                mu = [float(i) for i in nextRow]
                mus[bI,:] = mu
            self._model.set_mus(mus)
            
            for bI in range(self.numMix):
                for dI in range(self.dim):
                    nextRow = f.next().strip().split("\t")
                    sig_inv_row = [float(i) for i in nextRow]
                    sig_inv[bI][dI,:] = sig_inv_row
            self._model.set_sig_inv(sig_inv)

    def set_cFlag(self, cFlag): 
        if not cFlag == True:
            self._cFlag = False
        else:
            self._cFlag = True
    def get_cFlag(self):
        return self._cFlag
            
    def vec_2_posterior_c(self, vec=None):
        postVec  = matlib.zeros(shape=[1, self._model.numMix])
        return postVec

    def mat_2_nMerPosteriors(self, mat=None):
        numVecs,_ = mat.shape
        postMat = matlib.zeros(shape=(numVecs,self.numMix))
        for i in range(numVecs):
                postMat[i,:] = self._model.vec_2_posterior(mat[i,:])        
        return postMat
    
    def mat_2_cdrPosteriors(self, mat = None):
        postMatrix = self.mat_2_nMerPosteriors(mat=mat)
        numVecs,numBins = postMatrix.shape
        postVec = np.zeros(shape=(1,numBins))
        for i in range(numVecs):
            postVec += postMatrix[i,:]
        return (postVec, postMatrix)
    
    def num_mix(self):
        return self._model.numMix