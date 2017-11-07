import inspect
import os
import itertools
import numpy.matlib as matlib
import numpy as np

class Gmm_model:
    def __init__(self, dim=None, numMix=None):
        if not isinstance(dim, (int, long)) or not isinstance(numMix, (int, long)):
            raise IOError("Both dim and numMix must be specified")
        self.dim = dim
        self.numMix = numMix
        self.mixPriors = matlib.zeros(shape=(1, numMix))
        self.gaussCoeffs = matlib.zeros(shape=(1, numMix))
        self.mus = matlib.zeros(shape=(numMix, dim))
        self.sig_inv = []
        for _ in range(numMix):
            self.sig_inv.append(matlib.zeros(shape=(dim,dim)))
        
    def set_mixPriors(self, mixPriors=None):
        if len(mixPriors) != self.numMix:
            raise IndexError("Length of mixPriors doesn't match the number of mixtures") 
        self.mixPriors[:] = mixPriors
            
    def set_gaussCoeffs(self, gaussCoeffs=None):
        if len(gaussCoeffs) != self.numMix:
            raise IndexError("Length of gaussCoeffs doesn't match the number of mixtures") 
        self.gaussCoeffs[:] = gaussCoeffs
    
    def set_mu(self, mixI=None, mu=None):
        if len(mu) != self.dim:
            raise IndexError("Mu vector does not have the proper dimension")
        if mixI < 0 or mixI >= self.numMix:
            raise IndexError("mixI out of bound")
        self.mus[mixI,:] = mu
        
    def set_sig_inv_row(self, mixI=None, rowI=None, sig_inv_row=None):
        if len(sig_inv_row) != self.dim:
            raise IndexError("Covariance matrix row does not have the proper dimension")
        if mixI < 0 or mixI >= self.numMix:
            raise IndexError("mixI out of bound")
        if rowI < 0 or rowI >= self.dim:
            raise IndexError("Row index out of bound")
        self.sig_inv[mixI][rowI,:] = sig_inv_row

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
        
    def get_cdr_nmer_matrix(self, cdrSeq=None):
        
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
                
        nMerNum = len(nMerStrList)
#         print nMerNum, self._nMerVecDim

        nMerVecMat = matlib.zeros((nMerNum, self._nMerVecDim))
        for nMerI in range(nMerNum):
            nMerVecMat[nMerI,:] = self.get_nmer_vec(nMerStrList[nMerI], nMerPosList[nMerI], cdrSeqLen)
        
        return nMerVecMat


class Gmm:
    def __init__(self, modelFile = None, cFlag=None):
        classFileName = inspect.getfile(self.__class__)
        classDir = os.path.dirname(classFileName)
        
        
        self.set_cFlag(cFlag)
         
        self._dataDir = os.path.abspath(os.path.join(classDir, "..", "..", "data"))
        if modelFile == None:
            modelFile = os.path.join(self._dataDir, 'gmm_nMer03_nSpan04_nMix_600.tsv')
        if not os.path.isfile(modelFile):
            raise IOError("GMM model file file not found") 
        
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
            
            self._model = Gmm_model(dim=self.dim, numMix=self.numMix)
            
            nextRow = f.next().strip().split("\t")
            priors = [float(i) for i in nextRow]
            self._model.set_mixPriors(priors)
            nextRow = f.next().strip().split("\t")
            gaussCoeffs = [float(i) for i in nextRow]
            self._model.set_gaussCoeffs(gaussCoeffs)
            
            for bI in range(self.numMix):
                nextRow = f.next().strip().split("\t")
                mu = [float(i) for i in nextRow]
                self._model.set_mu(bI, mu)
            for bI in range(self.numMix):
                for dI in range(self.dim):
                    nextRow = f.next().strip().split("\t")
                    sig_inv_row = [float(i) for i in nextRow]
                    self._model.set_sig_inv_row(bI, dI, sig_inv_row)

    def set_cFlag(self, cFlag): 
        if not cFlag == True:
            self._cFlag = False
        else:
            self._cFlag = True
    def get_cFlag(self):
        return self._cFlag
            

    def vec_2_posterior(self, vec=None):
        _, vecLen = vec.shape
        if vecLen != self.dim:
            raise IndexError("nMer vector dimension doesn't agree with gmm dimension")
#         postVec = matlib.zeros(shape=[1, self._model.numMix])
        gaussP  = matlib.zeros(shape=[1, self._model.numMix])
        for bI in range(self._model.numMix):
            x = vec - self._model.mus[bI,:]
            x_sqr = x*self._model.sig_inv[bI]*x.transpose()
            pwr = -0.5*x_sqr
            gaussP[0,bI] = self._model.gaussCoeffs[0,bI] * np.exp(pwr)
        
        postNorm = 1/(gaussP*self._model.mixPriors.transpose())
        
        postVec = postNorm[0,0]*matlib.multiply(gaussP, self._model.mixPriors) 
#         print postVec
        return postVec

    def vec_2_posterior_c(self, vec=None):
        postVec  = matlib.zeros(shape=[1, self._model.numMix])
        return postVec

    def mat_2_nMerPosteriors(self, mat=None):
        numVecs,_ = mat.shape
        postMat = matlib.zeros(shape=(numVecs,self.numMix))
        for i in range(numVecs):
            if self._cFlag:
                postMat[i,:] = self.vec_2_posterior_c(mat[i,:])
            else:
                postMat[i,:] = self.vec_2_posterior(mat[i,:])        
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