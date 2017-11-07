import inspect
import os
# import itertools
import numpy.matlib as matlib
import numpy as np
import io_evd.tsv
# import csv
# import math
# import matplotlib.pyplot as plt

class EpiModel:
    def __init__(self, epiNames=None, epiWeights=None):
        self.epiNames = epiNames
        self.epiWeights = epiWeights
        

class Model:
    def __init__(self, modelsFile, numMix):
        classFileName = inspect.getfile(self.__class__)
        classDir = os.path.dirname(classFileName)
        self._dataDir = os.path.abspath(os.path.join(classDir, "..", "..", "data"))
        if modelsFile == None:
            modelsFile = os.path.join(self._dataDir, 'epiModels_30_2017_08_30.tsv')
        if not os.path.isfile(modelsFile):
            raise IOError("Epitope modelsfile not found")
        self._read_models(modelsFile, numMix)
        
        
    def _read_models(self, modelsFile, numMix):
        ioStream = io_evd.tsv.Data(inFile=modelsFile,outNAsOnly=False)
        ioStream.set_dialect('excel-tab')
        beta_fields = ['beta_%03d'%(i+1) for i in range(numMix)]
        inFields = ["model_name", "bias"] + beta_fields
        
        ioStream.set_in_fields(inFields)
        epiNames = []
        epiWeights = [];
        for modelInfo in ioStream:    
            epiNames.append(modelInfo[0])
            weights = []
            weights = [modelInfo[1]] + modelInfo[2:]
            weightsInt = [float(weight) for weight in weights]
            epiWeights.append(weightsInt)
        epiWeightMatrix = matlib.matrix(epiWeights).transpose()
        self._models = EpiModel(epiNames, epiWeightMatrix)
        
    def get_posteriors(self, feat_vec=None):
        scoresM = np.insert(feat_vec, 0, 1) * self._models.epiWeights
        scores = list(np.array(scoresM).reshape(-1,)) 
        sortI = sorted(range(len(scores)), key=lambda k: scores[k])
        epiNames = [self._models.epiNames[i] for i in reversed(sortI)]
        epiScores =  np.array([scores[i] for i in reversed(sortI)])
#         print epiScores
        epiPost = 1/(1 + np.exp(-epiScores))
        
        return (epiNames, epiPost)
    
    def get_aa_weights(self, cdrSeq=None, nMerMat=None, postMat= None, gmmModel=None, epiName=None):
        epiI = self._models.epiNames.index(epiName)

        pos_matrix = nMerMat[:,-gmmModel.nMerLen:]
        nRow,nCol = pos_matrix.shape
        pos_vec = pos_matrix.reshape((1, nRow*nCol))
        uniPos = set([pos_vec[0,i] for i in range(nRow*nCol)])
        cdrLen = len(uniPos)
        pos_matrix = np.round((cdrLen-1)*pos_matrix)
#         print pos_matrix
        if postMat == None:
            postMat = gmmModel.mat_2_nMerPosteriors(nMerMat)
        
        (_,numMix) = postMat.shape
        aaWeights = np.zeros(shape=cdrLen)
        for aaI in range(cdrLen):
            aaPost = matlib.zeros(shape=(1,numMix))
            occurI = [aaI in row for row in pos_matrix]
            for nMerI in range(nRow):
                if occurI[nMerI]:
                    aaPost += postMat[nMerI,:]
            aaWeight = aaPost * self._models.epiWeights[1:,epiI]
            aaWeight /= gmmModel.nMerLen 
            aaWeight += self._models.epiWeights[0,epiI]/cdrLen
            aaWeights[aaI] = aaWeight
        
#             print aaWeight
#             plt.figure(aaI)
#             plt.plot(np.asarray(aaPost)[0])
#             plt.show()
            
        
        return aaWeights
            
            