import os.path
import io_evd.tsv
import data_transform.nMer
import classify.per_epitope

class ReadOptions:
    def __init__(self, \
                 cdr3b_seq_field=None, \
                 cdr3a_seq_field=None, \
                 vb_gene_field=None, \
                 db_gene_field=None, \
                 jb_gene_field=None, \
                 va_gene_field=None, \
                 ja_gene_field=None):
        self._cdr3b_seq_field = cdr3b_seq_field
        self._cdr3a_seq_field = cdr3a_seq_field
        self._vb_gene_field = vb_gene_field
        self._db_gene_field = db_gene_field
        self._jb_gene_field = jb_gene_field
        self._va_gene_field = va_gene_field
        self._ja_gene_field = ja_gene_field
          
    def set_cdr3b_seq_field(self, cdr3b_seq_field):
        self._cdr3b_seq_field = cdr3b_seq_field
    def set_cdr3a_seq_field(self, cdr3a_seq_field):
        self._cdr3a_seq_field = cdr3a_seq_field
    def set_vb_gene_field(self, vb_gene_field):
        self._vb_gene_field = vb_gene_field
    def set_db_gene_field(self, db_gene_field):
        self._db_gene_field = db_gene_field
    def set_jb_gene_field(self, jb_gene_field):
        self._jb_gene_field = jb_gene_field
    def set_va_gene_field(self, va_gene_field):
        self._va_gene_field = va_gene_field
    def set_ja_gene_field(self, ja_gene_field):
        self._ja_gene_field = ja_gene_field
        
    def get_cdr3b_seq_field(self):
        return self._cdr3b_seq_field
    def get_cdr3a_seq_field(self):
        return self._cdr3a_seq_field
    def get_vb_gene_field(self):
        return self._vb_gene_field
    def get_db_gene_field(self):
        return self._db_gene_field
    def get_jb_gene_field(self):
        return self._jb_gene_field
    def get_va_gene_field(self):
        return self._va_gene_field
    def get_ja_gene_field(self):
        return self._ja_gene_field
    
class WriteOptions:
    def __init__(self, numEpi=10, showWeights=True, numWeightsMax=25):
        self._numEpi = numEpi
        self._showWeights=showWeights
        self._numWeightsMax=numWeightsMax
        
    def set_numEpi(self, numEpi):
        self._numEpi=numEpi
    def set_aaWeights(self, showWeights):
        self._showWeights=showWeights
    def set_numWeightsMax(self, numWeightsMax):
        self._numWeightsMax=numWeightsMax

    def get_numEpi(self):
        return self._numEpi
    def get_showWeights(self):
        return self._showWeights
    def get_numWeightsMax(self):
        return self._numWeightsMax
   
class ModelParams:
    def __init__(self, aaVecFile=None, gmmFile=None, epitopeFile=None):
        self._aaVec     = aaVecFile
        self._gmm       = gmmFile
        self._epitopes  = epitopeFile
        
    def set_aaVec_file(self, aaVecFile):
        self._aaVec = aaVecFile
    def set_gmm_file(self, gmmFile):
        self._gmm = gmmFile
    def set_epitopes_file(self, epiFile):
        self._epitopes = epiFile
          
    def get_aaVec_file(self):
        return self._aaVec
    def get_gmm_file(self):
        return self._gmm
    def get_epitopes_file(self):
        return self._epitopes
         
class Stream:
    epiNamePrefix = "Epitope_"
    epiPostPrefix = "Posterior_"
    aaWeightsPrefix = "aaWeights_"
    def __init__(self, \
                 inFile=None, \
                 outFile=None, \
                 readOptions=ReadOptions(), \
                 writeOptions=WriteOptions(), \
                 modelParams=ModelParams()):
        
        self.set_inFile(inFile)
        self.set_outFile(outFile)
        self.set_readOptions(readOptions)
        self.set_writeOptions(writeOptions)
        self.set_modelParams(modelParams)
        
        self._load_model()
             
    def set_inFile(self, inFile):
        if inFile == None:
            self._inFile = None
        else:
            if os.path.isfile(inFile):
                self._inFile = inFile
            else:
                raise IOError("inFile not found")
    def set_outFile(self, outFile):
        self._outFile = outFile
    def set_readOptions(self, readOptions):
        if isinstance(readOptions, ReadOptions):
            self._readOptions = readOptions
        else:
            raise IOError("readOptions is not the proper class type")
    def set_writeOptions(self, writeOptions):
        if isinstance(writeOptions, WriteOptions):
            self._writeOptions = writeOptions
        else:
            raise IOError("writeOptions is not the proper class type")
    def set_modelParams(self, modelParams):
        if isinstance(modelParams, ModelParams):
            self._modelParams = modelParams
        else:
            raise IOError("modelParams is not the proper class type")
    
    def add(self):
        if self._model == None:
            self._load_model()
        ioStream = io_evd.tsv.Data(inFile=self._inFile, outFile=self._outFile)
        ioStream.set_in_fields(self._get_infield_list())
        
        outFieldNames = self._get_outfield_list()
        ioStream.set_out_fields(outFieldNames)
        for entry in ioStream:
            cdrSeq = entry[0]
            print cdrSeq
            (epiNames, epiPosts, aaWeights) = self._model.get_epitopes(cdrSeq)
            outEntries = self._epiModelOut_2_out_fields(epiName=epiNames, epiPost=epiPosts, aaWeights=aaWeights)
            ioStream.write_out_fields(outEntries)
        ioStream.write()    
    
    def _load_model(self):
        if self._modelParams.get_aaVec_file() == None or \
           self._modelParams.get_gmm_file() == None or \
           self._modelParams.get_epitopes_file() == None:
            self._model = None
        else:
            self._model = Model(modelParams=self._modelParams)
        
    def _get_infield_list(self):
        inFieldList = []
        inFieldList.append(self._readOptions.get_cdr3b_seq_field())
        return inFieldList
    
    def _get_outfield_list(self):
        outFields = []
        for i in range(self._writeOptions.get_numEpi()):
            outFields.append("%s%02d"%(Stream.epiNamePrefix, i+1))
        for i in range(self._writeOptions.get_numEpi()):
            outFields.append("%s%02d"%(Stream.epiPostPrefix, i+1))
        if self._writeOptions.get_showWeights():
            for i in range(self._writeOptions.get_numWeightsMax()):
                outFields.append("%s%02d"%(Stream.aaWeightsPrefix, i+1))
        return outFields
    
    def _epiModelOut_2_out_fields(self, epiName=None, epiPost=None, aaWeights=None):
        outFields = []
        for i in range(self._writeOptions.get_numEpi()):
            outFields.append(epiName[i])
        for i in range(self._writeOptions.get_numEpi()):
            outFields.append("%.4f"%epiPost[i])
        if self._writeOptions.get_showWeights():
            numAAshow = min(self._writeOptions.get_numWeightsMax(), len(aaWeights))
            numAApad = self._writeOptions.get_numWeightsMax() - numAAshow
            for i in range(numAAshow):
                outFields.append("%.4f"%aaWeights[i])
            for i in range(numAApad):
                outFields.append("%.4f"%0)
        return outFields        
        
class Model:    
    def __init__(self, \
                 modelParams=ModelParams()):
        self.set_symVecFile(modelParams.get_aaVec_file())
        self.set_gmmModelFile(modelParams.get_gmm_file())
        self.set_epiModelFile(modelParams.get_epitopes_file())
        self._models_loaded = False
          
    def _load_models(self):
        if self._symVecFile == None or self._gmmModelFile == None or self._epiModelFile == None:
            raise IOError("Not all model files specified properly")
        
        self._gmmModel = data_transform.nMer.Gmm(modelFile=self._gmmModelFile)
        self._gmmModel.set_cFlag(False)
        self._nMerLen = self._gmmModel.nMerLen
        self._nMerSpan = self._gmmModel.nMerSpan
        
        self._aaVecSpace = data_transform.nMer.AaVecSpace( \
                    symVecFile=self._symVecFile, \
                    nMerLen=self._nMerLen, \
                    nMerSpan=self._nMerSpan)
        self._epiModel = classify.per_epitope.Model(self._epiModelFile, self._gmmModel.num_mix())
        self._models_loaded = True
    
    def set_symVecFile(self, symVecFile=None):
        if symVecFile == None:
            self._symVecFile = None
        else:
            if os.path.isfile(symVecFile):
                self._symVecFile = symVecFile
            else:
                self._symVecFile = None
                raise IOError("symVecFile not found.")
                  
    def set_gmmModelFile(self, gmmModelFile=None):
        if gmmModelFile == None:
            self._gmmModelFile = None
        else:
            if os.path.isfile(gmmModelFile):
                self._gmmModelFile = gmmModelFile
            else:
                self._gmmModelFile = None
                raise IOError("gmmModelFile not found.")
            
    def set_epiModelFile(self, epiModelFile=None):
        if epiModelFile == None:
            self._epiModelFile = None
        else:
            if os.path.isfile(epiModelFile):
                self._epiModelFile = epiModelFile
            else:
                self._epiModelFile = None
                raise IOError("epiModelFile not found.")
            
    def get_epitopes(self, cdrSeq):
        if not self._models_loaded:
            self._load_models()
        nMerMat = self._aaVecSpace.get_cdr_nmer_matrix(cdrSeq)
        (postVec, postMat) = self._gmmModel.mat_2_cdrPosteriors(nMerMat)
        (epiNames, epiPosts) = self._epiModel.get_posteriors(postVec)
        aaWeights = self._epiModel.get_aa_weights(cdrSeq=cdrSeq, nMerMat=nMerMat, postMat=postMat, gmmModel=self._gmmModel, epiName=epiNames[0])
        
        return (epiNames, epiPosts, aaWeights)
        