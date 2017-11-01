import os
import io_evd.tsv
import data_transform.nMer
import classify.per_epitope
import numpy as np

# User parameters
cFlag = False # Set cFlag to true to use Renske's c module nMer.vec_2_posterior_c


# Setup paths
script_path = os.path.dirname(os.path.realpath(__file__))
base_path = os.path.abspath(os.path.join(script_path, '..'))
data_path = os.path.join(base_path, "data")
test_path = os.path.join(base_path, 'test')

# Setup model files
symVecFile      = os.path.join(data_path, "pmbecAAvecs.tsv")
gmmModelFile    = os.path.join(data_path, "gmm_nMer03_nSpan04_nMix_600.tsv")
epiModelFile    = os.path.join(data_path, "epiModels_30_2017_08_30.tsv")

# Setup input and output file
inFile  = os.path.join(test_path, "rand_cdr3_in.csv")
if cFlag:
    outFile = os.path.join(test_path, "rand_cdr3_gmmPost_out_c.csv")
else:
    outFile = os.path.join(test_path, "rand_cdr3_gmmPost_out.csv")

# Initialize models
aaVecSpace = data_transform.nMer.AaVecSpace( \
            symVecFile=symVecFile, \
            nMerLen=3, \
            nMerSpan=4)
gmmModel = data_transform.nMer.Gmm( \
            modelFile=gmmModelFile)
gmmModel.set_cFlag(cFlag)
numMix = gmmModel.num_mix()

epiModel = classify.per_epitope.Model(epiModelFile, numMix)
#Setup IO
ioStream = io_evd.tsv.Data(inFile=inFile, outFile=outFile)
ioStream.set_dialect('excel')
ioStream.set_in_fields("CDR3 sequences")

outFields = ["posterior_%03d"%mixI for mixI in range(numMix)]
ioStream.set_out_fields(outFields)

for entry in ioStream:
    cdrSeq = entry[0]
    print cdrSeq
    nMerMat = aaVecSpace.get_cdr_nmer_matrix(cdrSeq)
    postVec = gmmModel.mat_2_cdrPosteriors(nMerMat)
    (epiNames, epiPosts) = epiModel.get_posteriors(postVec)
    aaWeights = epiModel.get_aa_weights(cdrSeq=cdrSeq, nMerMat=nMerMat, gmmModel=gmmModel, epiName=epiNames[0])
    print len(cdrSeq)
    print aaWeights.shape

    postVecStr = ["%0.5e"%postVec[0,i] for i in range(numMix)]
#     print postVecStr
    ioStream.write_out_fields(postVecStr)
ioStream.write()

