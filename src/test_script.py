import os
import io_evd.tsv
import data_transform.nMer

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
    postVecStr = ["%0.5e"%postVec[0,i] for i in range(numMix)]
#     print postVecStr
    ioStream.write_out_fields(postVecStr)
ioStream.write()

