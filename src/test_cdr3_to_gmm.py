import os
import data_transform.nMer

# Setup paths
script_path = os.path.dirname(os.path.realpath(__file__))
base_path = os.path.abspath(os.path.join(script_path, '..'))
data_path = os.path.join(base_path, "data")

# Setup model files
symVecFile      = os.path.join(data_path, "pmbecAAvecs.tsv")                    # Text file describing how amino acids are embedded in a space
gmmModelFile    = os.path.join(data_path, "gmm_nMer03_nSpan04_nMix_600.tsv")    # Gaussian mixture model that clusers similar n-mers together
epiModelFile    = os.path.join(data_path, "epiModels_30_2017_08_30.tsv")        # Text file with specificity linear models for multiple epitopes

# Load models (Only load the models once)
gmmModel = data_transform.nMer.Gmm(modelFile=gmmModelFile)
aaVecSpace = data_transform.nMer.AaVecSpace( \
                    symVecFile=symVecFile, \
                    nMerLen=gmmModel.nMerLen, \
                    nMerSpan=gmmModel.nMerSpan)


########### PART 1: Unpacking n-mers from CDR3 sequence and mapping them to bag of words ##############
#Unpacking cdr3 sequences into their nmers
cdr3_seq = "CASSLLVGTSGRWGGSYEQYF"                      #Input cdr3 beta chain sequence
print cdr3_seq
nmers = aaVecSpace.get_cdr_nmers(cdr3_seq)  # Unpack all combinations of n-mers
num_nmers = len(nmers.strList)
print "Number of nmers = %d"%num_nmers
print nmers.strList                         #aa string sequence of nMers
print nmers.posList                         #aa position sequence of nMers (start at index 0)

# Get bag-of words representation of nMers
nMerMat = aaVecSpace.get_cdr_nmer_matrix(cdr3_seq)      #Unpack nMers and represent them in aa vector space
(_, gmmMatrix) = gmmModel.mat_2_cdrPosteriors(nMerMat)  #Get bag counts for each nMer 
                                                        #Each row represents a unique n-mer. Sorted the same as "nmer" above
                                                        #Each column represent a bag (of similar n-mer)
                                                        #Intuitively, matrix entry is "1.0" if n-mer belongs to the bag. "0.0" otherwise
                                                        #In reality, the clustering is soft, so the matrix entries are weights representing the degree of membership                                


########### PART 2: Mapping entire CDR3 sequences to frequency of n-mer bags ##############
# Illustration on how to convert a cdr3 sequence (beta chain) to its gmm vector representations
cdr3_seq = "CASSLLVGTSGRWGGSYEQYF"                      #Input cdr3 beta chain sequence
nMerMat = aaVecSpace.get_cdr_nmer_matrix(cdr3_seq)      #Unpack nMers and represent them in aa vector space
(postVec, _) = gmmModel.mat_2_cdrPosteriors(nMerMat)    #postVec is the GMM representation of the cdr3seq
# Another example
cdr3_seq = "CASSHKQDSPLHF"                              #Input cdr3 beta chain sequence
nMerMat = aaVecSpace.get_cdr_nmer_matrix(cdr3_seq)      #Unpack nMers and represent them in aa vector space
(postVec, _) = gmmModel.mat_2_cdrPosteriors(nMerMat)    #postVec is the GMM representation of the cdr3seq





