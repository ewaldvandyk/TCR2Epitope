import os
import TCR2epitopePipe

# Setup paths
script_path = os.path.dirname(os.path.realpath(__file__))
base_path = os.path.abspath(os.path.join(script_path, '..'))
data_path = os.path.join(base_path, "data")
test_path = os.path.join(base_path, 'test')


# Setup input and output file. Output fields are appended as columns to the output file
# Input must be in comma/tab separated text format. Output format is the same   
inFile  = os.path.join(test_path, "rand_cdr3_in.tsv")
outFile = os.path.join(test_path, "rand_cdr3_gmmPost_out.tsv")


# Set input options, i.e. specify the relevant input field names as specified in the header 
readOptions = TCR2epitopePipe.ReadOptions()
readOptions.set_cdr3b_seq_field("CDR3 sequences")

# Set output options (used by TCR2epitope)
writeOptions = TCR2epitopePipe.WriteOptions()
writeOptions.set_numEpi(10)         # The number of top scoring epitopes to report
writeOptions.set_aaWeights(True)    # Flag indicating if amino acid weights should
                                    # be appended to the end of the file. Weights
                                    # are only reported for the top scoring epitope
writeOptions.set_numWeightsMax(25)  # The number of cdr3 amino acids for which 
                                    # weights are reported. aas beyond this number
                                    # are truncated. If there are fewer aas, 
                                    # zeros are padded 
                                    
                                    
# Setup model files
symVecFile      = os.path.join(data_path, "pmbecAAvecs.tsv")                    # Text file describing how amino acids are embedded in a space
gmmModelFile    = os.path.join(data_path, "gmm_nMer03_nSpan04_nMix_600.tsv")    # Gaussian mixture model that clusers similar n-mers together
epiModelFile    = os.path.join(data_path, "epiModels_30_2017_08_30.tsv")        # Text file with specificity linear models for multiple epitopes

# Load model file names into model parameter class (used by TCR2epitope)
modelParams = TCR2epitopePipe.ModelParams()
modelParams.set_aaVec_file(symVecFile) 
modelParams.set_gmm_file(gmmModelFile)
modelParams.set_epitopes_file(epiModelFile)


epiStream = TCR2epitopePipe.Stream(inFile = inFile, \
                                   outFile = outFile, \
                                   readOptions=readOptions, \
                                   writeOptions=writeOptions, \
                                   modelParams=modelParams)
epiStream.add()


