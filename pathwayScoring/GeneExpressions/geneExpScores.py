import scanpy as sc
import numpy as np

from GeneSets.geneSetScores import GeneSetScore 
from GeneSets.geneSetObjects import GeneSet

from jsonParser import pjson

def score(adata, geneSetScore):
    
    geneset = np.array([geneSetScore.get(var, 0.0) for var in adata.var_names]) # Create a 1D array for the given Gene Set 
                                                                                # with respect to gene names and their indexes
                                                                                # from annData object. If the gene name from adata.var_names
                                                                                # is absent in the geneset, put 0 as the value of the index.
    
    return adata.X.dot(geneset) # Return the dot product of gene names and our array. Since the indexes of both arrays point to the same gene names, 
                                # we can simply return the dot product.

if __name__ == "__main__" : # test

    newjson = pjson("/home/sadigungor/pathwayScoring/big_genesets_relations.json")

    adata = sc.read_h5ad("/home/sadigungor/Desktop/pathwayScoring/test/test_data/pbmc3k.h5ad")

    for i in newjson:

        newGeneSet = GeneSet(f"{i}",newjson[i])

        print(newGeneSet.getID)

        _geneNames = newGeneSet.getGeneNames

        newGeneSetScore = GeneSetScore(newGeneSet.getMatrix,_geneNames)


        print(score(adata,newGeneSetScore))
