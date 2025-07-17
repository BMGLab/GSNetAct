import numpy as np


def expScore(adata, geneSetScore):
    
    #TODO: adata.var_names kontrol et ve ENSEMBL id si varsa uyari mesaji dondur. 

    geneset = np.array([geneSetScore.get(var, 0.0) for var in adata.var_names])
    # Create a 1D array for the given Gene Set
    # with respect to gene names and their indexes
    # from annData object. If the gene name from adata.var_names
    # is absent in the geneset, put 0 as the value of the index.
    
    adata.X += 10**-6
    # To ensure that all expression scores have a non-zero contribution 
    # to the overall computation,add a negligible constant 
    # to all entries in adata.X.

    return adata.X.dot(geneset)
    # Return the dot product of gene names and our array. Since the indexes of
    # both arrays point to the same gene names,
    # we can simply return the dot product.
