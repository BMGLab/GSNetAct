from gsnetact import getGSNA, GeneSetScore, score

import scanpy as sc
import numpy as np


#Get the files.
jsonFilePath = "/home/sadigungor/Desktop/GSNetAct/test/test_data/deneme.json"

h5adFilePath = "/home/sadigungor/Desktop/GSNetAct/test/test_data/pbmc3k.h5ad"
adata = sc.read_h5ad(h5adFilePath)
geneSetList = getGSNA(jsonFilePath)

arr = []
for geneSet in geneSetList:
    # Calculate for each gene set.

    _geneNames = geneSet.getGeneNames
    # Get the gene names necessary for geneExpScores.score() from
    # GeneSet object.

    newGeneSetScore = GeneSetScore(geneSet.matrix, _geneNames)
    print(newGeneSetScore)
    # Create gene set score without
    # the expressions.
    
    newScore = score(adata, newGeneSetScore) 
    # Merge expression scores with gene set scores and print

"""
    # Print geneset names, scores and the percentage of scores that are not 0 in the given geneset.
    print(f"GENESET : {geneSet.getID}")
    print(f"{newScore} \n NONZEROS : %{round(((np.count_nonzero(newScore)/len(newScore)) * 100),2)}")
    arr.append(newScore)

  
    print("######################################################################################################")
"""

arr = np.array(arr)

np.savetxt("output.txt",arr)
