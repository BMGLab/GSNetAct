from GeneSets.geneSetObjects import GeneSet
from GeneSets.geneSetScores import GeneSetScore

from GeneExpressions.geneExpScores import score

from jsonParser import pjson

import scanpy as sc
import numpy as np


#Get the files.
jsonFilePath = "/home/sadigungor/Desktop/pathway_scorers/test/test_data/big_genesets_relations.json"
jsonFile = pjson(jsonFilePath)

h5adFilePath = "/home/sadigungor/Desktop/pathway_scorers/test/test_data/pbmc3k.h5ad"
adata = sc.read_h5ad(h5adFilePath)


for i in jsonFile : # Calculate for each gene set. 

    newGeneSet = GeneSet(f"{i}",jsonFile[i]) # Create GeneSet Object

    print(newGeneSet.getID)
 
    _geneNames = newGeneSet.getGeneNames # Get the gene names necessary for geneExpScores.score() from 
                                         # GeneSet object.

    newGeneSetScore = GeneSetScore(newGeneSet.getMatrix, _geneNames) # Create gene set score without 
                                                                     # the expressions.

    print(score(adata,newGeneSetScore)) # Merge expression scores with gene set scores and print
