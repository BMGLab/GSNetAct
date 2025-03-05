from GeneSets.geneSetObjects import GeneSet
from GeneSets.geneSetScores import GeneSetScore

from GeneExpressions.geneExpScores import score

from jsonParser import pjson

import scanpy as sc
import numpy as np


jsonFile = pjson("/home/sadigungor/Desktop/pathway_scorers/test/test_data/big_genesets_relations.json")

adata = sc.read_h5ad("/home/sadigungor/Desktop/pathway_scorers/test/test_data/pbmc3k.h5ad") 

for i in jsonFile : 

    newGeneSet = GeneSet(f"{i}",jsonFile[i])

    print(newGeneSet.getID)

    _geneNames = newGeneSet.getGeneNames

    newGeneSetScore = GeneSetScore(newGeneSet.getMatrix, _geneNames)

    print(score(adata,newGeneSetScore))
