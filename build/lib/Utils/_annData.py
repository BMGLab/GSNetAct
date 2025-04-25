from Scoring.GeneSets.geneSetObjects import createGeneSets
from Scoring.GeneSets.geneSetScores import GeneSetScore

from Scoring.GeneExpressions.geneExpScores import score

import numpy as np
import scanpy as sc


def createObject(adataPath, jsonPath, normalized=False):
    # TODO : Print normalization process info.

    adata = sc.read_h5ad(adataPath)
    scoresArray = []
    geneSetNamesArray = []

    geneSetList = createGeneSets(jsonPath)
    
    for geneset in geneSetList:

        geneSetNamesArray.append(geneset.getID)

        _geneNames = geneset.getGeneNames

        newGeneSetScore = GeneSetScore(geneset.matrix, _geneNames)

        newExpScore = score(adata, newGeneSetScore)

        scoresArray.append(newExpScore)

    scoresArray = np.array(scoresArray).T

    adataScores = sc.AnnData(X=scoresArray, var=geneSetNamesArray, 
                             obs=adata.obs)

    if normalized:

        from sklearn.preprocessing import quantile_transform, StandardScaler

        adataScores.X = quantile_transform(adataScores.X, axis=1, 
                                           output_distribution="normal")

        adataScores.X = StandardScaler().fit_transform(adataScores.X)

    return adataScores


