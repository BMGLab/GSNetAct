from ..GeneSets.geneSetObjects import getGSNA
from ..GeneSets.geneSetScores import GeneSetScore

from ..GeneExpressions.geneExpScores import score

import numpy as np
import scanpy as sc


def createAdataObject(adataPath, jsonPath, normalized=False):
    # TODO : Print normalization process info.

    adata = sc.read_h5ad(adataPath)

    scoresArray = []
    geneSetNamesArray = []
    # Create arrays to store scores and names of gene sets

    geneSetList = getGSNA(jsonPath)
    # Create GeneSet objects from the JSON file.

    for geneset in geneSetList:
        # Calculate scores for each gene set.

        geneSetNamesArray.append(geneset.getID)

        _geneNames = geneset.getGeneNames

        newGeneSetScore = GeneSetScore(geneset.matrix, _geneNames)

        newExpScore = score(adata, newGeneSetScore)

        scoresArray.append(newExpScore)

    scoresArray = np.array(scoresArray).T

    adataScores = sc.AnnData(X=scoresArray, var=geneSetNamesArray,
                             obs=adata.obs)

    if normalized:
        # If the normalized option is on, normalize the score data with
        # Quantile normalization and Z-Score normalization.

        from sklearn.preprocessing import quantile_transform, StandardScaler

        adataScores.X = quantile_transform(adataScores.X, axis=1,
                                           output_distribution="normal")

        adataScores.X = StandardScaler().fit_transform(adataScores.X)

    return adataScores
