from ..GeneSets.geneSetObjects import getGSNA
from ..GeneSets.geneSetScores import GeneSetScore

from ..GeneExpressions.geneExpScores import expScore

import numpy as np

from anndata import AnnData


def runGSNA(adata, jsonFile, normalized=False):
    # TODO : Print normalization process info.

    scoresArray = []
    geneSetNamesArray = []
    # Create arrays to store scores and names of gene sets

    geneSetList = getGSNA(jsonFile)
    # Create GeneSet objects from the JSON file.

    for geneset in geneSetList:
        # Calculate scores for each gene set.

        geneSetNamesArray.append(geneset.getID)

        _geneNames = geneset.getGeneNames

        newGeneSetScore = GeneSetScore(geneset.matrix, _geneNames)

        newExpScore = expScore(adata, newGeneSetScore)

        scoresArray.append(newExpScore)

    scoresArray = np.array(scoresArray).T

    adataScores = AnnData(X=scoresArray, var=geneSetNamesArray,
                          obs=adata.obs)

    if normalized:
        # If the normalized option is on, normalize the score data with
        # Quantile normalization and Z-Score normalization.

        from sklearn.preprocessing import quantile_transform, StandardScaler

        adataScores.X = quantile_transform(adataScores.X, axis=1,
                                           output_distribution="normal")

        adataScores.X = StandardScaler().fit_transform(adataScores.X)

    return adataScores
