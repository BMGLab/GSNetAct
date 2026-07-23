"""GSNetAct: network-weighted gene-set activity for single-cell data.

Quick start
-----------
>>> from gsnetact import pjson, runGSNA
>>> activity = runGSNA(adata, pjson("geneSets.json"), layer="lognorm")

The default weighting changed in 0.1.0.  ``weights="legacy"`` reproduces
pre-0.1.0 numbers; see ``CHANGELOG.md`` and
:mod:`gsnetact.Network.weights` for the reasoning and the benchmark evidence.
"""

__version__ = "0.1.0"

from .GeneExpressions.geneExpScores import buildWeightMatrix, expScore, scoreMatrix
from .GeneSets.geneSetObjects import GeneSet, GeneSetMatrix, getGSNA
from .GeneSets.geneSetScores import GeneSetScore, NetworkGeneWeights
from .Network import (
    GeneSetGraph,
    available_weightings,
    detection_rate,
    diffusion_weights,
    gene_weights,
    legacy_degree_weights,
    pagerank_weights,
    register_weighting,
    uniform_weights,
    weight_concentration,
)
from .Utils._annData import normalizeScores, runGSNA
from .Utils.jsonParser import pjson
from .Utils.makeJsonFile import makeJson
from .Utils.MatrixItem import MatrixItem

__all__ = [
    "__version__",
    # scoring
    "runGSNA",
    "expScore",
    "scoreMatrix",
    "buildWeightMatrix",
    "normalizeScores",
    # gene sets
    "GeneSet",
    "GeneSetMatrix",
    "GeneSetScore",
    "NetworkGeneWeights",
    "getGSNA",
    # network
    "GeneSetGraph",
    "gene_weights",
    "diffusion_weights",
    "legacy_degree_weights",
    "uniform_weights",
    "pagerank_weights",
    "available_weightings",
    "register_weighting",
    "weight_concentration",
    "detection_rate",
    # io
    "pjson",
    "makeJson",
    "MatrixItem",
]
