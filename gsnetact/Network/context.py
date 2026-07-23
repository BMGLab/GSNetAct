"""Restricting a context-agnostic prior to the context actually measured.

The argument
------------
STRING is assembled across tissues, cell types and conditions.  An edge in it
asserts that two proteins are functionally associated *somewhere* in the
organism.  It does not assert that the interaction is available in the cells
being scored: if one partner is not transcribed in this tissue, the edge cannot
be carrying signal here, and letting it shape the weights imports connectivity
from a context that was never sampled.

This matters most for the C8 cell-type collections, which are built from other
tissues.  A hepatocyte-derived programme scored on PBMCs will still have a dense
STRING subgraph, because STRING knows those genes interact in liver - but almost
none of those genes are expressed in a monocyte, so the prior is describing a
network that is not switched on.

The remedy is standard practice for prior-guided single-cell methods: prune the
prior with the data before using it.  Regulon inference does exactly this
(Margolin et al. 2006, ARACNe, *BMC Bioinformatics* 7:S7; Aibar et al. 2017,
SCENIC, *Nat Methods* 14:1083), as does context-specific metabolic model
extraction.  Here the pruning rule is deliberately the weakest one that is
defensible - detection, not co-expression:

    an edge is retained only if **both** partners are detected in at least
    ``min_detection`` of the cells in the dataset.

Detection is used rather than correlation for two reasons.  It is a statement
about whether the gene is transcribed at all, which is what "is this interaction
available here" actually asks; and it does not read the covariance structure
that the downstream clustering is supposed to discover, so it cannot manufacture
the co-expression signal it is later credited with finding.

Nodes are never removed, only their edges - gene-set membership is a property of
the gene set, not of the dataset.  A gene that falls below the threshold simply
stops contributing network evidence and, under the default weighting, drops to
weight zero.

Empirically, on the GSNetAct benchmark, restriction at the default
``min_detection = 0.05`` improved median NMI from 0.574 to 0.593 on CITE-seq PBMC
(protein-defined identity labels, 529 C8 sets) and from 0.149 to 0.196 on Kang
IFN-beta PBMC (condition labels, 48 Hallmark sets).  Every threshold from 0.01
to 0.20 beat no restriction on both datasets, but the optimum is
dataset-dependent - CITE-seq plateaus by 0.05 while IFN-beta keeps improving to
0.20 (0.241).  The default is the conservative end of that range, discarding the
least prior information, rather than the per-dataset optimum; raise it for
shallow data or for gene sets drawn from a different tissue.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import sparse

__all__ = ["detection_rate", "detection_mask", "restrict_to_context"]


def detection_rate(matrix: Any) -> np.ndarray:
    """Fraction of cells in which each gene is detected (non-zero).

    Accepts a dense array or any SciPy sparse matrix with cells in rows and
    genes in columns.  Computed on non-zero counts, which is invariant to
    library-size normalisation and to the log transform, so the same threshold
    means the same thing whichever layer is scored.
    """
    if sparse.issparse(matrix):
        n_cells = matrix.shape[0]
        detected = np.asarray((matrix != 0).sum(axis=0)).ravel()
    else:
        dense = np.asarray(matrix)
        n_cells = dense.shape[0]
        detected = np.count_nonzero(dense, axis=0)
    if n_cells == 0:
        return np.zeros(matrix.shape[1], dtype=np.float64)
    return detected.astype(np.float64) / float(n_cells)


def detection_mask(rates: np.ndarray, min_detection: float) -> np.ndarray:
    """Boolean mask of genes detected in at least ``min_detection`` of cells."""
    return np.asarray(rates, dtype=np.float64) >= float(min_detection)


def restrict_to_context(graph, keep_genes: set[str]):
    """Confine ``graph``'s edges to genes in ``keep_genes``.

    Returns the graph unchanged when every gene passes, so the common case costs
    nothing.
    """
    mask = np.array([gene in keep_genes for gene in graph.genes], dtype=bool)
    if mask.all():
        return graph
    return graph.restrict(mask)
