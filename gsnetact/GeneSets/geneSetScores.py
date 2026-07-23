"""Per-gene weights for one gene set.

:class:`GeneSetScore` keeps the historical ``{gene: weight}`` dict interface and
the historical formula, but is now vectorised and crash-free.  The old
implementation looped over incidence-matrix columns in Python and unpacked
``i, j = np.nonzero(column)``, which raised ``ValueError`` on any column that did
not have exactly two non-zeros - the case for self-loops and for hash-collided
edges.  The same quantity is a single matrix-vector product:

    ``w = A @ degree``

:class:`NetworkGeneWeights` is the new-style equivalent and is what
:func:`gsnetact.runGSNA` uses; it accepts any registered weighting scheme and
applies context restriction.  See :mod:`gsnetact.Network.weights` for why the
default is no longer the degree-based formula.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np

from ..Network.graph import GeneSetGraph
from ..Network.weights import LEGACY_EPSILON, gene_weights

__all__ = ["GeneSetScore", "NetworkGeneWeights"]


class GeneSetScore(dict):
    """Legacy per-gene weights: ``w_i = sum_j A_ij * degree(j)``.

    Parameters
    ----------
    geneID:
        Gene-set identifier (kept for signature compatibility).
    matrix:
        Incidence matrix from :class:`~gsnetact.GeneSets.geneSetObjects.GeneSetMatrix`,
        or ``None`` when a :class:`~gsnetact.Network.graph.GeneSetGraph` is
        passed as ``graph``.
    geneNamesList:
        Gene names in incidence-matrix row order.
    epsilon:
        Weight given to genes with no within-set partner.

    Notes
    -----
    Retained so published GSNetAct results stay reproducible.  The formula is
    dominated by network degree - see :mod:`gsnetact.Network.weights` - and
    :class:`NetworkGeneWeights` should be preferred for new work.
    """

    def __init__(
        self,
        geneID: str,
        matrix: np.ndarray | None,
        geneNamesList: Sequence[str],
        epsilon: float = LEGACY_EPSILON,
        graph: GeneSetGraph | None = None,
    ):
        super().__init__()
        self.geneID = str(geneID)
        self.epsilon = float(epsilon)
        self.matrix = matrix
        self.geneNamesList = [str(name) for name in geneNamesList]

        if graph is not None:
            adjacency = graph.adjacency
            degree = graph.degree
        else:
            if matrix is None:
                for gene in self.geneNamesList:
                    self[gene] = 0.0
                return
            dense = np.asarray(matrix, dtype=np.float64)
            if dense.shape[0] != len(self.geneNamesList):
                raise ValueError(
                    f"matrix has {dense.shape[0]} rows but "
                    f"{len(self.geneNamesList)} gene names were given"
                )
            # Rebuild the adjacency from the incidence matrix: each column holds
            # one edge, carrying its weight in both endpoint rows.
            size = dense.shape[0]
            adjacency = np.zeros((size, size), dtype=np.float64)
            for column in range(dense.shape[1]):
                endpoints = np.nonzero(dense[:, column])[0]
                if endpoints.size != 2:
                    # Degenerate column (self-loop or merged edges). Skip it
                    # rather than crash; the old code raised ValueError here.
                    continue
                left, right = endpoints
                weight = float(dense[left, column])
                adjacency[left, right] = weight
                adjacency[right, left] = weight
            degree = (adjacency > 0).sum(axis=1).astype(np.float64)

        weights = np.asarray(adjacency @ degree, dtype=np.float64).ravel()
        weights[np.asarray(degree).ravel() == 0] = self.epsilon
        for gene, weight in zip(self.geneNamesList, weights):
            self[gene] = float(weight)


class NetworkGeneWeights(dict):
    """Per-gene weights under any registered weighting scheme.

    Parameters
    ----------
    graph:
        The gene set's interaction subgraph.
    method:
        Weighting scheme name; see
        :func:`gsnetact.Network.weights.available_weightings`.
    normalize:
        Scale weights to sum to one, so the score is a weighted mean of member
        expression and gene-set size drops out of the score's scale.
    options:
        Extra keyword arguments for the weighting function.

    Attributes
    ----------
    stats:
        Per-set diagnostics: member count, edge count, isolated-gene count and
        weight share, weight Gini, largest-component fraction, and whether the
        uniform fallback was used.  These are what :func:`gsnetact.runGSNA` puts
        into ``var``, so every set's network support can be checked instead of
        assumed.
    """

    def __init__(
        self,
        graph: GeneSetGraph,
        method: str = "diffusion",
        normalize: bool = True,
        options: Mapping[str, Any] | None = None,
    ):
        super().__init__()
        from ..Network.weights import weight_concentration

        self.graph = graph
        self.method = str(method)
        weights = gene_weights(graph, method=method, normalize=normalize, options=options)

        degree = graph.degree
        isolated = degree == 0
        total = float(weights.sum())
        self.stats = {
            "n_genes": len(graph.genes),
            "n_edges": graph.n_edges,
            "n_connected_genes": int((~isolated).sum()),
            "n_isolated_genes": int(isolated.sum()),
            "isolated_weight_fraction": float(weights[isolated].sum() / total) if total > 0 else 0.0,
            "network_density": graph.density,
            "largest_component_fraction": graph.largest_component_fraction(),
            "weight_gini": weight_concentration(weights),
            "uniform_fallback": bool(graph.n_edges == 0),
            "weighting": self.method,
        }
        for gene, weight in zip(graph.genes, weights):
            self[gene] = float(weight)

    def top(self, n: int = 10) -> list[tuple[str, float]]:
        """The ``n`` genes carrying most of this set's weight.

        The interpretable counterpart of the score: it names which genes a
        gene-set activity value is actually reading.
        """
        return sorted(self.items(), key=lambda item: item[1], reverse=True)[: int(n)]
