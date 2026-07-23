"""Gene weighting schemes: how a network prior is turned into per-gene weights.

The scientific problem
----------------------
GSNetAct scores a gene set in a cell as a weighted sum of its members'
expression.  Everything the network contributes is in that weight vector, so the
weighting function *is* the method.

Releases up to 0.0.20 used

.. math::  w_i = \\sum_{j \\in N(i)} A_{ij} \\, \\deg(j)

- a single power-iteration step away from the degree vector.  Across the 1176
STRING subgraphs shipped with the GSNetAct benchmark this gives a Spearman-scale
correlation between weight and degree of 0.96, a median Gini coefficient of 0.61,
and a median 32% of a set's total weight on its top 10% of genes.  The weight
vector is, to a very good approximation, the degree vector.

That is a problem for two reasons, one statistical and one biological.

*Biological.*  Degree in STRING is dominated by how much a gene has been
studied, not by how central it is to any particular programme.  A small number
of genes accumulate interactions across the whole literature (``TP53``, ``AKT1``,
``EGFR``, ``ACTB``), and they appear in many unrelated gene sets.  Weighting by
neighbour degree therefore ranks a gene by *how well studied its neighbours are*.
This "multifunctionality" bias is a documented confounder of guilt-by-association
network analysis (Gillis & Pavlidis 2011, *PLoS ONE* 6:e17258; Ballouz et al.
2015, *Bioinformatics* 31:2123) and of gene-level attention generally (Stoeger et
al. 2018, *PLoS Biol* 16:e2006643).  Meanwhile the genes that actually mark a
pathway's *state* are usually its regulated, peripheral effectors, not its
constitutively expressed hubs.

*Statistical.*  Because hub genes are shared across gene sets, hub-dominated
weights make different sets' scores collinear.  On the benchmark's CITE-seq PBMC
data, 392 gene-set features scored this way have a participation-ratio effective
dimension of only ~15, with 23% of variance on PC1.  Downstream PCA then has
almost nothing left to separate cell states with.  The damage scales with the
size of the library: on 50 pre-selected high-variance sets the hub-weighted and
unweighted representations are indistinguishable (identity NMI 0.591 vs 0.580),
but on the full 529-set collection the hub-weighted representation falls behind
unweighted aggregation (0.511 vs 0.587) - exactly the discovery regime the method
is meant for.

The fix
-------
Damp the influence of connectivity instead of amplifying it.  The default scheme
is a symmetrically normalised ``t``-step diffusion,

.. math::  S = D^{-\\alpha} A D^{-\\alpha}, \\qquad w = S^{t} \\mathbf{1}

with ``alpha = 0.5`` and ``t = 2``, where ``D`` holds the weighted degree
(strength).  Dividing each edge by :math:`\\sqrt{d_i d_j}` is the normalisation
used in spectral clustering (Ng, Jordan & Weiss 2001) and in network propagation
for gene prioritisation (Vanunu et al. 2010, *PLoS Comput Biol* 6:e1000641;
Cowen et al. 2017, *Nat Rev Genet* 18:551).

Two properties make this the right correction.

*It converges to a square-root damping of strength.*  The leading eigenvector of
``S`` is :math:`\\sqrt{d}`, so as ``t`` grows ``w`` approaches
:math:`\\sqrt{\\mathrm{strength}}`.  Measured on a real STRING subgraph, the
correlation between ``w`` and :math:`\\sqrt{d}` is 0.88 at ``t = 1``, 0.97 at
``t = 2`` and 0.9999 by ``t = 50``.  Connectivity still counts - a gene wired
into a programme should matter more than one that is not - but its influence
grows with the square root of the evidence rather than with the evidence times
its neighbours' evidence, so a doubling of literature attention no longer
doubles a gene's say.

The step count is not a sensitive parameter: ``t`` of 1, 2, 3 and 5 give NMI
0.578 / 0.574 / 0.579 / 0.576 on CITE-seq and 0.144 / 0.149 / 0.147 / 0.145 on
IFN-beta, equivalent within seed noise on both.  ``t = 2`` is the default
because it is already close to the limiting behaviour while still cheap; nothing
here rests on the particular value.  ``alpha`` is the parameter that matters.

*It measurably de-concentrates the weights.*  Across the 752 STRING gene-set
subgraphs in the benchmark with at least 20 genes and 20 edges, going from the
legacy scheme to ``alpha = 0.5, t = 2`` moves the median weight Gini from 0.666
to 0.403, the median share of weight held by the top 10% of genes from 39% to
21%, and the median Spearman correlation with degree from 0.93 to 0.79.  Plain
weighted degree (``"strength"``) sits at Gini 0.618 and rho 0.97, i.e. no better
than legacy; PageRank lands in between (0.452, 0.86).

The exponent is not a free parameter chosen for convenience.  ``alpha = 0``
recovers weighted degree - the biased regime - and ``alpha = 1`` over-corrects,
moving weight onto peripheral singletons.  Both benchmark datasets peak at 0.5
and fall away on either side (NMI, CITE-seq / IFN-beta): 0.493 / 0.105 at
``alpha = 0``, 0.543 / 0.142 at 0.25, **0.574 / 0.149 at 0.5**, 0.433 / 0.091 at
0.75, 0.459 / 0.017 at 1.0.

One caveat, stated because it is easy to misread the mechanism: on a *pure star*
the legacy formula happens to give the hub and its leaves identical weights,
because a leaf's single edge is multiplied by the hub's large degree.  The
concentration it produces on real data comes from STRING subgraphs being
assortative - hubs connect to hubs - which turns "weight by your neighbours'
degree" into "weight by your own degree".  The numbers above are measured on the
real subgraphs, not inferred from toy topologies.

Weights are normalised to sum to one per gene set, so a score is a weighted
*mean* of member expression.  Under the old unnormalised weights a score scaled
with the set's edge count and mean degree, which made scores incomparable
between gene sets - invisible if every set is z-scored across cells, but wrong
for any within-cell or cross-dataset comparison.

The 0.0.20 scheme remains available as ``"legacy"`` so published results can be
reproduced exactly.
"""

from __future__ import annotations

from typing import Callable, Mapping

import numpy as np
from scipy import sparse

from .graph import GeneSetGraph

__all__ = [
    "diffusion_weights",
    "legacy_degree_weights",
    "uniform_weights",
    "strength_weights",
    "pagerank_weights",
    "gene_weights",
    "available_weightings",
    "register_weighting",
]

LEGACY_EPSILON = 1e-6


# --------------------------------------------------------------------- schemes
def _normalised_affinity(adjacency: sparse.csr_matrix, alpha: float) -> sparse.csr_matrix:
    """``D^-alpha A D^-alpha`` with zero-degree rows left at zero."""
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    if alpha == 0.0:
        return adjacency
    scale = np.zeros_like(degree)
    positive = degree > 0
    scale[positive] = np.power(degree[positive], -float(alpha))
    diagonal = sparse.diags(scale, format="csr")
    return diagonal @ adjacency @ diagonal


def diffusion_weights(
    graph: GeneSetGraph,
    alpha: float = 0.5,
    steps: int = 2,
) -> np.ndarray:
    """Symmetrically normalised ``steps``-step diffusion weights (the default).

    Computed as ``(D^-alpha A D^-alpha)^steps @ 1`` by repeated sparse
    matrix-vector products, so the cost is ``O(steps * n_edges)`` and no dense
    matrix is ever formed.

    Isolated genes receive weight zero.  Keeping them at zero is deliberate: a
    set member with no within-set STRING partner contributes no *network*
    evidence, and reserving a floor share of the set for such genes measurably
    degrades the representation - NMI falls from 0.593 to 0.461 (CITE-seq) and
    from 0.196 to 0.102 (IFN-beta) with a 5% floor, and no further with 20%.
    The fraction of members this affects is reported per set by
    :func:`gsnetact.runGSNA` so the exclusion is auditable rather than hidden -
    it is a real limitation for poorly annotated genes, not a free choice.
    """
    if steps < 1:
        raise ValueError("steps must be >= 1")
    affinity = _normalised_affinity(graph.adjacency, alpha)
    vector = np.ones(affinity.shape[0], dtype=np.float64)
    for _ in range(int(steps)):
        vector = affinity @ vector
    return np.clip(np.asarray(vector).ravel(), 0.0, None)


def legacy_degree_weights(graph: GeneSetGraph, epsilon: float = LEGACY_EPSILON) -> np.ndarray:
    """The GSNetAct <= 0.0.20 scheme: ``w_i = sum_j A_ij * degree(j)``.

    Retained verbatim for reproducibility of published results.  Isolated genes
    receive ``epsilon``.  Note that release 0.0.17 additionally divided the
    incidence matrix by its total sum; that rescales every gene set by a
    constant and therefore cancels under per-set standardisation, but it does
    change raw score magnitudes.  This function reproduces the 0.0.20 form,
    which applies no such division.

    See the module docstring for why this is no longer the default.
    """
    weights = graph.adjacency @ graph.degree
    weights = np.asarray(weights).ravel()
    weights[graph.degree == 0] = float(epsilon)
    return weights


def uniform_weights(graph: GeneSetGraph) -> np.ndarray:
    """Equal weight per member - the unweighted-aggregation control.

    Any claim that the network contributes information has to beat this.
    """
    return np.ones(len(graph.genes), dtype=np.float64)


def strength_weights(graph: GeneSetGraph) -> np.ndarray:
    """Weighted degree. Exposed as the un-normalised end of the ``alpha`` family."""
    return graph.strength


def pagerank_weights(
    graph: GeneSetGraph,
    damping: float = 0.85,
    max_iter: int = 200,
    tol: float = 1e-10,
) -> np.ndarray:
    """Personalised PageRank with a uniform restart over the set's own members.

    Offered as an alternative propagation prior.  It is degree-normalised on the
    outgoing side only, so it corrects hub bias less completely than
    :func:`diffusion_weights` and scored between the legacy and diffusion
    schemes on both benchmark datasets.
    """
    adjacency = graph.adjacency
    size = adjacency.shape[0]
    if size == 0:
        return np.zeros(0, dtype=np.float64)
    column_sums = np.asarray(adjacency.sum(axis=0)).ravel()
    scale = np.zeros_like(column_sums)
    positive = column_sums > 0
    scale[positive] = 1.0 / column_sums[positive]
    transition = adjacency @ sparse.diags(scale, format="csr")
    restart = np.full(size, 1.0 / size, dtype=np.float64)
    vector = restart.copy()
    for _ in range(int(max_iter)):
        updated = damping * (transition @ vector) + (1.0 - damping) * restart
        if np.abs(updated - vector).sum() < tol:
            vector = updated
            break
        vector = updated
    return np.clip(vector, 0.0, None)


_WEIGHTINGS: dict[str, Callable[..., np.ndarray]] = {
    "diffusion": diffusion_weights,
    "legacy": legacy_degree_weights,
    "uniform": uniform_weights,
    "strength": strength_weights,
    "pagerank": pagerank_weights,
}


def register_weighting(name: str, function: Callable[..., np.ndarray]) -> None:
    """Register a custom weighting under ``name``.

    The callable receives a :class:`~gsnetact.Network.graph.GeneSetGraph` and
    returns one non-negative weight per gene, in ``graph.genes`` order.
    """
    key = str(name)
    if key in _WEIGHTINGS:
        raise ValueError(f"weighting {key!r} is already registered")
    _WEIGHTINGS[key] = function


def available_weightings() -> list[str]:
    return sorted(_WEIGHTINGS)


# ---------------------------------------------------------------- entry point
def gene_weights(
    graph: GeneSetGraph,
    method: str = "diffusion",
    normalize: bool = True,
    options: Mapping[str, object] | None = None,
) -> np.ndarray:
    """Per-gene weights for one gene set.

    Parameters
    ----------
    graph:
        The gene set's interaction subgraph.
    method:
        One of :func:`available_weightings`.  Default ``"diffusion"``.
    normalize:
        Scale the weights to sum to one, making the score a weighted mean of
        member expression and removing gene-set size and network density from
        the score's scale.  ``normalize=False`` reproduces the historical
        un-normalised magnitudes.
    options:
        Extra keyword arguments forwarded to the weighting function
        (e.g. ``{"alpha": 0.5, "steps": 2}``).

    Notes
    -----
    A set whose graph carries no usable weight - no edges at all, or every edge
    removed by context restriction - falls back to uniform weights rather than
    to an all-zero score.  A gene set with no network evidence is then treated
    as exactly what it is: an unweighted gene set.  The fallback is reported in
    the per-set QC returned by :func:`gsnetact.runGSNA`.
    """
    try:
        function = _WEIGHTINGS[str(method)]
    except KeyError:
        raise ValueError(
            f"unknown weighting {method!r}; available: {available_weightings()}"
        ) from None

    weights = np.asarray(function(graph, **dict(options or {})), dtype=np.float64)
    if weights.shape != (len(graph.genes),):
        raise ValueError(
            f"weighting {method!r} returned {weights.shape}, "
            f"expected ({len(graph.genes)},)"
        )
    weights = np.where(np.isfinite(weights), weights, 0.0)
    weights = np.clip(weights, 0.0, None)

    if weights.sum() <= 0.0 and len(weights):
        weights = np.ones_like(weights)
    if normalize and weights.sum() > 0:
        weights = weights / weights.sum()
    return weights


def weight_concentration(weights: np.ndarray) -> float:
    """Gini coefficient of a weight vector: 0 = uniform, 1 = one gene takes all.

    Reported per gene set so that hub domination is a visible, checkable number
    rather than a property one has to take on trust.
    """
    values = np.sort(np.asarray(weights, dtype=np.float64))
    total = values.sum()
    size = values.size
    if size == 0 or total <= 0:
        return 0.0
    index = np.arange(1, size + 1)
    return float(((2 * index - size - 1) * values).sum() / (size * total))
