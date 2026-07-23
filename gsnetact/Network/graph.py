"""Sparse graph representation of one gene set's interaction network.

Why this replaces the incidence matrix
--------------------------------------
Releases up to 0.0.20 materialised a dense ``(n_genes x n_edges)`` incidence
matrix per gene set and then looped over its columns in Python.  Every quantity
that construction was used for is a function of the weighted adjacency matrix
alone, which is what this module builds instead:

* memory drops from ``O(V*E)`` to ``O(E)`` -- a 300-gene set with a dense STRING
  neighbourhood needs ~45k edge columns x 300 rows of float64 (~108 MB) as an
  incidence matrix, and ~0.7 MB as a sparse adjacency;
* edges are keyed by the ordered gene pair rather than by
  ``hash(a) * hash(b) * hash(score)``.  A product of hashes is not injective, so
  two distinct edges could land in the same column, and the resulting column had
  four non-zeros instead of two -- which the old scorer's
  ``i, j = np.nonzero(column)`` unpack turned into a ``ValueError``.  Self-loops
  produced the same crash from the other side (one non-zero);
* the node order is the JSON insertion order and no longer depends on
  ``PYTHONHASHSEED``.

Biological note on symmetry: STRING functional links are undirected evidence
statements, so the adjacency is symmetrised explicitly and duplicate ``(a, b)``
/ ``(b, a)`` entries are merged with ``max``.  A gene is never its own
interaction partner, so self-loops are dropped rather than counted as
connectivity.
"""

from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import numpy as np
from scipy import sparse

__all__ = ["GeneSetGraph"]


class GeneSetGraph:
    """The STRING-induced subgraph of a single gene set.

    Parameters
    ----------
    nodes:
        ``{gene: {partner: score}}`` as produced by :func:`gsnetact.makeJson`.
        Genes mapping to an empty dict are set members that STRING reported no
        within-set partner for; they are kept as isolated nodes so that gene-set
        membership is never silently redefined by the network.
    name:
        Optional gene-set identifier, used only in messages.

    Attributes
    ----------
    genes:
        Node names in JSON insertion order.
    adjacency:
        Symmetric CSR matrix of interaction scores, zero diagonal.
    """

    __slots__ = ("name", "genes", "adjacency", "_index")

    def __init__(self, nodes: Mapping[str, Mapping[str, float]], name: str = "") -> None:
        self.name = str(name)
        self.genes = [str(gene) for gene in nodes]
        self._index = {gene: position for position, gene in enumerate(self.genes)}

        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []
        for gene, partners in nodes.items():
            left = self._index[str(gene)]
            for partner, score in (partners or {}).items():
                right = self._index.get(str(partner))
                if right is None or right == left:
                    # Partners outside the set, and self-loops, carry no
                    # within-set connectivity information.
                    continue
                weight = float(score)
                if not np.isfinite(weight) or weight <= 0.0:
                    continue
                rows.append(left)
                cols.append(right)
                data.append(weight)

        size = len(self.genes)
        directed = sparse.coo_matrix(
            (data, (rows, cols)), shape=(size, size), dtype=np.float64
        ).tocsr()
        # ``maximum`` merges the two directions of every reported edge without
        # doubling its weight, which summation would do.
        self.adjacency = directed.maximum(directed.T)
        self.adjacency.eliminate_zeros()

    # ------------------------------------------------------------------ views
    def __len__(self) -> int:
        return len(self.genes)

    def __repr__(self) -> str:  # pragma: no cover - display only
        return (
            f"GeneSetGraph({self.name!r}, genes={len(self.genes)}, "
            f"edges={self.n_edges}, isolated={int((self.degree == 0).sum())})"
        )

    @property
    def n_edges(self) -> int:
        """Number of distinct undirected edges."""
        return int(self.adjacency.nnz // 2)

    @property
    def degree(self) -> np.ndarray:
        """Unweighted degree (number of within-set partners) per gene."""
        return np.asarray((self.adjacency > 0).sum(axis=1)).ravel().astype(np.float64)

    @property
    def strength(self) -> np.ndarray:
        """Weighted degree: the sum of incident interaction scores per gene."""
        return np.asarray(self.adjacency.sum(axis=1)).ravel()

    @property
    def density(self) -> float:
        size = len(self.genes)
        if size < 2:
            return 0.0
        return 2.0 * self.n_edges / (size * (size - 1))

    def index_of(self, genes: Iterable[str]) -> np.ndarray:
        """Positions of ``genes`` in :attr:`genes`; ``-1`` where absent."""
        return np.array([self._index.get(str(gene), -1) for gene in genes], dtype=np.int64)

    # ------------------------------------------------------------- operations
    def restrict(self, keep: Sequence[bool] | np.ndarray) -> "GeneSetGraph":
        """Return a copy whose edges are confined to the ``keep`` nodes.

        Nodes are retained -- only their edges are removed -- so the gene set
        keeps its full membership while the *prior topology* is narrowed to the
        part of the graph the data can support.  See
        :mod:`gsnetact.Network.context` for the biological argument.
        """
        mask = np.asarray(keep, dtype=bool)
        if mask.shape != (len(self.genes),):
            raise ValueError(
                f"keep mask has shape {mask.shape}, expected ({len(self.genes)},)"
            )
        restricted = object.__new__(GeneSetGraph)
        restricted.name = self.name
        restricted.genes = list(self.genes)
        restricted._index = dict(self._index)
        diagonal = sparse.diags(mask.astype(np.float64), format="csr")
        adjacency = diagonal @ self.adjacency @ diagonal
        adjacency.eliminate_zeros()
        restricted.adjacency = adjacency.tocsr()
        return restricted

    def largest_component_fraction(self) -> float:
        """Share of connected genes that sit in the largest connected component.

        A gene set whose STRING subgraph shatters into many small components is
        not one coherent programme, and a single activity score for it is harder
        to defend.  Reported per set so this can be audited rather than assumed.
        """
        connected = self.degree > 0
        n_connected = int(connected.sum())
        if n_connected == 0:
            return 0.0
        n_components, labels = sparse.csgraph.connected_components(
            self.adjacency, directed=False
        )
        del n_components
        sizes = np.bincount(labels[connected])
        return float(sizes.max() / n_connected)
