"""Gene-set objects.

The public names (:class:`GeneSetMatrix`, :class:`GeneSet`, :func:`getGSNA`) are
unchanged so existing scripts keep working, but the machinery underneath is now
the sparse :class:`~gsnetact.Network.graph.GeneSetGraph`.  Three failure modes
of the old dense-incidence-matrix construction are gone as a result:

* a self-loop (``{"A": {"A": 0.5}}``) produced an incidence column with one
  non-zero, and the scorer's ``i, j = np.nonzero(column)`` raised
  ``ValueError: not enough values to unpack``;
* edges were keyed by ``hash(a) * hash(b) * hash(score)``.  A product of hashes
  is not injective, so a collision merged two distinct edges into one column and
  raised the same unpack error from the other direction.  Keying on the ordered
  gene pair removes both the collision risk and the dependence on
  ``PYTHONHASHSEED``;
* the incidence matrix was ``O(n_genes * n_edges)`` dense, which for a 300-gene
  set with a dense STRING neighbourhood is ~100 MB for one gene set.

:class:`GeneSet` also gains :attr:`GeneSet.graph`, which is what new code should
use.
"""

from __future__ import annotations

import warnings
from typing import Mapping

import numpy as np

from ..Network.graph import GeneSetGraph

__all__ = ["GeneSetMatrix", "GeneSet", "getGSNA"]


class GeneSetMatrix:
    """Weighted incidence matrix of a gene set, kept for backward compatibility.

    Returns ``(matrix, error_flag)`` exactly as before: ``matrix`` is
    ``(n_genes x n_edges)`` with the edge weight in both endpoint rows, and
    ``error_flag`` is ``1`` when the gene set could not be processed.

    New code should use :class:`~gsnetact.Network.graph.GeneSetGraph` instead -
    every quantity this matrix was used for is a function of the adjacency
    matrix alone, at ``O(E)`` rather than ``O(V*E)`` memory.

    Unlike release 0.0.17, the matrix is **not** divided by its total sum.  That
    division (removed upstream in 0.0.20) rescales every gene set by a constant,
    which cancels under per-set standardisation but changes raw score
    magnitudes.
    """

    def __new__(cls, rawGeneSet: Mapping[str, Mapping[str, float]], _id: str = ""):
        graph = GeneSetGraph(rawGeneSet, name=str(_id))

        if len(graph) == 0:
            warnings.warn(
                f"gene set {_id!r} is empty and was skipped.",
                RuntimeWarning,
                stacklevel=2,
            )
            return None, 1

        upper = np.triu(graph.adjacency.toarray(), k=1)
        rows, cols = np.nonzero(upper)
        matrix = np.zeros((len(graph), len(rows)), dtype=np.float64)
        for column, (left, right) in enumerate(zip(rows, cols)):
            weight = upper[left, right]
            matrix[left, column] = weight
            matrix[right, column] = weight
        return matrix, 0


class GeneSet:
    """One gene set: its identifier, its raw JSON entry, and its network."""

    def __init__(self, _id: str, _rawGeneSet: Mapping[str, Mapping[str, float]]):
        self.id = str(_id)
        self.asJson = _rawGeneSet
        self.graph = GeneSetGraph(_rawGeneSet, name=self.id)
        self.err = 0 if len(self.graph) else 1
        self._matrix = None

    @property
    def matrix(self) -> np.ndarray | None:
        """Legacy incidence matrix, built lazily.

        Only materialised if something asks for it, so the ``O(V*E)`` cost is no
        longer paid by every run.
        """
        if self._matrix is None and not self.err:
            self._matrix, _ = GeneSetMatrix(self.asJson, self.id)
        return self._matrix

    @property
    def getID(self) -> str:
        return self.id

    @property
    def getAsJson(self):
        return self.asJson

    @property
    def getGeneNames(self) -> list[str]:
        return list(self.graph.genes)

    def __repr__(self) -> str:  # pragma: no cover - display only
        return f"GeneSet({self.id!r}, genes={len(self.graph)}, edges={self.graph.n_edges})"


def getGSNA(jsonFile: Mapping[str, Mapping[str, Mapping[str, float]]]) -> list[GeneSet]:
    """Build :class:`GeneSet` objects for every entry in ``jsonFile``.

    Sets that cannot be processed (no members at all) are dropped, with a warning
    naming each one, so a silently shorter output can be traced back.
    """
    sets: list[GeneSet] = []
    skipped: list[str] = []
    for name in jsonFile:
        candidate = GeneSet(f"{name}", jsonFile[name])
        if candidate.err:
            skipped.append(str(name))
            continue
        sets.append(candidate)
    if skipped:
        warnings.warn(
            f"{len(skipped)} gene set(s) had no members and were skipped: "
            + ", ".join(skipped[:5])
            + (" ..." if len(skipped) > 5 else ""),
            RuntimeWarning,
            stacklevel=2,
        )
    return sets
