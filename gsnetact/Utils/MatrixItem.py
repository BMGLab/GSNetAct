"""Deprecated: one entry of the old dense incidence matrix.

Nothing in the package uses this any more. Edges live in a sparse adjacency
matrix (:class:`gsnetact.Network.graph.GeneSetGraph`), which is where the
weighting schemes read them from. Kept only so existing imports do not break.
"""

from dataclasses import dataclass


@dataclass
class MatrixItem:
    """A single incidence-matrix entry: an edge weight and its position.

    .. deprecated:: 0.1.0
       Use :class:`gsnetact.Network.graph.GeneSetGraph` instead.
    """

    weight: float
    row: int
    column: int
