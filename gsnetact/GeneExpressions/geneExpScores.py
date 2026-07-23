"""Turning gene weights into per-cell gene-set activity.

Two bugs in the pre-0.1.0 implementation are fixed here.

1. ``adata.X += 1e-6`` mutated the caller's expression matrix *in place*, once
   per gene set.  Scoring 500 gene sets therefore left the user's ``adata.X``
   permanently shifted by ``500e-6``, and the shift accumulated further on every
   re-run.  Verified: a zero matrix scored against 500 sets came back at
   ``5.0e-4``.  Scoring must never modify its input, so the offset is applied to
   a copy - and only when it is requested, since it exists purely to reproduce
   historical numbers.
2. Sparse input crashed.  ``expScore`` converted ``adata.X`` to dense but kept
   scoring the stale sparse reference, and ``sparse += scalar`` raises
   ``NotImplementedError: adding a nonzero scalar to a sparse array is not
   supported``.  Sparse matrices are now scored natively, which is also what
   makes whole-collection scoring tractable: the product becomes a single sparse
   ``(cells x genes) @ (genes x sets)`` multiplication instead of one dense dot
   product per gene set.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
from scipy import sparse

__all__ = ["expScore", "scoreMatrix", "buildWeightMatrix"]


def _expression(adata: Any, layer: str | None):
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(
            f"layer {layer!r} not found; available layers: {list(adata.layers)}"
        )
    return adata.layers[layer]


def expScore(
    adata: Any,
    geneSetScore: Mapping[str, float],
    layer: str | None = None,
    epsilon: float = 0.0,
) -> np.ndarray:
    """Activity of one gene set across all cells.

    Parameters
    ----------
    adata:
        AnnData whose ``var_names`` are gene symbols. **Never modified.**
    geneSetScore:
        ``{gene: weight}``.  Genes absent from ``adata.var_names`` are ignored;
        genes in ``adata`` but not in the mapping get weight zero.
    layer:
        Layer to score.  ``None`` uses ``adata.X``.  Pass the log-normalised
        layer explicitly when ``X`` holds raw counts - a weighted sum of raw
        counts is dominated by library size.
    epsilon:
        Historical per-entry offset.  Left at ``0.0``; set to ``1e-6`` only to
        reproduce pre-0.1.0 output.  The offset shifts every cell's score by the
        same constant, so it cannot affect any downstream analysis that
        standardises scores, but it made the function destructive.
    """
    weights = np.array(
        [float(geneSetScore.get(str(name), 0.0)) for name in adata.var_names],
        dtype=np.float64,
    )
    matrix = _expression(adata, layer)

    if sparse.issparse(matrix):
        scores = np.asarray(matrix @ weights).ravel()
    else:
        scores = np.asarray(matrix, dtype=np.float64) @ weights
    if epsilon:
        # Equivalent to adding `epsilon` to every expression entry, without
        # touching the caller's matrix.
        scores = scores + float(epsilon) * weights.sum()
    return scores


def buildWeightMatrix(
    perSetWeights: Mapping[str, Mapping[str, float]],
    varNames: Sequence[str],
) -> tuple[sparse.csr_matrix, list[str]]:
    """Stack per-set gene weights into a sparse ``(genes x sets)`` matrix.

    Scoring a whole collection is then one sparse matrix product rather than a
    Python loop over gene sets, which is where nearly all of the speed-up over
    the pre-0.1.0 implementation comes from.
    """
    position = {str(name): index for index, name in enumerate(varNames)}
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    names = list(perSetWeights)
    for column, name in enumerate(names):
        for gene, weight in perSetWeights[name].items():
            row = position.get(str(gene))
            if row is None or weight == 0.0:
                continue
            rows.append(row)
            cols.append(column)
            data.append(float(weight))
    matrix = sparse.csr_matrix(
        (data, (rows, cols)), shape=(len(position), len(names)), dtype=np.float64
    )
    return matrix, names


def scoreMatrix(
    adata: Any,
    weightMatrix: sparse.csr_matrix,
    layer: str | None = None,
) -> np.ndarray:
    """Score every gene set at once: ``(cells x genes) @ (genes x sets)``."""
    matrix = _expression(adata, layer)
    if sparse.issparse(matrix):
        scores = matrix @ weightMatrix
        if sparse.issparse(scores):
            scores = scores.todense()
        return np.asarray(scores, dtype=np.float64)
    dense = np.asarray(matrix, dtype=np.float64)
    return np.asarray(weightMatrix.T @ dense.T, dtype=np.float64).T
