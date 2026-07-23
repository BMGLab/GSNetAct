"""``runGSNA``: network-weighted gene-set activity for an AnnData object.

What changed in 0.1.0
---------------------
* The default gene weighting is a symmetrically normalised two-step diffusion
  instead of the degree-based formula.  See :mod:`gsnetact.Network.weights` for
  the argument and the benchmark numbers; ``weights="legacy"`` reproduces
  pre-0.1.0 output.
* The prior network is restricted to genes actually detected in the data before
  weights are computed (:mod:`gsnetact.Network.context`).
* Weights sum to one per gene set, so a score is a weighted mean of member
  expression rather than a quantity that scales with the set's edge count.
* Scoring is one sparse matrix product for the whole collection, and never
  modifies the input AnnData.
* The returned object carries per-set QC in ``var`` and the run's parameters in
  ``uns``, so a score can be traced to the network evidence behind it.
"""

from __future__ import annotations

import warnings
from typing import Any, Mapping

import numpy as np
from anndata import AnnData

from ..GeneExpressions.geneExpScores import buildWeightMatrix, scoreMatrix
from ..GeneSets.geneSetObjects import getGSNA
from ..GeneSets.geneSetScores import NetworkGeneWeights
from ..Network.context import detection_mask, detection_rate
from ..Network.weights import LEGACY_EPSILON

__all__ = ["runGSNA", "normalizeScores"]

_ENSEMBL_PREFIX = "ENSG"


def normalizeScores(scores: np.ndarray, method: str | None, seed: int = 42) -> np.ndarray:
    """Standardise a ``(cells x sets)`` score matrix.

    ``"cell"``
        The historical transform: quantile-normalise each *cell* across gene
        sets, then standardise each gene set across cells.  Makes activities
        comparable *within* a cell, at the cost of imposing a compositional
        constraint - every cell is forced to the same distribution of pathway
        activity, so a globally quiescent cell cannot look quiescent.  Use it
        for within-cell comparisons, not for clustering.
    ``"set"``
        Robust z-score down each gene set across cells,
        ``(x - median) / (1.4826 * MAD)``, clipped at +/-10.  No cell-level
        constraint, so this is the right default for clustering and for
        differential-activity testing.  The MAD is used rather than the standard
        deviation because single-cell activity distributions are heavy-tailed: a
        handful of high-activity cells would otherwise inflate the scale and
        compress exactly the signal being looked for.
    ``None``
        Raw weighted means, on the expression layer's own scale.
    """
    values = np.asarray(scores, dtype=np.float64)
    if method in (None, "none"):
        return values.astype(np.float32)

    if method == "set":
        median = np.nanmedian(values, axis=0)
        mad = np.nanmedian(np.abs(values - median), axis=0) * 1.4826
        std = np.nanstd(values, axis=0)
        scale = np.where(mad > 1e-12, mad, np.where(std > 1e-12, std, 1.0))
        standardized = np.clip((values - median) / scale, -10.0, 10.0)
        standardized[~np.isfinite(standardized)] = 0.0
        return standardized.astype(np.float32)

    if method == "cell":
        from sklearn.preprocessing import StandardScaler, quantile_transform

        n_quantiles = max(2, min(1000, values.shape[1]))
        transformed = quantile_transform(
            values,
            axis=1,
            n_quantiles=n_quantiles,
            output_distribution="normal",
            random_state=int(seed),
            copy=True,
        )
        return StandardScaler().fit_transform(transformed).astype(np.float32)

    raise ValueError(f"unknown normalization {method!r}; use 'set', 'cell' or None")


def _check_var_names(adata: AnnData) -> None:
    sample = [str(name) for name in list(adata.var_names)[:50]]
    if sample and sum(name.startswith(_ENSEMBL_PREFIX) for name in sample) > len(sample) // 2:
        warnings.warn(
            "adata.var_names look like Ensembl gene IDs, but gene-set networks "
            "are keyed on HGNC symbols. Convert var_names to symbols first, or "
            "every set will score zero.",
            RuntimeWarning,
            stacklevel=3,
        )


def runGSNA(
    adata: AnnData,
    jsonFile: Mapping[str, Mapping[str, Mapping[str, float]]],
    weights: str = "diffusion",
    layer: str | None = None,
    min_detection: float = 0.05,
    restrict_to_context: bool = True,
    normalize: str | None = "set",
    weight_options: Mapping[str, Any] | None = None,
    normalized: bool | None = None,
    seed: int = 42,
) -> AnnData:
    """Score every gene set in ``jsonFile`` across every cell in ``adata``.

    Parameters
    ----------
    adata:
        Cells x genes, ``var_names`` as HGNC symbols. **Never modified.**
    jsonFile:
        ``{set: {gene: {partner: score}}}`` - a :class:`~gsnetact.pjson` or any
        mapping with that shape.
    weights:
        Weighting scheme. ``"diffusion"`` (default), ``"legacy"``, ``"uniform"``,
        ``"strength"``, ``"pagerank"``. See
        :func:`gsnetact.Network.weights.available_weightings`.
    layer:
        Layer to score; ``None`` uses ``adata.X``.  Pass a log-normalised layer
        when ``X`` holds raw counts.
    min_detection:
        A gene must be detected in at least this fraction of cells for its edges
        to count.  ``0.05`` by default; the benchmark optimum is broad over
        0.02-0.20.
    restrict_to_context:
        Apply that restriction. ``False`` uses the prior network as published.
    normalize:
        ``"set"`` (default), ``"cell"``, or ``None``; see :func:`normalizeScores`.
    weight_options:
        Extra arguments for the weighting function, e.g.
        ``{"alpha": 0.5, "steps": 2}``.
    normalized:
        Deprecated pre-0.1.0 flag. ``True`` maps to ``normalize="cell"``,
        ``False`` to ``normalize=None``.
    seed:
        Seed for the ``"cell"`` normalisation's quantile transform.

    Returns
    -------
    AnnData
        ``(cells x gene sets)`` activity, ``obs`` carried over from ``adata``.
        ``var`` holds per-set QC (member and edge counts, isolated-gene count and
        weight share, weight Gini, largest-component fraction, uniform-fallback
        flag, and the set's top-weighted gene); ``uns["gsnetact"]`` records the
        parameters, and ``uns["gsnetact_gene_weights"]`` the full per-set weight
        vectors for interpretation.
    """
    if normalized is not None:
        warnings.warn(
            "`normalized` is deprecated; use normalize='cell' (was True) or "
            "normalize=None (was False).",
            DeprecationWarning,
            stacklevel=2,
        )
        normalize = "cell" if normalized else None

    _check_var_names(adata)

    expression = adata.X if layer is None else adata.layers[layer]
    var_names = [str(name) for name in adata.var_names]

    keep_genes: set[str] | None = None
    detection = None
    if restrict_to_context:
        detection = detection_rate(expression)
        mask = detection_mask(detection, min_detection)
        keep_genes = {gene for gene, ok in zip(var_names, mask) if ok}
        if not keep_genes:
            warnings.warn(
                f"no gene reaches min_detection={min_detection}; context "
                "restriction disabled for this run.",
                RuntimeWarning,
                stacklevel=2,
            )
            keep_genes = None

    measured = set(var_names)
    gene_sets = getGSNA(jsonFile)
    if not gene_sets:
        raise ValueError("no usable gene set in jsonFile")

    per_set_weights: dict[str, dict[str, float]] = {}
    records: list[dict[str, Any]] = []
    options = dict(weight_options or {})
    if weights == "legacy" and "epsilon" not in options:
        options["epsilon"] = LEGACY_EPSILON

    for gene_set in gene_sets:
        graph = gene_set.graph
        if keep_genes is not None:
            # Edges are kept only between genes that are both measured and
            # detected often enough for the interaction to be available here.
            graph = graph.restrict(
                np.array([gene in keep_genes for gene in graph.genes], dtype=bool)
            )
        computed = NetworkGeneWeights(
            graph, method=weights, normalize=True, options=options
        )
        per_set_weights[gene_set.id] = dict(computed)

        stats = dict(computed.stats)
        stats["gene_set"] = gene_set.id
        stats["n_genes_in_data"] = int(sum(gene in measured for gene in graph.genes))
        top = computed.top(1)
        stats["top_weighted_gene"] = top[0][0] if top else ""
        stats["top_weight"] = float(top[0][1]) if top else 0.0
        records.append(stats)

    weight_matrix, set_names = buildWeightMatrix(per_set_weights, var_names)
    raw = scoreMatrix(adata, weight_matrix, layer=layer)
    values = normalizeScores(raw, normalize, seed=seed)

    scored = AnnData(X=np.asarray(values, dtype=np.float32), obs=adata.obs.copy())
    scored.var_names = [str(name) for name in set_names]
    for key in records[0]:
        if key == "gene_set":
            continue
        scored.var[key] = [record[key] for record in records]
    scored.layers["raw_activity"] = np.asarray(raw, dtype=np.float32)

    matched = int((np.asarray(weight_matrix.sum(axis=1)).ravel() > 0).sum())
    scored.uns["gsnetact"] = {
        "version": _package_version(),
        "weighting": str(weights),
        "weight_options": options,
        "layer": layer,
        "normalize": normalize,
        "restrict_to_context": bool(restrict_to_context),
        "min_detection": float(min_detection) if restrict_to_context else None,
        "n_gene_sets": len(set_names),
        "n_genes_used": matched,
        "n_genes_in_data": len(var_names),
        "median_detection_rate": float(np.median(detection)) if detection is not None else None,
    }
    scored.uns["gsnetact_gene_weights"] = per_set_weights
    return scored


def _package_version() -> str:
    from .. import __version__

    return __version__
