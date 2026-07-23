"""Unit tests for GSNetAct 0.1.0.

Each test is either a hand-computable numeric check, or a regression test for a
bug that was verified to exist in the pre-0.1.0 code.
"""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData
from scipy import sparse

from gsnetact import (
    GeneSetGraph,
    GeneSetMatrix,
    GeneSetScore,
    diffusion_weights,
    expScore,
    gene_weights,
    getGSNA,
    legacy_degree_weights,
    pjson,
    runGSNA,
    uniform_weights,
    weight_concentration,
)


# --------------------------------------------------------------------- fixtures
def path_graph():
    """A-B-C: B is the only gene with two partners."""
    return {"A": {"B": 0.5}, "B": {"A": 0.5, "C": 0.4}, "C": {"B": 0.4}}


def star_graph(n_leaves=6, score=0.8):
    """One hub connected to ``n_leaves`` mutually unconnected leaves."""
    leaves = [f"L{i}" for i in range(n_leaves)]
    nodes = {"HUB": {leaf: score for leaf in leaves}}
    nodes.update({leaf: {"HUB": score} for leaf in leaves})
    return nodes


def core_periphery_graph(n_hubs=5, n_leaves=10):
    """Assortative graph: mutually connected hubs, each carrying its own leaves.

    This is the regime real STRING gene-set subgraphs sit in - a pure star is
    maximally *dis*assortative and behaves differently under both schemes.
    """
    hubs = [f"H{i}" for i in range(n_hubs)]
    nodes = {hub: {other: 0.9 for other in hubs if other != hub} for hub in hubs}
    for index, hub in enumerate(hubs):
        for leaf_index in range(n_leaves):
            leaf = f"P{index}_{leaf_index}"
            nodes[leaf] = {hub: 0.7}
            nodes[hub][leaf] = 0.7
    return nodes


def toy_adata(n_cells=8, genes=("A", "B", "C"), sparse_x=False, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.random((n_cells, len(genes))) + 0.1
    adata = AnnData(X=sparse.csr_matrix(X) if sparse_x else X)
    adata.var_names = list(genes)
    return adata


# ----------------------------------------------------------------------- graph
def test_graph_is_symmetric_and_deduplicated():
    graph = GeneSetGraph(path_graph(), "path")
    dense = graph.adjacency.toarray()
    assert np.allclose(dense, dense.T)
    assert graph.n_edges == 2
    assert dense[0, 1] == pytest.approx(0.5)
    assert dense[1, 2] == pytest.approx(0.4)
    assert np.allclose(np.diag(dense), 0.0)


def test_asymmetric_input_is_merged_with_max_not_summed():
    graph = GeneSetGraph({"A": {"B": 0.9}, "B": {"A": 0.3}}, "asym")
    assert graph.adjacency.toarray()[0, 1] == pytest.approx(0.9)


def test_self_loops_are_dropped():
    """Pre-0.1.0 this raised ValueError from `i, j = np.nonzero(column)`."""
    graph = GeneSetGraph({"A": {"A": 0.7, "B": 0.3}, "B": {"A": 0.3}}, "selfloop")
    assert graph.n_edges == 1
    assert np.allclose(np.diag(graph.adjacency.toarray()), 0.0)
    matrix, err = GeneSetMatrix({"A": {"A": 0.7, "B": 0.3}, "B": {"A": 0.3}}, "selfloop")
    assert err == 0
    GeneSetScore("selfloop", matrix, ["A", "B"])  # must not raise


def test_isolated_members_are_kept_as_nodes():
    graph = GeneSetGraph({"A": {"B": 0.5}, "B": {"A": 0.5}, "ORPHAN": {}}, "orphan")
    assert graph.genes == ["A", "B", "ORPHAN"]
    assert graph.degree[2] == 0


def test_partners_outside_the_set_are_ignored():
    graph = GeneSetGraph({"A": {"B": 0.5, "OUTSIDE": 0.9}, "B": {"A": 0.5}}, "outside")
    assert graph.n_edges == 1
    assert graph.genes == ["A", "B"]


def test_edge_recovery_is_hash_seed_independent():
    """Edges are keyed on the gene pair, not on a product of hashes."""
    size = 40
    clique = {
        f"G{i}": {f"G{j}": 0.5 for j in range(size) if j != i} for i in range(size)
    }
    graph = GeneSetGraph(clique, "clique")
    assert graph.n_edges == size * (size - 1) // 2


def test_largest_component_fraction():
    two_components = {
        "A": {"B": 0.9}, "B": {"A": 0.9},
        "C": {"D": 0.9, "E": 0.9}, "D": {"C": 0.9}, "E": {"C": 0.9},
    }
    graph = GeneSetGraph(two_components, "split")
    assert graph.largest_component_fraction() == pytest.approx(3 / 5)


# --------------------------------------------------------------------- weights
def test_legacy_weights_match_the_published_formula_by_hand():
    """w_i = sum_j A_ij * degree(j) on A-B-C.

    degrees: A=1, B=2, C=1
    w_A = 0.5 * deg(B) = 1.0
    w_B = 0.5 * deg(A) + 0.4 * deg(C) = 0.9
    w_C = 0.4 * deg(B) = 0.8
    """
    graph = GeneSetGraph(path_graph(), "path")
    assert legacy_degree_weights(graph) == pytest.approx([1.0, 0.9, 0.8])


def test_legacy_weights_give_isolated_genes_epsilon():
    graph = GeneSetGraph({"A": {"B": 0.5}, "B": {"A": 0.5}, "ORPHAN": {}}, "orphan")
    assert legacy_degree_weights(graph)[2] == pytest.approx(1e-6)


def test_diffusion_weights_match_a_hand_computation():
    """S = D^-1/2 A D^-1/2 on A-B-C, then w = S^2 @ 1.

    strengths: A=0.5, B=0.9, C=0.4
    S_AB = 0.5/sqrt(0.5*0.9), S_BC = 0.4/sqrt(0.9*0.4)
    """
    graph = GeneSetGraph(path_graph(), "path")
    A = graph.adjacency.toarray()
    d = A.sum(1)
    S = A / np.sqrt(np.outer(d, d))
    expected = (S @ S) @ np.ones(3)
    assert diffusion_weights(graph) == pytest.approx(expected)


def test_diffusion_de_concentrates_weight_on_an_assortative_graph():
    """The central claim, on the topology STRING subgraphs actually have.

    Five mutually connected hubs, each carrying ten leaves: hubs connect to
    hubs, which is the assortative regime that turns the legacy scheme's
    "weight by your neighbours' degree" into "weight by your own degree".
    """
    graph = GeneSetGraph(core_periphery_graph(), "core-periphery")
    legacy = legacy_degree_weights(graph)
    diffusion = diffusion_weights(graph)
    hubs = [graph.genes.index(f"H{i}") for i in range(5)]

    assert legacy[hubs].sum() / legacy.sum() > 0.35
    assert diffusion[hubs].sum() / diffusion.sum() < 0.20
    assert weight_concentration(diffusion) < weight_concentration(legacy)


def test_diffusion_converges_to_a_square_root_damping_of_strength():
    """sqrt(strength) is the leading eigenvector of D^-1/2 A D^-1/2.

    This is what makes the scheme a *damping* of connectivity: influence grows
    with the square root of the interaction evidence, not with the evidence
    times the neighbours' evidence.
    """
    graph = GeneSetGraph(core_periphery_graph(), "core-periphery")
    root_strength = np.sqrt(graph.strength)
    far = diffusion_weights(graph, alpha=0.5, steps=40)
    assert np.corrcoef(far, root_strength)[0, 1] > 0.999
    # ...while the default t=2 stays short of the limit, keeping local structure
    near = diffusion_weights(graph, alpha=0.5, steps=2)
    assert np.corrcoef(near, root_strength)[0, 1] > 0.9


def test_alpha_zero_recovers_weighted_degree():
    """The unnormalised end of the family, and the biased regime."""
    graph = GeneSetGraph(core_periphery_graph(), "core-periphery")
    assert diffusion_weights(graph, alpha=0.0, steps=1) == pytest.approx(graph.strength)


def test_gene_weights_normalise_to_one():
    graph = GeneSetGraph(path_graph(), "path")
    for method in ("diffusion", "legacy", "uniform", "strength", "pagerank"):
        assert gene_weights(graph, method).sum() == pytest.approx(1.0)


def test_edgeless_set_falls_back_to_uniform_not_to_zero():
    graph = GeneSetGraph({"A": {}, "B": {}, "C": {}}, "edgeless")
    assert gene_weights(graph, "diffusion") == pytest.approx([1 / 3, 1 / 3, 1 / 3])


def test_unknown_weighting_is_rejected():
    graph = GeneSetGraph(path_graph(), "path")
    with pytest.raises(ValueError, match="unknown weighting"):
        gene_weights(graph, "not_a_method")


def test_weight_concentration_bounds():
    assert weight_concentration(np.ones(10)) == pytest.approx(0.0)
    winner_takes_all = np.zeros(10)
    winner_takes_all[0] = 1.0
    assert weight_concentration(winner_takes_all) > 0.85


def test_uniform_weights_are_the_meanset_control():
    graph = GeneSetGraph(star_graph(), "star")
    assert np.allclose(uniform_weights(graph), 1.0)


# --------------------------------------------------------- context restriction
def test_context_restriction_removes_edges_but_never_members():
    graph = GeneSetGraph(path_graph(), "path")
    restricted = graph.restrict(np.array([True, True, False]))
    assert restricted.genes == ["A", "B", "C"]
    assert restricted.n_edges == 1                 # B-C is gone
    assert restricted.degree[2] == 0


def test_undetected_genes_lose_their_weight():
    adata = toy_adata(n_cells=20, genes=("A", "B", "C"))
    adata.X[:, 2] = 0.0                            # C detected in 0% of cells
    adata.X[0, 2] = 1.0                            # ... except one, i.e. 5%
    network = {"S": path_graph()}
    scored = runGSNA(adata, network, min_detection=0.5, normalize=None)
    weights = scored.uns["gsnetact_gene_weights"]["S"]
    assert weights["C"] == pytest.approx(0.0)
    assert weights["A"] > 0 and weights["B"] > 0


# --------------------------------------------------------------------- scoring
def test_expscore_does_not_modify_the_input():
    """Pre-0.1.0 this shifted adata.X by 1e-6 per gene set, cumulatively."""
    adata = toy_adata()
    before = adata.X.copy()
    expScore(adata, {"A": 1.0, "B": 2.0, "C": 3.0})
    assert np.array_equal(before, adata.X)


def test_rungsna_does_not_modify_the_input_over_many_sets():
    adata = toy_adata()
    before = adata.X.copy()
    network = {f"S{i}": path_graph() for i in range(200)}
    runGSNA(adata, network, normalize=None)
    assert np.array_equal(before, adata.X)


def test_sparse_input_is_scored_natively():
    """Pre-0.1.0 this raised NotImplementedError on `sparse += scalar`."""
    dense = toy_adata(sparse_x=False)
    sparse_adata = AnnData(X=sparse.csr_matrix(dense.X))
    sparse_adata.var_names = list(dense.var_names)
    network = {"S": path_graph()}
    a = runGSNA(dense, network, restrict_to_context=False, normalize=None)
    b = runGSNA(sparse_adata, network, restrict_to_context=False, normalize=None)
    assert np.allclose(a.X, b.X, atol=1e-6)


def test_expscore_epsilon_reproduces_the_historical_offset():
    adata = toy_adata()
    weights = {"A": 1.0, "B": 2.0, "C": 3.0}
    plain = expScore(adata, weights)
    offset = expScore(adata, weights, epsilon=1e-6)
    assert offset - plain == pytest.approx(np.full(adata.n_obs, 1e-6 * 6.0))


def test_score_equals_a_weighted_mean_of_member_expression():
    adata = toy_adata(n_cells=5)
    network = {"S": path_graph()}
    scored = runGSNA(adata, network, restrict_to_context=False, normalize=None)
    weights = scored.uns["gsnetact_gene_weights"]["S"]
    manual = adata.X @ np.array([weights[g] for g in adata.var_names])
    assert np.allclose(scored.X.ravel(), manual, atol=1e-6)


def test_uniform_weighting_reproduces_the_plain_set_mean():
    adata = toy_adata(n_cells=5)
    network = {"S": path_graph()}
    scored = runGSNA(adata, network, weights="uniform",
                     restrict_to_context=False, normalize=None)
    assert np.allclose(scored.X.ravel(), adata.X.mean(axis=1), atol=1e-6)


def test_layer_selection():
    adata = toy_adata(n_cells=6)
    adata.layers["lognorm"] = np.log1p(adata.X)
    network = {"S": path_graph()}
    on_x = runGSNA(adata, network, restrict_to_context=False, normalize=None)
    on_layer = runGSNA(adata, network, layer="lognorm",
                       restrict_to_context=False, normalize=None)
    assert not np.allclose(on_x.X, on_layer.X)
    assert on_layer.uns["gsnetact"]["layer"] == "lognorm"


def test_missing_layer_is_a_clear_error():
    adata = toy_adata()
    with pytest.raises(KeyError):
        runGSNA(adata, {"S": path_graph()}, layer="nope")


# ---------------------------------------------------------------- output shape
def test_output_carries_per_set_qc_and_provenance():
    adata = toy_adata(n_cells=10)
    network = {"path": path_graph(), "star": star_graph()}
    adata = toy_adata(n_cells=10, genes=("A", "B", "C", "HUB", "L0", "L1",
                                         "L2", "L3", "L4", "L5"))
    scored = runGSNA(adata, network, restrict_to_context=False)
    assert scored.shape == (10, 2)
    assert list(scored.var_names) == ["path", "star"]
    for column in ("n_genes", "n_edges", "n_isolated_genes", "weight_gini",
                   "largest_component_fraction", "top_weighted_gene"):
        assert column in scored.var.columns
    assert scored.var.loc["star", "top_weighted_gene"] in {"HUB", "L0"}
    assert scored.uns["gsnetact"]["weighting"] == "diffusion"
    assert scored.uns["gsnetact"]["n_gene_sets"] == 2
    assert "raw_activity" in scored.layers


def test_obs_is_carried_through():
    adata = toy_adata(n_cells=6)
    adata.obs["celltype"] = ["a", "b"] * 3
    scored = runGSNA(adata, {"S": path_graph()}, restrict_to_context=False)
    assert list(scored.obs["celltype"]) == ["a", "b"] * 3


def test_normalisation_modes():
    adata = toy_adata(n_cells=50, genes=("A", "B", "C"), seed=3)
    network = {f"S{i}": path_graph() for i in range(6)}
    by_set = runGSNA(adata, network, restrict_to_context=False, normalize="set")
    assert np.allclose(np.median(by_set.X, axis=0), 0.0, atol=1e-5)
    by_cell = runGSNA(adata, network, restrict_to_context=False, normalize="cell")
    assert by_cell.X.shape == (50, 6)
    raw = runGSNA(adata, network, restrict_to_context=False, normalize=None)
    assert raw.X.min() > 0


def test_deprecated_normalized_flag_still_works():
    adata = toy_adata(n_cells=12)
    network = {f"S{i}": path_graph() for i in range(4)}
    with pytest.deprecated_call():
        scored = runGSNA(adata, network, normalized=False, restrict_to_context=False)
    assert scored.uns["gsnetact"]["normalize"] is None


def test_ensembl_ids_warn():
    adata = toy_adata(genes=("ENSG00000141510", "ENSG00000012048", "ENSG00000139618"))
    with pytest.warns(RuntimeWarning, match="Ensembl"):
        runGSNA(adata, {"S": path_graph()}, restrict_to_context=False)


def test_empty_gene_set_is_skipped_with_a_warning():
    with pytest.warns(RuntimeWarning):
        sets = getGSNA({"good": path_graph(), "empty": {}})
    assert [s.id for s in sets] == ["good"]


def test_pjson_rejects_a_malformed_file(tmp_path):
    import json as _json

    bad = tmp_path / "bad.json"
    bad.write_text(_json.dumps({"S": {"A": "not-a-dict"}}))
    with pytest.raises(ValueError, match="must map partners"):
        pjson(bad)


def test_pjson_round_trip(tmp_path):
    import json as _json

    good = tmp_path / "good.json"
    good.write_text(_json.dumps({"S": path_graph()}))
    parsed = pjson(good)
    assert parsed.getGeneSetCount == 1
    assert parsed.getUniqueGeneNames == ["A", "B", "C"]


# ------------------------------------------------------------- reproducibility
def test_scoring_is_deterministic():
    adata = toy_adata(n_cells=30, seed=11)
    network = {f"S{i}": path_graph() for i in range(5)}
    first = runGSNA(adata, network, restrict_to_context=False)
    second = runGSNA(adata, network, restrict_to_context=False)
    assert np.array_equal(first.X, second.X)


def test_legacy_weighting_is_still_reachable_end_to_end():
    adata = toy_adata(n_cells=10)
    scored = runGSNA(adata, {"S": path_graph()}, weights="legacy",
                     restrict_to_context=False, normalize=None)
    weights = scored.uns["gsnetact_gene_weights"]["S"]
    # normalised version of the hand-computed [1.0, 0.9, 0.8]
    assert weights["A"] == pytest.approx(1.0 / 2.7)
    assert weights["B"] == pytest.approx(0.9 / 2.7)
    assert weights["C"] == pytest.approx(0.8 / 2.7)
