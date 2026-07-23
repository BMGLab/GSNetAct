#!/usr/bin/env python
"""Lower-level example: inspect one gene set's network and its gene weights.

Shows what `runGSNA` does per set, and compares the 0.1.0 default weighting with
the pre-0.1.0 one on the same graph.

Run from this directory:  python basicTest.py
"""

from gsnetact import (
    GeneSetGraph,
    NetworkGeneWeights,
    pjson,
    weight_concentration,
)
from gsnetact.Network.weights import diffusion_weights, legacy_degree_weights

network = pjson("test_data/deneme.json")

for name in network:
    graph = GeneSetGraph(network[name], name)
    print(f"\n=== {name} ===")
    print(f"  {len(graph)} genes, {graph.n_edges} edges, "
          f"density {graph.density:.3f}, "
          f"largest component {graph.largest_component_fraction():.2f}")

    if graph.n_edges == 0:
        print("  no within-set interactions; weighting falls back to uniform")
        continue

    legacy = legacy_degree_weights(graph)
    diffusion = diffusion_weights(graph)
    print(f"  weight Gini   legacy {weight_concentration(legacy):.3f}"
          f"   diffusion {weight_concentration(diffusion):.3f}")

    weights = NetworkGeneWeights(graph)          # 0.1.0 default
    print("  gene weights (0.1.0 default, normalised to sum 1):")
    for gene, weight in weights.top(5):
        print(f"    {gene:12s} {weight:.4f}   degree={int(graph.degree[graph.genes.index(gene)])}")
    print(f"  isolated members: {weights.stats['n_isolated_genes']}"
          f"/{weights.stats['n_genes']}")
