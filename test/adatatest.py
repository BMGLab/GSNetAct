#!/usr/bin/env python
"""End-to-end example: score PBMC3k against the bundled toy network.

Run from this directory:  python adatatest.py
"""

import pandas as pd
import scanpy as sc

from gsnetact import pjson, runGSNA

adata = sc.read_h5ad("test_data/pbmc3k.h5ad")
network = pjson("test_data/deneme.json")

# pbmc3k ships raw counts in .X, so normalise before scoring: a weighted sum of
# raw counts measures library size more than it measures pathway activity.
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

activity = runGSNA(adata, network, layer="lognorm")

print(activity)
print("\nparameters:")
for key, value in activity.uns["gsnetact"].items():
    print(f"  {key}: {value}")

print("\nper-set network support:")
print(activity.var[["n_genes", "n_genes_in_data", "n_edges", "n_isolated_genes",
                    "weight_gini", "uniform_fallback", "top_weighted_gene"]])

first = activity.var_names[0]
weights = activity.uns["gsnetact_gene_weights"][first]
top = sorted(weights.items(), key=lambda item: item[1], reverse=True)[:5]
print(f"\ntop-weighted genes in {first}:")
for gene, weight in top:
    print(f"  {gene:12s} {weight:.4f}")

pd.DataFrame(activity.X, index=activity.obs_names, columns=activity.var_names).to_csv(
    "output.tsv", sep="\t"
)
print("\nwrote output.tsv")
