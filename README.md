# GSNetAct

Network-weighted gene-set activity for single-cell data.

GSNetAct scores a gene set in a cell as a weighted sum of its members'
expression, where the weights come from the interactions between those genes in
STRING. The idea is that a gene wired into the programme should count for more
than a gene that merely appears in the list.

**Version 0.1.0 changes the default weighting.** The old scheme weighted a gene
by its neighbours' degrees, which in practice reproduced STRING degree — a proxy
for how much a gene has been studied — and made gene-set scores collinear. The
new default damps connectivity instead of amplifying it, and restricts the prior
network to genes actually detected in the data. On the GSNetAct benchmark this is
the difference between falling behind unweighted aggregation and beating it. See
[CHANGELOG.md](CHANGELOG.md) for the argument, the numbers, and the bugs fixed
along the way. `weights="legacy"` reproduces 0.0.x exactly.

## Install

```bash
pip install gsnetact
```

or from source:

```bash
git clone https://github.com/BMGLab/GSNetAct
cd GSNetAct
pip install -e .
```

Requires Python 3.9+, NumPy, SciPy, AnnData, scikit-learn and Requests.

## Usage

```python
import scanpy as sc
from gsnetact import pjson, runGSNA

adata = sc.read_h5ad("pbmc.h5ad")     # var_names must be HGNC symbols
network = pjson("geneSets.json")      # see "Building the network" below

activity = runGSNA(adata, network, layer="lognorm")
```

`activity` is an AnnData of cells × gene sets. Pass the log-normalised layer
explicitly — a weighted sum of raw counts is dominated by library size.

The result carries its own provenance:

```python
activity.var                       # per-set QC, one row per gene set
activity.uns["gsnetact"]           # the parameters this run used
activity.uns["gsnetact_gene_weights"]["HALLMARK_INFLAMMATORY_RESPONSE"]
activity.layers["raw_activity"]    # scores before standardisation
```

`var` reports, for every set: how many members it has, how many are in your data,
how many edges survived, how many members ended up isolated and what share of the
weight they hold, the weight Gini, the largest-component fraction, whether the
set fell back to uniform weighting, and its top-weighted gene. A score you cannot
trace to network evidence is a score you should not trust.

### Which genes is a score reading?

```python
from gsnetact import GeneSetGraph, NetworkGeneWeights

graph = GeneSetGraph(network["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
weights = NetworkGeneWeights(graph)
weights.top(10)          # the ten genes carrying most of the set's weight
weights.stats            # the same QC as in var
```

### Options

```python
runGSNA(
    adata, network,
    weights="diffusion",        # "legacy" | "uniform" | "strength" | "pagerank"
    layer="lognorm",
    min_detection=0.05,         # edges need both partners detected this often
    restrict_to_context=True,
    normalize="set",            # "cell" | None
    weight_options={"alpha": 0.5, "steps": 2},
)
```

**Reproducing the GSNetAct paper.** The default (`alpha=0.5`, `min_detection=0.05`)
was fit for cell-identity scoring. For **coordinated condition/state programmes**
— the regime the paper focuses on (e.g. interferon stimulation) — the tuned
configuration is stronger hub damping and a stricter context filter:

```python
runGSNA(adata, network, weights="diffusion", layer="lognorm",
        min_detection=0.2, weight_options={"alpha": 1.0, "steps": 2})
```

`alpha` is a genuine tuning knob, not a universal constant; sweep it for your
signal rather than assuming the identity-fit default (see CHANGELOG).

**`weights`** — `"diffusion"` (default) weights a gene by
`(D^-0.5 A D^-0.5)^2 · 1`. `sqrt(strength)` is the leading eigenvector of that
operator, so influence grows with the square root of the interaction evidence
rather than with the evidence times the neighbours' evidence. `"uniform"` ignores
the network entirely and is the control any network claim has to beat.
`"legacy"` is the pre-0.1.0 formula.

**`min_detection`** — STRING asserts that two proteins associate *somewhere* in
the organism, not that the interaction is available in the cells you are
scoring. Edges between genes not transcribed here are dropped. Members are never
dropped, only their edges. A set that loses every edge falls back to uniform
weighting, which is the honest treatment of a gene set with no context-supported
network evidence.

**`normalize`** — `"set"` (default) robust-z-scores each gene set across cells.
`"cell"` is the pre-0.1.0 within-cell transform: it makes activities comparable
*within* a cell but forces every cell to the same distribution of pathway
activity, so a globally quiescent cell cannot look quiescent. Use `"cell"` for
within-cell comparisons, `"set"` for clustering and differential testing.

## Building the network

```bash
makeGeneSets --geneSymbols h.all.v2026.1.Hs.symbols.gmt --fileType gmt -o geneSets.json
```

or

```python
from gsnetact import makeJson

makeJson("h.all.v2026.1.Hs.symbols.gmt", "gmt", "geneSets.json",
         required_score=400,
         channels=["experiments", "database"])   # optional
```

`required_score` is STRING's confidence floor on a 0–1000 scale: 150 low, 400
medium (the default), 700 high, 900 highest. The pre-0.1.0 default was 100,
STRING's lowest possible cut, which admits links more likely wrong than right.

`channels` restricts the score to chosen evidence channels, recombined with
STRING's own independence formula. `["experiments", "database"]` gives a network
free of expression-derived and text-mined evidence — worth using when the
downstream analysis is itself expression-based, since STRING's `coexpression`
and `textmining` channels partly encode the co-expression structure the method
is then credited with finding.

## Network JSON format

```json
{
  "GeneSet1": {
    "Gene1": {"Gene2": 0.35, "Gene3": 0.77},
    "Gene2": {"Gene1": 0.35, "Gene3": 0.51},
    "Gene3": {"Gene1": 0.77, "Gene2": 0.51},
    "Gene4": {}
  }
}
```

Scores are STRING combined scores in `[0, 1]`. The graph is undirected; both
directions may be listed and are merged with `max`. A member with no within-set
partner is written with an empty dict — it stays in the gene set and receives no
network weight.

![Graph for GeneSet1](/genesets.png)

## Choosing a weighting

Measured across 752 STRING gene-set subgraphs (≥20 genes, ≥20 edges):

| scheme | weight Gini | ρ(weight, degree) | top-10% weight share |
|---|---|---|---|
| `legacy` | 0.666 | 0.93 | 39% |
| `strength` | 0.618 | 0.97 | 37% |
| `pagerank` | 0.452 | 0.86 | 28% |
| **`diffusion`** | **0.403** | **0.79** | **21%** |

Lower is less hub-dominated. Clustering benchmarks for each are in
[CHANGELOG.md](CHANGELOG.md).

## Tests

```bash
pip install -e ".[test]"
pytest
```

## Limitations

* Weights are unsigned. A gene set mixes activators and repressors, and a
  negative-feedback gene induced *by* pathway activity counts toward activity as
  if it were an effector.
* Genes with no within-set STRING partner get weight zero, so poorly annotated
  genes are invisible to the weighting. `var["n_isolated_genes"]` reports how
  many members that affects per set.
* Context restriction uses detection rate, which is confounded with expression
  level. It is a deliberately weak filter, not a claim about which interactions
  are physically occurring.
* STRING is context-agnostic even after restriction: a benefit from the topology
  is evidence for informative prior structure, not for tissue-specific causality.

## Citation

Please cite this repository if you use this code in your work.
