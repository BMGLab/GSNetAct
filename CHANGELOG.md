# Changelog

## 0.1.0

**The default gene weighting changed.** `runGSNA` no longer weights a gene by
its neighbours' degrees. Scores from 0.1.0 are not comparable with scores from
0.0.x unless you pass `weights="legacy"`, which reproduces the old formula
exactly.

### Why

GSNetAct scores a gene set in a cell as a weighted sum of its members'
expression, so the weighting function is the method. Releases up to 0.0.20 used

```
w_i = sum_j  A_ij * degree(j)
```

Measured across the 752 STRING gene-set subgraphs in the GSNetAct benchmark
with at least 20 genes and 20 edges, that vector is close to the degree vector:
median Spearman correlation with degree **0.93**, median Gini **0.666**, and a
median **39%** of a set's weight on its top 10% of genes.

Degree in STRING tracks how much a gene has been studied more than how central
it is to any one programme. A few genes accumulate interactions across the whole
literature and recur in unrelated gene sets, a bias documented for
guilt-by-association network analysis (Gillis & Pavlidis 2011, *PLoS ONE*
6:e17258; Ballouz et al. 2015, *Bioinformatics* 31:2123) and for gene-level
research attention generally (Stoeger et al. 2018, *PLoS Biol* 16:e2006643).
The genes that mark a pathway's *state* are usually its regulated peripheral
effectors, not its constitutively expressed hubs.

Because hub genes are shared between gene sets, hub-dominated weights also make
different sets' scores collinear. On the benchmark's CITE-seq PBMC data, 392
gene-set features scored the old way have a participation-ratio effective
dimension of ~15, with 23% of variance on PC1 — so a 30-component PCA has almost
nothing left to separate cell states with.

The damage grows with the size of the gene-set collection, which is why it went
unnoticed. On 50 pre-selected high-variance gene sets the old and unweighted
representations are indistinguishable (identity NMI 0.591 vs 0.580). On the full
529-set collection the old weighting falls *behind* unweighted aggregation
(0.511 vs 0.592) — and scoring a whole collection, rather than a hand-picked
shortlist, is the discovery regime the method exists for.

### What replaced it

A symmetrically normalised `t`-step diffusion,

```
S = D^-alpha  A  D^-alpha        w = S^t · 1        alpha = 0.5, t = 2
```

where `D` is the weighted degree (strength). Dividing each edge by
`sqrt(d_i · d_j)` is the normalisation used in spectral clustering (Ng, Jordan &
Weiss 2001) and in network propagation for gene prioritisation (Vanunu et al.
2010, *PLoS Comput Biol* 6:e1000641; Cowen et al. 2017, *Nat Rev Genet* 18:551).

`sqrt(d)` is the leading eigenvector of `S`, so as `t` grows the weight
approaches a **square-root damping of strength** — connectivity still counts, but
its influence grows with the square root of the interaction evidence instead of
with that evidence times the neighbours' evidence. A doubling of literature
attention no longer doubles a gene's say. On a real STRING subgraph the
correlation between `w` and `sqrt(d)` is 0.88 at `t=1`, 0.97 at `t=2` and 0.9999
by `t=50`.

`t` itself is not a sensitive parameter — `t` of 1, 2, 3 and 5 give CITE-seq NMI
0.578 / 0.574 / 0.579 / 0.576 and IFN-β NMI 0.144 / 0.149 / 0.147 / 0.145,
equivalent within seed noise on both. `t=2` is the default because it is already
close to the limiting behaviour and cheap; nothing rests on the particular
value. `alpha` is the parameter that matters.

The effect on the weights is measurable: median Gini **0.666 → 0.403**, median
top-10% weight share **39% → 21%**, median correlation with degree
**0.93 → 0.79**. Plain weighted degree (`"strength"`) does not help (0.618,
0.97); PageRank lands in between (0.452, 0.86).

`alpha` is not a free parameter. `alpha=0` recovers weighted degree — the biased
regime — and `alpha=1` over-corrects onto peripheral singletons. Both benchmark
datasets peak at 0.5 and fall away on either side:

| `alpha` | CITE-seq NMI | IFN-β NMI |
|---|---|---|
| 0.00 (= weighted degree) | 0.493 | 0.105 |
| 0.25 | 0.543 | 0.142 |
| **0.50 (default)** | **0.574** | **0.149** |
| 0.75 | 0.433 | 0.091 |
| 1.00 | 0.459 | 0.017 |

> One caveat, because the mechanism is easy to misread: on a *pure star* the old
> formula happens to give hub and leaves equal weight, since a leaf's single edge
> is multiplied by the hub's large degree. The concentration it produces on real
> data comes from STRING subgraphs being assortative — hubs connect to hubs —
> which turns "weight by your neighbours' degree" into "weight by your own
> degree". Every number above is measured on the real subgraphs.

### Context restriction (new, on by default)

STRING is assembled across tissues. An edge asserts that two proteins associate
*somewhere* in the organism, not that the interaction is available in the cells
being scored. Before computing weights, 0.1.0 keeps an edge only if **both**
partners are detected in at least `min_detection` (default 0.05) of the cells in
the dataset. Nodes are never removed — gene-set membership is a property of the
gene set, not of the dataset — only their edges.

Pruning a prior with the data is standard for prior-guided single-cell methods
(Margolin et al. 2006, ARACNe; Aibar et al. 2017, SCENIC). Detection is used
rather than co-expression deliberately: it answers "is this gene transcribed
here at all", and it does not read the covariance structure that the downstream
clustering is supposed to discover, so it cannot manufacture the signal it is
later credited with finding.

When a set loses every edge, weighting falls back to uniform — a gene set with no
context-supported network evidence is then treated as exactly what it is, an
unweighted gene set. `var["uniform_fallback"]` records this per set.

Every threshold tested beats no restriction on both datasets, but where the
optimum sits is dataset-dependent (NMI, with the number of sets that fall back
to uniform):

| `min_detection` | CITE-seq NMI | fallback / 529 | IFN-β NMI | fallback / 48 |
|---|---|---|---|---|
| off | 0.574 | 0 | 0.149 | 0 |
| 0.01 | 0.572 | 84 | 0.159 | 0 |
| 0.02 | 0.588 | 114 | 0.159 | 1 |
| **0.05 (default)** | **0.593** | 153 | **0.196** | 2 |
| 0.10 | 0.597 | 182 | 0.233 | 3 |
| 0.20 | 0.592 | 225 | 0.241 | 9 |

CITE-seq plateaus by 0.05 while IFN-β keeps improving to 0.20. The default is
0.05 — the conservative end, discarding the least prior information — rather than
the per-dataset optimum, because tuning a default on two datasets where they
disagree would be fitting to them. Raise it if your data are shallow or your
gene sets come from a different tissue.

That 153 of 529 C8 sets lose all network evidence at the default is the expected
result, not a failure: the C8 collections are built from other tissues, so a
hepatocyte programme scored on PBMCs *should* end up with no context-supported
topology. On the tissue-matched Hallmark collection only 2 of 48 fall back.

### Evidence

All numbers from the shipped `runGSNA`; downstream stack held fixed (robust z per
set → PCA 30 → 15-NN → Leiden, 5 resolutions × 3 seeds, median reported).

CITE-seq PBMC, full 529-set C8 label-blind collection, protein-defined labels:

| configuration | NMI | ARI |
|---|---|---|
| legacy weighting | 0.511 | 0.359 |
| unweighted aggregation (MeanSet control) | 0.592 | 0.443 |
| diffusion, no context restriction | 0.574 | 0.435 |
| **diffusion + context restriction (0.1.0 default)** | **0.593** | **0.464** |
| pagerank + context restriction | 0.557 | 0.411 |

Kang IFN-β PBMC, 48 Hallmark sets, condition labels:

| configuration | NMI | ARI |
|---|---|---|
| legacy weighting | 0.109 | 0.039 |
| unweighted aggregation (MeanSet control) | 0.157 | 0.079 |
| diffusion, no context restriction | 0.149 | 0.081 |
| **diffusion + context restriction (0.1.0 default)** | **0.196** | **0.092** |
| pagerank + context restriction | 0.120 | 0.045 |

On CITE-seq, 153 of 529 C8 sets fall back to uniform because they come from
other tissues. Restricting to the **376 sets that retain network evidence**, so
no fallback is involved, isolates the weighting itself:

| configuration | NMI | ARI |
|---|---|---|
| legacy weighting | 0.484 | 0.327 |
| unweighted aggregation (MeanSet control) | 0.570 | 0.387 |
| **diffusion + context restriction** | **0.603** | **0.468** |

This is the comparison the GSNetAct benchmark could not previously win: network
weighting exceeding unweighted aggregation on the same retained sets.

**Sequencing depth.** Binomially thinning the counts and rescoring (one thinning
replicate per level, median NMI):

| retained counts | \| | CITE-seq legacy | unweighted | 0.1.0 | \| | IFN-β legacy | unweighted | 0.1.0 |
|---|---|---|---|---|---|---|---|---|
| 100% | | 0.511 | 0.592 | 0.593 | | 0.109 | 0.157 | 0.196 |
| 50% | | 0.504 | 0.569 | 0.576 | | 0.005 | 0.131 | 0.160 |
| 25% | | 0.501 | 0.572 | 0.573 | | 0.003 | 0.088 | 0.101 |
| 10% | | 0.378 | 0.514 | 0.541 | | 0.003 | 0.075 | 0.101 |

The 0.1.0 default is at or above unweighted aggregation at every depth on both
datasets, and its margin is largest where counts are lowest — the direction a
useful prior should move in, since when dropout removes a gene its network
neighbours still carry the programme.

The legacy weighting collapses to near-chance on IFN-β below full depth (0.109 →
0.005 at 50%). Concentrating a set's weight on a handful of genes means that
when those genes drop out the score becomes noise, whereas a de-concentrated
weight vector degrades gracefully. This is one thinning replicate per level, so
treat the ordering as robust and the exact values as indicative.

### Bugs fixed

Each was reproduced against the pre-0.1.0 code before being fixed.

* **`expScore` mutated the caller's expression matrix.** `adata.X += 1e-6` ran
  once per gene set, in place, and accumulated across sets and across re-runs. A
  zero matrix scored against 500 gene sets came back at `5.0e-4`. Scoring now
  never touches its input; `expScore(..., epsilon=1e-6)` reproduces the offset
  without the mutation.
* **Sparse input crashed.** `expScore` replaced `adata.X` with a dense array but
  kept scoring the stale sparse reference, and `sparse += scalar` raises
  `NotImplementedError: adding a nonzero scalar to a sparse array is not
  supported`. Sparse matrices are now scored natively.
* **Self-loops crashed.** `{"A": {"A": 0.5}}` produced an incidence column with
  one non-zero, and `i, j = np.nonzero(column)` raised `ValueError: not enough
  values to unpack`. Self-loops are now dropped: a gene is not its own
  interaction partner.
* **Edges were keyed by `hash(a) * hash(b) * hash(score)`.** A product of hashes
  is not injective, so a collision would merge two distinct edges into one
  column and raise the same unpack error from the other side. Edges are now
  keyed on the gene pair, which also removes the `PYTHONHASHSEED` dependence.
* **`pip install .` failed on a clean clone** — `setup.py` read
  `README_PACKAGE.md`, which is git-ignored and not in the repository. Packaging
  moved to `pyproject.toml`.
* **Asymmetric and duplicate edges** are merged with `max` rather than summed,
  so reporting an edge in both directions no longer doubles its weight.
* **Edgeless gene sets** no longer divide by zero.
* Gene sets dropped because they have no members are now named in a warning
  instead of silently disappearing from the output.

### STRING retrieval

* **`required_score` default 100 → 400.** 100 is STRING's lowest possible cut; on
  its scale 150 is "low", 400 "medium", 700 "high", 900 "highest", and the
  combined score is a posterior probability that the association is real. A
  threshold of 100 admits links more likely wrong than right. This is the most
  likely explanation for the benchmark's topology-permutation result: across
  nine within-set weight-permutation tests the smallest BH-adjusted q-value was
  0.985, and in six of nine the real weights scored *below* the permuted mean. A
  prior indistinguishable from its own permutation is not carrying prior
  information.
* **`channels=` (new).** STRING's combined score merges seven evidence channels,
  two of which — `coexpression` and `textmining` — are themselves derived largely
  from expression compendia and the literature those produced. Using them to
  weight an expression-based score is partly circular. `channels=["experiments",
  "database"]` rebuilds the score from experimental and curated evidence only,
  using STRING's own independence formula. Default behaviour is unchanged.
* **`limit` removed.** It belongs to STRING's `interaction_partners` endpoint,
  not to `network`, and was silently ignored. The equivalent parameter is
  `add_nodes`, now explicit and defaulting to 0 — adding partners from outside
  the gene set would change what the gene set is.
* Requests are rate-limited and carry a `caller_identity`, as STRING's access
  terms ask. Results are collected with `as_completed` instead of by blocking on
  submission order, and written in input order so the JSON is reproducible.

### Performance and API

* Scoring a collection is now one sparse `(cells × genes) @ (genes × sets)`
  product instead of a Python loop with a dense dot product per set.
* The dense `(n_genes × n_edges)` incidence matrix is gone. A 300-gene set with a
  dense STRING neighbourhood needed ~100 MB as an incidence matrix and needs
  <1 MB as a sparse adjacency. `GeneSetMatrix` still exists and still returns
  `(matrix, err)`, but is built lazily and only if something asks for it.
* `runGSNA` gained `weights`, `layer`, `min_detection`, `restrict_to_context`,
  `normalize`, `weight_options` and `seed`. `normalized=` is deprecated but
  still works (`True` → `normalize="cell"`, `False` → `normalize=None`).
* `normalize="set"` is the new default: robust z per gene set across cells. The
  old within-cell transform (`normalize="cell"`) imposes a compositional
  constraint — every cell forced to the same distribution of pathway activity,
  so a globally quiescent cell cannot look quiescent. Right for within-cell
  comparisons, wrong for clustering.
* Weights are normalised to sum to one per set, so a score is a weighted *mean*
  of member expression. Previously a score scaled with the set's edge count and
  mean degree, making sets incomparable — invisible under per-set z-scoring, but
  wrong for any within-cell or cross-dataset comparison.
* The returned AnnData now carries per-set QC in `var` (member and edge counts,
  isolated-gene count and weight share, weight Gini, largest-component fraction,
  uniform-fallback flag, top-weighted gene), the run's parameters in
  `uns["gsnetact"]`, the full per-set weight vectors in
  `uns["gsnetact_gene_weights"]`, and unnormalised scores in
  `layers["raw_activity"]`.
* New: `GeneSetGraph`, `gene_weights`, `register_weighting`,
  `weight_concentration`, `detection_rate`, `NetworkGeneWeights`,
  `scoreMatrix`, `buildWeightMatrix`, `normalizeScores`.
* A warning fires if `var_names` look like Ensembl IDs, which would otherwise
  score every set as zero.
* `pjson` validates structure on load and accepts an already-parsed dict.
* 38 unit tests, including a regression test for every bug above.

### Known limitations

* Isolated genes get weight zero. A set member with no within-set STRING partner
  contributes no *network* evidence, and reserving a floor share of the set for
  such genes measurably degrades the representation — NMI falls 0.593 → 0.461
  (CITE-seq) and 0.196 → 0.102 (IFN-β) with a 5% floor, and no further with 20%.
  This does mean poorly annotated genes are invisible to the weighting;
  `var["n_isolated_genes"]` and `var["isolated_weight_fraction"]` report how many
  members that affects, so it can be audited rather than assumed.
* The weights are unsigned. A gene set mixes activators and repressors, and a
  negative-feedback gene induced *by* pathway activity (NFKBIA in NF-κB
  signalling, for instance) counts toward activity as if it were an effector.
  Signed membership is not yet supported.
* Context restriction uses detection rate, which is confounded with expression
  level. It is a deliberately weak filter and is not a claim about which
  interactions are physically occurring.
* STRING is context-agnostic even after restriction, so a benefit from the
  topology is evidence for informative prior structure, not for tissue-specific
  causality.
