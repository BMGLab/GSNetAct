"""Build gene-set interaction networks from STRING.

Two choices here are scientific, not cosmetic.

**Confidence threshold.**  The pre-0.1.0 default was ``required_score=100``,
STRING's lowest possible cut.  On STRING's own scale 150 is "low", 400 "medium",
700 "high" and 900 "highest" confidence, and the combined score is a posterior
probability that the association is real - so a 100 threshold admits links that
are more likely wrong than right.  The topology of a graph built that way is
close to noise, which is what the GSNetAct benchmark observed: across nine
within-set weight-permutation tests the smallest BH-adjusted q-value was 0.985,
and in six of nine the *real* weights scored below the permuted mean.  A prior
indistinguishable from its own permutation is not carrying prior information.
The default is now 400.

**Evidence channels.**  STRING's combined score merges seven channels, two of
which - ``coexpression`` and ``textmining`` - are themselves derived largely from
expression compendia and from the literature those compendia produced.  Using
them to weight an expression-based score is partly circular: the network
"knows" the co-expression structure the method is then credited with finding.
``channels`` lets a caller build the network from experimental and curated
evidence only (``["experiments", "database"]``), which is the conservative
choice for a methods claim.  The default keeps the combined score so existing
behaviour is unchanged, but the option is there and the trade-off is documented.

Other corrections: ``limit`` belongs to STRING's ``interaction_partners``
endpoint, not to ``network``, and was silently ignored - the equivalent
parameter is ``add_nodes``, now explicit and defaulting to 0 (never add genes
the gene set does not contain).  Requests are rate-limited and carry a
``caller_identity``, as STRING's access terms ask.  Results are collected with
``as_completed`` rather than by blocking on submission order, and the output is
written in the input gene sets' order so the JSON is reproducible.
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Iterable, Mapping, Sequence

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

__all__ = ["makeJson", "makeJson_console", "parse_gmt", "parse_tsv", "create_session"]

STRING_API_URL = "https://string-db.org/api/json/network"
CALLER_IDENTITY = "gsnetact"

#: STRING's own confidence bands, for reference.
CONFIDENCE_BANDS = {"low": 150, "medium": 400, "high": 700, "highest": 900}

#: Per-channel score fields returned by the STRING JSON API.
CHANNEL_FIELDS = {
    "neighborhood": "nscore",
    "fusion": "fscore",
    "cooccurence": "pscore",
    "coexpression": "ascore",
    "experiments": "escore",
    "database": "dscore",
    "textmining": "tscore",
}


class _RateLimiter:
    """At most ``per_second`` requests, shared across worker threads."""

    def __init__(self, per_second: float):
        self.interval = 1.0 / max(float(per_second), 1e-6)
        self.lock = threading.Lock()
        self.last = 0.0

    def wait(self) -> None:
        with self.lock:
            delay = self.interval - (time.monotonic() - self.last)
            if delay > 0:
                time.sleep(delay)
            self.last = time.monotonic()


def create_session(retries: int = 5, backoff: float = 0.5) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=retries,
        connect=retries,
        read=retries,
        status=retries,
        backoff_factor=backoff,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def _combined_from_channels(edge: Mapping[str, Any], channels: Sequence[str]) -> float:
    """Recombine selected channels using STRING's own independence formula.

    STRING combines channel probabilities as ``1 - prod(1 - s_c)`` after removing
    the random background ``p = 0.041``; the same formula restricted to a subset
    of channels gives the confidence one would have from that evidence alone.
    """
    prior = 0.041
    product = 1.0
    for channel in channels:
        field = CHANNEL_FIELDS[channel]
        raw = float(edge.get(field, 0.0) or 0.0)
        if raw > 1.0:
            raw /= 1000.0
        corrected = max(0.0, (raw - prior) / (1.0 - prior))
        product *= 1.0 - corrected
    combined = 1.0 - product
    return combined + prior * (1.0 - combined)


def fetch_filtered_string_network(
    genes: Sequence[str],
    session: requests.Session,
    species: int = 9606,
    required_score: int = 400,
    add_nodes: int = 0,
    network_type: str = "functional",
    channels: Sequence[str] | None = None,
    timeout: tuple[int, int] = (30, 600),
) -> list[dict[str, Any]] | None:
    """Fetch the STRING subnetwork induced by ``genes``.

    Only edges whose *both* endpoints are in ``genes`` are kept: the network is a
    property of the gene set, and adding partners from outside it would change
    which genes the set contains.
    """
    payload = {
        "identifiers": "\r".join(genes),
        "species": int(species),
        "required_score": int(required_score),
        "network_type": network_type,
        "add_nodes": int(add_nodes),
        "caller_identity": CALLER_IDENTITY,
    }
    try:
        response = session.post(STRING_API_URL, data=payload, timeout=timeout)
        response.raise_for_status()
        edges = response.json()
    except (requests.exceptions.RequestException, ValueError) as exc:
        print(f"STRING request failed: {exc!r}", file=sys.stderr)
        return None

    if not edges:
        return []

    wanted = {str(gene) for gene in genes}
    kept: list[dict[str, Any]] = []
    for edge in edges:
        left = str(edge.get("preferredName_A", ""))
        right = str(edge.get("preferredName_B", ""))
        if left not in wanted or right not in wanted or left == right:
            continue
        if channels:
            score = _combined_from_channels(edge, channels)
            if score * 1000.0 < required_score:
                continue
        else:
            score = float(edge.get("score", 0.0))
            if score > 1.0:
                score /= 1000.0
        if score <= 0:
            continue
        kept.append({"a": left, "b": right, "score": score})
    return kept


def _fetch_all(
    gene_sets: Mapping[str, Mapping[str, Sequence[str]]],
    max_workers: int,
    requests_per_second: float,
    **kwargs: Any,
) -> dict[str, list[dict[str, Any]] | None]:
    limiter = _RateLimiter(requests_per_second)
    results: dict[str, list[dict[str, Any]] | None] = {}

    def task(name: str) -> tuple[str, list[dict[str, Any]] | None]:
        limiter.wait()
        with create_session() as session:
            return name, fetch_filtered_string_network(
                gene_sets[name]["geneSymbols"], session, **kwargs
            )

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(task, name): name for name in gene_sets}
        for done, future in enumerate(as_completed(futures), start=1):
            name = futures[future]
            try:
                name, edges = future.result()
            except Exception as exc:  # noqa: BLE001 - one set must not kill the run
                print(f"gene set {name!r} failed: {exc!r}", file=sys.stderr)
                edges = None
            results[name] = edges
            if done % 25 == 0 or done == len(futures):
                print(f"STRING networks: {done}/{len(futures)}", file=sys.stderr)
    return results


def parse_gmt(file_path: str) -> dict[str, dict[str, list[str]]]:
    """Read an MSigDB GMT: ``name <tab> description <tab> gene ...``."""
    gene_sets: dict[str, dict[str, list[str]]] = {}
    with open(file_path) as handle:
        for line in handle:
            parts = [field.strip() for field in line.rstrip("\n").split("\t")]
            if len(parts) < 3:
                continue
            genes = [gene for gene in parts[2:] if gene]
            if genes:
                gene_sets[parts[0]] = {"geneSymbols": sorted(set(genes))}
    return gene_sets


def parse_tsv(file_path: str) -> dict[str, dict[str, list[str]]]:
    """Read a headerless TSV: ``name <tab> gene ...``."""
    gene_sets: dict[str, dict[str, list[str]]] = {}
    with open(file_path) as handle:
        for line in handle:
            parts = [field.strip() for field in line.rstrip("\n").split("\t")]
            if len(parts) < 2:
                continue
            genes = [gene for gene in parts[1:] if gene]
            if genes:
                gene_sets[parts[0]] = {"geneSymbols": sorted(set(genes))}
    return gene_sets


def makeJson(
    msigdbFile: str,
    fileType: str = "gmt",
    jsonFileName: str = "geneSets.json",
    required_score: int = 400,
    species: int = 9606,
    network_type: str = "functional",
    channels: Iterable[str] | None = None,
    add_nodes: int = 0,
    max_workers: int = 8,
    requests_per_second: float = 8.0,
) -> dict[str, dict[str, dict[str, float]]]:
    """Build ``{set: {gene: {partner: score}}}`` for every set in ``msigdbFile``.

    Parameters
    ----------
    required_score:
        STRING confidence floor on the 0-1000 scale.  Default 400 ("medium").
        See the module docstring for why the old default of 100 was a problem.
    channels:
        Restrict the score to these evidence channels, recombined with STRING's
        own formula.  ``["experiments", "database"]`` gives a network free of
        expression-derived and text-mined evidence, the conservative choice when
        the downstream analysis is itself expression-based.  ``None`` (default)
        uses STRING's combined score.
    add_nodes:
        Interaction partners to add from outside the gene set.  Kept at 0:
        adding them would change what the gene set is.

    Notes
    -----
    Genes for which STRING reports no within-set partner are written with an
    empty dict rather than dropped, so gene-set membership survives the network
    step intact.
    """
    if fileType == "json":
        with open(msigdbFile) as handle:
            gene_sets = json.load(handle)
    elif fileType == "gmt":
        gene_sets = parse_gmt(msigdbFile)
    elif fileType == "tsv":
        gene_sets = parse_tsv(msigdbFile)
    else:
        raise ValueError(f"unsupported file type: {fileType!r} (use gmt, tsv or json)")

    if not gene_sets:
        raise ValueError(f"no gene set parsed from {msigdbFile!r}")

    channel_list = list(channels) if channels else None
    if channel_list:
        unknown = set(channel_list) - set(CHANNEL_FIELDS)
        if unknown:
            raise ValueError(
                f"unknown STRING channel(s): {sorted(unknown)}; "
                f"available: {sorted(CHANNEL_FIELDS)}"
            )

    fetched = _fetch_all(
        gene_sets,
        max_workers=max_workers,
        requests_per_second=requests_per_second,
        species=species,
        required_score=required_score,
        add_nodes=add_nodes,
        network_type=network_type,
        channels=channel_list,
    )

    relations: dict[str, dict[str, dict[str, float]]] = {}
    failed: list[str] = []
    # Iterate over the input order, not the completion order, so the JSON is
    # byte-reproducible across runs.
    for name in gene_sets:
        edges = fetched.get(name)
        if edges is None:
            failed.append(name)
            continue
        nodes: dict[str, dict[str, float]] = {
            gene: {} for gene in gene_sets[name]["geneSymbols"]
        }
        for edge in edges:
            left, right, score = edge["a"], edge["b"], float(edge["score"])
            nodes.setdefault(left, {})
            nodes.setdefault(right, {})
            nodes[left][right] = max(score, nodes[left].get(right, 0.0))
            nodes[right][left] = max(score, nodes[right].get(left, 0.0))
        relations[name] = {gene: dict(sorted(nodes[gene].items())) for gene in sorted(nodes)}

    if failed:
        print(
            f"WARNING: {len(failed)} gene set(s) could not be retrieved and are "
            f"absent from the output: {', '.join(failed[:5])}"
            + (" ..." if len(failed) > 5 else ""),
            file=sys.stderr,
        )

    with open(jsonFileName, "w") as handle:
        json.dump(relations, handle, indent=2, sort_keys=False)
    return relations


def makeJson_console() -> None:
    parser = argparse.ArgumentParser(
        description="Build gene-set interaction networks from STRING."
    )
    parser.add_argument("--geneSymbols", required=True, help="Path to the gene-set file.")
    parser.add_argument("--fileType", default="gmt", choices=["gmt", "tsv", "json"])
    parser.add_argument("--output", "-o", default="geneSets.json")
    parser.add_argument(
        "--required-score",
        type=int,
        default=400,
        help="STRING confidence floor, 0-1000 (150 low, 400 medium, 700 high, "
             "900 highest). Default 400.",
    )
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument(
        "--network-type", default="functional", choices=["functional", "physical"]
    )
    parser.add_argument(
        "--channels",
        nargs="*",
        default=None,
        help="Restrict evidence to these STRING channels, e.g. "
             "--channels experiments database (avoids expression-derived and "
             "text-mined evidence). Default: STRING's combined score.",
    )
    parser.add_argument("--max-workers", type=int, default=8)
    parser.add_argument("--requests-per-second", type=float, default=8.0)
    args = parser.parse_args()

    makeJson(
        args.geneSymbols,
        args.fileType,
        args.output,
        required_score=args.required_score,
        species=args.species,
        network_type=args.network_type,
        channels=args.channels,
        max_workers=args.max_workers,
        requests_per_second=args.requests_per_second,
    )
