"""Reader for the gene-set network JSON produced by :func:`gsnetact.makeJson`.

Behaviour is unchanged apart from validation: a malformed file now fails at load
with a message naming the offending gene set, rather than at scoring time with
an ``AttributeError`` from somewhere in the network code.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class pjson(dict):
    """``{set: {gene: {partner: score}}}`` loaded from ``jsonFile``.

    Parameters
    ----------
    jsonFile:
        Path to the JSON file, or an already-parsed mapping.
    validate:
        Check the nesting depth and value types on load. Default ``True``.
    """

    def __init__(self, jsonFile: str | Path | dict, validate: bool = True):
        super().__init__()
        if isinstance(jsonFile, dict):
            self.jsonFile = "<dict>"
            self.js: dict[str, Any] = dict(jsonFile)
        else:
            self.jsonFile = str(jsonFile)
            with open(self.jsonFile) as handle:
                self.js = json.load(handle)

        if validate:
            self._validate()

        self.geneNamesList: list[str] = []
        self.geneSayi: list[int] = []
        for nodes in self.js.values():
            self.geneNamesList.extend(nodes.keys())
            self.geneSayi.append(len(nodes))

        self.iter_count = 0
        self.update(self.js)

    def _validate(self) -> None:
        if not isinstance(self.js, dict):
            raise ValueError(f"{self.jsonFile}: top level must be an object of gene sets")
        for name, nodes in self.js.items():
            if not isinstance(nodes, dict):
                raise ValueError(
                    f"{self.jsonFile}: gene set {name!r} must map genes to partner "
                    f"dicts, got {type(nodes).__name__}"
                )
            for gene, partners in nodes.items():
                if not isinstance(partners, dict):
                    raise ValueError(
                        f"{self.jsonFile}: gene set {name!r}, gene {gene!r} must map "
                        f"partners to scores, got {type(partners).__name__}"
                    )
                for partner, score in partners.items():
                    if not isinstance(score, (int, float)):
                        raise ValueError(
                            f"{self.jsonFile}: gene set {name!r}, edge "
                            f"{gene!r}-{partner!r} has non-numeric score {score!r}"
                        )

    def genesets(self):
        """The gene sets themselves, for iteration."""
        return self.values()

    @property
    def getAsDict(self) -> dict[str, Any]:
        return self.js

    @property
    def getGeneNames(self) -> list[str]:
        """Every gene name, with repeats across sets."""
        return self.geneNamesList

    @property
    def getUniqueGeneNames(self) -> list[str]:
        """Every gene name once, in first-appearance order."""
        return list(dict.fromkeys(self.geneNamesList))

    @property
    def getGeneCounts(self) -> list[int]:
        return self.geneSayi

    @property
    def getGeneSetCount(self) -> int:
        return len(self.geneSayi)

    @property
    def getFileInfo(self) -> str:
        return "".join(
            f"GeneSet : {name}, Gene Count : {count} \n"
            for name, count in zip(self.js, self.geneSayi)
        )

    def getFileName(self) -> str:
        return self.jsonFile

    def __repr__(self) -> str:  # pragma: no cover - display only
        return f"pjson({self.jsonFile!r}, gene_sets={len(self.js)})"
