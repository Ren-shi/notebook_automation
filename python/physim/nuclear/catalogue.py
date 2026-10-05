"""The detector catalogue: named models of real detectors, with their dimensions and segmentation (backlog
item 51).

A setup file names a model and may override any of its fields::

    [[detectors]]
    model = "S3"
    theta = "180 deg"
    distance = "30 mm"

    [[gamma_detectors]]
    model = "clover"
    theta = "135 deg"
    distance = "200 mm"

The models are in ``data/detector_models.toml``; each records where its numbers come from and which of them are
typical values instead of data-sheet values::

    from physim.nuclear import catalogue

    catalogue.names("particle")           # ["S3"]
    s3 = catalogue.model("S3")
    s3.fields["rings"], s3.typical        # 24, ("thickness", "dead_layer", "board")
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional

from . import _toml

FILE = Path(__file__).resolve().parent / "data" / "detector_models.toml"
KINDS = ("particle", "gamma")


@dataclass(frozen=True)
class Model:
    name: str
    #: "particle" for ``[[detectors]]``, "gamma" for ``[[gamma_detectors]]``.
    kind: str
    title: str
    sources: tuple
    #: Fields whose values are typical, not taken from a data sheet.
    typical: tuple
    notes: tuple
    #: Default values of setup-file fields.
    fields: dict
    #: Data-sheet values recorded for reference.
    facts: dict
    #: Dead material in the plane of the face that stops particles: annuli with a name and two radii.
    blocking: tuple


@lru_cache(maxsize=None)
def _models() -> dict:
    raw = _toml.loads(FILE.read_text(encoding="utf-8"))
    out = {}
    for name, m in raw.items():
        if m.get("kind") not in KINDS:
            raise ValueError(f"detector model '{name}': kind must be one of {', '.join(KINDS)}")
        out[name] = Model(name, m["kind"], m["title"], tuple(m.get("sources", ())), tuple(m.get("typical", ())),
                          tuple(m.get("notes", ())), dict(m.get("fields", {})), dict(m.get("facts", {})),
                          tuple(dict(b) for b in m.get("blocking", ())))
    return out


def names(kind: Optional[str] = None) -> list:
    """Names of the models, of one kind ("particle" or "gamma") or of both."""
    return [n for n, m in _models().items() if kind is None or m.kind == kind]


def model(name: str, kind: Optional[str] = None) -> Model:
    """A model by name. Raises ``ValueError`` naming the models that exist."""
    m = _models().get(name)
    if m is None or (kind is not None and m.kind != kind):
        what = {"particle": "particle detector", "gamma": "γ-ray detector", None: "detector"}[kind]
        raise ValueError(f"'{name}' is not a {what} model; the models are {', '.join(names(kind)) or 'none'}")
    return m


def with_defaults(raw: dict, kind: str) -> dict:
    """The fields of one detector of a setup file, with its model's defaults filled in where the file gives none.
    Without a ``model`` the fields are returned as they are."""
    if not isinstance(raw, dict) or "model" not in raw:
        return raw
    if not isinstance(raw["model"], str):
        raise ValueError(f"model must be text, one of {', '.join(names(kind))}")
    return {**model(raw["model"], kind).fields, **raw}


__all__ = ["FILE", "KINDS", "Model", "model", "names", "with_defaults"]
