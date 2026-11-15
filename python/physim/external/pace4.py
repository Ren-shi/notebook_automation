"""PACE4 (from LISE++), out of sight.

``PACE4.exe -s deck.pace`` runs without asking anything and writes ``deck.html``, ``deck.rtf`` and ``deck.cs4``
next to the deck. Each run gets a fresh scratch folder. Two PACE4 runs at the same time share a temporary file and
spoil each other's output, so runs on one machine wait for each other.

    from physim.external import pace4
    out = pace4.run("D:/Program Files/LISEcute/LISEcute/PACE4.exe", deck_text)
    out.html                   # the HTML output
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

from .runner import PathLike, ProgramFailed, Run
from .runner import run as _run


@dataclass
class Pace4Run:
    outputs: Dict[str, str] = field(default_factory=dict)  # by extension: "html", "rtf", "cs4"
    run: Optional[Run] = None

    @property
    def html(self) -> str:
        return self.outputs["html"]


def run(exe: PathLike, deck_text: str, timeout_s: float = 1800) -> Pace4Run:
    """Run one PACE4 deck and return its output files' text (read as Latin-1, as PACE4 writes them)."""
    with tempfile.TemporaryDirectory(prefix="pace4_") as tmp:
        work = Path(tmp)
        (work / "run.pace").write_text(deck_text, encoding="ascii", newline="\r\n")
        r = _run([Path(exe).resolve(), "-s", "run.pace"], work, timeout_s=timeout_s, exclusive="pace4")
        outputs = {ext: (work / f"run.{ext}").read_text(encoding="latin-1")
                   for ext in ("html", "rtf", "cs4") if (work / f"run.{ext}").exists()}
        if not r.ok or "html" not in outputs:
            raise ProgramFailed(f"PACE4 failed for the deck\n{deck_text}\n" + r.describe(), r)
        return Pace4Run(outputs=outputs, run=r)
