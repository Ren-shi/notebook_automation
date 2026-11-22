"""GOSIA, out of sight.

GOSIA is a console program: it reads its input on standard input and writes its results to ``fort.22`` (and its
working files, ``fort.99`` and others) in the folder it runs in. A run happens in a scratch folder unless one is
given.

    from physim.external import gosia
    out = gosia.run("tools/gosia/gosia.exe", Path("coulex_gosia.inp").read_text())
    print(out.output)          # the text of fort.22
"""

from __future__ import annotations

import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

from .runner import PathLike, ProgramFailed, Run
from .runner import run as _run


@dataclass
class GosiaRun:
    output: str  # fort.22
    files: Dict[str, str] = field(default_factory=dict)  # every file the run wrote, by name (fort.22 among them)
    run: Optional[Run] = None


def run(exe: PathLike, input_text: str, folder: Optional[PathLike] = None, timeout_s: float = 600) -> GosiaRun:
    """Run GOSIA on ``input_text`` in ``folder`` (a scratch folder if None) and return what it wrote."""
    if folder is None:
        with tempfile.TemporaryDirectory(prefix="gosia_") as tmp:
            return run(exe, input_text, tmp, timeout_s)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    before = {p.name: p.stat().st_mtime_ns for p in folder.iterdir() if p.is_file()}
    r = _run([Path(exe).resolve()], folder, timeout_s=timeout_s, stdin_text=input_text)
    files = {p.name: p.read_text(encoding="latin-1") for p in folder.iterdir()
             if p.is_file() and before.get(p.name) != p.stat().st_mtime_ns}
    if not r.ok or "fort.22" not in files:
        raise ProgramFailed("GOSIA failed: " + r.describe(), r)
    return GosiaRun(output=files["fort.22"], files=files, run=r)


def integrated(output: str) -> Dict[int, dict]:
    """OP,INTI's results from ``fort.22``, by experiment (numbered from 1): ``rutherford``, GOSIA's integrated
    Rutherford cross section (mb·mg/cm² over the energy range, per radian of azimuth times the azimuth covered —
    as printed, to four digits; it uses the kinematics of the inelastic scattering, not the elastic), and
    ``yields``, the integrated γ yield of each transition per detector: {detector: {(initial, final): yield}},
    with GOSIA's level numbers (from 1)."""
    res: Dict[int, dict] = {}
    for m in re.finditer(r"INTEGRATED RUTHERFORD CROSS SECTION=\s*([\d.E+-]+)\s+FOR EXP\.\s*(\d+)", output):
        res.setdefault(int(m.group(2)), {"yields": {}})["rutherford"] = float(m.group(1))
    for m in re.finditer(r"INTEGRATED YIELDS\s+EXPERIMENT\s+(\d+)\s+DETECTOR\s+(\d+)(.*?)(?=INTEGRATED|\Z)", output,
                         re.S):
        y = {(int(ni), int(nf)): float(val) for ni, nf, val in re.findall(
            r"^\s+(\d+)\s+(\d+)\s+[\d.]+\s+[\d.]+\s+([\d.E+-]+)\s+[\d.E+-]+\s*$", m.group(3), re.M)}
        res.setdefault(int(m.group(1)), {"yields": {}})["yields"][int(m.group(2))] = y
    return res
