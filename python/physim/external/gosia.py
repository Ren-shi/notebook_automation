"""GOSIA, out of sight.

GOSIA is a console program: it reads its input on standard input and writes its results to ``fort.22`` (and its
working files, ``fort.99`` and others) in the folder it runs in. A run happens in a scratch folder unless one is
given.

    from physim.external import gosia
    out = gosia.run("tools/gosia/gosia.exe", Path("coulex_gosia.inp").read_text())
    print(out.output)          # the text of fort.22
"""

from __future__ import annotations

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
