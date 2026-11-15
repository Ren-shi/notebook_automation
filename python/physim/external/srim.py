"""SRIM 2013, out of sight.

- **SR Module** (``SR Module/SRModule.exe``, a console program): reads ``SR.IN`` and writes the stopping table to
  the output file named in it.
- **TRIM** (``TRIM.exe``, a window program): with ``TRIMAUTO`` = 1 it needs no keyboard or mouse, runs ``TRIM.IN``
  and exits when all ions are done. It still opens its window, which goes to the hidden desktop.

The install folder is never written to: the parts each program needs are copied to a scratch folder the first time
it runs, and the copy is reused for later runs until the ``Srim`` is closed.

    from physim.external.srim import Srim
    with Srim("D:/Projects/srim/SRIM2013") as srim:
        table = srim.stopping(sr_in_text)                  # the SR Module output file's text
        outputs = srim.trim(trim_in_text)                  # files TRIM wrote in "SRIM Outputs", by name
        transmitted = outputs.get("TRANSMIT.txt", "")
"""

from __future__ import annotations

import re
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Optional

from .runner import PathLike, ProgramFailed, Run
from .runner import run as _run


class Srim:
    def __init__(self, install: PathLike):
        self.install = Path(install)
        if not (self.install / "TRIM.exe").exists():
            raise FileNotFoundError(f"no TRIM.exe in {self.install}: give the SRIM 2013 folder")
        self._tmp: Optional[tempfile.TemporaryDirectory] = None
        self._sr: Optional[Path] = None
        self._trim: Optional[Path] = None
        self.last_run: Optional[Run] = None

    # -- scratch copies ----------------------------------------------------------------------------------------------

    def _scratch(self) -> Path:
        if self._tmp is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="srim_")
        return Path(self._tmp.name)

    def _sr_folder(self) -> Path:
        if self._sr is None:
            self._sr = self._scratch() / "SR Module"
            shutil.copytree(self.install / "SR Module", self._sr)
        return self._sr

    def _trim_folder(self) -> Path:
        if self._trim is None:
            work = self._scratch() / "SRIM"
            work.mkdir()
            for name in ("Data", "SRIM Outputs", "SRIM Restore"):
                shutil.copytree(self.install / name, work / name)
            shutil.copyfile(self.install / "TRIM.exe", work / "TRIM.exe")
            auto = (self.install / "TRIMAUTO").read_text(encoding="latin-1").splitlines()
            (work / "TRIMAUTO").write_text("\r\n".join(["1", *auto[1:]]) + "\r\n", encoding="latin-1", newline="")
            self._trim = work
        return self._trim

    # -- runs --------------------------------------------------------------------------------------------------------

    def stopping(self, sr_in_text: str, timeout_s: float = 60) -> str:
        """Run SR Module on the text of an ``SR.IN`` file; the text of the output file it names."""
        work = self._sr_folder()
        lines = sr_in_text.splitlines()
        try:
            name = lines[lines.index("---Output File Name") + 1].strip().strip('"')
        except (ValueError, IndexError):
            raise ValueError("SR.IN has no '---Output File Name' line followed by a file name") from None
        out = work / name
        out.unlink(missing_ok=True)
        (work / "SR.IN").write_text(sr_in_text, encoding="ascii", newline="")
        self.last_run = r = _run([work / "SRModule.exe"], work, timeout_s=timeout_s)
        if not r.ok or not out.exists():
            raise ProgramFailed("SR Module failed: " + r.describe(), r)
        return out.read_text(encoding="latin-1")

    def trim(self, trim_in_text: str, timeout_s: float = 3 * 3600) -> Dict[str, str]:
        """Run TRIM on the text of a ``TRIM.IN`` file (written byte for byte); the files it wrote in
        ``SRIM Outputs``, by name. TRIM's exit code carries no meaning, so a run that wrote nothing is not an
        error here: the caller knows which file it expects (``TRANSMIT.txt`` only when transmitted ions are kept)."""
        work = self._trim_folder()
        outputs = work / "SRIM Outputs"
        before = {p.name: p.stat().st_mtime_ns for p in outputs.iterdir() if p.is_file()}
        for p in outputs.iterdir():  # stale results of an earlier run must not pass for this one's
            if p.is_file() and re.fullmatch(r"(?i)[a-z_]+\.txt", p.name):
                p.unlink()
        (work / "TRIM.IN").write_bytes(trim_in_text.encode("latin-1"))
        self.last_run = r = _run([work / "TRIM.exe"], work, timeout_s=timeout_s)
        if r.timed_out:
            raise ProgramFailed("TRIM failed: " + r.describe(), r)
        return {p.name: p.read_text(encoding="latin-1") for p in outputs.iterdir()
                if p.is_file() and before.get(p.name) != p.stat().st_mtime_ns}

    # -- closing -----------------------------------------------------------------------------------------------------

    def close(self) -> None:
        if self._tmp is not None:
            self._tmp.cleanup()
        self._tmp = self._sr = self._trim = None

    def __enter__(self) -> "Srim":
        return self

    def __exit__(self, *exc) -> None:
        self.close()
