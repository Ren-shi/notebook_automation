"""GEMINI++ (``Gemini.exe`` from LISE++), out of sight.

Gemini.exe has no silent mode: given a ``.gemini`` file on the command line it opens the file in its window and
waits for Execute. ``gemini_drive.ps1`` (beside this module) presses Execute through Windows UI Automation, waits for
the results, saves the text of every results tab as ``<name>.tab<i>.txt`` and closes the program. Both run on the
runner's hidden desktop, so GEMINI++'s windows never show or take the focus. GEMINI++ itself writes ``<name>.html``
(the summary) and ``<name>.cs4`` next to the input.

GEMINI++ runs from a scratch copy of the LISE++ install (``scratch_install``), so nothing is written into the
install; runs on one machine wait for each other.

    from physim.external import gemini
    exe = gemini.scratch_install("D:/Program Files/LISEcute/LISEcute", Path("work"))
    out = gemini.run(exe, Path("work/runs"), "16O_144Sm", text)
    out.html, out.tabs[1]
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from .runner import PathLike, ProgramFailed, Run
from .runner import run as _run

DRIVER = Path(__file__).with_name("gemini_drive.ps1")


class GeminiFailed(ProgramFailed):
    pass


@dataclass
class GeminiRun:
    html: str
    tabs: List[str] = field(default_factory=list)  # the results window's tabs, in order
    log: str = ""  # the driver's result line
    run: Optional[Run] = None


def scratch_install(lise: PathLike, work: PathLike) -> Path:
    """A copy of the LISE++ folder in ``work`` (made once, reused later); the path of its Gemini.exe."""
    copy = Path(work) / "LISEcute"
    if not (copy / "Gemini.exe").exists():
        shutil.copytree(Path(lise), copy, dirs_exist_ok=True)
    return copy / "Gemini.exe"


def run(exe: PathLike, folder: PathLike, name: str, text: str, tab: str = "Fusion reaction",
        timeout_s: float = 1800, stall_s: float = 45, attempts: int = 3, visible: bool = False) -> GeminiRun:
    """Run GEMINI++ on the input ``text``, saved as ``folder/<name>.gemini``, with the main window's ``tab``
    selected (which calculation runs follows the tab, not the file's mode= key). ``visible`` shows the windows, for
    watching the driver work. A run whose window did not come up or that crashed is tried again, up to
    ``attempts`` times; one that hung is not (GEMINI++'s seed is fixed: the same input fails the same way)."""
    if sys.platform != "win32":
        raise OSError("GEMINI++ from LISE++ runs only on Windows")
    exe, folder = Path(exe).resolve(), Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.glob(f"{name}.*"):
        old.unlink()
    (folder / f"{name}.gemini").write_text(text, encoding="ascii")
    log = folder / f"{name}.driver.log"
    args = ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(DRIVER), "-Exe", str(exe),
            "-InputFile", str(folder / f"{name}.gemini"), "-Tab", tab, "-TimeoutSeconds", str(int(timeout_s)),
            "-StallSeconds", str(int(stall_s)), "-LogFile", str(log)]
    message, r = "", None
    for attempt in range(attempts):
        log.unlink(missing_ok=True)
        if visible:
            code = subprocess.run(args, capture_output=True, text=True).returncode
        else:
            r = _run(args, folder, timeout_s=timeout_s + 200, exclusive="gemini")
            code = r.returncode
        message = log.read_text(encoding="utf-8-sig").strip() if log.exists() else "(no driver log)"
        html = folder / f"{name}.html"
        if code == 0 and html.exists():
            tabs = sorted(folder.glob(f"{name}.tab*.txt"), key=lambda p: int(p.suffixes[-2][4:]))
            return GeminiRun(html=html.read_text(encoding="utf-8", errors="replace"),
                             tabs=[p.read_text(encoding="utf-8-sig") for p in tabs], log=message, run=r)
        if "hung" in message or "exited during the run" in message:
            break
    raise GeminiFailed(f"GEMINI++ run {name} failed after {attempt + 1} attempt(s): {message}", r)
