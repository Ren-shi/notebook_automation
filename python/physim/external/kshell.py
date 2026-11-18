"""KSHELL (the shell-model code of Shimizu et al., Comput. Phys. Commun. 244, 372 (2019)), out of sight.

The low level only: this module writes KSHELL's namelist inputs from the fields it is given, runs ``kshell.exe``
(Lanczos diagonalisation, one 2J± at a time) and ``transit.exe`` (E2/M1 between two sets of states) through
:func:`physim.external.runner.run`, and returns their logs and the files they wrote. It knows nothing about nuclides
or interactions: partition (``.ptn``) and interaction (``.snt``) files are taken as given, and reading numbers out
of the logs is left to :mod:`physim.nuclear.shellmodel`.

    from physim.external import kshell
    install = kshell.find_install()            # PHYSIM_KSHELL, ~/.physim/programs.toml, bundled, <repo>/tools
    text = kshell.kshell_input(fn_int="usdb.snt", fn_ptn="Si28_usdb_p.ptn", fn_save_wave="Si28_usdb_j0p.wav",
                               is_double_j=True, mtot=0, n_eigen=3, eff_charge=(1.5, 0.5))
    out = kshell.run_kshell(install, work, text, name="Si28_usdb_j0p",
                            inputs=[install.snt_dir / "usdb.snt", "Si28_usdb_p.ptn"])
    out.log                                    # what kshell.exe printed (also written to log_Si28_usdb_j0p.txt)
    out.files["Si28_usdb_j0p.wav"]             # the wave functions, for transit.exe
    kshell.remove_wavefunctions(work)          # they can be GB; drop them once the transitions are done

An install folder holds ``kshell.exe``, ``transit.exe``, ``count_dim.exe`` (or the same names without ``.exe``;
also accepted in a ``bin/`` subfolder, as in KSHELL's own distribution), the interaction files in ``snt/``, and a
``VERSION`` file naming the build (``"unknown"`` without it, as for a user's own copy). On Windows the runtime DLLs
(OpenBLAS, the MinGW Fortran and OpenMP runtimes) sit next to the exes, where Windows looks first.

Runs set ``OMP_NUM_THREADS`` (all logical cores unless ``threads`` says otherwise) and ``OMP_STACKSIZE=64m``, and
only one KSHELL program runs at a time on the machine (``exclusive="kshell"``). A run that exits with an error, or
whose log has a Fortran ``STOP``, runtime error or ``error`` line, raises :class:`ProgramFailed` with the log's tail.
"""

from __future__ import annotations

import os
import re
import shutil
from dataclasses import dataclass, field
from numbers import Integral, Real
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Tuple

from . import programs
from .runner import PathLike, ProgramFailed, Run
from .runner import run as _run

__all__ = ["KshellInstall", "KshellRun", "find_install", "install_from", "kshell_input", "transit_input",
           "parse_namelist", "run_kshell", "run_transit", "remove_wavefunctions", "failure_line"]

#: Default time limit of one kshell.exe or transit.exe run.
DEFAULT_TIMEOUT_S = 6 * 3600.0

#: Lines that tell a KSHELL run failed even when the program's exit code is 0 (gfortran's ``STOP 'message'``
#: exits with 0 and prints the message, KSHELL prints ``error [...]`` or ``ERROR: ...`` before giving up).
_FAILURE = re.compile(r"^\s*(?:(?:ERROR\s+)?STOP\b|Fortran runtime error|Error termination|"
                      r"Program aborted|error\b)", re.IGNORECASE | re.MULTILINE)

#: The namelist fields that name files KSHELL reads from its working folder.
_READS = ("fn_int", "fn_ptn", "fn_load_wave", "fn_ptn_l", "fn_ptn_r", "fn_load_wave_l", "fn_load_wave_r")

#: Temporary Lanczos files kshell.exe leaves behind with ``mode_lv_hdd = 1`` (``tmp_lv_<ptn>_<mtot>_...``).
_SCRATCH = ("tmp_lv_*", "tmp_snapshot_*")

_FIELD = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


# --- where KSHELL is -------------------------------------------------------------------------------------------------


@dataclass
class KshellInstall:
    """A KSHELL build: its folder, the three programs, the interaction files and the build's name."""

    folder: Path
    kshell: Path
    transit: Path
    count_dim: Optional[Path]  # None when the build has no count_dim (a user's own copy)
    snt_dir: Optional[Path]  # None when the folder has no snt/
    version: str  # the VERSION file's text, "unknown" without one


def _program(folder: Path, name: str) -> Optional[Path]:
    for where in (folder, folder / "bin"):
        for exe in (f"{name}.exe", name):
            p = where / exe
            if p.is_file():
                return p
    return None


def install_from(folder: Optional[PathLike]) -> Optional[KshellInstall]:
    """The KSHELL build in ``folder``, or None when it has no ``kshell`` and ``transit`` program."""
    if folder is None:
        return None
    folder = Path(folder)
    kshell, transit = _program(folder, "kshell"), _program(folder, "transit")
    if kshell is None or transit is None:
        return None
    snt = folder / "snt"
    version_file = folder / "VERSION"
    version = version_file.read_text(encoding="utf-8", errors="replace").strip() if version_file.is_file() else ""
    return KshellInstall(folder=folder, kshell=kshell, transit=transit, count_dim=_program(folder, "count_dim"),
                         snt_dir=snt if snt.is_dir() else None, version=version or "unknown")


def find_install() -> Optional[KshellInstall]:
    """The KSHELL build physim uses (see :func:`physim.external.programs.find`), or None when there is none."""
    return install_from(programs.find("kshell"))


# --- namelists -------------------------------------------------------------------------------------------------------


def _scalar(v) -> str:
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, Integral):
        return str(int(v))
    if isinstance(v, Real):
        x = float(v)
        if x != x or x in (float("inf"), float("-inf")):
            raise ValueError(f"cannot write {x} to a namelist")
        return repr(x)
    if isinstance(v, (str, os.PathLike)):
        s = os.fspath(v)
        return '"' + s.replace('"', '""') + '"'
    raise TypeError(f"cannot write {type(v).__name__} to a namelist")


def _value(v, trailing_comma: bool) -> str:
    if isinstance(v, (list, tuple)):
        if not v:
            raise ValueError("an empty list cannot be written to a namelist")
        return ", ".join(_scalar(x) for x in v) + (", " if trailing_comma else "")
    return _scalar(v)


def _namelist(fields: Mapping[str, object], trailing_comma: bool) -> str:
    lines = ["&input"]
    for name, v in fields.items():
        if not _FIELD.match(name):
            raise ValueError(f"not a namelist field name: {name!r}")
        if v is None:  # left out: KSHELL's own default applies
            continue
        lines.append(f"  {name} = {_value(v, trailing_comma)}")
    lines.append("&end")
    return "\n".join(lines) + "\n"


def kshell_input(**fields) -> str:
    """kshell.exe's ``&input ... &end`` namelist with ``fields`` in the order given (None leaves a field out).

    Booleans become ``.true.``/``.false.``, strings and paths are quoted, sequences are written as KSHELL's
    ``kshell_ui.py`` writes them (``eff_charge = 1.5, 0.5, ``)."""
    return _namelist(fields, trailing_comma=True)


def transit_input(**fields) -> str:
    """transit.exe's ``&input ... &end`` namelist, as :func:`kshell_input` (sequences without the trailing comma,
    ``eff_charge = 1.5, 0.5``, as ``kshell_ui.py`` writes them for transit.exe)."""
    return _namelist(fields, trailing_comma=False)


def _split_values(text: str) -> list:
    out, cur, quote = [], "", None
    for ch in text:
        if quote:
            cur += ch
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
            cur += ch
        elif ch == ",":
            out.append(cur.strip())
            cur = ""
        elif ch == "!":
            break
        else:
            cur += ch
    out.append(cur.strip())
    while out and out[-1] == "":
        out.pop()
    return out


def _parse_scalar(s: str):
    low = s.lower()
    if low in (".true.", ".t.", "t", "true"):
        return True
    if low in (".false.", ".f.", "f", "false"):
        return False
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        q = s[0]
        return s[1:-1].replace(q + q, q)
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(low.replace("d", "e"))
    except ValueError:
        raise ValueError(f"cannot read namelist value {s!r}") from None


def parse_namelist(text: str) -> Dict[str, object]:
    """The fields of the first ``&input ... &end`` (or ``/``) group of a namelist: name (lower case) → value, a tuple
    for several values. For checking inputs; one ``name = value`` per line, as KSHELL's inputs are written."""
    fields: Dict[str, object] = {}
    inside, last = False, None
    for raw in text.splitlines():
        line = raw.strip()
        if not inside:
            if line.lower().startswith("&input"):
                inside = True
                line = line[len("&input"):].strip()
            else:
                continue
        if line.lower() in ("&end", "/", "$end"):
            break
        if not line or line.startswith("!"):
            continue
        if "=" in line.split("!", 1)[0] and not line.lstrip().startswith(('"', "'")):
            name, rest = line.split("=", 1)
            last = name.strip().lower()
            fields[last] = [_parse_scalar(x) for x in _split_values(rest)]
        elif last is not None:  # a continuation line
            fields[last] = list(fields[last]) + [_parse_scalar(x) for x in _split_values(line)]
    return {k: (v[0] if len(v) == 1 else tuple(v)) for k, v in fields.items()}


# --- running ---------------------------------------------------------------------------------------------------------


@dataclass
class KshellRun:
    """One kshell.exe or transit.exe run: its log, the files it wrote (by name) and how the run went."""

    log: str
    files: Dict[str, Path] = field(default_factory=dict)
    run: Optional[Run] = None


def failure_line(log: str) -> Optional[str]:
    """The first line of ``log`` that tells a KSHELL run failed, or None."""
    m = _FAILURE.search(log)
    if m is None:
        return None
    end = log.find("\n", m.start())
    return log[m.start():end if end >= 0 else None].strip()


def _tail(log: str, lines: int = 30) -> str:
    return "\n".join(log.rstrip().splitlines()[-lines:])


def _snapshot(folder: Path) -> Dict[str, Tuple[int, int]]:
    out = {}
    for p in folder.iterdir():
        if p.is_file():
            st = p.stat()
            out[p.name] = (st.st_mtime_ns, st.st_size)
    return out


def _remove_scratch(folder: Path, names: Iterable[str]) -> None:
    for name in names:
        if any(Path(name).match(pat) for pat in _SCRATCH):
            try:
                (folder / name).unlink()
            except OSError:
                pass


def _environment(threads: Optional[int]) -> Dict[str, str]:
    n = int(threads) if threads else (os.cpu_count() or 1)
    if n < 1:
        raise ValueError(f"threads must be at least 1, not {threads}")
    return {"OMP_NUM_THREADS": str(n), "OMP_STACKSIZE": "64m", "GFORTRAN_UNBUFFERED_PRECONNECTED": "y"}


def _execute(exe: Path, what: str, folder: PathLike, input_text: str, *, name: str, threads: Optional[int],
             timeout_s: float, inputs: Iterable[PathLike], log_name: Optional[str],
             writes: Iterable[str] = ()) -> KshellRun:
    if not _FIELD.match(name.replace("-", "_").replace(".", "_")):
        raise ValueError(f"not a usable run name: {name!r}")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for src in inputs:
        src = Path(src)
        if not src.is_absolute():
            src = folder / src
        dest = folder / src.name
        if not src.is_file():
            raise FileNotFoundError(f"{what} input file not found: {src}")
        if src.resolve() != dest.resolve():
            shutil.copyfile(src, dest)
    fields = parse_namelist(input_text)
    missing = [f"{k} = {fields[k]}" for k in _READS if isinstance(fields.get(k), str)
               and not (folder / str(fields[k])).is_file()]
    if missing:
        raise FileNotFoundError(f"{what} for {name}: not in {folder}: " + ", ".join(missing))

    input_file = folder / f"{name}.input"
    input_file.write_text(input_text, encoding="ascii", newline="\n")
    log_file = folder / (log_name or f"log_{name}.txt")
    before = _snapshot(folder)
    r = _run([exe.resolve(), input_file.name], folder, timeout_s=timeout_s, exclusive="kshell",
             env=_environment(threads))
    log = r.stdout + (("\n" if r.stdout and not r.stdout.endswith("\n") else "") + r.stderr if r.stderr else "")
    log_file.write_text(log, encoding="latin-1", errors="replace", newline="")
    after = _snapshot(folder)
    changed = [n for n, s in after.items() if before.get(n) != s]
    _remove_scratch(folder, changed)
    files = {n: folder / n for n in sorted(changed) if (folder / n).is_file()}
    files[input_file.name] = input_file
    files[log_file.name] = log_file

    reason = None
    if r.timed_out:
        reason = f"stopped after {r.elapsed_s:.0f} s (timeout)"
    elif r.returncode != 0:
        reason = f"exited with code {r.returncode}"
    elif failure_line(log):
        reason = f"stopped: {failure_line(log)}"
    else:
        absent = [w for w in writes if not (folder / w).is_file()]
        if absent:
            reason = "did not write " + ", ".join(absent)
    if reason:
        raise ProgramFailed(f"{what} {reason} ({name}, in {folder})\n{_tail(log)}", r)
    return KshellRun(log=log, files=files, run=r)


def run_kshell(install: KshellInstall, folder: PathLike, input_text: str, *, name: str,
               threads: Optional[int] = None, timeout_s: float = DEFAULT_TIMEOUT_S, inputs: Iterable[PathLike] = (),
               log_name: Optional[str] = None) -> KshellRun:
    """Run kshell.exe in ``folder`` on ``input_text`` (from :func:`kshell_input`), written to ``<name>.input``.

    ``inputs`` are files copied into ``folder`` first (the ``.snt`` and ``.ptn`` files, say; relative paths are
    taken in ``folder``); every file the namelist reads must be in ``folder`` by then. The log is returned and
    written to ``log_<name>.txt`` (or ``log_name``). ``files`` holds the files the run wrote or changed (the wave
    functions ``fn_save_wave`` names, the input and the log); the temporary Lanczos files are removed."""
    fields = parse_namelist(input_text)
    wave = fields.get("fn_save_wave")
    return _execute(install.kshell, "kshell.exe", folder, input_text, name=name, threads=threads,
                    timeout_s=timeout_s, inputs=inputs, log_name=log_name,
                    writes=[wave] if isinstance(wave, str) else [])


def run_transit(install: KshellInstall, folder: PathLike, input_text: str, *, name: str,
                threads: Optional[int] = None, timeout_s: float = DEFAULT_TIMEOUT_S, inputs: Iterable[PathLike] = (),
                log_name: Optional[str] = None) -> KshellRun:
    """Run transit.exe in ``folder`` on ``input_text`` (from :func:`transit_input`), as :func:`run_kshell`: the
    wave-function files ``fn_load_wave_l``/``_r`` name must be in ``folder`` (from earlier kshell.exe runs)."""
    return _execute(install.transit, "transit.exe", folder, input_text, name=name, threads=threads,
                    timeout_s=timeout_s, inputs=inputs, log_name=log_name)


def remove_wavefunctions(folder: PathLike) -> int:
    """Delete the wave-function files (``*.wav``) and any temporary Lanczos files kshell.exe left in ``folder``;
    returns the bytes freed. Logs, inputs and partition files stay."""
    folder = Path(folder)
    if not folder.is_dir():
        return 0
    freed = 0
    for pattern in ("*.wav",) + _SCRATCH:
        for p in folder.glob(pattern):
            if p.is_file():
                size = p.stat().st_size
                try:
                    p.unlink()
                except OSError:
                    continue
                freed += size
    return freed
