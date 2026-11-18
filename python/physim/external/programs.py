"""Where an outside program lives on this machine.

``find("kshell")`` gives the folder of a program, looking in this order and taking the first folder that exists:

1. the environment variable ``PHYSIM_<NAME>`` (``PHYSIM_KSHELL``);
2. ``~/.physim/programs.toml``, a line ``kshell = "C:/path/to/kshell"`` (written by :func:`set_path`);
3. the copy bundled with the installed planner, ``<install>/tools/<name>/`` (the installer puts the bundled
   Python in ``<install>/python``, so physim itself sits in ``<install>/python/Lib/site-packages/physim``);
4. ``<repo>/tools/<name>/`` when physim runs from a checkout (``<repo>/python/physim``), for development.

A setting that names a file rather than a folder (the program's exe, say) stands for the folder that holds it. A
setting that names nothing that exists is passed over, so a stale setting falls back to the bundled copy.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

try:
    import tomllib as _toml
except ImportError:  # Python < 3.11
    try:
        import tomli as _toml  # type: ignore[no-redef]
    except ImportError:
        _toml = None  # type: ignore[assignment]

PathLike = Union[str, Path]

_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def settings_file() -> Path:
    """``~/.physim/programs.toml``."""
    return Path.home() / ".physim" / "programs.toml"


def _check(name: str) -> str:
    if not _NAME.match(name):
        raise ValueError(f"not a program name: {name!r} (letters, digits and _, starting with a letter)")
    return name.lower()


def _read_settings() -> Dict[str, object]:
    path = settings_file()
    if not path.is_file() or _toml is None:
        return {}
    try:
        return _toml.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):  # a broken settings file is passed over, not fatal
        return {}


def _package_folder() -> Path:
    return Path(__file__).resolve().parent.parent  # .../physim


def install_folder() -> Optional[Path]:
    """The installed planner's folder (``<install>`` with ``<install>/python/Lib/site-packages/physim``), or None
    when physim is not running from the installed planner."""
    site = _package_folder().parent
    lib = site.parent
    python = lib.parent
    if site.name.lower() == "site-packages" and lib.name.lower() == "lib" and python.name.lower() == "python":
        return python.parent
    return None


def repo_folder() -> Optional[Path]:
    """The checkout physim runs from (``<repo>`` with ``<repo>/python/physim`` and ``<repo>/pyproject.toml``),
    or None."""
    python = _package_folder().parent
    if python.name.lower() == "python" and (python.parent / "pyproject.toml").is_file():
        return python.parent
    return None


def candidates(name: str) -> List[Tuple[str, Path]]:
    """Every place :func:`find` looks for ``name``, in order, as (where the setting comes from, folder)."""
    key = _check(name)
    out: List[Tuple[str, Path]] = []
    env = os.environ.get(f"PHYSIM_{key.upper()}", "").strip()
    if env:
        out.append((f"PHYSIM_{key.upper()}", Path(env).expanduser()))
    setting = _read_settings().get(key)
    if isinstance(setting, str) and setting.strip():
        out.append((str(settings_file()), Path(setting.strip()).expanduser()))
    install = install_folder()
    if install is not None:
        out.append(("bundled", install / "tools" / key))
    repo = repo_folder()
    if repo is not None:
        out.append(("repository", repo / "tools" / key))
    return out


def find(name: str) -> Optional[Path]:
    """The folder of the outside program ``name`` (``"kshell"``), or None when it is nowhere to be found."""
    for _, folder in candidates(name):
        if folder.is_file():
            folder = folder.parent
        if folder.is_dir():
            return folder.resolve()
    return None


def set_path(name: str, folder: Optional[PathLike]) -> None:
    """Remember ``folder`` as where ``name`` lives, in ``~/.physim/programs.toml`` (other programs' lines are kept).
    ``None`` forgets the setting, so the bundled copy is used again."""
    key = _check(name)
    settings = {k: v for k, v in _read_settings().items() if isinstance(v, str)}
    if folder is None:
        settings.pop(key, None)
    else:
        settings[key] = Path(folder).expanduser().resolve().as_posix()
    path = settings_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# Folders of outside programs physim runs (physim.external.programs)."]
    lines += [f"{k} = {json.dumps(v)}" for k, v in sorted(settings.items())]
    tmp = path.with_suffix(".toml.tmp")
    tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(tmp, path)
