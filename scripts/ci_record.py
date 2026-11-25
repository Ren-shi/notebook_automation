"""A pytest plugin: what of this repository a test file used, for ``scripts/ci_tests.py`` (which loads it).

``PHYSIM_CI_RECORD`` names the JSON file it writes when the session ends: the code that ran ("code": [path, line, name] of
every function, class body and module of the repository that was started, in any thread), the files that were opened
("files"), the directories that were listed ("lists") and the programs that were started ("spawned"; what those used cannot
be seen). It takes ``sys.monitoring`` (Python 3.12), which costs next to nothing: each piece of code reports
once. One test file a process: code that ran for an earlier one, and filled a cache, would be missed for the next.
"""

import json
import os
import sys
from pathlib import Path

ROOT = str(Path(__file__).resolve().parent.parent) + os.sep
OUT = os.environ.get("PHYSIM_CI_RECORD")
TOOL = 4

_code, _files, _lists, _spawned = set(), set(), set(), set()
#: Between the reader and whoever wanted the file or the listing.
_THROUGH = {"pathlib", "os", "glob", "fnmatch", "codecs", "io", "tokenize", "contextlib", "_local", "_abc", "shutil",
            "_collections_abc", "__init__"}
_SPAWN = {"subprocess.Popen", "os.system", "os.posix_spawn", "os.fork", "os.forkpty", "os.exec"}


def _started(code, offset):
    if code.co_filename.startswith(ROOT) and code.co_filename != __file__:
        _code.add((code.co_filename, code.co_firstlineno, code.co_name))
    return sys.monitoring.DISABLE


def _by_machinery() -> bool:
    """Whether the import system, pytest or a traceback read it, rather than the tests or the code they test."""
    frame = sys._getframe(2)
    while frame is not None:
        name = frame.f_code.co_filename.replace("\\", "/")
        if name.startswith("<frozen importlib") or "/importlib/metadata/" in name:
            return True
        stem = name.strip("<>").rpartition("/")[2].replace("frozen ", "").removesuffix(".py")
        if not name.startswith(ROOT.replace("\\", "/")) and "/_pytest/" not in name and stem in _THROUGH:
            frame = frame.f_back
            continue
        return ("/_pytest/" in name or "/pluggy/" in name or name.endswith("/tests/python/conftest.py")
                or stem in ("linecache", "traceback", "warnings"))
    return True


def _audit(event, args):
    if _SPAWN is None:  # the interpreter is closing
        return
    if event in _SPAWN:
        program = args[0] if event in ("subprocess.Popen", "os.posix_spawn", "os.exec") and args else event
        if event == "subprocess.Popen" and not isinstance(args[1], (str, bytes)) and len(args[1]):
            program = program or args[1][0]
        _spawned.add(os.path.basename(os.fsdecode(program)) if isinstance(program, (str, bytes, os.PathLike)) else event)
        return
    if event not in ("open", "os.listdir", "os.scandir", "sqlite3.connect"):
        return
    path = args[0]
    if path is None:
        path = "."
    if isinstance(path, int):
        return
    try:
        path = os.path.abspath(os.fsdecode(path))
    except (TypeError, ValueError):
        return
    if not path.startswith(ROOT) and path + os.sep != ROOT or _by_machinery():
        return
    (_lists if event in ("os.listdir", "os.scandir") else _files).add(path)


def _relative(path: str) -> str:
    return os.path.relpath(path, ROOT).replace(os.sep, "/")


def pytest_unconfigure(config):
    if OUT:
        sys.monitoring.set_events(TOOL, 0)
        Path(OUT).write_text(json.dumps({
            "code": sorted([_relative(f), line, name] for f, line, name in list(_code)),
            "files": sorted(map(_relative, list(_files))), "lists": sorted(map(_relative, list(_lists))),
            "spawned": sorted(_spawned)}), encoding="utf-8")


if OUT:
    sys.monitoring.use_tool_id(TOOL, "physim-ci")
    sys.monitoring.register_callback(TOOL, sys.monitoring.events.PY_START, _started)
    sys.monitoring.set_events(TOOL, sys.monitoring.events.PY_START)
    sys.addaudithook(_audit)
