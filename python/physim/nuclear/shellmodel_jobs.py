"""Shell-model calculations as background jobs, with a cache of their results.

A request goes in, a job ID comes back at once, and the calculation runs in a separate Python process that outlives
the caller (an MCP call or an app click returns straight away)::

    from physim.nuclear import shellmodel_jobs as jobs
    from physim.nuclear.shellmodel import ShellModelRequest

    job = jobs.submit(ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3), (8, 1, 2))))
    jobs.status(job)          # {"state": "running", "step": "kshell 2+ (3 states)", "step_index": 3, ...}
    jobs.wait(job)            # or poll status() until the state is done, failed or cancelled
    result = jobs.result(job) # a ShellModelResult (None while the job is queued or running)

Each job lives in ``~/.physim/shellmodel/<key>/`` (the root can be moved with the environment variable
``PHYSIM_SHELLMODEL``), where ``key`` is :meth:`ShellModelRequest.key` of the request, with the effective charges
and g-factors resolved to the interaction's defaults, and the KSHELL build's ``VERSION`` text as ``build``. The
job ID is that key, so the same request on the same build is never run twice: a second :func:`submit` finds
``results.json`` and is done at once (a cache hit).

The folder holds

- ``request.toml``: the request, the KSHELL folder and build it runs with;
- ``status.json``: ``state`` (queued, running, done, failed, cancelled), the current ``step`` with
  ``step_index`` of ``steps``, ``started``/``updated``/``finished`` times (UTC), a plain ``message``, the memory
  ``estimate``, and on failure the ``log_tail`` (``refused`` is ``"memory"`` when the job was refused for its
  size);
- ``results.json``: :meth:`ShellModelResult.to_dict`, once done;
- ``run/``: KSHELL's inputs, partition files and logs (wave functions are deleted when the job ends, whether it
  succeeded or not); ``worker.out``: what the background process printed.

What the worker does, step by step: estimate the size (:func:`shellmodel_spaces.estimate`) and refuse a request
over the memory limit (half of the RAM free when it starts) unless ``allow_over_limit``; copy the interaction file
and write the partition files; run ``kshell.exe`` for each requested J±, checking the free memory again before each
one (and keeping the Lanczos vectors on disk for that step if memory has become short); run ``transit.exe`` for
every pair of J± sets the requested multipoles connect; read the logs (:func:`shellmodel.read_folder`) into
``results.json``; delete the wave functions. Only one job runs at a time per cache folder (the others wait,
queued), and only one KSHELL program at a time on the machine (``exclusive="kshell"`` in the runner).

:func:`cancel` stops the worker. Its KSHELL program goes with it: the runner keeps every program in a Windows job
object that is closed, and its programs stopped, when the worker process ends (elsewhere the worker stops the
program's process group before exiting).

The cache is kept under a size limit (``PHYSIM_SHELLMODEL_CACHE_MB``, default 256 MB): when a job ends, the jobs
used longest ago are removed until the rest fit (:func:`clean_cache`).
"""

from __future__ import annotations

import contextlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Union

from ..external import kshell as _kshell
from ..external.runner import ProgramFailed
from . import shellmodel_spaces as _spaces
from .names import SYMBOLS, parse_nuclide
from .shellmodel import KINDS, ShellModelRequest, ShellModelResult, label, read_folder

__all__ = ["submit", "status", "result", "cancel", "list_jobs", "wait", "clean_cache", "run_job", "job_folder",
           "cache_root", "resolved", "job_key", "JobFailed", "QUEUED", "RUNNING", "DONE", "FAILED", "CANCELLED"]

QUEUED, RUNNING, DONE, FAILED, CANCELLED = "queued", "running", "done", "failed", "cancelled"
ACTIVE = (QUEUED, RUNNING)
FINISHED = (DONE, FAILED, CANCELLED)

#: Default size limit of the cache folder (all jobs together), in MB.
DEFAULT_CACHE_MB = 256.0
#: Lines of a failed run's log kept in ``status.json``.
LOG_TAIL_LINES = 40
#: KSHELL's iteration limit, as its helper writes it.
MAXITER = 300
#: A queued or running job whose worker is gone, and whose status has not changed for this long, has failed.
STALE_AFTER_S = 15.0
#: How long a rerun waits for the process of the job's earlier, finished run to exit.
WORKER_EXIT_S = 30.0

_KEY = re.compile(r"^[0-9a-f]{16}$")


class JobFailed(RuntimeError):
    """The job failed or was cancelled; ``status`` is its ``status.json``."""

    def __init__(self, message: str, status: Optional[dict] = None):
        super().__init__(message)
        self.status = status or {}


class _Cancelled(BaseException):
    """Raised inside the worker when it is asked to stop (a SIGTERM, or the ``cancel`` marker file)."""


# ---------------------------------------------------------------------------------------------------------------
# Where things are


def cache_root() -> Path:
    """``PHYSIM_SHELLMODEL`` if set, else ``~/.physim/shellmodel``."""
    env = os.environ.get("PHYSIM_SHELLMODEL", "").strip()
    return Path(env).expanduser() if env else Path.home() / ".physim" / "shellmodel"


def _check_key(job_id: str) -> str:
    if not isinstance(job_id, str) or not _KEY.match(job_id):
        raise KeyError(f"not a shell-model job ID: {job_id!r} (16 hexadecimal digits)")
    return job_id


def job_folder(job_id: str) -> Path:
    """The folder of job ``job_id`` (whether it exists or not)."""
    return cache_root() / _check_key(job_id)


def resolved(request: ShellModelRequest) -> ShellModelRequest:
    """The request with its effective charges and g-factors set to the interaction's defaults where it leaves them
    out, so that asking for the defaults and leaving them out are the same job."""
    spec = _spaces.interaction(request.interaction)
    return replace(request, eff_charge=request.eff_charge or spec.eff_charge, gl=request.gl or spec.gl,
                   gs=request.gs or spec.gs)


def job_key(request: ShellModelRequest, build: str) -> str:
    """The job ID of ``request`` run with the KSHELL build ``build`` (its ``VERSION`` text)."""
    return resolved(request).key(build)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _gb(n: float) -> str:
    gb = n / 1024**3
    return f"{gb:.2g} GB" if gb < 10 else f"{gb:.0f} GB"


# ---------------------------------------------------------------------------------------------------------------
# Files: written whole (so a reader never sees half a file), read with a few retries (Windows refuses to replace a
# file another process has open for a moment)


def _write_json(path: Path, data: dict) -> None:
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(data, indent=1, ensure_ascii=False), encoding="utf-8")
    for attempt in range(50):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.02 * (attempt + 1))
    os.replace(tmp, path)


def _read_json(path: Path) -> Optional[dict]:
    for attempt in range(20):
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except FileNotFoundError:
            return None
        except (PermissionError, json.JSONDecodeError):
            time.sleep(0.02 * (attempt + 1))
    return None


def _toml_value(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v)  # a JSON string with \uXXXX escapes is a valid TOML basic string
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    raise TypeError(f"cannot write {type(v).__name__} to request.toml")


def _write_request(folder: Path, request: ShellModelRequest, install: _kshell.KshellInstall, key: str) -> None:
    lines = ["# A physim shell-model job (physim.nuclear.shellmodel_jobs). The worker reads this file.",
             "", "[request]"]
    for name, value in request.to_dict().items():
        if value is not None:
            lines.append(f"{name} = {_toml_value(value)}")
    lines += ["", "[kshell]", f"folder = {_toml_value(Path(install.folder).resolve().as_posix())}",
              f"build = {_toml_value(install.version)}", "", "[job]", f"key = {_toml_value(key)}",
              f"submitted = {_toml_value(_now())}"]
    (folder / "request.toml").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_request(folder: Path) -> dict:
    try:
        import tomllib as toml
    except ImportError:  # Python < 3.11
        import tomli as toml  # type: ignore[no-redef]
    return toml.loads((folder / "request.toml").read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------------------------------------------
# Processes


def _process_created(pid: int) -> Optional[float]:
    """When process ``pid`` started (an opaque number that tells a reused PID apart), or None if it is not
    running."""
    if not pid or pid <= 0:
        return None
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes as w
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.OpenProcess.restype = w.HANDLE
        k32.OpenProcess.argtypes = [w.DWORD, w.BOOL, w.DWORD]
        k32.GetExitCodeProcess.argtypes = [w.HANDLE, ctypes.POINTER(w.DWORD)]
        k32.GetProcessTimes.argtypes = [w.HANDLE] + [ctypes.POINTER(w.FILETIME)] * 4
        k32.CloseHandle.argtypes = [w.HANDLE]
        h = k32.OpenProcess(0x1000, False, int(pid))  # PROCESS_QUERY_LIMITED_INFORMATION
        if not h:
            return None
        try:
            code = w.DWORD()
            if not k32.GetExitCodeProcess(h, ctypes.byref(code)) or code.value != 259:  # STILL_ACTIVE
                return None
            times = [w.FILETIME() for _ in range(4)]
            if not k32.GetProcessTimes(h, *[ctypes.byref(t) for t in times]):
                return 0.0
            return float((times[0].dwHighDateTime << 32) | times[0].dwLowDateTime)
        finally:
            k32.CloseHandle(h)
    for child in list(_children):  # reap our own finished workers, so they are not taken for running zombies
        if child.poll() is not None:
            _children.remove(child)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return None
    except PermissionError:
        pass
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
        fields = stat.rsplit(")", 1)[1].split()
        if fields[0] == "Z":
            return None
        return float(fields[19])
    except (OSError, IndexError, ValueError):
        return 0.0


def _alive(pid: Optional[int], created: Optional[float]) -> bool:
    if not pid:
        return False
    now = _process_created(int(pid))
    if now is None:
        return False
    return created is None or not now or not created or abs(now - created) < 1.0


def _terminate(pid: Optional[int], created: Optional[float], wait_s: float = 15.0) -> None:
    """Stop process ``pid`` (if it is still the process that started at ``created``) and wait for it to end."""
    if not _alive(pid, created) or int(pid) == os.getpid():  # a job run in this process is not stopped this way
        return
    pid = int(pid)
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes as w
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.OpenProcess.restype = w.HANDLE
        k32.OpenProcess.argtypes = [w.DWORD, w.BOOL, w.DWORD]
        k32.TerminateProcess.argtypes = [w.HANDLE, w.UINT]
        k32.WaitForSingleObject.argtypes = [w.HANDLE, w.DWORD]
        k32.CloseHandle.argtypes = [w.HANDLE]
        h = k32.OpenProcess(0x0001 | 0x00100000 | 0x1000, False, pid)  # TERMINATE | SYNCHRONIZE | QUERY_LIMITED
        if not h:
            return
        try:
            k32.TerminateProcess(h, 1)
            k32.WaitForSingleObject(h, int(wait_s * 1000))
        finally:
            k32.CloseHandle(h)
        return
    import signal
    with contextlib.suppress(ProcessLookupError):
        os.kill(pid, signal.SIGTERM)  # the worker stops its KSHELL program's process group, then exits
    deadline = time.monotonic() + wait_s
    while time.monotonic() < deadline and _alive(pid, None):
        time.sleep(0.1)
    if _alive(pid, None):
        with contextlib.suppress(ProcessLookupError):
            os.kill(pid, signal.SIGKILL)


_children: List[subprocess.Popen] = []


def _spawn(folder: Path) -> subprocess.Popen:
    """Start the worker for ``folder`` as a separate Python process, out of sight, that outlives this one."""
    env = dict(os.environ)
    package_parent = str(Path(__file__).resolve().parent.parent.parent)  # the folder holding physim/
    env["PYTHONPATH"] = os.pathsep.join([package_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep)
                                                            if p and p != package_parent])
    env["PYTHONIOENCODING"] = "utf-8"
    args = [sys.executable, "-m", "physim.nuclear.shellmodel_jobs", str(folder)]
    out = open(folder / "worker.out", "wb")
    try:
        common = dict(cwd=str(folder), env=env, stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT,
                      close_fds=True)
        if sys.platform == "win32":
            # No console window; its own process group (a Ctrl+C in the caller's console does not reach it); out of
            # the caller's job object where allowed, so that closing the caller does not stop the job.
            flags = 0x08000000 | 0x00000200  # CREATE_NO_WINDOW | CREATE_NEW_PROCESS_GROUP
            try:
                p = subprocess.Popen(args, creationflags=flags | 0x01000000, **common)  # CREATE_BREAKAWAY_FROM_JOB
            except OSError:  # the caller's job forbids breaking away
                p = subprocess.Popen(args, creationflags=flags, **common)
        else:
            p = subprocess.Popen(args, start_new_session=True, **common)
    finally:
        out.close()
    _children.append(p)
    return p


# ---------------------------------------------------------------------------------------------------------------
# Status


def _job_status(folder: Path) -> Optional[dict]:
    """``status.json`` of the job in ``folder``, corrected when it says queued or running but the worker is gone."""
    st = _read_json(folder / "status.json")
    has_results = (folder / "results.json").is_file()
    if st is None:
        if not has_results:
            return None
        st = {"job": folder.name, "state": DONE, "message": "Done (results from the cache)."}
    if st.get("state") in ACTIVE:
        if has_results:
            st.update(state=DONE, message="Done.")
            return st
        worker = _read_json(folder / "worker.json") or {}
        alive = _alive(st.get("pid"), st.get("pid_created")) if st.get("pid") else False
        alive = alive or _alive(worker.get("pid"), worker.get("created"))
        if not alive:
            # Give a worker that has just been started the chance to write its own status.
            last = float(worker.get("spawned_at") or 0.0)
            with contextlib.suppress(TypeError, ValueError):
                last = max(last, datetime.fromisoformat(st.get("updated")).timestamp())
            if time.time() - last > STALE_AFTER_S:
                tail = _tail_of(folder / "worker.out")
                st.update(state=FAILED, updated=_now(), finished=_now(),
                          message="The background process stopped without finishing (it may have been closed or "
                                  "the computer restarted). Submit the request again to rerun it.",
                          log_tail=tail or st.get("log_tail", ""))
                _write_json(folder / "status.json", st)
    return st


def _tail_of(path: Path, lines: int = LOG_TAIL_LINES) -> str:
    try:
        text = path.read_bytes().decode("utf-8", errors="replace")
    except OSError:
        return ""
    return "\n".join(text.rstrip().splitlines()[-lines:])


def status(job_id: str) -> dict:
    """The job's ``status.json`` (see the module's description): ``state``, ``step``, ``step_index``, ``steps``,
    ``message``, times, ``estimate``, and ``log_tail`` on failure. KeyError for a job that does not exist."""
    folder = job_folder(job_id)
    st = _job_status(folder)
    if st is None:
        raise KeyError(f"no shell-model job {job_id} in {cache_root()}")
    st.setdefault("job", job_id)
    st["folder"] = str(folder)
    return st


def result(job_id: str) -> Optional[ShellModelResult]:
    """The result of a finished job, None while it is queued or running. Raises :class:`JobFailed` for a failed
    or cancelled job and KeyError for one that does not exist."""
    folder = job_folder(job_id)
    data = _read_json(folder / "results.json")
    if data is not None:
        return ShellModelResult.from_dict(data)
    st = status(job_id)
    if st.get("state") in ACTIVE:
        return None
    raise JobFailed(f"shell-model job {job_id} {st.get('state')}: {st.get('message', '')}", st)


def wait(job_id: str, timeout_s: Optional[float] = None, poll_s: float = 0.5) -> dict:
    """Wait until the job is done, failed or cancelled (at most ``timeout_s``); returns its status.
    TimeoutError if it is still queued or running then."""
    deadline = None if timeout_s is None else time.monotonic() + timeout_s
    while True:
        st = status(job_id)
        if st.get("state") in FINISHED:
            return st
        if deadline is not None and time.monotonic() > deadline:
            raise TimeoutError(f"shell-model job {job_id} is still {st.get('state')} after {timeout_s} s "
                               f"({st.get('step', '')})")
        time.sleep(poll_s)


def list_jobs() -> List[dict]:
    """Every job in the cache, most recently updated first: its status with ``job``, ``state``, ``message`` and
    the request (``request``)."""
    root = cache_root()
    if not root.is_dir():
        return []
    out = []
    for folder in root.iterdir():
        if not folder.is_dir() or not _KEY.match(folder.name):
            continue
        st = _job_status(folder)
        if st is None:
            continue
        st.setdefault("job", folder.name)
        st["folder"] = str(folder)
        if "request" not in st:
            with contextlib.suppress(Exception):
                st["request"] = _read_request(folder).get("request")
        out.append(st)
    out.sort(key=lambda s: s.get("updated") or "", reverse=True)
    return out


# ---------------------------------------------------------------------------------------------------------------
# Submit and cancel


def _wait_for_workers(folder: Path, st: dict, timeout_s: float = WORKER_EXIT_S) -> None:
    """Wait until the processes of an earlier run of the job in ``folder`` have ended.

    A worker writes its final status just before it exits, so for a moment after a job has failed its process
    still holds ``worker.out`` (Windows then refuses to delete the file) and may still be tidying up. One that
    has not ended after ``timeout_s`` is stopped: its job has finished, so nothing is lost."""
    worker = _read_json(folder / "worker.json") or {}
    processes = [(st.get("pid"), st.get("pid_created")), (worker.get("pid"), worker.get("created"))]
    deadline = time.monotonic() + timeout_s
    while any(_alive(pid, created) and int(pid) != os.getpid() for pid, created in processes):
        if time.monotonic() > deadline:
            for pid, created in processes:
                _terminate(pid, created)
            break
        time.sleep(0.05)


def _remove(path: Path, attempts: int = 40) -> None:
    """Delete ``path``, waiting a little while another process (an antivirus scan, an indexer) holds it."""
    for attempt in range(attempts):
        try:
            path.unlink()
            return
        except FileNotFoundError:
            return
        except PermissionError:
            if attempt == attempts - 1:
                raise
            time.sleep(0.05)


def _clear(folder: Path, st: dict) -> None:
    """Remove what an earlier, unfinished run left in the job folder, once its processes have ended."""
    _wait_for_workers(folder, st)
    for name in ("results.json", "status.json", "worker.json", "worker.out", "cancel"):
        _remove(folder / name)
    run = folder / "run"
    if run.exists():
        shutil.rmtree(run, ignore_errors=True)


def submit(request: ShellModelRequest, *, background: bool = True,
           install: Optional[_kshell.KshellInstall] = None) -> str:
    """Start ``request`` and return its job ID at once.

    A request whose results are in the cache is done already (a cache hit); one that is queued or running is not
    started twice; one that failed or was cancelled is run again. ``background=False`` runs the job in this
    process and returns when it ends (for scripts and tests). ``install`` is the KSHELL build to use (default
    :func:`physim.external.kshell.find_install`). FileNotFoundError when there is no KSHELL; ValueError when the
    request does not fit its interaction's space."""
    install = install or _kshell.find_install()
    if install is None:
        raise FileNotFoundError("KSHELL was not found: set PHYSIM_KSHELL to the folder holding kshell.exe and "
                                "transit.exe, or install physim with its bundled copy")
    request = resolved(request)
    _spaces.nuclide_valence(request.nuclide, request.interaction)  # ValueError now, not in the background
    key = request.key(install.version)
    folder = cache_root() / key
    if (folder / "results.json").is_file():
        with contextlib.suppress(OSError):
            os.utime(folder / "results.json")  # most recently used: the last the clean-up removes
        return key
    if folder.is_dir():
        st = _job_status(folder)
        if st is not None and st.get("state") in ACTIVE:
            return key
        _clear(folder, st or {})
    folder.mkdir(parents=True, exist_ok=True)
    _write_request(folder, request, install, key)
    _write_json(folder / "status.json", {
        "job": key, "state": QUEUED, "step": "waiting to start", "step_index": 0, "steps": 0,
        "message": "Queued.", "submitted": _now(), "started": None, "updated": _now(), "finished": None,
        "request": request.to_dict(), "build": install.version.splitlines()[0] if install.version else ""})
    if not background:
        run_job(folder)
        return key
    p = _spawn(folder)
    _write_json(folder / "worker.json", {"pid": p.pid, "created": _process_created(p.pid),
                                         "spawned_at": time.time()})
    return key


def cancel(job_id: str) -> dict:
    """Stop a queued or running job (the worker and its KSHELL program) and mark it cancelled; returns its
    status. A finished job is left as it is."""
    folder = job_folder(job_id)
    st = status(job_id)
    if st.get("state") not in ACTIVE:
        return st
    (folder / "cancel").write_text(_now(), encoding="utf-8")  # a worker not stopped below stops at its next step
    worker = _read_json(folder / "worker.json") or {}
    _terminate(st.get("pid"), st.get("pid_created"))
    _terminate(worker.get("pid"), worker.get("created"))
    st = _read_json(folder / "status.json") or st
    if st.get("state") in ACTIVE:
        st.update(state=CANCELLED, message="Cancelled.", updated=_now(), finished=_now())
        _write_json(folder / "status.json", st)
    for attempt in range(20):  # the KSHELL program may take a moment to let go of its files
        _kshell.remove_wavefunctions(folder / "run")
        if not list((folder / "run").glob("*.wav")):
            break
        time.sleep(0.25)
    st["folder"] = str(folder)
    return st


# ---------------------------------------------------------------------------------------------------------------
# The cache's size


def _folder_bytes(folder: Path) -> int:
    total = 0
    for p in folder.rglob("*"):
        with contextlib.suppress(OSError):
            if p.is_file():
                total += p.stat().st_size
    return total


def _last_used(folder: Path) -> float:
    times = []
    for name in ("results.json", "status.json", "request.toml"):
        with contextlib.suppress(OSError):
            times.append((folder / name).stat().st_mtime)
    return max(times) if times else 0.0


def cache_limit_bytes() -> int:
    """``PHYSIM_SHELLMODEL_CACHE_MB`` (default 256 MB) in bytes."""
    try:
        mb = float(os.environ.get("PHYSIM_SHELLMODEL_CACHE_MB", "") or DEFAULT_CACHE_MB)
    except ValueError:
        mb = DEFAULT_CACHE_MB
    return int(mb * 1024 * 1024)


def clean_cache(limit_bytes: Optional[int] = None, keep: tuple = ()) -> List[str]:
    """Remove finished jobs, the one used longest ago first, until the cache is under ``limit_bytes`` (default
    :func:`cache_limit_bytes`). Queued and running jobs and the jobs in ``keep`` stay. Returns the removed IDs."""
    limit = cache_limit_bytes() if limit_bytes is None else int(limit_bytes)
    root = cache_root()
    if not root.is_dir():
        return []
    folders = [f for f in root.iterdir() if f.is_dir() and _KEY.match(f.name)]
    sizes = {f: _folder_bytes(f) for f in folders}
    total = sum(sizes.values())
    removed = []
    for folder in sorted(folders, key=_last_used):
        if total <= limit:
            break
        if folder.name in keep:
            continue
        st = _job_status(folder)
        if st is not None and st.get("state") in ACTIVE:
            continue
        shutil.rmtree(folder, ignore_errors=True)
        if not folder.exists():
            total -= sizes[folder]
            removed.append(folder.name)
    return removed


# ---------------------------------------------------------------------------------------------------------------
# The worker


@contextlib.contextmanager
def _job_lock(root: Path, waiting) -> Iterator[None]:
    """One job at a time per cache folder; ``waiting()`` is called once if another job holds the lock."""
    root.mkdir(parents=True, exist_ok=True)
    with open(root / ".worker.lock", "a+b") as f:
        told = False
        while True:
            try:
                if sys.platform == "win32":
                    import msvcrt
                    f.seek(0)
                    msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError:
                if not told:
                    waiting()
                    told = True
                time.sleep(0.5)
        try:
            yield
        finally:
            if sys.platform == "win32":
                import msvcrt
                f.seek(0)
                with contextlib.suppress(OSError):
                    msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def _parity_letter(parity: int) -> str:
    return "p" if parity > 0 else "n"


def _pairs(sets: List[tuple], kinds: tuple) -> List[tuple]:
    """The pairs of J± sets (left, right) that transit.exe runs for the multipoles ``kinds``, as KSHELL's helper
    chooses them: same parity and |ΔJ| ≤ 2 for E2 (|ΔJ| ≤ 1 for M1 alone), opposite parity and |ΔJ| ≤ 1 for E1;
    each set with itself too (moments and transitions within it), except 0 with 0."""
    reach_same = 4 if "E2" in kinds else 2 if "M1" in kinds else -1
    reach_other = 2 if "E1" in kinds else -1
    out = []
    for i, (j2a, pa) in enumerate(sets):
        for j2b, pb in sets[i:]:
            if j2a == 0 and j2b == 0:
                continue
            reach = reach_same if pa == pb else reach_other
            if abs(j2a - j2b) <= reach:
                out.append(((j2a, pa), (j2b, pb)))
    return out


class _Worker:
    def __init__(self, folder: Path):
        self.folder = folder
        self.run_dir = folder / "run"
        self.status: Dict[str, Any] = _read_json(folder / "status.json") or {"job": folder.name}
        self.status.update(pid=os.getpid(), pid_created=_process_created(os.getpid()))

    def write(self, **fields) -> None:
        if (self.folder / "cancel").exists():
            raise _Cancelled()
        self.status.update(fields, updated=_now())
        _write_json(self.folder / "status.json", self.status)

    def finish(self, state: str, message: str, **fields) -> None:
        self.status.update(fields, state=state, message=message, updated=_now(), finished=_now())
        _write_json(self.folder / "status.json", self.status)

    def outcome(self, error: BaseException) -> tuple:
        """(state, message, extra status fields) for a job that ended with ``error``."""
        step = self.status.get("step", "a step")
        if isinstance(error, _Cancelled):
            return CANCELLED, "Cancelled.", {}
        if isinstance(error, _Refused):
            return FAILED, str(error), {"refused": "memory", "log_tail": ""}
        if isinstance(error, ProgramFailed):
            return FAILED, f"{step} failed: {str(error).splitlines()[0]}", {
                "log_tail": "\n".join(str(error).splitlines()[-LOG_TAIL_LINES:])}
        tb = "".join(traceback.format_exception(type(error), error, error.__traceback__))  # the traceback is the log
        return FAILED, f"{step} failed: {type(error).__name__}: {error}", {
            "log_tail": "\n".join(tb.rstrip().splitlines()[-LOG_TAIL_LINES:])}

    def run(self) -> None:
        """Run the job. The final status is written last, after the wave functions are gone and while the job lock
        is still held, so that once a reader sees the job finished, the worker no longer touches its folder (it
        only tidies the rest of the cache and exits)."""
        outcome = None
        try:
            self.write()
            with _job_lock(self.folder.parent, lambda: self.write(
                    state=QUEUED, step="waiting for another shell-model job", message="Queued behind another job.")):
                try:
                    self.write(state=RUNNING, started=_now(), step="estimating the size", message="Running.")
                    outcome = self.calculate()
                except (_Cancelled, Exception) as e:
                    outcome = self.outcome(e)
                finally:
                    if self.run_dir.is_dir():
                        _kshell.remove_wavefunctions(self.run_dir)
                    if outcome is not None:
                        state, message, fields = outcome
                        self.finish(state, message, **fields)
        except (_Cancelled, Exception) as e:  # while waiting for the lock
            if outcome is None:
                state, message, fields = self.outcome(e)
                self.finish(state, message, **fields)
        if self.status.get("state") in FINISHED:
            with contextlib.suppress(Exception):
                clean_cache(keep=(self.folder.name,))

    def calculate(self) -> tuple:
        spec_text = _read_request(self.folder)
        request = ShellModelRequest.from_dict(spec_text["request"])
        install = _kshell.install_from(spec_text["kshell"]["folder"])
        if install is None:
            raise FileNotFoundError(f"KSHELL is no longer in {spec_text['kshell']['folder']}")
        if install.version != spec_text["kshell"].get("build", install.version):
            raise RuntimeError(f"the KSHELL build in {install.folder} changed since the job was submitted; "
                               f"submit the request again")
        spec = _spaces.interaction(request.interaction)
        z_val, n_val = _spaces.nuclide_valence(request.nuclide, spec)

        est = _spaces.estimate(request, free_bytes=_spaces.available_memory_bytes())
        self.status["estimate"] = est.to_dict()
        if not est.fits and not request.allow_over_limit:
            raise _Refused(f"Refused: {est.reason} (allow_over_limit runs it anyway.)")

        notes: List[str] = []
        if not est.fits:
            notes.append(f"Run over the memory limit at the user's request: {est.reason}")
        sets, counts = [], {}
        for j2, parity, n in request.states:
            jdim = _spaces.j_dimension(spec.space, z_val, n_val, parity, j2)
            if jdim == 0:
                notes.append(f"No {label(j2, parity)} state exists for {request.nuclide} in the {spec.space.name} "
                             f"space; skipped.")
                continue
            if jdim < n:
                notes.append(f"Only {jdim} {label(j2, parity)} state{'s' if jdim > 1 else ''} exist; "
                             f"{n} were asked for.")
            sets.append((j2, parity))
            counts[(j2, parity)] = min(n, jdim)
        if not sets:
            raise ValueError("none of the requested states exists in this space")
        kinds = tuple(k for k in request.transitions if k in KINDS)
        pairs = _pairs(sets, kinds) if kinds else []
        steps = 1 + len(sets) + len(pairs) + 1
        self.write(steps=steps, step_index=1, step="writing the inputs")

        z, a = parse_nuclide(request.nuclide)
        base = f"{SYMBOLS[z]}{a}_{spec.name}"
        self.run_dir.mkdir(parents=True, exist_ok=True)
        if install.snt_dir is None or not (install.snt_dir / spec.snt).is_file():
            raise FileNotFoundError(f"{spec.snt} is not in {install.snt_dir or install.folder}")
        shutil.copyfile(install.snt_dir / spec.snt, self.run_dir / spec.snt)
        ptn = {}
        for parity in sorted({p for _, p in sets}, reverse=True):
            ptn[parity] = f"{base}_{_parity_letter(parity)}.ptn"
            (self.run_dir / ptn[parity]).write_text(
                _spaces.partitions(spec.space, z_val, n_val, parity, snt_name=spec.snt), encoding="ascii")

        waves = {}
        for i, (j2, parity) in enumerate(sets):
            n = counts[(j2, parity)]
            what = f"{label(j2, parity)} ({n} state{'s' if n > 1 else ''})"
            mode = self.check_memory(request, j2, parity, n, what)
            self.write(step_index=2 + i, step=f"kshell {what}")
            name = f"{base}_j{j2}{_parity_letter(parity)}"
            waves[(j2, parity)] = f"{name}.wav"
            n_restart, max_lanc = _spaces.lanczos_parameters(n)
            text = _kshell.kshell_input(
                beta_cm=float(spec.beta_cm), eff_charge=tuple(request.eff_charge), fn_int=spec.snt,
                fn_ptn=ptn[parity], fn_save_wave=waves[(j2, parity)], gl=tuple(request.gl), gs=tuple(request.gs),
                hw_type=int(spec.hw_type), is_double_j=True, max_lanc_vec=max_lanc, maxiter=MAXITER,
                mode_lv_hdd=mode, mtot=j2, n_block=0, n_eigen=n, n_restart_vec=n_restart)
            _kshell.run_kshell(install, self.run_dir, text, name=name, threads=request.threads)

        for i, (left, right) in enumerate(pairs):
            what = f"{label(*left)} × {label(*right)}"
            self.write(step_index=2 + len(sets) + i, step=f"transit {what}")
            name = f"{base}_tr_j{left[0]}{_parity_letter(left[1])}_j{right[0]}{_parity_letter(right[1])}"
            text = _kshell.transit_input(
                fn_int=spec.snt, fn_ptn_l=ptn[left[1]], fn_ptn_r=ptn[right[1]], fn_load_wave_l=waves[left],
                fn_load_wave_r=waves[right], hw_type=int(spec.hw_type), eff_charge=tuple(request.eff_charge),
                gl=tuple(request.gl), gs=tuple(request.gs))
            _kshell.run_transit(install, self.run_dir, text, name=name, threads=request.threads)

        self.write(step_index=steps, step="reading the results")
        res = read_folder(self.run_dir, request, kshell_version=install.version)
        res.transitions = [t for t in res.transitions if t.kind in kinds]
        res.notes = notes + list(res.notes)
        payload = res.to_dict()
        payload["job"] = self.folder.name
        _kshell.remove_wavefunctions(self.run_dir)
        with contextlib.suppress(OSError):
            (self.run_dir / spec.snt).unlink()  # the install has it; the logs name it
        _write_json(self.folder / "results.json", payload)  # last: results.json alone means done (a cache hit)
        return DONE, f"Done: {len(res.states)} states, {len(res.transitions)} transitions.", {
            "step": "done", "step_index": steps}

    def check_memory(self, request: ShellModelRequest, j2: int, parity: int, n: int, what: str) -> int:
        """Check the free memory again before a kshell.exe run; returns ``mode_lv_hdd`` for it (0, vectors in
        memory, when that fits the limit now, else 1). Raises :class:`_Refused` when the run no longer fits even
        with the vectors on disk, unless the request allows going over the limit."""
        free = _spaces.available_memory_bytes()
        one = _spaces.estimate(replace(request, states=((j2, parity, n),)), free_bytes=free)
        self.status.setdefault("memory_checks", []).append(
            {"step": f"kshell {what}", "free_bytes": one.free_bytes, "limit_bytes": one.limit_bytes,
             "needed_bytes": one.memory_bytes, "mode_lv_hdd": one.mode_lv_hdd, "fits": one.fits})
        if one.fits:
            return one.mode_lv_hdd
        if request.allow_over_limit:
            return 1
        raise _Refused(f"Stopped before kshell {what}: the free memory fell to {_gb(one.free_bytes)}, so the limit "
                       f"is now {_gb(one.limit_bytes)} (half of it), and this step needs about "
                       f"{_gb(one.memory_bytes_on_disk_mode)} even with the Lanczos vectors on disk. Close other "
                       f"programs and submit the request again (allow_over_limit runs it anyway).")


class _Refused(Exception):
    """The job does not fit the memory limit."""


def run_job(folder: Union[str, Path]) -> dict:
    """Run the job in ``folder`` (with ``request.toml`` and ``status.json`` written by :func:`submit`) in this
    process; returns its final status. This is what the background process does."""
    folder = Path(folder)
    worker = _Worker(folder)
    worker.run()
    return worker.status


def main(argv: Optional[List[str]] = None) -> int:
    """``python -m physim.nuclear.shellmodel_jobs <job folder>``: the background process."""
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print("usage: python -m physim.nuclear.shellmodel_jobs <job folder>", file=sys.stderr)
        return 2
    if sys.platform != "win32":
        import signal

        def stop(signum, frame):
            raise _Cancelled()
        signal.signal(signal.SIGTERM, stop)
    st = run_job(Path(argv[0]))
    return 0 if st.get("state") == DONE else 1


if __name__ == "__main__":
    sys.exit(main())
