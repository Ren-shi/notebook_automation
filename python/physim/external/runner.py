"""Run an outside program so that nothing of it reaches the user's screen.

On Windows every program starts

- on a separate desktop of this window station (``physim_hidden``) that is never shown: its windows, and those of
  the programs it starts (they inherit the desktop), never appear, take the focus or get a taskbar button. UI
  Automation run on that desktop still sees them, which is how window programs without a silent mode are driven;
- with no console window (``CREATE_NO_WINDOW``), so a console program started from the app (which runs without a
  console) does not flash one;
- in a job object: when the program ends, times out, or this Python process goes away, everything it started is
  stopped with it. A dialog no one can see (an error box on the hidden desktop) cannot leave a program hanging
  past the timeout.

Elsewhere the program runs in its own process group, which is killed on a timeout.

Standard input comes from text or a file, standard output and error go to files and come back as text: no pipes,
so a chatty program cannot block on a full pipe.
"""

from __future__ import annotations

import contextlib
import os
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Mapping, Optional, Sequence, Union

HIDDEN_DESKTOP = "physim_hidden"

PathLike = Union[str, Path]


@dataclass
class Run:
    """How one run of an outside program went."""

    args: list
    returncode: Optional[int]  # None when it was stopped at the timeout
    timed_out: bool
    elapsed_s: float
    stdout: str
    stderr: str
    hidden_desktop: bool  # False when the hidden desktop could not be made (or off Windows)

    @property
    def ok(self) -> bool:
        return not self.timed_out and self.returncode == 0

    def describe(self) -> str:
        name = Path(str(self.args[0])).name
        if self.timed_out:
            what = f"{name} stopped after {self.elapsed_s:.0f} s (timeout)"
        else:
            what = f"{name} exited with code {self.returncode} after {self.elapsed_s:.1f} s"
        tail = (self.stderr.strip() or self.stdout.strip())[-600:]
        return what + (f"\n{tail}" if tail else "")


class ProgramFailed(RuntimeError):
    """An outside program did not finish, exited with an error, or left out the output it should write."""

    def __init__(self, message: str, run: Optional[Run] = None):
        super().__init__(message)
        self.run = run


def run(args: Sequence[PathLike], cwd: PathLike, *, timeout_s: float, stdin_text: Optional[str] = None,
        stdin_file: Optional[PathLike] = None, encoding: str = "latin-1", exclusive: Optional[str] = None,
        env: Optional[Mapping[str, str]] = None) -> Run:
    """Run ``args`` in ``cwd`` out of sight and wait for it, at most ``timeout_s`` seconds.

    ``stdin_text`` or ``stdin_file`` is the program's standard input (else it reads an empty stream). ``exclusive``
    names a program of which only one copy may run at a time on this machine (across processes): the call waits for
    the copy already running. ``env`` holds environment variables set for the program on top of this process's own
    environment (the program otherwise inherits this process's environment unchanged). The return code is the
    program's; this does not raise for a failed run."""
    if stdin_text is not None and stdin_file is not None:
        raise ValueError("give stdin_text or stdin_file, not both")
    args = [str(a) for a in args]
    environment = _environment(env)
    with contextlib.ExitStack() as stack:
        if exclusive:
            stack.enter_context(_machine_lock(exclusive))
        tmp = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix="physim_run_")))
        if stdin_text is not None:
            stdin_file = tmp / "stdin.txt"
            Path(stdin_file).write_text(stdin_text, encoding=encoding, newline="")
        out, err = tmp / "stdout.txt", tmp / "stderr.txt"
        start = time.monotonic()
        if sys.platform == "win32":
            code, hidden = _run_windows(args, Path(cwd), stdin_file, out, err, timeout_s, environment)
        else:
            code, hidden = _run_posix(args, Path(cwd), stdin_file, out, err, timeout_s, environment), False
        elapsed = time.monotonic() - start
        return Run(args=args, returncode=code, timed_out=code is None, elapsed_s=elapsed,
                   stdout=out.read_text(encoding=encoding, errors="replace"),
                   stderr=err.read_text(encoding=encoding, errors="replace"), hidden_desktop=hidden)


def _environment(extra: Optional[Mapping[str, str]]) -> Optional[dict]:
    """This process's environment with ``extra`` set on top (None: inherit it unchanged). Names compare without
    case on Windows, as Windows does."""
    if not extra:
        return None
    merged = dict(os.environ)
    for name, value in extra.items():
        if sys.platform == "win32":
            for k in [k for k in merged if k.upper() == name.upper()]:
                del merged[k]
        merged[str(name)] = str(value)
    return merged


# --- one copy at a time ---------------------------------------------------------------------------------------------

_thread_locks: dict = {}
_thread_locks_guard = threading.Lock()


@contextlib.contextmanager
def _machine_lock(name: str) -> Iterator[None]:
    """Hold the lock ``name`` against other threads of this process and other processes on this machine."""
    with _thread_locks_guard:
        tl = _thread_locks.setdefault(name, threading.Lock())
    with tl:
        path = Path(tempfile.gettempdir()) / f"physim_{name}.lock"
        with open(path, "a+b") as f:
            if sys.platform == "win32":
                import msvcrt
                f.seek(0)
                while True:
                    try:
                        msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
                        break
                    except OSError:
                        time.sleep(0.2)
                try:
                    yield
                finally:
                    f.seek(0)
                    msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                try:
                    yield
                finally:
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)


# --- POSIX ----------------------------------------------------------------------------------------------------------


def _run_posix(args: list, cwd: Path, stdin_file, out: Path, err: Path, timeout_s: float,
               env: Optional[dict] = None) -> Optional[int]:
    import signal
    import subprocess
    with contextlib.ExitStack() as stack:
        fin = stack.enter_context(open(stdin_file, "rb")) if stdin_file else subprocess.DEVNULL
        fout, ferr = stack.enter_context(open(out, "wb")), stack.enter_context(open(err, "wb"))
        p = subprocess.Popen(args, cwd=cwd, stdin=fin, stdout=fout, stderr=ferr, start_new_session=True, env=env)
        try:
            return p.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            return None
        finally:
            with contextlib.suppress(ProcessLookupError):
                os.killpg(p.pid, signal.SIGKILL)
            p.wait()


# --- Windows --------------------------------------------------------------------------------------------------------

_win = None  # the Windows API bindings, made on first use
_desktop = None  # handle of the hidden desktop, kept open while this process runs (0 if it could not be made)
_create_guard = threading.Lock()  # std handles are made inheritable only around CreateProcess


def _api():
    global _win
    if _win is not None:
        return _win
    import ctypes
    from ctypes import wintypes as w

    class STARTUPINFOW(ctypes.Structure):
        _fields_ = [("cb", w.DWORD), ("lpReserved", w.LPWSTR), ("lpDesktop", w.LPWSTR), ("lpTitle", w.LPWSTR),
                    ("dwX", w.DWORD), ("dwY", w.DWORD), ("dwXSize", w.DWORD), ("dwYSize", w.DWORD),
                    ("dwXCountChars", w.DWORD), ("dwYCountChars", w.DWORD), ("dwFillAttribute", w.DWORD),
                    ("dwFlags", w.DWORD), ("wShowWindow", w.WORD), ("cbReserved2", w.WORD),
                    ("lpReserved2", ctypes.c_void_p), ("hStdInput", w.HANDLE), ("hStdOutput", w.HANDLE),
                    ("hStdError", w.HANDLE)]

    class PROCESS_INFORMATION(ctypes.Structure):
        _fields_ = [("hProcess", w.HANDLE), ("hThread", w.HANDLE), ("dwProcessId", w.DWORD),
                    ("dwThreadId", w.DWORD)]

    class BASIC_LIMIT(ctypes.Structure):
        _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64), ("PerJobUserTimeLimit", ctypes.c_int64),
                    ("LimitFlags", w.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                    ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", w.DWORD),
                    ("Affinity", ctypes.c_size_t), ("PriorityClass", w.DWORD), ("SchedulingClass", w.DWORD)]

    class EXTENDED_LIMIT(ctypes.Structure):
        _fields_ = [("BasicLimitInformation", BASIC_LIMIT), ("IoInfo", ctypes.c_uint64 * 6),
                    ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                    ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]

    user32 = ctypes.WinDLL("user32", use_last_error=True)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    user32.CreateDesktopW.restype = w.HANDLE
    user32.CreateDesktopW.argtypes = [w.LPCWSTR, w.LPCWSTR, ctypes.c_void_p, w.DWORD, w.DWORD, ctypes.c_void_p]
    kernel32.CreateProcessW.restype = w.BOOL
    kernel32.CreateProcessW.argtypes = [w.LPCWSTR, w.LPWSTR, ctypes.c_void_p, ctypes.c_void_p, w.BOOL, w.DWORD,
                                        ctypes.c_void_p, w.LPCWSTR, ctypes.POINTER(STARTUPINFOW),
                                        ctypes.POINTER(PROCESS_INFORMATION)]
    kernel32.CreateJobObjectW.restype = w.HANDLE
    kernel32.CreateJobObjectW.argtypes = [ctypes.c_void_p, w.LPCWSTR]
    kernel32.SetInformationJobObject.argtypes = [w.HANDLE, ctypes.c_int, ctypes.c_void_p, w.DWORD]
    kernel32.AssignProcessToJobObject.argtypes = [w.HANDLE, w.HANDLE]
    kernel32.TerminateJobObject.argtypes = [w.HANDLE, w.UINT]
    kernel32.ResumeThread.argtypes = [w.HANDLE]
    kernel32.ResumeThread.restype = w.DWORD
    kernel32.WaitForSingleObject.argtypes = [w.HANDLE, w.DWORD]
    kernel32.WaitForSingleObject.restype = w.DWORD
    kernel32.GetExitCodeProcess.argtypes = [w.HANDLE, ctypes.POINTER(w.DWORD)]
    kernel32.TerminateProcess.argtypes = [w.HANDLE, w.UINT]
    kernel32.CloseHandle.argtypes = [w.HANDLE]
    _win = {"ctypes": ctypes, "w": w, "user32": user32, "kernel32": kernel32, "STARTUPINFOW": STARTUPINFOW,
            "PROCESS_INFORMATION": PROCESS_INFORMATION, "EXTENDED_LIMIT": EXTENDED_LIMIT}
    return _win


def _hidden_desktop() -> bool:
    """Make (or open) the hidden desktop once; False if Windows refused (the program then runs on the user's)."""
    global _desktop
    if _desktop is None:
        api = _api()
        h = api["user32"].CreateDesktopW(HIDDEN_DESKTOP, None, None, 0, 0x01FF, None)  # DESKTOP_ALL_ACCESS
        _desktop = h or 0
    return bool(_desktop)


def _run_windows(args: list, cwd: Path, stdin_file, out: Path, err: Path, timeout_s: float,
                 env: Optional[dict] = None):
    import msvcrt
    import subprocess
    api = _api()
    ctypes, w, k32 = api["ctypes"], api["w"], api["kernel32"]
    hidden = _hidden_desktop()

    job = k32.CreateJobObjectW(None, None)
    if not job:
        raise OSError(ctypes.get_last_error(), "CreateJobObject failed")
    try:
        limit = api["EXTENDED_LIMIT"]()
        limit.BasicLimitInformation.LimitFlags = 0x2000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not k32.SetInformationJobObject(job, 9, ctypes.byref(limit), ctypes.sizeof(limit)):  # extended limits
            raise OSError(ctypes.get_last_error(), "SetInformationJobObject failed")

        with contextlib.ExitStack() as files:
            fin = files.enter_context(open(stdin_file if stdin_file else os.devnull, "rb"))
            fout, ferr = files.enter_context(open(out, "wb")), files.enter_context(open(err, "wb"))
            handles = [msvcrt.get_osfhandle(f.fileno()) for f in (fin, fout, ferr)]
            si = api["STARTUPINFOW"](cb=ctypes.sizeof(api["STARTUPINFOW"]))
            si.dwFlags = 0x100  # STARTF_USESTDHANDLES
            si.hStdInput, si.hStdOutput, si.hStdError = handles
            if hidden:
                si.lpDesktop = HIDDEN_DESKTOP  # a desktop of this process's window station
            else:
                si.dwFlags |= 0x1  # STARTF_USESHOWWINDOW, SW_HIDE: the best that can be done on the user's desktop
                si.wShowWindow = 0
            pi = api["PROCESS_INFORMATION"]()
            cmdline = ctypes.create_unicode_buffer(subprocess.list2cmdline(args))
            flags = 0x08000000 | 0x00000004  # CREATE_NO_WINDOW | CREATE_SUSPENDED (joins the job before it runs)
            block = None
            if env is not None:  # a Unicode environment block: NAME=value entries, sorted, each ending in a NUL
                flags |= 0x00000400  # CREATE_UNICODE_ENVIRONMENT
                text = "".join(f"{k}={v}\0" for k, v in sorted(env.items(), key=lambda kv: kv[0].upper())) + "\0"
                block = ctypes.create_unicode_buffer(text, len(text))
            with _create_guard:
                for h in handles:
                    os.set_handle_inheritable(h, True)
                try:
                    ok = k32.CreateProcessW(None, cmdline, None, None, True, flags, block, str(cwd),
                                            ctypes.byref(si), ctypes.byref(pi))
                    error = ctypes.get_last_error()
                finally:
                    for h in handles:
                        os.set_handle_inheritable(h, False)
            if not ok:
                raise OSError(error, f"CreateProcess failed for {args[0]}")
            try:
                if not k32.AssignProcessToJobObject(job, pi.hProcess):
                    error = ctypes.get_last_error()
                    k32.TerminateProcess(pi.hProcess, 1)
                    raise OSError(error, "AssignProcessToJobObject failed")
                k32.ResumeThread(pi.hThread)
                ms = min(int(timeout_s * 1000), 0xFFFFFFFE)
                if k32.WaitForSingleObject(pi.hProcess, ms) != 0:  # WAIT_OBJECT_0
                    k32.TerminateJobObject(job, 1)
                    k32.WaitForSingleObject(pi.hProcess, 10000)
                    return None, hidden
                code = w.DWORD()
                k32.GetExitCodeProcess(pi.hProcess, ctypes.byref(code))
                return code.value, hidden
            finally:
                k32.CloseHandle(pi.hThread)
                k32.CloseHandle(pi.hProcess)
    finally:
        k32.TerminateJobObject(job, 1)  # anything the program left running ends with the run
        k32.CloseHandle(job)
