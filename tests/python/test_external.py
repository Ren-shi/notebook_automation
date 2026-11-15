"""physim.external: outside programs run out of sight (no windows, no console), stopped at a timeout with
everything they started. Python itself stands in for the outside programs."""

import sys
import threading
import time
from pathlib import Path

import pytest

from physim.external import HIDDEN_DESKTOP, ProgramFailed, gosia, run

PY = sys.executable
windows = pytest.mark.skipif(sys.platform != "win32", reason="Windows desktops and job objects")


def _alive(pid: int) -> bool:
    if sys.platform != "win32":
        import os
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        return True
    import ctypes
    k32 = ctypes.WinDLL("kernel32")
    k32.OpenProcess.restype = ctypes.c_void_p
    h = k32.OpenProcess(0x00100000, False, pid)  # SYNCHRONIZE
    if not h:
        return False
    try:
        return k32.WaitForSingleObject(ctypes.c_void_p(h), 0) != 0  # not signalled: still running
    finally:
        k32.CloseHandle(ctypes.c_void_p(h))


def test_input_output_and_exit_code(tmp_path):
    r = run([PY, "-c", "import sys; print(sys.stdin.read().upper()); sys.stderr.write('warn'); sys.exit(3)"],
            tmp_path, timeout_s=60, stdin_text="deck")
    assert (r.returncode, r.stdout.strip(), r.stderr, r.timed_out, r.ok) == (3, "DECK", "warn", False, False)
    assert "exited with code 3" in r.describe()


def test_runs_in_the_given_folder_with_a_file_as_input(tmp_path):
    (tmp_path / "in.txt").write_text("from a file")
    r = run([PY, "-c", "import os, sys; print(os.getcwd()); print(sys.stdin.read())"], tmp_path, timeout_s=60,
            stdin_file=tmp_path / "in.txt")
    cwd, text = r.stdout.splitlines()
    assert Path(cwd).resolve() == tmp_path.resolve() and text == "from a file" and r.ok


@windows
def test_program_runs_on_the_hidden_desktop(tmp_path):
    code = ("import ctypes; from ctypes import wintypes as w; u = ctypes.windll.user32; "
            "h = u.GetThreadDesktop(ctypes.windll.kernel32.GetCurrentThreadId()); "
            "b = ctypes.create_unicode_buffer(256); n = w.DWORD(); "
            "u.GetUserObjectInformationW(h, 2, b, 512, ctypes.byref(n)); print(b.value)")  # UOI_NAME
    r = run([PY, "-c", code], tmp_path, timeout_s=60)
    assert r.hidden_desktop and r.stdout.strip() == HIDDEN_DESKTOP


@windows
def test_a_window_never_reaches_the_users_desktop(tmp_path):
    tk = pytest.importorskip("tkinter")
    del tk
    title = f"physim hidden window test {time.time_ns()}"
    shown = tmp_path / "shown"
    code = ("import tkinter, time; r = tkinter.Tk(); r.title(%r); r.geometry('300x200'); r.update(); "
            "open(%r, 'w').close(); time.sleep(2); r.destroy()") % (title, str(shown))
    result = {}
    t = threading.Thread(target=lambda: result.setdefault("r", run([PY, "-c", code], tmp_path, timeout_s=60)))
    t.start()
    deadline = time.monotonic() + 30
    while not shown.exists() and t.is_alive() and time.monotonic() < deadline:
        time.sleep(0.05)
    import ctypes
    from ctypes import wintypes as w
    u32 = ctypes.WinDLL("user32")
    found = u32.FindWindowW(None, title)  # searches the desktop this test runs on
    on_hidden = []  # and the window is there, on the hidden desktop
    u32.OpenDesktopW.restype = ctypes.c_void_p
    desk = u32.OpenDesktopW(HIDDEN_DESKTOP, 0, False, 0x0040)  # DESKTOP_ENUMERATE
    proc = ctypes.WINFUNCTYPE(w.BOOL, w.HWND, w.LPARAM)

    def each(hwnd, _):
        buf = ctypes.create_unicode_buffer(256)
        u32.GetWindowTextW(hwnd, buf, 256)
        on_hidden.append(buf.value)
        return True

    u32.EnumDesktopWindows(ctypes.c_void_p(desk), proc(each), 0)
    u32.CloseDesktop(ctypes.c_void_p(desk))
    t.join()
    assert shown.exists(), result["r"].describe()
    assert result["r"].ok and not found and title in on_hidden


def test_timeout_stops_the_program_and_what_it_started(tmp_path):
    child = tmp_path / "child.txt"
    code = ("import subprocess, sys, time; "
            "p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)']); "
            f"open({str(child)!r}, 'w').write(str(p.pid)); time.sleep(120)")
    start = time.monotonic()
    r = run([PY, "-c", code], tmp_path, timeout_s=3)
    assert r.timed_out and r.returncode is None and not r.ok and time.monotonic() - start < 30
    time.sleep(0.5)
    assert not _alive(int(child.read_text()))


@windows
def test_what_a_program_leaves_running_ends_with_it(tmp_path):
    child = tmp_path / "child.txt"
    code = ("import subprocess, sys; "
            "p = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(120)']); "
            f"open({str(child)!r}, 'w').write(str(p.pid))")
    assert run([PY, "-c", code], tmp_path, timeout_s=60).ok
    time.sleep(0.5)
    assert not _alive(int(child.read_text()))


def test_exclusive_runs_wait_for_each_other(tmp_path):
    code = "import time; print(time.time()); time.sleep(1.0); print(time.time())"
    spans = []

    def one():
        r = run([PY, "-c", code], tmp_path, timeout_s=60, exclusive="physim_test")
        spans.append(tuple(map(float, r.stdout.split())))

    threads = [threading.Thread(target=one) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    (a0, a1), (b0, b1) = sorted(spans)
    assert b0 >= a1


def test_missing_program_raises(tmp_path):
    with pytest.raises(OSError):
        run([tmp_path / "no_such_program.exe"], tmp_path, timeout_s=10)


def test_gosia_wrapper_collects_fort22_and_reports_failures(tmp_path):
    # python with no arguments runs its standard input as a script: a stand-in for GOSIA reading its input there
    out = gosia.run(PY, "open('fort.22', 'w').write('results')\nopen('fort.99', 'w').write('scratch')\n", tmp_path)
    assert out.output == "results" and set(out.files) == {"fort.22", "fort.99"}
    with pytest.raises(ProgramFailed, match="GOSIA failed"):
        gosia.run(PY, "import sys; sys.exit(2)\n")
