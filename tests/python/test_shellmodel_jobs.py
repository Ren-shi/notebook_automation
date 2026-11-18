"""physim.nuclear.shellmodel_jobs: shell-model calculations as background jobs with a result cache.

Most tests use a stand-in for KSHELL: a small Python program, made into ``kshell.exe`` and ``transit.exe`` (with
pip's script launcher on Windows; a script with a shebang elsewhere), that prints the wave 0 fixture logs of
²⁸Si/USDB and writes empty wave-function files. It goes through the real runner, so the background process, the
job object and the one-KSHELL-at-a-time lock are the real ones. When a KSHELL build is found, one test runs ²⁸Si
for real and compares with the fixtures.
"""

import os
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from physim.external import kshell
from physim.nuclear import shellmodel, shellmodel_jobs as jobs, shellmodel_spaces as spaces
from physim.nuclear.shellmodel import ShellModelRequest

FIXTURES = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "shellmodel"
SI28 = FIXTURES / "kshell" / "si28_usdb"
INPUTS = SI28 / "inputs"
LOGS = SI28 / "logs"

#: The wave 0 fixture set: three 0+, three 2+, two 4+, with E2 and M1.
SI28_REQUEST = ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3), (8, 1, 2)), ("E2", "M1"))
GB = 1024**3

FAKE = textwrap.dedent('''\
    #!python
    """A stand-in for kshell.exe and transit.exe: prints the 28Si/USDB fixture logs."""
    import os, re, sys, time
    from pathlib import Path

    fixtures = Path(os.environ["FAKE_KSHELL_LOGS"])
    mode = os.environ.get("FAKE_KSHELL_MODE", "ok")
    text = Path(sys.argv[1]).read_text()
    fields = dict(re.findall(r"^\\s*(\\w+)\\s*=\\s*(.*?)\\s*$", text, re.M))
    transit = "fn_load_wave_l" in fields
    name = Path(sys.argv[1]).stem
    calls = os.environ.get("FAKE_KSHELL_CALLS")
    if calls:
        with open(calls, "a") as f:
            f.write(("transit " if transit else "kshell ") + name + "\\n")
    if mode == "sleep":
        Path(os.environ["FAKE_KSHELL_PIDFILE"]).write_text(str(os.getpid()))
        time.sleep(120)
    if mode == "slow":
        time.sleep(float(os.environ.get("FAKE_KSHELL_SLEEP", "1")))
    if mode == "fail" or (mode == "fail_transit" and transit):
        for i in range(60):
            print(f"fake iteration {i}")
        print("STOP 'fake: possible JJ collapse'")
        sys.exit(0)
    if transit:
        j = [re.search(r"_j(\\d+)p", fields[k]).group(1) for k in ("fn_load_wave_l", "fn_load_wave_r")]
        log = fixtures / f"log_Si28_usdb_tr_j{j[0]}p_j{j[1]}p.txt"
    else:
        log = fixtures / f"log_Si28_usdb_j{fields['mtot']}p.txt"
        Path(fields["fn_save_wave"].strip('"')).write_bytes(b"\\0" * 1000)
    sys.stdout.buffer.write(log.read_bytes())
''')


def _make_program(folder: Path, name: str) -> None:
    if sys.platform == "win32":
        distlib = pytest.importorskip("pip._vendor.distlib.scripts")
        source = folder.parent / f"{folder.name}_source"
        source.mkdir(exist_ok=True)
        (source / f"{name}.py").write_text(FAKE, encoding="utf-8")
        maker = distlib.ScriptMaker(str(source), str(folder), add_launchers=True)
        maker.executable = sys.executable
        maker.make(f"{name}.py")  # kshell.exe: pip's launcher with the script inside
        assert (folder / f"{name}.exe").is_file()
    else:
        exe = folder / name
        exe.write_text(FAKE.replace("#!python", f"#!{sys.executable}", 1), encoding="utf-8")
        exe.chmod(0o755)


@pytest.fixture(scope="session")
def fake_folder(tmp_path_factory):
    folder = tmp_path_factory.mktemp("fake_kshell")
    (folder / "snt").mkdir()
    shutil.copyfile(FIXTURES / "snt" / "usdb.snt", folder / "snt" / "usdb.snt")
    (folder / "VERSION").write_text("fake KSHELL for tests\nbuilt: never\n")
    for name in ("kshell", "transit"):
        _make_program(folder, name)
    return folder


@pytest.fixture
def cache(tmp_path, monkeypatch):
    root = tmp_path / "cache"
    monkeypatch.setenv("PHYSIM_SHELLMODEL", str(root))
    yield root
    for st in jobs.list_jobs():  # never leave a background job behind
        if st["state"] in jobs.ACTIVE:
            jobs.cancel(st["job"])


@pytest.fixture
def fake(fake_folder, cache, tmp_path, monkeypatch):
    """The fake install, with its call log; FAKE_KSHELL_MODE selects what it does."""
    calls = tmp_path / "calls.txt"
    monkeypatch.setenv("FAKE_KSHELL_LOGS", str(LOGS))
    monkeypatch.setenv("FAKE_KSHELL_CALLS", str(calls))
    monkeypatch.setenv("FAKE_KSHELL_MODE", "ok")
    install = kshell.install_from(fake_folder)
    install.calls = calls
    return install


def _calls(install) -> list:
    return install.calls.read_text().splitlines() if install.calls.is_file() else []


def _fixture_result():
    return shellmodel.read_folder(LOGS, SI28_REQUEST)


def _namelist(path: Path) -> dict:
    return kshell.parse_namelist(path.read_bytes().decode("latin-1"))


# ---------------------------------------------------------------------------------------------------------------
# Keys and transit pairs


def test_defaults_and_explicit_defaults_are_the_same_job():
    spec = spaces.interaction("usdb")
    explicit = ShellModelRequest("28Si", "usdb", ((0, 1, 3),), eff_charge=spec.eff_charge, gl=spec.gl, gs=spec.gs)
    implicit = ShellModelRequest("28Si", "usdb", ((0, 1, 3),))
    assert jobs.job_key(implicit, "b") == jobs.job_key(explicit, "b")
    assert jobs.job_key(implicit, "b") != jobs.job_key(implicit, "c")  # another build, another job
    assert jobs.job_key(implicit, "b") != jobs.job_key(ShellModelRequest("28Si", "usdb", ((0, 1, 3),),
                                                                         eff_charge=(1.31, 0.46)), "b")
    over = ShellModelRequest("28Si", "usdb", ((0, 1, 3),), threads=2, allow_over_limit=True)
    assert jobs.job_key(over, "b") == jobs.job_key(implicit, "b")  # how it runs is not what it computes


def test_transit_pairs_are_the_ones_kshells_helper_runs():
    sets = [(0, 1), (4, 1), (8, 1)]
    # The fixture's transit inputs, less 0 x 0 (no E2 or M1 between two 0+ states; KSHELL 2019's helper skips it)
    assert jobs._pairs(sets, ("E2", "M1")) == [((0, 1), (4, 1)), ((4, 1), (4, 1)), ((4, 1), (8, 1)),
                                                ((8, 1), (8, 1))]
    assert jobs._pairs(sets, ("M1",)) == [((4, 1), (4, 1)), ((8, 1), (8, 1))]
    assert jobs._pairs([(1, 1), (3, -1), (7, -1)], ("E1",)) == [((1, 1), (3, -1))]
    assert jobs._pairs([(1, 1), (3, -1), (7, -1)], ("E2",)) == [((1, 1), (1, 1)), ((3, -1), (3, -1)),
                                                                 ((3, -1), (7, -1)), ((7, -1), (7, -1))]


# ---------------------------------------------------------------------------------------------------------------
# A job from submit to result in the background


def test_background_job_from_submit_to_result(fake, cache):
    t = time.monotonic()
    job = jobs.submit(SI28_REQUEST, install=fake)
    assert time.monotonic() - t < 5  # returns at once; the work is in another process
    assert job == jobs.job_key(SI28_REQUEST, fake.version)
    folder = cache / job
    assert (folder / "request.toml").is_file()
    assert jobs.status(job)["state"] in (jobs.QUEUED, jobs.RUNNING, jobs.DONE)

    st = jobs.wait(job, timeout_s=120)
    assert st["state"] == jobs.DONE, st
    assert st["step_index"] == st["steps"] == 1 + 3 + 4 + 1
    assert st["started"] and st["finished"] and st["pid"] != os.getpid()
    assert [c["step"] for c in st["memory_checks"]] == ["kshell 0+ (3 states)", "kshell 2+ (3 states)",
                                                       "kshell 4+ (2 states)"]

    res = jobs.result(job)
    want = _fixture_result()
    assert [s.key for s in res.states] == [s.key for s in want.states]
    for s in res.states:
        assert s.energy_mev == pytest.approx(want.state(*s.key).energy_mev, abs=1e-6)
    assert res.transition("E2", (4, 1, 1), (0, 1, 1)).b_down == pytest.approx(99.8442, abs=1e-4)
    assert res.kshell_version == fake.version and res.interaction == "usdb"

    run = folder / "run"
    assert not list(run.glob("*.wav")), "wave functions are never kept"
    assert not (run / "usdb.snt").exists()
    assert sorted(p.name for p in run.glob("log_*.txt")) == sorted(
        ["log_Si28_usdb_j0p.txt", "log_Si28_usdb_j4p.txt", "log_Si28_usdb_j8p.txt",
         "log_Si28_usdb_tr_j0p_j4p.txt", "log_Si28_usdb_tr_j4p_j4p.txt", "log_Si28_usdb_tr_j4p_j8p.txt",
         "log_Si28_usdb_tr_j8p_j8p.txt"])
    assert (run / "Si28_usdb_p.ptn").read_text().replace("\r\n", "\n").splitlines()[1:] == \
        (INPUTS / "Si28_usdb_p.ptn").read_text().replace("\r\n", "\n").splitlines()[1:]
    # The inputs are the ones KSHELL's helper writes, apart from keeping the Lanczos vectors in memory when that fits
    for j2 in (0, 4, 8):
        got, want_fields = _namelist(run / f"Si28_usdb_j{j2}p.input"), _namelist(INPUTS / f"Si28_usdb_{j2}.input")
        assert got.pop("mode_lv_hdd") in (0, 1)
        want_fields.pop("mode_lv_hdd")
        assert got == want_fields, j2
    for a, b in ((0, 4), (4, 4), (4, 8), (8, 8)):
        assert _namelist(run / f"Si28_usdb_tr_j{a}p_j{b}p.input") == _namelist(INPUTS / f"Si28_usdb_{a}_{b}.input")
    assert [s["job"] for s in jobs.list_jobs()] == [job]


def test_a_second_submit_is_a_cache_hit(fake, cache):
    job = jobs.submit(SI28_REQUEST, install=fake, background=False)
    assert jobs.status(job)["state"] == jobs.DONE
    n = len(_calls(fake))
    assert n == 7
    t = time.monotonic()
    again = jobs.submit(ShellModelRequest("28Si", "USDB", ((8, 1, 2), (0, 1, 3), (4, 1, 3)),
                                          eff_charge=(1.5, 0.5)), install=fake)
    assert again == job and time.monotonic() - t < 1
    assert len(_calls(fake)) == n  # nothing ran
    assert not (cache / job / "worker.json").exists()  # no process was started
    assert jobs.result(job).state(0, 1, 1) is not None


def test_the_job_outlives_the_process_that_submitted_it(fake, cache, tmp_path):
    script = tmp_path / "submit_and_go.py"
    script.write_text(textwrap.dedent(f'''
        from physim.external import kshell
        from physim.nuclear import shellmodel_jobs as jobs
        from physim.nuclear.shellmodel import ShellModelRequest
        req = ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3)), ("E2",))
        print(jobs.submit(req, install=kshell.install_from({str(fake.folder)!r})))
    '''), encoding="utf-8")
    out = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, timeout=60,
                         env=dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path)))
    assert out.returncode == 0, out.stderr
    job = out.stdout.strip().splitlines()[-1]
    st = jobs.wait(job, timeout_s=120)  # the submitting process is gone; the worker carries on
    assert st["state"] == jobs.DONE, st
    res = jobs.result(job)
    assert {t.kind for t in res.transitions} == {"E2"}  # only what was asked for
    assert _calls(fake) == ["kshell Si28_usdb_j0p", "kshell Si28_usdb_j4p", "transit Si28_usdb_tr_j0p_j4p",
                            "transit Si28_usdb_tr_j4p_j4p"]


def test_jobs_run_one_at_a_time_and_wait_queued(fake, cache, monkeypatch):
    monkeypatch.setenv("FAKE_KSHELL_MODE", "slow")
    monkeypatch.setenv("FAKE_KSHELL_SLEEP", "0.7")
    first = jobs.submit(ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3))), install=fake)
    deadline = time.monotonic() + 60
    while jobs.status(first)["state"] != jobs.RUNNING and time.monotonic() < deadline:
        time.sleep(0.1)
    second = jobs.submit(ShellModelRequest("28Si", "usdb", ((4, 1, 3),)), install=fake)
    assert second != first
    seen = set()
    while jobs.status(first)["state"] == jobs.RUNNING and time.monotonic() < deadline:
        seen.add((jobs.status(second)["state"], jobs.status(second).get("step")))
        time.sleep(0.1)
    assert (jobs.QUEUED, "waiting for another shell-model job") in seen
    assert jobs.wait(first, 60)["state"] == jobs.DONE
    assert jobs.wait(second, 60)["state"] == jobs.DONE
    calls = _calls(fake)
    assert calls.index("kshell Si28_usdb_j4p") < calls.index("transit Si28_usdb_tr_j4p_j4p")
    assert calls[-2:] == ["kshell Si28_usdb_j4p", "transit Si28_usdb_tr_j4p_j4p"]  # the second job, after


# ---------------------------------------------------------------------------------------------------------------
# Failure and cancel


def test_a_failed_step_puts_the_log_tail_in_the_status(fake, cache, monkeypatch):
    monkeypatch.setenv("FAKE_KSHELL_MODE", "fail_transit")
    job = jobs.submit(SI28_REQUEST, install=fake)
    st = jobs.wait(job, timeout_s=120)
    assert st["state"] == jobs.FAILED
    assert st["step"] == "transit 0+ × 2+"
    assert "transit 0+ × 2+ failed" in st["message"] and "possible JJ collapse" in st["message"]
    assert "fake iteration 59" in st["log_tail"] and "STOP 'fake: possible JJ collapse'" in st["log_tail"]
    assert not list((cache / job / "run").glob("*.wav"))
    with pytest.raises(jobs.JobFailed, match="failed"):
        jobs.result(job)

    monkeypatch.setenv("FAKE_KSHELL_MODE", "ok")  # submitting again reruns a failed job
    assert jobs.submit(SI28_REQUEST, install=fake) == job
    assert jobs.wait(job, timeout_s=120)["state"] == jobs.DONE
    assert "log_tail" not in jobs.status(job)


def test_cancel_stops_the_worker_and_its_kshell(fake, cache, tmp_path, monkeypatch):
    pidfile = tmp_path / "kshell.pid"
    monkeypatch.setenv("FAKE_KSHELL_MODE", "sleep")
    monkeypatch.setenv("FAKE_KSHELL_PIDFILE", str(pidfile))
    job = jobs.submit(SI28_REQUEST, install=fake)
    deadline = time.monotonic() + 60
    while not pidfile.is_file() and time.monotonic() < deadline:
        time.sleep(0.1)
    kshell_pid = int(pidfile.read_text())
    st = jobs.status(job)
    assert st["state"] == jobs.RUNNING and st["step"] == "kshell 0+ (3 states)"
    worker_pid = st["pid"]
    assert jobs._alive(kshell_pid, None) and jobs._alive(worker_pid, None)

    st = jobs.cancel(job)
    assert st["state"] == jobs.CANCELLED and st["finished"]
    deadline = time.monotonic() + 15
    while (jobs._alive(kshell_pid, None) or jobs._alive(worker_pid, None)) and time.monotonic() < deadline:
        time.sleep(0.1)
    assert not jobs._alive(worker_pid, None)
    assert not jobs._alive(kshell_pid, None), "the KSHELL program goes with its worker"
    time.sleep(0.5)
    assert jobs.status(job)["state"] == jobs.CANCELLED
    with pytest.raises(jobs.JobFailed, match="cancelled"):
        jobs.result(job)
    assert jobs.cancel(job)["state"] == jobs.CANCELLED  # cancelling twice is harmless


def test_a_worker_that_vanished_is_reported_as_failed(cache):
    folder = cache / "0123456789abcdef"
    folder.mkdir(parents=True)
    gone = subprocess.Popen([sys.executable, "-c", "pass"])
    gone.wait()
    jobs._write_json(folder / "status.json", {"job": folder.name, "state": jobs.RUNNING, "step": "kshell 0+",
                                              "pid": gone.pid, "updated": "2026-01-01T00:00:00+00:00"})
    (folder / "worker.out").write_text("Traceback (most recent call last):\nMemoryError\n")
    st = jobs.status(folder.name)
    assert st["state"] == jobs.FAILED
    assert "stopped without finishing" in st["message"] and "MemoryError" in st["log_tail"]


# ---------------------------------------------------------------------------------------------------------------
# Memory


def _free(monkeypatch, values):
    """Make the free memory read ``values`` one after the other (the last one from then on)."""
    seq = list(values)

    def free():
        return seq.pop(0) if len(seq) > 1 else seq[0]

    monkeypatch.setattr(spaces, "available_memory_bytes", free)


def test_a_request_over_the_memory_limit_is_refused(fake, cache, monkeypatch):
    _free(monkeypatch, [20 * 1024**2])  # 20 MB free: the limit is 10 MB
    job = jobs.submit(SI28_REQUEST, install=fake, background=False)
    st = jobs.status(job)
    assert st["state"] == jobs.FAILED and st["refused"] == "memory"
    assert st["message"].startswith("Refused: Too big for this computer now")
    assert "allow_over_limit" in st["message"]
    assert _calls(fake) == []  # nothing ran

    over = ShellModelRequest(*[getattr(SI28_REQUEST, f) for f in ("nuclide", "interaction", "states",
                                                                  "transitions")], allow_over_limit=True)
    assert jobs.submit(over, install=fake, background=False) == job
    st = jobs.status(job)
    assert st["state"] == jobs.DONE
    assert any("over the memory limit" in n for n in jobs.result(job).notes)
    assert all(_namelist(cache / job / "run" / f"Si28_usdb_j{j}p.input")["mode_lv_hdd"] == 1 for j in (0, 4, 8))


def test_free_memory_is_checked_again_before_each_step(fake, cache, monkeypatch):
    # Plenty at the start and for the 0+ run, then too little before the 2+ run.
    _free(monkeypatch, [16 * GB, 16 * GB, 20 * 1024**2])
    job = jobs.submit(SI28_REQUEST, install=fake, background=False)
    st = jobs.status(job)
    assert st["state"] == jobs.FAILED and st["refused"] == "memory"
    assert st["message"].startswith("Stopped before kshell 2+ (3 states): the free memory fell to")
    assert [c["fits"] for c in st["memory_checks"]] == [True, False]
    assert _calls(fake) == ["kshell Si28_usdb_j0p"]
    assert not list((cache / job / "run").glob("*.wav"))


def test_lanczos_vectors_go_to_disk_when_memory_becomes_short(fake, cache, monkeypatch):
    one = spaces.estimate(ShellModelRequest("28Si", "usdb", ((4, 1, 3),)), free_bytes=64 * GB)
    between = one.memory_bytes_on_disk_mode + one.memory_bytes_in_memory_mode  # limit between the two modes
    _free(monkeypatch, [16 * GB, 16 * GB, between, 16 * GB])
    job = jobs.submit(SI28_REQUEST, install=fake, background=False)
    assert jobs.status(job)["state"] == jobs.DONE
    run = cache / job / "run"
    assert [_namelist(run / f"Si28_usdb_j{j}p.input")["mode_lv_hdd"] for j in (0, 4, 8)] == [0, 1, 0]


# ---------------------------------------------------------------------------------------------------------------
# The cache's size, and odd requests


def test_clean_cache_removes_the_oldest_finished_jobs_first(cache):
    now = time.time()
    names = [f"{i:016x}" for i in range(5)]
    for age, name in enumerate(names):  # names[0] is the newest
        folder = cache / name
        folder.mkdir(parents=True)
        jobs._write_json(folder / "status.json", {"job": name, "state": jobs.DONE})
        (folder / "results.json").write_bytes(b"x" * 1000)
        for p in folder.iterdir():
            os.utime(p, (now - 100 * age, now - 100 * age))
    busy = cache / names[4]  # the oldest, but running
    jobs._write_json(busy / "status.json", {"job": names[4], "state": jobs.RUNNING, "pid": os.getpid()})
    (busy / "results.json").unlink()
    os.utime(busy / "status.json", (now - 1000, now - 1000))
    removed = jobs.clean_cache(limit_bytes=2500, keep=(names[3],))
    assert removed == [names[2], names[1]]
    assert sorted(p.name for p in cache.iterdir()) == [names[0], names[3], names[4]]
    assert jobs.clean_cache(limit_bytes=10**9) == []


def test_unknown_jobs_and_bad_ids(cache):
    with pytest.raises(KeyError):
        jobs.status("0123456789abcdef")
    with pytest.raises(KeyError):
        jobs.result("../../etc")
    assert jobs.list_jobs() == []


def test_submit_refuses_what_cannot_run(fake, cache, monkeypatch):
    with pytest.raises(ValueError, match="does not fit"):
        jobs.submit(ShellModelRequest("48Ca", "usdb", ((0, 1, 1),)), install=fake)
    monkeypatch.setattr(kshell, "find_install", lambda: None)
    with pytest.raises(FileNotFoundError, match="PHYSIM_KSHELL"):
        jobs.submit(SI28_REQUEST)
    assert not cache.exists() or not list(cache.iterdir())


# ---------------------------------------------------------------------------------------------------------------
# The real thing


INSTALL = kshell.find_install()


@pytest.mark.skipif(INSTALL is None, reason="no KSHELL build found (PHYSIM_KSHELL, programs.toml, tools/kshell)")
def test_si28_usdb_for_real_in_the_background(cache):
    job = jobs.submit(SI28_REQUEST, install=INSTALL)
    st = jobs.wait(job, timeout_s=900)
    assert st["state"] == jobs.DONE, st
    res, want = jobs.result(job), _fixture_result()
    assert sorted(s.key for s in res.states) == sorted(s.key for s in want.states)
    for s in res.states:
        assert s.energy_mev == pytest.approx(want.state(*s.key).energy_mev, abs=1e-3)  # 1 keV
    b, b_want = (r.transition("E2", (4, 1, 1), (0, 1, 1)).b_down for r in (res, want))
    assert b == pytest.approx(b_want, rel=1e-3)  # 0.1 %
    assert not list((cache / job / "run").glob("*.wav"))

    t = time.monotonic()
    assert jobs.submit(SI28_REQUEST, install=INSTALL) == job
    assert time.monotonic() - t < 1 and jobs.status(job)["state"] == jobs.DONE  # a cache hit
