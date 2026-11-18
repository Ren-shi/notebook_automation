"""physim.external.kshell: KSHELL's namelist inputs (checked field by field against the wave 0 fixture inputs), the
run wrapper's failure detection and file handling (with a stand-in for the program), and, when a KSHELL build is
found, a real ²⁸Si/USDB run compared with the fixture logs."""

import re
import shutil
import sys
from pathlib import Path

import pytest

from physim.external import ProgramFailed, Run, kshell, run

FIXTURES = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "shellmodel"
SI28 = FIXTURES / "kshell" / "si28_usdb"
INPUTS = SI28 / "inputs"
LOGS = SI28 / "logs"

KSHELL_INPUTS = sorted(p.name for p in INPUTS.glob("Si28_usdb_?.input"))
TRANSIT_INPUTS = sorted(p.name for p in INPUTS.glob("Si28_usdb_?_?.input"))


def _text(path: Path) -> str:
    return path.read_bytes().decode("latin-1")


# --- namelists -------------------------------------------------------------------------------------------------------


def test_fixture_inputs_are_all_there():
    assert KSHELL_INPUTS == ["Si28_usdb_0.input", "Si28_usdb_4.input", "Si28_usdb_8.input"]
    assert len(TRANSIT_INPUTS) == 5


@pytest.mark.parametrize("name", KSHELL_INPUTS)
def test_kshell_input_reproduces_the_fixture(name):
    fixture = _text(INPUTS / name)
    fields = kshell.parse_namelist(fixture)
    written = kshell.kshell_input(**fields)
    assert kshell.parse_namelist(written) == fields
    assert written == fixture  # same order as kshell_ui.py (alphabetical), so the very same bytes


@pytest.mark.parametrize("name", TRANSIT_INPUTS + ["../logs_eff_1.31_0.46/transit_template.input"])
def test_transit_input_reproduces_the_fixture(name):
    fixture = _text(INPUTS / name)
    fields = kshell.parse_namelist(fixture)
    written = kshell.transit_input(**fields)
    assert kshell.parse_namelist(written) == fields
    assert "eff_charge = " + ", ".join(repr(x) for x in fields["eff_charge"]) + "\n" in written


def test_fields_written_from_python_values():
    """The ²⁸Si 2⁺ input from Python values, field for field."""
    text = kshell.kshell_input(beta_cm=0.0, eff_charge=(1.5, 0.5), fn_int="usdb.snt", fn_ptn="Si28_usdb_p.ptn",
                               fn_save_wave="Si28_usdb_j4p.wav", gl=[1.0, 0.0], gs=(5.585, -3.826), hw_type=2,
                               is_double_j=True, max_lanc_vec=200, maxiter=300, mode_lv_hdd=1, mtot=4, n_block=0,
                               n_eigen=3, n_restart_vec=10)
    assert kshell.parse_namelist(text) == kshell.parse_namelist(_text(INPUTS / "Si28_usdb_4.input"))
    assert "  eff_charge = 1.5, 0.5, \n" in text and "  is_double_j = .true.\n" in text
    assert '  fn_int = "usdb.snt"\n' in text and "  mtot = 4\n" in text


def test_value_formatting():
    text = kshell.transit_input(a=False, b=Path("x.ptn"), c='say "hi"', d=None, e=1e-5, f=-2, g=(1, 2, 3))
    assert text.splitlines() == ["&input", "  a = .false.", '  b = "x.ptn"', '  c = "say ""hi"""', "  e = 1e-05",
                                 "  f = -2", "  g = 1, 2, 3", "&end"]
    assert kshell.parse_namelist(text) == {"a": False, "b": "x.ptn", "c": 'say "hi"', "e": 1e-5, "f": -2,
                                           "g": (1, 2, 3)}
    for bad in ({"x y": 1}, {"v": object()}, {"v": float("nan")}, {"v": ()}):
        with pytest.raises((ValueError, TypeError)):
            kshell.kshell_input(**bad)


def test_numpy_numbers_are_written_as_numbers():
    np = pytest.importorskip("numpy")
    text = kshell.kshell_input(mtot=np.int64(4), eff_charge=np.array([1.5, 0.5]).tolist(), x=np.float64(0.25))
    assert kshell.parse_namelist(text) == {"mtot": 4, "eff_charge": (1.5, 0.5), "x": 0.25}


# --- failure detection -----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("line", [
    "STOP possible JJ collapse ",
    "STOP ERROR: maxiter/=1 in s.f. ",
    "ERROR STOP 1",
    "Fortran runtime error: End of file",
    "Error termination. Backtrace:",
    " error [read_sps]: index of  3 is wrong",
    " ERROR: local_dim > max_int4",
])
def test_failure_lines(line):
    log = "  some output\n   1  <H>:  -135.86116\n" + line + "\n#0  0x7ff6\n"
    assert kshell.failure_line(log) == line.strip()


@pytest.mark.parametrize("log", sorted(p.name for p in LOGS.glob("log_*.txt")))
def test_successful_fixture_logs_are_not_failures(log):
    assert kshell.failure_line(_text(LOGS / log)) is None


# --- running, with a stand-in for the program ------------------------------------------------------------------------


@pytest.fixture
def fake_install(tmp_path):
    folder = tmp_path / "kshell"
    (folder / "snt").mkdir(parents=True)
    for exe in ("kshell.exe", "transit.exe", "count_dim.exe"):
        (folder / exe).write_bytes(b"")
    (folder / "snt" / "usdb.snt").write_text("! interaction\n")
    (folder / "VERSION").write_text("KSHELL 2019 (test)\n")
    return kshell.install_from(folder)


def test_install_layout(fake_install, tmp_path):
    inst = fake_install
    assert inst.kshell.name == "kshell.exe" and inst.transit.name == "transit.exe"
    assert inst.count_dim.name == "count_dim.exe" and inst.snt_dir == inst.folder / "snt"
    assert inst.version == "KSHELL 2019 (test)"

    own = tmp_path / "own"  # KSHELL's own distribution: programs in bin/, no VERSION, no count_dim
    (own / "bin").mkdir(parents=True)
    for exe in ("kshell.exe", "transit.exe"):
        (own / "bin" / exe).write_bytes(b"")
    inst = kshell.install_from(own)
    assert inst.kshell == own / "bin" / "kshell.exe" and inst.count_dim is None and inst.snt_dir is None
    assert inst.version == "unknown"

    assert kshell.install_from(tmp_path / "empty") is None and kshell.install_from(None) is None


def test_find_install_follows_the_program_lookup(fake_install, monkeypatch):
    monkeypatch.setenv("PHYSIM_KSHELL", str(fake_install.folder))
    assert kshell.find_install().folder == fake_install.folder.resolve()


def _fake_run(writes=None, stdout="", stderr="", code=0, seen=None):
    """A stand-in for runner.run: records how it was called and writes ``writes`` (name -> bytes) in the folder."""
    def fake(args, cwd, *, timeout_s, exclusive=None, env=None, **kw):
        if seen is not None:
            seen.update(args=args, cwd=Path(cwd), timeout_s=timeout_s, exclusive=exclusive, env=env,
                        input=(Path(cwd) / args[1]).read_text())
        for name, data in (writes or {}).items():
            (Path(cwd) / name).write_bytes(data)
        return Run(args=list(args), returncode=code, timed_out=code is None, elapsed_s=0.1, stdout=stdout,
                   stderr=stderr, hidden_desktop=True)
    return fake


def _si28_0() -> str:
    return _text(INPUTS / "Si28_usdb_0.input")


def test_run_kshell_writes_input_and_log_and_collects_files(fake_install, tmp_path, monkeypatch):
    work = tmp_path / "work"
    seen = {}
    monkeypatch.setattr(kshell, "_run", _fake_run(
        writes={"Si28_usdb_j0p.wav": b"x" * 100, "tmp_lv_Si28_usdb_p.ptn_0_H_0_1": b"y" * 50,
                "tmp_snapshot_Si28_usdb_p.ptn_0_1": b"z"},
        stdout="   1  <H>:  -135.86116\n", seen=seen))
    out = kshell.run_kshell(fake_install, work, _si28_0(), name="Si28_usdb_j0p", threads=3, timeout_s=99,
                            inputs=[fake_install.snt_dir / "usdb.snt", INPUTS / "Si28_usdb_p.ptn"])
    assert seen["args"] == [fake_install.kshell.resolve(), "Si28_usdb_j0p.input"]
    assert seen["cwd"] == work and seen["timeout_s"] == 99 and seen["exclusive"] == "kshell"
    assert seen["env"] == {"OMP_NUM_THREADS": "3", "OMP_STACKSIZE": "64m", "GFORTRAN_UNBUFFERED_PRECONNECTED": "y"}
    assert seen["input"] == _si28_0()
    assert out.log == "   1  <H>:  -135.86116\n" and (work / "log_Si28_usdb_j0p.txt").read_text() == out.log
    assert set(out.files) == {"Si28_usdb_j0p.wav", "Si28_usdb_j0p.input", "log_Si28_usdb_j0p.txt"}
    assert not list(work.glob("tmp_*"))  # the temporary Lanczos files are gone
    assert (work / "usdb.snt").is_file() and (work / "Si28_usdb_p.ptn").is_file()


def test_threads_default_to_all_logical_cores(fake_install, tmp_path, monkeypatch):
    import os
    seen = {}
    monkeypatch.setattr(kshell, "_run", _fake_run(writes={"Si28_usdb_j0p.wav": b"x"}, seen=seen))
    kshell.run_kshell(fake_install, tmp_path, _si28_0(), name="j0",
                      inputs=[fake_install.snt_dir / "usdb.snt", INPUTS / "Si28_usdb_p.ptn"])
    assert seen["env"]["OMP_NUM_THREADS"] == str(os.cpu_count())


def test_missing_input_files_are_reported_before_running(fake_install, tmp_path, monkeypatch):
    monkeypatch.setattr(kshell, "_run", lambda *a, **k: pytest.fail("ran without its input files"))
    with pytest.raises(FileNotFoundError, match="fn_ptn = Si28_usdb_p.ptn"):
        kshell.run_kshell(fake_install, tmp_path, _si28_0(), name="j0", inputs=[fake_install.snt_dir / "usdb.snt"])
    with pytest.raises(FileNotFoundError):
        kshell.run_kshell(fake_install, tmp_path, _si28_0(), name="j0", inputs=[tmp_path / "nothing.snt"])
    with pytest.raises(FileNotFoundError, match="fn_load_wave_l"):
        kshell.run_transit(fake_install, tmp_path, _text(INPUTS / "Si28_usdb_0_4.input"), name="tr",
                           inputs=[fake_install.snt_dir / "usdb.snt", INPUTS / "Si28_usdb_p.ptn"])


@pytest.mark.parametrize("kw, expect", [
    ({"stdout": "start\n", "stderr": "STOP possible JJ collapse \n", "writes": {"Si28_usdb_j0p.wav": b"x"}},
     "stopped: STOP possible JJ collapse"),
    ({"stdout": "start\nFortran runtime error: End of file\n", "code": 2}, "exited with code 2"),
    ({"stdout": "start\n", "code": None}, "timeout"),
    ({"stdout": "all fine, but no wave functions\n"}, "did not write Si28_usdb_j0p.wav"),
])
def test_failures_raise_with_the_log_tail(fake_install, tmp_path, monkeypatch, kw, expect):
    monkeypatch.setattr(kshell, "_run", _fake_run(**kw))
    with pytest.raises(ProgramFailed) as e:
        kshell.run_kshell(fake_install, tmp_path, _si28_0(), name="Si28_usdb_j0p",
                          inputs=[fake_install.snt_dir / "usdb.snt", INPUTS / "Si28_usdb_p.ptn"])
    message = str(e.value)
    assert expect in message and ("start" in message or "all fine" in message)
    assert e.value.run is not None
    assert (tmp_path / "log_Si28_usdb_j0p.txt").is_file()  # the log is kept for a look


def test_remove_wavefunctions(tmp_path):
    (tmp_path / "a.wav").write_bytes(b"x" * 1000)
    (tmp_path / "b.wav").write_bytes(b"x" * 24)
    (tmp_path / "tmp_lv_x.ptn_0_J_0_1").write_bytes(b"x" * 10)
    (tmp_path / "log_a.txt").write_text("log")
    (tmp_path / "a.ptn").write_text("ptn")
    assert kshell.remove_wavefunctions(tmp_path) == 1034
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.ptn", "log_a.txt"]
    assert kshell.remove_wavefunctions(tmp_path) == 0
    assert kshell.remove_wavefunctions(tmp_path / "missing") == 0


def test_runner_sets_environment_variables_on_top_of_its_own(tmp_path, monkeypatch):
    monkeypatch.setenv("PHYSIM_TEST_INHERITED", "kept")
    r = run([sys.executable, "-c", "import os; print(os.environ['PHYSIM_TEST_INHERITED'], "
             "os.environ['PHYSIM_TEST_SET'], os.environ.get('PATH') is not None)"],
            tmp_path, timeout_s=60, env={"PHYSIM_TEST_SET": "new value"})
    assert r.ok and r.stdout.split() == ["kept", "new", "value", "True"]
    r = run([sys.executable, "-c", "import os; print(os.environ.get('PHYSIM_TEST_SET'))"], tmp_path, timeout_s=60)
    assert r.stdout.strip() == "None"


# --- a real run ------------------------------------------------------------------------------------------------------

INSTALL = kshell.find_install()
needs_kshell = pytest.mark.skipif(INSTALL is None, reason="no KSHELL build found (PHYSIM_KSHELL, programs.toml, "
                                                          "tools/kshell)")

_ENERGY = re.compile(r"^\s*(\d+)\s+<H>:\s+(-?\d+\.\d+)", re.MULTILINE)
_E2 = re.compile(r"^\s*(\d+)\(\s*(\d+)\)\s+(-?\d+\.\d+)\s+(\d+)\(\s*(\d+)\)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+"
                 r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", re.MULTILINE)


def _energies(log: str) -> dict:
    """State number -> energy (MeV), the last value printed for each (the converged one)."""
    return {int(n): float(e) for n, e in _ENERGY.findall(log)}


def _e2_lines(log: str) -> list:
    block = log[log.index("E2 transition"):]
    block = block[:block.index("M1 transition")] if "M1 transition" in block else block
    return [tuple(float(x) for x in m) for m in _E2.findall(block)]


@needs_kshell
def test_si28_usdb_against_the_fixture_logs(tmp_path):
    snt = INSTALL.snt_dir / "usdb.snt" if INSTALL.snt_dir else FIXTURES / "snt" / "usdb.snt"
    shutil.copyfile(INPUTS / "Si28_usdb_p.ptn", tmp_path / "Si28_usdb_p.ptn")
    for j2 in (0, 4, 8):  # 0+, 2+, 4+
        fields = kshell.parse_namelist(_text(INPUTS / f"Si28_usdb_{j2}.input"))
        out = kshell.run_kshell(INSTALL, tmp_path, kshell.kshell_input(**fields), name=f"Si28_usdb_j{j2}p",
                                inputs=[snt], timeout_s=600)
        assert f"Si28_usdb_j{j2}p.wav" in out.files
        want = _energies(_text(LOGS / f"log_Si28_usdb_j{j2}p.txt"))
        got = _energies(out.log)
        assert sorted(got) == sorted(want) == list(range(1, fields["n_eigen"] + 1))
        for n in want:
            assert got[n] == pytest.approx(want[n], abs=1e-3), (j2, n)  # 1 keV

    fields = kshell.parse_namelist(_text(INPUTS / "Si28_usdb_0_4.input"))  # 2+ -> 0+
    out = kshell.run_transit(INSTALL, tmp_path, kshell.transit_input(**fields), name="Si28_usdb_tr_j0p_j4p",
                             timeout_s=600)
    want = _e2_lines(_text(LOGS / "log_Si28_usdb_tr_j0p_j4p.txt"))
    got = _e2_lines(out.log)
    assert len(got) == len(want) == 9
    for g, w in zip(got, want):
        assert g[:2] == w[:2] and g[3:5] == w[3:5]  # the same pair of states
        assert g[2] == pytest.approx(w[2], abs=1e-3) and g[5] == pytest.approx(w[5], abs=1e-3)  # energies, 1 keV
        assert abs(g[7]) == pytest.approx(abs(w[7]), abs=2e-4)  # |Mred| (a state's overall sign is arbitrary)
        assert g[8] == pytest.approx(w[8], abs=2e-3) and g[9] == pytest.approx(w[9], abs=2e-3)  # B(E2) both ways

    assert kshell.remove_wavefunctions(tmp_path) > 0
    assert not list(tmp_path.glob("*.wav"))
