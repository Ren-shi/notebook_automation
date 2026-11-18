"""The shell-model entries of the validation suite: KSHELL's results, read by physim, against published calculations
with the same Hamiltonians; physim's dimension count against published dimensions; ²⁸Si/USDB against ENSDF as
context. The checks pass on the fixtures and fail when a reference number is perturbed."""

import shutil
from pathlib import Path

import pytest

from physim.nuclear import ensdf, validation
from physim.nuclear.levels import LevelScheme

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "tests" / "reference" / "nuclear"
BENCHMARKS = REFERENCE / validation.SHELLMODEL_BENCHMARKS
CAPABILITY = "Shell-model calculations (KSHELL)"
needs_reference = pytest.mark.skipif(not BENCHMARKS.exists(), reason="reference files are not shipped with the wheel")


def run_one(name, reference):
    (r,) = validation.run(reference=reference, names={name})
    return r


def perturbed_reference(tmp_path, replace=None, ensdf_shift_kev=0.0):
    """A copy of the shell-model references with one benchmark value replaced: ``replace`` is (case, state, new
    value text). ``ensdf_shift_kev`` moves every excited ENSDF level."""
    ref = tmp_path / "nuclear"
    shutil.copytree(REFERENCE / "shellmodel" / "kshell", ref / "shellmodel" / "kshell")
    lines = BENCHMARKS.read_text(encoding="utf-8").splitlines(keepends=True)
    if replace is not None:
        case, state, new = replace
        hits = [i for i, line in enumerate(lines) if line.startswith(case + ",") and f",{state}," in line]
        assert len(hits) == 1, (case, state)
        cells = lines[hits[0]].split(",")
        cells[4] = new  # value
        lines[hits[0]] = ",".join(cells)
    (ref / validation.SHELLMODEL_BENCHMARKS).write_text("".join(lines), encoding="utf-8")
    snap = (REFERENCE / validation.SHELLMODEL_ENSDF_SI28).read_text(encoding="utf-8").splitlines(keepends=True)
    if ensdf_shift_kev:
        for i, line in enumerate(snap):
            cells = line.split(",")
            if cells[0] == "level" and cells[1] != "0":
                cells[4] = f"{float(cells[4]) + ensdf_shift_kev:g}"
                snap[i] = ",".join(cells)
    (ref / validation.SHELLMODEL_ENSDF_SI28).write_text("".join(snap), encoding="utf-8")
    return ref


def test_shell_model_checks_are_registered():
    checks = {c.name: c for c in validation.CHECKS if c.capability == CAPABILITY}
    assert checks["shellmodel_vs_published"].kind == "tool"
    assert checks["shellmodel_dimensions"].kind == "literature"
    for c in checks.values():  # every check names its sources
        assert "arXiv:" in c.reference, c.name


@needs_reference
def test_kshell_agrees_with_the_published_calculations():
    r = run_one("shellmodel_vs_published", REFERENCE)
    assert r.status == "pass", r.detail
    assert r.points == 11
    assert 0 < r.worst <= 1.0
    assert len(r.plot["x"]) == 11


@needs_reference
@pytest.mark.parametrize("case, state, new", [
    ("28Si USDB", "0+1", "-136.0"),      # E(0⁺₁): KSHELL −135.861
    ("28Si USDB", "2+1>0+1", "78.4"),    # B(E2): 78.20
    ("28Si USDB", "4+1>2+1", "112"),     # B(E2): 110.02
    ("28Si USDB", "2+1", "18.7"),        # Q_s(2⁺₁): 18.50
    ("28Si USDB", "M=0+", "93711"),      # dimensions are exact
    ("48Cr KB3", "0+2", "-28.566"),      # E(0⁺₂): −28.56458
    ("48Cr KB3", "M=0+", "1963460"),
])
def test_a_perturbed_benchmark_fails(tmp_path, case, state, new):
    r = run_one("shellmodel_vs_published", perturbed_reference(tmp_path, (case, state, new)))
    assert r.status == "fail", r.detail
    assert "outside" in r.detail


@needs_reference
def test_the_tolerance_follows_the_printed_digits(tmp_path):
    # 78.3 allows 0.1 (KSHELL gives 78.20); printed as 78.30 it would allow 0.01, and the check must fail.
    r = run_one("shellmodel_vs_published", perturbed_reference(tmp_path, ("28Si USDB", "2+1>0+1", "78.30")))
    assert r.status == "fail"
    assert validation._printed_step("78.3") == pytest.approx(0.1)
    assert validation._printed_step("110") == 1.0
    assert validation._printed_step("-32.953") == pytest.approx(0.001)


@needs_reference
def test_benchmark_dimensions_match_the_published_list():
    rows = validation.read_reference(BENCHMARKS)[1]
    in_csv = sorted(int(r["value"]) for r in rows if r["quantity"] == "dimension")
    listed = sorted(v for _, _, _, _, j2, v, _ in validation.PUBLISHED_SHELL_DIMENSIONS if j2 is None)
    assert set(in_csv) <= set(listed)


def test_physim_counts_the_published_dimensions(tmp_path):
    # Runs without reference files (the dimensions are in the code); the ENSDF context then says it is missing.
    r = run_one("shellmodel_dimensions", tmp_path)
    assert r.status == "pass", r.detail
    assert r.points == len(validation.PUBLISHED_SHELL_DIMENSIONS)
    assert "not available" in r.detail


def test_a_wrong_dimension_fails(monkeypatch, tmp_path):
    dims = list(validation.PUBLISHED_SHELL_DIMENSIONS)
    nucleus, space, z, n, j2, value, source = dims[0]
    dims[0] = (nucleus, space, z, n, j2, value + 1, source)
    monkeypatch.setattr(validation, "PUBLISHED_SHELL_DIMENSIONS", dims)
    r = run_one("shellmodel_dimensions", tmp_path)
    assert r.status == "fail"
    assert "different" in r.detail


@needs_reference
def test_ensdf_comparison_is_reported_not_judged(tmp_path):
    ctx = validation.shellmodel_si28_vs_ensdf(REFERENCE)
    assert ctx["excited"] == 7  # 2⁺₁ … 2⁺₃ of the fixture, each matched to a known level
    assert 300 < ctx["rms_kev"] < 340  # 320 keV, dominated by 0⁺₃ and 2⁺₃ (the prolate band USDB misses)
    pairs = {(a, b) for a, b, _, _ in ctx["be2"]}
    assert ("2+1", "0+1") in pairs and ("4+1", "2+1") in pairs
    r = run_one("shellmodel_dimensions", REFERENCE)
    assert r.status == "pass"
    assert "Context, not judged" in r.detail and "rms 320 keV" in r.detail
    # Moving every known level by 1 MeV changes the reported rms but not the outcome.
    shifted = run_one("shellmodel_dimensions", perturbed_reference(tmp_path, ensdf_shift_kev=1000.0))
    assert shifted.status == "pass"
    assert "rms 320 keV" not in shifted.detail


@needs_reference
@pytest.mark.skipif(not ensdf.available(), reason="no local ENSDF copy")
def test_ensdf_snapshot_matches_the_local_copy():
    scheme = LevelScheme.from_ensdf("28Si", max_energy_kev=10000.0)
    snap = validation._si28_ensdf_scheme(REFERENCE / validation.SHELLMODEL_ENSDF_SI28)
    assert [lv.jpi for lv in snap.levels] == [lv.jpi for lv in scheme.levels]
    for a, b in zip(snap.levels, scheme.levels):
        assert a.energy.value == pytest.approx(b.energy.value, abs=1e-6)
    for i in range(len(scheme.levels)):
        for f in range(i):
            ours, theirs = snap.b(i, f), scheme.b(i, f)
            if theirs is not None and scheme.levels[i].energy.value > scheme.levels[f].energy.value:
                assert ours == pytest.approx(theirs, rel=1e-3), (i, f)


@needs_reference
def test_report_lists_the_shell_model_rows():
    text = validation.report(validation.run(capability=CAPABILITY, reference=REFERENCE))
    assert text.count("Shell-model calculations (KSHELL)") == 2
