"""The validation suite (backlog item 40): every check passes, the physics register says no more than the checks
show, and deliberately broken formulas are caught."""

import dataclasses
import re
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import Experiment, data, export, validation
from physim.nuclear import coulex, gamma, kinematics, rates, rutherford, stopping
from physim.nuclear.detectors import Geometry

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "tests" / "reference" / "nuclear"
REGISTER = ROOT / "docs" / "physics-register" / "README.md"
needs_reference = pytest.mark.skipif(not REFERENCE.is_dir(), reason="reference files are not shipped with the wheel")


def statuses(**kw):
    return {(r.capability, r.kind, r.name): r for r in validation.run(reference=REFERENCE, **kw)}


# -- the checks themselves --------------------------------------------------------------------------------------


def test_every_capability_has_a_tool_and_a_literature_check():
    for cap in validation.CAPABILITIES:
        kinds = {c.kind for c in validation.CHECKS if c.capability == cap}
        if cap in ("Coulomb trajectories", "Detector solid angles and response", "γ-ray Doppler shift and broadening"):
            # Exact closed forms are a stronger test than any tool; no tool computes these for arbitrary setups.
            assert "literature" in kinds, cap
        else:
            assert kinds == {"tool", "literature"}, cap


@needs_reference
def test_every_check_passes_or_waits_for_its_reference():
    results = validation.run(reference=REFERENCE)
    failed = [(r.name, r.detail) for r in results if r.status == "fail"]
    assert not failed, failed
    pending = {r.name for r in results if r.status == "pending"}
    # Exactly the references still requested from the user (tests/reference/nuclear/pending/).
    assert pending == {"closest_approach_above_barrier", "rates_vs_lise", "spectra_vs_trim", "coulex_vs_gosia"}
    assert "| Capability |" in validation.report(results)


def test_literature_checks_run_without_reference_files(tmp_path):
    for r in validation.run(kind="literature", reference=tmp_path):
        if r.name != "closest_approach_above_barrier":
            assert r.status == "pass", (r.name, r.detail)


def test_reference_files_record_their_provenance():
    if not REFERENCE.is_dir():
        pytest.skip("no reference folder")
    for path in sorted(REFERENCE.glob("*.csv")):
        header, rows = validation.read_reference(path)
        assert {"tool", "produced by", "settings"} <= set(header), path.name
        assert re.search(r"20\d\d-\d\d-\d\d", header["produced by"]), path.name
        assert rows, path.name


# -- the register says no more than the checks show ---------------------------------------------------------------


def _register_rows():
    text = REGISTER.read_text(encoding="utf-8")
    table = text.split("## Nuclear experiment planner", 1)[1].split("\n## ", 1)[0]
    rows = {}
    for line in table.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 7 and cells[0] not in ("Capability",) and not cells[0].startswith("---"):
            rows[cells[0]] = {"tool": cells[3], "literature": cells[4]}
    return rows


@needs_reference
@pytest.mark.skipif(not REGISTER.exists(), reason="no docs in this checkout")
def test_register_ticks_are_backed_by_passing_checks():
    rows = _register_rows()
    results = validation.run(reference=REFERENCE)
    for cap in validation.CAPABILITIES:
        assert cap in rows, f"register has no row for {cap}"
        for kind in ("tool", "literature"):
            mine = [r for r in results if r.capability == cap and r.kind == kind]
            cell = rows[cap][kind]
            if cell.startswith("✅"):
                assert mine and all(r.passed for r in mine), (cap, kind, [(r.name, r.status) for r in mine])
            if cell.startswith("🟡"):
                assert mine, (cap, kind, "a 🟡 entry needs a check waiting for its reference")


# -- deliberately broken formulas are caught ----------------------------------------------------------------------


def _scale_method(monkeypatch, cls, name, factor):
    orig = getattr(cls, name)

    def broken(self, *a, **k):
        return np.asarray(orig(self, *a, **k)) * factor

    monkeypatch.setattr(cls, name, broken)


def _break_masses(monkeypatch):
    orig = data.Nuclide.nuclear_mass_mev
    monkeypatch.setattr(data.Nuclide, "nuclear_mass_mev", property(lambda self: orig.fget(self) + 0.005))
    monkeypatch.setattr(data, "q_value", lambda *a, **k: 17.589 + 0.005)


def _break_kinematics(monkeypatch):
    orig = kinematics.TwoBody.at_cm

    def broken(self, theta_cm, particle="ejectile"):
        p = orig(self, theta_cm, particle)
        return dataclasses.replace(p, energy=np.asarray(p.energy) * 1.002 if np.ndim(p.energy) else p.energy * 1.002)

    monkeypatch.setattr(kinematics.TwoBody, "at_cm", broken)


def _break_ranges(monkeypatch):
    orig = stopping.Stopping.range
    monkeypatch.setattr(stopping.Stopping, "range", lambda self, e: orig(self, e) * 1.03)


def _break_rates(monkeypatch):
    orig = rates.Rates._depth_nodes
    monkeypatch.setattr(rates.Rates, "_depth_nodes", lambda self, layer: (orig(self, layer)[0],
                                                                          orig(self, layer)[1] * 1.02))


def _break_spectra(monkeypatch):
    orig = stopping.Stopping.transport_table

    def broken(self):
        t = orig(self)
        t["stopping"] = [2 * s for s in t["stopping"]]  # straggling added in the generator ×4 in variance
        return t

    monkeypatch.setattr(stopping.Stopping, "transport_table", broken)
    rates._w_table.cache_clear()


BREAKS = {
    "masses": (_break_masses, ["masses_vs_lise", "q_values"]),
    "kinematics": (_break_kinematics, ["kinematics_vs_lise", "kinematics_classical_limit"]),
    "stopping power": (lambda m: _scale_method(m, stopping.Stopping, "stopping_power", 1.10), ["stopping_vs_lise"]),
    "range": (_break_ranges, ["ranges_vs_nist"]),
    "straggling": (lambda m: _scale_method(m, stopping.Stopping, "straggling", 1.3),
                   ["straggling_vs_lise", "straggling_bohr"]),
    "rutherford": (lambda m: _scale_method(m, rutherford.Rutherford, "cross_section_cm", 1.01),
                   ["rutherford_vs_lise"]),
    "rutherford shape": (lambda m: m.setattr(rutherford.Rutherford, "cross_section_lab", lambda self, th, p="ejectile": (
        np.asarray(rutherford.Rutherford.cross_section_cm(self, th)) * np.sin(np.radians(th) / 2) ** 2, None)),
        ["rutherford_geiger_marsden"]),
    "closest approach": (lambda m: _scale_method(m, rutherford.Rutherford, "impact_parameter", 1.05),
                         ["closest_approach_vs_lise"]),
    "deflection": (lambda m: m.setattr(rutherford.Orbit, "deflection",
                                       lambda self, _o=rutherford.Orbit.deflection: _o(self) + 1e-4),
                   ["trajectories_analytic"]),
    "solid angle": (lambda m: _scale_method(m, Geometry, "solid_angle", 1 + 1e-8), ["solid_angles_closed_form"]),
    "rates": (_break_rates, ["rates_small_detector"]),
    "spectra": (_break_spectra, ["spectra_vs_analytic"]),
    "coulex": (lambda m: _scale_method(m, coulex.Coulex, "_p", 1.01), ["coulex_closed_form"]),
    "doppler": (lambda m: m.setattr(gamma, "doppler_energy", lambda e0, b, c: e0 * np.sqrt(1 - np.asarray(b) ** 2)
                                    / (1 - 0.999 * np.asarray(b) * np.asarray(c))), ["doppler_lorentz"]),
}


@pytest.mark.parametrize("what", list(BREAKS))
def test_a_broken_formula_fails_its_check(monkeypatch, what):
    if not REFERENCE.is_dir():
        pytest.skip("tool checks need the reference files")
    breaker, expected = BREAKS[what]
    before = {r.name: r.status for r in validation.run(reference=REFERENCE) if r.name in expected}
    assert all(s == "pass" for s in before.values()), before
    breaker(monkeypatch)
    after = {}
    for c in validation.CHECKS:
        if c.name in expected:
            r = validation.run(capability=c.capability, kind=c.kind, reference=REFERENCE)
            after.update({x.name: x.status for x in r if x.name == c.name})
    assert after == {name: "fail" for name in expected}, after


# -- exports ------------------------------------------------------------------------------------------------------


def test_trim_input_file():
    text = export.trim_in("4He", 5.5, [("Mylar", "10 um")], ions=500)
    lines = text.split("\r\n")
    assert lines[0].startswith("==> SRIM-2013.00")
    assert lines[2].split()[:5] == ["2", "4.0026", "5500.000", "0.00", "500"]
    assert lines[6].split()[2] == "1"  # transmitted ions saved
    assert lines[8].split()[-2:] == ["3", "1"]  # 3 elements, 1 layer
    layer = [ln for ln in lines if ln.startswith(" 1      ")][0].split()
    assert float(layer[2]) == pytest.approx(1e5)  # 10 µm in Å
    assert sum(map(float, layer[4:7])) == pytest.approx(1.0, abs=1e-5)
    # H : C : O = 8 : 10 : 4 in Mylar.
    assert float(layer[5]) / float(layer[4]) == pytest.approx(10 / 8, rel=1e-4)
    two = export.trim_in_for(Experiment.example("oxygen_on_lead_array"))
    assert two.split("\r\n")[8].split()[-2:] == ["2", "2"]
    assert "207.977" in two  # the enriched isotope's mass, not natural lead's
    with pytest.raises(ValueError, match="density"):
        export.trim_in("4He", 5.5, [("CD2", "1 um")])


def test_trim_files_in_the_request_are_current():
    trim = REFERENCE / "pending" / "trim"
    if not trim.is_dir():
        pytest.skip("no request folder")

    def lf(text):  # git may change the line endings on checkout
        return text.replace("\r\n", "\n")

    assert lf((trim / "TRIM_2.IN").read_bytes().decode()) == lf(export.trim_in("16O", 64.0, [("Au", "0.5 um")]))
    assert lf((trim / "TRIM_oxygen_on_lead_array.IN").read_bytes().decode()) == lf(export.trim_in_for(
        Experiment.example("oxygen_on_lead_array")))


def test_lise_settings():
    text = export.lise_settings(Experiment.example("alpha_on_gold"))
    assert "4He, 1.3750 MeV/u" in text and "Target: Au, 0.5 mg/cm2." in text
    a30 = [ln for ln in text.splitlines() if ln.startswith("A30")][0].split()
    assert float(a30[2]) == pytest.approx(5.4701, abs=1e-4)


def test_pictures():
    import matplotlib

    matplotlib.use("Agg")
    for r in validation.run(capability="Rutherford cross section", reference=REFERENCE):
        assert validation.plot(r) is not None
