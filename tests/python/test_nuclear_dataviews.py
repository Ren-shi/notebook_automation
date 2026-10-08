"""The Data tab's model (backlog item 67): every spectrum of a run, the Doppler correction as the user's choice,
gates as named objects used by the analysis, the 2D views, run overlays and exports."""

import csv
import io

import numpy as np
import pytest

from physim.nuclear import dataviews as dv
from physim.nuclear.planner import Planner


def _cd_and_close_germaniums(thickness="0.5 mg/cm2", current="2 pnA") -> Planner:
    """The ⁵⁸Ni example with the CD only and the germanium detectors at 60 mm: many coincidences per minute."""
    p = Planner.example("coulex_ni58")
    d = p.draft
    d["detectors"] = [x for x in d["detectors"] if x["name"] == "CD"]
    for gd in d["gamma_detectors"]:
        gd["distance"] = "60 mm"
    d["gamma_detectors"] = [gd for gd in d["gamma_detectors"] if gd["name"] != "Ge135"]  # it would overlap Ge90
    d["beam"]["current"] = current
    d["target"]["thickness"] = thickness
    assert p._apply(d), p.problems
    return p


@pytest.fixture(scope="module")
def taken(tmp_path_factory):
    p = _cd_and_close_germaniums()
    p.create_experiment(root=tmp_path_factory.mktemp("exp"))
    p.start_run("beam", duration="10 min", budget="5 min")
    return p


def test_the_gamma_spectrum_holds_the_whole_run_and_a_view_only_zooms(taken):
    pytest.importorskip("plotly")
    from physim.nuclear.app.figures import figure_histogram

    run = taken.run
    g = run.gammas()
    lo, hi = dv.gamma_range(g)
    assert lo == 0.0 and hi >= 1.3 * g.energy_mev and hi >= float(g["measured"].max())
    name = g.detector_names()[0]
    whole = dv.gamma_spectrum(run, name, randoms="none")
    assert whole["edges"][-1] == pytest.approx(hi)
    wide = dv.gamma_spectrum(run, name, randoms="none", range=(0.0, 10.0))
    assert whole["counts"].sum() == pytest.approx(wide["counts"].sum())   # nothing of the run is cut off
    spec = dict(whole, scale=1e3, xlabel="γ-ray energy (keV)", view=(400.0, 800.0))
    fig = figure_histogram([spec])
    assert list(fig.layout.xaxis.range) == [400.0, 800.0]
    assert fig.data[0].y.sum() == pytest.approx(2 * whole["counts"].sum())  # the data are all still there
    assert figure_histogram([dict(spec, view=None)]).layout.xaxis.range is None


def test_every_spectrum_is_there_and_the_recoil_correction_finds_the_line(taken):
    run = taken.run
    grid = dv.grid(run)
    assert set(grid["particles"]) == set(run.events().detectors)
    assert set(grid["gammas"]) == set(run.gammas().detector_names())
    assert all(s["counts"].sum() > 0 for s in grid["particles"].values())
    e0 = 1e3 * run.experiment.excitation.energy_mev
    peaks = {}
    for corr in ("off", "recoil"):
        total = None
        for name in grid["gammas"]:
            h = dv.gamma_spectrum(run, name, correction=corr, randoms="none", bins=75, range=(1.3, 1.6))
            total = h["counts"] if total is None else total + h["counts"]
        centres = 1e3 * (h["edges"][:-1] + h["edges"][1:]) / 2
        peaks[corr] = centres[np.argmax(total)]
    assert abs(peaks["recoil"] - e0) <= 4.0, peaks  # within a 4 keV bin
    assert abs(peaks["off"] - e0) > abs(peaks["recoil"] - e0)
    assert "Right for this setup" in dv.correction_reason(run.experiment, "recoil")
    assert "wrong velocity" in dv.correction_reason(run.experiment, "projectile")


def test_randoms_and_singles_and_the_scaled_run(taken):
    run = taken.run
    name = run.gammas().detector_names()[0]
    true = dv.gamma_spectrum(run, name, randoms="none")["counts"].sum()
    shown = dv.gamma_spectrum(run, name, randoms="shown")["counts"].sum()
    assert shown >= true
    singles = dv.gamma_spectrum(run, name, mode="singles")
    assert singles["counts"].sum() > true and "from the rates" in singles["label"]
    whole = dv.particle_spectrum(run, "CD", scaled=True)["counts"].sum()
    real = dv.particle_spectrum(run, "CD")["counts"].sum()
    assert whole == pytest.approx(real * run.scale, rel=0.01)


def test_a_gate_on_the_energy_against_ring_view_selects_the_excited_group(tmp_path):
    # A thin target, so the elastic and excited groups separate in every ring.
    p = _cd_and_close_germaniums(thickness="0.05 mg/cm2", current="20 pnA")
    p.create_experiment(root=tmp_path)
    run = p.start_run("beam", duration="1 min")
    view = dv.energy_vs_ring(run, "CD")
    elastic = next(x for x in view["lines"] if x["group"] == "elastic")
    excited = next(x for x in view["lines"] if x["group"] == "excited")
    # A box around the excited line over a few middle rings, as one would drag it on the view.
    rings = (6, 9)
    el = dict(zip(elastic["rings"], elastic["energy"]))
    ex = dict(zip(excited["rings"], excited["energy"]))
    gap = min(el[r] - ex[r] for r in range(rings[0], rings[1] + 1))
    lo = min(ex[r] for r in range(rings[0], rings[1] + 1)) - gap / 3
    hi = max(ex[r] for r in range(rings[0], rings[1] + 1)) + gap / 3
    gate = dv.Gate.from_region("excited CD", "CD", rings, (lo, hi))
    ev = run.events()
    m = gate.mask(ev, run.experiment)
    truly = np.isin(ev["channel"], [i for i, c in enumerate(ev.channels) if "excited" in c])
    assert truly[m].mean() > 0.5, "the gate should keep mostly excited particles"
    in_box = ev.select("CD") & np.isin(ev["segment_i"], gate.rings)
    assert m[truly & in_box].mean() > 0.5, "and most of those in its rings"
    by_energy = dv.group_of(ev, run.experiment) == "excited"
    assert truly[by_energy].mean() > 0.5
    # The analysis uses the gate by name.
    p.put_gate(gate)
    out = p.analyse(gate="excited CD")
    if out["available"]:
        assert out["settings"]["particle_energy"] == tuple(gate.energy)
        assert out["settings"]["detectors"] == ["CD"]
    # Gates are kept with the experiment.
    q = Planner.open_experiment(p.experiment_folder)
    assert [g.name for g in q.gates()] == ["excited CD"] and q.gate("excited CD").energy == gate.energy


def test_two_runs_overlay_with_their_setups_named(tmp_path):
    p = Planner.example("alpha_on_gold")
    p.create_experiment(root=tmp_path)
    first = p.start_run("beam", duration="5 s")
    p.set("detector 3", "distance", "60 mm")
    second = p.start_run("beam", duration="5 s")
    ov = dv.overlay(first, second, "A45")
    labels = [s["label"] for s in ov["spectra"]]
    assert labels[0].startswith("run 1:") and labels[1].startswith("run 2:")
    assert any("distance" in d for d in ov["differences"])
    a, b = (s["counts"].sum() for s in ov["spectra"])
    assert b > 2 * a  # closer, so a larger solid angle


def test_exports_have_the_panels_counts(taken, tmp_path):
    pytest.importorskip("uproot")
    from physim.nuclear.rootio import read_root

    run = taken.run
    gate = dv.Gate("ring 3-5", ["CD"], [2, 3, 4])
    path = dv.write_root(run, tmp_path / "run.root", [gate])
    back = read_root(path)
    spec = dv.particle_spectrum(run, "CD")
    assert (back["events"]["detector"] == 0).sum() == spec["counts"].sum()
    rows = list(csv.reader(io.StringIO(dv.spectrum_csv(spec))))
    assert rows[0] == ["low", "high", "counts"] and sum(float(r[2]) for r in rows[1:]) == spec["counts"].sum()
    import uproot

    with uproot.open(path) as f:
        assert "segment_i == 2" in str(f["cut_ring_3_5"])
