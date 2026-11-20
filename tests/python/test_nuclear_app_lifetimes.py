"""The planner app for lifetime measurements (plunger RDDS and DSAM, wave 3B of the lifetime work): the lifetime
card in Setup, the "lifetimes" kind, the Plan's lifetimes panel, runs at a plunger distance and fusion runs in Run
and Data, and the lifetime analysis.

The planner's lifetime methods (``Planner.lifetimes``, ``lifetime_spectra``, ``analyse_lifetimes``,
``dataviews.distance_spectra``) are written in parallel (wave 3A): the panels are tested here against fixture dicts
of their frozen shapes, the page must say "not available yet" without them, and one end-to-end test with the real
calls runs once they are there. No number here hangs on the lifetime examples' values (wave 3C may change them),
and nothing needs a local ENSDF copy."""

import asyncio
import copy

import numpy as np
import pytest

from physim.nuclear import workbench as wb
from physim.nuclear.experiment import Experiment
from physim.nuclear.planner import Planner

pytest.importorskip("plotly")

# -- fixtures of the frozen shapes (made-up numbers, not the examples') ---------------------------------------------

LIFETIMES = {
    "available": True, "reason": None, "reaction": "coulex", "nucleus": "98Ru", "method": "rdds",
    "distances_um": [3.0, 10.0, 30.0, 100.0, 300.0], "precision": 0.1,
    "levels": [
        {"level": "2+", "energy_kev": 650.0, "gamma_kev": 650.0, "tau_ps": 7.5, "tau_source": "ENSDF",
         "beta": 0.0444, "method": "rdds", "recommended_distances_um": [11.0, 22.0, 44.0, 88.0, 176.0],
         "sensitive_um": [15.0, 240.0],
         "separation": [{"detector": "F1", "theta_deg": 40.0, "separation_kev": 21.9, "fwhm_kev": 2.5,
                         "resolved": True},
                        {"detector": "M1", "theta_deg": 90.0, "separation_kev": 0.0, "fwhm_kev": 2.5,
                         "resolved": False}],
         "attenuation": None, "hours_per_distance": 4.25},
        {"level": "4+", "energy_kev": 1400.0, "gamma_kev": 750.0, "tau_ps": None, "tau_source": None,
         "beta": 0.0444, "method": "outside", "recommended_distances_um": [], "sensitive_um": None,
         "separation": [], "attenuation": None, "hours_per_distance": None},
    ],
    "time_to_stop_ps": 1.7, "problems": [],
    "warnings": ["Fixture warning: the detectors at 90 degrees see no Doppler shift."],
}

SPECTRA = {
    "available": True, "reason": None, "method": "rdds", "edges_kev": list(np.linspace(600.0, 700.0, 101)),
    "spectra": [{"distance_um": d, "theta_deg": th, "counts": list(np.full(100, 3.0 + d / 100)),
                 "shifted_kev": 672.0 if th < 90 else 628.0, "unshifted_kev": 650.0}
                for d in (3.0, 100.0) for th in (40.0, 140.0)],
}

ANALYSIS = {
    "available": True, "reason": None, "level": "2+", "method": "ddcm", "tau_ps": 7.6, "stat_ps": 0.31,
    "budget": {"statistics": 0.31, "stopping": 0.12, "distance": 0.05},
    "points": [{"distance_um": d, "unshifted": u, "shifted": 1000.0 - u, "error": 25.0, "tau_ps": t,
                "tau_err_ps": 0.6 if t else None}
               for d, u, t in ((3.0, 900.0, None), (30.0, 640.0, 7.3), (100.0, 300.0, 7.8), (300.0, 40.0, None))],
    "curve": {"x": [1.0, 10.0, 100.0, 1000.0], "y": [0.97, 0.75, 0.32, 0.0], "kind": "decay_curve"},
    "b_down": [610.0, 30.0], "b_unit": "e2fm4",
    "compare": [{"source": "ENSDF", "tau_ps": 8.6, "b_down": [585.0, 20.0]},
                {"source": "Coulex B(E2) of the setup", "tau_ps": 8.4, "b_down": 600.0}],
    "notes": ["Fixture note: deorientation corrected."],
}

LINESHAPE = dict(ANALYSIS, method="lineshape", points=[], curve={"x": [640.0, 650.0, 660.0], "y": [1.0, 5.0, 2.0],
                                                               "kind": "lineshape"})


def _distance_view(distances=(3.0, 30.0, 300.0)) -> dict:
    edges = np.linspace(600.0, 700.0, 51)
    return {"edges_kev": edges, "by_distance": {d: {"theta_deg": [40.0, 140.0], "counts": np.ones((2, 50)) * (1 + i),
                                                    "scale": 2.0, "runs": ["x"]}
                                                for i, d in enumerate(distances)}}


def _patch(monkeypatch, lifetimes=LIFETIMES, spectra=SPECTRA, analysis=ANALYSIS, view=None):
    """The planner's lifetime methods and dataviews.distance_spectra as the fixtures (what wave 3A writes)."""
    from physim.nuclear import dataviews

    monkeypatch.setattr(Planner, "lifetimes", lambda self, level=None: copy.deepcopy(lifetimes), raising=False)
    monkeypatch.setattr(Planner, "lifetime_spectra", lambda self, distance_um=None, level=None: copy.deepcopy(spectra),
                        raising=False)
    monkeypatch.setattr(Planner, "analyse_lifetimes",
                        lambda self, level=None, series=None, method=None: copy.deepcopy(analysis), raising=False)
    monkeypatch.setattr(dataviews, "distance_spectra",
                        lambda series, edges_kev=None, doppler=True: view or _distance_view(series.distances),
                        raising=False)


def _without(monkeypatch):
    """A physim without the planner's lifetime methods (as before wave 3A merges)."""
    from physim.nuclear import dataviews

    for name in ("lifetimes", "lifetime_spectra", "analyse_lifetimes"):
        monkeypatch.delattr(Planner, name, raising=False)
    monkeypatch.delattr(dataviews, "distance_spectra", raising=False)


# -- workbench ------------------------------------------------------------------------------------------------------

def test_the_lifetimes_kind_is_live_with_its_templates():
    kinds = {k.key: k for k in wb.KINDS}
    k = kinds["lifetimes"]
    assert k.live and k.templates == tuple(t[0] for t in wb.LIFETIME_TEMPLATES) and k.measure == "lifetimes"
    assert set(k.examples) <= set(Planner.examples())
    for key, label, holds, example, reaction, method in wb.LIFETIME_TEMPLATES:
        exp = wb.template(key)
        assert exp.lifetime is not None and exp.lifetime.method == method and exp.reaction == reaction, key
        assert wb.template_label(key) == (label, holds)
        assert wb.kind_of(exp)[0] == "lifetimes"
    assert wb.kind_of(Experiment.example("coulex_ru98_plunger"))[1] == "Lifetimes · plunger after Coulex"
    assert wb.kind_of(Experiment.example("fusion_sm143_dsam"))[1] == "Lifetimes · DSAM after fusion"
    assert wb.kind_of(Experiment.example("fusion_dy158"))[0] == "fusion-evaporation"


def test_a_new_lifetime_experiment_from_a_template():
    # Fusion-evaporation: the template's γ-ray array, and a plunger whose stopper is made for this recoil.
    exp = wb.new_setup("32S", "100 MeV", "24Mg", "0.6 mg/cm2", "fusion-evaporation",
                       detectors_from="lifetimes-fusion-plunger")
    assert exp.fusion is not None and exp.lifetime.method == "rdds" and exp.lifetime.stopper.material == "Au"
    assert len(exp.gamma_detectors) == len(Experiment.example("fusion_cr50_plunger").gamma_detectors)
    assert not wb.new_setup_notes(exp)
    # DSAM: a target without a backing gets the template's.
    exp = wb.new_setup("32S", "100 MeV", "124Sn", "0.8 mg/cm2", "fusion-evaporation",
                       detectors_from="lifetimes-fusion-dsam")
    assert exp.lifetime.method == "dsam" and exp.target.backing is not None
    # Coulomb excitation needs the excited state from ENSDF: without a copy (as in CI) the setup stays elastic and
    # a note says the lifetime measurement waits; with one it has the plunger.
    exp = wb.new_setup("98Ru", "300 MeV", "24Mg", "0.7 mg/cm2", "coulex-beam", detectors_from="lifetimes-coulex-plunger")
    assert (exp.reaction == "coulex" and exp.lifetime is not None) or (
        exp.lifetime is None and any("lifetime measurement waits" in n for n in wb.new_setup_notes(exp)))


def test_parse_distances():
    assert wb.parse_distances("10, 20, 50 um") == ["10 um", "20 um", "50 um"]
    assert wb.parse_distances("3 um; 0.02 mm 7") == ["3 um", "0.02 mm", "7 um"]
    assert wb.parse_distances("5µm 1e2") == ["5 um", "1e2 um"]
    assert wb.parse_distances("x, 5") == ["x", "5 um"]  # kept for the setup's checks to report
    assert wb.parse_distances("") == []


def test_the_lifetime_card():
    p = Planner.example("coulex_ru98_plunger")
    cards = {d.key: d for d in wb.decisions(p)}
    keys = list(cards)
    assert keys.index("lifetime") == len(keys) - 2 and keys[-1] == "run"
    card = cards["lifetime"]
    assert card.group == "Lifetime measurement" and card.title == "Plunger (RDDS)"
    n = len(p.experiment.lifetime.distances_um)
    assert card.state.startswith(f"{n} distances") and "stopper" in card.state and "deorientation" in card.state
    assert "lifetime" in [d.key for d in wb.draft_decisions(p)]  # the draft's cards, nothing computed
    assert "Lifetime measurement" in wb.GROUPS and wb.GROUPS[-1] == "Run conditions"
    # DSAM's card reads the backing.
    dsam = {d.key: d for d in wb.decisions(Planner.example("fusion_sm143_dsam"))}["lifetime"]
    assert dsam.title == "Line shapes (DSAM)" and "backing" in dsam.state
    # No card for a setup without a lifetime measurement.
    assert "lifetime" not in [d.key for d in wb.decisions(Planner.example("coulex_ni62_pt194"))]


def test_card_setup_and_changes_of_the_lifetime_section():
    p = Planner.example("coulex_ru98_plunger")
    before = copy.deepcopy(p.draft)
    assert wb.card_setup(before, "lifetime")["lifetime"]["method"] == "rdds"
    assert p.set_stopper("differential", "Ta", "4 mg/cm2") is False  # a differential plunger needs its second foil
    assert p.set_degrader2("Au", "12 mg/cm2", "100 um")
    assert p.set_lifetime(offset="5 um", deorientation={})
    lines = wb.changes(before, p.draft)
    assert "Plunger foil: kind stopper → differential" in lines and "Second foil: added" in lines
    assert "Deorientation: removed" in lines and any(x.startswith("Lifetime: offset") for x in lines)
    assert wb.card_setup(before, "lifetime") != wb.card_setup(p.draft, "lifetime")
    assert wb.changes(before, {k: v for k, v in before.items() if k != "lifetime"}) == ["Lifetime measurement: removed"]
    # An ancillary's card is found by its section (it used to raise).
    fus = Planner(wb.template("fusion-all"))
    assert wb.card_setup(fus.draft, "anc:neutron_wall") == fus.draft["neutron_wall"]


def test_the_lifetime_line_reads_the_plan(monkeypatch):
    p = Planner.example("coulex_ru98_plunger")
    _without(monkeypatch)
    line = wb.lifetime_line(p)
    assert "98Ru" in line and "not in this version yet" in line
    _patch(monkeypatch)
    p = Planner.example("coulex_ru98_plunger")
    line = wb.lifetime_line(p)
    assert line.startswith("2+ → 650 keV: τ 7.5 ps (ENSDF)")
    assert "peaks apart at 1 of 2 detectors" in line and "11–176 µm" in line and "4.25 h per distance" in line
    assert wb.lifetime_plan(p) is wb.lifetime_plan(p)  # worked out once per setup
    # A lifetime problem of the setup comes first.
    assert p.set_stopper("stopper", "Au", "0.01 mg/cm2") is False
    assert wb.lifetime_line(p).startswith("lifetime")
    # A failing plan is a line, not an error.
    monkeypatch.setattr(Planner, "lifetimes", lambda self, level=None: 1 / 0, raising=False)
    assert "could not be worked out" in wb.lifetime_line(Planner.example("fusion_cr50_plunger"))


def test_steps_say_what_a_lifetime_measurement_needs():
    p = Planner.example("coulex_ru98_plunger")
    analysis = next(s for s in wb.steps(p) if s["key"] == "analysis")
    assert analysis["status"] == "lifetimes: take a run at each distance" and analysis["state"] == "ready"
    hint = wb.next_step(p, "run")["hint"]
    assert "plunger distance" in hint


# -- figures --------------------------------------------------------------------------------------------------------

def test_lifetime_figures_render_in_both_themes():
    from physim.nuclear.app import figures as f

    figs = [f.figure_lifetime_spectra(SPECTRA, 3.0), f.figure_lifetime_spectra(dict(SPECTRA, method="dsam")),
            f.figure_distance_spectra(_distance_view(), 1), f.figure_distance_spectra(_distance_view(), 0, normalise=False),
            f.figure_decay_curve(ANALYSIS), f.figure_tau_distance(ANALYSIS), f.figure_lineshape(LINESHAPE),
            f.figure_lineshape(LINESHAPE, counts=[1, 2], edges_kev=[640, 650, 660]),
            f.figure_gg_matrix(np.ones((4, 4)), np.linspace(0, 400, 5))]
    for fig in figs:
        for theme in f.THEMES:
            f.themed(fig, theme).to_json()
    # Two angles at one distance; the marks at the unshifted and shifted energies.
    one = f.figure_lifetime_spectra(SPECTRA, 3.0)
    assert len(one.data) == 2 and len(one.layout.shapes) == 2
    # The decay curve: the unshifted share with its error, and the fit.
    dc = f.figure_decay_curve(ANALYSIS)
    assert np.allclose(dc.data[0].y, [0.9, 0.64, 0.3, 0.04]) and len(dc.data) == 2
    assert len(f.figure_tau_distance(ANALYSIS).data[0].x) == 2  # the distances with a τ
    assert len(f.figure_distance_spectra(_distance_view()).data) == 3
    assert f._shade(0, 3) != f._shade(2, 3)


# -- the page -------------------------------------------------------------------------------------------------------

def _texts(user) -> str:
    from nicegui import ui

    return " | ".join(str(e.text) for e in user.find(kind=ui.label).elements)


def _on_page(scenario, **kwargs) -> None:
    """``scenario(user)`` on the planner page built whole (``prebuild``: every stage at once)."""
    pytest.importorskip("nicegui")
    from nicegui.testing import user_simulation

    from physim.nuclear.app.page import build_page

    async def main():
        async with user_simulation(lambda: build_page(prebuild=True, **kwargs)) as user:
            await user.open("/")
            await scenario(user)
            await asyncio.sleep(0.2)

    asyncio.run(main())


def test_the_page_without_the_planner_lifetime_methods(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    _without(monkeypatch)

    async def scenario(user):
        await user.should_see("Lifetimes · Plunger (RDDS)")
        texts = _texts(user)
        # A note instead of an error, with what the setup alone gives.
        assert "is not in this version of physim yet" in texts and "The recoil: 98Ru" in texts
        assert "could not be shown" not in texts and "could not be computed" not in texts
        # The lifetime card stands in the setup panel, its group with it.
        assert "Plunger (RDDS)" in texts and "Lifetime measurement" in texts
        # The Run step asks for the plunger distance; the Analysis step has its lifetimes part.
        assert "Run each distance not taken" in " ".join(
            str(b.text) for b in user.find(kind=__import__("nicegui").ui.button).elements)
        assert "Lifetimes" in texts

    _on_page(scenario, example="coulex_ru98_plunger")


def test_the_plan_and_analysis_panels_on_the_frozen_shapes(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    _patch(monkeypatch)

    async def scenario(user):
        await user.should_see("Lifetimes · Plunger (RDDS)")
        texts = _texts(user)
        for part in ("Studied level", "7.5 ps", "Peaks apart", "1 of 2", "Beam time per distance", "4.25 h",
                     "Levels", "Shifted against unshifted peak", "Fixture warning", "lie in the sensitive range"):
            assert part in texts, part
        # The predicted spectra, a distance at a time (a page built whole builds the closed panel too).
        assert "the shifted peak (dashed)" in texts
        # The analysis: τ with its budget, B(E2), the decay curve and τ(x), the comparison and the notes.
        for part in ("7.6 ± 0.31 ps", "Total uncertainty", "610 ± 30 e2fm4", "Decay curve",
                     "τ at each distance (DDCM)", "Uncertainty budget", "Beside τ from elsewhere",
                     "Fixture note: deorientation corrected."):
            assert part in texts, part
        assert "could not be shown" not in texts

    _on_page(scenario, example="coulex_ru98_plunger")


def test_the_dsam_analysis_and_a_failing_plan(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    plan = dict(LIFETIMES, method="dsam", levels=[dict(LIFETIMES["levels"][0], method="dsam", attenuation=0.42,
                                                       separation=[])])
    _patch(monkeypatch, lifetimes=plan, analysis=LINESHAPE,
           spectra=dict(SPECTRA, method="dsam", spectra=[dict(s, distance_um=None) for s in SPECTRA["spectra"][:2]]))

    async def scenario(user):
        await user.should_see("Lifetimes · Line shapes (DSAM)")
        texts = _texts(user)
        assert "Attenuation F(τ)" in texts and "0.420" in texts and "line shape between the full Doppler" in texts
        assert "Decay curve" not in texts  # a line-shape fit has no decay curve

    _on_page(scenario, example="fusion_sm143_dsam")
    monkeypatch.setattr(Planner, "analyse_lifetimes", lambda self, **k: 1 / 0, raising=False)

    async def failing(user):
        await user.should_see("Planner.analyse_lifetimes failed")

    _on_page(failing, example="fusion_sm143_dsam")


#: The bottom of 50Cr's yrast band typed in from ENSDF (adopted levels and γ rays of 50Cr, NDS 157, 1 (2019)), as in
#: test_nuclear_fusion_runs.py, so the fusion runs here have a cascade (hence γγ coincidences) without a local
#: ENSDF copy (CI has none).
CR50_LEVELS = {
    "nuclide": "50Cr",
    "level": [{"energy": "0 keV", "jpi": "0+", "source": "ensdf"},
              {"energy": "783.31 keV", "jpi": "2+", "half_life": "9.08 ps", "source": "ensdf"},
              {"energy": "1881.42 keV", "jpi": "4+", "half_life": "2.2 ps", "source": "ensdf"},
              {"energy": "3164.06 keV", "jpi": "6+", "half_life": "0.8 ps", "source": "ensdf"}],
    "transition": [{"from": 1, "to": 0, "energy": "783.3 keV", "intensity": 100, "multipolarity": "E2"},
                   {"from": 2, "to": 1, "energy": "1098.1 keV", "intensity": 100, "multipolarity": "E2"},
                   {"from": 3, "to": 2, "energy": "1282.5 keV", "intensity": 100, "multipolarity": "E2"}]}


def _fusion_run_folder(tmp_path, monkeypatch):
    """A fusion experiment (the 50Cr plunger, its channel alone, with :data:`CR50_LEVELS`) with one short run at a
    plunger distance (half a second of beam simulated). ENSDF is hidden, as in CI."""
    from physim.nuclear import ensdf

    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    d = Experiment.example("fusion_cr50_plunger").to_dict()
    d.setdefault("levels", {})["50Cr"] = copy.deepcopy(CR50_LEVELS)
    d["reaction"]["channels"] = [ch for ch in d["reaction"]["channels"] if ch["residue"] == "50Cr"]
    d["run"].pop("channel", None)
    p = Planner(Experiment.from_dict(d))
    p.create_experiment("cr50")
    d = p.experiment.lifetime.distances_um[0]
    run = p.start_run("fusion", "1 s", "0.5 s", distance_um=d)
    return p, run, d


def test_fusion_runs_in_run_and_data(monkeypatch, tmp_path):
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    p, run, d = _fusion_run_folder(tmp_path, monkeypatch)
    assert run.kind == "fusion" and run.distance_um == d
    from physim.nuclear.app import tab_run

    assert list(tab_run.kinds(p.experiment)) == ["fusion", "source"]
    assert list(tab_run.kinds(Experiment.example("coulex_ni62_pt194"))) == ["beam", "source", "alignment"]
    rows = tab_run.describe_runs(p)
    assert rows[0]["distance"] == f"{d:g}"
    _patch(monkeypatch)

    async def scenario(user):
        await user.should_see("γγ coincidences")
        texts = _texts(user)
        assert "Distance series" in texts and "Not taken yet" in texts  # one distance of the setup's nine
        assert "Give a γ-ray energy" in texts and "pairs in the real part" in texts
        assert "could not be shown" not in texts
        # The Run step: the fusion kind, the distance, the run's numbers.
        n = len(p.experiment.lifetime.distances_um)
        assert "Residues" in texts and "Events kept" in texts and f"1 of {n} distances taken" in texts

    _on_page(scenario, experiment=str(p.folder.path))


def test_the_fusion_data_views_and_the_distance_series(monkeypatch, tmp_path):
    """The fusion run's views as the Data step computes them on the worker (no page)."""
    from types import SimpleNamespace

    from physim.nuclear.app import tab_lifetimes as tl

    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    p, run, d = _fusion_run_folder(tmp_path, monkeypatch)
    kept = {}

    def keep(planner, key, make, of_run=False):
        if (of_run,) + key not in kept:
            kept[(of_run,) + key] = make()
        return kept[(of_run,) + key]

    ctx = SimpleNamespace(keep=keep, data_options={"bins": 100}, state={"theme": "light"})
    views = tl.fusion_views(ctx, p)
    lo, hi = views["range"]
    assert lo == 0 and hi >= 500 and len(views["singles"]) == len(p.experiment.gamma_detectors)
    counts, edges = views["matrix"]
    assert counts.shape[0] == counts.shape[1] and np.allclose(counts, counts.T) and counts.sum() > 0
    # Gate on the 2+ → 0+ line as measured: at this short distance most 2+ decays are in the stopper, unshifted.
    ctx.data_options.update(fusion_gate=783.3, fusion_corrected=False)
    gated = tl.fusion_views(ctx, p)["gated"]
    assert gated is not None and len(gated[0]) == 100 and gated[0].sum() > 0  # the 4+ → 2+ line, in coincidence
    # The series this run belongs to; its spectra per distance through distance_spectra (guarded).
    series = tl.series_of(p, p.run)
    assert series is not None and series.distances == [d]
    _without(monkeypatch)
    assert tl.distance_spectra(series)["available"] is False
    _patch(monkeypatch)
    out = tl.distance_spectra(series)
    assert out["available"] and list(out["by_distance"]) == [d]
    assert tl.taken_distances(p) == {d: 1}


def test_call_guards_the_planner():
    from physim.nuclear.app import tab_lifetimes as tl

    class Bare:
        _cache: dict = {}

    out = tl.call(Bare(), "analyse_lifetimes")
    assert out["available"] is False and out["missing"]

    class Broken:
        def lifetime_spectra(self):
            raise RuntimeError("boom")

    assert "boom" in tl.call(Broken(), "lifetime_spectra")["reason"]
    assert tl._tau(0.25) == "250 fs" and tl._tau(7.5, 0.3) == "7.5 ± 0.3 ps" and tl._tau(2500) == "2.5 ns"


@pytest.mark.skipif(not hasattr(Planner, "lifetimes"), reason="needs the planner's lifetime methods (wave 3A)")
def test_end_to_end_with_the_planner(tmp_path, monkeypatch):
    """The real planner calls on an example, through the page: the panels read what the planner gives."""
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    p = Planner.example("coulex_ru98_plunger")
    plan = p.lifetimes()
    assert plan["available"] and plan["levels"]

    async def scenario(user):
        await user.should_see("Lifetimes · Plunger (RDDS)")
        texts = _texts(user)
        assert "Studied level" in texts and "Levels" in texts
        assert "could not be shown" not in texts and "failed" not in texts

    _on_page(scenario, example="coulex_ru98_plunger")


# -- editing the lifetime card on the page ---------------------------------------------------------------------------

def _inputs(user, label):
    from nicegui import ui

    return [e for e in user.find(kind=ui.input).elements if e.props.get("label") == label]


def _type(user, box, text: str) -> None:
    from nicegui.testing import UserInteraction

    box.value = text
    UserInteraction(user, {box}, None).trigger("blur")


def _click(user, label: str) -> None:
    from nicegui import ui
    from nicegui.testing import UserInteraction

    found = [e for e in user.find(kind=ui.button).elements if e.text == label]
    assert found, label
    UserInteraction(user, {found[0]}, None).click()


async def _wait(predicate, tries: int = 400) -> None:
    for _ in range(tries):
        if predicate():
            return
        await asyncio.sleep(0.05)
    assert predicate()


def _with_planner(monkeypatch, planner) -> None:
    monkeypatch.setattr(Planner, "example", classmethod(lambda cls, name: planner))


def test_editing_the_plunger_card(monkeypatch, tmp_path):
    from nicegui import ui

    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    _without(monkeypatch)
    p = Planner.example("coulex_ru98_plunger")
    applied = p.experiment
    _with_planner(monkeypatch, p)

    async def scenario(user):
        await user.should_see("Plunger (RDDS)")
        box = _inputs(user, "Distances (target to stopper)")[0]
        _type(user, box, "5, 50, 500 um")
        await _wait(lambda: "1 change not applied" in _texts(user))
        assert p.experiment is applied  # waits for Apply
        assert "Lifetime: distances" in _texts(user)
        # A differential plunger: its second foil's fields come, and the setup has a problem until they are given.
        kind = next(e for e in user.find(kind=ui.select).elements if e.props.get("label") == "Kind")
        kind.value = "differential"
        await _wait(lambda: bool(_inputs(user, "Second foil material")))
        assert "needs its second foil" in _texts(user)
        for label, text in (("Second foil material", "Au"), ("Second foil thickness", "12 mg/cm2"),
                            ("Gap after the first foil", "300 um")):
            _type(user, _inputs(user, label)[0], text)
            await asyncio.sleep(0.05)
        # The first foil is a degrader now: as thick as the stopper was, it stops the recoil (the card says so,
        # following the draft as the fields are typed in).
        await _wait(lambda: "it must let the recoil through" in _texts(user))
        _type(user, _inputs(user, "Foil thickness")[0], "1 mg/cm2")
        await asyncio.sleep(0.05)
        # Deorientation as g with C₂ and C₄: α stays, the rates go.
        toggle = next(e for e in user.find(kind=ui.toggle).elements if "g" in (e.options or {}))
        toggle.value = "g"
        await _wait(lambda: bool(_inputs(user, "g factor")))
        for label, text in (("g factor", "0.45"), ("C₂", "3.1 ps"), ("C₄", "1.4 ps")):
            _type(user, _inputs(user, label)[0], text)
            await asyncio.sleep(0.05)
        texts = _texts(user)
        assert "missing; give alpha2" not in texts and "must let the recoil through" not in texts  # followed the draft
        _click(user, "Apply")
        await _wait(lambda: p.experiment is not applied)
        lt = p.experiment.lifetime
        assert lt.distances_um == [5.0, 50.0, 500.0] and lt.stopper.kind == "differential"
        assert lt.degrader2.material == "Au" and lt.degrader2.gap_um == 300.0
        de = lt.deorientation
        assert de.g == 0.45 and de.gamma2 is None and de.alpha2 == pytest.approx(0.18)
        assert de.gamma_per_ps[2] == pytest.approx(0.45 / 3.1)

    _on_page(scenario, example="coulex_ru98_plunger")


def test_adding_and_removing_a_lifetime_measurement(monkeypatch, tmp_path):
    from nicegui import ui
    from nicegui.testing import UserInteraction

    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    _without(monkeypatch)
    p = Planner.example("fusion_dy158")
    _with_planner(monkeypatch, p)

    async def scenario(user):
        await user.should_see("No lifetime measurement")
        # The fusion plan's stage 2 has its buttons on now.
        _click(user, "Add a plunger…")
        await _wait(lambda: "Lifetime measurement: added" in _texts(user))
        assert p.experiment.lifetime is None  # in the draft only, until Apply
        _click(user, "Apply")
        await _wait(lambda: p.experiment.lifetime is not None)
        assert p.experiment.lifetime.method == "rdds" and p.experiment.lifetime.stopper is not None
        # Switch it to DSAM: the plunger's fields are kept aside, and DSAM needs a backing (a problem until then).
        method = next(e for e in user.find(kind=ui.select).elements if e.props.get("label") == "Method")
        method.value = "dsam"
        await _wait(lambda: bool(_inputs(user, "Backing material")))
        _click(user, "Discard")
        await asyncio.sleep(0.1)
        # Remove it, with the card's own button (beside its Method; the detectors' cards have one too).
        method = next(e for e in user.find(kind=ui.select).elements if e.props.get("label") == "Method")
        card = method.parent_slot.parent
        remove = next(c for c in card.default_slot.children if isinstance(c, ui.button) and c.text == "Remove")
        UserInteraction(user, {remove}, None).click()
        await _wait(lambda: "Lifetime measurement: removed" in _texts(user))

    _on_page(scenario, example="fusion_dy158")
