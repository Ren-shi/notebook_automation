"""The run record: every number explains itself (backlog 57)."""

import math

import pytest

from physim.nuclear import Experiment
from physim.nuclear.analysis import Analysis, Settings
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.planner import Planner
from physim.nuclear.record import LEFT_OUT, analysis_explanations, explanations, provenance, record_html


@pytest.fixture(scope="module")
def ni():
    p = Planner.example("coulex_ni58")
    result = Analysis(p.experiment, simulate_gammas(p.experiment, 200_000, seed=1)).run(Settings(particle_gate="all"))
    return p, result


def test_every_explanation_recomputes_to_its_value(ni):
    p, result = ni
    xs = explanations(p) + analysis_explanations(p, result)
    assert [x.key for x in explanations(p)] == ["solid_angle", "rate", "counts", "beam_time",
                                                "excitation_probability", "gamma_efficiency", "coincidences",
                                                "doppler"]
    assert [x.key for x in analysis_explanations(p, result)] == ["area", "yield", "mean_probability", "b_e2",
                                                                 "weisskopf", "beta2", "q0", "lifetime",
                                                                 "uncertainty"]
    for x in xs:
        assert x.recompute() == pytest.approx(x.value, rel=1e-9), x.key
        for part in (x.formula, x.substituted, x.meaning, x.assumptions, x.reference):
            assert part.strip(), (x.key, part)
        # The numbers shown are the inputs: each appears in the substituted formula (to its printed precision).
        for name, v in x.inputs.items():
            if isinstance(v, (int, float)) and math.isfinite(v) and v != 0:
                assert any(f"{v:.{d}g}" in x.substituted for d in (3, 4, 6)) or f"{v:.1f}" in x.substituted \
                    or f"{v:.2f}" in x.substituted or f"{v:.3f}" in x.substituted or f"{v:.4f}" in x.substituted \
                    or str(v) in x.substituted, (x.key, name, v, x.substituted)


def test_the_numbers_are_those_of_the_calculation(ni):
    p, result = ni
    r = p.rates()
    by = {x.key: x for x in explanations(p)}
    first = r["rows"][0]
    assert by["solid_angle"].value == pytest.approx(first["solid_angle_msr"])
    assert by["counts"].value == pytest.approx(first["counts_in_run"])
    assert by["beam_time"].value == pytest.approx(first["beam_time_s"] / 3600)
    assert by["coincidences"].value == pytest.approx(first["coincidence_per_s"])
    assert by["gamma_efficiency"].value == pytest.approx(
        p.experiment.gamma_detectors[0].peak_efficiency(p.experiment.excitation.energy_mev, p.experiment))
    g = p.gamma()
    assert by["doppler"].value == pytest.approx(g["doppler"][0]["mean_kev"])
    ax = {x.key: x for x in analysis_explanations(p, result)}
    assert ax["b_e2"].value == result.b_e2fm4 and ax["weisskopf"].value == result.b_wu
    assert ax["beta2"].value == result.shape["beta2"] and ax["uncertainty"].value == result.total_unc
    assert ax["area"].value == pytest.approx(result.area)
    # The lifetime follows from B(E2↓) = B(E2↑)/5 and the transition energy.
    from physim.nuclear.levels import rate_per_b

    tau_ps = 1e12 / (rate_per_b(2, 1454.0) * result.b_e2fm4 / 5)
    assert ax["lifetime"].value == pytest.approx(tau_ps, rel=1e-3)


def test_explanations_pick_the_detector_asked_for(ni):
    p, _ = ni
    cd = {x.key: x for x in explanations(p, "CD", "Ge45")}
    dssd = {x.key: x for x in explanations(p, "DSSD-L")}
    assert "CD" in cd["rate"].title and "DSSD-L" in dssd["rate"].title
    assert cd["rate"].value != dssd["rate"].value
    assert "Ge45" in cd["gamma_efficiency"].title
    plain = explanations(Planner.example("alpha_on_gold"))
    assert [x.key for x in plain] == ["solid_angle", "rate", "counts", "beam_time"]


def test_the_record_holds_the_same_content_as_the_app_view(ni):
    p, result = ni
    page = record_html(p, result)
    assert page.startswith("<!doctype html>") and "run record" in page
    for heading in ("The setup", "The nuclear data", "The method", "Every number, explained",
                    "From the peak area to B(E2)", "The uncertainty budget", "What the simulation leaves out"):
        assert heading in page, heading
    for x in explanations(p) + analysis_explanations(p, result):
        assert f'id="{x.key}"' in page
    for line in LEFT_OUT:
        assert line[:40] in page
    # The app shows the body of the same page; the planner method gives it with the last analysis.
    body = page.split("<body>", 1)[1].rsplit("</body>", 1)[0]
    assert "Every number, explained" in body
    assert p.record_html(result) == page
    without = record_html(p)
    assert "From the peak area" not in without and "What the simulation leaves out" in without
    assert len(p.explanations(result=result)) == len(explanations(p)) + len(analysis_explanations(p, result))


def test_provenance_names_every_source():
    d = Experiment.example("coulex_ni58").to_dict()
    exp = Experiment.from_dict(d)
    rows = provenance(exp)
    assert any("B(E2↑)" in q and s == "user" for q, _, s in rows)
    assert any(q.startswith("Efficiency of") and "typical" in s for q, _, s in rows)
    d["gamma_detectors"][0]["efficiency"] = "1 %"
    rows = provenance(Experiment.from_dict(d))
    assert any(q == "Efficiency of Ge90" and s == "user" for q, _, s in rows)


def test_the_report_zip_includes_the_record(tmp_path):
    from physim.nuclear.report import build

    rep = build(Experiment.example("coulex_ni58"), seed=1, events=20_000, validation=False)
    paths = rep.write(tmp_path, figures=False, root=False)
    names = {p.name for p in paths}
    assert "record.html" in names and "report.html" in names
    text = (tmp_path / "record.html").read_text(encoding="utf-8")
    assert "Every number, explained" in text and 'id="coincidences"' in text
