"""The record of an experiment and its runs (backlog item 69): the report lists every run with its setup and the
analysis of the beam run, and the run record explains the counters and the gates."""

import io
import zipfile

import pytest

from physim.nuclear import dataviews as dv
from physim.nuclear import logbook
from physim.nuclear.planner import Planner


@pytest.fixture(scope="module")
def experiment(tmp_path_factory):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path_factory.mktemp("exp"))
    p.start_run("source", duration="1 h", source="152Eu")
    p.set("target", "thickness", "0.6 mg/cm2")
    p.start_run("beam", duration="24 h", budget="2 min")
    p.put_gate(dv.Gate("CD rings 3-10", ["CD"], list(range(2, 10))))
    out = p.analyse(settings={"particle_gate": "all"})
    result = p.analysis().history[-1] if out["available"] else None
    return p, result


def test_the_report_lists_both_runs_with_their_setups_and_the_analysis(experiment):
    p, result = experiment
    page = logbook.experiment_html(p, result)
    assert "Run 1: 152Eu source run" in page and "Run 2: Beam run, 24h" in page
    assert "the setup as it started" in page and "target.thickness: 0.5 mg/cm2 → 0.6 mg/cm2" in page
    assert "predicted against measured" in page
    assert "CD rings 3-10" in page
    if result is not None:
        assert "B(E2↑) =" in page and "run's real events" in page
    assert "physics-register" in page


def test_every_run_has_its_spectra_in_the_report_with_the_doppler_correction(experiment):
    pytest.importorskip("matplotlib")
    p, result = experiment
    page = logbook.experiment_html(p, result)
    assert page.count("<figure><img alt=\"Spectra of run") == 2
    assert "Doppler-corrected for the recoil" in page and "the counts each crystal recorded from the source" in page
    beam = p.run
    v = logbook.run_views(p, beam, result)
    assert set(v["particles"]) == set(beam.events().detectors)
    assert set(v["gammas"]) == set(beam.gammas().detector_names()) and v["correction"] == "recoil"
    for s in v["gammas"].values():
        assert s["raw"]["counts"].sum() == pytest.approx(s["corrected"]["counts"].sum(), rel=0.02)
    png, caption = logbook.run_spectra_png(p, beam, result)
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and "corrected for the recoil" in caption
    source = p.folder.load_run(1)
    assert logbook.run_views(p, source)["source"] and logbook.run_spectra_png(p, source)[0][:4] == b"\x89PNG"
    # An analysis run with a named gate: the report's spectra take that gate and say so.
    assert v["gate"] is None and "every counted particle" in caption
    out = p.analyse(gate="CD rings 3-10")
    if out["available"]:
        gated = p.analysis().history[-1]
        assert logbook.run_views(p, beam, gated)["gate"].name == "CD rings 3-10"
        assert "gate “CD rings 3-10”" in logbook.run_spectra_png(p, beam, gated)[1]


def test_the_record_explains_the_counters_and_the_gate(experiment):
    p, result = experiment
    xs = logbook.run_explanations(p)
    keys = {x.key for x in xs}
    assert {"run_rate", "run_counts", "run_busy", "run_coincidences", "run_live", "gate_CD rings 3-10"} <= keys
    assert any(k.startswith("compare_Rate") for k in keys)
    for x in xs:  # the reader's check of item 57: every number recomputes from the numbers it shows
        assert x.recompute() == pytest.approx(x.value, rel=1e-9, abs=1e-12), x.key
        assert x.formula and x.substituted and x.meaning and x.assumptions
    page = p.record_html(result)
    assert "Measured rate of CD" in page and "Gate “CD rings 3-10”" in page


def test_the_zip_holds_the_record_and_the_experiment_folder(experiment):
    p, result = experiment
    z = zipfile.ZipFile(io.BytesIO(logbook.experiment_zip(p, result)))
    names = z.namelist()
    for want in ("experiment.html", "record.html", "setup.toml", "experiment/experiment.toml", "experiment/gates.json"):
        assert want in names, want
    assert any(n.startswith("report/") for n in names)
    assert sum(n.endswith("summary.json") for n in names) == 2 and sum(n.endswith("events.npz") for n in names) == 2
