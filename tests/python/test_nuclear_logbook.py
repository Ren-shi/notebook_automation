"""The report of an experiment (backlog item 69): one document with the setup, the plan, every run with its
spectra, the gates, the Doppler correction, the analysis and every number explained."""

import io
import zipfile

import numpy as np
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
    png, caption = logbook.run_gamma_png(p, beam, result)
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and "corrected for the recoil" in caption
    assert logbook.run_spectra_png(p, beam, result)[0][:4] == b"\x89PNG"
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
    for want in ("report.html", "setup.toml", "tables/detectors.csv", "experiment/experiment.toml",
                 "experiment/gates.json"):
        assert want in names, want
    # One report: the pages that used to stand beside it are its sections.
    assert not {"experiment.html", "record.html"} & set(names)
    if logbook.pdf_browser():
        assert z.read("report.pdf")[:5] == b"%PDF-"
    assert sum(n.endswith("summary.json") for n in names) == 2 and sum(n.endswith("events.npz") for n in names) == 2


def test_the_analysis_is_kept_with_the_run_and_outlives_the_session(tmp_path):
    """The report used to lose the analysis on any edit of the setup, on reopening the experiment, and when a
    source run became the current one: it lived only in the planner's memory."""
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="24 h", budget="30 min")
    assert p.run_analysis() is None and "No run has been analysed" in logbook.experiment_html(p)
    assert p.analyse()["available"]
    assert (p.run.folder / "analysis.json").is_file()
    kept = p.run_analysis()
    assert kept is p.analysis().history[-1] and kept.b_e2fm4 > 0
    q = Planner.open_experiment(p.experiment_folder)  # as after closing the app
    again = q.run_analysis()
    assert again.b_e2fm4 == pytest.approx(kept.b_e2fm4) and again.fit.area == pytest.approx(kept.fit.area)
    assert again.steps == kept.steps and again.spectrum["counts"] == kept.spectrum["counts"]
    assert again.total_unc == pytest.approx(kept.total_unc) and again.budget == pytest.approx(kept.budget)
    assert q.set("run", "counts_wanted", 700) and q.run_analysis() is not None  # an edit clears the memory only
    q.start_run("source", duration="10 min", source="60Co")  # the current run is now a source run
    assert q.run_analysis() is None and list(logbook.run_analyses(q)) == [1]
    page = logbook.experiment_html(q)
    assert "Run 1: from the peak area to B(E2)" in page and "No run has been analysed" not in page
    for want in ("The chain, step by step", "The uncertainty budget", "The shape of the nucleus", "Weisskopf units"):
        assert want in page, want
    q2 = Planner.open_experiment(p.experiment_folder)
    q2.run = q2.folder.load_run(1)
    assert "From the peak area to B(E2) and the shape" in q2.record_html()  # the run record finds it too


def test_the_report_draws_the_level_scheme_and_the_fitted_peak(experiment):
    pytest.importorskip("matplotlib")
    p, result = experiment
    schemes = logbook.level_schemes(p)
    assert schemes and schemes[0]["table"]["nuclide"] == "58Ni"
    assert logbook.level_scheme_png(schemes[0]["table"], 1454.0)[:4] == b"\x89PNG"
    page = logbook.experiment_html(p, result)
    assert 'alt="Level scheme of 58Ni"' in page
    if result is not None and result.b_e2fm4 == result.b_e2fm4:
        assert logbook.analysis_png(result)[:4] == b"\x89PNG" and 'alt="The fitted peak"' in page
        assert "The shape of the nucleus" in page and "Weisskopf units" in page
        fit = result.fit  # the curve drawn is the fit: at its centre, the height over the line
        at_centre = fit.height + fit.level + fit.slope * (fit.centroid - fit.origin)
        assert fit.curve(fit.centroid) == pytest.approx(at_centre)
    from physim.nuclear.report import build

    assert 'alt="Level scheme of 58Ni"' in build(p.experiment, events=5_000, validation=False).to_html()


def test_the_app_download_is_the_one_report_with_the_runs_in_it(experiment):
    pytest.importorskip("matplotlib")
    from physim.nuclear.app.figures import report_zip

    p, _ = experiment
    z = zipfile.ZipFile(io.BytesIO(report_zip(p, events=5_000)))
    names = set(z.namelist())
    assert {"report.html", "tables/detectors.csv", "figures/spectra.pdf", "setup.toml"} <= names
    assert not {"experiment.html", "record.html"} & names
    assert "Spectra of run" in z.read("report.html").decode("utf-8")


def test_the_report_is_one_document_with_every_section(experiment):
    """There were three pages (the beam-time report, the record of the experiment and the run record); they are
    the sections of one."""
    pytest.importorskip("matplotlib")
    p, result = experiment
    page = logbook.experiment_html(p, result, plan_events=5_000)
    at = [page.index(f'<h2 id="{key}">{i}. ') for i, (key, _) in enumerate(logbook.SECTIONS, 1)]
    assert at == sorted(at) and len(at) == 9
    # What each of the three held is in it: the plan's tables and figures, the runs, the numbers' explanations.
    for want in ("Detectors: coverage, cross sections, rates and beam time", 'alt="Kinematics"', "Expected peaks",
                 "Energy loss", 'alt="Spectra"', "Run 1: 152Eu source run", 'id="coincidences"',
                 "The nuclear data and where it comes from", "What the simulation leaves out",
                 "Models and their validation", "AME2020", "@media print"):
        assert want in page, want
    # Every beam run with its particles, the gates on the energy against ring, and its γ rays, both ways.
    for alt in ("Spectra of run 2", "Energy against ring of run 2", "γ-ray spectra of run 2", "Run 2 under each gate"):
        assert f'<img alt="{alt}"' in page, alt
    assert "CD rings 3-10" in page and "Particles that pass" in page
    # The Doppler correction and what it works with.
    for want in ("E₀ = E_measured × (1 − β cos α) / √(1 − β²)", "Opening half-angle", "The speeds and the shifts",
                 "Corrected line, FWHM per crystal"):
        assert want in page, want
    png, caption, rows = logbook.run_gates_png(p, p.run, result)
    assert png[:4] == b"\x89PNG" and [r["gate"] for r in rows] == ["every counted particle", "CD rings 3-10"]
    assert 0 < rows[1]["particles"] < rows[0]["particles"]  # eight of the CD's rings, of three detectors
    assert rows[1]["gammas"] <= rows[0]["gammas"] and rows[1]["in_peak"] <= rows[0]["in_peak"]
    png, caption = logbook.run_rings_png(p, p.run)
    assert png[:4] == b"\x89PNG" and "“CD rings 3-10”" in caption


def test_the_gamma_ray_spectra_reach_the_top_of_the_level_scheme(tmp_path):
    """They stopped at 1.3 × the transition energy (1525 keV for 62Ni's 1173 keV), short of the 2302 keV γ ray of
    its second 2+ state."""
    from physim.nuclear.gamma_events import spectrum_top

    p = Planner.example("coulex_ni62_pt194")
    assert spectrum_top(p.experiment) == pytest.approx(1.1 * 2.30178, rel=1e-6)
    assert spectrum_top(Planner.example("coulex_ni58").experiment) == pytest.approx(1.3 * 1.454)
    predicted = p.gamma_spectra(events=20_000)["spectra"]
    edges = next(iter(predicted.values()))["edges"]
    assert edges[-1] == pytest.approx(2.532, abs=0.001) and 1e3 * (edges[1] - edges[0]) == pytest.approx(5.08, abs=0.05)
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="24 h", budget="1 min")
    v = logbook.run_views(p, p.run)
    assert v["correction"] == "projectile" and v["total"]["corrected"]["edges"][-1] >= 2.5
    width = 1e3 * float(np.diff(v["total"]["raw"]["edges"])[0])
    assert width == pytest.approx(logbook.GAMMA_BIN_KEV, rel=0.02)
    assert dv.gamma_range(p.run.gammas(), p.run.experiment)[1] >= 2.5
    lines = {round(x["energy_kev"]): x["label"] for x in logbook.gamma_lines(p.experiment)}
    assert lines[1173] == "2+ → 0+" and lines[1163] == "4+ → 2+" and lines[2302] == "2+ → 0+" and 2048 not in lines
    # The events populate every level of this scheme, to all orders: the report has no note of levels left out.
    assert logbook._not_reached(p.experiment) == []
    page = logbook.experiment_html(p, plan_events=5_000)
    assert "Their γ rays are not in these spectra" not in page
    # First-order events excite each level in one step: the report says which levels they never reach.
    assert p.set("reaction", "order", "first")
    assert [label for label, _ in logbook._not_reached(p.experiment)] == ["0+", "4+", "0+"]
    page = logbook.experiment_html(p, plan_events=5_000)
    assert "4+ at 2336 keV" in page and "one step from the ground state" in page


def test_the_report_prints_to_pdf_where_a_browser_is_found():
    if logbook.pdf_browser() is None:
        assert logbook.html_to_pdf("<html><body><p>x</p></body></html>") is None
        pytest.skip("no Edge, Chrome or Chromium here")
    data = logbook.html_to_pdf("<!doctype html><html><body><h1>A report</h1><p>One page.</p></body></html>")
    assert data is not None and data[:5] == b"%PDF-" and len(data) > 500
