"""The report leads with the results (the results ledger, wave 1): the order of its sections, a projected row
shown with its closure and never pulled against the reference, the comparisons still to run, the numbers'
formatting and the direction of every B(E2). Synthetic rows and setups without runs: seconds."""

from types import SimpleNamespace

import pytest

from physim.nuclear import compare, ensdf, logbook, nucleus, results
from physim.nuclear.planner import Planner

KEY = "62Ni:2+1>0+1:E2"


def _rows() -> list:
    """62Ni's B(E2) of the first 2+ by four routes: the evaluation, a measurement given in Weisskopf units, a
    shell-model value without an uncertainty, and a simulated run analysed back."""
    def row(value, unit, method, evidence, origin, **more):
        return results.Row("62Ni", "B(E2)", KEY, value, unit, method, evidence, origin=origin, **more)

    return [row(884.0, "e2fm4", "evaluation", "evaluated", results.reference("ENSDF", 2012), direction="up",
                total=13.0),
            row(12.3, "W.u.", "coulex-fit", "measured", results.reference("Allmond et al.", 2014),
                direction="down", total=0.4),
            row(408.0, "e2fm4", "shell-model", "theory", results.shell_model_job("JUN45", "f5pg9"), direction="up"),
            row(918.6, "e2fm4", "coulex-yield", "projected", results.from_run("62Ni-on-194Pt", "03-beam-24h"),
                direction="up", stat=13.0, sys=63.7, budget={"efficiency": 45.9, "beam energy": 31.0}, truth=884.1)]


def test_numbers_have_matching_decimals_and_no_exponent_form():
    assert logbook.plus_minus(745, 130) == "745 ± 130"  # it was "745 ± 1.3e+02"
    assert logbook.plus_minus(918.6259, 12.49, 63.66) == "919 ± 12 ± 64"
    assert logbook.plus_minus(8.399, 0.32) == "8.40 ± 0.32"
    assert logbook.number(437100.0) == "437,100" and logbook.number(0.002265) == "0.002265"
    assert logbook.number(3.2e12) == "3.2 × 10¹²" and logbook.number(None) == "—"
    assert logbook.number(9.6e12, 1) == "1 × 10¹³" and logbook.number(9.99996e12) == "1 × 10¹³"
    # The power of ten of a value out of the ordinary is carried, never dropped.
    assert logbook.matched(5e-5, 1e-5, 2e-5) == ("5.0", "1.0", "2.0", "× 10⁻⁵")
    assert logbook.plus_minus(2.5e-12, 1e-13) == "(2.50 ± 0.10) × 10⁻¹²"
    assert logbook.asymmetric(2.5e-12, 1e-13, 2e-13) == "(2.50 +0.10 −0.20) × 10⁻¹²"
    tau = results.Row("62Ni", "tau", "62Ni:2+1", 2.5e-12, "s", "rdds", "measured", stat=1e-13, sys=2e-13,
                      origin=results.reference("A. Author", 2020))
    assert logbook.row_value(tau) == "(2.50 ± 0.10 (stat) ± 0.20 (sys)) × 10⁻¹² s"
    # An uncertainty that is not zero is never printed as 0; one that is no number is not "± 0" either.
    assert logbook.plus_minus(919, 0.3, 64) == "919.0 ± 0.3 ± 64.0"
    assert logbook.plus_minus(919, float("nan"), 64) == "919 ± n/a ± 64" == logbook.plus_minus(919, -1.0, 64)
    assert logbook.plain_numbers("2.293e+06 particles, ⟨P⟩ = 1.910e-01") == "2,293,000 particles, ⟨P⟩ = 0.191"
    # Settings an analysis wrote into a note as a Python dict are words.
    assert logbook.note_words("Changed: {'correction': 'projectile'}; done") == "Changed: correction = projectile; done"


def test_every_b_says_up_or_down_and_weisskopf_units_are_down():
    evaluated, measured, theory, projected = _rows()
    assert logbook.quantity_words("B(E2)", KEY, "up") == "B(E2; 0⁺₁ → 2⁺₁)↑"
    assert logbook.quantity_words("B(E2)", KEY, "down") == "B(E2; 2⁺₁ → 0⁺₁)↓"
    assert logbook.row_value(projected) == "919 ± 13 (stat) ± 64 (sys) e²fm⁴"
    assert logbook.row_value(evaluated, direction=True) == "884 ± 13 e²fm⁴ ↑"
    assert logbook.row_value(measured) == "12.30 ± 0.40 W.u. ↓"
    # The excitation's value is never printed in Weisskopf units: the decay's is, marked so.
    assert "B(E2; 2⁺₁ → 0⁺₁)↓ = 183.7 e²fm⁴ = 12.6 W.u. ↓" in logbook.decay_words(projected)
    # Where the spins of the two levels are not known no ratio is assumed: no ↓ value, no Weisskopf units.
    assert logbook.decay_of(918.6, KEY)[2] == 5 and logbook.decay_of(918.6, "62Ni:?@1173>0+1:E2") is None
    assert logbook.decay_of(918.6, None) is None
    assert logbook.unit_words("e2fm2") == "e²fm²" and logbook.unit_words("e2fm6") == "e²fm⁶"
    assert logbook.wu_marked("from B(E2) = 12.13 W.u. in ENSDF") == "from B(E2)↓ = 12.13 W.u. in ENSDF"


def test_a_projected_row_is_closed_against_its_truth_and_not_pulled():
    page = logbook.comparison_html(compare.compare(_rows()), here="62Ni-on-194Pt")
    independent, projected = page.split("<h3>Projected")
    assert "the reference" in independent and "+1.4 %, pull +0.39" in independent and "agree" in independent
    assert "given as 12.30 ± 0.40 W.u. ↓" in independent  # converted to B(E2)↑, and said
    assert "× 0.46 (no uncertainty: no pull)" in independent  # the shell model: a ratio, no verdict
    assert "this experiment, run 3" not in independent
    table = projected.split("<tbody>", 1)[1].split("</tbody>")[0]  # the rows of the projected table
    assert "this experiment, run 3" in table and "884.1" in table and "+0.53σ: closes" in table
    assert "pull" not in table and "agree" not in table  # closure against the truth, not a pull
    assert ("This setup reaches ±7.1 %; the value is known to ±1.5 % (ENSDF 2012, the reference); the best single "
            "measurement reaches ±3.3 % (Allmond et al. 2014).") in projected
    assert logbook.RULE in page


def test_the_comparisons_still_to_run_for_a_nucleus_with_one_route_done():
    evaluated, measured, theory, projected = _rows()
    # The states are the nucleus page's: a simulated run is not a measurement, so its route is "projected", not
    # done. A measured value makes it done. A route that does not apply is crossed out with the reason.
    routes = nucleus.routes("62Ni", "B(E2)", KEY, [projected, measured], tests=False)
    assert {r["method"]: r["status"] for r in routes if r["kind"]} == {
        "coulex-yield": "projected", "coulex-fit": "done", "rdds": "in physim", "dsam": "in physim",
        "inelastic-scattering": "later"}
    off = next(r for r in routes if r["method"] == "dsam")
    off.update(applies=False, reason="τ = 5 ps is above the range")
    page = logbook.to_run_html(routes, [projected, measured], here="62Ni-on-194Pt", unanalysed=[1, 2], outside=True)
    assert "✗ <b>DSAM (line shape).</b>" in page and "Does not apply: τ = 5 ps is above the range." in page
    assert "◐ <b>Coulex from the yield (first order).</b>" in page
    assert "Projected: a simulated run, no measurement yet: this experiment, run 3" in page
    assert "☑ <b>Coulex fitted to all orders.</b>" in page and "Done: Allmond et al. 2014" in page
    assert "☐ <b>RDDS (decay curve, DDCM).</b>" in page
    assert "physim has this kind of experiment (Lifetimes (RDDS, DSAM))" in page
    assert "☐ <b>Inelastic scattering.</b>" in page and "does not have this kind of experiment yet" in page
    assert "☐ <b>Runs 1 and 2.</b> Taken and not analysed" in page
    for missing in ("No evaluated value", "No shell-model value", "GOSIA"):
        assert missing in page, missing
    # With the evaluation, a measurement and the shell model entered, nothing independent is missing.
    assert all(x["have"] for x in logbook.evidence_missing(_rows()))


def test_the_report_and_the_nucleus_page_show_the_same_rows_and_route_states(tmp_path, monkeypatch):
    """A τ entered for the level gave a B(E2) on the nucleus page and not in the report."""
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setattr(nucleus, "decay", lambda nuclide, level: (
        results.TransitionKey.parse(KEY), {"energy_kev": 1172.9, "branching": 1.0, "conversion": 0.0, "mixing": None},
        ""))
    evaluated, measured, theory, projected = _rows()
    results.add_points([evaluated, results.Row(
        "62Ni", "tau", "62Ni:2+1", 2.01, "ps", "dsam", "measured", total=0.12,
        origin=results.reference("Kenn et al.", 2001))])
    p = Planner.example("coulex_ni62_pt194")
    p.create_experiment(root=tmp_path / "experiments")
    here = p.experiment_folder.name
    own = results.Row("62Ni", "B(E2)", KEY, 918.6, "e2fm4", "coulex-yield", "projected", direction="up", stat=13.0,
                      sys=63.7, truth=884.1, origin=results.from_run(here, "03-beam-24h"))
    results.append(p.experiment_folder, own)
    before = {f: f.read_bytes() for f in tmp_path.rglob("*") if f.is_file()}
    sec = {name: [] for name, _ in logbook.RESULTS}
    c = logbook._routes_section(sec, p, "62Ni", "B(E2)", KEY, [own], own)
    page = nucleus.observable_view("62Ni", "B(E2)", KEY, root=tmp_path / "experiments")
    assert {f: f.read_bytes() for f in tmp_path.rglob("*") if f.is_file()} == before  # nothing is written
    ids = lambda x: ({e.id for e in x.entries}, {q.id for q in x.projected}, x.reference.id)  # noqa: E731
    assert ids(c) == ids(page["comparison"]) and len(c.entries) == 2 and len(c.projected) == 1
    view = logbook.observable_rows(p, "62Ni", "B(E2)", KEY, [own])[0]
    assert view["routes"] == page["routes"] and view["checks"] == page["checks"]
    made = "Kenn et al. 2001, DSAM (line shape), converted from τ = 2.01 ± 0.12 ps"
    beside, to_run = "".join(sec["routes"]), "".join(sec["to-run"])
    assert made in beside and "it is no second result" in beside
    assert "☑ <b>DSAM (line shape).</b>" in to_run and made in to_run
    assert "◐ <b>Coulex from the yield (first order).</b>" in to_run and "no measurement yet" in to_run
    assert "☐ <b>RDDS (decay curve, DDCM).</b>" in to_run and "☐ <b>Inelastic scattering.</b>" in to_run


def test_the_report_leads_with_the_results_and_keeps_the_rest_as_appendices(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path / "experiments")
    key = results.key_of_reaction(p.experiment)
    results.append(p.experiment_folder, results.Row(
        key.nuclide, "B(E2)", key, 700.0, "e2fm4", "coulex-yield", "projected", direction="up", stat=20.0, sys=40.0,
        budget={"efficiency": 35.0}, truth=695.0, origin=results.from_run(p.experiment_folder.name, "01-beam-24h")))
    before = {f: f.read_bytes() for f in p.experiment_folder.rglob("*") if f.is_file()}
    page = logbook.experiment_html(p, plan=object())  # a plan without figures: the order is what is tested
    # Building the report reads: the experiment's folder is as it was.
    assert {f: f.read_bytes() for f in p.experiment_folder.rglob("*") if f.is_file()} == before
    at = [page.index(f'<h2 id="{name}">{i}. ') for i, (name, _) in enumerate(logbook.RESULTS, 1)]
    at.append(page.index('<h2 id="appendices"'))
    at += [page.index(f'<h2 id="{name}">{chr(65 + i)}. ') for i, (name, _) in enumerate(logbook.APPENDICES)]
    assert at == sorted(at) and len(at) == 5 + 1 + 8
    lead = page[at[0]:at[1]]
    assert "PROJECTED" in lead and "700 ± 20 (stat) ± 40 (sys) e²fm⁴" in lead and "↑" in lead
    assert "not a measurement" in lead and "closes" in lead
    assert "No run yet" in page[at[5]:]  # the appendices are the report as it was
    # An elastic experiment has no result in the ledger: it opens with what it gives, then the appendices.
    page = logbook.experiment_html(Planner.example("alpha_on_gold"), plan=object())
    assert '<h2 id="result">1. What this setup gives</h2>' in page and '<h2 id="routes">' not in page
    assert page.index('id="result"') < page.index('id="appendices"') < page.index('<h2 id="setup">A. ')


def test_an_analysed_run_without_a_row_says_why():
    """It said "the analysis found no peak" for every analysed run that had no row."""
    p = Planner.example("coulex_ni58")
    nan = float("nan")
    row, why = logbook.analysis_row(p, None, SimpleNamespace(b_e2fm4=nan, truth_e2fm4=700.0))
    assert row is None and "found no peak" in why
    row, why = logbook.analysis_row(p, None, SimpleNamespace(b_e2fm4=690.0, truth_e2fm4=nan))
    assert row is None and "no value that the simulation was given" in why and "peak" not in why
    problems = {}
    assert logbook.experiment_rows(p, [], {}, SimpleNamespace(b_e2fm4=690.0, truth_e2fm4=nan), problems) == []
    assert list(problems) == [None] and "simulation was given" in problems[None]
