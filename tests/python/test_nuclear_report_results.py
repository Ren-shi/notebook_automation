"""The report leads with the results (the results ledger, wave 1): the order of its sections, a projected row
shown with its closure and never pulled against the reference, the comparisons still to run, the numbers'
formatting and the direction of every B(E2). Synthetic rows and setups without runs: seconds."""

import pytest

from physim.nuclear import compare, logbook, results
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
    projected = _rows()[3]
    done = {r["method"]: bool(r["done"]) for r in logbook.routes_to_run("B(E2)", [projected])}
    assert done == {"coulex-yield": True, "coulex-fit": False, "rdds": False, "dsam": False, "other": False}
    page = logbook.to_run_html("B(E2)", [projected], here="62Ni-on-194Pt", unanalysed=[1, 2], outside=True)
    assert "☑ <b>Coulex from the yield (first order).</b>" in page and "Done: this experiment, run 3" in page
    assert "☐ <b>RDDS (decay curve, DDCM).</b>" in page
    assert "physim has this kind of experiment (Lifetimes (RDDS, DSAM))" in page
    assert "☐ <b>Inelastic scattering.</b>" in page and "does not have this kind of experiment yet" in page
    assert "☐ <b>Runs 1 and 2.</b> Taken and not analysed" in page
    for missing in ("No evaluated value", "No measured point", "No shell-model value", "GOSIA"):
        assert missing in page, missing
    # With the evaluation, a measurement and the shell model entered, nothing independent is missing.
    assert all(x["have"] for x in logbook.evidence_missing(_rows()))


def test_the_report_leads_with_the_results_and_keeps_the_rest_as_appendices(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path / "experiments")
    key = results.key_of_reaction(p.experiment)
    results.append(p.experiment_folder, results.Row(
        key.nuclide, "B(E2)", key, 700.0, "e2fm4", "coulex-yield", "projected", direction="up", stat=20.0, sys=40.0,
        budget={"efficiency": 35.0}, truth=695.0, origin=results.from_run(p.experiment_folder.name, "01-beam-24h")))
    page = logbook.experiment_html(p, plan=object())  # a plan without figures: the order is what is tested
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
