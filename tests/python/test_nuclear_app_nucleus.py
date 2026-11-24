"""The nuclei screen and the page of one nucleus (``physim.nuclear.app.start`` and ``.tab_nucleus``) without a
browser, on temporary folders with a handful of synthetic ⁶²Ni rows. No run is simulated."""

import asyncio
from types import SimpleNamespace

import pytest

from physim.nuclear import Experiment, ensdf, nucleus, results
from physim.nuclear.app import start
from physim.nuclear.app import tab_nucleus as tn
from physim.nuclear.runs import ExperimentFolder

KEY = "62Ni:2+1>0+1:E2"
NI62 = ((0.0, 0, 1, 1), (1172.91, 4, 1, 1), (2336.46, 8, 1, 1))


@pytest.fixture
def lab(tmp_path, monkeypatch):
    """A ⁶²Ni Coulex experiment with one projected row, an elastic one, the evaluated B(E2), and a measured point in
    tension with it that the evaluation depends on; no local ENSDF copy (the levels to anchor to are here)."""
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path / "experiments"))
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setattr(results, "ensdf_levels", lambda nuclide: NI62 if nuclide == "62Ni" else ())
    monkeypatch.setattr(nucleus, "decay", lambda nuclide, level_key: (None, None, "no γ ray known"))
    root = tmp_path / "experiments"
    coulex = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=root)
    ExperimentFolder.create(Experiment.example("alpha_on_gold"), root=root)
    run = coulex.path / "runs" / "03-beam-24h"
    run.mkdir()
    (run / "summary.json").write_text("{}", encoding="utf-8")
    results.append(coulex.path, results.Row(
        "62Ni", "B(E2)", KEY, 918.6, "e2fm4", "coulex-yield", "projected", direction="up", stat=12.9, sys=63.4,
        truth=884.1, origin=results.from_run(coulex.path.name, run.name)))
    doe = nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 980.0, direction="up", uncertainty=37.0,
                                 method="coulex-fit", reference={"authors": "Doe et al.", "year": 2014})
    results.add_points([results.Row(
        "62Ni", "B(E2)", KEY, 884.1, "e2fm4", "evaluation", "evaluated", direction="up", total=13.1,
        depends_on=(doe.id,), origin=results.reference("ENSDF", 2025, where="NDS 204, 1 (2025)",
                                                       evaluation=nucleus.ENSDF))])
    return SimpleNamespace(coulex=coulex, doe=doe)


def _texts(user) -> str:
    """The text of the labels and buttons, in the order of the page."""
    from nicegui import ui

    found = {e for kind in (ui.label, ui.button) for e in user.find(kind=kind).elements}
    return "\n".join(str(e.text) for e in sorted(found, key=lambda e: e.id))


def _on_screen(scenario, **st) -> list:
    """``scenario(user, state)`` on the start page alone, built again as the page's ``refresh_start`` does; the
    experiments it opened."""
    pytest.importorskip("nicegui")
    from nicegui import ui
    from nicegui.testing import user_simulation

    opened, state = [], {"start": {"kind": None, **st}}

    def page() -> None:
        box = ui.column()

        def refresh() -> None:
            box.clear()
            with box:
                start.render(ctx)

        ctx = SimpleNamespace(ui=ui, state=state, open_experiment=opened.append, refresh_start=refresh,
                              create_experiment=None, open_example=None, import_setup=None, close_start=None)
        refresh()

    async def main():
        async with user_simulation(page) as user:
            await user.open("/")
            await scenario(user, state["start"])
            await asyncio.sleep(0.1)

    asyncio.run(main())
    return opened


def _press(user, text: str) -> None:
    from nicegui import ui
    from nicegui.testing import UserInteraction

    found = [e for e in user.find(kind=ui.button).elements if str(e.text).startswith(text)]
    UserInteraction(user, {min(found, key=lambda e: e.id)}, None).click()


def test_the_nuclei_screen_and_the_page_of_one_nucleus(lab):
    # The evaluation of the partner is here (its page was opened once): it stays in the line under the cards.
    results.add_points([results.Row("194Pt", "E_level", "194Pt:2+1", 328.5, "keV", "evaluation", "evaluated",
                                    total=0.1, origin=results.reference("ENSDF", 2025, evaluation=nucleus.ENSDF))])

    async def scenario(user, st):
        from nicegui import ui
        from nicegui.testing import UserInteraction

        # The nuclei come first: a card for the nucleus studied, the partner in a line, the unassigned in a table.
        await user.should_see("WHAT IS KNOWN AND PROJECTED", retries=100)
        texts = _texts(user)
        assert "⁶²Ni →" in texts and "62Ni on 194Pt" in texts and "Experiments tied to no nucleus" in texts
        assert ("884 ± 13 e²fm⁴↑ evaluated · 980 ± 37 e²fm⁴↑ measured (Coulex fitted to all orders) · projected: "
                "reaches ±7.0 %") in texts
        assert "¹⁹⁴Pt · in 1 experiment" in texts and "¹⁹⁴Pt →" not in texts
        [unassigned] = [r for t in user.find(kind=ui.table).elements for r in t.rows]
        assert unassigned["reason"] == "elastic scattering studies no level of a nucleus"
        [row] = [e for e in user.find(kind=ui.element).elements if "ps-nt-click" in e.classes]
        UserInteraction(user, {row}, None).click()  # a row of the card opens the experiment
        _press(user, "⁶²Ni →")
        # One nucleus: the two groups apart, each with its own comparison and the method in words.
        await user.should_see("Given to the simulation", retries=100)
        texts = _texts(user)
        assert "INDEPENDENT · against the reference: ENSDF adopted" in texts
        assert "△ tension" in texts and "◦ not independent of the reference" in texts  # the measured point
        assert "✓ closes (+0.53σ)" in texts and "884.1" in texts
        assert texts.index("△ tension") < texts.index("PROJECTED · simulated runs") < texts.index("✓ closes")
        assert "This setup reaches ±7.0 %" in texts and "the value is known to ±1.5 %" in texts
        # The evaluation is here already, so the page's own note of a missing ENSDF copy is not (a route may
        # still say that it could not test itself against ENSDF: that is another sentence).
        assert "There is no local copy of ENSDF, so the adopted values" not in texts
        # Plan it: the new-experiment form of that kind, with the nucleus as the one that is excited.
        _press(user, "Plan it for ⁶²Ni")
        await user.should_see("62Ni is filled in as the beam", retries=100)
        assert st["kind"] in ("coulex", "lifetimes") and st["values"]["beam"] == "62Ni"

    assert _on_screen(scenario) == [str(lab.coulex.path)]


def test_a_data_point_through_the_form(lab):
    form = tn.blank_point("B(E2)", KEY)
    assert (form["key"], form["unit"], form["direction"]) == ("2+1 > 0+1", "e2fm4", "up")

    async def scenario(user, st):
        await user.should_see("Add the point", retries=100)
        st["point"].update(value="many", method="dsam", authors="Roe", year="2020")
        _press(user, "Add the point")  # a bad point: everything that is wrong, beside the form
        await user.should_see("the value must be a number, got 'many'")
        await user.should_see("a measured point needs its uncertainty")
        st["point"].update(value="0.0906", unit="e2b2", uncertainty="0.0037")
        _press(user, "Add the point")
        await user.should_see("906 ± 37", retries=100)
        texts = _texts(user)
        assert "906 ± 37" in texts and "✓ agree" in texts and "2 data points" in texts

    _on_screen(scenario, view="nucleus", nucleus="62Ni", observable="B(E2)", key=KEY.replace("+", " "))
    # Corrected through the same form: the old row is superseded, never changed.
    assert tn.add_point("62Ni", {**tn.point_of(lab.doe), "value": "906"}).supersedes == lab.doe.id
    # The form's other shapes of an uncertainty and of where the point is from.
    args = tn.point_arguments({**form, "value": "900", "form": "asymmetric", "plus": "60", "minus": "50",
                               "from": "own", "date": "2027-06-14", "note": "runs 12 to 31", "method": "coulex-fit"})
    assert (args["plus"], args["minus"]) == (60.0, 50.0) and "uncertainty" not in args and "reference" not in args
    assert args["own"] == {"date": "2027-06-14", "note": "runs 12 to 31"}


def test_words_units_addresses_and_the_figure(lab, monkeypatch):
    assert tn.url("62Ni", "B(E2)", KEY) == "/?nucleus=62Ni&observable=B%28E2%29&key=62Ni%3A2%2B1%3E0%2B1%3AE2"
    assert tn.shown(12.13, 0.18, 0.18, "W.u.") == "12.13 ± 0.18 W.u.↓"  # always marked ↓
    factor, unit = tn.display_factor("B(E2)", "e2b2")
    assert tn.shown(884.1, 13.1, 13.1, unit, "up", factor) == "0.0884 ± 0.0013 e²b²↑"
    assert tn.display_factor("tau", "e2b2") == (1.0, "ps")
    figure = tn.figure_values(nucleus.observable_view("62Ni", "B(E2)", KEY), "dark", "e2b2")
    assert [t.name for t in figure.data] == ["evaluated", "measured", "projected (simulated runs)",
                                             "given to the simulation"]
    # What a value is used for: one string, or texts with their references.
    entry = SimpleNamespace(used_for=[{"text": "collectivity", "reference": "a book"}, {"text": "a text"}])
    monkeypatch.setattr(tn.observables, "get", lambda key: entry)
    assert tn.uses("B(E2)") == [{"text": "collectivity", "reference": "a book"}, {"text": "a text", "reference": ""}]
    entry.used_for = "one line"
    assert tn.uses("B(E2)") == [{"text": "one line", "reference": ""}]


def test_an_experiment_in_the_address_that_is_not_there(lab, tmp_path):
    from physim.nuclear.app import page

    assert page.missing_experiment(None) is None and page.missing_experiment(str(lab.coulex.path)) is None
    said = page.missing_experiment(str(tmp_path / "no-such-experiment"))
    assert "there is no experiment folder there" in said

    async def scenario(user, st):
        await user.should_see("there is no experiment folder there", retries=100)
        await user.should_see("What are you studying?")
        assert "note" not in st  # said once

    _on_screen(scenario, view="nuclei", note=said)


def test_an_import_and_the_report_from_the_page(lab, tmp_path):
    data = (b"observable,key,value,unit,uncertainty,method,authors,year\n"
            b"tau,2+1,2.01,ps,0.12,dsam,Kenn et al.,2001\ntau,2+1,2.2,ps,,,Roe,2020\n")
    said = tn.import_lines(tn.import_file("62Ni", "mine.csv", data))
    assert said["added"] == ["line 2 · tau of 2+1: 2.01 ± 0.12 ps, dsam, Kenn et al. 2001 → DSAM (line shape)"]
    assert not said["already"] and "pick the method" in said["refused"][0]
    assert tn.import_lines(tn.import_file("62Ni", "mine.csv", data))["already"][0].startswith("line 2 · tau")
    assert "B(E2; 0⁺₁ → 2⁺₁)↑" in tn.report_html("62Ni")


def test_the_b_values_of_a_level_scheme_carry_their_arrows():
    from physim.nuclear.app.page import b_words

    table = {"nuclide": "62Ni", "levels": [{"energy_kev": 0.0, "jpi": "0+"}, {"energy_kev": 1172.91, "jpi": "2+"},
                                           {"energy_kev": 2000.0, "jpi": ""}]}
    m = {"from": 0, "to": 1, "multipolarity": "E2", "b_up_e2fm": 884.1}
    assert b_words(table, m) == "B↑ = 884.1 e²fm^4 · B↓ = 176.8 e²fm^4 = 12.1 W.u. ↓"
    assert b_words(table, {**m, "to": 2}).endswith("(B↓ and W.u. need the spins of both levels)")
