"""The deorientation in vacuum in the app (wave 3B): the Setup card's two fields, greyed out beside a lifetime
measurement, the Physics page's part, the record's row, the report's sentence and the words of a run's changes.

The pages are NiceGUI's simulated browser, with the app's folders under the test's own temporary folder, so nothing
reads the machine's ~/.physim."""

from __future__ import annotations

import asyncio

import pytest

from physim.nuclear import record
from physim.nuclear import workbench as wb
from physim.nuclear.planner import Planner


@pytest.fixture(autouse=True)
def own_folders(tmp_path, monkeypatch):
    for name in ("experiments", "shellmodel", "nuclei"):
        monkeypatch.setenv("PHYSIM_" + name.upper(), str(tmp_path / name))


# -- what the card and the pages are made of ------------------------------------------------------------------------


def test_the_card_state_follows_the_setup():
    p = Planner.example("coulex_ni62_pt194")
    st = wb.deorientation_state(p.draft)
    assert (st["set"], st["g"], st["in_effect"], st["lifetime"], st["ignored"]) == (None, None, "vacuum", None, [])
    assert st["nuclide"] == "62Ni" and st["z_over_a"] == pytest.approx(28 / 62)
    assert "1.5 to 3%" in wb.deorientation_note(st) and "20 to 23%" in wb.deorientation_note(st)
    raw = p.draft
    raw["reaction"].update(deorientation="none", deorientation_g=0.3)
    st = wb.deorientation_state(raw)
    assert st["in_effect"] == "none" and st["ignored"] == ["deorientation_g"]
    assert wb.deorientation_note(st).startswith("Left out") and "not used" in wb.deorientation_note(st)
    # A natural target takes its most abundant isotope, as the excitation does; a compound has no Z/A to show.
    raw = {"beam": {"nuclide": "16O"}, "target": {"material": "Pt"}, "reaction": {"type": "coulex"}}
    assert wb.deorientation_state(raw)["nuclide"] == "195Pt"
    raw["target"]["material"] = "CD2"
    assert wb.deorientation_state(raw)["z_over_a"] is None


def test_a_lifetime_measurement_has_its_own():
    p = Planner.example("coulex_ru98_plunger")
    st = wb.deorientation_state(p.draft)
    assert st["in_effect"] == "none" and st["lifetime"] == "rdds" and st["ignored"] == []
    assert wb.deorientation_note(st).startswith("Not in effect: a plunger has its own")
    assert not any("does nothing" in w.text for w in p.warnings())
    assert p.set("reaction", "deorientation", "vacuum") and p.set("reaction", "deorientation_g", 0.4)
    st = wb.deorientation_state(p.draft)
    assert st["ignored"] == ["deorientation", "deorientation_g"] and "does nothing" in wb.deorientation_note(st)
    notes = [w for w in p.warnings() if "does nothing beside the lifetime measurement" in w.text]
    assert len(notes) == 1 and notes[0].level == "note"
    assert "reaction.deorientation and reaction.deorientation_g" in notes[0].text
    assert record.deorientation_row(p.experiment) == (
        "Deorientation in vacuum", "not applied (lifetime measurement: the plunger's own)", "setup")
    dsam = {"reaction": {"type": "coulex", "deorientation": "vacuum"}, "lifetime": {"method": "dsam"}}
    assert wb.deorientation_note(wb.deorientation_state(dsam)).startswith("Not in effect: in a DSAM")


def test_the_physics_page_levels_of_both_nuclei():
    p = Planner.example("coulex_ni62_pt194")
    d = wb.deorientation_levels(p)
    assert d["in_effect"] == "vacuum" and d["g_source"] == "Z/A" and d["g"] == pytest.approx(28 / 62)
    by = {(r["nucleus"], r["energy_kev"]): r for r in d["levels"]}
    assert {n for n, _ in by} == {"62Ni", "194Pt"} and len(by) == 8
    pt = by[("194Pt", 328.473)]["g"]
    ni = by[("62Ni", 1172.91)]["g"]
    # G_k only falls, faster with speed and life, faster for the higher rank; the 0⁺ levels have none.
    assert 1 > ni[0.02][0] > ni[0.05][0] > ni[0.08][0] > pt[0.08][0] > 0 and pt[0.05][1] < pt[0.05][0]
    assert pt[0.05][0] == pytest.approx(0.155, abs=0.002)  # GOSIA prints 0.1509 at β = 0.0491 (gosia_gk.csv)
    assert by[("194Pt", 1411.83)]["g"][0.05][1] is not None
    p.set("reaction", "deorientation_g", 0.2)
    q = wb.deorientation_levels(Planner(p.experiment))
    assert q["g"] == 0.2 and q["g_source"] == "setup"
    # The other nucleus keeps Z/A; the excited one deorients less with the smaller g.
    after = {(r["nucleus"], r["energy_kev"]): r["g"] for r in q["levels"]}
    assert after[("194Pt", 328.473)] == pt and after[("62Ni", 1172.91)][0.05][0] > ni[0.05][0]
    p.set("reaction", "deorientation", "none")
    q = wb.deorientation_levels(Planner(p.experiment))
    assert q["in_effect"] == "none" and q["levels"] and not any("g" in r for r in q["levels"])
    for name in Planner.examples():  # never an error, whatever the setup
        wb.deorientation_levels(Planner.example(name))


def test_the_record_row_and_the_report_sentence():
    p = Planner.example("coulex_ni62_pt194")
    assert record.deorientation_row(p.experiment) == (
        "Deorientation in vacuum", "GOSIA's two-state model, its parameters, g = Z/A = 0.452", "GOSIA's default")
    assert record.deorientation_row(p.experiment) in record.provenance(p.experiment)
    assert any("transient field" in x for x in record.LEFT_OUT)
    p.set("reaction", "deorientation_g", 0.25)
    assert record.deorientation_row(p.experiment)[1:] == ("GOSIA's two-state model, its parameters, g = 0.25", "user")
    p.set("reaction", "deorientation", "none")
    assert record.deorientation_row(p.experiment) == ("Deorientation in vacuum", "none", "user")
    assert wb.decisions(p)[1].state.endswith(", no deorientation")


def _report_section(p) -> str:
    """The report's Coulomb-excitation section of a planner's setup, without building the whole report."""
    from physim.nuclear.report import Report

    r = Report.__new__(Report)
    r.excitation, r.planner = p.gamma(), p
    return r._excitation_html(schemes=False)


def test_the_report_says_when_it_is_applied():
    p = Planner.example("coulex_ni58")
    section = _report_section(p)
    assert "is deoriented in vacuum (GOSIA's two-state model with its parameters, g = Z/A)" in section
    p.set("reaction", "deorientation_g", 0.5)
    assert "g = 0.5)" in _report_section(Planner(p.experiment))
    p.set("reaction", "deorientation", "none")
    assert "deoriented" not in _report_section(Planner(p.experiment))


def test_a_runs_changes_in_words():
    assert wb.change_words("reaction.deorientation: none → vacuum") == "deorientation in vacuum: none → vacuum"
    assert wb.change_words("reaction.e1_polarisation: 0.0 → 0.0038") == "E1 polarisation (GOSIA's d): 0.0 → 0.0038"
    assert wb.change_words("target.thickness: 0.5 mg/cm2 → 1.0 mg/cm2") == "target.thickness: 0.5 mg/cm2 → 1.0 mg/cm2"


# -- the page -------------------------------------------------------------------------------------------------------


def _on_page(planner, scenario, prebuild=False) -> None:
    """``scenario(user)`` on the planner page of ``planner`` (taken as the example it was made from)."""
    pytest.importorskip("nicegui")
    from nicegui.testing import user_simulation

    from physim.nuclear.app.page import build_page

    async def main():
        async with user_simulation(lambda: build_page(example="x", prebuild=prebuild)) as user:
            await user.open("/")
            await scenario(user)
            await asyncio.sleep(0.2)

    original = Planner.example
    Planner.example = classmethod(lambda cls, name: planner)
    try:
        asyncio.run(main())
    finally:
        Planner.example = original


def _texts(user) -> str:
    from nicegui import ui

    return " | ".join(str(e.text) for e in user.find(kind=ui.label).elements if e.visible)


def _elements(user, kind, label):
    try:
        found = user.find(kind=kind).elements
    except AssertionError:
        return []
    return [e for e in found if e.props.get("label") == label]


async def _card(user, title: str):
    """Open the setup panel's card ``title`` and wait for its deorientation field."""
    from nicegui import ui

    await user.should_see("Run conditions", retries=600)
    card = next(e for e in user.find(kind=ui.expansion).elements if "header" in e.slots and any(
        getattr(c, "text", None) == title for h in e.slots["header"].children for c in h.descendants()))
    card.value = True
    for _ in range(200):
        await asyncio.sleep(0.05)
        if _elements(user, ui.select, "Deorientation in vacuum"):
            break


def _select(user):
    from nicegui import ui

    return _elements(user, ui.select, "Deorientation in vacuum")[0]


def _g(user):
    from nicegui import ui

    found = user.find(kind=ui.input).elements
    return next(e for e in found if str(e.props.get("label", "")).startswith("g factor (empty: "))


def _press(user, label: str) -> None:
    from nicegui import ui
    from nicegui.testing import UserInteraction

    button = next(e for e in user.find(kind=ui.button).elements if e.text == label)
    UserInteraction(user, {button}, None).click()


async def _applied(user, p, before) -> None:
    _press(user, "Apply")
    for _ in range(600):
        if p.experiment is not before:
            return
        await asyncio.sleep(0.05)
    raise AssertionError("Apply did not take the change")


def test_the_card_sets_the_deorientation_and_its_g_factor():
    from nicegui.testing import UserInteraction

    p = Planner.example("coulex_ni58")

    async def scenario(user):
        await _card(user, "Target and reaction")
        sel, g = _select(user), _g(user)
        assert sel.value == "vacuum" and "disable" not in sel.props and "disable" not in g.props
        assert g.value == "" and g.props.get("label") == "g factor (empty: Z/A = 0.483)"
        assert "GOSIA's two-state model" in _texts(user) and "1.5 to 3%" in _texts(user)
        # None, applied: the setup says so, and the field for g is greyed out.
        before = p.experiment
        sel.value = "none"
        await asyncio.sleep(0.2)
        assert "1 change not applied" in _texts(user)
        await _applied(user, p, before)
        assert p.experiment.excitation.deorientation == "none"
        assert 'deorientation = "none"' in p.to_toml()
        assert "disable" in _g(user).props and "Left out" in _texts(user)
        # Back to in vacuum: the setup says nothing again (the default).
        before = p.experiment
        _select(user).value = "vacuum"
        await asyncio.sleep(0.2)
        await _applied(user, p, before)
        assert p.experiment.excitation.deorientation is None and "deorientation" not in p.to_toml()
        # A g factor typed in, and taken out again.
        before, box = p.experiment, _g(user)
        box.value = "0.33"
        UserInteraction(user, {box}, None).trigger("blur")
        await asyncio.sleep(0.2)
        await _applied(user, p, before)
        assert p.experiment.excitation.deorientation_g == pytest.approx(0.33)
        before, box = p.experiment, _g(user)
        box.value = ""
        UserInteraction(user, {box}, None).trigger("blur")
        await asyncio.sleep(0.2)
        await _applied(user, p, before)
        assert p.experiment.excitation.deorientation_g is None

    _on_page(p, scenario)


def test_beside_a_plunger_the_fields_are_greyed_out_and_what_is_written_can_go():
    p = Planner.example("coulex_ru98_plunger")
    assert p.set("reaction", "deorientation", "vacuum")
    p = Planner(p.experiment)

    async def scenario(user):
        await _card(user, "Target and reaction")
        sel, g = _select(user), _g(user)
        assert "disable" in sel.props and "disable" in g.props and sel.value == "none"
        texts = _texts(user)
        assert "Not in effect: a plunger has its own deorientation" in texts and "does nothing" in texts
        # The Physics page says why too.
        assert "Explain: Deorientation in vacuum" in " ".join(
            str(e.props.get("label")) for e in user.find(kind=__import__("nicegui").ui.expansion).elements)
        before = p.experiment
        _press(user, "Remove it")
        await asyncio.sleep(0.2)
        await _applied(user, p, before)
        assert p.experiment.excitation.deorientation is None
        assert not any("does nothing" in w.text for w in p.warnings())

    _on_page(p, scenario, prebuild=True)
