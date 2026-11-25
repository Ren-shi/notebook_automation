"""The comparison of one observable with the literature (``physim.nuclear.app.tab_compare``), without a browser, on
temporary folders with a handful of synthetic ⁶²Ni rows: the five tiles, the table grouped by method or ordered by
year or by pull, the unit switch, "What stands out", the window of a data point and a nucleus with almost nothing.
No run is simulated."""

import asyncio
from types import SimpleNamespace

import pytest

from physim.nuclear import Experiment, ensdf, nucleus, results
from physim.nuclear.app import start
from physim.nuclear.runs import ExperimentFolder

KEY = "62Ni:2+1>0+1:E2"
LEVELS = ((0.0, 0, 1, 1), (1172.91, 4, 1, 1), (2336.46, 8, 1, 1))


@pytest.fixture
def lab(tmp_path, monkeypatch):
    """A ⁶²Ni Coulex experiment with one projected row, the evaluated B(E2), a measured point in tension with it that
    the evaluation depends on, a second measured point that agrees, another evaluation and a theory value; no local
    ENSDF copy (the levels to anchor to are here, for ⁶⁰Ni too)."""
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path / "experiments"))
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setattr(results, "ensdf_levels", lambda nuclide: LEVELS if nuclide in ("62Ni", "60Ni") else ())
    monkeypatch.setattr(nucleus, "decay", lambda nuclide, level_key: (None, None, "no γ ray known"))
    root = tmp_path / "experiments"
    coulex = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=root)
    run = coulex.path / "runs" / "03-beam-24h"
    run.mkdir()
    (run / "summary.json").write_text("{}", encoding="utf-8")
    results.append(coulex.path, results.Row(
        "62Ni", "B(E2)", KEY, 918.6, "e2fm4", "coulex-yield", "projected", direction="up", stat=12.9, sys=63.4,
        truth=884.1, origin=results.from_run(coulex.path.name, run.name)))
    doe = nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 980.0, direction="up", uncertainty=37.0,
                                 method="coulex-fit", reference={"authors": "Doe et al.", "year": 2014})
    roe = nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 877.0, direction="up", uncertainty=11.0,
                                 method="electron-scattering", reference={"authors": "Roe et al.", "year": 1967})
    nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 890.0, direction="up", uncertainty=25.0,
                           evidence="evaluated", method="evaluation",
                           reference={"authors": "Raman et al.", "year": 2001})
    nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 1130.0, direction="up", evidence="theory",
                           method="shell-model", reference={"authors": "Poe et al.", "year": 2011})
    results.add_points([results.Row(
        "62Ni", "B(E2)", KEY, 884.1, "e2fm4", "evaluation", "evaluated", direction="up", total=13.1,
        depends_on=(doe.id,), origin=results.reference("ENSDF", 2025, where="NDS 204, 1 (2025)",
                                                       evaluation=nucleus.ENSDF))])
    return SimpleNamespace(coulex=coulex, doe=doe, roe=roe)


def _texts(user) -> list:
    """The text of the labels and buttons, in the order of the page."""
    from nicegui import ui

    found = {e for kind in (ui.label, ui.button) for e in user.find(kind=kind).elements}
    return [str(e.text) for e in sorted(found, key=lambda e: e.id)]


def _on_screen(scenario, **st) -> None:
    """``scenario(user, state)`` on the start page alone (the comparison view), built again as the page's
    ``refresh_start`` does."""
    pytest.importorskip("nicegui")
    from nicegui import ui
    from nicegui.testing import user_simulation

    state = {"start": {"kind": None, "view": "observable", "nucleus": "62Ni", "observable": "B(E2)", "key": KEY,
                       **st}}

    def page() -> None:
        box = ui.column()

        def refresh() -> None:
            box.clear()
            with box:
                start.render(ctx)

        ctx = SimpleNamespace(ui=ui, state=state, open_experiment=lambda path: None, refresh_start=refresh,
                              create_experiment=None, open_example=None, import_setup=None, close_start=None,
                              remember_nucleus=lambda name: None)
        refresh()

    async def main():
        async with user_simulation(page) as user:
            await user.open("/")
            await scenario(user, state["start"])
            await asyncio.sleep(0.1)

    asyncio.run(main())


def _press(user, text: str) -> None:
    from nicegui import ui
    from nicegui.testing import UserInteraction

    found = [e for e in user.find(kind=ui.button).elements if str(e.text).startswith(text)]
    UserInteraction(user, {min(found, key=lambda e: e.id)}, None).click()


def _toggle(user, option: str, value: str) -> None:
    """Set the toggle that offers ``option`` to ``value``, as a click on it does."""
    from nicegui import ui

    [toggle] = [t for t in user.find(kind=ui.toggle).elements if option in t.options]
    toggle.value = value


def _in_table(user, kind: str) -> list:
    """The texts of the table of one kind ("ps-cmp-src": the sources, "ps-cmp-group-label": the groups), in the order
    they stand."""
    from nicegui import ui

    [grid] = [e for e in user.find(kind=ui.element).elements if "ps-cmp-grid" in e.classes]
    return [str(e.text) for e in grid.descendants() if isinstance(e, ui.label) and kind in e.classes]


def _sources(user, names) -> list:
    """The sources of the table's rows, in the order they stand, of those named."""
    return [t for t in _in_table(user, "ps-cmp-src") if t in names]


def _window(user):
    from nicegui import ui

    return [d for d in user.find(kind=ui.dialog).elements if "ps-cmp-window" in
            " ".join(" ".join(c.classes) for c in d.descendants())]


def test_the_tiles_the_groups_and_what_stands_out(lab):
    async def scenario(user, st):
        await user.should_see("The reference", retries=100)
        texts = _texts(user)
        # The head and the tabs.
        assert "B(E2; 0⁺₁ → 2⁺₁)↑" in texts and " of ⁶²Ni, every value side by side" in texts
        assert "Add a data point" in texts and "Import points from a file" in texts
        # The five tiles, from the counts and the reference.
        for tile in ("The reference", "Measured", "Best single measurement", "Theory", "Your setup, projected"):
            assert tile in texts
        assert "884 ± 13" in texts and "ENSDF adopted, NDS 204, 1 (2025) · ±1.5 %" in texts
        assert "2 values" in texts and "1 agrees" in texts and "1 in tension" in texts  # measured
        assert "877 ± 11" in texts and "1 value" in texts  # the best one; theory
        assert "reaches ±7.0 %" in texts and "1 simulated run" in texts
        # The groups in their order, with their rows and verdicts in the data's words.
        groups = ["Evaluated", "Other evaluations", "Measured · Coulomb excitation", "Measured · electron scattering",
                  "Theory", "Projected · simulated runs, analysed back"]
        assert _in_table(user, "ps-cmp-group-label") == groups
        assert "△ tension" in texts and "✓ agree" in texts and "✓ closes" in texts and "+0.53σ" in texts
        assert "◦ not independent of the reference" in texts  # Doe: the evaluation depends on it
        assert "not judged: not independent" in texts  # the other evaluation has no verdict
        assert "not judged: no uncertainty" in texts and "× 1.28" in texts
        assert "Given to the simulation" in texts and "Reference ± its uncertainty" in texts
        # What stands out, and what the projected runs say.
        assert texts.index("What stands out") < texts.index("Doe et al. 2014 is in tension")
        assert "The theory value lies above the reference" in texts
        assert any(t.startswith("This setup reaches ±7.0 %") for t in texts)
        assert "Open 62Ni-on-194Pt →" in texts
        # Under them, today's sections.
        assert "Routes to B(E2; 0⁺₁ → 2⁺₁)↑ in ⁶²Ni" in texts and "History" in texts
        assert "What this value gives" in texts

    _on_screen(scenario)


def test_the_orders_and_the_unit(lab):
    names = ("ENSDF adopted, NDS 204, 1 (2025)", "Doe et al. 2014", "Roe et al. 1967", "Raman et al. 2001",
             "Poe et al. 2011", "62Ni-on-194Pt · 03-beam-24h")

    async def scenario(user, st):
        await user.should_see("The reference", retries=100)
        assert _sources(user, names) == list(names[:1]) + ["Raman et al. 2001", "Doe et al. 2014", "Roe et al. 1967",
                                                            "Poe et al. 2011", names[-1]]
        _toggle(user, "year", "year")  # one list, the reference first and the projected run last
        await user.should_see("By year", retries=100)
        await asyncio.sleep(0.1)
        by_year = _sources(user, names)
        assert st["order"] == "year" and by_year[0] == names[0] and by_year[-1] == names[-1]
        assert not _in_table(user, "ps-cmp-group-label")
        _toggle(user, "pull", "pull")
        await asyncio.sleep(0.1)
        by_pull = _sources(user, names)
        assert by_pull[0] == names[0] and by_pull[-1] == names[-1] and by_pull != by_year
        assert by_pull.index("Doe et al. 2014") < by_pull.index("Roe et al. 1967")  # the larger pull first
        # The unit: every value, the reference tile and the axis in e²b².
        _toggle(user, "e2b2", "e2b2")
        await user.should_see("0.0884 ± 0.0013", retries=100)
        texts = _texts(user)
        assert st["unit"] == "e2b2" and "0.0980 ± 0.0037" in texts
        assert "Value (e²b²↑)" in texts

    _on_screen(scenario)


def test_the_window_of_a_data_point(lab):
    async def scenario(user, st):
        await user.should_see("Add the point", retries=100)
        [window] = _window(user)
        assert window.value is True  # asked for by the landing page: open at once
        assert "to_point_form" not in st
        st["point"].update(value="0.0906", unit="e2b2", uncertainty="0.0037", method="coulex-fit", authors="Allmond",
                           year="2014")
        _press(user, "Add the point")
        await user.should_see("906 ± 37", retries=100)
        [window] = _window(user)
        assert window.value is False  # the page is read again, the window closed
        # Correct, from the row's menu: the window opens with the point in it.
        from nicegui import ui
        from nicegui.testing import UserInteraction

        def in_row_of(item, source: str) -> bool:
            row = next((a for a in item.ancestors() if "ps-cmp-row" in a.classes), None)
            return row is not None and any(getattr(x, "text", None) == source for x in row.descendants())

        [item] = [m for m in user.find(kind=ui.menu_item).elements if in_row_of(m, "Roe et al. 1967")
                  and any(getattr(x, "text", None) == "Correct…" for x in m.descendants())]
        UserInteraction(user, {item}, None).click()
        await user.should_see("Correct a data point", retries=100)
        [window] = _window(user)
        assert window.value is True and st["point"]["corrects"] == lab.roe.id

    _on_screen(scenario, to_point_form=True)


def test_a_nucleus_with_almost_nothing(lab):
    # ⁶⁰Ni: one calculated τ without an uncertainty, no evaluation, no measurement, no run.
    nucleus.add_data_point("60Ni", "tau", "2+1", 1.0, evidence="theory", method="shell-model",
                           reference={"authors": "Poe et al.", "year": 2011})

    async def scenario(user, st):
        await user.should_see("The reference", retries=100)
        texts = _texts(user)
        assert "τ(2⁺₁)" in texts and "Poe et al. 2011" in texts and "1 value" in texts
        assert "nothing evaluated yet: pick a row to compare against" in texts
        assert "no measured value" in texts and "no simulated run" in texts
        assert any(t.startswith("Nothing stands out") for t in texts)
        assert any(t.startswith("No simulated run of this is analysed yet") for t in texts)
        assert "Reference ± its uncertainty" not in texts and "no reference" in texts

    _on_screen(scenario, nucleus="60Ni", observable="tau", key="60Ni:2+1")
