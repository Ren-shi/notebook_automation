"""The Shell model panel of the planner app (``physim.nuclear.app.tab_shellmodel``): its helpers, and the panel on
the page of the Coulomb-excitation examples, with and without a level scheme. No KSHELL is needed: a fake install
stands in for it, and a result written into the job cache stands in for a run."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from physim.nuclear import shellmodel as sm
from physim.nuclear import shellmodel_jobs as jobs
from physim.nuclear import shellmodel_spaces as spaces
from physim.nuclear.app import tab_shellmodel as tsm
from physim.nuclear.planner import Planner

FAKE_VERSION = "fake KSHELL for the app tests\nbuilt: never\n"

#: The default request of the panel for the ⁶²Ni beam of coulex_ni62_pt194 (JUN45, the smallest space).
NI62_REQUEST = sm.ShellModelRequest("62Ni", "jun45", ((0, 1, 3), (4, 1, 3), (8, 1, 2)), ("E2", "M1"))


def ni62_result() -> sm.ShellModelResult:
    """A made-up ⁶²Ni result: five states that pair with the example's ENSDF levels and one (2⁺₃) that does not."""
    def state(j2, index, ex, q=None, mu=None):
        return sm.ShellState(j2, 1, index, -40.0 + ex, ex, q_s=q, mu=mu)

    states = [state(0, 1, 0.0), state(4, 1, 1.10, q=-12.0, mu=0.6), state(0, 2, 2.00),
              state(4, 2, 2.20, q=9.0, mu=0.9), state(8, 1, 2.40, q=-15.0, mu=2.1), state(4, 3, 2.95, q=3.0)]

    def tr(kind, i, f, mred, ji):
        return sm.ShellTransition(kind, (i[0], 1, i[1]), (f[0], 1, f[1]), mred, mred ** 2 / (ji + 1),
                                  e_initial_mev=None)

    transitions = [tr("E2", (4, 1), (0, 1), -25.0, 4), tr("E2", (8, 1), (4, 1), -40.0, 8),
                   tr("E2", (4, 2), (0, 1), 5.0, 4), tr("E2", (4, 2), (4, 1), 20.0, 4),
                   tr("M1", (4, 2), (4, 1), 0.3, 4)]
    return sm.ShellModelResult(NI62_REQUEST, "f5pg9", "jun45", "fake", (1.5, 1.1), (1.0, 0.0), (3.91, -2.678),
                               states, transitions, {(0, 1): 3103}, 12.0, [], "62Ni")


def fake_install(folder: Path):
    return SimpleNamespace(version=FAKE_VERSION, folder=folder)


def put_in_cache(root: Path, request: sm.ShellModelRequest, result: sm.ShellModelResult) -> str:
    """Write ``result`` where the job cache keeps the result of ``request`` run with the fake build."""
    key = jobs.job_key(request, FAKE_VERSION)
    folder = root / key
    folder.mkdir(parents=True)
    (folder / "results.json").write_text(json.dumps(dict(result.to_dict(), job=key)), encoding="utf-8")
    return key


# -- helpers ----------------------------------------------------------------------------------------------------------

def test_states_are_read_and_written_as_text():
    assert tsm.parse_states("0+ ×3, 2+ x3, 4+*2") == ((0, 1, 3), (4, 1, 3), (8, 1, 2))
    assert tsm.parse_states("3/2- ×2; 1/2-") == ((3, -1, 2), (1, -1, 1))
    assert tsm.states_text(((0, 1, 3), (4, 1, 3), (3, -1, 2))) == "0+ ×3, 2+ ×3, 3/2- ×2"
    for bad in ("", "two plus", "2 ×3"):
        with pytest.raises(ValueError):
            tsm.parse_states(bad)


def test_the_first_choice_is_the_smallest_calculation_with_sensible_states():
    ni58 = spaces.options("58Ni")
    assert tsm.default_option(ni58).interaction == "jun45"  # 2 neutrons above ⁵⁶Ni, not the full fp shell
    si28 = tsm.default_option(spaces.options("28Si"))
    assert si28.interaction == "usdb"  # USDA and USDB are the same size: the usual one first
    assert tsm.default_states(si28) == ((0, 1, 3), (4, 1, 3), (8, 1, 2))  # even-even: 3×0⁺, 3×2⁺, 2×4⁺
    assert tsm.default_states(tsm.default_option(spaces.options("16O"))) == ((0, 1, 1),)  # only one state exists
    odd = tsm.default_states(tsm.default_option(spaces.options("27Al")))
    assert {j2 for j2, _, _ in odd} == {1, 3, 5, 7} and all(p == 1 for _, p, _ in odd)
    assert tsm.default_option(spaces.options("194Pt")) is None  # outside every space


def test_what_to_apply_and_its_preview_keep_ensdf_unless_replace_is_ticked():
    assert tsm.which_of({"energy": True, "E2": True, "Q": True, "new": False, "states": []}) == ["energy", "E2", "Q"]
    assert tsm.which_of({"energy": False, "E2": True, "Q": False, "new": True, "states": ["2+1"]}) == [
        "E2", "new", "2+1"]
    assert tsm.which_of({"energy": False, "E2": False, "Q": False, "new": True}) == ["energy", "new"]
    assert tsm.which_of({"energy": False, "E2": False, "Q": False, "new": False, "states": ["2+1"]}) == []
    p = Planner.example("coulex_ni62_pt194")
    scheme = p.experiment.levels["beam"]
    before = scheme.to_dict()
    kept = tsm.preview(ni62_result(), scheme, ["energy", "E2", "Q"], replace=False)
    assert scheme.to_dict() == before  # a preview changes nothing
    assert any(line.startswith("kept ") and "Replace existing values" in line for line in kept)
    assert not any("energy 1172.91 keV" in line for line in kept)  # the ENSDF energy stays
    replaced = tsm.preview(ni62_result(), scheme, ["energy", "E2", "Q"], replace=True)
    assert any("energy 1172.91 keV (ensdf) → 1100.0 keV" in line for line in replaced)
    assert not any(line.startswith("kept ") for line in replaced)
    again = tsm.preview(ni62_result(), scheme, ["E2"], replace=False)  # every E2 is ENSDF-derived and kept
    assert again[0].startswith("Nothing would change") and tsm.kept(again[-1])
    assert "starts one" in tsm.preview(ni62_result(), None, "all", False)[0]


def test_elapsed_time_and_the_cache_lookup(tmp_path, monkeypatch):
    st = {"started": "2026-11-18T10:00:00Z", "finished": "2026-11-18T10:01:30Z"}
    assert tsm.elapsed(st) == pytest.approx(90.0)
    assert tsm.elapsed({"started": None}) is None
    monkeypatch.setenv("PHYSIM_SHELLMODEL", str(tmp_path))
    install = fake_install(tmp_path)
    assert tsm.cached_job(NI62_REQUEST, install) is None and tsm.cached_job(NI62_REQUEST, None) is None
    key = put_in_cache(tmp_path, NI62_REQUEST, ni62_result())
    assert tsm.cached_job(NI62_REQUEST, install) == key
    assert jobs.result(key).states[1].ex_mev == pytest.approx(1.10)


# -- on the page ------------------------------------------------------------------------------------------------------

def _expansion(user, nuclide: str):
    from nicegui import ui

    label = f"Shell model: calculate or import the levels of {nuclide}"
    found = [e for e in user.find(kind=ui.expansion).elements if e.props.get("label") == label]
    return found[0] if found else None


async def _open_panel(user, nuclide: str):
    from test_nuclear_app import asyncio_sleep

    for _ in range(600):  # the Physics tab, which holds the level schemes, is built ahead of time
        exp = _expansion(user, nuclide)
        if exp is not None:
            break
        await asyncio_sleep(0.05)
    assert exp is not None, nuclide
    exp.value = True
    await asyncio_sleep(0.2)
    return exp


def _button(user, text: str):
    """The button ``text`` of a Shell model panel (the page has other buttons of the same name, e.g. Cancel)."""
    from nicegui import ui

    def in_panel(e) -> bool:
        return any(isinstance(a, ui.expansion) and str(a.props.get("label", "")).startswith("Shell model:")
                   for a in e.ancestors())

    return next(e for e in user.find(kind=ui.button).elements if e.text == text and in_panel(e))


def _click(user, text: str) -> None:
    from nicegui.testing import UserInteraction

    UserInteraction(user, {_button(user, text)}, None).click()


def _checkbox(user, text: str):
    from nicegui import ui

    return next(e for e in user.find(kind=ui.checkbox).elements if e.text == text)


def test_the_panel_without_a_level_scheme_estimates_and_gates_the_run(tmp_path, monkeypatch):
    pytest.importorskip("nicegui")
    from nicegui import ui

    from test_nuclear_app import _page, _texts, asyncio_sleep

    monkeypatch.setenv("PHYSIM_SHELLMODEL", str(tmp_path / "cache"))
    monkeypatch.setattr(tsm, "_install", lambda: fake_install(tmp_path))

    async def scenario(user, p, on_loop):
        await user.open("/")
        await user.should_see("Beam time needed", retries=600)
        assert p.experiment.levels == {}  # coulex_ni58 has no level scheme
        await _open_panel(user, "58Ni")
        await user.should_see("Before running", retries=100)
        texts = _texts(user)
        assert "JUN45" in " ".join(str(e.options) for e in user.find(kind=ui.select).elements)
        assert "of RAM free now" in texts and "Rough time" in texts and "Fits" in texts
        run = _button(user, "Run KSHELL")
        assert run.enabled
        # The full fp shell is far too big: refused unless the override is ticked.
        select = next(e for e in user.find(kind=ui.select).elements if e.props.get("label") == "Interaction")
        select.value = "kb3g"
        await asyncio_sleep(0.2)
        assert "Too big for this computer now" in _texts(user)
        assert not _button(user, "Run KSHELL").enabled
        _checkbox(user, "Run anyway, over the limit").value = True
        await asyncio_sleep(0.1)
        assert _button(user, "Run KSHELL").enabled
        # A state list that is not one says so, and nothing can run.
        states = next(e for e in user.find(kind=ui.input).elements if e.props.get("label") == "States (J± ×how many)")
        from nicegui.testing import UserInteraction

        states.value = "two plus"
        UserInteraction(user, {states}, None).trigger("blur")
        await asyncio_sleep(0.1)
        assert "Not a calculation yet" in _texts(user) and not _button(user, "Run KSHELL").enabled
        # The beam, ¹⁶O, has its own panel.
        assert _expansion(user, "16O") is not None

    _page("coulex_ni58", scenario)


def test_the_panel_without_kshell_says_where_to_point_it(tmp_path, monkeypatch):
    pytest.importorskip("nicegui")
    from test_nuclear_app import _page, _texts

    monkeypatch.setattr(tsm, "_install", lambda: None)

    async def scenario(user, p, on_loop):
        await user.open("/")
        await user.should_see("Beam time needed", retries=600)
        await _open_panel(user, "58Ni")
        await user.should_see("Before running", retries=100)
        assert "KSHELL was not found" in _texts(user) and not _button(user, "Run KSHELL").enabled

    _page("coulex_ni58", scenario)


def test_a_cached_result_shows_at_once_and_applies_without_replacing_ensdf(tmp_path, monkeypatch):
    pytest.importorskip("nicegui")
    from test_nuclear_app import _page, _texts, asyncio_sleep

    root = tmp_path / "cache"
    monkeypatch.setenv("PHYSIM_SHELLMODEL", str(root))
    monkeypatch.setattr(tsm, "_install", lambda: fake_install(tmp_path))
    put_in_cache(root, NI62_REQUEST, ni62_result())

    def beam(p):
        return p.experiment.levels["beam"]

    async def scenario(user, p, on_loop):
        await user.open("/")
        await user.should_see("Beam time needed", retries=600)
        energy, e2 = beam(p).levels[1].energy, beam(p).matrix_element(0, 1, "E2")
        assert (energy.source, e2.source) == ("ensdf", "derived")
        await _open_panel(user, "62Ni")
        await user.should_see("from the cache", retries=100)  # no run needed
        texts = _texts(user)
        assert "States, beside the known levels" in texts and "Applying would:" in texts
        assert "Replace existing values" in texts and "kept " in texts
        # The target, ¹⁹⁴Pt, is outside every space physim knows: said so, and import is offered.
        await _open_panel(user, "194Pt")
        await user.should_see("No interaction physim knows covers 194Pt", retries=100)
        # Apply, not replacing: the ENSDF energy and the derived matrix element stay, gaps are filled.
        assert not _checkbox(user, "Replace existing values").value
        _click(user, "Apply to the level scheme")
        await asyncio_sleep(0.3)
        s = beam(p)
        assert s.levels[1].energy.value == pytest.approx(1172.91) and s.levels[1].energy.source == "ensdf"
        assert s.matrix_element(0, 1, "E2").value == pytest.approx(e2.value)
        assert s.matrix_element(0, 1, "E2").source == "derived"
        q22 = s.matrix_element(3, 3, "E2")  # Q of the 2⁺₂, which the scheme did not have
        assert q22 is not None and q22.source == "shell model"
        assert len(s.levels) == 6  # no new level unless asked
        # With "Replace existing values" ticked, the shell model's values go in.
        await _open_panel(user, "62Ni")  # the level scheme was drawn again
        await user.should_see("Applying would:", retries=100)
        _checkbox(user, "Replace existing values").value = True
        await asyncio_sleep(0.1)
        _click(user, "Apply to the level scheme")
        await asyncio_sleep(0.3)
        s = beam(p)
        assert s.levels[1].energy.value == pytest.approx(1100.0) and s.levels[1].energy.source == "shell model"
        assert s.matrix_element(0, 1, "E2").value == pytest.approx(-25.0)
        assert s.matrix_element(0, 1, "E2").source == "shell model"

    _page("coulex_ni62_pt194", scenario)


def test_a_run_shows_its_progress_can_be_cancelled_and_says_why_it_failed(tmp_path, monkeypatch):
    pytest.importorskip("nicegui")
    from test_nuclear_app import _page, _texts, asyncio_sleep

    monkeypatch.setenv("PHYSIM_SHELLMODEL", str(tmp_path / "cache"))
    monkeypatch.setattr(tsm, "_install", lambda: fake_install(tmp_path))
    job = {"id": "0123456789abcdef", "state": "running", "runs": 0}

    def submit(request):
        job["runs"] += 1
        job["state"] = "running" if job["runs"] == 1 else "failed"
        return job["id"]

    def status(job_id):
        assert job_id == job["id"]
        st = {"job": job_id, "state": job["state"], "started": "2026-11-18T10:00:00Z", "steps": 5,
              "step_index": 2, "step": "kshell 2+ (3 states)"}
        if job["state"] == "failed":
            st.update(message="kshell.exe stopped: possible JJ collapse", log_tail="... the last lines ...",
                      finished="2026-11-18T10:00:09Z")
        return st

    def cancel(job_id):
        job["state"] = "cancelled"
        return status(job_id)

    monkeypatch.setattr(jobs, "submit", submit)
    monkeypatch.setattr(jobs, "status", status)
    monkeypatch.setattr(jobs, "cancel", cancel)

    async def scenario(user, p, on_loop):
        await user.open("/")
        await user.should_see("Beam time needed", retries=600)
        await _open_panel(user, "58Ni")
        await user.should_see("Before running", retries=100)
        _click(user, "Run KSHELL")
        await user.should_see("Step 2 of 5: kshell 2+ (3 states)", retries=100)
        assert not _button(user, "Run KSHELL").enabled  # one run at a time
        _click(user, "Cancel")
        await user.should_see("Cancelled.", retries=600)
        assert _button(user, "Run KSHELL").enabled
        _click(user, "Run KSHELL")
        await user.should_see("kshell.exe stopped: possible JJ collapse", retries=600)
        assert "The end of the log" in " ".join(str(e.props.get("label")) for e in user.find(
            kind=__import__("nicegui").ui.expansion).elements)
        await asyncio_sleep(0.1)

    _page("coulex_ni58", scenario)
