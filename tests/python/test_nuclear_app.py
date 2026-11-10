"""The planner web app (backlog items 41 and 64): every figure renders, and the app serves every example page in
one mode, with the seven stages as tabs."""

import os
import socket
import subprocess
import sys
import time
import urllib.request

import pytest

from physim.nuclear.planner import Planner

pytest.importorskip("plotly")


@pytest.mark.parametrize("name", Planner.examples())
def test_every_figure_renders(name):
    from physim.nuclear import app

    p = Planner.example(name)
    for tab, fig in app.FIGURES.items():
        f = fig(p, 20_000) if tab == "spectra" else fig(p)
        # The excitation figure has data only for a Coulomb-excitation setup.
        assert len(f.data) > 0 or (tab == "gamma" and p.experiment.excitation is None), tab
        f.to_json()  # serialisable, as the browser needs it
    first = p.experiment.detectors[0].name
    assert app.figure_strips(p, first).data[0].z.sum() > 0
    s = p.sweep("beam energy", ["4 MeV", "5 MeV"], detector=first)
    assert len(app.figure_sweep(s).data[0].x) == 2


def test_figures_take_the_page_theme():
    from physim.nuclear import app

    p = Planner.example("alpha_on_gold")
    for theme, colours in app.THEMES.items():
        for fig in (app.figure_kinematics(p), app.figure_geometry(p)):
            f = app.themed(fig, theme)
            assert f.layout.paper_bgcolor == colours["paper"] and f.layout.font.color == colours["ink"]
            f.to_json()
    assert app.themed(app.figure_kinematics(p), "dark").layout.xaxis.linecolor == app.THEMES["dark"]["ink"]


def test_field_values():
    from physim.nuclear.app import _value

    assert _value("energy", " 5 MeV ") == "5 MeV"
    assert _value("strips_x", "16") == 16
    assert _value("strips_x", "sixteen") == "sixteen"  # the setup checks report it
    assert _value("tilt", "") is None


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def test_the_app_serves_every_example(tmp_path):
    pytest.importorskip("nicegui")
    from physim.nuclear import workbench

    port = _free_port()
    # NiceGUI switches to its own test mode when it sees pytest's variables; the app here is a normal run.
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTEST")}
    env["PYTHONIOENCODING"] = "utf-8"
    env["PHYSIM_EXPERIMENTS"] = str(tmp_path)
    proc = subprocess.Popen([sys.executable, "-m", "physim", "app", "--no-browser", "--port", str(port)],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    try:
        deadline = time.time() + 60
        while True:
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=30)
                break
            except OSError:
                if time.time() > deadline or proc.poll() is not None:
                    out = proc.stdout.read().decode(errors="replace") if proc.poll() is not None else ""
                    pytest.fail(f"the app did not start: {out}")
                time.sleep(0.5)
        steps = [label for _, label in workbench.STEPS]
        pages = [f"?example={name}" for name in Planner.examples()]
        pages += [f"?template={key}" for key, _, _ in workbench.TEMPLATES] + [""]
        for query in pages:
            # Only the tab on screen is built when the page is served; ``prebuild=all`` builds every tab now, so
            # an error in any of them gives a 500 here.
            query += ("&" if query else "?") + "prebuild=all"
            with urllib.request.urlopen(f"http://127.0.0.1:{port}/{query}", timeout=120) as r:
                assert r.status == 200
                body = r.read().decode()
            assert "physim" in body and "Count rates and beam time" in body, query
            # One mode: the six steps, in order, and the Physics reference (its button in the header comes
            # first in the page); no guided steps, no mode switch, no example drop-down.
            at = [body.find(f'"label":"{label}"') for label in steps]
            assert all(i >= 0 for i in at) and at == sorted(at), (query, at)
            assert '"label":"Physics"' in body, query
            for gone in ("What do you want to measure?", "Expert view", "Start from example"):
                assert gone not in body, (query, gone)
            assert "New experiment" in body and "Run conditions" in body
            # The stages are a step bar; each step has a foot that leads on, and Physics is in the header.
            assert "ps-steps" in body and "ps-step-status" in body and "ps-step-hidden" in body, query
            assert "Next: check the plan" in body and "Next: take a run" in body, query
            assert "Open an experiment" in body and "menu_book" in body, query
            # The predicted spectra are a check of the Plan; the Data step without a run says what it needs.
            assert "Predicted spectra: what each detector will record" in body and "No run yet" in body, query
            # The Plan tab answers without a button pressed, and its checks are there (backlog 65).
            assert "Beam time needed" in body and "Kinematics: where the particles go" in body, query
            # The Run tab: the button, the real part, and what is not simulated (backlog 66).
            assert "Event by event" in body and "Not simulated yet" in body, query
            # The Physics tab: every section of the register, with the panels of the old tabs (backlog 70).
            for _, title, _, _ in workbench.PHYSICS:
                if title.startswith("Coulomb excitation") and "coulex" not in query and "ni58" not in query:
                    continue
                assert title in body, (query, title)
            for panel in ("Explain: Two-body kinematics", "Explain: Energy loss and straggling",
                          "Explain: Coulomb trajectories"):
                assert panel in body, (query, panel)
    finally:
        proc.terminate()
        proc.wait(timeout=30)


def _start(args, env_extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTEST")}
    env.update(env_extra or {})
    return subprocess.Popen([sys.executable, "-m", "physim", "app", *args], stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, env=env)


def _wait_for(predicate, seconds):
    deadline = time.time() + seconds
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(0.5)
    return False


def test_ports():
    from physim.nuclear import app

    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        taken = s.getsockname()[1]
        assert not app.port_free(taken)
        other = app.pick_port(taken)
        assert other != taken and app.port_free(other)
        assert not app.planner_at(taken)  # something else holds the port
    free = _free_port()
    assert app.pick_port(free) == free


def test_desktop_launch(tmp_path):
    """What the installer's shortcut runs: a log file instead of a console, a second launch reuses the running
    planner, and the server stops once no browser tab has been open for a while."""
    pytest.importorskip("nicegui")
    from physim.nuclear import app

    port = _free_port()
    env = {"LOCALAPPDATA": str(tmp_path), "PYTHONIOENCODING": "utf-8"}
    first = _start(["--desktop", "--no-browser", "--port", str(port), "--idle-exit", "15"], env)
    try:
        assert _wait_for(lambda: app.planner_at(port), 60), "the planner did not start"
        second = _start(["--desktop", "--no-browser", "--port", str(port)], env)
        assert second.wait(timeout=60) == 0
        log = (tmp_path / "physim" / "planner.log").read_text(encoding="utf-8")
        assert "already running" in log
        # Nothing ever connected: the server stops by itself after --idle-exit.
        assert first.wait(timeout=60) == 0
        assert "no browser tab open for 15 s" in (tmp_path / "physim" / "planner.log").read_text(encoding="utf-8")
    finally:
        if first.poll() is None:
            first.terminate()
            first.wait(timeout=30)


def test_installer_bundle_helpers(tmp_path):
    """The parts of installer/windows/build.py that run without downloading anything."""
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "installer" / "windows" / "build.py"
    spec = importlib.util.spec_from_file_location("installer_build", path)
    build = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build)

    (tmp_path / "python312._pth").write_text("python312.zip\n.\n\n# Uncomment to run site.main() automatically\n"
                                             "#import site\n")
    build.enable_site_packages(tmp_path)
    lines = (tmp_path / "python312._pth").read_text().splitlines()
    assert lines[-2:] == ["Lib\\site-packages", "import site"] and "#import site" not in lines
    assert (tmp_path / "Lib" / "site-packages").is_dir()
    major, minor = build.PYTHON_VERSION.split(".")[:2]
    assert (int(major), int(minor)) >= (3, 11)
    # build.py precompiles with the host Python, so the workflow must run the same minor version.
    workflow = (path.parents[2] / ".github" / "workflows" / "installer.yml").read_text(encoding="utf-8")
    assert f'python-version: "{major}.{minor}"' in workflow

    pytest.importorskip("PIL")
    from PIL import Image

    assert build.make_icon(tmp_path / "planner.ico")
    assert (256, 256) in Image.open(tmp_path / "planner.ico").info["sizes"]


def test_the_window_is_chosen_from_the_flags():
    """The installed shortcut (--desktop) opens the planner in its own window when pywebview is there; --browser
    or --no-browser keep the browser; without pywebview the caller is told to use the browser (None)."""
    from physim.nuclear import app

    assert app.wants_window(True, False, False, False, available=True) is True
    assert app.wants_window(False, True, False, False, available=True) is True
    assert app.wants_window(False, True, True, False, available=True) is False
    assert app.wants_window(False, True, False, True, available=True) is False
    assert app.wants_window(False, False, False, False, available=True) is False
    assert app.wants_window(True, False, False, False, available=False) is None
    assert isinstance(app.window_available(), bool)


# -- the page stays responsive (the numbers on a worker thread, the tabs built from what is kept) -------------------

def _page(example: str, scenario):
    """Run ``scenario(user, planner, on_loop)`` on the planner page of an example, in NiceGUI's simulated browser.
    ``on_loop`` collects, for every rate integral computed, whether it ran on the event loop's thread (it must
    not: the page would stop answering for as long as it takes)."""
    import asyncio
    import threading

    pytest.importorskip("nicegui")
    from nicegui.testing import user_simulation

    from physim.nuclear import planner as planner_module
    from physim.nuclear.app.page import build_page

    planner = Planner.example(example)
    on_loop = []
    rates = planner_module.Rates

    def watched(*args, **kwargs):
        on_loop.append(threading.current_thread() is threading.main_thread())
        return rates(*args, **kwargs)

    async def main():
        async with user_simulation(lambda: build_page(example=example)) as user:
            await scenario(user, planner, on_loop)
            await asyncio.sleep(0.2)  # let what the last click started end with the page still there

    original = Planner.example
    Planner.example = classmethod(lambda cls, name: planner)
    planner_module.Rates = watched
    try:
        asyncio.run(main())
    finally:
        Planner.example = original
        planner_module.Rates = rates


def _kept(planner, *key) -> bool:
    return ("page",) + key in planner._cache


def test_the_page_is_served_before_its_numbers_and_builds_the_tabs_ahead():
    pytest.importorskip("nicegui")

    async def scenario(user, p, on_loop):
        from nicegui import ui

        await user.open("/")
        # Served at once: nothing was computed to build it, and a spinner stands where the numbers will be.
        assert "rates" not in p._cache
        await user.should_see("Computing the setup's numbers")
        await user.should_see("Run conditions", retries=600)  # the setup panel, once the rates are there
        await user.should_see("Drag changes the angle", retries=600)  # the Setup tab
        assert "rates" in p._cache and "warnings" in p._cache
        # The other tabs follow without being asked for, and the Plan's checks are computed though still closed.
        await user.should_see("Beam time needed", retries=600)
        await user.should_see("Not simulated yet", retries=600)
        for _ in range(600):
            if _kept(p, "figure", "light", "strips", "A20"):
                break
            await asyncio_sleep(0.05)
        assert _kept(p, "figure", "light", "kinematics") and _kept(p, "figure", "light", "strips", "A20")
        assert on_loop == [False]  # the rates: once, on the worker thread
        # A check is built when it is opened, from the figure that is waiting for it.
        assert "Sweep one parameter" not in _texts(user)
        check = next(e for e in user.find(kind=ui.expansion).elements
                     if e.props.get("label") == "Rates per strip, and the beam-time sweep")
        check.value = True
        await user.should_see("Sweep one parameter", retries=100)

    _page("alpha_on_gold", scenario)


async def asyncio_sleep(seconds: float) -> None:
    import asyncio

    await asyncio.sleep(seconds)


def _texts(user) -> str:
    """Every label's text on the simulated page."""
    from nicegui import ui

    return " | ".join(str(e.text) for e in user.find(kind=ui.label).elements if e.visible)


def _fields(user, label):
    from nicegui import ui

    return [e for e in user.find(kind=ui.input).elements if e.props.get("label") == label]


async def _open_card(user, title: str, field: str):
    """Open a card of the setup panel (its fields are built then) and return its field ``field``."""
    from nicegui import ui

    assert not _fields(user, field)
    card = next(e for e in user.find(kind=ui.expansion).elements if "header" in e.slots and any(
        getattr(c, "text", None) == title for h in e.slots["header"].children for c in h.descendants()))
    card.value = True
    for _ in range(100):
        await asyncio_sleep(0.05)
        if _fields(user, field):
            break
    return _fields(user, field)[0]


def _type(user, box, text: str) -> None:
    """Type into a field and leave it, as a user does."""
    from nicegui.testing import UserInteraction

    box.value = text
    UserInteraction(user, {box}, None).trigger("blur")


def _press(user, label: str) -> None:
    from nicegui import ui
    from nicegui.testing import UserInteraction

    button = next(e for e in user.find(kind=ui.button).elements if e.text == label)
    UserInteraction(user, {button}, None).click()


def test_changes_wait_for_apply_and_are_computed_off_the_loop():
    pytest.importorskip("nicegui")

    async def scenario(user, p, on_loop):
        await user.open("/")
        await user.should_see("Beam time needed", retries=600)
        first, setup, beam = p._cache, p.experiment, p.experiment.beam.particles_per_second
        before = p.rates()["rows"][0]["rate_per_s"]
        current = await _open_card(user, "Beam", "Current")
        new = "3 pnA" if str(current.value) != "3 pnA" else "2 pnA"
        _type(user, current, new)
        await asyncio_sleep(0.1)
        # The change waits: the applied setup and its numbers are as they were, the field is marked, and the bar
        # at the top of the panel names the change.
        assert p.experiment is setup and p._cache is first
        assert "ps-pending" in current.classes
        texts = _texts(user)
        assert "1 change not applied" in texts and f"→ {new}" in texts and "not applied" in texts
        counts = (await _open_card(user, "Run conditions", "Counts wanted"))
        _type(user, counts, "4321")
        await asyncio_sleep(0.1)
        assert "2 changes not applied" in _texts(user) and p.experiment is setup
        # The draft's numbers are worked out ahead, on the worker, on the draft's own planner.
        await asyncio_sleep(2.0)
        assert p.experiment is setup and p._cache is first
        _press(user, "Apply")
        await asyncio_sleep(0.1)
        assert p.experiment is not setup and p._cache is not first
        assert "rates" in p._cache  # adopted from the draft: nothing to wait for
        assert p.experiment.beam.particles_per_second != beam and p.experiment.run.counts_wanted == 4321
        for _ in range(600):
            if _kept(p, "figure", "light", "strips", "A20"):  # the last thing the pass after a change computes
                break
            await asyncio_sleep(0.05)
        assert _kept(p, "figure", "light", "strips", "A20")
        assert p.rates()["rows"][0]["rate_per_s"] == pytest.approx(
            before * p.experiment.beam.particles_per_second / beam, rel=1e-9)
        texts = _texts(user)
        assert "not applied" not in texts and new in texts  # the Beam card's line; no bar, no mark
        assert "ps-pending" not in current.classes
        assert on_loop and not any(on_loop)

    _page("alpha_on_gold", scenario)


def test_discard_and_applying_each_change_as_it_is_made():
    pytest.importorskip("nicegui")
    from nicegui import ui

    async def scenario(user, p, on_loop):
        await user.open("/")
        await user.should_see("Beam time needed", retries=600)
        setup = p.experiment
        await _open_card(user, "Beam", "Current")
        energy = next(e for e in _fields(user, "Energy") if e.value != "64 MeV")  # not the New experiment form's
        was = energy.value
        _type(user, energy, "6.25 MeV")
        await asyncio_sleep(0.1)
        assert "1 change not applied" in _texts(user) and p.experiment is setup
        _type(user, energy, "fast")  # not a setup: said at the top of the panel, and Apply would not take it
        await asyncio_sleep(0.1)
        assert p.experiment is setup
        _press(user, "Apply")
        await asyncio_sleep(0.1)
        assert p.experiment is setup
        _press(user, "Discard")
        await asyncio_sleep(0.2)
        assert "not applied" not in _texts(user) and p.experiment is setup
        energy = next(e for e in _fields(user, "Energy") if str(e.value) == str(was))  # the panel was built again
        # With the switch on, a change is applied as it is made.
        switch = next(e for e in user.find(kind=ui.switch).elements if e.text == "Apply each change as I make it")
        switch.value = True
        await asyncio_sleep(0.1)
        _type(user, energy, "6.25 MeV")
        await asyncio_sleep(0.2)
        assert p.experiment is not setup and p.experiment.beam.energy_mev == pytest.approx(6.25)
        assert "not applied" not in _texts(user)
        for _ in range(600):
            if "rates" in p._cache:
                break
            await asyncio_sleep(0.05)
        assert "rates" in p._cache and not any(on_loop)

    _page("alpha_on_gold", scenario)
