"""The planner web app (backlog item 41): every figure renders, and the app serves every example page."""

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


def test_the_app_serves_every_example():
    pytest.importorskip("nicegui")
    port = _free_port()
    # NiceGUI switches to its own test mode when it sees pytest's variables; the app here is a normal run.
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTEST")}
    env["PYTHONIOENCODING"] = "utf-8"
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
        for name in Planner.examples():
            # Building the page computes every tab on the server; an error there gives a 500.
            for mode, marker in (("expert", "Rates and beam time"), ("guided", "What do you want to measure?")):
                with urllib.request.urlopen(f"http://127.0.0.1:{port}/?example={name}&mode={mode}",
                                            timeout=120) as r:
                    assert r.status == 200
                    body = r.read().decode()
                assert "physim" in body and marker in body, (name, mode)
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
