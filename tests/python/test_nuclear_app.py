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
            with urllib.request.urlopen(f"http://127.0.0.1:{port}/?example={name}", timeout=120) as r:
                assert r.status == 200
                body = r.read().decode()
            assert "physim" in body
    finally:
        proc.terminate()
        proc.wait(timeout=30)
