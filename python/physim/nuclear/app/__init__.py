"""The experiment planner web app: an experiment's setup as a list of decisions, its runs and their data, in the
stages of the experiment (Setup · Plan · Run · Data · Analysis · Report · Physics).

Start it with ``physim app`` (or ``python -m physim app``); it opens in the browser at http://localhost:8080. It runs
on your computer only: no accounts, nothing uploaded. Needs the ``app`` extra: ``pip install physim-engine[app]``.

The app is a thin layer over :class:`physim.nuclear.planner.Planner` and :mod:`physim.nuclear.workbench`: every
number it shows can be had from Python too. The figure functions (``figure_geometry``, ``figure_kinematics``, ...)
return Plotly figures and work in a notebook as well.

Modules: :mod:`.figures` (figures, style), :mod:`.page` (the page for one client), :mod:`.panel` (the setup panel),
:mod:`.strip` (the status strip), one module per stage tab (:mod:`.tab_setup` ... :mod:`.tab_physics`), and
:mod:`.server` (starting the server).
"""

from .figures import *  # noqa: F401,F403 - the figures, style and field lists, as before the split
from .figures import (FIGURES, THEMES, _crystals_of, _shown, _value, figure_energy_loss,  # noqa: F401
                      figure_excitation, figure_geometry, figure_kinematics, figure_spectra, figure_strips,
                      figure_sweep, figure_trajectories, report_zip, themed)
from .server import MARKER, log_path, main, pick_port, planner_at, port_free  # noqa: F401


def build_page(*args, **kwargs):
    """Build the planner page for the current client (call inside a NiceGUI page function); see
    :func:`physim.nuclear.app.page.build_page`."""
    from .page import build_page as _build

    return _build(*args, **kwargs)


__all__ = ["FIGURES", "MARKER", "build_page", "log_path", "main", "pick_port", "planner_at", "figure_energy_loss",
           "figure_geometry", "figure_kinematics", "figure_excitation", "figure_spectra", "figure_strips",
           "figure_sweep", "figure_trajectories", "report_zip", "themed"]
