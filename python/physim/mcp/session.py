"""The MCP server's state: the experiments folder and the one experiment open in it.

As on the app's page, one experiment is open at a time (a :class:`~physim.nuclear.planner.Planner` on an
experiment folder). The tools that read or change the setup act on it; opening or creating another switches it.
Every valid edit is saved to the experiment folder by the planner itself, so nothing is lost if the client stops.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from ..nuclear import runs as _runs
from ..nuclear.planner import Planner


class ToolProblem(ValueError):
    """A request the tools refuse, with a message meant for the model (what was wrong, and what is allowed)."""


class Session:
    """The experiments folder (``root``; by default ``$PHYSIM_EXPERIMENTS`` or ``~/.physim/experiments``) and the
    experiment open in it."""

    def __init__(self, root: Union[str, Path, None] = None):
        self._root = Path(root).expanduser() if root is not None else None
        self.planner: Optional[Planner] = None

    @property
    def root(self) -> Path:
        return self._root if self._root is not None else _runs.home()

    def require(self) -> Planner:
        """The open experiment's planner; refuses when none is open."""
        if self.planner is None:
            raise ToolProblem("No experiment is open. Open one with open_experiment (list_experiments shows them) "
                              "or make one with new_experiment (list_examples shows the starting setups).")
        return self.planner

    def switch(self, planner: Planner) -> Planner:
        """Make ``planner`` the open experiment; refuses while the open one is taking a run."""
        if self.planner is not None and self.planner.running:
            raise ToolProblem("A run is being taken in the open experiment; stop it before opening another.")
        self.planner = planner
        return planner

    def find(self, name_or_path: str) -> Path:
        """The folder of an experiment given by its path, its folder name or its name (case does not matter)."""
        p = Path(name_or_path).expanduser()
        if (p / "experiment.toml").is_file():
            return p
        found = _runs.list_experiments(self.root)
        key = name_or_path.strip().lower()
        for e in found:
            if key in (e["name"].lower(), Path(e["path"]).name.lower()):
                return Path(e["path"])
        names = ", ".join(repr(e["name"]) for e in found) or "none"
        raise ToolProblem(f"No experiment {name_or_path!r} in {self.root}. The experiments there: {names}.")
