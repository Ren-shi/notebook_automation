"""The experiment and setup tools: list, open and make experiments; read, replace and edit the setup.

Each tool is a method of :class:`SetupTools`; its docstring is what the model reads, so it says what the tool
does, what the arguments look like, and what comes back. The methods are plain Python: the tests call them
directly, and :mod:`physim.mcp.server` registers them as MCP tools.

Every tool that changes the setup returns the experiment's :func:`summary`, so the model sees the state after
each change without another call.
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Optional, Union

from ..nuclear import runs as _runs
from ..nuclear.experiment import SHAPES, Experiment, SetupError
from ..nuclear.planner import Planner
from .session import Session, ToolProblem
from .shape import plain

#: A value in the setup: text with a unit ("30 MeV", "100 mm"), a number (rings = 16), or a flag.
Value = Union[str, int, float, bool, None]

#: The fields a detector line in the summary shows, in order (the rest are in get_setup's TOML).
_LINE = ("shape", "model", "theta", "phi", "distance")


def _line(number: int, det: dict) -> str:
    where = ", ".join(f"{k} {det[k]}" for k in _LINE if k in det)
    return f"{number}. {det.get('name', '(no name)')}: {where}"


def summary(p: Planner) -> dict:
    """The open experiment in brief: name and folder, beam, target, reaction, run, detectors (one line each),
    level schemes, the current run, problems, and warnings (most important first)."""
    d = p.draft
    levels = {role: s.get("nuclide") for role, s in d.get("levels", {}).items()}
    return plain({
        "experiment": {"name": p.name, "folder": p.experiment_folder},
        "title": d.get("title"),
        "beam": d.get("beam"),
        "target": d.get("target"),
        "reaction": d.get("reaction", {"type": "elastic"}),
        "run": d.get("run"),
        "detectors": [_line(i + 1, det) for i, det in enumerate(d.get("detectors", []))],
        "gamma_detectors": [_line(i + 1, det) for i, det in enumerate(d.get("gamma_detectors", []))],
        "level_schemes": levels,
        "current_run": p.run_status(),
        "problems": list(p.problems),
        "warnings": [{"level": w.level, "source": w.source, "text": w.text} for w in p.warnings()],
    })


def _changed(p: Planner, ok: bool) -> dict:
    out = {"ok": ok}
    if not ok:
        out["note"] = ("The edit was kept in the setup being edited, but the setup is not valid: see problems. "
                       "Results and the saved setup stay at the last valid setup until the problems are fixed.")
    out.update(summary(p))
    return out


def _which(which: Union[str, int]) -> Union[str, int]:
    """A detector given by name, or by number counting from 1 (as the summary lists them), for the planner, which
    counts from 0."""
    if isinstance(which, bool):
        raise ToolProblem("Give a detector by its name or its number in the summary (from 1).")
    if isinstance(which, int) or (isinstance(which, str) and which.strip().isdigit()):
        n = int(which)
        if n < 1:
            raise ToolProblem("Detector numbers count from 1, as in the summary.")
        return n - 1
    return which


def _section(p: Planner, section: str) -> str:
    """``section`` with a detector named instead of numbered ("detector CD") turned into its number."""
    s = section.strip()
    for prefix, kind in (("gamma detector ", "gamma_detectors"), ("detector ", "detectors")):
        if s.lower().startswith(prefix):
            rest = s[len(prefix):].strip()
            if rest.isdigit():
                return prefix + rest
            names = [det.get("name") for det in p.draft.get(kind, [])]
            if rest not in names:
                raise ToolProblem(f"No {prefix.strip()} named {rest!r}; they are {names}.")
            return f"{prefix}{names.index(rest) + 1}"
    return s


class SetupTools:
    """The experiment and setup tools, on one :class:`~physim.mcp.session.Session`."""

    def __init__(self, session: Session):
        self.session = session

    # -- experiments -----------------------------------------------------------------------------------------------

    def list_experiments(self) -> dict:
        """List the saved experiments in the experiments folder, newest first: name, folder, number of runs, and
        when each last changed. Open one with open_experiment."""
        root = self.session.root
        found = []
        for e in _runs.list_experiments(root):
            when = _dt.datetime.fromtimestamp(e["modified"]).isoformat(timespec="minutes")
            found.append({"name": e["name"], "folder": e["path"], "runs": e["runs"], "modified": when})
        return plain({"root": root, "experiments": found})

    def list_examples(self) -> dict:
        """List the example setups a new experiment can start from (new_experiment's from_example), each with
        its title and description."""
        out = []
        for name in Planner.examples():
            exp = Experiment.example(name)
            out.append({"name": name, "title": exp.title, "description": exp.description})
        return {"examples": out}

    def open_experiment(self, name_or_path: str) -> dict:
        """Open a saved experiment by its name, its folder name, or its folder path; it becomes the open
        experiment that the other tools act on, with its last run current. Returns its summary."""
        folder = self.session.find(name_or_path)
        p = self.session.switch(Planner.open_experiment(folder))
        return summary(p)

    def new_experiment(self, name: Optional[str] = None, from_example: Optional[str] = None,
                       from_setup_file: Optional[str] = None) -> dict:
        """Make a new experiment folder and open it. Start from one of the examples (from_example, a name from
        list_examples) or from a setup file on disk (from_setup_file, a path to a .toml setup); give exactly one.
        name is the experiment's name (by default the beam and target, e.g. "16O on 58Ni"); a number is added if
        the name is taken. Returns the new experiment's summary."""
        if (from_example is None) == (from_setup_file is None):
            raise ToolProblem("Give exactly one of from_example (see list_examples) or from_setup_file.")
        if from_example is not None:
            if from_example not in Planner.examples():
                raise ToolProblem(f"No example {from_example!r}; the examples are {Planner.examples()}.")
            p = Planner.example(from_example)
        else:
            path = Path(from_setup_file).expanduser()
            if not path.is_file():
                raise ToolProblem(f"No setup file at {path}.")
            p = Planner.load(path)
        self.session.switch(p)
        p.create_experiment(name, root=self.session.root)
        return summary(p)

    # -- the setup -------------------------------------------------------------------------------------------------

    def get_setup(self) -> dict:
        """The open experiment's setup as TOML text: every field of the beam, target, reaction, run, detectors,
        γ-ray detectors and level schemes. When the setup being edited has problems, this is the last valid setup
        and the problems are listed."""
        p = self.session.require()
        return {"toml": p.to_toml(), "problems": list(p.problems)}

    def replace_setup(self, toml_text: str) -> dict:
        """Replace the whole setup of the open experiment with a setup in TOML text (the format get_setup returns
        and setup files use). It is checked first: an invalid setup is refused with its problems, and the setup
        is left as it was. A valid one is saved to the experiment folder. Returns ok and the summary."""
        p = self.session.require()
        try:
            other = Planner(Experiment.from_toml(toml_text))
        except SetupError as e:
            return {"ok": False, "problems": list(e.problems),
                    "note": "The setup was not replaced; it is as it was."}
        except Exception as e:  # noqa: BLE001 - TOML syntax errors and the like: say what the parser said
            return {"ok": False, "problems": [str(e)], "note": "The setup was not replaced; it is as it was."}
        return _changed(p, p.adopt(other))

    def set_value(self, section: str, field: str, value: Value) -> dict:
        """Set one field of the setup. section is "beam", "target", "backing" (the target's backing layer),
        "reaction", "run", "title", "description", "detector N" or "gamma detector N" (N counting from 1, or the
        detector's name, e.g. "detector CD"). value is text with a unit for quantities ("6 MeV", "0.5 mg/cm2",
        "45 deg", "100 mm", "24 h"), a number for counts (rings = 16), or null to remove an optional field.
        Examples: ("beam", "energy", "32 MeV"); ("target", "material", "208Pb"); ("reaction", "type", "coulex");
        ("detector CD", "distance", "35 mm"). Returns ok (whether the setup is valid) and the summary."""
        p = self.session.require()
        return _changed(p, p.set(_section(p, section), field, value))

    def add_detector(self, fields: dict[str, Value]) -> dict:
        """Add a particle detector. fields are as in a setup file's [[detectors]]: name, shape, theta, phi,
        distance, thickness, and optionally dead_layer, resolution, threshold. shape is "rectangle" (width,
        height, strips_x, strips_y), "annular" (inner_radius, outer_radius, rings, sectors) or "circle" (radius).
        Quantities are text with a unit. Example: {"name": "A75", "shape": "circle", "theta": "75 deg",
        "distance": "100 mm", "radius": "2.5 mm", "thickness": "300 um"}. Returns ok and the summary."""
        p = self.session.require()
        shape = fields.get("shape")
        if shape is not None and shape not in SHAPES:
            raise ToolProblem(f"Unknown shape {shape!r}; the shapes are {sorted(SHAPES)}.")
        return _changed(p, p.add_detector(**fields))

    def remove_detector(self, which: Union[str, int]) -> dict:
        """Remove a particle detector, given by name or by its number in the summary (from 1). Returns ok and
        the summary."""
        p = self.session.require()
        return _changed(p, p.remove_detector(_which(which)))

    def duplicate_detector(self, which: Union[str, int], changes: Optional[dict[str, Value]] = None) -> dict:
        """Copy a particle detector (by name or number from 1), apply changes (fields as in add_detector, e.g.
        {"phi": "180 deg"}), and add the copy under a new name. Returns ok and the summary."""
        p = self.session.require()
        return _changed(p, p.duplicate_detector(_which(which), **(changes or {})))

    def add_gamma_detector(self, fields: dict[str, Value]) -> dict:
        """Add a γ-ray detector. fields are as in a setup file's [[gamma_detectors]]: name, model (e.g. "clover",
        or a coaxial crystal described by crystal_diameter and crystal_length), theta, phi, distance, resolution,
        material, and for clovers crystals, crystal_pitch, addback, shield. Quantities are text with a unit.
        Returns ok and the summary."""
        p = self.session.require()
        return _changed(p, p.add_gamma_detector(**fields))

    def remove_gamma_detector(self, which: Union[str, int]) -> dict:
        """Remove a γ-ray detector, given by name or by its number in the summary (from 1). Returns ok and the
        summary."""
        p = self.session.require()
        return _changed(p, p.remove_gamma_detector(_which(which)))

    def duplicate_gamma_detector(self, which: Union[str, int], changes: Optional[dict[str, Value]] = None) -> dict:
        """Copy a γ-ray detector (by name or number from 1), apply changes, and add the copy under a new name.
        Returns ok and the summary."""
        p = self.session.require()
        return _changed(p, p.duplicate_gamma_detector(_which(which), **(changes or {})))

    # -- level schemes ---------------------------------------------------------------------------------------------

    def _with_levels(self, p: Planner, ok: bool) -> dict:
        out = _changed(p, ok)
        out["levels"] = plain(p.levels())
        return out

    def lookup_levels(self, role: str, max_energy: str = "3 MeV", nuclide: Optional[str] = None) -> dict:
        """Read the level scheme of the beam (role "beam") or target (role "target") nucleus from the local
        ENSDF copy into the setup, up to max_energy. nuclide picks another nuclide of the target material (by
        default its most abundant). Matrix elements set by hand for the same nuclide are kept. Returns ok, the
        summary, and the level schemes (levels numbered from 0, the ground state)."""
        p = self.session.require()
        return self._with_levels(p, p.lookup_levels(role, max_energy, nuclide))

    def use_state(self, role: str, level: int, multipolarity: str = "E2") -> dict:
        """Plan Coulomb excitation of one state of a level scheme in the setup: role "beam" or "target", level its
        number in the scheme (0 is the ground state), multipolarity "E2" by default. Sets the reaction's excited
        nucleus, state energy and B(Eλ↑) from the scheme. Returns ok, the summary and the level schemes."""
        p = self.session.require()
        if role not in p.experiment.levels:
            raise ToolProblem(f"The setup has no level scheme for the {role}; read one with lookup_levels first.")
        return self._with_levels(p, p.use_state(role, level, multipolarity))

    def set_matrix_element(self, role: str, from_level: int, to_level: int, multipolarity: str, value: str) -> dict:
        """Set a matrix element of a level scheme in the setup, marked as set by hand: role "beam" or "target",
        the two levels by number (0 is the ground state), multipolarity ("E2", "M1", ...), and value as text with
        a unit (e.g. "1.28 eb"). Returns ok, the summary and the level schemes."""
        p = self.session.require()
        if role not in p.experiment.levels:
            raise ToolProblem(f"The setup has no level scheme for the {role}; read one with lookup_levels first.")
        return self._with_levels(p, p.set_matrix_element(role, from_level, to_level, multipolarity, value))

    #: The tools in the order they are listed to the model, and whether each writes files.
    TOOLS: tuple = (
        ("list_experiments", False), ("list_examples", False), ("open_experiment", False),
        ("new_experiment", True), ("get_setup", False), ("replace_setup", True), ("set_value", True),
        ("add_detector", True), ("remove_detector", True), ("duplicate_detector", True),
        ("add_gamma_detector", True), ("remove_gamma_detector", True), ("duplicate_gamma_detector", True),
        ("lookup_levels", True), ("use_state", True), ("set_matrix_element", True),
    )

