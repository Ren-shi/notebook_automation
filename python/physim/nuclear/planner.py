"""The experiment planner app's model: everything the web app shows, without the web app.

The planner app (backlog item 41) is a thin layer over this module. It edits a setup, keeps the last valid
:class:`~physim.nuclear.Experiment`, and computes the content of every result tab as plain data (arrays, tables,
text) that any front end, or a notebook, can draw::

    from physim.nuclear.planner import Planner

    p = Planner.example("alpha_on_gold")
    p.set("beam", "energy", "6 MeV")          # edits are checked; problems go to p.problems
    p.add_detector(name="A75", shape="circle", theta="75 deg", distance="100 mm", radius="2.5 mm",
                   thickness="300 um")
    p.warnings()                              # validity and rate warnings, most important first
    p.rates()                                 # one row per detector
    p.explain("rates")                        # the formula, assumptions, limits, register link
    p.sweep("beam energy", ["4 MeV", "5 MeV", "6 MeV"], "rate", detector="A45")

Tabs: :data:`TABS`. Each ``Planner.<tab>()`` returns a dict; :meth:`Planner.tab` calls one by name.

**Experiments and runs** (backlog item 63, :mod:`physim.nuclear.runs`): a planner can belong to an experiment folder,
where every valid edit is saved and runs are taken and kept::

    p.create_experiment(root="~/.physim/experiments")
    run = p.start_run("beam", duration="10 min")   # unweighted events for the real part, the rest scaled
    p.run_status()                                # which run is current, and whether the setup changed since
    p.gamma_spectra()                             # read from the current run: nothing simulates on its own

While a run is being taken the setup is locked (:class:`~physim.nuclear.runs.RunInProgress`). Without a current
run the spectra, γ rays, analysis, alignment and tracks simulate their own sample as before.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np

from . import data
from . import data as _data
from . import ensdf as _ensdf
from . import response as _response
from . import scene as _scene
from .detectors import Array, Response
from .experiment import Experiment, SetupError, example_names
from .kinematics import TwoBody
from .levels import LevelScheme
from .names import parse_nuclide
from .quantity import Quantity, format_number
from .rates import Rates, beam_energy_at, beam_ion, stack, stopping, tilt_deg
from .rutherford import Rutherford
from . import runs as _runs

#: The result tabs, in display order.
TABS = ("geometry", "kinematics", "rates", "energy_loss", "spectra", "trajectories", "gamma", "report")

#: Where each tab's physics is documented and validated (paths under the docs root).
REGISTER = {
    "geometry": "physics-register/detectors",
    "kinematics": "physics-register/kinematics",
    "rates": "physics-register/rates-and-events",
    "energy_loss": "physics-register/stopping",
    "spectra": "physics-register/rates-and-events",
    "trajectories": "physics-register/rutherford",
    "gamma": "theory/coulex",
    "report": "physics-register/README",
}

#: The "Explain" panel of each tab: formula, assumptions, where it stops being valid.
EXPLAIN = {
    "geometry": {
        "title": "Detector geometry",
        "formula": "Ω = ∫ (r̂ · n̂) / r² dA over each face (Gauss–Legendre quadrature).",
        "assumptions": "Flat faces; the beam runs along +z through the target at the origin; point source unless a "
                       "beam spot is given.",
        "limits": "No gaps between strips, no inter-strip effects. A detector in the beam path is flagged.",
    },
    "kinematics": {
        "title": "Two-body kinematics",
        "formula": "Relativistic: E and θ in the lab from the CM angle by a Lorentz boost with β_cm = p₁/(E₁ + m₂).",
        "assumptions": "Nuclear masses from AME2020; the target is at rest.",
        "limits": "Exact for two-body reactions. Where the beam is heavier than the target there are two energies "
                  "at each angle, up to a maximum angle.",
    },
    "rates": {
        "title": "Count rates and beam time",
        "formula": "counts/s = (beam particles/s) × (nuclei/cm²) × ∫ dσ/dΩ dΩ; Rutherford "
                   "dσ/dΩ = (Z₁Z₂e²/4E)² / sin⁴(θ/2). Beam time for N counts = N / rate; relative error 1/√N.",
        "assumptions": "Elastic Rutherford scattering, averaged over the target depth; ejectile and recoil, target "
                       "and backing.",
        "limits": "Rutherford fails above the Coulomb barrier at large angles, at very low energy (screening) and "
                  "for identical nuclei (Mott). The warnings say when.",
    },
    "energy_loss": {
        "title": "Energy loss and straggling",
        "formula": "E_out = R⁻¹(R(E_in) − path), with R the CSDA range from tabulated (ICRU 49/90) and modelled "
                   "stopping powers; straggling from Bohr's formula carried through the layer.",
        "assumptions": "Paths grow as 1/cos of the angle to the layer normal; fully stripped ions at high velocity, "
                       "effective charge below.",
        "limits": "Slow heavy ions (below ~1 MeV/u) are uncertain to 10–50%; straggling below 10 MeV/u to ~2×. "
                  "Ions heavier than Ar in layers with hydrogen or helium lose up to 20% too little energy below "
                  "10 MeV/u; the warnings say when.",
    },
    "spectra": {
        "title": "Simulated spectra",
        "formula": "Monte Carlo: reaction depth, beam spread and spot, straggling, Rutherford angle, kinematics, "
                   "exit energy loss and detector response, event by event; each event weighted so spectra are "
                   "counts in the planned beam time.",
        "assumptions": "Straight tracks; Gaussian straggling and resolution.",
        "limits": "No multiple scattering, beam divergence, background or pile-up.",
    },
    "trajectories": {
        "title": "Coulomb trajectories",
        "formula": "Classical orbits integrated by physim's engine in the CM frame; deflection 2 arctan(d₀/2b), "
                   "closest approach d₀/2 (1 + 1/sin(θ/2)).",
        "assumptions": "Point charges, classical mechanics (valid when the Sommerfeld parameter η ≫ 1).",
        "limits": "Inside the interaction radius nuclear forces act and the orbit is no longer Coulomb.",
    },
    "gamma": {
        "title": "Coulomb excitation and γ rays",
        "formula": "First-order semiclassical Coulomb excitation: P(θ) from the multipole field integrated along the "
                   "Rutherford orbit (Alder and Winther); dσ/dΩ = P dσ_R/dΩ. γ energy E₀ √(1 − β²)/(1 − β cos α) "
                   "for a nucleus moving at β, emitting at α to its velocity.",
        "assumptions": "One state reached from a 0⁺ ground state; B(Eλ) as given in the setup; decay in flight after "
                       "leaving the target, with the particle–γ angular correlation of first-order excitation "
                       "(or isotropic, by the switch).",
        "limits": "Closer than Cline's safe distance nuclear forces interfere; lifetimes and deorientation are not "
                  "modelled. The rates and the table of all levels are first order; 'Solve with all orders' "
                  "includes multi-step excitation and reorientation for the level scheme. With a level scheme "
                  "the simulated events, and the excitation probability of the rates, take the populations of "
                  "all orders too (unless the setup says order = \"first\").",
    },
    "report": {
        "title": "The report",
        "formula": "One document: the setup, the plan, every run with its spectra, the gates, the Doppler "
                   "correction, the analysis, and every number with its explanation.",
        "assumptions": "As for each tab.",
        "limits": "Only what the tabs compute and what the runs recorded; check the warnings first.",
    },
}


@dataclass
class Warning:
    """A message for the warnings banner. ``level`` is "error" (the setup cannot be used), "warning" or "note"."""

    level: str
    source: str
    text: str


def _q(x) -> Quantity:
    return Quantity.parse(x) if isinstance(x, str) else x


#: Reactions a lifetime measurement goes with.
LIFETIME_REACTIONS = ("coulex", "fusion-evaporation")

#: The fields of a lifetime section that belong to a plunger (RDDS) only.
PLUNGER_FIELDS = ("distances", "offset", "deorientation", "stopper", "degrader2")


def first_state(nuclide: str, multipolarity: str = "E2", max_energy_kev: float = 3000.0) -> tuple:
    """The lowest state of ``nuclide`` reached from the ground state by ``multipolarity`` with a known B(Eλ↑), from
    the local ENSDF copy: (scheme, level index, B(Eλ↑) in e²fm^2λ). Raises when there is no copy, no data for the
    nuclide or no such state."""
    scheme = LevelScheme.from_ensdf(nuclide, max_energy_kev=max_energy_kev)
    lam = int(multipolarity[1])
    for n in range(1, len(scheme.levels)):
        if scheme.b(0, n, multipolarity):
            return scheme, n, scheme.b(0, n, multipolarity)
    raise ValueError(f"no state with a known B({multipolarity}↑) from the ground state of {nuclide}")


class Planner:
    """An editable setup and the results computed from it.

    Edits go to a working copy of the setup file's structure. When the copy is valid it becomes :attr:`experiment`
    and the results are recomputed (lazily); otherwise :attr:`problems` lists what is wrong and the results keep
    showing the last valid setup.
    """

    def __init__(self, experiment: Experiment):
        experiment.validate()
        self._draft: dict = experiment.to_dict()
        self.experiment: Experiment = experiment
        #: Problems with the current draft (empty when it is valid).
        self.problems: list = []
        self._cache: dict = {}
        #: The rates before the last edit: detectors the edit left alone keep their rows (see :class:`Rates`).
        self._earlier_rates: Optional[Rates] = None
        #: The Coulomb-excitation settings while the reaction is switched to elastic, so switching back restores them.
        self._stashed_excitation: dict = {}
        #: The fusion-evaporation settings (the channel table) while another reaction is chosen.
        self._stashed_fusion: dict = {}
        #: The lifetime section while the reaction is one without lifetimes (elastic).
        self._stashed_lifetime: dict = {}
        #: A plunger's distances and foils while the lifetime method is DSAM, and a differential plunger's second
        #: foil while the stopper is of another kind, so switching back restores them.
        self._stashed_plunger: dict = {}
        self._stashed_degrader2: dict = {}
        #: The experiment folder this setup belongs to (None until one is created or opened).
        self.folder: Optional[_runs.ExperimentFolder] = None
        #: The current run, which the spectra, analysis and the other views read (None: they simulate their own).
        self.run: Optional[_runs.RunData] = None
        #: A run being taken on a thread.
        self._background: Optional[_runs.Background] = None
        self._taking: Optional[_runs.RunTaker] = None
        self._gates = None

    # -- loading and saving -----------------------------------------------------------------------------------------

    @classmethod
    def example(cls, name: str) -> Planner:
        return cls(Experiment.example(name))

    @staticmethod
    def examples() -> list:
        return example_names()

    @classmethod
    def load(cls, path: Union[str, Path]) -> Planner:
        return cls(Experiment.load(path))

    def save(self, path: Union[str, Path]) -> None:
        """Write the current (valid) setup to a setup file."""
        if self.problems:
            raise SetupError(self.problems)
        self.experiment.save(path)

    def to_toml(self) -> str:
        return self.experiment.to_toml()

    # -- editing ----------------------------------------------------------------------------------------------------

    @property
    def draft(self) -> dict:
        """A copy of the setup being edited, in the structure of a setup file."""
        return copy.deepcopy(self._draft)

    def _apply(self, draft: dict) -> bool:
        if self.running:
            raise _runs.RunInProgress("the setup cannot change while a run is being taken; stop the run first")
        self._draft = draft
        try:
            exp = Experiment.from_dict(copy.deepcopy(draft))
        except SetupError as e:
            self.problems = list(e.problems)
            return False
        self.problems = []
        self.experiment = exp
        self._earlier_rates = self._cache.get("rates", self._earlier_rates)
        # A new place for the new setup's results, not the old one emptied: a :meth:`snapshot` still computing
        # the old setup's results on another thread keeps writing where they belong.
        self._cache = {}
        if self.folder is not None:
            self.folder.save_setup(exp)
        return True

    def adopt(self, other: Planner) -> bool:
        """Take the setup ``other`` holds, and the results it has kept for it. This is how a draft edited on a
        planner of its own is applied (the planner adopts the draft) or thrown away (the draft adopts the
        planner). An invalid setup is not taken: returns whether it was. In an experiment folder the setup is
        saved, as after any edit."""
        if other.problems:
            return False
        if not self._apply(other.draft):
            return False
        # The same setup, so the same results: what the other has computed (or is still computing on a worker
        # thread) is this planner's too.
        self._cache = other._cache
        return True

    def snapshot(self) -> Planner:
        """This planner as it is now, for a computation on another thread. It has the same setup, run and kept
        results (what it computes, the planner has too) and does not follow later edits: a result still being
        computed when the setup changes is kept with the setup it belongs to, never with the new one."""
        return copy.copy(self)

    # -- experiments and runs ---------------------------------------------------------------------------------------

    @property
    def experiment_folder(self) -> Optional[Path]:
        """The folder of the experiment this setup belongs to, or None."""
        return self.folder.path if self.folder is not None else None

    @property
    def name(self) -> str:
        """The experiment's name ("58Ni on 208Pb" by default)."""
        return self.folder.name if self.folder is not None else _runs.default_name(self.experiment)

    @staticmethod
    def experiments(root=None) -> list:
        """The experiments on disk, newest first (:func:`physim.nuclear.runs.list_experiments`)."""
        return _runs.list_experiments(root)

    def create_experiment(self, name: Optional[str] = None, root=None) -> Path:
        """Make an experiment folder for the current setup (under ``root``, by default ``~/.physim/experiments``)
        and work in it from now on. Returns its path."""
        if self.problems:
            raise SetupError(self.problems)
        self.folder = _runs.ExperimentFolder.create(self.experiment, name, root)
        self.run = None
        return self.folder.path

    @classmethod
    def open_experiment(cls, path) -> Planner:
        """A planner on an experiment folder's current setup, with its last run current."""
        folder = _runs.ExperimentFolder(path)
        p = cls(folder.setup())
        p.folder = folder
        from . import result_rows

        result_rows.backfill_if_new(folder.path)  # an experiment from before the results ledger
        runs = folder.runs()
        if runs:
            p.run = folder.load_run(runs[-1]["number"])
        return p

    def rename_experiment(self, name: str) -> None:
        if self.folder is None:
            raise ValueError("there is no experiment folder yet (create_experiment)")
        self.folder.rename(name)

    # -- gates ------------------------------------------------------------------------------------------------------

    @property
    def gatebook(self):
        """The experiment's named gates (:class:`physim.nuclear.dataviews.GateBook`), saved in its folder."""
        from .dataviews import GateBook

        folder = self.folder.path if self.folder is not None else None
        if self._gates is None or self._gates.folder != folder:
            old = self._gates
            self._gates = GateBook.load(folder)
            if old is not None and old.folder is None and folder is not None:  # gates made before the folder
                for g in old.gates.values():
                    self._gates.put(g)
        return self._gates

    def gates(self) -> list:
        """The experiment's gates, in the order they were made."""
        return list(self.gatebook.gates.values())

    def put_gate(self, gate) -> None:
        """Add or replace a gate (a :class:`~physim.nuclear.dataviews.Gate`, or a dict of its fields)."""
        from .dataviews import Gate

        self.gatebook.put(gate if isinstance(gate, Gate) else Gate.from_dict(gate))

    def remove_gate(self, name: str) -> None:
        self.gatebook.remove(name)

    def gate(self, name: str):
        return self.gatebook[name]

    def runs(self) -> list:
        """The summaries of the experiment's runs, in order, each with a one-line description and whether the
        setup's physics changed since (``stale``)."""
        if self.folder is None:
            return []
        out = []
        current = self.experiment.to_dict()
        for s in self.folder.runs():
            r = _runs.RunData(Path(s["_folder"]), s, _runs.taken_setup(s["_folder"], s.get("kind")).to_dict())
            s = dict(s, describe=r.describe(), stale=bool(r.changes(current)), label=s["_label"],
                     one_nucleus=r.one_nucleus)
            out.append(s)
        return out

    @property
    def running(self) -> bool:
        """Whether a run is being taken (the setup is locked meanwhile)."""
        return self._taking is not None

    def _plan_numbers(self) -> dict:
        """The Plan's predictions at this moment, kept in the next run's summary (to compare with what it
        measured)."""
        try:
            return self.plan()
        except Exception:  # noqa: BLE001 - a setup whose plan fails still runs; it just has no predictions
            return {}

    def start_run(self, kind: str = "beam", duration=None, budget=_runs.DEFAULT_BUDGET_S, seed: Optional[int] = None,
                  background: bool = False, progress=None, **options):
        """Take a run with the current setup: ``kind`` "beam", "source" or "alignment"; ``duration`` (text with a
        unit, by default the setup's beam time); ``budget``, the part simulated event by event (10 min by default,
        at most an hour). Source runs take ``source``, ``activity``, ``position``; alignment runs ``offset_mm``.

        Creates an experiment folder first if there is none. The setup is locked until the run ends. With
        ``background`` the run is taken on a thread and this returns at once (see :meth:`run_progress`,
        :meth:`stop_run`); otherwise it returns the stored run, which becomes current."""
        if self.running:
            raise _runs.RunInProgress("a run is already being taken")
        if self.problems:
            raise SetupError(self.problems)
        if self.folder is None:
            self.create_experiment()
        out = _runs.start(self.folder, kind, duration, budget, seed, self.experiment, self._plan_numbers(), **options)
        if isinstance(out, _runs.RunData):
            self._set_run(out)
            return out
        return self._take(out, background, progress)

    def _take(self, taker: _runs.RunTaker, background: bool, progress):
        self._taking = taker
        if background:
            def done(bg):
                self._taking = None
                if bg.result is not None:
                    self._set_run(bg.result)

            self._background = _runs.Background(taker, done)
            return self._background
        try:
            while not taker.done:
                taker.step()
                if progress is not None:
                    progress(taker.progress())
            run = taker.write()
        finally:
            self._taking = None
        self._set_run(run)
        return run

    def _set_run(self, run: Optional[_runs.RunData]) -> None:
        self.run = run
        for key in [k for k in list(self._cache) if isinstance(k, tuple) and k and k[0] in ("run",)]:
            self._cache.pop(key, None)

    def stop_run(self) -> Optional[_runs.RunData]:
        """Stop the run being taken: what is accumulated is kept (the run's duration becomes the beam time taken).
        Returns the stored run once it is written."""
        if self._background is None or self._taking is None:
            return self.run
        self._background.stop()
        return self._background.join()

    def extend_run(self, more, background: bool = False, progress=None):
        """Add ``more`` beam time to the current run, with its own setup and seed stream; the real part grows up to
        the run's budget."""
        if self.run is None:
            raise ValueError("there is no current run to extend")
        if self.running:
            raise _runs.RunInProgress("a run is already being taken")
        return self._take(_runs.extend(self.run, more), background, progress)

    def retake_run(self, number: Optional[int] = None, background: bool = False, progress=None):
        """Take run ``number`` (the current run by default) again as a new run: its own setup, seed, duration and
        real part, with today's events. A run taken before the events held both nuclei
        (:attr:`~physim.nuclear.runs.RunData.one_nucleus`) comes back with the γ rays of both; the old run
        stays. The new run becomes current."""
        if self.folder is None:
            raise ValueError("there is no experiment folder")
        if self.running:
            raise _runs.RunInProgress("a run is already being taken")
        run = self.run if number is None else self.folder.load_run(number)
        if run is None:
            raise ValueError("there is no run to take again")
        try:  # the plan of that run's own setup, for predicted against measured
            plan = Planner(run.experiment).plan()
        except Exception:  # noqa: BLE001 - as for any run, a plan that fails leaves it without predictions
            plan = {}
        return self._take(_runs.retake(self.folder, run, plan), background, progress)

    def load_run(self, number: int) -> _runs.RunData:
        """Make run ``number`` current (its events are loaded when a view asks for them; one run at a time)."""
        if self.folder is None:
            raise ValueError("there is no experiment folder")
        self._set_run(self.folder.load_run(number))
        return self.run

    def estimate_run(self, duration=None, budget=_runs.DEFAULT_BUDGET_S) -> dict:
        """What a beam run of ``duration`` with a real part of ``budget`` will take: its real and scaled parts,
        about how many particles, the CPU time and the disk space (:func:`physim.nuclear.runs.estimate`)."""
        return _runs.estimate(self.experiment, duration, budget, self._rates())

    def run_sample(self, n: int = 12) -> list:
        """The last particles the run being taken counted (empty when none is being taken)."""
        return self._taking.sample(n) if self._taking is not None else []

    def run_progress(self) -> Optional[dict]:
        """The live counters of the run being taken (None when none is)."""
        return self._taking.progress() if self._taking is not None else None

    def run_status(self, run: Optional[_runs.RunData] = None) -> Optional[dict]:
        """About the current run (or ``run``): number, label, one-line description, and whether the setup's
        physics changed since it was taken (``stale``, with the ``changes``). None without a run."""
        run = run or self.run
        if run is None:
            return None
        changes = run.changes(self.experiment.to_dict())
        return {"number": run.number, "label": run.label, "kind": run.kind, "describe": run.describe(),
                "duration_s": run.duration_s, "real_s": run.real_s, "scaled": run.scale > 1.000001,
                "stale": bool(changes), "changes": changes}

    def _data_experiment(self) -> Experiment:
        """The setup the data views read: the current beam run's own, or the current setup without one."""
        return self.run.experiment if self._beam_run() else self.experiment

    def _beam_run(self) -> Optional[_runs.RunData]:
        return self.run if self.run is not None and self.run.kind != "source" else None

    def place(self, key: str, **fields) -> bool:
        """Set several fields of one detector at once; ``key`` is "detector:N" or "gamma:N", counting from 0 as
        the scene does (:mod:`physim.nuclear.scene`). ``None`` removes a field."""
        d = self.draft
        kind, i = key.split(":")
        entry = d["detectors" if kind == "detector" else "gamma_detectors"][int(i)]
        for k, v in fields.items():
            if v is None:
                entry.pop(k, None)
            else:
                entry[k] = v
        return self._apply(d)

    def move(self, key: str, position, mode: str = "angle") -> bool:
        """Move a detector's front face to ``position`` (x, y, z in mm), as a drag in the scene does: in ``mode``
        "angle" it keeps its distance and takes the new direction, in "distance" it keeps its direction. The
        result is the same as typing the new angles or distance."""
        return self.place(key, **_scene.move_fields(self.experiment, key, position, mode))

    def set(self, section: str, field: str, value: Any) -> bool:
        """Set ``field`` of ``section`` ("beam", "target", "backing", "run", "reaction", "detector N" or
        "gamma detector N", counting from 1); ``None`` removes an optional field. Returns whether the setup is valid
        afterwards.

        Setting the reaction ``type`` to ``"coulex"`` fills in E2 excitation of the target, with the first state
        that has a known B(E2↑) from the local ENSDF copy when there is one; changing the excited nucleus or the
        multipolarity looks the state up again. Switching back to ``"elastic"`` keeps the excitation settings for
        later."""
        d = self.draft
        if section == "reaction" and field == "type":
            return self._set_reaction_type(d, value)
        if section.startswith("reaction channel"):
            sec = d.get("reaction", {}).get("channels", [])[self._channel_index(section)]
            if value is None:
                sec.pop(field, None)
            else:
                sec[field] = value
            self._drop_orphan_levels(d)
            return self._apply(d)
        if section == "reaction" and field in ("excite", "multipolarity") and value is not None:
            d.setdefault("reaction", {})[field] = value
            if d["reaction"].get("type") == "coulex":
                d["reaction"].pop("energy", None)
                d["reaction"].pop("b_up", None)
                self._fill_state(d)
            return self._apply(d)
        if section == "backing":
            sec = d["target"].setdefault("backing", {})
        elif section.startswith("gamma detector"):
            sec = d["gamma_detectors"][self._index(section[len("gamma "):], "gamma_detectors")]
        elif section.startswith("detector"):
            sec = d["detectors"][self._index(section)]
        elif section in ("title", "description"):
            d[section] = value
            return self._apply(d)
        else:
            sec = d.setdefault(section, {})
        if value is None:
            sec.pop(field, None)
        else:
            sec[field] = value
        if section == "backing" and not sec:
            d["target"].pop("backing")
        if (section, field) in (("beam", "nuclide"), ("target", "material")):
            self._new_nucleus(d, section)
            self._drop_orphan_levels(d)
        return self._apply(d)

    def _new_nucleus(self, d: dict, role: str) -> None:
        """After the beam or target changed: its old level scheme goes, and when it is the excited nucleus its state
        is looked up again."""
        nuclide = self.level_nuclides(d)[role]
        if nuclide == self.level_nuclides()[role]:
            return
        scheme = d.get("levels", {}).get(role)
        if scheme is not None and scheme.get("nuclide") != nuclide:
            d["levels"].pop(role)
            if not d["levels"]:
                d.pop("levels")
        reaction = d.get("reaction", {})
        excited = "target" if reaction.get("excite", "target") == "target" else "beam"
        if reaction.get("type") == "coulex" and excited == role:
            reaction.pop("energy", None)
            reaction.pop("b_up", None)
            self._fill_state(d)

    def _set_reaction_type(self, d: dict, value: str) -> bool:
        reaction = d.setdefault("reaction", {"type": "elastic"})
        old = reaction.get("type", "elastic")
        kept = {k: v for k, v in reaction.items() if k != "type"}
        if kept and old == "fusion-evaporation":
            # The residues' level schemes go with the channel table, and come back with it.
            residue_levels = {role: s for role, s in d.get("levels", {}).items() if role not in ("beam", "target")}
            from .experiment import ANCILLARIES

            self._stashed_fusion = dict(kept, _levels=residue_levels, _wall=d.pop("neutron_wall", None),
                                        _ancillaries={a[0]: d.pop(a[0]) for a in ANCILLARIES if a[0] in d})
        elif kept and old != value:
            self._stashed_excitation = kept
        if value == old:
            return self._apply(d)
        # A lifetime measurement goes with Coulex and fusion-evaporation; it waits while another reaction is chosen.
        if value not in LIFETIME_REACTIONS and "lifetime" in d:
            self._stashed_lifetime = d.pop("lifetime")
        elif value in LIFETIME_REACTIONS and "lifetime" not in d and self._stashed_lifetime:
            d["lifetime"] = copy.deepcopy(self._stashed_lifetime)
        if value == "coulex":
            reaction = d["reaction"] = {"type": "coulex"}
            reaction.update(self._stashed_excitation or {"excite": "target", "multipolarity": "E2"})
            if not reaction.get("energy") or not reaction.get("b_up"):
                self._fill_state(d)
        elif value == "fusion-evaporation":
            stash = copy.deepcopy(self._stashed_fusion) or self._starter_channels(d)
            levels = stash.pop("_levels", {})
            wall = stash.pop("_wall", None)
            if wall:
                d["neutron_wall"] = wall
            d.update(stash.pop("_ancillaries", {}))
            d["reaction"] = {"type": value, **stash}
            if levels:
                d.setdefault("levels", {}).update(levels)
        else:
            d["reaction"] = {"type": value}
        self._drop_orphan_levels(d)
        return self._apply(d)

    @staticmethod
    def _drop_orphan_levels(d: dict) -> None:
        """Level schemes of residues the channel table no longer makes leave the draft (a scheme must belong to a
        nucleus of the setup). Left alone when the residues cannot be worked out (the setup's problems say why)."""
        levels = d.get("levels")
        if not levels:
            return
        keep = {"beam", "target"}
        reaction = d.get("reaction", {})
        if reaction.get("type") == "fusion-evaporation":
            from .residues import parse_particles, residue_of

            try:
                b = parse_nuclide(d["beam"]["nuclide"])
                z, a, _ = max(data.material(d["target"]["material"]).atoms, key=lambda t: t[2])
                cn = (b[0] + z, b[1] + a)
                for ch in reaction.get("channels", []):
                    keep.add(data.nuclide(residue_of(cn, parse_particles(ch["particles"]))).name)
            except Exception:  # noqa: BLE001 - an unfinished setup: keep every scheme for now
                return
        for role in [r for r in levels if r not in keep]:
            levels.pop(role)
        if not levels:
            d.pop("levels")

    def _starter_channels(self, d: dict) -> dict:
        """A first channel table for a setup switched to fusion-evaporation: the xn channel the excitation energy
        favours (:func:`physim.nuclear.residues.likely_xn`), with a placeholder cross section the planner warns
        about. Without a valid beam and target, an empty table (the setup says what is missing)."""
        from .residues import likely_xn

        try:
            exp = self.experiment
            from .experiment import main_nuclide

            k = likely_xn(beam_ion(exp), main_nuclide(exp.target), exp.beam.energy_mev)
        except Exception:  # noqa: BLE001 - no valid setup to estimate from
            k = 0
        if k < 1:
            return {"channels": []}
        return {"channels": [{"particles": f"{k}n", "cross_section": "100 mb", "source": "placeholder"}]}

    def _channel_index(self, which) -> int:
        """A channel of the fusion-evaporation table: "reaction channel N" (from 1), its particles ("4n") or
        residue ("158Dy"), or an index from 0."""
        chans = self._draft.get("reaction", {}).get("channels", [])
        if isinstance(which, str):
            key = which[len("reaction channel"):].strip() if which.startswith("reaction channel") else which
            if key.isdigit():
                i = int(key) - 1
            else:
                found = [k for k, c in enumerate(chans) if key in (c.get("particles"), c.get("residue"))]
                if not found:
                    raise KeyError(f"no channel {key!r}")
                i = found[0]
        else:
            i = int(which)
        if not 0 <= i < len(chans):
            raise IndexError(f"there are {len(chans)} channels")
        return i

    def add_channel(self, **fields) -> bool:
        """Append an exit channel to the fusion-evaporation table (fields as in ``[[reaction.channels]]``)."""
        d = self.draft
        if d.get("reaction", {}).get("type") != "fusion-evaporation":
            raise ValueError("channels belong to a fusion-evaporation reaction: set the reaction type first")
        d["reaction"].setdefault("channels", []).append(dict(fields))
        return self._apply(d)

    def remove_channel(self, which) -> bool:
        """Remove an exit channel, by index (from 0), particles or residue."""
        d = self.draft
        d["reaction"]["channels"].pop(self._channel_index(which))
        self._drop_orphan_levels(d)
        return self._apply(d)

    def physim_channels(self, events: Optional[int] = None, seed: int = 1):
        """The exit channels from physim's own fusion and evaporation
        (:func:`physim.nuclear.exit_channels.physim_channels`) for the current setup. Kept until the setup
        changes."""
        from .exit_channels import EVENTS, physim_channels

        if self.experiment.fusion is None:
            raise ValueError("the setup is not a fusion-evaporation reaction")
        key = ("physim_channels", events or EVENTS, seed)
        if key not in self._cache:
            self._cache[key] = physim_channels(self.experiment, events or EVENTS, seed)
        return self._cache[key]

    def fill_channels(self, events: Optional[int] = None, seed: int = 1) -> bool:
        """Replace the channel table with the one physim's own fusion and evaporation give (each channel's
        ``source`` says so). The setup keeps its shape: the channels can still be edited, or replaced by numbers
        from PACE4 or a measurement. Level schemes of residues no longer made leave the setup."""
        if self.experiment.fusion is None:
            raise ValueError("channels belong to a fusion-evaporation reaction: set the reaction type first")
        made = self.physim_channels(events, seed)
        d = self.draft
        d["reaction"]["channels"] = made.table.to_setup()
        self._drop_orphan_levels(d)
        return self._apply(d)

    def _fill_state(self, d: dict) -> None:
        """The Coulomb-excitation reaction of ``d`` gets the excited nucleus's first state with a known B(Eλ↑) from
        the local ENSDF copy, and its level scheme; without a copy or data it is left for the user to enter."""
        reaction = d["reaction"]
        role = "target" if reaction.get("excite", "target") == "target" else "beam"
        multipolarity = reaction.get("multipolarity", "E2")
        nuclide = self.level_nuclides(d)[role]
        if nuclide is None:
            return
        try:
            scheme, level, b_up = first_state(nuclide, multipolarity)
        except Exception:  # noqa: BLE001 - no ENSDF copy, or no such state: the user enters it
            return
        reaction["energy"] = f"{scheme.levels[level].energy.value:.10g} keV"
        reaction["b_up"] = f"{b_up:.6g} e2fm{2 * int(multipolarity[1])}"
        d.setdefault("levels", {})[role] = scheme.to_dict()

    def _index(self, which, kind: str = "detectors") -> int:
        n = len(self._draft.get(kind, []))
        if isinstance(which, str):
            if which.startswith("detector "):
                i = int(which.split()[1]) - 1
            else:
                names = [det.get("name") for det in self._draft.get(kind, [])]
                if which not in names:
                    raise KeyError(f"no detector named {which!r}")
                i = names.index(which)
        else:
            i = int(which)
        if not 0 <= i < n:
            raise IndexError(f"there are {n} {kind.replace('_', ' ')}")
        return i

    def add_gamma_detector(self, **fields) -> bool:
        """Append a γ-ray detector with the given fields (as in a setup file's ``[[gamma_detectors]]``)."""
        d = self.draft
        d.setdefault("gamma_detectors", []).append(dict(fields))
        return self._apply(d)

    def remove_gamma_detector(self, which) -> bool:
        """Remove a γ-ray detector, by index (from 0) or name."""
        d = self.draft
        d["gamma_detectors"].pop(self._index(which, "gamma_detectors"))
        if not d["gamma_detectors"]:
            d.pop("gamma_detectors")
        return self._apply(d)

    def duplicate_gamma_detector(self, which, **changes) -> bool:
        """Copy a γ-ray detector (by index or name), apply ``changes``, and append it under a new name."""
        return self._duplicate("gamma_detectors", which, changes)

    def add_detector(self, **fields) -> bool:
        """Append a detector with the given fields (as in a setup file)."""
        d = self.draft
        d.setdefault("detectors", []).append(dict(fields))
        return self._apply(d)

    def remove_detector(self, which) -> bool:
        """Remove a detector, by index (from 0) or name."""
        d = self.draft
        d["detectors"].pop(self._index(which))
        return self._apply(d)

    def duplicate_detector(self, which, **changes) -> bool:
        """Copy a detector (by index or name), apply ``changes``, and append it; a copy gets a new name."""
        return self._duplicate("detectors", which, changes)

    def _duplicate(self, kind: str, which, changes: dict) -> bool:
        d = self.draft
        det = copy.deepcopy(d[kind][self._index(which, kind)])
        names = {x.get("name") for x in d[kind]}
        base = det.get("name") or "D"
        k = 2
        while f"{base}-{k}" in names:
            k += 1
        det["name"] = f"{base}-{k}"
        det.update(changes)
        d[kind].append(det)
        return self._apply(d)

    # -- shared results -----------------------------------------------------------------------------------------

    def _rates(self) -> Rates:
        if "rates" not in self._cache:
            self._cache["rates"] = Rates(self.experiment, previous=self._earlier_rates)
            self._earlier_rates = None
        return self._cache["rates"]

    # -- the scene ----------------------------------------------------------------------------------------------

    def live(self, key: str, fields: Optional[dict] = None) -> dict:
        """The quantities that follow a detector while it is dragged, computed quickly: its angles, solid angle,
        how much of it other detectors hide, and (for a particle detector) a coarse estimate of its rate.
        ``fields`` are the setup fields of the trial position (see :func:`physim.nuclear.scene.move_fields`); the
        setup itself is not changed. "problem" says why the position is not allowed, if it is not."""
        try:
            exp = _scene.moved(self.experiment, key, fields) if fields else self.experiment
        except SetupError as e:
            return {"problem": "; ".join(e.problems)}
        found = _scene.problems(exp, only=key)
        i = int(key.split(":")[1])
        if key.startswith("gamma"):
            gd = exp.gamma_detectors[i]
            return {"name": gd.name or f"γ{i + 1}", "theta": _q(gd.theta).to("deg"),
                    "phi": _q(gd.phi).to("deg") if gd.phi is not None else 0.0,
                    "distance_mm": gd.distance_mm(), "half_angle_deg": gd.half_angle_deg(),
                    "geometric_efficiency": gd.geometric_efficiency(),
                    "peak_efficiency": gd.peak_efficiency(self.gamma_energy_mev(), exp),
                    "problem": found[0].text if found else None}
        array = Array.from_experiment(exp)
        g = array.geometries[i]
        quick = Rates(exp, depth_points=2, order=3, only=g.name)
        hidden = sum(frac for (_, behind), frac in array.shadowing().items() if behind == g.name)
        x, y, z = g.centre
        return {"name": g.name, "theta": math.degrees(math.atan2(math.hypot(x, y), z)),
                "phi": math.degrees(math.atan2(y, x)), "distance_mm": float(np.linalg.norm(g.centre)),
                "theta_range": g.theta_range(), "solid_angle_msr": g.solid_angle(), "hidden": min(hidden, 1.0),
                "rate_per_s": quick.rate(g.name, counted=False), "problem": found[0].text if found else None}

    def gamma_energy_mev(self) -> float:
        """The γ-ray energy the efficiencies are quoted at: that of the excited state, or 1332 keV without one."""
        exc = self.experiment.excitation
        return exc.energy_mev if exc is not None else _response.REFERENCE_MEV

    def efficiency(self, key: str, points: int = 120) -> dict:
        """The efficiency curve of γ-ray detector ``key`` ("gamma:N", from 0): full-energy-peak and total
        efficiency (fractions of all γ rays emitted at the target) and the resolution against energy, with what
        the response is built from (see :mod:`physim.nuclear.response`)."""
        gd = self.experiment.gamma_detectors[int(key.split(":")[1])]
        r = _response.Response(self.experiment, gd)
        e = np.geomspace(0.03, 3.0, points)
        return {"name": gd.name or f"γ{int(key.split(':')[1]) + 1}", "energy_mev": e,
                "peak": r.peak_efficiency(e), "total": r.total_efficiency(e), "fwhm_kev": 1e3 * r.fwhm(e),
                "transmission": r.transmission(e), **r.describe()}

    def source_run(self, nuclide: str = "152Eu", activity: str = "37 kBq", time: str = "1 h", seed: int = 1):
        """A simulated run with a calibration source at the target position, without beam
        (:func:`physim.nuclear.response.source_run`). Kept until the setup changes."""
        key = ("source", nuclide, activity, time, seed)
        if key not in self._cache:
            self._cache[key] = _response.source_run(self.experiment, nuclide, activity, time, seed)
        return self._cache[key]

    def safety(self) -> dict:
        """For Coulomb excitation, which rings (or strips) of each particle detector see collisions closer than
        Cline's safe distance: {detector: set of first segment indices}. Empty for elastic scattering.

        A ring is marked if the scattered beam or the recoil can reach it from a collision at a CM angle above the
        largest safe one."""
        if "safety" in self._cache:
            return self._cache["safety"]
        from .rates import coulex_for, layer_energy

        r = self._rates()
        out: dict = {}
        ion = beam_ion(self.experiment)
        th_cm = np.linspace(0.5, 180.0, 720)
        for ch in r.channels:
            if ch.excitation is None or not layer_energy(self.experiment, ch.layer, r.layers) > 0:
                continue  # the beam stops before this layer: no collisions there
            safe = coulex_for(self.experiment, ch, r.layers).max_safe_angle()
            if safe >= 180.0:
                continue
            tb = TwoBody(ion, ch.nuclide.name, self.experiment.beam.energy_mev)
            close = th_cm > safe
            lab = np.r_[tb.at_cm(th_cm).theta_lab[close], tb.recoil_for(th_cm).theta_lab[close]]
            for g in r.array:
                for k in sorted({seg[0] for seg in g.segments}):
                    ranges = [g.theta_range(seg) for seg in g.segments if seg[0] == k]
                    lo, hi = min(a for a, _ in ranges), max(b for _, b in ranges)
                    if np.any((lab >= lo) & (lab <= hi)):
                        out.setdefault(g.name, set()).add(k)
        self._cache["safety"] = out
        return out

    def selection(self, key: Optional[str] = None, element: Optional[int] = None) -> dict:
        """What the scene's side panel shows for a selected object: a detector ("detector:N"), one of its rings or
        strips (``element``), a γ-ray detector ("gamma:N"), or with no key (or "target") the whole experiment."""
        r = self._rates()
        what = r.measured
        if key is None or key == "target":
            names = [g.name for g in r.array]
            return {"kind": "experiment", "title": self.experiment.title, "detectors": len(names),
                    "gamma_detectors": len(self.experiment.gamma_detectors),
                    "solid_angle_msr": sum(r.array[n].solid_angle() for n in names),
                    "rate_per_s": sum(r.rate(n) for n in names), "measured": what,
                    "measured_rate_per_s": sum(r.rate(n, what=what) for n in names),
                    "gamma_efficiency": r.gamma_efficiency()[0] if self.experiment.gamma_detectors else None,
                    "target": {"material": self.experiment.target.material,
                               "thickness_um": 1e3 * _scene.target_thickness_mm(self.experiment),
                               "tilt_deg": tilt_deg(self.experiment)},
                    "advice": _scene.advice(self.experiment)}
        i = int(key.split(":")[1])
        if key.startswith("gamma"):
            gd = self.experiment.gamma_detectors[i]
            out = dict(self.live(key), kind="gamma", model=gd.model)
            out["crystals"] = [label or "crystal" for label, _, _ in gd.elements()]
            out["gamma_energy_kev"] = 1e3 * self.gamma_energy_mev()
            out["response"] = _response.Response(self.experiment, gd).describe()
            out["coincidence_rate_per_s"] = None if what == "all" else out["peak_efficiency"] * sum(
                r.rate(g.name, what="excitations") for g in r.array)
            if element is not None and len(out["crystals"]) > 1:
                _, c, rad = gd.elements()[element]
                d = math.hypot(*c)
                out["element"] = {"label": f"crystal {out['crystals'][element]}",
                                  "theta": math.degrees(math.acos(c[2] / d)),
                                  "half_angle_deg": math.degrees(math.atan2(rad, d))}
            return out
        g = r.array.geometries[i]
        hidden = sum(frac for (_, behind), frac in r.array.shadowing().items() if behind == g.name)
        x, y, z = g.centre
        out = {"kind": "detector", "name": g.name, "model": self.experiment.detectors[i].model, "shape": g.shape,
               "theta": math.degrees(math.atan2(math.hypot(x, y), z)), "phi": math.degrees(math.atan2(y, x)),
               "distance_mm": float(np.linalg.norm(g.centre)), "theta_range": g.theta_range(),
               "phi_range": g.phi_range(), "solid_angle_msr": g.solid_angle(), "hidden": min(hidden, 1.0),
               "rate_per_s": r.rate(g.name), "measured": what, "measured_rate_per_s": r.rate(g.name, what=what),
               "counts_in_run": r.counts_in_run(g.name, what=what), "segments": len(g.segments),
               "unsafe": sorted(self.safety().get(g.name, ())), "on_axis": _scene.on_axis(self.experiment, key)}
        if element is not None and len({seg[0] for seg in g.segments}) > 1:
            segs = [seg for seg in g.segments if seg[0] == element]
            ranges = [g.theta_range(seg) for seg in segs]
            omega = g.segment_solid_angles()
            out["element"] = {"label": f"{'strip' if g.shape == 'rectangle' else 'ring'} {element + 1}",
                              "theta_range": (min(a for a, _ in ranges), max(b for _, b in ranges)),
                              "solid_angle_msr": sum(omega[seg] for seg in segs),
                              "rate_per_s": sum(r.rate(g.name, seg) for seg in segs),
                              "measured_rate_per_s": sum(r.rate(g.name, seg, what=what) for seg in segs),
                              "safe": element not in out["unsafe"]}
        return out

    def warnings(self) -> list:
        """Setup problems, then validity and rate warnings: everything the banner should show."""
        cache = self._cache
        if "warnings" not in cache:
            # The valid setup's own warnings are worked out once: the panel, the strip and the read-outs all ask.
            cache["warnings"] = self._setup_warnings()
        order = {"error": 0, "warning": 1, "note": 2}
        return sorted([Warning("error", "setup", p) for p in self.problems] + cache["warnings"],
                      key=lambda w: order[w.level])

    def _setup_warnings(self) -> list:
        out = []
        if self.experiment.fusion is not None:
            from .residues import warnings as fusion_warnings

            plain = self._plain_fusion()
            for text in fusion_warnings(self.experiment, self.residues()):
                if plain is not None and "below the Coulomb barrier" in text:
                    continue  # the fusion model's own (Bass) barrier warning below says it
                level = "warning" if ("below the Coulomb barrier" in text or "needs E*" in text
                                      or "placeholder" in text) else "note"
                out.append(Warning(level, "physics", text))
            for text in (plain or {}).get("warnings", []):
                level = "warning" if ("below the fusion barrier" in text or "quasi-fission" in text) else "note"
                out.append(Warning(level, "physics", text))
        for text in self._rates().warnings():
            level = "warning" if any(k in text for k in ("pile-up", "not Rutherford", "Mott", "stops inside",
                                                         "records nothing", "beam path", "blocks", "hides")) else "note"
            if self.experiment.fusion is not None and ("not Rutherford" in text or "Coulomb barrier" in text):
                # Fusion is meant to happen: for a beam monitor's elastic counts this is a note, not a fault.
                level, text = "note", "For the beam monitor: " + text[:1].lower() + text[1:]
            source = "rates" if ("counts" in text or "collects" in text or "records" in text) else "physics"
            out.append(Warning(level, source, text))
        out += [Warning("warning", "physics", t) for t in self._stopping_shortfalls()]
        if self.experiment.lifetime is not None:
            from .lifetime_plan import deorientation_warning, stopper_warnings

            text = deorientation_warning(self.experiment)
            if text:
                out.append(Warning("warning", "physics", text))
            # The beam in the plunger's foils: its heating and their own Coulomb excitation (kept for the setup,
            # shared with Planner.lifetimes).
            for text in stopper_warnings(self):
                out.append(Warning("note" if text.startswith("No level data") else "warning", "physics", text))
        for p in _scene.problems(self.experiment):
            if p.kind == "overlap":  # a detector in the beam is already reported above
                out.append(Warning("warning", "setup", p.text.replace("would overlap", "overlaps")
                                   + " They cannot both stand there."))
        return out

    def _stopping_shortfalls(self) -> list:
        """One warning per heavy ion (beam or target nuclide, Z > 18) that crosses a layer or detector containing
        hydrogen or helium, where physim's stopping is known to fall short of measured data (see
        :func:`physim.nuclear.stopping.known_shortfall`)."""
        from .stopping import known_shortfall

        e_u = self.experiment.beam.energy_mev / data.nuclide(beam_ion(self.experiment)).A
        ions = {data.nuclide(beam_ion(self.experiment))}
        layers = stack(self.experiment)
        ions |= {data.nuclide(key) for lay in layers for key, n in lay.nuclides.items() if n > 0}
        places = [(lay.material, lay.name) for lay in layers]
        places += [(d.material_data(), d.name) for d in self.experiment.detectors]
        out = []
        for ion in sorted(ions, key=lambda n: (n.Z, n.A)):
            where = [f"{name} ({mat.name})" for mat, name in places if known_shortfall(ion.Z, mat, e_u)]
            if where:
                out.append(f"Energy loss of {ion.name} in {', '.join(dict.fromkeys(where))} may be up to 20% too "
                           f"low: for ions heavier than Ar below 10 MeV/u in hydrogen or helium, physim's stopping "
                           f"falls short of measured data (physics register, stopping). Check the energies after "
                           f"these layers against SRIM or data.")
        return out

    def explain(self, tab: str) -> dict:
        """The "Explain" panel of a tab, with a link to its physics-register page."""
        if tab not in EXPLAIN:
            raise KeyError(f"unknown tab {tab!r}; tabs are {TABS}")
        return dict(EXPLAIN[tab], register=REGISTER[tab])

    def tab(self, name: str, **kw) -> dict:
        """The content of a result tab by name (see :data:`TABS`)."""
        if name not in TABS:
            raise KeyError(f"unknown tab {name!r}; tabs are {TABS}")
        return getattr(self, name)(**kw)

    # -- tabs -------------------------------------------------------------------------------------------------------

    def geometry(self) -> dict:
        """Beam axis, target and detector outlines in 3D (mm), with each detector's solid angle and angles."""
        array = Array.from_experiment(self.experiment)
        extent = max(float(np.linalg.norm(g.centre)) for g in array) * 1.25 if array.geometries else 100.0
        dets = []
        for g in array:
            lo, hi = g.theta_range()
            dets.append({"name": g.name, "outline": g.outline(), "centre": g.centre, "normal": g.n,
                         "solid_angle_msr": g.solid_angle(), "theta_range": (lo, hi), "phi_range": g.phi_range(),
                         "segments": len(g.segments)})
        gammas = []
        gap = np.full((1, 3), np.nan)
        for i, gd in enumerate(self.experiment.gamma_detectors):
            u = np.array(gd.direction())
            dist = gd.distance_mm()
            a, b = (np.array(x) for x in gd.face_axes())
            t = np.linspace(0.0, 2 * np.pi, 49)[:, None]
            # The front face of each crystal, and the crystal's length behind it where the setup gives one.
            length = _q(gd.crystal_length).to("mm") if gd.crystal_length is not None else 0.0
            parts = []
            for _, centre, rad in gd.elements():
                ring = np.array(centre) + rad * (np.cos(t) * a + np.sin(t) * b)
                parts += [ring, gap]
                if length:
                    parts += [ring + length * u, gap]
                    for k in (0, 12, 24, 36):
                        parts += [np.array([ring[k], ring[k] + length * u]), gap]
            housing = gd.housing()
            if housing is not None:
                c, side = np.array(housing[0]), housing[1] / 2
                corners = [c + sa * side * a + sb * side * b for sa, sb in ((1, 1), (-1, 1), (-1, -1), (1, -1), (1, 1))]
                parts += [np.array(corners), gap]
            gammas.append({"name": gd.name or f"γ{i + 1}", "outline": np.concatenate(parts[:-1]),
                           "centre": dist * u, "theta": _q(gd.theta).to("deg"),
                           "half_angle_deg": gd.half_angle_deg(), "distance_mm": dist, "model": gd.model,
                           "crystals": len(gd.elements())})
        if gammas:
            extent = max(extent, 1.25 * max(g["distance_mm"] for g in gammas))
        blocking = [{"name": b.name, "outline": b.outline()} for b in array.blockers if "(housing)" not in b.name]
        return {"beam": np.array([[0.0, 0.0, -extent], [0.0, 0.0, extent]]), "extent": extent,
                "target_size_mm": 0.04 * extent, "detectors": dets, "gamma_detectors": gammas,
                "blocking": blocking}

    def kinematics(self, points: int = 361) -> dict:
        """Lab energy against lab angle for the scattered beam and the recoil of every target nuclide, with the
        angular range each detector covers."""
        beam = self.experiment.beam
        ion = beam_ion(self.experiment)
        th_cm = np.linspace(0.0, 180.0, points)
        curves = []
        for lay in stack(self.experiment):
            for key in lay.nuclides:
                nuc = data.nuclide(key).name
                tb = TwoBody(ion, nuc, beam.energy_mev)
                ej = tb.at_cm(th_cm)
                rec = tb.recoil_for(th_cm)
                label = f"{nuc} ({lay.name})"
                curves.append({"label": f"{ion} on {label}", "particle": "ejectile", "theta": ej.theta_lab,
                               "energy": ej.energy, "max_angle": tb.max_angle()})
                curves.append({"label": f"{nuc} recoil ({lay.name})", "particle": "recoil", "theta": rec.theta_lab,
                               "energy": rec.energy, "max_angle": tb.max_angle("recoil")})
        for ch in self._rates().channels:
            if ch.excitation is None:
                continue
            try:
                tb = TwoBody(ion, ch.nuclide.name, beam.energy_mev, **ch.kinematics_args)
            except ValueError:  # below the excitation threshold
                continue
            ej, rec = tb.at_cm(th_cm), tb.recoil_for(th_cm)
            star = f"{ch.excitation.energy_mev * 1e3:g} keV"
            who = ch.nuclide.name if ch.excitation.excite == "target" else ion
            curves.append({"label": f"{ion} on {ch.nuclide.name}, {who}* {star}", "particle": "ejectile",
                           "theta": ej.theta_lab, "energy": ej.energy, "max_angle": tb.max_angle(),
                           "inelastic": True})
            curves.append({"label": f"{ch.nuclide.name} recoil, {who}* {star}", "particle": "recoil",
                           "theta": rec.theta_lab, "energy": rec.energy, "max_angle": tb.max_angle("recoil"),
                           "inelastic": True})
        cover = [{"name": g.name, "theta_range": g.theta_range()} for g in Array.from_experiment(self.experiment)]
        return {"beam_energy_mev": beam.energy_mev, "curves": curves, "detectors": cover}

    def rates(self) -> dict:
        """One row per detector: solid angle, angles, rate, counts in the run, beam time for the counts wanted,
        relative error; and per-strip rates for each detector."""
        r = self._rates()
        what = r.measured
        rows, strips = [], {}
        for g in r.array:
            rate = r.rate(g.name)
            counts = r.counts_in_run(g.name, what=what)
            row = {"detector": g.name, "solid_angle_msr": g.solid_angle(), "theta_range": g.theta_range(),
                   "rate_per_s": rate, "rate_all_per_s": r.rate(g.name, counted=False),
                   "counts_in_run": counts,
                   "relative_error": r.relative_error(g.name, what=what) if counts > 0 else math.inf,
                   "beam_time_s": (r.beam_time_for(g.name, what=what) if r.counts_wanted is not None else None)}
            if what != "all":
                row["excitation_per_s"] = r.rate(g.name, what="excitations")
                row["coincidence_per_s"] = r.rate(g.name, what="coincidences") if what == "coincidences" else None
                # The γ-ray efficiency in coincidence with this detector, with the angular correlation.
                row["gamma_efficiency"] = r.gamma_efficiency(g.name)[0] if what == "coincidences" else None
            rows.append(row)
            strips[g.name] = r.per_segment(g.name)
        eff, typical = r.gamma_efficiency() if what == "coincidences" else (None, False)
        return {"rows": rows, "strips": strips, "beam_time_s": r.beam_time_s, "counts_wanted": r.counts_wanted,
                "particles_per_second": r.particles_per_second, "measured": what, "gamma_efficiency": eff,
                "gamma_efficiency_typical": typical}

    def plan(self) -> dict:
        """What one checks before asking for beam, answers first (backlog item 65): the counts wanted and the beam
        time they need against the beam time planned; per particle detector its rate, the share of its time the
        acquisition is busy, its share of the excitations and its safe and unsafe rings; per γ-ray detector its
        full-energy efficiency at the transition and its coincidence rate; and the particle × γ coincidence rates.
        Plain numbers (no arrays), so a run's summary can keep them."""
        r = self._rates()
        t = self.rates()
        exp = self.experiment
        dead = _q(exp.run.dead_time).to("s") if exp.run.dead_time is not None else 0.0
        total_rate = sum(row["rate_per_s"] for row in t["rows"])
        excitations = {row["detector"]: row.get("excitation_per_s") or 0.0 for row in t["rows"]}
        all_exc = sum(excitations.values())
        unsafe = self.safety()
        detectors = []
        other = r.other_peak()
        for row, g in zip(t["rows"], r.array):
            rings = sorted({seg[0] for seg in g.segments})
            bad = sorted(unsafe.get(g.name, ()))
            busy = row["rate_per_s"] * dead
            detectors.append({
                "other_coincidence_per_s": other["coincidence_per_s"].get(g.name) if other else None,
                "detector": row["detector"], "theta_range": [float(x) for x in row["theta_range"]],
                "solid_angle_msr": row["solid_angle_msr"], "rate_per_s": row["rate_per_s"],
                "dead_time_fraction": busy / (1 + busy), "excitation_per_s": excitations[row["detector"]],
                "excitation_share": excitations[row["detector"]] / all_exc if all_exc > 0 else None,
                "coincidence_per_s": row.get("coincidence_per_s"), "counts_in_run": row["counts_in_run"],
                "beam_time_s": row["beam_time_s"], "relative_error": row["relative_error"],
                "rings": len(rings), "unsafe_rings": bad, "safe_rings": [k for k in rings if k not in bad]})
        gammas, matrix = [], {}
        if exp.excitation is not None and exp.gamma_detectors:
            r.gamma_efficiency()
            effs = r._gamma_effs
            factors = r.correlation()
            for i, (gd, eff) in enumerate(zip(exp.gamma_detectors, effs)):
                name = gd.name or f"G{i + 1}"
                pairs = {d: excitations[d] * eff * factors.get(d, {}).get(name, 1.0) for d in excitations}
                for d, v in pairs.items():
                    matrix.setdefault(d, {})[name] = v
                gammas.append({"detector": name, "efficiency": eff, "energy_kev": 1e3 * exp.excitation.energy_mev,
                               "coincidence_per_s": sum(pairs.values())})
        timed = [d for d in detectors if d["beam_time_s"] is not None]
        needed = max(timed, key=lambda d: d["beam_time_s"]) if timed else None
        residues = self.residues() if exp.fusion is not None else None
        fusion = self._plain_fusion() if exp.fusion is not None else None
        needed_s = needed["beam_time_s"] if needed else None
        study = None
        if residues is not None and residues["rows"]:
            # Fusion-evaporation: the beam time is the one the studied channel's strongest line needs.
            i = self.studied_channel()
            row = residues["rows"][i if i is not None and i < len(residues["rows"]) else 0]
            line = row["lines"][0] if row["lines"] else None
            study = {"index": residues["rows"].index(row), "channel": row["channel"], "residue": row["residue"],
                     "cross_section_mb": row["cross_section_mb"], "share": row["share"],
                     "rate_per_s": row["rate_per_s"], "beam_time_s": row["beam_time_s"],
                     "line": (None if line is None else
                              {"energy_kev": line["energy_kev"], "label": line["label"],
                               "detected_per_s": line["detected_per_s"]}),
                     "chosen": exp.run.channel is not None}
            needed, needed_s = None, row["beam_time_s"]
        stops_in = next((lay.name for j, lay in enumerate(r.layers)
                         if float(beam_energy_at(exp, j, [lay.thickness], r.layers)[0]) <= 0), None)
        return {"measured": t["measured"], "counts_wanted": t["counts_wanted"], "beam_time_s": t["beam_time_s"],
                # The layer the beam stops in (a thick target or backing), or None when it goes through.
                "beam_stops_in": stops_in,
                "beam_time_needed_s": needed_s,
                "limiting_detector": needed["detector"] if needed else None,
                "enough": (needed_s <= t["beam_time_s"]) if needed_s is not None else None,
                # Fusion-evaporation: the channel studied, whose strongest line sets the beam time.
                "study": study,
                "particles_per_second": t["particles_per_second"], "dead_time_s": dead,
                "live_fraction": 1 / (1 + dead * total_rate), "detectors": detectors, "gamma_detectors": gammas,
                "coincidences": matrix, "gamma_efficiency": t["gamma_efficiency"],
                # The other nucleus's peak, which the events hold when the setup has its level scheme.
                "other_peak": ({"label": other["label"], "nuclide": other["nuclide"],
                                "energy_kev": other["energy_kev"]} if other else None),
                # Fusion-evaporation: the residues of the channel table (backlog item 75).
                "residues": residues, "fusion": fusion}

    def fusion(self) -> dict:
        """The fusion of the setup's beam and target (:func:`physim.nuclear.exit_channels.fusion_summary`): the
        barrier, σ_fus at the beam energy and through the target, the excitation function and the model's
        warnings. Kept until the setup changes."""
        if "fusion" not in self._cache:
            from .exit_channels import fusion_summary

            self._cache["fusion"] = fusion_summary(self.experiment)
        return self._cache["fusion"]

    def _plain_fusion(self) -> Optional[dict]:
        """:meth:`fusion` without the excitation function's arrays, for the plan (None if it cannot be worked
        out; the warnings then say why the setup is incomplete)."""
        try:
            f = self.fusion()
        except Exception:  # noqa: BLE001 - the plan still comes back
            return None
        return {k: (tuple(float(x) for x in v) if isinstance(v, tuple) else
                    (float(v) if isinstance(v, (float, np.floating)) else v))
                for k, v in f.items() if k != "excitation_function"}

    def residues(self) -> dict:
        """The fusion-evaporation residues (:func:`physim.nuclear.residues.residue_table`): per channel its cross
        section, residues per second and in the run, and its strongest lines. Kept until the setup changes."""
        if self.experiment.fusion is None:
            raise ValueError("the setup is not a fusion-evaporation reaction")
        if "residues" not in self._cache:
            from .residues import residue_table

            self._cache["residues"] = residue_table(self.experiment)
        return self._cache["residues"]

    def energy_loss(self) -> dict:
        """The beam through each layer (energy in and out, loss, straggling), the beam energy through the target,
        and each detector's dead-layer loss and punch-through energy for the scattered beam."""
        exp = self.experiment
        ion = beam_ion(exp)
        layers = stack(exp)
        rows = []
        e = exp.beam.energy_mev
        c = math.cos(math.radians(tilt_deg(exp)))
        profile_z, profile_e, offset = [], [], 0.0
        for i, lay in enumerate(layers):
            st = stopping(ion, lay.material)
            path = lay.thickness / c
            out = st.energy_after(e, path)
            rows.append({"layer": lay.name, "material": lay.material.name, "thickness_mg_cm2": lay.thickness,
                         "energy_in_mev": e, "energy_out_mev": out, "loss_mev": e - out,
                         "straggling_fwhm_mev": 2.3548 * st.straggling(e, path, steps=50) if out > 0 else 0.0})
            z = np.linspace(0.0, lay.thickness, 41)
            profile_z.append(offset + z)
            profile_e.append(beam_energy_at(exp, i, z, layers))
            offset += lay.thickness
            e = out
        main = data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        dets = []
        for det, g in zip(exp.detectors, Array.from_experiment(exp)):
            resp = Response(ion, det)
            k = float(TwoBody(ion, main, exp.beam.energy_mev).at_lab(g.mean_theta())[0].energy)
            dets.append({"detector": g.name, "ejectile_energy_mev": k,
                         "after_dead_layer_mev": float(resp.after_dead_layer(k)) if np.isfinite(k) else math.nan,
                         "punch_through_mev": resp.punch_through_energy()})
        return {"layers": rows, "depth_mg_cm2": np.concatenate(profile_z), "beam_energy_mev": np.concatenate(profile_e),
                "detectors": dets}

    def spectra(self, events: int = 200_000, seed: int = 1, bins: int = 200) -> dict:
        """Measured-energy spectra per detector: of the current run's real part (counts per bin), or without a
        run simulated (counts in the planned beam time per bin). A fusion run holds the residues' γ rays and no
        particle events: the spectra are then the setup's simulated ones, and ``note`` says so."""
        fusion = self._beam_run() is not None and self.run.kind == "fusion"
        if self._beam_run() and not fusion:
            ev = self.run.events()
        else:
            key = ("events", events, seed)
            if key not in self._cache:
                from .events import simulate

                self._cache[key] = simulate(self.experiment, events, seed)
            ev = self._cache[key]
        out = {}
        for name in ev.detectors:
            m = ev.select(name)
            if not m.any():
                out[name] = {"counts": np.zeros(bins), "edges": np.linspace(0, 1, bins + 1)}
                continue
            hi = float(ev["measured"][m].max()) * 1.03
            h, edges = ev.spectrum(name, bins=bins, range=(0.0, hi))
            out[name] = {"counts": h, "edges": edges}
        result = {"events": ev, "spectra": out, "run": self.run_status() if self._beam_run() and not fusion else None}
        if fusion:
            result["note"] = ("The current run is a fusion run: it holds the residues' γ rays and no particle events, "
                              "so these are the setup's simulated particle spectra (counts in the planned beam time).")
        return result

    def gamma_events(self, events: int = 400_000, seed: int = 1, plain: bool = False):
        """γ rays in coincidence with the detected particles (:func:`physim.nuclear.gamma_events.simulate_gammas`),
        built on the same particle events as :meth:`spectra`; ``plain`` leaves out the add-back and the shields.
        Kept until the setup changes. With a current beam run, its γ rays (nothing is simulated)."""
        from .events import simulate
        from .gamma_events import simulate_gammas

        if self._beam_run():
            g = self.run.gammas(plain)
            if g is None:
                raise ValueError("the current run has no γ rays")
            return g
        key = ("gammas", events, seed, plain)
        if key not in self._cache:
            ev_key = ("events", events, seed)
            if ev_key not in self._cache:
                self._cache[ev_key] = simulate(self.experiment, events, seed)
            self._cache[key] = simulate_gammas(self.experiment, events, seed, particle_events=self._cache[ev_key],
                                               plain=plain)
        return self._cache[key]

    def gamma_modes(self) -> dict:
        """Per γ-ray detector, {"addback": bool, "shield": material or None}: which run with add-back or a
        Compton-suppression shield."""
        return {gd.name or f"G{i + 1}": {"addback": bool(gd.addback), "shield": gd.shield}
                for i, gd in enumerate(self._data_experiment().gamma_detectors)}

    def gamma_spectra(self, events: int = 400_000, seed: int = 1, bins: int = 300) -> dict:
        """The γ-ray side of the Monte Carlo: for each γ-ray detector the raw and Doppler-corrected spectra in
        coincidence with all particle detectors (counts in the run per bin, random coincidences included), the
        particle × γ matrix of true and random coincidences, and the live fraction.
        ``{"available": False}`` without an excited state or γ-ray detectors. With a current beam run, its real part
        (``"run"`` says which run, and whether the setup changed since)."""
        exp = self._data_experiment()
        if exp.excitation is None or not exp.gamma_detectors:
            return {"available": False, "reason": "Coulomb excitation with γ-ray detectors is needed."}
        from .gamma_events import spectrum_top

        g = self.gamma_events(events, seed)
        emitter = "recoil" if exp.excitation.excite == "target" else "projectile"
        # To the top of the level scheme, at the bin width ``bins`` gives over 1.3 × the transition energy.
        top = max(spectrum_top(exp), 1.3 * g.energy_mev)
        n_bins = int(round(bins * top / (1.3 * g.energy_mev)))

        def spectrum_set(gg, name: str) -> dict:
            out = {}
            for key in (None, "projectile", "recoil"):
                parts = [gg.spectrum(name, d, corrected=key, bins=n_bins, range=(0.0, top))
                         for d in gg.events.detectors]
                out[key or "measured"] = np.sum([h for h, _ in parts], axis=0)
                out["edges"] = parts[0][1]
            return {**out, "counts": gg.counts(name)}

        modes = self.gamma_modes()
        spectra = {name: spectrum_set(g, name) for name in g.detector_names()}
        # A detector with add-back or a shield also shows its spectrum without them, on the same γ rays.
        if any(m["addback"] or m["shield"] for m in modes.values()):
            plain = self.gamma_events(events, seed, plain=True)
            for name, m in modes.items():
                if (m["addback"] or m["shield"]) and name in spectra:
                    spectra[name]["plain"] = spectrum_set(plain, name)
        return {"available": True, "spectra": spectra, "emitter": emitter, "energy_kev": 1e3 * g.energy_mev,
                "coincidences": g.coincidences(), "live_fraction": g.live_fraction, "window_s": g.window_s,
                "gamma_detectors": g.detector_names(), "particle_detectors": list(g.events.detectors),
                "notes": g.notes, "n_gammas": len(g), "modes": modes,
                "run": self.run_status() if self._beam_run() else None}

    def channel_selection(self, events: Optional[int] = None, seed: int = 1) -> dict:
        """The channel-selection matrix of the neutron wall and the particle detectors
        (:func:`physim.nuclear.exit_channels.channel_selection`). Kept until the setup changes."""
        from .exit_channels import SELECTION_EVENTS, channel_selection

        key = ("channel_selection", events or SELECTION_EVENTS, seed)
        if key not in self._cache:
            self._cache[key] = channel_selection(self.experiment, events or SELECTION_EVENTS, seed)
        return self._cache[key]

    def add_neutron_wall(self, **fields) -> bool:
        """Give a fusion-evaporation setup a neutron wall (``[neutron_wall]``); without fields, 50 typical EJ-301
        cells 510 mm from the target."""
        d = self.draft
        d["neutron_wall"] = dict(fields) or {"cells": 50, "distance": "510 mm"}
        return self._apply(d)

    def remove_neutron_wall(self) -> bool:
        d = self.draft
        d.pop("neutron_wall", None)
        return self._apply(d)

    #: What an ancillary detector of a fusion-evaporation setup starts with when it is added.
    ANCILLARY_DEFAULTS = {
        "neutron_wall": {"cells": 50, "distance": "510 mm", "theta_min": "4 deg"},
        "charged_particle_ball": {"elements": 40, "radius": "60 mm", "material": "CsI", "threshold": "2 MeV"},
        "recoil_separator": {"theta_acceptance": "3 deg", "kind": "gas-filled", "distance": "600 mm",
                             "momentum_acceptance": "10 %"},
        "gamma_calorimeter": {"elements": 60, "radius": "120 mm", "material": "BGO", "thickness": "75 mm",
                              "fold_min": 6},
    }

    def add_ancillary(self, name: str, **fields) -> bool:
        """Give a fusion-evaporation setup an ancillary detector: "neutron_wall", "charged_particle_ball",
        "recoil_separator" or "gamma_calorimeter" (:data:`ANCILLARY_DEFAULTS` where no fields are given).
        "plunger" gives the setup (Coulex too) a lifetime section: a plunger with the default distances and a
        stopper (:meth:`set_lifetime`), or the ``[lifetime]`` fields given."""
        if name == "plunger":
            d = self.draft
            d["lifetime"] = copy.deepcopy(fields) or self._default_lifetime("rdds")
            return self._apply(d)
        if name not in self.ANCILLARY_DEFAULTS:
            known = ", ".join(list(self.ANCILLARY_DEFAULTS) + ["plunger"])
            raise KeyError(f"unknown ancillary {name!r}; known: {known}")
        d = self.draft
        d[name] = dict(fields) or dict(self.ANCILLARY_DEFAULTS[name])
        return self._apply(d)

    # -- the lifetime section (plunger or DSAM) ---------------------------------------------------------------------

    #: The distances a new plunger starts with: log-spaced, for lifetimes of about 1 to 100 ps at a few % of c.
    PLUNGER_DISTANCES = ("10 um", "20 um", "50 um", "100 um", "200 um", "500 um")

    def _default_lifetime(self, method: str) -> dict:
        """A new lifetime section: DSAM has no fields of its own; a plunger gets :data:`PLUNGER_DISTANCES` and a
        gold stopper thick enough to stop the setup's recoil (half as thick again as its range, at least
        10 mg/cm²)."""
        if method != "rdds":
            return {"method": method}
        return {"method": "rdds", "distances": list(self.PLUNGER_DISTANCES),
                "stopper": {"kind": "stopper", "material": "Au", "thickness": self._stopper_thickness("Au")}}

    def _stopper_thickness(self, material: str) -> str:
        from .experiment import _stopping

        recoil = self.experiment.lifetime_recoil()
        need = 0.0
        if recoil is not None and recoil[2] > 0:
            try:
                need = 1.5 * float(_stopping(recoil[0], material, None).range(recoil[2]))
            except ValueError:
                need = 0.0
        return f"{max(10, math.ceil(need))} mg/cm2"

    @staticmethod
    def _text(value, unit: str = "um") -> str:
        """A quantity as setup text: text and quantities as they are, a bare number in ``unit``."""
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return f"{format_number(float(value))} {unit}"
        return str(value)

    def set_lifetime(self, method: Optional[str] = None, distances=None, offset=None, deorientation=None) -> bool:
        """Set the lifetime measurement, making the section if there is none (a plunger unless ``method`` says
        "dsam"). ``distances`` are the target-to-stopper distances and ``offset`` the distance at electrical
        contact (text with a unit, or numbers in µm); ``deorientation`` a dict of ``alpha2``, ``alpha4`` and either
        ``gamma2``, ``gamma4`` or ``g``, ``c2``, ``c4`` (an empty dict removes it). Switching a plunger to DSAM
        keeps its distances and foils aside, and switching back brings them back. Returns whether the setup is
        valid afterwards."""
        d = self.draft
        lt = d.get("lifetime")
        if lt is None:
            lt = d["lifetime"] = self._default_lifetime(method or "rdds")
        elif method is not None and method != lt.get("method"):
            if method == "dsam":
                self._stashed_plunger = {k: lt.pop(k) for k in PLUNGER_FIELDS if k in lt}
                lt["method"] = method
            else:
                lt["method"] = method
                if method == "rdds" and not any(k in lt for k in PLUNGER_FIELDS):
                    back = copy.deepcopy(self._stashed_plunger) or self._default_lifetime("rdds")
                    lt.update({k: v for k, v in back.items() if k != "method"})
        if distances is not None:
            lt["distances"] = [self._text(x) for x in distances]
        if offset is not None:
            lt["offset"] = self._text(offset)
        if deorientation is not None:
            if deorientation:
                lt["deorientation"] = {k: v if isinstance(v, (int, float)) else str(v)
                                       for k, v in deorientation.items()}
            else:
                lt.pop("deorientation", None)
        return self._apply(d)

    def set_stopper(self, kind: str, material: str, thickness, density=None) -> bool:
        """Set the plunger's foil after the gap: ``kind`` "stopper" (stops the recoil), "degrader" (lets it through,
        slower) or "differential" (the first of two foils; see :meth:`set_degrader2`), its material and thickness
        ("10 mg/cm2", or a length with a known density). Makes a plunger section if there is none. A differential
        plunger's second foil is kept aside while the kind is another, and comes back with "differential"."""
        d = self.draft
        lt = d.setdefault("lifetime", {"method": "rdds", "distances": list(self.PLUNGER_DISTANCES)})
        stopper = {"kind": kind, "material": material, "thickness": self._text(thickness, "mg/cm2")}
        if density is not None:
            stopper["density"] = self._text(density, "g/cm3")
        lt["stopper"] = stopper
        if kind != "differential" and "degrader2" in lt:
            self._stashed_degrader2 = lt.pop("degrader2")
        elif kind == "differential" and "degrader2" not in lt and self._stashed_degrader2:
            lt["degrader2"] = copy.deepcopy(self._stashed_degrader2)
        return self._apply(d)

    def set_degrader2(self, material: str, thickness, gap, density=None) -> bool:
        """Set the second foil of a differential plunger, ``gap`` behind the first (text with a unit, or µm); it
        stops the recoil. The plunger's first foil becomes of kind "differential"."""
        d = self.draft
        lt = d.setdefault("lifetime", {"method": "rdds", "distances": list(self.PLUNGER_DISTANCES)})
        foil = {"material": material, "thickness": self._text(thickness, "mg/cm2"), "gap": self._text(gap)}
        if density is not None:
            foil["density"] = self._text(density, "g/cm3")
        lt["degrader2"] = foil
        if isinstance(lt.get("stopper"), dict):
            lt["stopper"]["kind"] = "differential"
        return self._apply(d)

    def remove_lifetime(self) -> bool:
        """Take the lifetime measurement out of the setup."""
        d = self.draft
        d.pop("lifetime", None)
        return self._apply(d)

    # -- planning, predicting and analysing a lifetime measurement (lifetime wave 3A) ------------------------------

    def lifetimes(self, level=None, precision: float = 0.1) -> dict:
        """How to measure the lifetimes of the studied nucleus's levels (:mod:`physim.nuclear.lifetime_plan`):
        Coulomb excitation's excited nucleus (the projectile or the target), or fusion-evaporation's studied residue
        (:meth:`studied_channel`; its levels as :meth:`band` lists them).

        Returns ``available`` (and ``reason`` when not), ``reaction``, ``nucleus``, ``method`` (the setup's, or
        None), ``distances_um`` (the setup's, without the offset), ``precision`` (the relative precision on τ
        wanted, default 10%), ``time_to_stop_ps`` (DSAM: the recoil's time to stop in the target and backing; a
        plunger: in its stopper, the DSAM-like tail under the unshifted peak), ``problems``, ``warnings`` and
        ``levels``, the studied level first (``level`` asks for one level only: its name, number or energy), each
        with:

        - ``level``, ``energy_kev``, ``gamma_kev`` (its strongest line), ``tau_ps`` and ``tau_source``
          (:func:`physim.nuclear.lifetimes.level_taus`: the half-life, else the B(Eλ↓) of the scheme or the
          setup's Coulex B(E2), with its source), ``beta`` (the recoil's speed: in the plunger's gap, or where it
          is made for DSAM);
        - ``method``: "rdds" for τ of about 1 ps to 1 ns, "dsam" for 10 fs to 1 ps, else "outside" (also for a
          level without a known τ); "none" for every level when the recoil does not fly forward out of the target
          (β = 0: no Doppler-shift method fits; the warnings say why);
        - ``recommended_distances_um`` and ``sensitive_um`` (:func:`~physim.nuclear.lifetimes.recommend_distances`,
          :func:`~physim.nuclear.lifetimes.sensitive_range`) for a plunger level; empty in DSAM mode (the setup's
          DSAM, or a recoil that stops in the target);
        - ``separation``: per γ-ray detector its angle, the shifted − unshifted separation against its FWHM, and
          whether they are ``resolved`` (at least two FWHM apart);
        - ``attenuation``: F(τ) on the setup's backing for DSAM (:func:`~physim.nuclear.lifetimes.dsam_attenuation`
          on the slowing-down kernel of the typical recoil), else None;
        - ``hours_per_distance``: the beam time at each distance (the setup's, else the recommended ones) for the
          wanted precision on τ, from the line's full-energy rate in all γ-ray detectors (``rate_per_s``) and the
          one-level binomial precision of the unshifted share at each distance, combined over the distances;
        - ``backing_advice`` (levels in the DSAM range, and every level with a τ in a DSAM setup; else None): the
          usual backings ranked by the DSAM sensitivity ``|dF/d ln τ|``, each with F(τ), the slowing-down time and the
          thickness that stops the recoil, and the setup's own backing among them
          (:func:`physim.nuclear.lifetime_plan.backing_advice`).

        Also ``stopper``: for a plunger, the beam in each foil it reaches (energy in and out, power deposited
        against the foil's limit, and the foil's own Coulex lines near the studied peaks); its warnings are in
        ``warnings``. The recoil is the setup's typical one (:meth:`~physim.nuclear.experiment.Experiment.lifetime_recoil`).
        Kept until the setup changes."""
        key = ("lifetimes", level, float(precision))
        if key not in self._cache:
            from . import lifetime_plan

            self._cache[key] = lifetime_plan.plan(self, level, precision)
        return self._cache[key]

    def lifetime_spectra(self, distance_um=None, level=None) -> dict:
        """The predicted γ-ray spectra of a level's line (``level``; the studied one by default) at each plunger
        distance (``distance_um``: one distance, else the setup's) and each ring of γ-ray detectors (their distinct
        polar angles), or for DSAM the line shape at each ring.

        Counts are expected full-energy counts in the beam time per distance (the setup's beam time shared over its
        distances; DSAM: all of it). How they are made (``model``):

        - fusion-evaporation: :func:`physim.nuclear.residues.plunger_spectra` or
          :func:`~physim.nuclear.residues.dsam_lineshapes` (residues drawn from the evaporation, the cascade with
          side feeding), scaled by the residue rate and the ring's full-energy efficiency, folded with its
          resolution;
        - Coulomb excitation, a plunger: analytic, fast enough for a UI call. The unshifted share is
          :func:`~physim.nuclear.lifetimes.unshifted` of the level's cascade at the flight time of the typical
          recoil (one speed, head-on); the shifted peak is a Gaussian at its Doppler energy, as wide as the
          detector's resolution and the spread of angles over a crystal's face. The real recoils' spread of speeds
          and angles (forward-scattered beam caught by a 0° detector, say) broadens and moves the shifted peak: a
          run shows that;
        - Coulomb excitation, DSAM: :func:`~physim.nuclear.lifetimes.lineshape` on the slowing-down kernel of the
          typical recoil (ZBL scattering), folded with the resolution.

        Returns ``available``, ``reason``, ``method``, ``edges_kev``, ``spectra`` = [{``distance_um`` (None for
        DSAM), ``theta_deg``, ``counts``, ``shifted_kev``, ``unshifted_kev``}], and ``level``, ``gamma_kev``,
        ``model``, ``seconds_per_distance``, ``notes``. Kept until the setup changes."""
        key = ("lifetime_spectra", distance_um, level)
        if key not in self._cache:
            from . import lifetime_plan

            self._cache[key] = lifetime_plan.spectra(self, distance_um, level)
        return self._cache[key]

    def analyse_lifetimes(self, level=None, series=None, method: Optional[str] = None) -> dict:
        """τ of a level (the studied one by default) from the runs, as an experimenter would get it, then B(λ↓) and
        the comparison with what is known.

        A plunger: ``series`` (a :class:`~physim.nuclear.runs.DistanceSeries`; by default this setup's, from
        :meth:`~physim.nuclear.runs.ExperimentFolder.distance_series`) → the spectra per distance and ring
        (:func:`physim.nuclear.dataviews.distance_spectra`, Doppler shifts kept) → the shifted and unshifted peak
        areas in the rings away from 90°, added (:func:`physim.nuclear.lifetime_analysis.peak_areas`; the recoil's
        speed from where the shifted peaks are) → ``method`` "ddcm" (default; with the strongest feeding line's
        areas for the feeding correction), "ddcm_three_foil" (the default for a differential plunger) or
        "decay_curve". DSAM: the runs of this setup → the line shape at each ring away from 90° →
        :func:`~physim.nuclear.lifetime_analysis.fit_lineshape` on the transport kernel ("lineshape"), τ the
        weighted mean of the rings. Then B(λ↓) by :func:`~physim.nuclear.lifetime_analysis.b_from_result` and
        ``compare``: the level scheme's own τ (ENSDF, typed, derived), the local ENSDF copy, the setup's Coulex
        B(E2) and shell-model matrix elements in the scheme.

        Returns ``available``, ``reason``, ``level``, ``method``, ``tau_ps``, ``stat_ps``, ``budget``, ``points`` =
        [{``distance_um``, ``unshifted``, ``shifted`` (shares of the line), ``error``, ``tau_ps``, ``tau_err_ps``}]
        (empty for DSAM), ``curve`` = {``x``, ``y``, ``kind``} (the fitted curve, or the line shape), ``b_down`` =
        [value, uncertainty], ``b_unit``, ``compare`` = [{``source``, ``tau_ps``, ``b_down``}], ``notes``, and
        ``total_ps``, ``beta``, ``rings_deg``, and ``transition`` (the nucleus and the two levels of the analysed
        line) and ``runs`` (the names of the runs it read).

        An analysis that gave a τ is written to the experiment's results ledger (τ, and the B(E2) that depends on
        it: :func:`physim.nuclear.result_rows.write_lifetimes`), and, when it is of this setup's own runs, kept in
        the experiment's folder (``lifetime.json``, :meth:`kept_lifetimes`)."""
        from . import lifetime_plan, result_rows

        out = lifetime_plan.analyse(self, level, series, method)
        if out.get("available"):
            if series is None:
                self._keep("lifetime.json", out, self.experiment_folder)
            result_rows.write_lifetimes(self, out)  # never fails the analysis
        return out

    def kept_lifetimes(self) -> Optional[dict]:
        """The last lifetime analysis of the current setup's runs (the dict of :meth:`analyse_lifetimes`), as
        kept in the experiment's folder (``lifetime.json``): there after the app is opened again. None if there
        is none, or if the setup's physics changed since. Never starts the analysis."""
        return self._kept("lifetime.json", self.experiment_folder)

    def remove_ancillary(self, name: str) -> bool:
        d = self.draft
        d.pop(name, None)
        return self._apply(d)

    def studied_channel(self) -> Optional[int]:
        """The index of the exit channel a fusion-evaporation setup studies (``[run] channel``, else the largest
        cross section), or None for another reaction."""
        exp = self.experiment
        if exp.fusion is None:
            return None
        from .experiment import studied_channel_index

        return studied_channel_index(exp.fusion.channels, exp.run.channel)

    def band(self, index: Optional[int] = None, lines: int = 12) -> dict:
        """The studied residue (``index``: a channel of the table, else the studied one) as the run sees it: its
        strongest γ-ray transitions with their rates in the γ-ray detectors (from
        :func:`physim.nuclear.residues.residue_table`), the lifetime of the level each leaves, and B(E2)↓ in
        Weisskopf units where the level scheme determines it (its half-life, branching and conversion, by
        :meth:`physim.nuclear.levels.LevelScheme.derive`). Known values, beside the rates: nothing here is
        measured by the run. Kept until the setup changes."""
        exp = self.experiment
        if exp.fusion is None:
            raise ValueError("the setup is not a fusion-evaporation reaction")
        i = self.studied_channel() if index is None else index
        key = ("band", i, lines)
        if key in self._cache:
            return self._cache[key]
        from .residues import residue_scheme, residue_table

        ch = exp.fusion.channels[i]
        residue = ch.residue_name
        scheme = residue_scheme(exp, residue)
        if scheme is None:
            out = {"available": False, "residue": residue, "channel": ch.particles,
                   "reason": f"No level scheme of {residue}: neither the setup nor the local ENSDF copy has one."}
            self._cache[key] = out
            return out
        table = residue_table(exp, lines=lines)
        row = table["rows"][i]
        rates = {round(x["energy_kev"], 1): x["detected_per_s"] for x in row["lines"]}
        derived = copy.deepcopy(scheme)
        derived.derive()
        rows = []
        for t in derived.transitions:
            e = round(t.energy.value, 1)
            if e not in rates:
                continue
            up, low = derived.levels[t.initial], derived.levels[t.final]
            b = me = None
            if t.multipolarity in ("", "E2") or "E2" in t.multipolarity:
                b = derived.b_weisskopf(t.initial, t.final, "E2")
                me = derived.matrix_element(t.initial, t.final, "E2")
            rel = (2 * me.unc / me.value) if (me is not None and me.unc and me.value) else None
            tau = up.half_life.value / math.log(2) if up.half_life is not None else None
            rows.append({"initial": up.label or up.jpi, "final": low.label or low.jpi,
                         "level_kev": up.energy.value, "energy_kev": t.energy.value, "rate_per_s": rates[e],
                         "multipolarity": t.multipolarity, "tau_ps": None if tau is None else 1e12 * tau,
                         "b_e2_wu": b, "b_e2_unc_wu": None if (b is None or rel is None) else b * rel})
        rows = sorted(sorted(rows, key=lambda r: -r["rate_per_s"])[:10], key=lambda r: r["level_kev"])
        rates = {round(r["energy_kev"], 1): r["rate_per_s"] for r in rows}
        levels = sorted({(r["level_kev"], r["initial"]) for r in rows} | {
            (derived.levels[t.final].energy.value, derived.levels[t.final].label or derived.levels[t.final].jpi)
            for t in derived.transitions if round(t.energy.value, 1) in rates})
        out = {"available": True, "residue": residue, "channel": row["channel"], "rows": rows,
               "levels": [{"energy_kev": e, "label": lab} for e, lab in levels],
               "reference": scheme.reference, "from_setup": residue in exp.levels}
        self._cache[key] = out
        return out

    def residue_events(self, events: int = 20_000, seed: int = 1):
        """A sample of the fusion-evaporation events (:func:`physim.nuclear.residues.simulate_residues`): the
        residues, their recoil and the particles they evaporated. Kept until the setup changes."""
        if self.experiment.fusion is None:
            raise ValueError("the setup is not a fusion-evaporation reaction")
        key = ("residue_events", events, seed)
        if key not in self._cache:
            from .residues import simulate_residues

            self._cache[key] = simulate_residues(self.experiment, events, seed)
        return self._cache[key]

    def recoil_cone(self) -> dict:
        """Where the residues go: the half-angle holding 90 % of them, their mean speed and its spread."""
        ev = self.residue_events()
        theta, beta = np.asarray(ev["theta"]), np.asarray(ev["beta"])
        return {"half_angle_deg": float(np.percentile(theta, 90)), "beta": float(np.mean(beta)),
                "beta_spread": float(np.std(beta))}

    def evaporation(self, points: int = 1500) -> dict:
        """What the compound nucleus boils off, from :meth:`residue_events`: each light particle's energy
        spectrum in the compound nucleus's frame, the number of each kind per event, and for the studied channel
        the entry states (E*, J) and the residues' recoil (θ, β), as samples of at most ``points`` events. With a
        neutron wall, the neutrons' flight times to it."""
        from .residues import LIGHT_LABEL

        ev = self.residue_events()
        kinds = list(LIGHT_LABEL)
        part = ev.particles
        energies = {}
        edges = np.linspace(0.0, 30.0, 61)
        for k, name in enumerate(kinds):
            e = np.asarray(part["energy_cm"])[np.asarray(part["kind"]) == k]
            if len(e):
                counts, _ = np.histogram(e, bins=edges)
                energies[LIGHT_LABEL[name]] = {"edges": edges.tolist(), "counts": counts.tolist(),
                                               "mean_mev": float(np.mean(e))}
        n_events = len(ev)
        multiplicity = {}
        for k, name in enumerate(kinds):
            per = np.bincount(np.asarray(part["event"])[np.asarray(part["kind"]) == k], minlength=n_events)
            if per.any():
                multiplicity[LIGHT_LABEL[name]] = {"mean": float(per.mean()),
                                                   "shares": (np.bincount(per) / n_events).tolist()}
        i = self.studied_channel()
        label = ev.channels[i] if i is not None and i < len(ev.channels) else None
        mask = ev.select(label)
        rng = np.random.default_rng(1)
        idx = np.flatnonzero(mask)
        if len(idx) > points:
            idx = rng.choice(idx, points, replace=False)
        out = {"energies": energies, "multiplicity": multiplicity, "channel": label,
               "entry": {"spin": np.asarray(ev["entry_spin"])[idx].tolist(),
                         "excitation_mev": np.asarray(ev["entry_excitation"])[idx].tolist()},
               "recoil": {"theta_deg": np.asarray(ev["theta"])[idx].tolist(),
                          "beta": np.asarray(ev["beta"])[idx].tolist()},
               "cone": self.recoil_cone()}
        wall = self.experiment.neutron_wall
        if wall is not None:
            from .neutrons import time_of_flight

            cell = wall.model().cell
            e_n = np.asarray(part["energy"])[np.asarray(part["kind"]) == kinds.index("n")]
            if len(e_n):
                t = time_of_flight(e_n, cell.flight_path_mm())
                counts, t_edges = np.histogram(t, bins=np.linspace(0.0, max(80.0, float(np.percentile(t, 99))), 81))
                out["tof"] = {"edges": t_edges.tolist(), "counts": counts.tolist(),
                              "gamma_flash_ns": float(cell.gamma_flash()), "path_mm": float(cell.flight_path_mm())}
        return out

    def residue_gamma_spectra(self, events: int = 100_000, seed: int = 1, bin_kev: float = 1.0) -> dict:
        """The fusion-evaporation residues' γ rays (:func:`physim.nuclear.residues.simulate_residue_gammas`): for
        each γ-ray detector (and under "all", all of them together) its measured and Doppler-corrected singles spectrum
        (counts in the run per bin, also per exit channel under "channels"), the strongest lines of the residue table, and the generator's notes. ``{"available": False}`` without a
        fusion-evaporation reaction or γ-ray detectors. Kept until the setup changes."""
        exp = self.experiment
        if exp.fusion is None or not exp.gamma_detectors:
            return {"available": False, "reason": "A fusion-evaporation reaction with γ-ray detectors is needed."}
        key = ("residue_gammas", events, seed, bin_kev)
        if key not in self._cache:
            from .residues import simulate_residue_gammas

            g = simulate_residue_gammas(exp, events, seed)
            lines = [(x["energy_kev"], f"{r['residue']} {x['label']}") for r in self.residues()["rows"]
                     for x in r["lines"]]
            top = 1.15 * max([e for e, _ in lines] + [500.0]) * 1e-3
            edges = np.arange(0.0, top + bin_kev * 1e-3, bin_kev * 1e-3)
            names = list(dict.fromkeys(c.name.rsplit(" ", 1)[0] if c.element is not None else c.name
                                       for c in g.crystals))
            spectra = {}
            labels = list(g.events.channels)
            for name in names + [None]:
                measured, _ = g.spectrum(name, corrected=False, bins=edges)
                corrected, _ = g.spectrum(name, corrected=True, bins=edges)
                per = {label: {"measured": g.spectrum(name, corrected=False, channel=label, bins=edges)[0],
                               "corrected": g.spectrum(name, corrected=True, channel=label, bins=edges)[0]}
                       for label in labels}
                spectra[name or "all"] = {"measured": measured, "corrected": corrected, "edges": edges,
                                          "counts": g.counts(name) if name else float(measured.sum()),
                                          "channels": per}
            self._cache[key] = {"available": True, "spectra": spectra, "gamma_detectors": names, "lines": lines,
                                "notes": g.notes, "n_gammas": len(g), "beta": g.beta_nominal,
                                "channels": labels}
        return self._cache[key]

    def tracks(self, n: int = 30, select: str = "all", weighted: bool = True, seed: int = 1,
               events: int = 400_000) -> dict:
        """A sample of simulated events as tracks for the scene (:func:`physim.nuclear.tracks.sample_tracks`):
        {"tracks": [...], "description": how the sample was chosen, "channels": the channel labels}. With
        Coulomb excitation and γ-ray detectors the γ events are used, else the particle events alone."""
        from .tracks import describe, sample_tracks

        exp = self._data_experiment()
        if exp.excitation is not None and exp.gamma_detectors:
            g = self.gamma_events(events, seed)
            ev = g.events
        else:
            g = None
            ev = self.spectra(events=events, seed=seed)["events"]
        tracks = sample_tracks(exp, g, ev, n=n, select=select, weighted=weighted, seed=seed)
        return {"tracks": tracks, "description": describe(tracks, select, weighted, ev.n_events),
                "channels": list(ev.channels), "gammas": g is not None}

    def alignment(self, offset_mm: Optional[float] = None, events: int = 400_000, seed: int = 1,
                  fit: bool = True) -> dict:
        """What a misplaced target does to the analysis: the corrected peak with the true geometry and with the
        target assumed ``offset_mm`` along the beam from its true place (an overlay with the shift and the
        broadening), the diagnostic plot (centroid against ring, per crystal) with the assumed geometry, and the
        offset the plot gives back when fitted (:mod:`physim.nuclear.alignment`). ``{"available": False}``
        without Coulomb excitation and γ-ray detectors. ``offset_mm`` defaults to the current alignment run's
        assumed offset, or 0."""
        from .alignment import diagnostic, fit_offset, overlay, with_offset

        if offset_mm is None:
            offset_mm = float(self.run.summary.get("assumed_offset_mm", 0.0)) if self._beam_run() else 0.0
        exp = self._data_experiment()
        if exp.excitation is None or not exp.gamma_detectors:
            return {"available": False, "reason": "Coulomb excitation with γ-ray detectors is needed."}
        g = self.gamma_events(events, seed)
        assumed = with_offset(exp, offset_mm) if offset_mm else exp
        out = {"available": True, "offset_mm": offset_mm, "overlay": overlay(g, exp, assumed),
               "diagnostic": diagnostic(g, assumed), "energy_kev": 1e3 * g.energy_mev}
        if fit:
            found = fit_offset(g, assumed)
            # The fit says where the target is, from where the analysis assumed it: the truth is −offset.
            out["fit"] = found
            out["fit"]["truth_mm"] = -offset_mm
        if self._beam_run() is not None:
            o = out["overlay"]
            self._keep("alignment.json", {"offset_mm": offset_mm, "energy_kev": out["energy_kev"],
                                          "overlay": {k: o[k] for k in ("shift_kev", "broadening_kev",
                                                                        "fwhm_true_kev", "fwhm_assumed_kev")},
                                          "fit": out.get("fit")}, self.run.folder)
        return out

    def analysis(self, events: int = 400_000, seed: int = 1):
        """The :class:`~physim.nuclear.analysis.Analysis` of the simulated γ rays (:meth:`gamma_events`), kept
        with its history until the setup changes."""
        from .analysis import Analysis

        key = ("analysis", events, seed) if not self._beam_run() else ("run", "analysis", self.run.label)
        if key not in self._cache:
            self._cache[key] = Analysis(self._data_experiment(), self.gamma_events(events, seed),
                                        run=self._beam_run())
        return self._cache[key]

    def analyse(self, settings=None, events: int = 400_000, seed: int = 1, gate: Optional[str] = None) -> dict:
        """Run the automatic analysis with ``settings`` (a :class:`~physim.nuclear.analysis.Settings`, or a dict
        of its fields) and return the result as a dict for display. ``{"available": False}`` without Coulomb
        excitation and γ-ray detectors."""
        from dataclasses import asdict

        from .analysis import Settings

        exp = self._data_experiment()
        if exp.excitation is None or not exp.gamma_detectors:
            return {"available": False, "reason": "Coulomb excitation with γ-ray detectors is needed."}
        if gate is not None:  # a named gate of the Data tab, with any other settings on top
            extra = settings.as_dict() if isinstance(settings, Settings) else dict(settings or {})
            for k in ("detectors", "rings", "particle_gate", "particle_energy"):
                extra.pop(k, None)
            settings = self.gate(gate).settings(**extra)
        elif isinstance(settings, dict):
            settings = Settings(**settings)
        try:
            r = self.analysis(events, seed).run(settings)
        except ValueError as e:
            return {"available": False, "reason": str(e)}
        if self._beam_run() is not None:  # kept with the run: the reports read it from there
            from .analysis import save_result

            try:
                save_result(self.run.folder, r)
            except OSError:
                pass
            from . import result_rows

            result_rows.write_coulex(self, r)  # a row of the results ledger; never fails the analysis
        out = asdict(r)
        out.update(available=True, total_unc=r.total_unc, runs=len(self.analysis(events, seed).history))
        return out

    def run_analysis(self, run=None):
        """The last analysis of a run (the current one by default), a
        :class:`~physim.nuclear.analysis.Result`: the one made in this session, or the one kept in the run's
        folder (``analysis.json``), which outlives a change of the setup and closing the app. None for a source
        run or a run that was never analysed."""
        from .analysis import load_result

        run = run if run is not None else self.run
        if run is None or run.kind == "source":
            return None
        kept = self._cache.get(("run", "analysis", run.label)) if run is self.run else None
        if kept is not None and kept.history:
            return kept.history[-1]
        return load_result(run.folder)

    def compare(self, run=None) -> dict:
        """Predicted against measured (backlog item 68): the Plan's numbers kept in the run's summary beside what
        the run measured, with the measured value's statistical uncertainty and whether the two agree (within three
        standard deviations). Rates per detector, full-energy-peak coincidences per detector, the beam time for
        the counts wanted, and each γ-ray detector's efficiency at the transition (measured from the experiment's
        latest source run, if there is one). {"rows": [...], "run": number, "source_run": number or None}."""
        run = run or self.run
        if run is None or run.kind == "source":
            return {"rows": [], "run": None, "source_run": None, "notes": []}
        plan = run.summary.get("plan") or {}
        real = run.real_s
        rows = []

        def row(quantity, what, predicted, measured, unc, unit, model=0.0):
            agree = None
            if predicted is not None and measured is not None and unc is not None:
                # ``model``: how well the plan itself knows the number, relative, beside the run's statistics.
                agree = (abs(measured - predicted) <= math.hypot(3 * max(unc, 1e-300), model * abs(predicted))
                         or (measured == predicted))
            rows.append({"quantity": quantity, "what": what, "predicted": predicted, "measured": measured,
                         "uncertainty": unc, "unit": unit, "agree": agree})

        if run.kind == "fusion":
            return self._compare_fusion(run, plan, row, rows)
        counts = run.summary.get("counts", {})
        peak = run.summary.get("peak_coincidences", {})
        for d in plan.get("detectors", []):
            n = counts.get(d["detector"], 0)
            row("Rate", d["detector"], d["rate_per_s"], n / real, math.sqrt(max(n, 1)) / real, "1/s")
        for d in plan.get("detectors", []):
            if d.get("coincidence_per_s") is None:
                continue
            k = sum(v for key, v in peak.items() if key.split("|")[0] == d["detector"])
            row("Coincidences in the peak", d["detector"], d["coincidence_per_s"], k / real,
                math.sqrt(max(k, 1)) / real, "1/s")
        other, other_peak = plan.get("other_peak"), run.summary.get("other_peak_coincidences")
        if other and other_peak is not None:
            for d in plan.get("detectors", []):
                if d.get("other_coincidence_per_s") is None:
                    continue
                k = sum(v for key, v in other_peak.items() if key.split("|")[0] == d["detector"])
                # The plan leaves out the summing of this nucleus's cascade in a crystal and the Doppler shift
                # of the efficiency: a few percent, more than the statistics of so strong a peak.
                row(f"Coincidences in the {other['label']} peak", d["detector"], d["other_coincidence_per_s"],
                    k / real, math.sqrt(max(k, 1)) / real, "1/s", model=0.05)
        limiting = plan.get("limiting_detector")
        if limiting and plan.get("counts_wanted"):
            what = plan.get("measured", "all")
            if what == "coincidences":
                k = sum(v for key, v in peak.items() if key.split("|")[0] == limiting)
            else:
                k = counts.get(limiting, 0)
            rate = k / real
            t = plan["counts_wanted"] / rate if rate > 0 else math.inf
            row("Beam time for the counts wanted", limiting, plan.get("beam_time_needed_s"), t,
                t / math.sqrt(k) if k > 0 else None, "s")
        source = next((r for r in reversed(self.runs()) if r["kind"] == "source"), None) if self.folder else None
        if plan.get("gamma_detectors"):
            sr = self.folder.load_run(source["number"]).source() if source else None
            for g in plan["gamma_detectors"]:
                measured = unc = None
                if sr is not None:
                    measured, unc = _efficiency_at(sr, g["detector"], g["energy_kev"] * 1e-3)
                row("Efficiency at the transition", g["detector"], g["efficiency"], measured, unc, "")
        return {"rows": rows, "run": run.number, "source_run": source["number"] if source else None, "notes": []}

    def _compare_fusion(self, run, plan: dict, row, rows: list) -> dict:
        """:meth:`compare` for a fusion run, which has no particle detectors: the studied channel's residues per
        second, and its strongest line's full-energy counts per second (shown, not judged: the run keeps only the
        events its γ-ray trigger takes, the plan's line rate has no trigger)."""
        real = run.real_s
        notes = ["A fusion run has no particle detectors: the plan's particle rates, coincidences and beam time for "
                 "the counts wanted do not apply to it."]
        study = plan.get("study")
        if not study:
            notes.append("The run's plan has no studied channel (no residue table when it was taken): nothing to "
                         "compare.")
            return {"rows": rows, "run": run.number, "source_run": None, "notes": notes}
        label = study["channel"]
        per = (run.summary.get("per_channel") or {}).get(label)
        if per is not None:
            n = per.get("residues", 0)
            row("Residues made", label, study.get("rate_per_s"), n / real, math.sqrt(max(n, 1)) / real, "1/s")
        line = study.get("line")
        channels = list(run.summary.get("channels") or [])
        if line and label in channels:
            g = run.gammas()
            m = (np.asarray(g["channel"]) == channels.index(label)) & np.asarray(g["full"], dtype=bool)
            m &= np.abs(1e3 * np.asarray(g["energy0"], dtype=float) - line["energy_kev"]) < 0.05
            k = int(m.sum())
            rows.append({"quantity": "Full-energy counts in the line", "what": f"{label} {line['energy_kev']:.1f} keV",
                         "predicted": line.get("detected_per_s"), "measured": k / real,
                         "uncertainty": math.sqrt(max(k, 1)) / real, "unit": "1/s", "agree": None})
            trig = run.summary.get("trigger", {})
            notes.append(f"The line's counts are those of the events the trigger kept (γ-ray fold ≥ "
                         f"{trig.get('fold', '?')}); the plan's line rate has no trigger, so the two are shown side by "
                         "side, not compared.")
        return {"rows": rows, "run": run.number, "source_run": None, "notes": notes}

    def run_yields(self, gate: Optional[str] = None, settings=None) -> dict:
        """The γ-ray yields the current run measured, per particle detector, for the multi-step fit
        (:meth:`fit_matrix_elements`): the analysis of each detector's gated, corrected peak, its area over the
        efficiency and the correlation, as {detector: {(1, 0): (γ rays emitted, uncertainty)}} in the run's real
        part, and the beam time that is (s)."""
        from .analysis import Settings

        run = self._beam_run()
        if run is None:
            raise ValueError("take a beam run first")
        if run.kind == "fusion":
            raise ValueError("the current run is a fusion run: it has no particle detectors, so there are no "
                             "particle-gated γ-ray yields to fit (take a Coulex beam run)")
        a = self.analysis()
        out = {}
        for d in run.events().detectors:
            base = self.gate(gate).settings() if gate else (settings or Settings())
            s = Settings(**{**base.as_dict(), "detectors": [d]})
            try:
                r = a.run(s)
            except ValueError:
                continue
            if not r.area > 0:
                continue
            # The inelastic gate keeps a share of the excitations; the analysis knows it from the simulation.
            gate_eff = getattr(a, "_gate_efficiency", 1.0) if s.particle_gate == "inelastic" else 1.0
            y = r.area / (r.efficiency * r.correlation * r.branch * gate_eff)
            out[d] = {(1, 0): (y, y * math.hypot(r.fit.area_unc / r.area, 0.0))}
        return {"measured": out, "beam_time_s": run.real_s}

    def explanations(self, detector: Optional[str] = None, gamma_detector: Optional[str] = None,
                     result=None) -> list:
        """Every number of the plan (and of the analysis ``result``, if given) with its explanation, as dicts
        (:mod:`physim.nuclear.record`)."""
        from .record import analysis_explanations, explanations

        out = [x.as_dict() for x in explanations(self, detector, gamma_detector)]
        if result is not None:
            out += [x.as_dict() for x in analysis_explanations(self, result)]
        return out

    def record_html(self, result=None, scene_png: Optional[str] = None) -> str:
        """The run record as one HTML page (:func:`physim.nuclear.record.record_html`); ``result`` is the last
        analysis result if None and one has been run."""
        from .record import record_html

        if result is None:
            result = self.run_analysis()
        if result is None and self._beam_run() is None:
            for key, value in list(self._cache.items()):  # a worker thread may be adding to it
                if isinstance(key, tuple) and key[:1] == ("analysis",) and value.history:
                    result = value.history[-1]
        return record_html(self, result, scene_png=scene_png)

    def multistep(self, role: Optional[str] = None, angle_step: float = 4.0, energies: int = 2,
                  shapes: bool = True) -> dict:
        """Multi-step Coulomb excitation of the setup's level scheme (:mod:`physim.nuclear.multistep`): the γ
        yields of every transition in every particle detector with all orders and with first order, the
        populated levels, and (``shapes``) the prolate–zero–oblate comparison of the first 2⁺ state. Uses the
        scheme of ``role`` ("target" or "beam"; the excited nucleus's by default). Kept until the setup changes.
        ``{"available": False}`` without a level scheme."""
        from .multistep import Multistep

        exp = self.experiment
        if role is None:
            wanted = None if exp.excitation is None else ("target" if exp.excitation.excite == "target" else "beam")
            role = wanted if wanted in exp.levels else next(iter(exp.levels), None)
        if role not in exp.levels:
            return {"available": False, "reason": "Look up a level scheme first (Excitation and γ rays)."}
        key = ("multistep", role, angle_step, energies, shapes)
        if key in self._cache:
            return self._cache[key]
        scheme = exp.levels[role]
        ms = Multistep(exp, scheme, role, angle_step=angle_step, energies=energies)
        try:
            full = ms.yields()
            first = Multistep(exp, scheme, role, angle_step=angle_step, energies=energies, first_order=True).yields()
        except ValueError as e:
            return {"available": False, "reason": str(e)}
        labels = {t: f"{scheme.levels[t[0]].label} ({scheme.levels[t[0]].energy.value:g} keV) → "
                     f"{scheme.levels[t[1]].label}" for t in full.transitions}
        out = {"available": True, "role": role, "nuclide": scheme.nuclide, "detectors": list(full.detectors),
               "transitions": [{"transition": t, "label": labels[t], "energy_kev": scheme.levels[t[0]].energy.value
                                - scheme.levels[t[1]].energy.value,
                                "all_orders": {d: full.detectors[d].get(t, 0.0) for d in full.detectors},
                                "first_order": {d: first.detectors[d].get(t, 0.0) for d in first.detectors}}
                               for t in full.transitions],
               "levels": [{"index": n, "label": scheme.levels[n].label, "energy_kev": scheme.levels[n].energy.value,
                           "all_orders": {d: full.levels[d].get(n, 0.0) for d in full.levels},
                           "first_order": {d: first.levels[d].get(n, 0.0) for d in first.levels}}
                          for n in range(1, len(scheme.levels)) if any(full.levels[d].get(n, 0.0) > 0
                                                                       for d in full.levels)],
               "notes": full.notes, "gosia_input": None, "beam_time_s": _q(exp.run.beam_time).to("s")}
        try:
            out["gosia_input"] = ms.gosia_input(gamma="array" if exp.gamma_detectors else "4pi")
        except ValueError as e:
            out["notes"] = list(out["notes"]) + [f"No GOSIA input: {e}."]
        if shapes:
            twos = [n for n, lev in enumerate(scheme.levels) if n and lev.spin == 2 and scheme.b(0, n, "E2")]
            if twos:
                sh = ms.shapes(twos[0])
                out["shapes"] = {"level": twos[0], "q_efm2": sh["q_efm2"], "rings": sh["rings"],
                                 "totals": sh["totals"], "separable": sh["separable"],
                                 "total_difference_sigma": sh["total_difference_sigma"],
                                 "transition": sh["transition"]}
        self._cache[key] = out
        if (role, angle_step, energies, shapes) == (self._default_role(), 4.0, 2, True):
            self._keep("multistep.json", {k: v for k, v in out.items() if k != "gosia_input"}, self.experiment_folder)
        return out

    def _default_role(self) -> Optional[str]:
        exp = self.experiment
        wanted = None if exp.excitation is None else ("target" if exp.excitation.excite == "target" else "beam")
        return wanted if wanted in exp.levels else next(iter(exp.levels), None)

    def _keep(self, name: str, what: dict, folder) -> None:
        """Write a result beside the data it belongs to (the experiment's folder or a run's), with the setup it
        was computed for and the E1 polarisation of its collisions, so that the reports find it later."""
        import json

        from .gamma import e1_polarisation

        if folder is None:
            return
        try:
            (Path(folder) / name).write_text(json.dumps(
                {"setup": self.experiment.to_dict(), "e1_polarisation": e1_polarisation(self.experiment),
                 "result": what}, default=_runs._json_default), encoding="utf-8")
        except (OSError, TypeError, ValueError):
            pass

    def _kept(self, name: str, folder, same_setup: bool = True) -> Optional[dict]:
        """What :meth:`_keep` wrote, or None. With ``same_setup``, only when it was computed for the current setup's
        physics; a result without an E1 polarisation beside it is from before there was one, and has none."""
        import json

        from .gamma import e1_polarisation

        if folder is None or not (Path(folder) / name).is_file():
            return None
        try:
            d = json.loads((Path(folder) / name).read_text(encoding="utf-8"))
            if same_setup and _runs.physics_changes(d["setup"], json.loads(json.dumps(
                    self.experiment.to_dict(), default=_runs._json_default))):
                return None
            if same_setup and d.get("e1_polarisation", 0.0) != e1_polarisation(self.experiment):
                return None
            return d["result"]
        except (OSError, ValueError, KeyError, TypeError):
            return None

    def kept_multistep(self) -> Optional[dict]:
        """The multi-step solution of the current setup (:meth:`multistep` with its defaults) if it has been
        solved, in this session or earlier (kept in the experiment's folder as ``multistep.json``); None if not,
        or if the setup's physics changed since. Never starts the calculation."""
        solved = self._cache.get(("multistep", self._default_role(), 4.0, 2, True))
        if solved is not None:
            return solved
        return self._kept("multistep.json", self.experiment_folder)

    def kept_alignment(self, run=None) -> Optional[dict]:
        """The last alignment check of a run (the current one by default), as kept in its folder
        (``alignment.json``): the assumed offset, the peak's shift and broadening, and the fitted offset."""
        run = run if run is not None else self.run
        if run is None or run.kind == "source":
            return None
        return self._kept("alignment.json", run.folder, same_setup=False)

    def fit_matrix_elements(self, measured: Optional[dict] = None, free: Optional[list] = None,
                            role: Optional[str] = None, beam_time_s: Optional[float] = None,
                            gate: Optional[str] = None) -> dict:
        """Fit up to three matrix elements of the level scheme to measured counts (see
        :meth:`physim.nuclear.multistep.Multistep.fit`). Without ``measured``, to the current run's yields
        (:meth:`run_yields`, in its real part); ``free`` defaults to the ground state's E2 to the first excited
        state."""
        from .multistep import Multistep

        exp = self._data_experiment()
        if role is None:
            role = "target" if exp.excitation is None or exp.excitation.excite == "target" else "beam"
        of_run = measured is None
        if measured is None:
            y = self.run_yields(gate)
            measured, beam_time_s = y["measured"], y["beam_time_s"]
            if not measured:
                raise ValueError("the run has no peak to fit in any particle detector")
        free = free or [(0, 1, "E2")]
        fit = Multistep(exp, exp.levels[role], role).fit(measured, free, beam_time_s=beam_time_s)
        if of_run:  # the run's own yields: rows of the results ledger; counts given from outside have no run
            from . import result_rows

            result_rows.write_fit(self, fit, role, gate)
        return fit

    def trajectories(self, impact_parameters: Optional[list] = None, nuclide: Optional[str] = None) -> dict:
        """Coulomb orbits (CM frame, fm) for a range of impact parameters, from physim's engine."""
        exp = self.experiment
        ion = beam_ion(exp)
        target = nuclide or data.nuclide(max(stack(exp)[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        r = Rutherford(ion, target, exp.beam.energy_mev)
        bs = np.asarray(impact_parameters if impact_parameters is not None
                        else r.d0 * np.array([0.1, 0.25, 0.5, 1.0, 2.0, 4.0]), dtype=float)
        orbits = r.trajectories(bs, distance=40 * r.d0)
        return {"target": target, "d0_fm": r.d0, "interaction_radius_fm": r.interaction_radius,
                "orbits": [{"b_fm": o.b, "xy": o.pos, "deflection_deg": o.deflection(),
                            "closest_fm": o.closest_approach()} for o in orbits],
                "grazing_angle_deg": r.grazing_angle()}

    def gamma(self) -> dict:
        """Coulomb excitation: the excited state, its cross section, the excitation probability against angle, and
        the Doppler-shifted γ-ray energy for every particle × γ detector pair. ``{"available": False}`` without an
        excited state."""
        exc = self.experiment.excitation
        if exc is None:
            return {"available": False, "reason": "The setup has no excited state ([reaction] type = \"coulex\")."}
        from .gamma import doppler_table
        from .rates import coulex_for

        r = self._rates()
        ch = next(c for c in r.channels if c.excitation is not None)
        cx = coulex_for(self.experiment, ch, r.layers)
        th = np.linspace(1.0, 180.0, 180)
        rates = {g.name: sum(v for (label, _), v in r.by_channel(g.name).items() if label == ch.label)
                 for g in r.array}
        if "gamma" not in self._cache:
            self._cache["gamma"] = doppler_table(self.experiment) if self.experiment.gamma_detectors else []
        return {"available": True, "state": {"excite": exc.excite, "energy_kev": exc.energy_mev * 1e3,
                                             "multipolarity": exc.multipolarity, "b_up_e2fm": exc.b_up_e2fm,
                                             "e1_polarisation": cx.e1_polarisation},
                "xi": cx.xi, "eta": cx.eta, "safe_distance_fm": cx.safe_distance, "max_safe_angle": cx.max_safe_angle(),
                "total_mb": cx.total(), "theta_cm": th, "probability": cx.probability(th),
                "rates": rates, "particles": self._particle_energies(ch), "doppler": self._cache["gamma"],
                "emission": exc.emission or "correlated", "correlation": self.correlation()}

    def correlation(self) -> list:
        """For every particle detector × γ-ray crystal, the γ rays seen in coincidence relative to isotropic
        emission (:func:`physim.nuclear.gamma.correlation_table`). Empty without γ-ray detectors."""
        if "correlation" not in self._cache:
            from .gamma import correlation_table

            exp = self.experiment
            self._cache["correlation"] = (correlation_table(exp) if exp.excitation is not None
                                          and exp.gamma_detectors else [])
        return self._cache["correlation"]

    def populations(self, role: Optional[str] = None) -> dict:
        """First-order Coulomb excitation of a whole level scheme of the setup (``role`` "target" or "beam"; by
        default the nucleus the reaction excites, if it has a scheme, else whichever has one): for each level the
        cross section for exciting it directly and with feeding from above, and for each transition the cross
        section of its γ ray, integrated over all scattering angles at the mid-target energy.
        ``{"available": False}`` without a level scheme."""
        from .gamma import e1_polarisation
        from .orientation import Excitation

        exp = self.experiment
        if role is None:
            wanted = None if exp.excitation is None else ("target" if exp.excitation.excite == "target" else "beam")
            role = wanted if wanted in exp.levels else next(iter(exp.levels), None)
        if role not in exp.levels:
            return {"available": False, "reason": "Look up a level scheme first."}
        key = ("populations", role)
        if key not in self._cache:
            scheme = exp.levels[role]
            layers = stack(exp)
            e_mid = round(float(beam_energy_at(exp, 0, [layers[0].thickness / 2], layers)[0]), 9)
            target = data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
            if role == "target":
                target = scheme.nuclide
            try:
                ex = Excitation(beam_ion(exp), target, e_mid, scheme,
                                excite="target" if role == "target" else "projectile",
                                e1_polarisation=e1_polarisation(exp))
            except ValueError as e:
                return {"available": False, "reason": str(e)}
            cs = ex.cross_sections()
            levels = [{"index": n, "label": lev.label, "energy_kev": lev.energy.value,
                       "direct_mb": cs["direct"].get(n, 0.0), "populated_mb": cs["populated"].get(n, 0.0)}
                      for n, lev in enumerate(scheme.levels) if n and cs["populated"].get(n, 0.0) > 0]
            gammas = [{"initial": i, "final": f, "energy_kev": scheme.levels[i].energy.value
                       - scheme.levels[f].energy.value, "label": f"{scheme.levels[i].label} → "
                       f"{scheme.levels[f].label}", "gamma_mb": v}
                      for (i, f), v in cs["gamma"].items() if v > 0]
            self._cache[key] = {"available": True, "role": role, "nuclide": scheme.nuclide, "levels": levels,
                                "gammas": sorted(gammas, key=lambda g: -g["gamma_mb"]), "notes": list(ex.notes),
                                "beam_energy_mev": e_mid}
        return self._cache[key]

    # -- level schemes ----------------------------------------------------------------------------------------------

    def level_nuclides(self, draft: Optional[dict] = None) -> dict:
        """The nuclide whose level scheme each role takes: the beam, and the target's most abundant nuclide (all of
        them are listed under "target_choices"); of the setup, or of ``draft`` (a setup being edited)."""
        d = draft if draft is not None else self._draft
        out = {"beam": None, "target": None, "target_choices": []}
        try:
            z, a = parse_nuclide(d["beam"]["nuclide"])
            out["beam"] = _data.nuclide((z, a)).name
        except (KeyError, ValueError):
            pass
        try:
            atoms = sorted(_data.material(d["target"]["material"]).atoms, key=lambda x: -x[2])
            out["target_choices"] = [_data.nuclide((z, a)).name for z, a, _ in atoms]
            out["target"] = out["target_choices"][0]
        except (KeyError, ValueError, IndexError):
            pass
        return out

    def levels(self) -> dict:
        """The level schemes in the setup, as tables for display (:meth:`LevelScheme.table`), with whether a local
        ENSDF copy is present and which nuclide each role would look up."""
        schemes = self.experiment.levels if not self.problems else {}
        return {"ensdf": _ensdf.available(), "ensdf_folder": str(_ensdf.folder()), "nuclides": self.level_nuclides(),
                "beam": schemes["beam"].table() if "beam" in schemes else None,
                "target": schemes["target"].table() if "target" in schemes else None}

    def lookup_levels(self, role: str, max_energy: str = "3 MeV", nuclide: Optional[str] = None) -> bool:
        """Read the level scheme of the beam or the target nuclide from the local ENSDF copy into the setup, up to
        ``max_energy``. Matrix elements the user had set for the same nuclide are kept. Raises
        :class:`~physim.nuclear.ensdf.EnsdfMissing` without a copy or without data for the nuclide."""
        if role not in ("beam", "target"):
            raise ValueError("role must be 'beam' or 'target'")
        nuclide = nuclide or self.level_nuclides()[role]
        if nuclide is None:
            raise ValueError(f"the setup has no {role} nuclide yet")
        scheme = LevelScheme.from_ensdf(nuclide, max_energy_kev=_q(max_energy).to("keV"))
        d = self.draft
        old = d.get("levels", {}).get(role)
        if old and old.get("nuclide") == scheme.nuclide:
            energies = [lv.energy.value for lv in scheme.levels]
            old_energy = [_q(lv["energy"]).to("keV") for lv in old.get("level", [])]
            for m in old.get("matrix_element", []):
                if m.get("source") not in ("user", "shell model"):  # values someone chose to put in
                    continue
                try:  # the same pair of levels, found by energy
                    a, b = (energies.index(old_energy[m[k]]) for k in ("from", "to"))
                except (ValueError, IndexError, KeyError):
                    continue
                scheme.set_matrix_element(a, b, m.get("multipolarity", "E2"), m["value"], m.get("value_unc"),
                                          source=m["source"], note=m.get("note", ""))
        d.setdefault("levels", {})[role] = scheme.to_dict()
        return self._apply(d)

    def apply_shell_model(self, role: str, result, which: Any = "all", replace: bool = False) -> list:
        """Put chosen values of a shell-model calculation (a :class:`~physim.nuclear.shellmodel.ShellModelResult`)
        into the level scheme of the beam or the target, with the source ``"shell model"``, through
        :func:`physim.nuclear.shellmodel.apply` (``which``: ``"all"``, or quantities and states such as
        ``["E2", "2+1"]``). ENSDF, derived and user values are kept unless ``replace`` is true, which only the user's
        explicit choice should set. Without a scheme for the role in the setup, one is started from the result: its
        states become the levels. Returns what changed, one line each."""
        from . import shellmodel as _shellmodel

        if role not in ("beam", "target"):
            raise ValueError("role must be 'beam' or 'target'")
        nuclide = self.level_nuclides()[role]
        if nuclide is None:
            raise ValueError(f"the setup has no {role} nuclide yet")
        if result.nuclide and parse_nuclide(result.nuclide) != parse_nuclide(nuclide):
            raise ValueError(f"the shell-model result is for {result.nuclide}, but the {role} nuclide is {nuclide}")
        scheme = self.experiment.levels.get(role)
        if scheme is None or not scheme.levels:
            if not result.states:
                raise ValueError("the shell-model result has no states")
            # An empty scheme takes the calculated states as its levels, starting from the lowest as the ground state.
            from .levels import Level, Value

            lowest = min(result.states, key=lambda s: s.energy_mev)
            scheme = LevelScheme(nuclide, reference=f"shell model ({result.interaction})")
            scheme.levels = [Level(Value(0.0, None, "shell model", f"{result.interaction}: {lowest.label}"),
                                   lowest.j2 / 2, lowest.parity, _shellmodel.label(lowest.j2, lowest.parity))]
            items = [which] if isinstance(which, str) else list(which)
            which = [*items, "new"]
        changes = _shellmodel.apply(result, scheme, which, replace=replace)
        d = self.draft
        d.setdefault("levels", {})[role] = scheme.to_dict()
        self._apply(d)
        from . import result_rows

        result_rows.copy_shell_model(result)  # the whole calculation, as theory rows of the nucleus
        return changes

    def remove_levels(self, role: str) -> bool:
        """Take a level scheme out of the setup."""
        d = self.draft
        d.get("levels", {}).pop(role, None)
        if not d.get("levels"):
            d.pop("levels", None)
        return self._apply(d)

    def set_matrix_element(self, role: str, a: int, b: int, multipolarity: str, value: str) -> bool:
        """Set a matrix element of a level scheme (text with a unit, e.g. ``"1.28 eb"``); it is marked as the
        user's."""
        scheme = self.experiment.levels[role]
        scheme.set_matrix_element(a, b, multipolarity, value)
        d = self.draft
        d["levels"][role] = scheme.to_dict()
        return self._apply(d)

    def use_state(self, role: str, level: int, multipolarity: str = "E2") -> bool:
        """Plan Coulomb excitation of one state of a level scheme: sets the reaction's excited nucleus, energy,
        multipolarity and B(Eλ↑) from the scheme."""
        scheme = self.experiment.levels[role]
        b_up = scheme.b(0, level, multipolarity)
        if not b_up:
            raise ValueError(f"the scheme has no {multipolarity} matrix element between the ground state and "
                             f"level {level}")
        lam = int(multipolarity[1])
        d = self.draft
        d.setdefault("reaction", {}).update(
            {"type": "coulex", "excite": "target" if role == "target" else "projectile",
             "energy": f"{scheme.levels[level].energy.value:.10g} keV", "multipolarity": multipolarity,
             "b_up": f"{b_up:.6g} e2fm{2 * lam}"})
        return self._apply(d)

    def _particle_energies(self, channel) -> list:
        """What each particle detector measures in a Coulomb-excitation run: the energy of the beam particle or
        recoil at the detector's smallest, central and largest lab angle, after elastic scattering and after exciting
        the state, and the speed β of the excited nucleus in that event (for the Doppler correction). Energies are at
        the reaction point, before losses in the target."""
        beam = self.experiment.beam
        ion = beam_ion(self.experiment)
        elastic = TwoBody(ion, channel.nuclide.name, beam.energy_mev)
        try:
            inelastic = TwoBody(ion, channel.nuclide.name, beam.energy_mev, **channel.kinematics_args)
        except ValueError:  # below the excitation threshold
            return []
        excited = channel.kinematics_args["excite"]  # "recoil" or "ejectile"
        m_exc = inelastic.m4 if excited == "recoil" else inelastic.m3

        def beta(t, m):
            return float(np.sqrt(t * (t + 2 * m)) / (t + m))

        rows = []
        for g in Array.from_experiment(self.experiment):
            lo, hi = g.theta_range()
            for particle, name in (("ejectile", ion), ("recoil", channel.nuclide.name)):
                reach = min(elastic.max_angle(particle), inelastic.max_angle(particle))
                for where, th in (("min", lo), ("centre", 0.5 * (lo + hi)), ("max", hi)):
                    if not 0 < th <= reach:
                        continue
                    el = elastic.at_lab(th, particle)[0]
                    ex = inelastic.at_lab(th, particle)[0]
                    if not (np.isfinite(el.energy) and np.isfinite(ex.energy)):
                        continue
                    # The excited nucleus is the detected particle itself or its partner in the same event.
                    if particle == excited:
                        t_exc = float(ex.energy)
                    else:
                        partner = "recoil" if particle == "ejectile" else "ejectile"
                        t_exc = float(inelastic.at_cm(180.0 - float(ex.theta_cm), partner).energy)
                    rows.append({"detector": g.name, "particle": particle, "nuclide": name, "where": where,
                                 "theta_lab": float(th), "elastic_mev": float(el.energy),
                                 "excited_mev": float(ex.energy),
                                 "difference_mev": float(el.energy) - float(ex.energy),
                                 "beta_excited": beta(t_exc, m_exc)})
        return rows

    def report(self) -> dict:
        """What the beam-time report will contain (the export itself is item 42)."""
        return {"title": self.experiment.title, "setup": self.experiment.to_dict(), "rates": self.rates()["rows"],
                "peaks": {g.name: self._rates().peaks(g.name) for g in self._rates().array},
                "warnings": self.warnings()}

    # -- sweeps -----------------------------------------------------------------------------------------------------

    def sweep(self, parameter: str, values: list, quantity: str = "rate", detector: Optional[str] = None) -> dict:
        """Vary one parameter and recompute one quantity for a detector.

        ``parameter``: "beam energy", "target thickness", or "detector angle" (of ``detector``). ``quantity``:
        "rate" (counts/s), "beam time" (s for the counts wanted), "peak energy" (MeV, the strongest peak) or
        "peak width" (FWHM, MeV). Values are given with units, as in a setup file.
        """
        detector = detector or self.experiment.detectors[0].name
        out = []
        for v in values:
            d = self.draft
            # The rates do not depend on a lifetime measurement, and a swept value may make the plunger impossible
            # (a recoil too slow to leave the target) without changing what is swept.
            d.pop("lifetime", None)
            if parameter == "beam energy":
                d["beam"]["energy"] = v
            elif parameter == "target thickness":
                d["target"]["thickness"] = v
            elif parameter == "detector angle":
                det = d["detectors"][self._index(detector)]
                det["theta"] = v
                det.pop("position", None)
            else:
                raise ValueError("parameter must be 'beam energy', 'target thickness' or 'detector angle'")
            r = Rates(Experiment.from_dict(d))
            if quantity == "rate":
                y = r.rate(detector)
            elif quantity == "beam time":
                y = r.beam_time_for(detector, what="measured")
            elif quantity in ("peak energy", "peak width"):
                pk = r.peaks(detector)
                y = (pk[0].mean if quantity == "peak energy" else pk[0].fwhm) if pk else math.nan
            else:
                raise ValueError("quantity must be 'rate', 'beam time', 'peak energy' or 'peak width'")
            out.append(y)
        return {"parameter": parameter, "values": list(values), "x": [_q(v).value for v in values],
                "quantity": quantity, "detector": detector, "y": np.array(out, dtype=float)}


#: The relative uncertainty allowed for reading an efficiency between (or beyond) a source's lines.
INTERPOLATION_UNC = 0.03


def _efficiency_at(source_run, name: str, energy_mev: float) -> tuple:
    """A γ-ray detector's (or crystal's) full-energy efficiency at ``energy_mev`` from a source run's efficiency
    points, interpolated in log–log between the lines: (efficiency, uncertainty), or (None, None). The uncertainty
    is the nearest line's, with 3 % added for the interpolation (or the extrapolation beyond the last line)."""
    names = [n for n in source_run.names() if n == name or n.rsplit(" ", 1)[0] == name]
    if not names:
        return None, None
    pts = {}
    for n in names:  # a detector of several crystals: their efficiencies add
        for p in source_run.efficiency_points(n):
            e, u = pts.get(p["energy_mev"], (0.0, 0.0))
            pts[p["energy_mev"]] = (e + p["efficiency"], math.hypot(u, p["uncertainty"]))
    good = sorted((e, v) for e, v in pts.items() if v[0] > 0)
    if len(good) < 2:
        return None, None
    x = np.log([e for e, _ in good])
    y = np.log([v[0] for _, v in good])
    eff = float(np.exp(np.interp(math.log(energy_mev), x, y)))
    k = int(np.argmin([abs(e - energy_mev) for e, _ in good]))
    return eff, eff * math.hypot(good[k][1][1] / good[k][1][0], INTERPOLATION_UNC)


__all__ = ["EXPLAIN", "REGISTER", "TABS", "Planner", "Warning", "first_state"]
