"""Experiments and runs: the lab's model under the workbench (backlog item 63).

An **experiment** is the unit of work: a name, the current setup, and the runs taken with it. On disk it is one
folder::

    58Ni-on-208Pb/
      experiment.toml             the current setup (a setup file)
      runs/
        01-152Eu-source/   setup.toml  events.npz  summary.json
        02-beam-10min/     ...

A **run** is a beam run of a duration, a source run (beam off, a calibration source at a position) or an alignment
run (a beam run analysed with a deliberately wrong assumed geometry). ``setup.toml`` is the exact setup the run was
taken with, ``summary.json`` its counters, and ``events.npz`` its data as compact named columns readable from any
Python::

    from physim.nuclear import Experiment
    from physim.nuclear.runs import ExperimentFolder

    lab = ExperimentFolder.create(Experiment.example("coulex_ni58"), root="~/.physim/experiments")
    run = lab.take_run("beam", duration="10 min")
    run.counts("CD")                  # particles counted in the real part of the run
    run.events()                      # the particle events, one row per counted particle, each counting once
    run.gammas()                      # the γ rays, with the Doppler corrections derived from the stored columns

**Real statistics.** A beam run generates *unweighted* events — one event is one event — for the first part of
the beam time, up to a budget (10 min by default, at most an hour), and scales the rest: histograms beyond the
real part are the real part scaled to the full duration and redrawn with the full run's Poisson fluctuations
(:meth:`RunData.scaled`). The generator's weighted events are turned into unweighted ones by keeping each with a
probability proportional to its weight (hit-or-miss, :class:`BeamRunner`); the γ chain then emits one γ ray per
excited event. In the stored events each row counts once: its ``weight`` is 1 / (real seconds), so the existing
rates, spectra and analysis code reads them unchanged.

**What is stored.** Per counted particle: event, detector, ring/strip i, sector/strip j, channel (16-bit integers),
whether it is the recoil (8 bits), measured energy, θ_cm and depth (32-bit floats). Per γ ray: the particle row,
the crystal, the deposited and measured energy, and the lab θ and φ. Everything else is derived on loading: the
particle's direction is the centre of the segment it hit (what the detector knows), the Doppler corrections come
from :func:`~physim.nuclear.gamma_events.doppler_correct`, the singles and randoms from the rates. Particles below
threshold are not recorded, as an acquisition would not record them.
"""

from __future__ import annotations

import copy
import json
import math
import os
import re
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional, Union

import numpy as np

from .._core import nuclear_events
from .detectors import Array
from .events import Events, generator_config
from .experiment import Experiment
from .quantity import Quantity
from .rates import Rates

#: The kinds of run.
KINDS = ("beam", "source", "alignment")
#: The real part of a beam run by default, and the most that is offered, s.
DEFAULT_BUDGET_S = 600.0
MAX_BUDGET_S = 3600.0
#: The most reactions generated at once (a chunk); the real part is built chunk by chunk.
MAX_CHUNK = 8_000_000
#: The pilot run that sets the channel mixture and the largest weight of the unweighting.
PILOT_EVENTS = 400_000
#: Setup fields that are labels: changing them does not make a run stale.
LABELS = ("title", "description", "name")

#: The stored particle columns: {key in events.npz: (column of Events, type)}. ``p_event_step`` is the step from
#: the previous row's event number (rows are in event order; the numbers are its running sum), and
#: ``p_measured_kev`` the measured energy in whole keV, far finer than a silicon detector's resolution: both
#: compress to a fraction of a byte or two per row.
PARTICLE_COLUMNS = {
    "p_event_step": ("event", np.int32), "p_detector": ("detector", np.int16),
    "p_segment_i": ("segment_i", np.int16), "p_segment_j": ("segment_j", np.int16),
    "p_channel": ("channel", np.int16), "p_recoil": ("recoil", np.int8),
    "p_measured_kev": ("measured", np.uint32), "p_theta_cm": ("theta_cm", np.float16),
    "p_depth": ("depth", np.float16),
}


def _store_particles(cols: dict, last_event: int) -> dict:
    """Particle columns as stored (see :data:`PARTICLE_COLUMNS`)."""
    out = {}
    for key, (name, t) in PARTICLE_COLUMNS.items():
        x = np.asarray(cols[name])
        if key == "p_event_step":
            x = np.diff(x.astype(np.int64), prepend=last_event)
        elif key == "p_measured_kev":
            x = np.round(np.maximum(x, 0.0) * 1e3)
        out[key] = x.astype(t)
    return out


def _load_particles(a: dict) -> dict:
    cols = {}
    for key, (name, t) in PARTICLE_COLUMNS.items():
        x = a[key]
        if key == "p_event_step":
            x = np.cumsum(x.astype(np.int64))
        elif key == "p_measured_kev":
            x = x.astype(float) * 1e-3
        cols[name] = x.astype(np.int64) if np.dtype(t).kind in "iu" and key != "p_measured_kev" else x.astype(float)
    return cols
GAMMA_COLUMNS = {
    "particle": np.int32, "crystal": np.int16, "deposited": np.float32, "measured": np.float32,
    "theta": np.float32, "phi": np.float32,
}


class RunInProgress(RuntimeError):
    """The setup cannot change while a run is being taken."""


def _seconds(x) -> float:
    if x is None:
        return math.nan
    if isinstance(x, (int, float)):
        return float(x)
    return Quantity.parse(x).to("s") if isinstance(x, str) else x.to("s")


def duration_label(seconds: float) -> str:
    """"10min", "2h", "90s": a short label for a run's folder name."""
    if seconds >= 3600 and seconds % 3600 == 0:
        return f"{int(seconds // 3600)}h"
    if seconds >= 60 and seconds % 60 == 0:
        return f"{int(seconds // 60)}min"
    return f"{seconds:g}s"


def slug(text: str) -> str:
    """A folder name from an experiment's name: "⁵⁸Ni on ²⁰⁸Pb" → "58Ni-on-208Pb"."""
    sup = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")
    out = re.sub(r"[^A-Za-z0-9.+_-]+", "-", text.translate(sup)).strip("-")
    return out or "experiment"


def default_name(experiment: Experiment) -> str:
    """The name an experiment gets from its beam and target, e.g. "58Ni on 208Pb"."""
    from . import data
    from .rates import beam_ion, stack

    try:
        target = data.nuclide(max(stack(experiment)[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        return f"{beam_ion(experiment)} on {target}"
    except (KeyError, ValueError, IndexError):
        return experiment.title or "Experiment"


def home() -> Path:
    """Where experiments live: ``$PHYSIM_EXPERIMENTS``, or ``~/.physim/experiments``."""
    return Path(os.environ.get("PHYSIM_EXPERIMENTS") or Path.home() / ".physim" / "experiments").expanduser()


def list_experiments(root: Union[str, Path, None] = None) -> list:
    """The experiments under ``root`` (default :func:`home`), newest first: dicts with name, path, number of runs
    and when the experiment last changed."""
    root = Path(root).expanduser() if root is not None else home()
    out = []
    if not root.is_dir():
        return out
    for p in root.iterdir():
        if not (p / "experiment.toml").is_file():
            continue
        f = ExperimentFolder(p)
        runs = f.runs()
        stamp = max([(p / "experiment.toml").stat().st_mtime] + [r["_mtime"] for r in runs])
        out.append({"name": f.name, "path": str(p), "runs": len(runs), "modified": stamp})
    return sorted(out, key=lambda x: -x["modified"])


def _strip_labels(obj):
    if isinstance(obj, dict):
        return {k: _strip_labels(v) for k, v in obj.items() if k not in LABELS}
    if isinstance(obj, list):
        return [_strip_labels(v) for v in obj]
    return obj


def physics_changes(taken: dict, current: dict) -> list:
    """What differs in the physics between the setup a run was taken with and the current one, as
    ["target.thickness: 0.5 mg/cm2 → 1.0 mg/cm2", ...]. Labels (titles, descriptions, names) are ignored."""
    out: list = []

    def walk(a, b, path):
        if isinstance(a, dict) and isinstance(b, dict):
            for k in list(a) + [k for k in b if k not in a]:
                walk(a.get(k), b.get(k), f"{path}.{k}" if path else k)
        elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
            for i, (x, y) in enumerate(zip(a, b)):
                walk(x, y, f"{path}[{i + 1}]")
        elif a != b:
            def show(x):
                if x is None:
                    return "(none)"
                if isinstance(x, list):
                    return f"{len(x)} entries"
                return str(x)
            out.append(f"{path}: {show(a)} → {show(b)}")

    walk(_strip_labels(taken), _strip_labels(current), "")
    return out


# -- generating unweighted events ----------------------------------------------------------------------------------

class BeamRunner:
    """Unweighted particle events from the weighted generator, chunk by chunk.

    A pilot run (:data:`PILOT_EVENTS`, its own seed) sets how often each reaction channel is picked, so that the
    events of every channel carry similar weights, and the largest weight q (its 99.9 % quantile). A chunk of N
    reactions then covers t = N / (q × N_pilot) seconds of beam: each event, of weight w (per second), occurs
    w × t times on average, and is kept that many times — the whole part of w × t, and once more with
    probability its fraction. The kept events are exactly as frequent as in a beam of that duration; one of
    weight above q may appear twice. The reactions continue one seed stream (``first`` counts on), so a run that
    is extended carries on where it stopped."""

    def __init__(self, experiment: Experiment, seed: int, state: Optional[dict] = None):
        self.experiment = experiment
        self.seed = int(seed)
        self.config, used = generator_config(experiment)
        self.config["skip_misses"] = True
        self.channels = [ch.label for ch in used]
        self.detectors = [g.name for g in Array.from_experiment(experiment)]
        if state:
            for c, p in zip(self.config["channels"], state["probabilities"]):
                c["probability"] = p
            self.q, self.cursor, self.chunk, self.real_s = state["q"], state["cursor"], state["chunk"], state["real_s"]
            self.generated = state["generated"]
        else:
            self._pilot()
            self.cursor, self.chunk, self.real_s, self.generated = 0, 0, 0.0, 0

    def state(self) -> dict:
        return {"probabilities": [c["probability"] for c in self.config["channels"]], "q": self.q,
                "cursor": self.cursor, "chunk": self.chunk, "real_s": self.real_s, "generated": self.generated}

    def _pilot(self) -> None:
        cols = nuclear_events(self.config, PILOT_EVENTS, self.seed + 0x9E3779B9, PILOT_EVENTS, 0)
        event, first = np.unique(np.asarray(cols["event"]), return_index=True)
        w = np.asarray(cols["weight"])[first]
        ch = np.asarray(cols["channel"])[first]
        chans = self.config["channels"]
        old = np.array([c["probability"] for c in chans])
        mean = np.array([w[ch == k].mean() if np.any(ch == k) else 0.0 for k in range(len(chans))])
        floor = mean[mean > 0].min() * 1e-3 if np.any(mean > 0) else 1.0
        new = old * np.where(mean > 0, mean, floor)
        new = new / new.sum()
        for c, p in zip(chans, new):
            c["probability"] = float(p)
        # The pilot's weights as they would be with the new mixture.
        scale = (old / old.sum()) / new
        w = w * scale[ch]
        self.q = float(np.quantile(w, 0.999)) if len(w) else 1.0

    def seconds_per_reaction(self) -> float:
        return 1.0 / (self.q * PILOT_EVENTS)

    def next_chunk(self, seconds_left: float) -> tuple:
        """The next chunk of at most ``seconds_left`` of beam: (particle columns of the kept events, counted
        particles only, with event numbers that go on from the previous chunks; seconds covered)."""
        per = self.seconds_per_reaction()
        n = int(min(MAX_CHUNK, math.ceil(seconds_left / per)))
        t = min(n * per, seconds_left)
        cols = {k: np.asarray(v) for k, v in nuclear_events(self.config, n, self.seed, n, self.cursor).items()}
        self.cursor += n
        self.generated += n
        rng = np.random.default_rng([self.seed, self.chunk])
        self.chunk += 1
        if not len(cols["event"]):
            self.real_s += t
            return {}, t
        _, first, inverse = np.unique(cols["event"], return_index=True, return_inverse=True)
        p = cols["weight"][first] * t
        times = np.floor(p).astype(np.int64) + (rng.uniform(0.0, 1.0, len(p)) < p - np.floor(p))
        start = np.cumsum(times) - times
        per_row = times[inverse]
        rows = np.repeat(np.arange(len(per_row)), per_row)
        copy_k = np.arange(len(rows)) - np.repeat(np.cumsum(per_row) - per_row, per_row)
        number = start[inverse[rows]] + copy_k
        order = np.argsort(number, kind="stable")
        rows, number = rows[order], number[order]
        keep = (cols["detector"][rows] < len(self.detectors)) & cols["counted"][rows]
        rows, number = rows[keep], number[keep]
        out = {k: v[rows] for k, v in cols.items()}
        base = getattr(self, "_next_event", 0)
        out["event"] = number + base
        self._next_event = base + int(times.sum())
        self.real_s += t
        return out, t


# -- runs on disk -------------------------------------------------------------------------------------------------

def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class RunData:
    """One run, as stored: its setup snapshot, its summary (counters), and its events, loaded when asked for."""

    folder: Path
    summary: dict
    setup: dict
    _events: Optional[Events] = field(default=None, repr=False)
    _gammas: dict = field(default_factory=dict, repr=False)

    @classmethod
    def load(cls, folder: Union[str, Path]) -> RunData:
        folder = Path(folder)
        summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
        setup = Experiment.load(folder / "setup.toml").to_dict()
        return cls(folder, summary, setup)

    # -- what the run is ------------------------------------------------------------------------------------------

    @property
    def number(self) -> int:
        return int(self.summary["number"])

    @property
    def kind(self) -> str:
        return self.summary["kind"]

    @property
    def label(self) -> str:
        return self.folder.name

    @property
    def experiment(self) -> Experiment:
        """The setup the run was taken with."""
        return Experiment.from_dict(copy.deepcopy(self.setup))

    @property
    def duration_s(self) -> float:
        return float(self.summary["duration_s"])

    @property
    def real_s(self) -> float:
        """The part of the run simulated event by event, s."""
        return float(self.summary["real_s"])

    @property
    def scale(self) -> float:
        """Duration over the real part: what histograms of the real part are scaled by for the whole run."""
        return self.duration_s / self.real_s if self.real_s > 0 else 1.0

    def describe(self) -> str:
        """One line: kind, duration, what was real."""
        s = self.summary
        what = {"beam": "Beam run", "source": f"{s.get('source', {}).get('nuclide', 'Source')} source run",
                "alignment": "Alignment run"}[self.kind]
        line = f"{what}, {duration_label(self.duration_s)}"
        if self.kind != "source":
            line += (f"; the first {duration_label(self.real_s)} simulated event by event, the rest scaled"
                     if self.scale > 1.000001 else "; every event simulated")
        if s.get("status") == "stopped":
            line += " (stopped)"
        return line

    def changes(self, current: dict) -> list:
        """The physics changes from this run's setup to ``current`` (a setup dict); empty if the run is current."""
        return physics_changes(self.setup, current)

    def counts(self, detector: Optional[str] = None) -> float:
        """Particles counted in the real part (all detectors, or one)."""
        c = self.summary.get("counts", {})
        return float(sum(c.values()) if detector is None else c.get(detector, 0.0))

    def counts_in_run(self, detector: Optional[str] = None) -> float:
        """Counts over the whole duration: the real part's, scaled."""
        return self.counts(detector) * self.scale

    def scaled(self, counts: np.ndarray, seed: Optional[int] = None) -> np.ndarray:
        """A histogram of the real part as it would be for the whole run: scaled to the duration and redrawn with
        the full run's Poisson fluctuations. Unchanged when the whole run is real."""
        counts = np.asarray(counts, dtype=float)
        if self.scale <= 1.000001:
            return counts
        rng = np.random.default_rng(self.summary["seed"] + 104_729 if seed is None else seed)
        return rng.poisson(np.maximum(counts, 0.0) * self.scale).astype(float)

    # -- the data ---------------------------------------------------------------------------------------------------

    def _arrays(self) -> dict:
        with np.load(self.folder / "events.npz") as z:
            return {k: z[k] for k in z.files}

    def events(self) -> Optional[Events]:
        """The particle events of the real part (see :class:`~physim.nuclear.events.Events`): one row per counted
        particle, each standing for one count (weight 1 / real seconds). None for a source run."""
        if self.kind == "source":
            return None
        if self._events is None:
            a = self._arrays()
            cols = _load_particles(a)
            cols["recoil"] = cols["recoil"].astype(bool)
            cols["counted"] = np.ones(len(cols["event"]), dtype=bool)
            cols["weight"] = np.full(len(cols["event"]), 1.0 / self.real_s if self.real_s > 0 else 0.0)
            theta, phi = _segment_directions(self.experiment, cols)
            cols["theta"], cols["phi"] = theta, phi
            self._events = Events(cols, int(self.summary.get("generated", 0)), int(self.summary["seed"]),
                                  list(self.summary["detectors"]), list(self.summary["channels"]), self.real_s, 0.0)
        return self._events

    def gammas(self, plain: bool = False):
        """The γ rays of the real part, in coincidence with the counted particles
        (:class:`~physim.nuclear.gamma_events.GammaEvents`), with add-back and shields as the setup has them, or
        without (``plain``). None without γ rays (an elastic setup, a source run)."""
        if plain not in self._gammas:
            prefix = "gp_" if plain else "g_"
            a = self._arrays()
            if prefix + "particle" not in a:
                if plain and "g_particle" in a:  # no add-back or shields: the plain chain is the same
                    return self.gammas(False)
                return None
            self._gammas[plain] = _gamma_events(self, {k: a[prefix + k] for k in GAMMA_COLUMNS}, plain)
        return self._gammas[plain]

    def source(self):
        """The calibration-source spectra of a source run (:class:`~physim.nuclear.response.SourceRun`)."""
        if self.kind != "source":
            return None
        from .response import Response, SourceRun
        from .response import source as source_lines

        a = self._arrays()
        s = self.summary["source"]
        exp = self.experiment
        names = s["crystals"]
        responses = {}
        position = s.get("position")
        for i, gd in enumerate(exp.gamma_detectors):
            r = Response(exp, gd, source=position)
            for k, (label, _, _) in enumerate(r.elements):
                responses[(gd.name or f"G{i + 1}") + (f" {label}" if label else "")] = (r, k if len(r.elements) > 1
                                                                                       else None)
        return SourceRun(source_lines(s["nuclide"]), s["decays"], self.duration_s, a["s_edges"],
                         {n: a[f"s_spectrum_{i}"] for i, n in enumerate(names)},
                         {n: a[f"s_expected_{i}"] for i, n in enumerate(names)}, responses)


def _segment_directions(experiment: Experiment, cols: dict) -> tuple:
    """The lab θ and φ (degrees) of each particle row from the centre of the segment it hit."""
    array = Array.from_experiment(experiment)
    theta = np.zeros(len(cols["detector"]))
    phi = np.zeros(len(cols["detector"]))
    if not len(theta):
        return theta, phi
    keys = np.stack([cols["detector"], cols["segment_i"], cols["segment_j"]], axis=1)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    inverse = np.asarray(inverse).reshape(-1)
    th_u, ph_u = np.zeros(len(unique)), np.zeros(len(unique))
    for n, (d, i, j) in enumerate(unique):
        c = array.geometries[int(d)].segment_centre((int(i), int(j)), weighted=True)
        th_u[n] = math.degrees(math.atan2(math.hypot(c[0], c[1]), c[2]))
        ph_u[n] = math.degrees(math.atan2(c[1], c[0]))
    return th_u[inverse], ph_u[inverse]


def _gamma_events(run: RunData, a: dict, plain: bool):
    from .gamma_events import GammaEvents, backgrounds, crystals_of, doppler_correct

    exp = run.experiment
    ev = run.events()
    crystals, _ = crystals_of(exp, plain)
    bg = backgrounds(exp, crystals)
    rows = a["particle"].astype(np.int64)
    which = a["crystal"].astype(np.int64)
    measured = a["measured"].astype(float)
    thresholds = np.array([c.threshold for c in crystals])
    corrected = doppler_correct(exp, ev.columns, rows, which, measured)
    e0 = exp.excitation.energy_mev
    cols = {"event": ev.columns["event"][rows], "particle": rows, "crystal": which,
            "energy0": np.full(len(rows), e0), "deposited": a["deposited"].astype(float), "measured": measured,
            "counted": measured >= thresholds[which] if len(rows) else np.zeros(0, dtype=bool),
            "theta": a["theta"].astype(float), "phi": a["phi"].astype(float),
            "corrected_projectile": corrected["ejectile"], "corrected_recoil": corrected["recoil"],
            "weight": np.full(len(rows), 1.0 / run.real_s if run.real_s > 0 else 0.0)}
    notes = [f"{len(rows)} γ rays in the real part of run {run.number} ({duration_label(run.real_s)}); each counts "
             "once."]
    from .gamma import excitation_of

    return GammaEvents(cols, ev, crystals, int(run.summary["seed"]), e0, excitation_of(exp).scheme.nuclide,
                       bg["window_s"], bg["dead_time_s"], bg["live_fraction"], bg["singles_rate"],
                       bg["singles_spectrum"], bg["singles_edges"], bg["particle_rate"], notes)


# -- taking a run -------------------------------------------------------------------------------------------------

class RunTaker:
    """A beam (or alignment) run being taken: the generator, the accumulated compact columns and the live
    counters. :meth:`step` adds one chunk; :meth:`write` stores the run. Stop by setting :attr:`stop`."""

    def __init__(self, folder: Path, number: int, kind: str, experiment: Experiment, duration_s: float,
                 budget_s: float, seed: int, extra: Optional[dict] = None, plan: Optional[dict] = None,
                 previous: Optional[RunData] = None):
        self.folder, self.number, self.kind = folder, number, kind
        self.experiment = experiment
        self.duration_s = float(duration_s)
        self.budget_s = float(min(budget_s, MAX_BUDGET_S))
        self.seed = int(seed)
        self.extra = dict(extra or {})
        self.plan = plan
        self.stop = False
        self.started = time.time()
        self.created = _now()
        self.rates = Rates(experiment)
        self.gamma = experiment.excitation is not None and bool(experiment.gamma_detectors)
        self.plain_too = self.gamma and any(gd.addback or gd.shield for gd in experiment.gamma_detectors)
        self.parts: dict = {k: [] for k in PARTICLE_COLUMNS}
        self.gparts: dict = {p + k: [] for p in (("g_", "gp_") if self.plain_too else ("g_",)) for k in GAMMA_COLUMNS}
        self.n_rows = 0
        if previous is not None:  # an extension: continue the same stream and append
            self.runner = BeamRunner(experiment, self.seed, previous.summary["generator"])
            self.runner._next_event = previous.summary.get("next_event", 0)
            a = previous._arrays()
            for k in self.parts:
                self.parts[k].append(a[k])
            for k in self.gparts:
                if k in a:
                    self.gparts[k].append(a[k])
            self.n_rows = len(a["p_event_step"])
            self.last_event = int(np.sum(a["p_event_step"].astype(np.int64)))
            self.created = previous.summary.get("created", self.created)
            self.counters = copy.deepcopy(previous.summary.get("counters_raw") or self._empty_counters())
            self.extended = previous.summary.get("extended", 0) + 1
        else:
            self.runner = BeamRunner(experiment, self.seed)
            self.last_event = 0
            self.counters = self._empty_counters()
            self.extended = 0
        self.real_target = min(self.duration_s, self.budget_s)

    def _empty_counters(self) -> dict:
        return {"counts": {d: 0 for d in self.runner.detectors}, "gamma": {}, "coincidences": {}}

    @property
    def real_s(self) -> float:
        return self.runner.real_s

    @property
    def done(self) -> bool:
        return self.stop or self.real_s >= self.real_target * (1 - 1e-12)

    def cpu_estimate_s(self) -> float:
        """A rough CPU time for the real part, from the generator's speed (about 6 million reactions per second
        on all cores) and the γ chain."""
        reactions = self.real_target / self.runner.seconds_per_reaction()
        return reactions / 6e6 + 0.7 * max(1, reactions / MAX_CHUNK)

    def step(self) -> None:
        """Take one more chunk of the real part."""
        if self.done:
            return
        cols, _ = self.runner.next_chunk(self.real_target - self.real_s)
        if not cols or not len(cols.get("event", ())):
            return
        n = len(cols["event"])
        for key, x in _store_particles(cols, self.last_event).items():
            self.parts[key].append(x)
        self.last_event = int(cols["event"][-1])
        names = self.runner.detectors
        counts = np.bincount(cols["detector"].astype(np.int64), minlength=len(names))
        for d, k in zip(names, counts):
            self.counters["counts"][d] += int(k)
        if self.gamma:
            ev = Events(cols, 0, self.seed, names, self.runner.channels, 1.0, 0.0)
            ev.columns["weight"] = np.ones(n)
            chunk_seed = self.seed * 100_003 + self.runner.chunk
            for prefix, plain in (("g_", False), ("gp_", True)) if self.plain_too else (("g_", False),):
                from .gamma_events import simulate_gammas

                g = simulate_gammas(self.experiment, seed=chunk_seed, particle_events=ev, gammas_per_event=1,
                                    plain=plain, rates=self.rates)
                rows = g["particle"]
                part = {"particle": rows + self.n_rows, "crystal": g["crystal"], "deposited": g["deposited"],
                        "measured": g["measured"], "theta": g["theta"], "phi": g["phi"]}
                for k, t in GAMMA_COLUMNS.items():
                    self.gparts[prefix + k].append(np.asarray(part[k]).astype(t))
                if not plain:
                    gdet = [g.crystals[c].name for c in g["crystal"]]
                    counted = g["counted"]
                    for c, ok in zip(gdet, counted):
                        if ok:
                            self.counters["gamma"][c] = self.counters["gamma"].get(c, 0) + 1
                    parent = [c.rsplit(" ", 1)[0] if g.crystals[k].element is not None else c
                              for c, k in zip(gdet, g["crystal"])]
                    pdet = cols["detector"][rows]
                    for gd, pd, ok in zip(parent, pdet, counted):
                        if ok:
                            key = f"{names[int(pd)]}|{gd}"
                            self.counters["coincidences"][key] = self.counters["coincidences"].get(key, 0) + 1
        self.n_rows += n

    def progress(self) -> dict:
        """The live counters: experiment time, real part done, counts and rates per detector, γ singles per
        crystal, coincidences."""
        real = self.real_s
        return {"number": self.number, "kind": self.kind, "duration_s": self.duration_s, "real_target_s": self.real_target,
                "real_s": real, "fraction": real / self.real_target if self.real_target else 1.0,
                "counts": dict(self.counters["counts"]),
                "rates_per_s": {d: c / real if real else 0.0 for d, c in self.counters["counts"].items()},
                "gamma": dict(self.counters["gamma"]), "coincidences": dict(self.counters["coincidences"]),
                "elapsed_s": time.time() - self.started, "events": self.n_rows, "stopping": self.stop}

    def write(self) -> RunData:
        """Store the run (setup snapshot, compact events, summary) and return it."""
        self.folder.mkdir(parents=True, exist_ok=True)
        self.experiment.save(self.folder / "setup.toml")
        arrays = {k: (np.concatenate(v) if v else np.zeros(0, dtype=PARTICLE_COLUMNS[k][1]))
                  for k, v in self.parts.items()}
        for k, v in self.gparts.items():
            arrays[k] = np.concatenate(v) if v else np.zeros(0, dtype=GAMMA_COLUMNS[k.split("_", 1)[1]])
        np.savez_compressed(self.folder / "events.npz", **arrays)
        real = self.real_s
        stopped = self.stop and real < self.real_target * (1 - 1e-9)
        duration = real if stopped else self.duration_s
        live = None
        if self.gamma:
            from .gamma_events import backgrounds, crystals_of

            live = backgrounds(self.experiment, crystals_of(self.experiment)[0], self.rates)["live_fraction"]
        counts = self.counters["counts"]
        summary = {
            "number": self.number, "kind": self.kind, "status": "stopped" if stopped else "complete",
            "created": self.created, "finished": _now(), "seed": self.seed, "duration_s": duration,
            "real_s": real, "budget_s": self.budget_s, "scaled": duration > real * (1 + 1e-9),
            "generated": self.runner.generated, "generator": self.runner.state(),
            "next_event": getattr(self.runner, "_next_event", 0), "extended": self.extended,
            "particles": int(len(arrays["p_event_step"])), "gammas": int(len(arrays.get("g_particle", ()))),
            "detectors": list(self.runner.detectors), "channels": list(self.runner.channels),
            "counts": counts, "rates_per_s": {d: c / real if real else 0.0 for d, c in counts.items()},
            "counts_in_run": {d: c * duration / real if real else 0.0 for d, c in counts.items()},
            "segments": _segment_counts(arrays, self.runner.detectors),
            "gamma_counts": self.counters["gamma"], "coincidences": self.counters["coincidences"],
            "live_fraction": live, "cpu_s": time.time() - self.started,
            "counters_raw": self.counters,
            "plan": self.plan,
        }
        summary.update(self.extra)
        (self.folder / "summary.json").write_text(json.dumps(summary, indent=1, default=_json_default),
                                                  encoding="utf-8")
        return RunData.load(self.folder)


def _segment_counts(arrays: dict, detectors: list) -> dict:
    out = {}
    det, i = arrays["p_detector"].astype(np.int64), arrays["p_segment_i"].astype(np.int64)
    for k, name in enumerate(detectors):
        m = det == k
        if m.any():
            out[name] = {str(int(a)): int(b) for a, b in zip(*np.unique(i[m], return_counts=True))}
    return out


def _json_default(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (set, tuple)):
        return list(x)
    return str(x)


class ExperimentFolder:
    """An experiment on disk: its current setup (``experiment.toml``) and its runs (``runs/NN-...``)."""

    def __init__(self, path: Union[str, Path]):
        self.path = Path(path).expanduser()

    @classmethod
    def create(cls, experiment: Experiment, name: Optional[str] = None,
               root: Union[str, Path, None] = None) -> ExperimentFolder:
        """A new experiment folder under ``root`` (default :func:`home`), named after ``name`` (by default the beam
        and target); a number is added if the name is taken."""
        root = Path(root).expanduser() if root is not None else home()
        name = name or default_name(experiment)
        base = slug(name)
        path, k = root / base, 2
        while path.exists():
            path, k = root / f"{base}-{k}", k + 1
        path.mkdir(parents=True)
        (path / "runs").mkdir()
        f = cls(path)
        f.save_setup(experiment)
        f._write_meta({"name": name, "created": _now()})
        return f

    def _meta(self) -> dict:
        p = self.path / "experiment.json"
        return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {}

    def _write_meta(self, meta: dict) -> None:
        (self.path / "experiment.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")

    @property
    def name(self) -> str:
        return self._meta().get("name") or self.path.name

    def rename(self, name: str) -> None:
        meta = self._meta()
        meta["name"] = name
        self._write_meta(meta)

    def setup(self) -> Experiment:
        return Experiment.load(self.path / "experiment.toml")

    def save_setup(self, experiment: Experiment) -> None:
        experiment.save(self.path / "experiment.toml")

    def run_folders(self) -> list:
        d = self.path / "runs"
        if not d.is_dir():
            return []
        return sorted((p for p in d.iterdir() if (p / "summary.json").is_file()), key=lambda p: p.name)

    def runs(self) -> list:
        """The summaries of all runs, in order (small; the events stay on disk)."""
        out = []
        for p in self.run_folders():
            s = json.loads((p / "summary.json").read_text(encoding="utf-8"))
            s["_folder"] = str(p)
            s["_label"] = p.name
            s["_mtime"] = (p / "summary.json").stat().st_mtime
            out.append(s)
        return out

    def next_number(self) -> int:
        nums = [int(p.name.split("-")[0]) for p in self.run_folders() if p.name.split("-")[0].isdigit()]
        return max(nums, default=0) + 1

    def load_run(self, number: int) -> RunData:
        for p in self.run_folders():
            if p.name.split("-")[0].isdigit() and int(p.name.split("-")[0]) == int(number):
                return RunData.load(p)
        raise KeyError(f"no run {number} in {self.path}")

    def run_folder(self, number: int, kind: str, duration_s: float, extra: Optional[dict] = None) -> Path:
        if kind == "source":
            what = f"{slug((extra or {}).get('source', {}).get('nuclide', 'source'))}-source"
        else:
            what = f"{kind}-{duration_label(duration_s)}"
        return self.path / "runs" / f"{number:02d}-{what}"

    def take_run(self, kind: str = "beam", duration=None, budget=DEFAULT_BUDGET_S, seed: Optional[int] = None,
                 experiment: Optional[Experiment] = None, plan: Optional[dict] = None,
                 progress: Optional[Callable] = None, **options) -> RunData:
        """Take a run with the folder's setup (or ``experiment``) and store it. See :func:`start`."""
        taker = start(self, kind, duration, budget, seed, experiment, plan, **options)
        if isinstance(taker, RunData):
            return taker
        while not taker.done:
            taker.step()
            if progress is not None:
                progress(taker.progress())
        return taker.write()


def start(folder: ExperimentFolder, kind: str = "beam", duration=None, budget=DEFAULT_BUDGET_S,
          seed: Optional[int] = None, experiment: Optional[Experiment] = None, plan: Optional[dict] = None,
          source: str = "152Eu", activity: str = "37 kBq", position=None, offset_mm: float = 0.0):
    """Begin a run in ``folder``: a :class:`RunTaker` to step for a beam or alignment run, or the stored
    :class:`RunData` of a source run (which is computed at once). ``duration`` defaults to the setup's beam time;
    ``budget`` is the real part (at most :data:`MAX_BUDGET_S`); ``seed`` defaults to the run's number. A source
    run takes ``source``, ``activity`` and ``position`` (mm from the target); an alignment run ``offset_mm``, the
    target's assumed offset along the beam."""
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    exp = experiment if experiment is not None else folder.setup()
    exp.validate()
    duration_s = _seconds(duration) if duration is not None else _seconds(exp.run.beam_time)
    if not duration_s > 0:
        raise ValueError("the run needs a duration")
    budget_s = min(_seconds(budget), MAX_BUDGET_S)
    number = folder.next_number()
    seed = number if seed is None else int(seed)
    if kind == "source":
        return _source_run(folder, number, exp, duration_s, seed, source, activity, position, plan)
    extra = {"assumed_offset_mm": float(offset_mm)} if kind == "alignment" else {}
    path = folder.run_folder(number, kind, duration_s)
    return RunTaker(path, number, kind, exp, duration_s, budget_s, seed, extra, plan)


def extend(run: RunData, more) -> RunTaker:
    """Continue ``run`` for ``more`` beam time (text with a unit, or seconds) with its own setup and seed stream;
    the real part grows up to the run's budget, the rest is scaled. Step the result and :meth:`RunTaker.write`
    it, which rewrites the run's folder."""
    if run.kind == "source":
        raise ValueError("a source run cannot be extended; take another")
    extra = {k: run.summary[k] for k in ("assumed_offset_mm",) if k in run.summary}
    taker = RunTaker(run.folder, run.number, run.kind, run.experiment, run.duration_s + _seconds(more),
                     run.summary["budget_s"], run.summary["seed"], extra, run.summary.get("plan"), previous=run)
    return taker


def _source_run(folder, number, exp, duration_s, seed, nuclide, activity, position, plan) -> RunData:
    from .response import source_run

    sr = source_run(exp, nuclide, activity, f"{duration_s} s", seed, position=position)
    names = sr.names()
    path = folder.run_folder(number, "source", duration_s, {"source": {"nuclide": nuclide}})
    path.mkdir(parents=True, exist_ok=True)
    exp.save(path / "setup.toml")
    arrays = {"s_edges": sr.edges}
    for i, n in enumerate(names):
        arrays[f"s_spectrum_{i}"] = sr.spectra[n].astype(np.float32)
        arrays[f"s_expected_{i}"] = np.asarray(sr.expected[n], dtype=np.float32)
    np.savez_compressed(path / "events.npz", **arrays)
    summary = {"number": number, "kind": "source", "status": "complete", "created": _now(), "finished": _now(),
               "seed": seed, "duration_s": duration_s, "real_s": duration_s, "scaled": False,
               "source": {"nuclide": sr.source.name if hasattr(sr.source, "name") else nuclide,
                          "activity": activity, "position": list(position) if position is not None else None,
                          "decays": sr.decays, "crystals": names},
               "detectors": [], "channels": [], "counts": {},
               "gamma_counts": {n: float(sr.spectra[n].sum()) for n in names},
               "rates_per_s": {n: float(sr.spectra[n].sum()) / duration_s for n in names}, "plan": plan}
    (path / "summary.json").write_text(json.dumps(summary, indent=1, default=_json_default), encoding="utf-8")
    return RunData.load(path)


class Background:
    """A run taken on a thread, so a front end stays responsive: :attr:`taker` has the live counters,
    :meth:`stop` asks it to stop after the current chunk, :attr:`result` is the stored run when done."""

    def __init__(self, taker: RunTaker, done: Optional[Callable] = None):
        self.taker = taker
        self.result: Optional[RunData] = None
        self.error: Optional[BaseException] = None
        self._done = done
        self.thread = threading.Thread(target=self._work, daemon=True)
        self.thread.start()

    def _work(self):
        try:
            while not self.taker.done:
                self.taker.step()
            self.result = self.taker.write()
        except BaseException as e:  # noqa: BLE001 - reported to the front end
            self.error = e
        finally:
            if self._done is not None:
                self._done(self)

    def stop(self) -> None:
        self.taker.stop = True

    def join(self, timeout: Optional[float] = None) -> Optional[RunData]:
        self.thread.join(timeout)
        if self.error is not None:
            raise self.error
        return self.result

    @property
    def running(self) -> bool:
        return self.thread.is_alive()


__all__ = ["DEFAULT_BUDGET_S", "KINDS", "MAX_BUDGET_S", "Background", "BeamRunner", "ExperimentFolder", "RunData",
           "RunInProgress", "RunTaker", "default_name", "duration_label", "extend", "home", "list_experiments",
           "physics_changes", "slug", "start"]
