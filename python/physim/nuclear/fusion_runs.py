"""Simulated runs of a fusion-evaporation setup: residue events, their γ rays, a trigger and γγ coincidences.

A fusion run is stored like a particle run (:mod:`physim.nuclear.runs`), in a run folder of the experiment::

    runs/03-fusion-10min/   setup.toml  events.npz  summary.json

and it is started, stopped, extended and loaded the same way::

    from physim.nuclear import Experiment, fusion_runs
    from physim.nuclear.runs import ExperimentFolder

    lab = ExperimentFolder.create(Experiment.example("fusion_dy158"))
    taker = fusion_runs.start(lab, duration="1 h", budget="20 s")   # fold ≥ 2 by default
    while not taker.done:
        taker.step()                       # one block of residues; taker.stop = True ends the run early
    run = taker.write()                    # a FusionRunData
    run.gammas().spectrum()                # the Doppler-corrected γ rays of the triggered events
    counts, edges = run.gated(218, 3)      # γ rays in coincidence with the 218 keV 4⁺ → 2⁺ line of ¹⁵⁸Dy
    fusion_runs.extend(run, "1 h")         # a taker that continues the same stream

**Events.** The residues come from :func:`~physim.nuclear.residues.simulate_residues` and their γ rays from
:func:`~physim.nuclear.residues.simulate_residue_gammas`, with one cascade per residue (``gammas_per_event=1``), so
the γ rays of one event really are in coincidence. Every residue stands for one residue: for a real part of t
seconds the run generates a Poisson number of residues of mean σ_total × luminosity × t, in blocks of
:data:`BLOCK_EVENTS` (each block with its own seed, so a run that is extended carries on where it stopped).
The rest of the duration is scaled, as in a particle run (:meth:`~physim.nuclear.runs.RunData.scaled`).

**Trigger.** An event is recorded when at least ``fold`` crystals count a γ ray (γ-ray fold, 2 by default: γγ)
and the setup's ancillaries accept it:

- a recoil separator: the residue flies within its angular acceptance (by angle only, as the plan counts it;
  its transmission is not computed);
- a γ calorimeter with a ``fold_min``: at least that many of its elements fire. Each γ ray of the residue's
  cascade (drawn again from the same level scheme, so the fold has the right distribution but is not the same
  draw as the germanium's γ rays) and :data:`STATISTICAL_GAMMAS` statistical γ rays fire an element with the
  calorimeter's ``efficiency`` (:data:`CALORIMETER_EFFICIENCY` where the setup gives none), each at an element
  drawn at random; the fold is the number of distinct elements;
- a neutron wall and a charged-particle ball are recorded per event (detected neutrons; protons and α particles
  in the ball's polar range above its threshold) and gate the trigger only when asked for (``neutrons_min``,
  ``charged_min``). The wall's efficiency is its :meth:`~physim.nuclear.neutrons.NeutronWall.detection_probability`;
  cross-talk and γ rays taken for neutrons are left out (the channel selection on the Plan has them).

``gates=False`` turns the ancillaries' gates off (they are still recorded).

**What is stored.** Per recorded event: event number (as a step from the previous one), channel, γ-ray fold,
detected neutrons, charged particles, calorimeter fold, whether the separator accepted the residue, its speed,
direction and depth, its entry spin. Per counted γ ray of a recorded event: the event's row, crystal, initial
and final level, energy in the residue's frame, deposited and measured energy, whether it left its full energy,
lab θ and φ (and, from a residue generator that has it, whether the residue had stopped). γ rays below their
crystal's threshold are not recorded. Two γ rays in one crystal are kept apart (the residues' γ rays have no
summing); the fold counts crystals. The Doppler correction (to β_CN along the beam, as
:func:`~physim.nuclear.residues.doppler_correct`) is derived on loading.

**Memory.** A run uses at most half of the RAM free when it starts: the block size is cut to fit if need be
(about :data:`BYTES_PER_EVENT` per residue while a block is made), and a run whose stored events would grow past
that limit stops, saying so in its summary.

**Plunger.** ``distance_um`` (µm, the setup's distance without the offset) is stored in ``summary.json`` and passed
to :func:`~physim.nuclear.residues.simulate_residue_gammas` when it is set: one run per distance makes a distance
series.
"""

from __future__ import annotations

import copy
import inspect
import json
import math
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import numpy as np

from .experiment import Experiment
from .runs import (DEFAULT_BUDGET_S, MAX_BUDGET_S, ExperimentFolder, RunData, _json_default, _now, _seconds,
                   duration_label)

#: The kind of run, as written in ``summary.json``.
KIND = "fusion"
#: The γ-ray fold the trigger asks for by default (γγ coincidences).
DEFAULT_FOLD = 2
#: Residues per block: each block has its own seed, and one step of a run makes one block.
BLOCK_EVENTS = 100_000
#: The fewest residues per block, however little RAM is free.
MIN_BLOCK_EVENTS = 5_000
#: Peak memory while a block is made, bytes per residue (measured: about 2.0–2.4 kB for ``fusion_dy158``).
BYTES_PER_EVENT = 4_000
#: CPU time per residue, s (measured: about 12 µs for ``fusion_dy158``, residues and γ rays).
CPU_S_PER_EVENT = 1.5e-5
#: Statistical γ rays per event a γ calorimeter sees besides the discrete cascade (as the channel selection).
STATISTICAL_GAMMAS = 4.0
#: The chance that a γ ray fires a calorimeter element, where the setup gives no ``efficiency``.
CALORIMETER_EFFICIENCY = 0.7

#: The stored columns of the recorded events: {key in events.npz: type}.
EVENT_COLUMNS = {
    "r_event_step": np.int32, "r_channel": np.int8, "r_fold": np.int8, "r_neutrons": np.int8,
    "r_charged": np.int8, "r_calorimeter": np.int16, "r_separator": np.int8, "r_beta": np.float32,
    "r_theta": np.float32, "r_phi": np.float32, "r_depth": np.float32, "r_entry_spin": np.float16,
}
#: The stored columns of the γ rays: {key in events.npz: type}. ``g_row`` is the row of the event.
GAMMA_COLUMNS = {
    "g_row": np.int32, "g_crystal": np.int16, "g_initial": np.int16, "g_final": np.int16, "g_energy0": np.float32,
    "g_deposited": np.float32, "g_measured": np.float32, "g_full": np.int8, "g_theta": np.float32,
    "g_phi": np.float32,
}
#: γ-ray columns stored when the residue generator gives them: {key: (column of ResidueGammas, type)}.
OPTIONAL_GAMMA_COLUMNS = {"g_stopped": ("stopped", np.int8)}


def _memory_limit() -> int:
    """Half of the RAM free now, bytes (the user's rule for a run)."""
    from .shellmodel_spaces import available_memory_bytes

    try:
        return int(available_memory_bytes()) // 2
    except OSError:
        return 2 * 1024**3


def _takes_distance() -> bool:
    from .residues import simulate_residue_gammas

    return "distance_um" in inspect.signature(simulate_residue_gammas).parameters


# -- the generator ------------------------------------------------------------------------------------------------

class FusionRunner:
    """Residue events with their γ rays, block by block, and the trigger.

    Block k has seed (seed, k): its number of residues is Poisson, of mean the residue rate × the beam time it
    covers (at most :attr:`block` / rate), and its residues and γ rays come from the residue generator with a seed
    of its own. ``state`` (from :meth:`state`) continues a run."""

    def __init__(self, experiment: Experiment, seed: int, trigger: dict, distance_um: Optional[float] = None,
                 state: Optional[dict] = None, block: Optional[int] = None):
        from .residues import ChannelTable, simulate_residues

        self.experiment = experiment
        self.seed = int(seed)
        self.trigger = dict(trigger)
        self.distance_um = distance_um
        probe = simulate_residues(experiment, 1, self.seed)
        self.table = ChannelTable.from_experiment(experiment)
        #: Residues per second (σ_total × luminosity) and the Doppler correction's β.
        self.rate = float(probe.expected())
        self.beta_nominal = float(probe.beta_nominal)
        from .gamma_events import crystals_of

        self.crystals = crystals_of(experiment, plain=True)[0]
        self.block = int(state["block"]) if state else int(block or BLOCK_EVENTS)
        if state:
            self.chunk, self.real_s, self.generated = state["chunk"], state["real_s"], state["generated"]
            self.next_event = state["next_event"]
        else:
            self.chunk, self.real_s, self.generated, self.next_event = 0, 0.0, 0, 0
        self._cascades: dict = {}
        self._wall = experiment.neutron_wall.model() if experiment.neutron_wall is not None else None

    def state(self) -> dict:
        return {"block": self.block, "chunk": self.chunk, "real_s": self.real_s, "generated": self.generated,
                "next_event": self.next_event}

    def seconds_per_block(self) -> float:
        return self.block / self.rate if self.rate > 0 else math.inf

    def next_chunk(self, seconds_left: float) -> tuple:
        """The next block of at most ``seconds_left`` of beam: (the block's events and γ rays, seconds covered).
        The events are a dict of per-residue arrays (all residues, the trigger's verdict in ``"trigger"``), the γ
        rays a dict of per-γ-ray arrays (counted γ rays only, ``"event"`` the index into the block's residues)."""
        from .residues import simulate_residue_gammas, simulate_residues

        t = min(self.seconds_per_block(), seconds_left)
        k = self.chunk
        self.chunk += 1
        rng = np.random.default_rng([self.seed, k, 1])
        n = int(rng.poisson(self.rate * t)) if self.rate > 0 and t > 0 else 0
        self.real_s += t
        if n < 1:
            return {}, {}, t
        block_seed = self.seed * 1_000_003 + k
        ev = simulate_residues(self.experiment, n, block_seed)
        kw = {"distance_um": self.distance_um} if self.distance_um is not None else {}
        g = simulate_residue_gammas(self.experiment, seed=block_seed, residue_events=ev, gammas_per_event=1, **kw)
        c = ev.columns
        gc = g.columns
        counted = np.asarray(gc["counted"], dtype=bool)
        gam = {"event": np.asarray(gc["event"])[counted].astype(np.int64)}
        for name in ("crystal", "initial", "final", "energy0", "deposited", "measured", "theta", "phi"):
            gam[name] = np.asarray(gc[name])[counted]
        gam["full"] = (np.abs(np.asarray(gc["deposited"]) - np.asarray(gc["energy_lab"])) < 1e-9)[counted]
        for _, (name, _) in OPTIONAL_GAMMA_COLUMNS.items():
            if name in gc:
                gam[name] = np.asarray(gc[name])[counted]
        # The γ-ray fold: crystals with a counted γ ray, per residue.
        fold = np.zeros(n, dtype=np.int64)
        if len(gam["event"]):
            nc = max(len(self.crystals), 1)
            pairs = np.unique(gam["event"] * nc + gam["crystal"].astype(np.int64))
            fold = np.bincount(pairs // nc, minlength=n)
        arng = np.random.default_rng([self.seed, k, 7])
        neutrons, charged = self._particles(ev, arng)
        calorimeter = self._calorimeter(ev, arng)
        separator = self._separator(ev)
        tr = self.trigger
        trig = fold >= tr["fold"]
        if tr["neutrons_min"] is not None:
            trig &= neutrons >= tr["neutrons_min"]
        if tr["charged_min"] is not None:
            trig &= charged >= tr["charged_min"]
        if tr["gates"]:
            if tr["calorimeter_fold"] is not None:
                trig &= calorimeter >= tr["calorimeter_fold"]
            if tr["separator"]:
                trig &= separator
        events = {"event": self.next_event + np.arange(n), "channel": np.asarray(c["channel"]), "fold": fold,
                  "neutrons": neutrons, "charged": charged, "calorimeter": calorimeter, "separator": separator,
                  "beta": np.asarray(c["beta"]), "theta": np.asarray(c["theta"]), "phi": np.asarray(c["phi"]),
                  "depth": np.asarray(c["depth"]), "entry_spin": np.asarray(c["entry_spin"]), "trigger": trig}
        self.next_event += n
        self.generated += n
        return events, gam, t

    # -- the ancillaries ------------------------------------------------------------------------------------------

    def _particles(self, ev, rng) -> tuple:
        """Detected neutrons (the wall) and charged particles (the ball) per residue; zeros without them."""
        from .residues import LIGHT

        n = len(ev)
        neutrons, charged = np.zeros(n, dtype=np.int64), np.zeros(n, dtype=np.int64)
        p = ev.particles
        if not len(p["event"]):
            return neutrons, charged
        kinds = np.array(list(LIGHT))[np.asarray(p["kind"], dtype=np.int64)]
        th, ph = np.radians(p["theta"]), np.radians(p["phi"])
        event = np.asarray(p["event"], dtype=np.int64)
        if self._wall is not None:
            sel = kinds == "n"
            if sel.any():
                dirs = np.stack([np.sin(th[sel]) * np.cos(ph[sel]), np.sin(th[sel]) * np.sin(ph[sel]),
                                 np.cos(th[sel])], axis=1)
                prob = self._wall.detection_probability(np.asarray(p["energy"])[sel], dirs)
                hit = rng.uniform(0.0, 1.0, int(sel.sum())) < prob
                neutrons = np.bincount(event[sel][hit], minlength=n)
        ball = self.experiment.charged_particle_ball
        if ball is not None:
            from .exit_channels import PARTICLE_THRESHOLD_MEV
            from .experiment import _q

            lo = _q(ball.theta_min).to("deg") if ball.theta_min is not None else 0.0
            hi = _q(ball.theta_max).to("deg") if ball.theta_max is not None else 180.0
            thr = _q(ball.threshold).to("MeV") if ball.threshold is not None else PARTICLE_THRESHOLD_MEV
            theta = np.asarray(p["theta"])
            sel = (kinds != "n") & (theta >= lo) & (theta <= hi) & (np.asarray(p["energy"]) >= thr)
            charged = np.bincount(event[sel], minlength=n)
        return neutrons, charged

    def _calorimeter(self, ev, rng) -> np.ndarray:
        """Elements of the γ calorimeter that fire, per residue; zeros without one."""
        from .residues import ResidueCascade, residue_scheme

        cal = self.experiment.gamma_calorimeter
        n = len(ev)
        if cal is None:
            return np.zeros(n, dtype=np.int64)
        from .experiment import _q

        eff = _q(cal.efficiency).to("%") / 100 if cal.efficiency is not None else CALORIMETER_EFFICIENCY
        mult = np.zeros(n, dtype=np.int64)
        channel = np.asarray(ev["channel"])
        spin = np.asarray(ev["entry_spin"], dtype=float)
        for i, ch in enumerate(self.table.channels):
            rows = np.flatnonzero(channel == i)
            if not len(rows):
                continue
            if i not in self._cascades:
                scheme = residue_scheme(self.experiment, ch.residue)
                self._cascades[i] = (ResidueCascade(scheme, self.table.alignment, self.table.side_feeding)
                                     if scheme is not None and scheme.transitions else None)
            cas = self._cascades[i]
            if cas is not None:
                src = cas.run(spin[rows], rng)["source"]
                mult[rows] += np.bincount(src, minlength=len(rows))
        low = math.floor(STATISTICAL_GAMMAS)
        mult += low + (rng.uniform(0.0, 1.0, n) < STATISTICAL_GAMMAS - low)
        total = int(mult.sum())
        owner = np.repeat(np.arange(n), mult)
        fires = rng.uniform(0.0, 1.0, total) < eff
        element = rng.integers(0, cal.elements, total)
        code = np.unique(owner[fires] * cal.elements + element[fires])
        return np.bincount(code // cal.elements, minlength=n)

    def _separator(self, ev) -> np.ndarray:
        """Whether the recoil separator accepts each residue (by angle); all True without one."""
        sep = self.experiment.recoil_separator
        if sep is None:
            return np.ones(len(ev), dtype=bool)
        from .experiment import _q

        return np.asarray(ev["theta"]) <= _q(sep.theta_acceptance).to("deg")


# -- taking a run -------------------------------------------------------------------------------------------------

def _trigger(experiment: Experiment, fold, neutrons_min, charged_min, gates: bool) -> dict:
    fold = DEFAULT_FOLD if fold is None else int(fold)
    if fold < 1:
        raise ValueError("the trigger's γ-ray fold must be at least 1")
    for name, x in (("neutrons_min", neutrons_min), ("charged_min", charged_min)):
        if x is not None and int(x) < 0:
            raise ValueError(f"{name} must not be negative")
    if neutrons_min is not None and experiment.neutron_wall is None:
        raise ValueError("neutrons_min needs a [neutron_wall] in the setup")
    if charged_min is not None and experiment.charged_particle_ball is None:
        raise ValueError("charged_min needs a [charged_particle_ball] in the setup")
    cal = experiment.gamma_calorimeter
    return {"fold": fold, "neutrons_min": None if neutrons_min is None else int(neutrons_min),
            "charged_min": None if charged_min is None else int(charged_min), "gates": bool(gates),
            "calorimeter_fold": cal.fold_min if cal is not None else None,
            "separator": experiment.recoil_separator is not None}


class FusionRunTaker:
    """A fusion run being taken, with the surface of :class:`~physim.nuclear.runs.RunTaker`: :meth:`step` adds one
    block, :attr:`done` says when the real part is complete, setting :attr:`stop` ends it early, :meth:`write`
    stores the run, :meth:`progress` and :meth:`sample` are the live counters. They may be called from another thread
    while :meth:`step` runs: a block is made aside, and shown (its beam time, counters and events together) only when
    all of it is there."""

    kind = KIND

    def __init__(self, folder: Path, number: int, experiment: Experiment, duration_s: float, budget_s: float,
                 seed: int, trigger: dict, distance_um: Optional[float] = None, extra: Optional[dict] = None,
                 plan: Optional[dict] = None, previous: Optional[RunData] = None):
        self.folder, self.number = Path(folder), int(number)
        self.experiment = experiment
        self.duration_s = float(duration_s)
        self.budget_s = float(min(budget_s, MAX_BUDGET_S))
        self.seed = int(seed)
        self.trigger = dict(trigger)
        self.distance_um = distance_um
        self.extra = dict(extra or {})
        self.plan = plan
        self.stop = False
        self.started = time.time()
        self.created = _now()
        self.memory_limit = _memory_limit()
        self.stopped_for: Optional[str] = None
        self.parts: dict = {k: [] for k in list(EVENT_COLUMNS) + list(GAMMA_COLUMNS)}
        self.optional: set = set()
        if previous is not None:
            s = previous.summary
            self.runner = FusionRunner(experiment, self.seed, self.trigger, distance_um, state=s["generator"])
            a = previous._arrays()
            for k in self.parts:
                if k in a:
                    self.parts[k].append(a[k])
            for k in OPTIONAL_GAMMA_COLUMNS:
                if k in a:
                    self.optional.add(k)
                    self.parts[k] = [a[k]]
            self.n_rows = len(a["r_event_step"])
            self.last_event = int(np.sum(a["r_event_step"].astype(np.int64)))
            self.created = s.get("created", self.created)
            self.counters = copy.deepcopy(s["counters_raw"])
            self.extended = s.get("extended", 0) + 1
            self.stored_bytes = sum(v.nbytes for v in a.values())
        else:
            block = BLOCK_EVENTS
            fits = self.memory_limit // BYTES_PER_EVENT
            if fits < block:
                block = max(MIN_BLOCK_EVENTS, int(fits))
            self.runner = FusionRunner(experiment, self.seed, self.trigger, distance_um, block=block)
            self.n_rows = 0
            self.last_event = 0
            self.counters = self._empty_counters()
            self.extended = 0
            self.stored_bytes = 0
        #: Held while a finished block is shown, and while the front end reads what there is.
        self._lock = threading.Lock()
        #: The beam time (s) of the blocks that are stored: the runner's own is ahead of it during a step.
        self._real_s = self.runner.real_s
        self.real_target = min(self.duration_s, self.budget_s)
        self._last_sample: list = []

    def _empty_counters(self) -> dict:
        names = [c.name for c in self.runner.crystals]
        return {"residues": 0, "triggered": 0, "gamma": {n: 0 for n in names}, "singles": {n: 0 for n in names},
                "folds": [], "channels": {lab: {"residues": 0, "triggered": 0} for lab in self.runner.table.labels},
                "neutrons": [], "charged": [], "calorimeter": [], "separator": 0}

    @property
    def real_s(self) -> float:
        """The beam time simulated so far (s): that of the events stored."""
        return self._real_s

    @property
    def done(self) -> bool:
        return self._done_at(self.real_s)

    def _done_at(self, real_s: float) -> bool:
        return self.stop or real_s >= self.real_target * (1 - 1e-12)

    def cpu_estimate_s(self) -> float:
        """A rough CPU time for the real part, s: :data:`CPU_S_PER_EVENT` per residue."""
        return 1.0 + self.runner.rate * max(self.real_target - self.real_s, 0.0) * CPU_S_PER_EVENT

    def step(self) -> None:
        """Take one more block of the real part."""
        if self.done:
            return
        if self.stored_bytes + self.runner.block * BYTES_PER_EVENT > self.memory_limit:
            self.stop, self.stopped_for = True, "memory"
            return
        events, gam, _ = self.runner.next_chunk(self.real_target - self.real_s)
        if not events:
            with self._lock:
                self._real_s = self.runner.real_s
            return
        cnt = copy.deepcopy(self.counters)
        n = len(events["event"])
        names = [c.name for c in self.runner.crystals]
        cnt["residues"] += n
        _add_hist(cnt, "folds", events["fold"])
        for key in ("neutrons", "charged", "calorimeter"):
            _add_hist(cnt, key, events[key])
        cnt["separator"] += int(np.count_nonzero(events["separator"]))
        for k, x in enumerate(np.bincount(gam["crystal"].astype(np.int64), minlength=len(names))):
            cnt["singles"][names[k]] += int(x)
        labels = self.runner.table.labels
        trig = events["trigger"]
        for i, (a, b) in enumerate(zip(np.bincount(events["channel"], minlength=len(labels)),
                                       np.bincount(events["channel"][trig], minlength=len(labels)))):
            cnt["channels"][labels[i]]["residues"] += int(a)
            cnt["channels"][labels[i]]["triggered"] += int(b)
        rows = np.flatnonzero(trig)
        if not len(rows):
            with self._lock:
                self.counters = cnt
                self._real_s = self.runner.real_s
            return
        cnt["triggered"] += len(rows)
        # Rows of the recorded events, and their γ rays.
        new_row = np.full(n, -1, dtype=np.int64)
        new_row[rows] = self.n_rows + np.arange(len(rows))
        gsel = trig[gam["event"]]
        numbers = events["event"][rows]
        out = {"r_event_step": np.diff(numbers.astype(np.int64), prepend=self.last_event),
               "r_channel": events["channel"][rows], "r_fold": np.minimum(events["fold"][rows], 127),
               "r_neutrons": np.minimum(events["neutrons"][rows], 127),
               "r_charged": np.minimum(events["charged"][rows], 127),
               "r_calorimeter": events["calorimeter"][rows], "r_separator": events["separator"][rows],
               "r_beta": events["beta"][rows], "r_theta": events["theta"][rows], "r_phi": events["phi"][rows],
               "r_depth": events["depth"][rows], "r_entry_spin": events["entry_spin"][rows],
               "g_row": new_row[gam["event"][gsel]], "g_crystal": gam["crystal"][gsel],
               "g_initial": gam["initial"][gsel], "g_final": gam["final"][gsel],
               "g_energy0": gam["energy0"][gsel], "g_deposited": gam["deposited"][gsel],
               "g_measured": gam["measured"][gsel], "g_full": gam["full"][gsel], "g_theta": gam["theta"][gsel],
               "g_phi": gam["phi"][gsel]}
        for key, (name, _) in OPTIONAL_GAMMA_COLUMNS.items():
            if name in gam:
                out[key] = gam[name][gsel]
        types = {**EVENT_COLUMNS, **GAMMA_COLUMNS, **{k: t for k, (_, t) in OPTIONAL_GAMMA_COLUMNS.items()}}
        stored = {key: np.asarray(x).astype(types[key]) for key, x in out.items()}
        for k, x in enumerate(np.bincount(out["g_crystal"].astype(np.int64), minlength=len(names))):
            cnt["gamma"][names[k]] += int(x)
        sample = self._sample_from(events, gam, rows)
        with self._lock:  # all of the block at once
            for key, arr in stored.items():
                if key in OPTIONAL_GAMMA_COLUMNS:
                    self.optional.add(key)
                self.parts.setdefault(key, []).append(arr)
                self.stored_bytes += arr.nbytes
            self.counters = cnt
            self.last_event = int(numbers[-1])
            self.n_rows += len(rows)
            self._last_sample = sample
            self._real_s = self.runner.real_s

    def _sample_from(self, events: dict, gam: dict, rows: np.ndarray, n: int = 12) -> list:
        labels = self.runner.table.labels
        out = []
        for r in rows[::-1][:n]:
            m = gam["event"] == r
            out.append({"event": int(events["event"][r]), "channel": labels[int(events["channel"][r])],
                        "fold": int(events["fold"][r]),
                        "energies_kev": [round(float(e) * 1e3, 1) for e in gam["measured"][m]]})
        return out

    def progress(self) -> dict:
        """The live counters, in experiment time: the beam time simulated so far, the residues made and the events
        the trigger kept (with their rate), per crystal the γ rays in the kept events and the singles of all
        residues, the γ-ray fold distribution of all residues, and per channel its residues and kept events."""
        with self._lock:
            real, c, recorded = self._real_s, self.counters, self.n_rows
        trig = c["triggered"]
        scale = self.duration_s / real if real else 0.0
        return {"number": self.number, "kind": KIND, "duration_s": self.duration_s,
                "real_target_s": self.real_target, "scaled_s": max(self.duration_s - self.real_target, 0.0),
                "real_s": real, "fraction": real / self.real_target if self.real_target else 1.0,
                "real_done": self._done_at(real), "residues": c["residues"], "residues_per_s": self.runner.rate,
                "counts": {"trigger": trig}, "rates_per_s": {"trigger": trig / real if real else 0.0},
                "projected": {"trigger": trig * scale}, "trigger": dict(self.trigger),
                "gamma": dict(c["gamma"]), "gamma_singles": dict(c["singles"]), "folds": list(c["folds"]),
                "channels": copy.deepcopy(c["channels"]), "distance_um": self.distance_um,
                "elapsed_s": time.time() - self.started, "events": recorded, "stopping": self.stop}

    def sample(self, n: int = 12) -> list:
        """The last ``n`` recorded events of the latest block, newest first: {event, channel, fold, measured γ-ray
        energies (keV)}, for "watch events"."""
        with self._lock:
            return self._last_sample[:n]

    def write(self) -> FusionRunData:
        """Store the run (setup snapshot, compact events, summary) and return it."""
        self.folder.mkdir(parents=True, exist_ok=True)
        self.experiment.save(self.folder / "setup.toml")
        types = {**EVENT_COLUMNS, **GAMMA_COLUMNS, **{k: t for k, (_, t) in OPTIONAL_GAMMA_COLUMNS.items()}}
        arrays = {k: (np.concatenate(v) if v else np.zeros(0, dtype=types[k])) for k, v in self.parts.items()
                  if k in EVENT_COLUMNS or k in GAMMA_COLUMNS or k in self.optional}
        np.savez_compressed(self.folder / "events.npz", **arrays)
        real = self.real_s
        stopped = self.stop and real < self.real_target * (1 - 1e-9)
        duration = real if stopped else self.duration_s
        c = self.counters
        trig = c["triggered"]
        notes = [f"{c['residues']} residues in {duration_label(real)} of beam simulated event by event "
                 f"({self.runner.rate:.4g} per second); the trigger kept {trig}."]
        if self.stopped_for == "memory":
            notes.append(f"Stopped before the end of the real part: its events would have grown past half of the "
                         f"RAM free when the run started ({self.memory_limit / 1e9:.2f} GB).")
        if self.experiment.gamma_calorimeter is not None and self.experiment.gamma_calorimeter.efficiency is None:
            notes.append(f"The γ calorimeter gives no efficiency: each γ ray fires an element with "
                         f"{100 * CALORIMETER_EFFICIENCY:.0f} %.")
        summary = {
            "number": self.number, "kind": KIND, "status": "stopped" if stopped else "complete",
            "stopped_for": self.stopped_for, "created": self.created, "finished": _now(), "seed": self.seed,
            "duration_s": duration, "real_s": real, "budget_s": self.budget_s,
            "scaled": duration > real * (1 + 1e-9), "generated": self.runner.generated,
            "generator": self.runner.state(), "extended": self.extended, "distance_um": self.distance_um,
            "trigger": dict(self.trigger), "residues": c["residues"], "residues_per_s": self.runner.rate,
            "triggered": trig, "events": int(len(arrays["r_event_step"])),
            "gammas": int(len(arrays["g_row"])), "beta_nominal": self.runner.beta_nominal,
            "channels": list(self.runner.table.labels), "crystals": [cr.name for cr in self.runner.crystals],
            "detectors": ["trigger"], "counts": {"trigger": trig},
            "rates_per_s": {"trigger": trig / real if real else 0.0},
            "counts_in_run": {"trigger": trig * duration / real if real else 0.0},
            "gamma_counts": dict(c["gamma"]), "gamma_singles": dict(c["singles"]), "folds": list(c["folds"]),
            "per_channel": copy.deepcopy(c["channels"]),
            "ancillaries": {"neutrons": list(c["neutrons"]), "charged": list(c["charged"]),
                            "calorimeter": list(c["calorimeter"]), "separator_accepted": c["separator"]},
            "memory_limit_bytes": self.memory_limit, "cpu_s": time.time() - self.started,
            "counters_raw": c, "plan": self.plan, "notes": notes,
        }
        summary.update(self.extra)
        (self.folder / "summary.json").write_text(json.dumps(summary, indent=1, default=_json_default),
                                                  encoding="utf-8")
        return FusionRunData.load(self.folder)


def _add_hist(counters: dict, key: str, values: np.ndarray) -> None:
    h = np.bincount(np.asarray(values, dtype=np.int64), minlength=len(counters[key]))
    old = counters[key] + [0] * (len(h) - len(counters[key]))
    counters[key] = [int(a + b) for a, b in zip(old, h)]


def run_folder(folder: ExperimentFolder, number: int, duration_s: float, distance_um: Optional[float]) -> Path:
    """``runs/NN-fusion-<duration>``, with ``-<d>um`` for a plunger distance."""
    what = f"{KIND}-{duration_label(duration_s)}" + (f"-{distance_um:g}um" if distance_um is not None else "")
    return folder.path / "runs" / f"{number:02d}-{what}"


def start(folder: ExperimentFolder, duration=None, budget=DEFAULT_BUDGET_S, seed: Optional[int] = None,
          experiment: Optional[Experiment] = None, plan: Optional[dict] = None, *,
          distance_um: Optional[float] = None, fold: Optional[int] = None, neutrons_min: Optional[int] = None,
          charged_min: Optional[int] = None, gates: bool = True) -> FusionRunTaker:
    """Begin a fusion run in ``folder``: a :class:`FusionRunTaker` to step, as :func:`physim.nuclear.runs.start`
    gives for a beam run. ``duration`` defaults to the setup's beam time; ``budget`` is the real part (at most
    :data:`~physim.nuclear.runs.MAX_BUDGET_S`); ``seed`` defaults to the run's number. ``distance_um`` is the
    plunger distance (µm, without the offset) of a distance series. The trigger: γ-ray ``fold`` (default
    :data:`DEFAULT_FOLD`), and the ancillaries' gates (see the module); ``neutrons_min`` and ``charged_min``
    ask for detected neutrons and charged particles, ``gates=False`` turns the separator's and calorimeter's
    gates off."""
    exp = experiment if experiment is not None else folder.setup()
    exp.validate()
    if exp.fusion is None:
        raise ValueError("a fusion run needs a fusion-evaporation setup ([reaction] type = \"fusion-evaporation\")")
    if not exp.gamma_detectors:
        raise ValueError("a fusion run needs γ-ray detectors: the trigger counts their γ rays")
    if distance_um is not None:
        distance_um = float(distance_um)
        if not distance_um >= 0:
            raise ValueError("distance_um must not be negative")
        if not _takes_distance():
            raise ValueError("this physim's residue generator cannot place the residues in a plunger yet "
                             "(simulate_residue_gammas takes no distance_um)")
    duration_s = _seconds(duration) if duration is not None else _seconds(exp.run.beam_time)
    if not duration_s > 0:
        raise ValueError("the run needs a duration")
    budget_s = min(_seconds(budget), MAX_BUDGET_S)
    if not budget_s > 0:
        raise ValueError("the run needs a real part (budget) above zero")
    trigger = _trigger(exp, fold, neutrons_min, charged_min, gates)
    number = folder.next_number()
    seed = number if seed is None else int(seed)
    path = run_folder(folder, number, duration_s, distance_um)
    return FusionRunTaker(path, number, exp, duration_s, budget_s, seed, trigger, distance_um, plan=plan)


def extend(run: RunData, more) -> FusionRunTaker:
    """Continue the fusion ``run`` for ``more`` beam time (text with a unit, or seconds) with its own setup, trigger
    and seed stream; the real part grows up to the run's budget, the rest is scaled. Step the result and
    :meth:`FusionRunTaker.write` it, which rewrites the run's folder."""
    if run.kind != KIND:
        raise ValueError(f"not a fusion run (kind {run.kind!r}); extend it with physim.nuclear.runs.extend")
    s = run.summary
    extra = {k: s[k] for k in ("retake_of",) if k in s}
    return FusionRunTaker(run.folder, run.number, run.experiment, run.duration_s + _seconds(more), s["budget_s"],
                          s["seed"], s["trigger"], s.get("distance_um"), extra, s.get("plan"), previous=run)


def estimate(experiment: Experiment, duration=None, budget=DEFAULT_BUDGET_S) -> dict:
    """What a fusion run will take, before it is taken: the real and scaled parts (s), the residues the real part
    makes, and a rough CPU time (s)."""
    from .residues import simulate_residues

    duration_s = _seconds(duration) if duration is not None else _seconds(experiment.run.beam_time)
    real = min(duration_s, min(_seconds(budget), MAX_BUDGET_S))
    rate = float(simulate_residues(experiment, 1, 1).expected())
    return {"duration_s": duration_s, "real_s": real, "scaled_s": max(duration_s - real, 0.0),
            "residues_per_s": rate, "residues": rate * real, "cpu_s": 1.0 + rate * real * CPU_S_PER_EVENT}


# -- a stored fusion run ------------------------------------------------------------------------------------------

@dataclass
class FusionGammas:
    """The counted γ rays of a fusion run's recorded events: columns ``row`` (the event's row), ``event``,
    ``channel``, ``crystal``, ``initial``, ``final``, ``energy0``, ``deposited``, ``measured`` and ``corrected``
    (MeV), ``full``, ``theta``, ``phi`` (degrees) and ``weight`` (1 / real seconds: each counts once)."""

    columns: dict
    crystals: list
    channels: list
    beta_nominal: float
    real_s: float

    def __len__(self) -> int:
        return len(self.columns["row"])

    def __getitem__(self, name: str) -> np.ndarray:
        return self.columns[name]

    def select(self, crystal: Optional[str] = None, channel: Optional[str] = None) -> np.ndarray:
        m = np.ones(len(self), dtype=bool)
        if crystal is not None:
            idx = [i for i, c in enumerate(self.crystals) if c == crystal or c.rsplit(" ", 1)[0] == crystal]
            if not idx:
                raise KeyError(f"no crystal {crystal!r}; the crystals are {', '.join(self.crystals)}")
            m &= np.isin(self.columns["crystal"], idx)
        if channel is not None:
            m &= self.columns["channel"] == self.channels.index(channel)
        return m

    def spectrum(self, crystal: Optional[str] = None, corrected: bool = True, channel: Optional[str] = None,
                 bins: int = 2000, range_kev: tuple = (0.0, 2000.0)) -> tuple:
        """Counts per bin in the real part, and the edges (keV): Doppler-corrected energies, or measured ones."""
        m = self.select(crystal, channel)
        x = self.columns["corrected" if corrected else "measured"][m] * 1e3
        return np.histogram(x, bins=bins, range=range_kev)

    def pairs(self, corrected: bool = True) -> tuple:
        """Every ordered pair of distinct γ rays of one event: (first energy, second energy), keV."""
        rows = self.columns["row"]
        e = self.columns["corrected" if corrected else "measured"] * 1e3
        order = np.argsort(rows, kind="stable")
        rows, e = rows[order], e[order]
        first, second = [], []
        if len(rows) > 1:
            longest = int(np.bincount(rows).max())
            for s in range(1, longest):
                same = np.flatnonzero(rows[s:] == rows[:-s])
                first += [e[same], e[same + s]]
                second += [e[same + s], e[same]]
        if not first:
            return np.zeros(0), np.zeros(0)
        return np.concatenate(first), np.concatenate(second)

    def matrix(self, bins: int = 1000, range_kev: tuple = (0.0, 2000.0), corrected: bool = True) -> tuple:
        """The γγ matrix of the real part, symmetric (each pair counted both ways), and its edges (keV)."""
        a, b = self.pairs(corrected)
        counts, edges, _ = np.histogram2d(a, b, bins=bins, range=[range_kev, range_kev])
        return counts, edges

    def gated(self, energy_kev: float, half_width_kev: float, bins: int = 2000,
              range_kev: tuple = (0.0, 2000.0), corrected: bool = True) -> tuple:
        """The γ rays in coincidence with one inside the gate energy ± half width (keV): counts and edges (keV).
        The gating γ ray itself is not in the spectrum."""
        a, b = self.pairs(corrected)
        inside = np.abs(a - energy_kev) <= half_width_kev
        return np.histogram(b[inside], bins=bins, range=range_kev)


class FusionRunData(RunData):
    """A stored fusion run: the :class:`~physim.nuclear.runs.RunData` of its folder, with the residue events and
    their γ rays (:meth:`events`, :meth:`gammas`) and the γγ coincidences (:meth:`matrix`, :meth:`gated`)."""

    @classmethod
    def load(cls, folder: Union[str, Path]) -> FusionRunData:
        run = RunData.load(folder)
        if run.kind != KIND:
            raise ValueError(f"{folder} holds a {run.kind} run, not a fusion run")
        return cls(run.folder, run.summary, run.setup)

    @classmethod
    def of(cls, run: RunData) -> FusionRunData:
        """The fusion run of a :class:`~physim.nuclear.runs.RunData` loaded by the generic loader."""
        return run if isinstance(run, cls) else cls.load(run.folder)

    @property
    def distance_um(self) -> Optional[float]:
        return self.summary.get("distance_um")

    @property
    def one_nucleus(self) -> bool:
        return False

    @property
    def channels(self) -> list:
        return list(self.summary["channels"])

    def describe(self) -> str:
        s = self.summary
        tr = s["trigger"]
        line = f"Fusion run, {duration_label(self.duration_s)}, γ-ray fold ≥ {tr['fold']}"
        if self.distance_um is not None:
            line += f", plunger at {self.distance_um:g} µm"
        line += (f"; the first {duration_label(self.real_s)} simulated event by event, the rest scaled"
                 if self.scale > 1.000001 else "; every event simulated")
        if s.get("status") == "stopped":
            line += " (stopped)"
        return line

    def counts(self, detector: Optional[str] = None) -> float:
        """Events the trigger kept in the real part."""
        return float(self.summary.get("triggered", 0))

    def events(self) -> dict:
        """The recorded events of the real part, one row each: ``event``, ``channel``, ``fold``, ``neutrons``,
        ``charged``, ``calorimeter``, ``separator``, ``beta``, ``theta``, ``phi``, ``depth``, ``entry_spin`` and
        ``weight`` (1 / real seconds)."""
        if self._events is None:
            a = self._arrays(lambda k: k in EVENT_COLUMNS)
            cols = {k[2:]: a[k] for k in EVENT_COLUMNS}
            cols["event"] = np.cumsum(cols.pop("event_step").astype(np.int64))
            for k in ("channel", "fold", "neutrons", "charged", "calorimeter"):
                cols[k] = cols[k].astype(np.int64)
            cols["separator"] = cols["separator"].astype(bool)
            for k in ("beta", "theta", "phi", "depth", "entry_spin"):
                cols[k] = cols[k].astype(float)
            cols["weight"] = np.full(len(cols["event"]), 1.0 / self.real_s if self.real_s > 0 else 0.0)
            self._events = cols
        return self._events

    def gammas(self, plain: bool = False) -> FusionGammas:
        """The counted γ rays of the recorded events (:class:`FusionGammas`), Doppler corrected for β_CN along the
        beam. ``plain`` is accepted for the particle runs' signature: fusion runs have no add-back."""
        if "g" not in self._gammas:
            from .gamma_events import _unit, crystals_of
            from .residues import doppler_correct

            a = self._arrays(lambda k: k.startswith("g_"))
            cols = {k[2:]: a[k] for k in a}
            rows = cols["row"].astype(np.int64)
            ev = self.events()
            out = {"row": rows, "event": ev["event"][rows], "channel": ev["channel"][rows],
                   "crystal": cols["crystal"].astype(np.int64), "initial": cols["initial"].astype(np.int64),
                   "final": cols["final"].astype(np.int64), "energy0": cols["energy0"].astype(float),
                   "deposited": cols["deposited"].astype(float), "measured": cols["measured"].astype(float),
                   "full": cols["full"].astype(bool), "theta": cols["theta"].astype(float),
                   "phi": cols["phi"].astype(float)}
            if "stopped" in cols:
                out["stopped"] = cols["stopped"].astype(bool)
            crystals = crystals_of(self.experiment, plain=True)[0]
            normals = _unit(np.array([c.centre for c in crystals]))
            beta = float(self.summary["beta_nominal"])
            out["corrected"] = (doppler_correct(beta, normals[out["crystal"]], out["measured"]) if len(rows)
                                else np.zeros(0))
            out["weight"] = np.full(len(rows), 1.0 / self.real_s if self.real_s > 0 else 0.0)
            self._gammas["g"] = FusionGammas(out, list(self.summary["crystals"]), self.channels, beta, self.real_s)
        return self._gammas["g"]

    def matrix(self, bins: int = 1000, range_kev: tuple = (0.0, 2000.0), corrected: bool = True) -> tuple:
        """The γγ matrix of the real part (:meth:`FusionGammas.matrix`)."""
        return self.gammas().matrix(bins, range_kev, corrected)

    def gated(self, energy_kev: float, half_width_kev: float, bins: int = 2000,
              range_kev: tuple = (0.0, 2000.0), corrected: bool = True) -> tuple:
        """The spectrum gated on a γ-ray energy (:meth:`FusionGammas.gated`)."""
        return self.gammas().gated(energy_kev, half_width_kev, bins, range_kev, corrected)


def load(folder: Union[str, Path, RunData]) -> FusionRunData:
    """A stored fusion run, from its folder or from the generic :class:`~physim.nuclear.runs.RunData`."""
    return FusionRunData.of(folder) if isinstance(folder, RunData) else FusionRunData.load(folder)


__all__ = ["BLOCK_EVENTS", "DEFAULT_FOLD", "KIND", "FusionGammas", "FusionRunData", "FusionRunTaker", "FusionRunner",
           "estimate", "extend", "load", "run_folder", "start"]
