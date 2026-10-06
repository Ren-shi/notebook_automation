"""Particle and γ-ray tracks of simulated events, for the scene: a sample of events as straight paths from the
beam line to the target and from the target to the segment or crystal that recorded each hit.

::

    from physim.nuclear.tracks import sample_tracks

    tracks = sample_tracks(experiment, gammas, n=40, select="coincidences")
    for t in tracks:
        t.paths        # [{"what": "beam" | "ejectile" | "recoil" | "gamma", "points": [[x, y, z], ...], "hit": ...}]
        t.event        # the event's numbers

**How the sample is chosen.** The generator draws scattering angles evenly over the detectors' range, so its
events over-represent the rare large-angle scatterings; each event carries a weight, the rate it stands for.
``weighted=True`` (the default) picks events with probability proportional to that weight, so the sample looks
like a run; ``weighted=False`` picks them as generated, and the rare events show. Either way the sample's
description says which (:func:`describe`).

A particle's path is drawn straight from the target to where its track crosses the detector's face (the segment
in the event record is the one that contains that point); an undetected partner is drawn to the edge of the
scene. γ rays start where the nucleus leaves the target and run to the crystal they hit, or to the edge. Paths
inside detectors are not drawn.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .detectors import Array
from .kinematics import TwoBody
from .rates import beam_ion, stack
from . import data


@dataclass
class Track:
    """One event's paths and numbers."""

    event: int
    paths: list
    #: The particle rows and γ rows of the event, as dicts of the columns.
    particles: list
    gammas: list = field(default_factory=list)
    weight: float = 0.0
    coincidence: bool = False
    channel: str = ""
    #: What the event record names: [(detector, segment), ...] and [crystal, ...].
    hits: list = field(default_factory=list)
    crystal_hits: list = field(default_factory=list)


def _direction(theta_deg: float, phi_deg: float) -> np.ndarray:
    th, ph = math.radians(theta_deg), math.radians(phi_deg)
    return np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])


def _crystal_centres(experiment) -> list:
    out = []
    for i, gd in enumerate(experiment.gamma_detectors):
        for label, centre, radius in gd.elements():
            name = (gd.name or f"G{i + 1}") + (f" {label}" if label else "")
            out.append((name, np.array(centre, dtype=float), radius))
    return out


def sample_tracks(experiment, gammas=None, events=None, n: int = 30, select: str = "all", weighted: bool = True,
                  seed: int = 1, extent_mm: Optional[float] = None) -> list:
    """A sample of ``n`` events as :class:`Track` s. ``gammas`` is a
    :class:`~physim.nuclear.gamma_events.GammaEvents` (its particle events are used), or ``events`` the particle
    events alone. ``select`` is "all", "coincidences" (a particle hit and a γ hit) or a channel label."""
    ev = gammas.events if gammas is not None else events
    if ev is None:
        raise ValueError("give the γ events or the particle events")
    c = ev.columns
    rng = np.random.default_rng(seed)
    by_event: dict = {}
    for row in range(len(ev)):
        by_event.setdefault(int(c["event"][row]), []).append(row)
    gamma_rows: dict = {}
    if gammas is not None:
        for g in range(len(gammas)):
            if gammas["counted"][g]:
                gamma_rows.setdefault(int(gammas["event"][g]), []).append(g)
    keys = sorted(by_event)
    if select == "coincidences":
        keys = [k for k in keys if k in gamma_rows]
    elif select != "all":
        if select not in ev.channels:
            raise ValueError(f"no channel {select!r}; the channels are {', '.join(ev.channels)}")
        idx = ev.channels.index(select)
        keys = [k for k in keys if any(c["channel"][r] == idx for r in by_event[k])]
    if not keys:
        return []
    weights = np.array([max(sum(float(c["weight"][r]) for r in by_event[k]), 1e-300) for k in keys])
    p = weights / weights.sum() if weighted else None
    chosen = rng.choice(len(keys), size=min(n, len(keys)), replace=False, p=p)
    array = Array.from_experiment(experiment)
    extent = extent_mm or 1.3 * max([float(np.linalg.norm(g.centre)) for g in array]
                                    + [float(np.linalg.norm(cc[1])) for cc in _crystal_centres(experiment)] + [100.0])
    crystals = _crystal_centres(experiment)
    layers = stack(experiment)
    ion = beam_ion(experiment)
    target = data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
    tracks = []
    for i in chosen:
        k = keys[i]
        rows = by_event[k]
        particles = [{name: (col[r].item() if hasattr(col[r], "item") else col[r]) for name, col in c.items()}
                     for r in rows]
        paths = [{"what": "beam", "points": [[0.0, 0.0, -extent], [0.0, 0.0, 0.0]], "hit": None}]
        hits = []
        seen = set()
        for pr in particles:
            d = _direction(pr["theta"], pr["phi"])
            g = array.geometries[pr["detector"]]
            h = g.hit(d)
            end = (d * h.distance).tolist() if h is not None else (d * extent).tolist()
            what = "recoil" if pr["recoil"] else "ejectile"
            seen.add(what)
            paths.append({"what": what, "points": [[0.0, 0.0, 0.0], end],
                          "hit": {"detector": g.name, "segment": [int(pr["segment_i"]), int(pr["segment_j"])],
                                  "index": int(pr["detector"])}})
            hits.append((g.name, (int(pr["segment_i"]), int(pr["segment_j"]))))
        # The partner that was not detected, from two-body kinematics at the event's beam energy.
        first = particles[0]
        for what in ("ejectile", "recoil"):
            if what in seen:
                continue
            try:
                tb = TwoBody(ion, target, float(first["beam_energy"]))
                point = tb.at_cm(float(first["theta_cm"]), what)
                th = float(point.theta_lab)
                ph = first["phi"] + (180.0 if (what == "recoil") != bool(first["recoil"]) else 0.0)
                d = _direction(th, ph)
                h = array.first_hit(d)
                end = (d * h.distance).tolist() if h is not None else (d * extent).tolist()
                paths.append({"what": what, "points": [[0.0, 0.0, 0.0], end], "hit": None})
            except (ValueError, AttributeError):
                pass
        gam = []
        crystal_hits = []
        for g in gamma_rows.get(k, []):
            row = {name: (col[g].item() if hasattr(col[g], "item") else col[g]) for name, col in gammas.columns.items()}
            d = _direction(row["theta"], row["phi"])
            name, centre, radius = crystals[int(row["crystal"])]
            u = centre / np.linalg.norm(centre)
            t = float(centre @ u) / max(float(d @ u), 1e-9)
            paths.append({"what": "gamma", "points": [[0.0, 0.0, 0.0], (d * t).tolist()],
                          "hit": {"crystal": name, "index": int(row["crystal"])}})
            crystal_hits.append(name)
            gam.append(row)
        tracks.append(Track(k, paths, particles, gam, float(weights[i]), bool(gam), ev.channels[int(first["channel"])],
                            hits, crystal_hits))
    return tracks


def describe(tracks: list, select: str, weighted: bool, n_events: int) -> str:
    """How the sample was chosen, in words, for the scene to state."""
    if not tracks:
        return "No event matches the filter."
    which = {"all": "all events", "coincidences": "events with a particle hit and a γ hit"}.get(select,
                                                                                               f"events of {select}")
    how = ("picked with probability proportional to the rate each stands for, so the sample looks like a run"
           if weighted else "picked as the generator made them: it draws angles evenly, so rare large-angle "
                            "scatterings are over-represented")
    return (f"{len(tracks)} of {n_events} generated reactions, {which}, {how}. Each track stands for a rate of its "
            f"own, from {min(t.weight for t in tracks):.2g} to {max(t.weight for t in tracks):.2g} per second.")


__all__ = ["Track", "describe", "sample_tracks"]
