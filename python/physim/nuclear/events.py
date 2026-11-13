"""Monte Carlo events: simulated spectra for an experiment before it gets beam.

Each event is one reaction in the target (or its backing). The generator samples where in the target it
happens, slows the beam to that depth with straggling, draws the scattering angle from the cross section,
works out the kinematics, carries the ejectile and the recoil out of the target to the first detector each one
crosses, and applies that detector's dead layer, punch-through, resolution and threshold. The per-event loop runs
in physim's Rust engine on all cores; the physics tables come from :mod:`physim.nuclear.stopping`::

    from physim.nuclear import Experiment
    from physim.nuclear.events import simulate

    ev = simulate(Experiment.example("alpha_on_gold"), events=1_000_000, seed=1)
    ev.rate("A45")                       # counts/s, and its statistical error
    counts, edges = ev.spectrum("A45")   # measured-energy spectrum, counts in the planned beam time
    ev.columns["theta"]                  # per-particle arrays: angles, energies, depth, detector, strip, ...

Events are weighted: each stands for a rate (per second), so rates and spectra are absolute. The scattering angle
is drawn only over the ranges the detectors can see, which wastes few events on particles that miss. The same
seed gives the same events on any machine and with any number of threads.

**How the events are shared out.** A channel's detectors may see separate ranges of the CM angle
(:func:`~physim.nuclear.rates.cm_ranges`): backscattered beam in one, and in another the slow recoils of distant
collisions that just clear a detector's edge near 90°. Rutherford's cross section in the second can be hundreds of
times that in the first, while none of its recoils leaves the target with enough energy to be counted. So when
some channel has more than one range, a short pilot run (:data:`PILOT_PER_RANGE` events per range, always the same
seed) finds the fraction of each range's particles that are counted, and the events are shared out by the rate
*counted* from each range rather than the rate into it (:func:`share_out`). A range in which nothing is counted
keeps a small share (:data:`FLOOR_SHARE`), so its rate is still there, with the larger statistical error its few
events give; no count is dropped, whatever the shares.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .._core import nuclear_events
from . import data
from .detectors import Array
from .rates import (
    MAX_PATH_FACTOR,
    Channel,
    _q,
    beam_ion,
    beam_energy_at,
    channels,
    cm_ranges,
    energy_cut,
    energy_sigma,
    excitation_probability,
    excited_partner,
    spot_sigma_mm,
    stack,
    stopping,
    tilt_deg,
)
from .rutherford import E2_MEV_FM

#: The particle-level columns of :class:`Events`, with units.
COLUMNS = {
    "event": "event number",
    "detector": "index into Events.detectors",
    "segment_i": "strip along u (rectangles) or ring (annular detectors)",
    "segment_j": "strip along v (rectangles) or sector (annular detectors)",
    "recoil": "True for the recoiling target nucleus, False for the scattered beam (ejectile)",
    "channel": "index into Events.channels",
    "depth": "depth of the reaction below the target's front face, along its normal, mg/cm2",
    "beam_energy": "beam energy at the reaction, MeV",
    "energy": "energy of the particle as it leaves the reaction, MeV",
    "energy_face": "energy reaching the detector face, MeV",
    "deposited": "energy deposited in the active volume, MeV",
    "measured": "measured energy (deposited, with resolution), MeV",
    "counted": "reached the detector and was measured above threshold",
    "theta": "lab polar angle of the track, degrees",
    "phi": "lab azimuth of the track, degrees",
    "theta_cm": "CM angle of the ejectile, degrees",
    "weight": "rate this particle stands for, per second",
    "x": "where on the target the reaction happened, across the beam (the spot or the halo), mm",
    "y": "the same, vertically, mm",
}


#: Events of the pilot run per sampled range (see :func:`share_out`), and its seed: always the same, so that an
#: experiment has one generator input whatever seed its events are drawn with.
PILOT_PER_RANGE = 4000
PILOT_SEED = 0x51A7ED
#: Share of a channel's events that its ranges keep whatever the pilot run finds in them, and the least a channel
#: is picked relative to the average channel.
FLOOR_SHARE = 0.05


def share_out(config: dict, counted: Optional[list] = None) -> None:
    """Sets how often each sampled range of ``config`` is picked (its ``probability``; any choice is unbiased, as
    the weights carry it).

    A channel is picked in proportion to its Rutherford rate, and its ranges half in proportion to their
    Rutherford rates and half to their widths in ln sin²(θ*/2), as the angle is drawn inside a range. With
    ``counted`` (the fraction of each range's particles that a detector counts, from a pilot run) every rate and
    width is first multiplied by that fraction, so the events go where particles are counted;
    :data:`FLOOR_SHARE` of a channel's events stay spread over its ranges by width, and no channel is picked less
    than :data:`FLOOR_SHARE` of the average."""
    chans = config["channels"]
    owner = np.asarray(config["channel_of"])
    r = np.array([1 / c["u_min"] - 1 / c["u_max"] for c in chans])
    width = np.array([math.log(c["u_max"] / c["u_min"]) for c in chans])
    amplitude = np.array([c["amplitude"] for c in chans])
    f = np.ones(len(chans)) if counted is None else np.asarray(counted, dtype=float)
    if not np.any(amplitude * r * f > 0):  # nothing counted anywhere: as without a pilot run
        counted, f = None, np.ones(len(chans))
    mass, share = np.zeros(len(chans)), np.ones(len(chans))
    for k in np.unique(owner):
        m = owner == k
        mass[m] = (amplitude[m] * r[m] * f[m]).sum()
        if m.sum() > 1:
            rf, wf = r[m] * f[m], width[m] * f[m]
            main = 0.5 * rf / rf.sum() + 0.5 * wf / wf.sum() if rf.sum() > 0 else width[m] / width[m].sum()
            share[m] = (1 - FLOOR_SHARE) * main + FLOOR_SHARE * width[m] / width[m].sum()
    if counted is not None:
        per_channel = np.array([mass[owner == k][0] for k in np.unique(owner)])
        mass = np.maximum(mass, FLOOR_SHARE * per_channel.mean())
    for c, p in zip(chans, mass * share):
        c["probability"] = max(float(p), 1e-300)


def _counted_fractions(config: dict) -> list:
    """The fraction of the particles reaching a face from each sampled range that the face counts, from a pilot
    run that picks every range equally often; 1 for a range that sent none anywhere."""
    chans = config["channels"]
    pilot = dict(config, channels=[dict(c, probability=1.0) for c in chans])
    n = PILOT_PER_RANGE * len(chans)
    cols = nuclear_events(pilot, n, PILOT_SEED, n, 0)
    reached = np.bincount(np.asarray(cols["channel"]), minlength=len(chans))
    counted = np.bincount(np.asarray(cols["channel"]), weights=np.asarray(cols["counted"], dtype=float),
                          minlength=len(chans))
    return np.where(reached > 0, counted / np.maximum(reached, 1), 1.0).tolist()


def channel_index(config: dict, column) -> np.ndarray:
    """The ``channel`` column of :func:`physim._core.nuclear_events`, which numbers the sampled ranges of
    ``config``, as indices into the channels :func:`generator_config` returned."""
    column = np.asarray(column)
    return np.asarray(config["channel_of"], dtype=column.dtype)[column]


def generator_config(experiment, theta_floor: float = 0.5, pilot: bool = True) -> tuple:
    """The event generator's input for an experiment, and the reaction channels it covers.

    Returns (config dict for :func:`physim._core.nuclear_events`, list of :class:`~physim.nuclear.rates.Channel`).
    ``config["channels"]`` has one entry per range of CM angle a channel is sampled over
    (:func:`~physim.nuclear.rates.cm_ranges`), and ``config["channel_of"]`` the channel each belongs to
    (:func:`channel_index` turns the generator's ``channel`` column into channel numbers). With ``pilot``, and a
    channel that has more than one range, a pilot run sets how the events are shared out (:func:`share_out`).
    """
    layers = stack(experiment)
    array = Array.from_experiment(experiment)
    beam = beam_ion(experiment)
    species = [beam]
    for ch in channels(layers, experiment):
        if ch.nuclide.name not in species:
            species.append(ch.nuclide.name)
    materials: list = []

    def material_index(mat) -> int:
        for i, m in enumerate(materials):
            if m == mat:
                return i
        materials.append(mat)
        return len(materials) - 1

    layer_cfg = [{"thickness": lay.thickness, "material": material_index(lay.material)} for lay in layers]
    faces = []
    for g, d in zip(array, experiment.detectors):
        mat = d.material_data()
        areal = lambda x: mat.areal_density_mg_cm2(_q(x)) if x is not None else 0.0  # noqa: E731
        face = {"centre": g.centre.tolist(), "n": g.n.tolist(), "u": g.u.tolist(), "v": g.v.tolist(),
                "shape": g.shape, "material": material_index(mat), "dead_layer": areal(d.dead_layer),
                "thickness": areal(d.thickness),
                "fwhm": _q(d.resolution).to("MeV") if d.resolution is not None else 0.0,
                "threshold": _q(d.threshold).to("MeV") if d.threshold is not None else 0.0}
        if g.shape == "rectangle":
            face.update(width=g.width, height=g.height, strips_x=g.strips_x, strips_y=g.strips_y)
        else:
            face.update(inner_radius=g.inner_radius or 0.0, outer_radius=g.outer_radius, rings=g.rings,
                        sectors=g.sectors)
        faces.append(face)
    # Dead material stops particles like a detector that never counts: a face with a threshold out of reach.
    dead = material_index(data.material("Si"))
    for b in array.blockers:
        face = {"centre": b.centre.tolist(), "n": b.n.tolist(), "u": b.u.tolist(), "v": b.v.tolist(),
                "shape": b.shape, "material": dead, "dead_layer": 0.0, "thickness": 1e6, "fwhm": 0.0,
                "threshold": 1e30}
        if b.shape == "rectangle":
            face.update(width=b.width, height=b.height, strips_x=1, strips_y=1)
        else:
            face.update(inner_radius=b.inner_radius or 0.0, outer_radius=b.outer_radius, rings=1, sectors=1)
        faces.append(face)

    z1 = experiment.beam.Z
    pps = experiment.beam.particles_per_second
    chan_cfg, used, channel_of = [], [], []
    for ch in channels(layers, experiment):
        ranges = cm_ranges(experiment, layers, ch, array, theta_floor)
        if not ranges:
            continue
        k = z1 * ch.nuclide.Z * E2_MEV_FM
        # Rate at the middle of the layer: sets how often the channel is picked (any choice is unbiased). The
        # rate into a range is ``amplitude`` × (1/u_min − 1/u_max).
        e_mid = float(beam_energy_at(experiment, ch.layer, [layers[ch.layer].thickness / 2], layers)[0])
        m1, m2 = data.nuclide(beam).atomic_mass_u, ch.nuclide.atomic_mass_u
        e_cm = max(e_mid * m2 / (m1 + m2), 1e-6)
        cfg = {"layer": ch.layer, "target": species.index(ch.nuclide.name), "atoms_per_cm2": ch.atoms_per_cm2,
               "k": k, "amplitude": pps * ch.atoms_per_cm2 * 4 * math.pi * (k / (4 * e_cm)) ** 2 * 1e-26}
        partner = excited_partner(experiment, ch, layers)
        if partner is not None:
            # One particle per collision: elastic scattering is what the excitation on this nuclide leaves,
            # 1 − P at each angle, so that elastic and excited together are Rutherford scattering.
            p = np.asarray(excitation_probability(experiment, partner, layers).probability(
                np.linspace(0.0, 180.0, 721)), dtype=float)
            cfg.update(p_table=(1.0 - np.clip(p, 0.0, 1.0)).tolist())
        if ch.excitation is not None:
            # Excitation probability every 0.25°, at mid-layer, the one the analytic rates use: of the
            # [reaction] state in first order, or with a level scheme that decays by cascades, of every level
            # its solution populates (to all orders unless the setup says first; backlog item 71).
            p = np.asarray(excitation_probability(experiment, ch, layers).probability(np.linspace(0.0, 180.0, 721)))
            # Picked as often as elastic scattering on the same nuclei (the weights carry P), so the rare
            # inelastic events get as many samples as the elastic ones.
            cfg.update(excitation=ch.excitation.energy_mev, excite_recoil=ch.excitation.excite == "target",
                       p_table=p.tolist())
        for lo, hi in ranges:
            u_min, u_max = (math.sin(math.radians(a) / 2) ** 2 for a in (lo, hi))
            chan_cfg.append(dict(cfg, u_min=u_min, u_max=u_max))
            channel_of.append(len(used))
        used.append(ch)
    if not chan_cfg:
        raise ValueError("no detector can see any scattered particle or recoil: nothing to simulate")
    tables = [stopping(s, m).transport_table() for s in species for m in materials]
    config = {
        "masses": [data.nuclide(s).nuclear_mass_mev for s in species],
        "beam": 0,
        "beam_energy": experiment.beam.energy_mev,
        "energy_sigma": energy_sigma(experiment),
        "spot_sigma": spot_sigma_mm(experiment),
        "particles_per_second": pps,
        "tilt": math.radians(tilt_deg(experiment)),
        "layers": layer_cfg,
        "channels": chan_cfg,
        "channel_of": channel_of,
        "materials": len(materials),
        "tables": tables,
        "faces": faces,
        "max_path_factor": MAX_PATH_FACTOR,
        "species": species,
        # The beam halo (backlog item 74): a share of the beam spread evenly over the frame's aperture.
        "halo_fraction": (_q(experiment.beam.halo_fraction).to("%") / 100
                          if experiment.beam.halo_fraction is not None else 0.0),
        "halo_radius": _q(experiment.beam.halo_radius).to("mm") if experiment.beam.halo_radius is not None else 0.0,
    }
    share_out(config)
    if pilot and len(chan_cfg) > len(used):
        share_out(config, _counted_fractions(config))
    return config, used


@dataclass
class Events:
    """Simulated particles that reached a detector, as columns (see :data:`COLUMNS`), with what they stand for."""

    columns: dict
    #: Events generated (reactions, including those whose particles missed every detector).
    n_events: int
    seed: int
    #: Detector names, indexed by the ``detector`` column.
    detectors: list
    #: Channel labels such as ``"197Au (target)"``, indexed by the ``channel`` column.
    channels: list
    #: Planned beam time, s: counts in the run are rate × beam time.
    beam_time_s: float
    #: Particles leaving the reaction below this energy (MeV) are outside the sampled range; :meth:`rate` with
    #: ``counted=False`` leaves them out, to compare with :class:`~physim.nuclear.rates.Rates`.
    min_energy: float

    def __len__(self) -> int:
        return len(self.columns["event"])

    def __getitem__(self, name: str) -> np.ndarray:
        return self.columns[name]

    def select(self, detector: Optional[str] = None, segment: Optional[tuple] = None, counted: bool = True,
               particle: Optional[str] = None, channel: Optional[str] = None) -> np.ndarray:
        """Boolean mask over the particles."""
        c = self.columns
        mask = np.ones(len(self), dtype=bool)
        if detector is not None:
            mask &= c["detector"] == self.detectors.index(detector)
        if segment is not None:
            mask &= (c["segment_i"] == segment[0]) & (c["segment_j"] == segment[1])
        if counted:
            mask &= c["counted"]
        else:
            mask &= c["energy"] >= self.min_energy
        if particle is not None:
            if particle not in ("ejectile", "recoil"):
                raise ValueError("particle must be 'ejectile' or 'recoil'")
            mask &= c["recoil"] == (particle == "recoil")
        if channel is not None:
            mask &= c["channel"] == self.channels.index(channel)
        return mask

    def rate(self, detector: Optional[str] = None, segment: Optional[tuple] = None, counted: bool = True,
             **kw) -> tuple:
        """(counts per second, statistical standard error). ``counted=False`` counts every particle reaching the
        face that left the reaction with at least :attr:`min_energy`, as :class:`~physim.nuclear.rates.Rates`
        does."""
        w = self.columns["weight"][self.select(detector, segment, counted, **kw)]
        return float(w.sum()), float(math.sqrt(np.sum(w**2)))

    def counts(self, detector: Optional[str] = None, segment: Optional[tuple] = None, **kw) -> float:
        """Counts expected in the planned beam time."""
        return self.rate(detector, segment, **kw)[0] * self.beam_time_s

    def spectrum(self, detector: Optional[str] = None, segment: Optional[tuple] = None, bins=200, range=None,
                 quantity: str = "measured", **kw) -> tuple:
        """Histogram of ``quantity`` (a column, by default the measured energy): (counts in the planned beam
        time per bin, bin edges)."""
        m = self.select(detector, segment, **kw)
        x = self.columns[quantity][m]
        if range is None and len(x):
            range = (float(x.min()), float(x.max()) + 1e-9)
        h, edges = np.histogram(x, bins=bins, range=range, weights=self.columns["weight"][m] * self.beam_time_s)
        return h, edges


def simulate(experiment, events: int = 1_000_000, seed: int = 1, theta_floor: float = 0.5) -> Events:
    """Generate ``events`` reactions for ``experiment`` with ``seed`` and return the particles that reached a
    detector."""
    if events < 1:
        raise ValueError("events must be at least 1")
    config, used = generator_config(experiment, theta_floor)
    cols = dict(nuclear_events(config, int(events), int(seed)))
    cols["channel"] = channel_index(config, cols["channel"])
    names = [d.name for d in Array.from_experiment(experiment)]
    reached = np.asarray(cols["detector"]) < len(names)  # particles stopped by dead material are left out
    if not reached.all():
        cols = {k: np.asarray(v)[reached] for k, v in cols.items()}
    return Events(cols, int(events), int(seed), names,
                  [ch.label for ch in used], _q(experiment.run.beam_time).to("s"), energy_cut(experiment))


__all__ = ["COLUMNS", "Channel", "Events", "channel_index", "generator_config", "share_out", "simulate"]
