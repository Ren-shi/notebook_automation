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
is drawn only over the range the detectors can see, which wastes few events on particles that miss. The same
seed gives the same events on any machine and with any number of threads.
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
    cm_acceptance,
    coulex_for,
    energy_cut,
    energy_sigma,
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
}


def generator_config(experiment, theta_floor: float = 0.5) -> tuple:
    """The event generator's input for an experiment, and the reaction channels it covers.

    Returns (config dict for :func:`physim._core.nuclear_events`, list of :class:`~physim.nuclear.rates.Channel`).
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

    z1 = experiment.beam.Z
    pps = experiment.beam.particles_per_second
    chan_cfg, used = [], []
    for ch in channels(layers, experiment):
        acc = cm_acceptance(experiment, layers, ch, array, theta_floor)
        if acc is None:
            continue
        u_min, u_max = (math.sin(math.radians(a) / 2) ** 2 for a in acc)
        k = z1 * ch.nuclide.Z * E2_MEV_FM
        # Rate at the middle of the layer: sets how often the channel is picked (any choice is unbiased).
        e_mid = float(beam_energy_at(experiment, ch.layer, [layers[ch.layer].thickness / 2], layers)[0])
        m1, m2 = data.nuclide(beam).atomic_mass_u, ch.nuclide.atomic_mass_u
        e_cm = max(e_mid * m2 / (m1 + m2), 1e-6)
        rate = pps * ch.atoms_per_cm2 * 4 * math.pi * (k / (4 * e_cm)) ** 2 * (1 / u_min - 1 / u_max) * 1e-26
        cfg = {"layer": ch.layer, "target": species.index(ch.nuclide.name), "atoms_per_cm2": ch.atoms_per_cm2,
               "k": k, "u_min": u_min, "u_max": u_max, "probability": max(rate, 1e-300)}
        if ch.excitation is not None:
            # Excitation probability every 0.25° (from the Coulex at mid-layer, as for the analytic rates).
            p = np.asarray(coulex_for(experiment, ch, layers).probability(np.linspace(0.0, 180.0, 721)))
            # Picked as often as elastic scattering on the same nuclei (the weights carry P), so the rare
            # inelastic events get as many samples as the elastic ones.
            cfg.update(excitation=ch.excitation.energy_mev, excite_recoil=ch.excitation.excite == "target",
                       p_table=p.tolist())
        chan_cfg.append(cfg)
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
        "materials": len(materials),
        "tables": tables,
        "faces": faces,
        "max_path_factor": MAX_PATH_FACTOR,
        "species": species,
    }
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
    cols = nuclear_events(config, int(events), int(seed))
    return Events(dict(cols), int(events), int(seed), [d.name for d in Array.from_experiment(experiment)],
                  [ch.label for ch in used], _q(experiment.run.beam_time).to("s"), energy_cut(experiment))


__all__ = ["COLUMNS", "Channel", "Events", "generator_config", "simulate"]
