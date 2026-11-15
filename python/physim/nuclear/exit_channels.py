"""Exit channels from physim's own fusion and evaporation codes (backlog item 75, the wiring of Stages 2 and 3).

Stage 1 (:mod:`physim.nuclear.residues`) takes each exit channel's cross section from a table. This module makes
that table from the app's own physics instead::

    from physim.nuclear import Experiment
    from physim.nuclear.exit_channels import physim_channels

    exp = Experiment.example("fusion_dy158")
    pc = physim_channels(exp, events=20_000, seed=1)
    pc.table                 # a ChannelTable, the same object the setup's [[reaction.channels]] fill
    pc.table.to_setup()      # the table written as [[reaction.channels]] sections
    pc.sigma_fusion_mb       # σ_fus averaged over the target, mb
    pc.warnings              # where the numbers stop being reliable

**How.** :class:`physim.nuclear.fusion.Fusion` gives the fusion cross section and its partial waves σ_ℓ averaged
over the beam's path through the target (Bass barrier, Wong's parabolic transmission). The compound nucleus is made
at the σ-weighted mean excitation energy through the target, and :func:`physim.nuclear.evaporation.evaporate`
follows one statistical decay chain per event from that E* and the σ_ℓ distribution. Each exit channel's cross
section is σ_fus times its share of the events; its entry state is the mean (E*, J) its events end in. Channels
with less than ``min_share`` of the events, events that fission and events that emit no particle are left out of
the table (Stage 1 needs at least one particle per channel); :attr:`PhysimChannels.left_out` says how much each
took.

**Limits.** One E* for the whole target (the spread of E* through the target is not followed into the
evaporation). The evaporation's known differences from PACE4 (``docs/theory/evaporation.md``): neutrons 10–25%
softer, so the xn distribution leans one neutron further for light projectiles; for heavy projectiles at high
spin it leaves too much spin in the residue and its channel shares can be off by a large factor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from . import data
from .residues import ChannelTable, ResidueChannel, parse_particles

#: Where the channel table's numbers come from, as written in each channel's ``source``.
SOURCE = "physim (Bass fusion + evaporation)"
#: Evaporation chains per table by default.
EVENTS = 20_000
#: Channels with a smaller share of the events are left out of the table.
MIN_SHARE = 0.002
#: Beams this heavy or heavier (mass number), or mean ℓ this high (ħ), are where the evaporation is known to
#: disagree strongly with PACE4 (docs/theory/evaporation.md, the ³²S and ⁴⁸Ca comparisons).
HEAVY_BEAM_A = 28
HIGH_SPIN_L = 30.0


@dataclass
class PhysimChannels:
    """The channel table made by :func:`physim_channels`, with what it was made from."""

    table: ChannelTable
    #: σ_fus averaged over the target, mb.
    sigma_fusion_mb: float
    #: The E* the evaporation started from (σ-weighted mean through the target), and its range through the target,
    #: MeV.
    excitation_mev: float
    excitation_range_mev: tuple
    #: ⟨ℓ⟩ and the highest ℓ that fuses, ħ.
    mean_l: float
    max_l: int
    #: The l = 0 fusion barrier: height (CM) and the lab beam energy it corresponds to, MeV; radius, fm; ħω, MeV.
    barrier: dict
    critical_l: int
    #: Shares of the events not in the table: "fission", "no particle" (γ rays only) and "small channels".
    left_out: dict
    #: Mean kinetic energies of n, p and α in the compound nucleus's frame, MeV, and mean multiplicities.
    mean_energy_mev: dict
    multiplicity: dict
    events: int
    seed: int
    warnings: list = field(default_factory=list)

    def summary(self) -> dict:
        """Plain numbers for the app, the MCP tools and the report."""
        return {"source": SOURCE, "sigma_fusion_mb": self.sigma_fusion_mb, "excitation_mev": self.excitation_mev,
                "excitation_range_mev": self.excitation_range_mev, "mean_l": self.mean_l, "max_l": self.max_l,
                "barrier": self.barrier, "critical_l": self.critical_l, "left_out": self.left_out,
                "mean_energy_mev": self.mean_energy_mev, "multiplicity": self.multiplicity,
                "events": self.events, "seed": self.seed, "warnings": list(self.warnings),
                "channels": [{"channel": ch.label, "residue": ch.residue, "cross_section_mb": ch.cross_section_mb,
                              "entry_spin": ch.entry_spin, "entry_energy_mev": ch.entry_energy_mev}
                             for ch in self.table.channels]}


def _target(experiment) -> tuple:
    """The beam ion, the fused nuclide, and the target's areal density along the beam (mg/cm²)."""
    from .experiment import main_nuclide
    from .rates import beam_ion, stack, tilt_deg

    layer = stack(experiment)[0]
    along = layer.thickness / math.cos(math.radians(tilt_deg(experiment)))
    return beam_ion(experiment), main_nuclide(experiment.target), along


def fusion_of(experiment):
    """The :class:`physim.nuclear.fusion.Fusion` of the setup's beam on the target's main nuclide."""
    from .fusion import Fusion

    ion, target, _ = _target(experiment)
    return Fusion(ion, target)


def thick_target(experiment):
    """The setup's fusion averaged over its target (:meth:`physim.nuclear.fusion.Fusion.thick_target`)."""
    _, _, along = _target(experiment)
    t = experiment.target
    density = t.density if t.density is not None else None
    return fusion_of(experiment).thick_target(experiment.beam.energy_mev, along, t.material, density)


def fusion_summary(experiment, points: int = 41) -> dict:
    """The fusion of the setup's beam and target, for the Plan tab: the barrier, σ_fus at the beam energy and
    averaged over the target, ⟨ℓ⟩ and the critical ℓ, the excitation function (thin target) from just below the
    barrier to 1.6 times it or the beam energy, and the model's warnings."""
    f = fusion_of(experiment)
    tt = thick_target(experiment)
    e = experiment.beam.energy_mev
    b = f.barrier()
    lab_barrier = float(f.barrier_lab_mev)
    cn = f.compound(e)
    hi = max(1.6 * lab_barrier, 1.1 * e)
    energies = np.linspace(0.85 * lab_barrier, hi, points)
    good = tt.sigma_depth_mb > 0
    ex = tt.excitation_mev[good] if good.any() else tt.excitation_mev[:1]
    weights = tt.sigma_depth_mb[good] if good.any() else np.ones(1)
    _, _, along = _target(experiment)
    return {"beam": f.beam.name, "target": f.target.name, "compound": f.compound_nuclide.name,
            "barrier_cm_mev": b.height_mev, "barrier_lab_mev": lab_barrier, "barrier_radius_fm": b.radius_fm,
            "hbar_omega_mev": b.hbar_omega_mev, "critical_l": f.critical_l,
            "sigma_mb": cn.sigma_mb, "sigma_thick_mb": tt.sigma_mb, "mean_l": cn.mean_l, "max_l": cn.max_l,
            "excitation_mev": float(np.average(ex, weights=weights)),
            "excitation_range_mev": (float(np.nanmin(tt.excitation_mev)), float(np.nanmax(tt.excitation_mev))),
            "e_in_mev": tt.e_in_mev, "e_out_mev": tt.e_out_mev,
            "excitation_function": {"beam_energy_mev": energies,
                                    "sigma_mb": np.atleast_1d(f.cross_section(energies))},
            "warnings": f.warnings(e, along, experiment.target.material, experiment.target.density)}


def _entry_spin(j: float, residue: str) -> float:
    """A mean spin rounded down to the residue's kind of spin (whole for even A, half-whole for odd A)."""
    odd = data.nuclide(residue).A % 2
    j = max(j, 0.5 * odd)
    return math.floor(j - 0.5 * odd) + 0.5 * odd


def physim_channels(experiment, events: int = EVENTS, seed: int = 1, min_share: float = MIN_SHARE) -> PhysimChannels:
    """The exit channels of the setup's fusion-evaporation reaction from physim's own fusion and evaporation (see
    the module docstring). The table keeps the setup's defaults (temperature, alignment, side feeding)."""
    from .evaporation import evaporate

    if experiment.fusion is None:
        raise ValueError("the setup is not a fusion-evaporation reaction ([reaction] type = \"fusion-evaporation\")")
    ion, target, along = _target(experiment)
    f = fusion_of(experiment)
    tt = thick_target(experiment)
    if tt.sigma_mb <= 0:
        raise ValueError(f"nothing fuses: the beam ({experiment.beam.energy_mev:.4g} MeV) is far below the fusion "
                         f"barrier ({f.barrier_lab_mev:.4g} MeV in the lab)")
    good = tt.sigma_depth_mb > 0
    e_star = float(np.average(tt.excitation_mev[good], weights=tt.sigma_depth_mb[good]))
    ev = evaporate(f.compound_nuclide.name, e_star, sigma_l_mb=tt.sigma_l_mb, events=events, seed=seed)

    base = ChannelTable.from_experiment(experiment)
    table = ChannelTable(base.beam, base.target, temperature_mev=base.temperature_mev, alignment=base.alignment,
                         side_feeding=base.side_feeding)
    left = {"fission": 0.0, "no particle": 0.0, "small channels": 0.0}
    labels = ev.label_of_event
    for row in ev.channel_table():
        label = row["label"]
        if label == "fission":
            left["fission"] += row["share"]
            continue
        if label == "g":
            left["no particle"] += row["share"]
            continue
        if row["share"] < min_share:
            left["small channels"] += row["share"]
            continue
        sel = labels == label
        residue = row["residue"]
        table.channels.append(ResidueChannel(
            parse_particles(label), float(row["sigma_mb"]), residue, label,
            entry_energy_mev=float(ev.entry_energy_mev[sel].mean()),
            entry_spin=_entry_spin(float(ev.entry_spin[sel].mean()), residue), source=SOURCE))
    if not table.channels:
        raise ValueError("the evaporation left no channel with particles above the cut: every event fissioned or "
                         "emitted only γ rays")

    p = tt.sigma_l_mb / tt.sigma_l_mb.sum()
    mean_l = float(np.dot(np.arange(p.size), p))
    nz = np.nonzero(tt.sigma_l_mb > 1e-6 * tt.sigma_l_mb.max())[0]
    b = f.barrier()
    warns = list(f.warnings(experiment.beam.energy_mev, along, experiment.target.material,
                            experiment.target.density))
    if data.nuclide(ion).A >= HEAVY_BEAM_A or mean_l >= HIGH_SPIN_L:
        warns.append(f"High spin (⟨ℓ⟩ = {mean_l:.0f} ħ, beam A = {data.nuclide(ion).A}): physim's evaporation leaves "
                     "too much spin in the residues here and its channel shares can differ from PACE4 by a large "
                     "factor (docs/theory/evaporation.md). Check the table against PACE4 before planning on it.")
    warns.append("physim's evaporation gives neutrons 10–25% softer than PACE4's, which can shift the xn channels "
                 "by about one neutron; the table is an estimate, PACE4 or a measurement is the better source.")
    if ev.truncated:
        warns.append(f"{ev.truncated} of {events} evaporation chains were cut short.")
    if left["fission"] > 0.05:
        warns.append(f"{100 * left['fission']:.0f}% of the compound nuclei fission; those events make no residue.")
    return PhysimChannels(
        table=table, sigma_fusion_mb=tt.sigma_mb, excitation_mev=e_star,
        excitation_range_mev=(float(np.nanmin(tt.excitation_mev[good])), float(np.nanmax(tt.excitation_mev[good]))),
        mean_l=mean_l, max_l=int(nz[-1]) if nz.size else 0,
        barrier={"height_cm_mev": b.height_mev, "lab_mev": float(f.barrier_lab_mev), "radius_fm": b.radius_fm,
                 "hbar_omega_mev": b.hbar_omega_mev},
        critical_l=f.critical_l, left_out=left,
        mean_energy_mev={k: ev.mean_energy_mev(k) for k in ("n", "p", "a")},
        multiplicity={k: ev.multiplicity(k) for k in ("n", "p", "a")}, events=events, seed=seed, warnings=warns)


# ---------------------------------------------------------------------------------------------------------------
# Channel selection: which channel each event looks like to the neutron wall and the particle detectors

#: Residue events the channel selection follows by default.
SELECTION_EVENTS = 50_000
#: Statistical γ rays per event on top of the stretched E2 cascade from the entry spin (typical: about 4 for
#: fusion-evaporation, from the usual conversion J ≈ 2 (M_γ − 4) between γ-ray multiplicity and spin).
STATISTICAL_GAMMAS = 4.0
#: Energy, MeV, a proton or α particle must deposit in a particle detector that sets no threshold of its own
#: (typical of silicon and CsI arrays for evaporated particles; ``[neutron_wall] particle_threshold``).
PARTICLE_THRESHOLD_MEV = 1.0
#: Ion of each charged particle kind, for its energy loss in the detectors.
_IONS = {"p": "1H", "a": "4He"}


class _Charged:
    """Which particle detector a proton or α particle leaving the target meets first (dead material in front
    stops it), and whether it deposits more than that detector's threshold: its own ``threshold``, or
    ``default_mev`` when it sets none. The deposit is after the dead layer and through the thickness (a particle
    that punches through deposits only part of its energy). Energy lost in the target on the way out is left out."""

    def __init__(self, experiment, default_mev: float):
        from .detectors import Array

        self.array = Array.from_experiment(experiment)
        self.setup = list(experiment.detectors)
        self.thresholds = [(d.name or g.name, _threshold_mev(d, default_mev), d.threshold is not None)
                           for d, g in zip(self.setup, self.array.geometries)]
        self._responses: dict = {}

    def first(self, direction) -> tuple:
        """(index of the detector met, or −1; incidence on its face, degrees) for each direction."""
        d = np.atleast_2d(np.asarray(direction, dtype=float))
        if not self.array.geometries:
            return np.full(len(d), -1), np.zeros(len(d))
        dist = np.array([g.distances(d) for g in self.array.geometries])
        which = np.argmin(dist, axis=0)
        own = dist[which, np.arange(len(d))]
        ok = np.isfinite(own)
        for b in self.array.blockers:
            ok &= ~(b.distances(d) < own - 1e-9)
        normals = np.array([g.n for g in self.array.geometries])[which]
        unit = d / np.linalg.norm(d, axis=1, keepdims=True)
        incidence = np.degrees(np.arccos(np.clip(-np.sum(unit * normals, axis=1), -1, 1)))
        return np.where(ok, which, -1), incidence

    def reached(self, direction) -> np.ndarray:
        return (self.first(direction)[0] >= 0).astype(float)

    def response(self, i: int, kind: str):
        from .detectors import Response

        key = (i, kind)
        if key not in self._responses:
            self._responses[key] = Response(_IONS[kind], self.setup[i])
        return self._responses[key]

    def counted(self, kind: str):
        """f(energy, direction) → 1 where a particle of this kind reaches a detector and deposits more than its
        threshold, else 0."""
        def f(energy, direction):
            which, incidence = self.first(direction)
            e = np.broadcast_to(np.asarray(energy, dtype=float), which.shape)
            out = np.zeros(len(which))
            for i in np.unique(which[which >= 0]):
                sel = which == i
                dep = self.response(int(i), kind).deposited_each(e[sel], incidence[sel])
                out[sel] = dep >= self.thresholds[i][1]
            return out
        return f


def _threshold_mev(detector, default_mev: float) -> float:
    from .detectors import _q

    return default_mev if detector.threshold is None else _q(detector.threshold).to("MeV")


def _gamma_numbers(entry_spin: np.ndarray, fixed: Optional[float], rng) -> np.ndarray:
    """γ rays per event: ``fixed`` if the setup gives it, else J/2 + STATISTICAL_GAMMAS from each event's entry spin
    J; a non-whole mean is rounded up or down at random so that the mean is kept."""
    mean = np.full(len(entry_spin), float(fixed)) if fixed is not None else entry_spin / 2 + STATISTICAL_GAMMAS
    low = np.floor(mean)
    return (low + (rng.random(len(mean)) < mean - low)).astype(int)


def channel_selection(experiment, events: int = SELECTION_EVENTS, seed: int = 1) -> dict:
    """The channel-selection matrix of the setup's neutron wall and particle detectors, for its channel table.

    The residue events of Stage 1 (:func:`physim.nuclear.residues.simulate_residues`) give each event's neutrons,
    protons and α particles with their lab energies and directions. A neutron is counted with the wall's
    efficiency at its energy where it flies into the wall's cone (:meth:`physim.nuclear.neutrons.NeutronWall.
    detection_probability`), fires a neighbour with the wall's cross-talk at the mean neutron energy, and is lost or
    taken for a γ ray as the cells' pulse-shape discrimination says. A proton or α particle is counted when it
    meets a particle detector of the setup and deposits more than the detector's threshold (``threshold``, or
    ``[neutron_wall] particle_threshold`` for detectors that set none; after the dead layer, and only part of its
    energy if it punches through). Each event also has γ rays (``[neutron_wall] gamma_multiplicity``, or J/2 +
    :data:`STATISTICAL_GAMMAS` from its entry spin J); each fires a cell with the wall's geometric efficiency times
    the cells' ``gamma_efficiency``, and adds a false neutron when the cell calls it a neutron
    (``gamma_as_neutron``).

    Each event is sorted by the counts it shows: its detected class, such as "1n" or "p2n"
    (:func:`physim.nuclear.neutrons.simulate_selection`, with every class up to the most particles a channel emits
    as a column; an event showing more, through cross-talk or false neutrons, is "unassigned"). With a real wall few
    events show every neutron, so the classes are what one gates on. Returns the matrix (rows: true channels,
    columns: the detected classes), its statistical errors, the analytic matrix from the same events' mean
    probabilities (:func:`physim.nuclear.neutrons.selection_matrix`, which the Monte Carlo converges to), for each
    class the share of its events from each true channel for the table's cross sections (what a gate on that class
    contains), and the efficiencies, thresholds and γ-ray numbers used."""
    from . import neutrons
    from .residues import LIGHT, simulate_residues

    if experiment.fusion is None or experiment.neutron_wall is None:
        raise ValueError("channel selection needs a fusion-evaporation reaction and a [neutron_wall]")
    setup = experiment.neutron_wall
    wall = setup.model()
    ev = simulate_residues(experiment, events, seed)
    table = ev.table
    kinds = list(LIGHT)
    usable, skipped = [], []
    for ch in table.channels:
        (skipped if any(k in ("d", "t") for k in ch.particles) else usable).append(ch)
    true = {(ch.particles.count("n"), ch.particles.count("p"), ch.particles.count("a")): ch.label for ch in usable}
    top = [max(m[k] for m in true) for k in range(3)]
    chans = [neutrons.Channel(true.get((n, pr, a), neutrons.channel_name(n, pr, a)), n, pr, a)
             for a in range(top[2] + 1) for pr in range(top[1] + 1) for n in range(top[0] + 1)]

    rng = np.random.default_rng([seed, 75])
    fixed = setup.gamma_multiplicity
    gammas = _gamma_numbers(np.asarray(ev["entry_spin"], dtype=float), fixed, rng)
    part = ev.particles
    th, ph = np.radians(part["theta"]), np.radians(part["phi"])
    dirs = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=1)
    letters = np.array(kinds)[part["kind"]]
    order = np.argsort(part["event"], kind="stable")
    bounds = np.searchsorted(part["event"][order], np.arange(len(ev) + 1))
    event_list = []
    for i in range(len(ev)):
        rows = order[bounds[i]:bounds[i + 1]]
        event_list.append({"particles": [(letters[r], float(part["energy"][r]), dirs[r]) for r in rows],
                           "gammas": int(gammas[i])})

    is_n = letters == "n"
    mean_n = float(part["energy"][is_n].mean()) if is_n.any() else 2.0
    crosstalk = wall.crosstalk(mean_n)
    from .detectors import _q

    default_mev = (PARTICLE_THRESHOLD_MEV if setup.particle_threshold is None
                   else _q(setup.particle_threshold).to("MeV"))
    charged = _Charged(experiment, default_mev)
    gamma_cell = neutrons._fraction(wall.cell.gamma_efficiency)
    detection = neutrons.Detection.from_detector(
        wall.cell, wall.detection_probability, crosstalk, proton_efficiency=charged.counted("p"),
        alpha_efficiency=charged.counted("a"), gamma_efficiency=gamma_cell * wall.geometric_efficiency())
    m = neutrons.simulate_selection(event_list, chans, detection, seed=seed)
    analytic = neutrons.selection_matrix(neutrons.gamma_distributions(event_list, chans),
                                         detection.averaged(event_list, chans))
    rows = [i for i, name in enumerate(m.true) if name in true.values()]
    names = [m.true[i] for i in rows]
    shown = np.nonzero(m.matrix[rows].max(axis=0) > 0)[0]          # the classes some event shows
    assigned = [m.assigned[j] for j in shown]
    matrix, err, counts = m.matrix[rows][:, shown], m.uncertainty()[rows][:, shown], m.events[rows]
    sigma = {ch.label: ch.cross_section_mb for ch in usable}
    w = np.array([sigma[n] for n in names])
    per_class = w @ matrix                                   # mb of events showing each class
    contents = {}
    for j, cls in enumerate(assigned):
        if per_class[j] > 0:
            contents[cls] = {"cross_section_mb": float(per_class[j]),
                             "from": {names[i]: float(w[i] * matrix[i, j] / per_class[j]) for i in range(len(names))
                                      if matrix[i, j] > 0}}
    reached, counted = {}, {}
    for k in ("p", "a"):
        sel = letters == k
        reached[k] = float(charged.reached(dirs[sel]).mean()) if sel.any() else None
        counted[k] = float(charged.counted(k)(part["energy"][sel], dirs[sel]).mean()) if sel.any() else None
    per_gamma = detection.gamma_fake()
    listed = ", ".join(f"{n} {t:.3g} MeV" + ("" if own else " (default)") for n, t, own in charged.thresholds)
    how_many = (f"a fixed {fixed:g}" if fixed is not None
                else f"J/2 + {STATISTICAL_GAMMAS:g} from the entry spin J")
    notes = [f"Protons and α particles count when they deposit more than the particle detector's threshold "
             f"({listed or 'no particle detectors'}); the energy they lose in the target, and the detectors' own "
             "particle identification, are left out.",
             f"γ rays: {how_many} per event (mean {gammas.mean():.1f}); each fires a cell with "
             f"{100 * detection.gamma_efficiency:.1f}% and is called a neutron with "
             f"{100 * detection.gamma_as_neutron:.1f}%: a false neutron for {100 * per_gamma:.3f}% of the γ rays. "
             "Time of flight is taken to be part of that γ-as-neutron share.",
             f"Cross-talk is taken at the mean neutron energy ({mean_n:.2f} MeV in the lab)."]
    if skipped:
        notes.append("Channels with deuterons or tritons are left out of the selection: "
                     + ", ".join(ch.label for ch in skipped) + ".")
    a_rows = [analytic.true.index(n) for n in names]
    a_cols = [analytic.assigned.index(c) for c in assigned]
    return {"true": names, "assigned": assigned, "matrix": matrix.tolist(), "errors": err.tolist(),
            "analytic": analytic.matrix[a_rows][:, a_cols].tolist(),
            "events": counts.tolist(), "contents": contents, "cross_sections_mb": sigma,
            "neutron": {"cells": wall.cells, "solid_angle_msr": wall.solid_angle(),
                        "geometric_efficiency": wall.geometric_efficiency(), "theta_max_deg": wall.theta_max(),
                        "mean_energy_mev": mean_n,
                        "intrinsic_at_mean": float(wall.cell.intrinsic_efficiency(mean_n)),
                        "efficiency_at_mean": float(wall.efficiency(mean_n)), "crosstalk": crosstalk,
                        "misidentification": wall.cell.misidentification},
            "charged_reached": reached, "charged_counted": counted,
            "thresholds_mev": {n: t for n, t, _ in charged.thresholds},
            "gamma": {"multiplicity": float(gammas.mean()), "from_entry_spin": fixed is None,
                      "fires_a_cell": detection.gamma_efficiency, "false_neutron_per_gamma": per_gamma},
            "notes": notes, "n_events": len(ev), "seed": seed}


__all__ = ["EVENTS", "MIN_SHARE", "PARTICLE_THRESHOLD_MEV", "SELECTION_EVENTS", "SOURCE", "STATISTICAL_GAMMAS",
           "PhysimChannels", "channel_selection", "fusion_of",
           "fusion_summary", "physim_channels", "thick_target"]
