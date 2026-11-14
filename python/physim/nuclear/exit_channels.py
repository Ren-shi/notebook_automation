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


def _charged_hit(experiment):
    """f(energy, directions) → 1 where a particle leaving the target along the direction meets a particle
    detector's face before any dead material, else 0 (every particle that arrives is counted and identified:
    thresholds, punch-through and the detector's own identification are left out)."""
    from .detectors import Array

    array = Array.from_experiment(experiment)

    def hit(energy, direction):
        d = np.atleast_2d(np.asarray(direction, dtype=float))
        if not array.geometries:
            return np.zeros(len(d))
        own = np.min([g.distances(d) for g in array.geometries], axis=0)
        ok = np.isfinite(own)
        for b in array.blockers:
            ok &= ~(b.distances(d) < own - 1e-9)
        return ok.astype(float)

    return hit


def channel_selection(experiment, events: int = SELECTION_EVENTS, seed: int = 1) -> dict:
    """The channel-selection matrix of the setup's neutron wall and particle detectors, for its channel table.

    The residue events of Stage 1 (:func:`physim.nuclear.residues.simulate_residues`) give each event's neutrons,
    protons and α particles with their lab energies and directions. A neutron is counted with the wall's
    efficiency at its energy where it flies into the wall's cone (:meth:`physim.nuclear.neutrons.NeutronWall.
    detection_probability`), fires a neighbour with the wall's cross-talk at the mean neutron energy, and is lost or
    taken for a γ ray as the cells' pulse-shape discrimination says; a proton or α particle is counted when it
    meets a particle detector of the setup.

    Each event is sorted by the counts it shows: its detected class, such as "1n" or "p2n"
    (:func:`physim.nuclear.neutrons.simulate_selection`, with every class up to the most particles a channel emits
    as a column; an event showing more, through cross-talk, is "unassigned"). With a real wall few events show
    every neutron, so the classes are what one gates on. Returns the matrix (rows: true channels, columns: the
    detected classes), its statistical errors, for each class the share of its events from each true channel for
    the table's cross sections (what a gate on that class contains), and the efficiencies used."""
    from . import neutrons
    from .residues import LIGHT, simulate_residues

    if experiment.fusion is None or experiment.neutron_wall is None:
        raise ValueError("channel selection needs a fusion-evaporation reaction and a [neutron_wall]")
    wall = experiment.neutron_wall.model()
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

    part = ev.particles
    th, ph = np.radians(part["theta"]), np.radians(part["phi"])
    dirs = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=1)
    letters = np.array(kinds)[part["kind"]]
    order = np.argsort(part["event"], kind="stable")
    bounds = np.searchsorted(part["event"][order], np.arange(len(ev) + 1))
    event_list = []
    for i in range(len(ev)):
        rows = order[bounds[i]:bounds[i + 1]]
        event_list.append([(letters[r], float(part["energy"][r]), dirs[r]) for r in rows])

    is_n = letters == "n"
    mean_n = float(part["energy"][is_n].mean()) if is_n.any() else 2.0
    crosstalk = wall.crosstalk(mean_n)
    charged = _charged_hit(experiment)
    detection = neutrons.Detection.from_detector(wall.cell, wall.detection_probability, crosstalk,
                                                 proton_efficiency=charged, alpha_efficiency=charged)
    m = neutrons.simulate_selection(event_list, chans, detection, seed=seed)
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
    charged_share = {k: float(charged(None, dirs[letters == k]).mean()) if (letters == k).any() else None
                     for k in ("p", "a")}
    notes = ["Every proton and α particle that meets a particle detector is counted and told apart; thresholds, "
             "punch-through and the detectors' own identification are left out.",
             "γ rays that fire a neutron cell are left out (the events do not carry the cascades' γ rays here).",
             f"Cross-talk is taken at the mean neutron energy ({mean_n:.2f} MeV in the lab)."]
    if skipped:
        notes.append("Channels with deuterons or tritons are left out of the selection: "
                     + ", ".join(ch.label for ch in skipped) + ".")
    return {"true": names, "assigned": assigned, "matrix": matrix.tolist(), "errors": err.tolist(),
            "events": counts.tolist(), "contents": contents, "cross_sections_mb": sigma,
            "neutron": {"cells": wall.cells, "solid_angle_msr": wall.solid_angle(),
                        "geometric_efficiency": wall.geometric_efficiency(), "theta_max_deg": wall.theta_max(),
                        "mean_energy_mev": mean_n,
                        "intrinsic_at_mean": float(wall.cell.intrinsic_efficiency(mean_n)),
                        "efficiency_at_mean": float(wall.efficiency(mean_n)), "crosstalk": crosstalk,
                        "misidentification": wall.cell.misidentification},
            "charged_reached": charged_share, "notes": notes, "n_events": len(ev), "seed": seed}


__all__ = ["EVENTS", "MIN_SHARE", "SELECTION_EVENTS", "SOURCE", "PhysimChannels", "channel_selection", "fusion_of",
           "fusion_summary", "physim_channels", "thick_target"]
