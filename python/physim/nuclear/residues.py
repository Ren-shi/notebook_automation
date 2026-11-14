"""Fusion-evaporation, Stage 1 (backlog item 75): the residues of a channel table, their kinematics and γ rays.

The beam fuses with the target's main nuclide into a compound nucleus, which boils off light particles (neutrons,
protons, α particles) and leaves a residue; the residue's γ-ray cascade is what a germanium array sees. Stage 1
takes the cross section of each exit channel from a table (typed into the setup, or later from
:mod:`physim.nuclear.fusion` and :mod:`physim.nuclear.evaporation`) and works out everything else::

    from physim.nuclear import Experiment
    from physim.nuclear.residues import ChannelTable, compound_profile, simulate_residues, simulate_residue_gammas

    exp = Experiment.example("fusion_dy158")
    table = ChannelTable.from_experiment(exp)       # the channels, their cross sections and residues
    cn = table.compound_at(exp.beam.energy_mev)     # E* = E_cm + Q and the compound nucleus's speed, exact
    ev = simulate_residues(exp, events=100_000, seed=1)
    ev.per_channel()                                # residues per second in each channel, with errors
    g = simulate_residue_gammas(exp, events=100_000, seed=1)
    g.spectrum("Ge90", corrected=True)              # the Doppler-corrected γ-ray spectrum, counts in the run

**The compound nucleus.** E* = √s − M_CN, with √s the invariant mass of beam and target and the nuclear masses of
AME2020: this is E_cm + Q with E_cm = √s − m_b − m_t, exact (relativistic). It moves along the beam with
β_CN = p_b / (E_b + m_t), which is v_beam m_b / (m_b + m_t) in the non-relativistic limit. Through the target the
beam slows down (:mod:`physim.nuclear.stopping`): each event happens at a depth drawn evenly through the target,
with the mean beam energy there and its spread (the beam's own spread carried through the layer, plus
straggling), so E* and β_CN spread as the target makes them.

**The evaporation.** Each event picks a channel with probability σ_i / Σσ and emits the channel's particles one
after the other. In the frame of the emitting nucleus a particle leaves isotropically, with a relative kinetic
energy ε drawn from the evaporation spectrum (ε − V) exp(−(ε − V)/T), V being the Coulomb barrier for a charged
particle (0 for a neutron) and T the temperature (the channel's, the reaction's, or √(E*/a) with a = A/8 MeV⁻¹ of
the emitting nucleus). Energy is conserved exactly: the daughter keeps E* − S − ε, S being the particle's
separation energy from the masses, and ε never takes more than leaves the particles still to come their
separation energies and barriers. The chain must end with a residue that cannot emit another particle (E* below
its lowest separation energy plus barrier, for n, p and α); a chain that does not is drawn again, up to
:data:`CHAIN_ATTEMPTS` times, and then its last particle is drawn from the tail of its spectrum above what it must
carry away (the notes say how often). Momentum is conserved particle by particle with a Lorentz boost from the
emitter's frame to the laboratory, so the residue's direction and speed (the Doppler cone) come from the sum of
the recoils.

**The γ rays.** The residue enters its known levels (``[levels.<residue>]`` in the setup, or the local ENSDF copy)
at the lowest level of the highest spin not above the entry spin (the yrast entry); a share ``side_feeding`` of
the cascades enters instead at one of the lower yrast levels. The cascade then follows the level scheme's
branching ratios (with internal conversion). The compound nucleus's spin is perpendicular to the beam, so the
levels are aligned: a Gaussian distribution of magnetic substates about the beam of width σ = ``alignment`` × J
at the entry level, carried down by the U_k coefficients through every γ ray, and each γ ray has the angular
distribution W(θ) = 1 + Σ B_k A_k P_k(cos θ) about the beam in the residue's frame. The γ rays are emitted
where the residue was made, at the speed it was made with (prompt decays; lifetimes are left out), followed to
the crystals with the response of :mod:`physim.nuclear.response`, and Doppler corrected as an experimentalist
does without a recoil detector: as if every residue moved along the beam at β_CN for the nominal beam energy at
the middle of the target.

**Not in Stage 1:** the fusion cross section (Stage 2) and the evaporation's own channel cross sections
(Stage 3); the energy dependence of the cross sections through the target (each channel's cross section is the
table's, at every depth); lifetimes and slowing down of the residues; summing and add-back of the residues'
γ rays; the statistical γ rays above the yrast line; and particle detection of the evaporated particles.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, replace
from typing import Iterable, Optional

import numpy as np

from . import angular, data
from .names import parse_nuclide
from .quantity import Quantity

#: The light particles a compound nucleus can emit: symbol as written in a channel, and the nuclide.
LIGHT = {"n": "n", "p": "1H", "d": "2H", "t": "3H", "a": "4He"}
#: How each light particle is printed.
LIGHT_LABEL = {"n": "n", "p": "p", "d": "d", "t": "t", "a": "α"}
#: Order of the light particles in a printed channel: "p2n", "α2n".
_ORDER = ("a", "t", "d", "p", "n")

#: The level-density parameter a = A / LEVEL_DENSITY_DIVISOR (MeV⁻¹) of the default temperature T = √(E*/a).
LEVEL_DENSITY_DIVISOR = 8.0
#: Radius parameter of the Coulomb barrier for charged-particle emission, fm: V = Z_x Z_d e² / (r0 (A_x^⅓ + A_d^⅓)).
EMISSION_R0 = 1.5
#: Default width σ/J of the compound nucleus's magnetic-substate distribution about the beam.
ALIGNMENT = 0.3
#: Default share of the cascades that enter below the top yrast level.
SIDE_FEEDING = 0.25
#: Angular momentum each kind of particle takes away on average, ħ, for the default entry spin.
SPIN_REMOVED = {"n": 1.0, "p": 1.0, "d": 1.5, "t": 2.0, "a": 3.0}
#: e² in MeV fm.
E2_MEV_FM = 1.439964547
#: The γ-ray columns of :class:`ResidueGammas`, with units.
GAMMA_COLUMNS = {
    "event": "event number, as in the residue events",
    "cascade": "which cascade the γ ray belongs to (γ rays of one cascade share it)",
    "channel": "index into the channel table",
    "initial": "level the γ ray leaves (index in the residue's level scheme)",
    "final": "level the γ ray goes to",
    "crystal": "index into ResidueGammas.crystals",
    "energy0": "γ-ray energy in the residue's frame, MeV",
    "energy_lab": "γ-ray energy in the laboratory (Doppler shifted), MeV",
    "deposited": "energy left in the crystal, MeV",
    "measured": "measured energy (deposited, with resolution), MeV",
    "counted": "above the crystal's threshold",
    "theta": "lab polar angle of the γ ray, degrees",
    "phi": "lab azimuth of the γ ray, degrees",
    "beta": "speed of the residue at emission, units of c",
    "corrected": "measured energy Doppler corrected for a residue along the beam at the nominal β_CN, MeV",
    "weight": "rate this γ ray stands for, per second",
}


# ---------------------------------------------------------------------------------------------------------------
# Channels: which particles, which residue


def parse_particles(text: str) -> tuple:
    """The particles of a channel, in the order written: ``"p2n"`` → ("p", "n", "n"); ``"a n"`` and ``"αn"`` →
    ("a", "n"). Each is one of n, p, d, t and a (or α). Raises ValueError with a readable message."""
    if not isinstance(text, str):
        raise ValueError(f"must be text such as \"4n\" or \"p2n\", got {text!r}")
    s = re.sub(r"[\s,+]", "", text.replace("α", "a"))
    if not s:
        raise ValueError("must name at least one particle, e.g. \"4n\", \"p2n\" or \"a n\"")
    if not re.fullmatch(r"(?:\d*[npdta])+", s):
        raise ValueError(f"'{text}' is not a list of particles: write each as a count and n, p, d, t or a (α), "
                         "e.g. \"4n\", \"p2n\" or \"a 2n\"")
    out = []
    for count, kind in re.findall(r"(\d*)([npdta])", s):
        k = int(count) if count else 1
        if k < 1:
            raise ValueError(f"'{text}': a count must be at least 1")
        out += [kind] * k
    if len(out) > 20:
        raise ValueError(f"'{text}' names {len(out)} particles; at most 20 are followed")
    return tuple(out)


def particles_label(particles: Iterable[str]) -> str:
    """The particles as a channel is usually written: ("n", "p", "n") → "p2n", ("a", "n") → "αn"."""
    parts = list(particles)
    out = ""
    for k in _ORDER:
        n = parts.count(k)
        if n:
            out += (str(n) if n > 1 else "") + LIGHT_LABEL[k]
    return out


def residue_of(compound: tuple, particles: Iterable[str]) -> tuple:
    """(Z, A) of the residue a compound nucleus (Z, A) leaves after emitting ``particles``. Raises ValueError if
    nothing would be left, or the mass table does not know it."""
    z, a = compound
    for p in particles:
        pz, pa = parse_nuclide(LIGHT[p])
        z, a = z - pz, a - pa
    if z < 1 or a <= z:
        raise ValueError("nothing would be left of the compound nucleus")
    data.nuclide((z, a))  # raises if AME does not know it
    return (z, a)


# ---------------------------------------------------------------------------------------------------------------
# The compound nucleus


@dataclass
class CompoundState:
    """The compound nucleus made by a beam of lab kinetic energy ``beam_energy`` (MeV; an array works too)."""

    nuclide: str
    beam_energy: object
    #: Q-value of the fusion, m_b + m_t − M_CN (nuclear masses), MeV.
    q_value: float
    #: Kinetic energy in the CM frame, √s − m_b − m_t, MeV.
    e_cm: object
    #: Excitation energy √s − M_CN = E_cm + Q, MeV.
    excitation: object
    #: Speed of the compound nucleus along the beam, units of c, and its Lorentz factor.
    beta: object
    gamma: object

    @property
    def kinetic_energy(self):
        """Kinetic energy of the compound nucleus in the laboratory, MeV."""
        m = data.nuclide(self.nuclide).nuclear_mass_mev + np.asarray(self.excitation)
        return m * (np.asarray(self.gamma) - 1.0)


def compound(beam: str, target: str, beam_energy) -> CompoundState:
    """The compound nucleus of ``beam`` on ``target`` at a lab beam energy (MeV, text with a unit, or an array of
    MeV), exactly: E* = √s − M_CN and β_CN = p_b / (E_b + m_t)."""
    if isinstance(beam_energy, str):
        beam_energy = Quantity.parse(beam_energy).to("MeV")
    b, t = data.nuclide(beam), data.nuclide(target)
    c = data.nuclide((b.Z + t.Z, b.A + t.A))
    mb, mt, mc = b.nuclear_mass_mev, t.nuclear_mass_mev, c.nuclear_mass_mev
    e = np.asarray(beam_energy, dtype=float)
    total = e + mb + mt
    p = np.sqrt(e * e + 2.0 * e * mb)
    roots = np.sqrt(total * total - p * p)
    out = CompoundState(c.name, beam_energy, mb + mt - mc, roots - mb - mt, roots - mc, p / total, total / roots)
    if np.ndim(e) == 0:
        for name in ("e_cm", "excitation", "beta", "gamma"):
            setattr(out, name, float(getattr(out, name)))
    return out


def compound_profile(experiment, points: int = 41) -> dict:
    """The compound nucleus through the target: at ``points`` depths (mg/cm² along the normal) the mean beam energy,
    its standard deviation (the beam's spread carried through the layer, and straggling), and E* and β_CN at the
    mean energy, with their standard deviations from the energy spread."""
    from .rates import _carry, beam_ion, energy_sigma, stack, stopping, tilt_deg
    from .experiment import main_nuclide

    layers = stack(experiment)
    ion = beam_ion(experiment)
    st = stopping(ion, layers[0].material)
    c = math.cos(math.radians(tilt_deg(experiment)))
    z = np.linspace(0.0, layers[0].thickness, points)
    mean, sigma = np.zeros(points), np.zeros(points)
    e0, var0 = experiment.beam.energy_mev, energy_sigma(experiment) ** 2
    for k, depth in enumerate(z):
        e, var = _carry(st, e0, var0, depth / c) if depth > 0 else (e0, var0)
        mean[k], sigma[k] = e, math.sqrt(max(var, 0.0))
    target = main_nuclide(experiment.target)
    cn = compound(ion, target, mean)
    de = 1e-3
    hi = compound(ion, target, mean + de)
    d_ex, d_beta = (np.asarray(hi.excitation) - cn.excitation) / de, (np.asarray(hi.beta) - cn.beta) / de
    return {"depth": z, "beam_energy": mean, "beam_sigma": sigma, "excitation": np.asarray(cn.excitation),
            "excitation_sigma": np.abs(d_ex) * sigma, "beta": np.asarray(cn.beta),
            "beta_sigma": np.abs(d_beta) * sigma, "compound": cn.nuclide, "q_value": cn.q_value}


def separation_energy(parent: tuple, particle: str) -> float:
    """Energy (MeV) needed to take ``particle`` out of the nucleus (Z, A) in its ground state, from the masses."""
    pz, pa = parse_nuclide(LIGHT[particle])
    d = data.nuclide((parent[0] - pz, parent[1] - pa)).nuclear_mass_mev
    return d + data.nuclide(LIGHT[particle]).nuclear_mass_mev - data.nuclide(parent).nuclear_mass_mev


def emission_barrier(parent: tuple, particle: str) -> float:
    """Coulomb barrier (MeV) for a charged particle leaving the nucleus (Z, A); 0 for a neutron."""
    pz, pa = parse_nuclide(LIGHT[particle])
    if pz == 0:
        return 0.0
    dz, da = parent[0] - pz, parent[1] - pa
    return pz * dz * E2_MEV_FM / (EMISSION_R0 * (pa ** (1 / 3) + da ** (1 / 3)))


def default_temperature(excitation_mev, a: int):
    """T = √(E*/a), a = A / 8 MeV⁻¹, MeV (at least 0.1 MeV)."""
    return np.sqrt(np.maximum(np.asarray(excitation_mev, dtype=float), 0.08 * a / LEVEL_DENSITY_DIVISOR)
                   * LEVEL_DENSITY_DIVISOR / a)


def likely_xn(beam: str, target: str, beam_energy: float) -> int:
    """How many neutrons the compound nucleus most likely boils off, as a first guess: neutrons leave while E*
    exceeds the next neutron's separation energy plus its mean kinetic energy 2T. 0 when not even one can."""
    cn = compound(beam, target, beam_energy)
    za, e, k = parse_nuclide(cn.nuclide), float(cn.excitation), 0
    while k < 12:
        need = separation_energy(za, "n") + 2 * float(default_temperature(e, za[1]))
        if e < need:
            break
        e, za, k = e - need, (za[0], za[1] - 1), k + 1
    return k


def grazing_l(beam: str, target: str, beam_energy: float) -> float:
    """The largest angular momentum (ħ) that fuses in a sharp cut-off picture: ℓ(ℓ+1)ħ² = 2μR²(E_cm − V_B) at the
    interaction radius R of :class:`~physim.nuclear.rutherford.Rutherford`; 0 below the barrier."""
    from .levels import HBARC_MEV_FM
    from .rutherford import Rutherford

    r = Rutherford(beam, target, beam_energy)
    over = r.e_cm - r.coulomb_barrier_cm
    if over <= 0:
        return 0.0
    mb, mt = data.nuclide(beam).nuclear_mass_mev, data.nuclide(target).nuclear_mass_mev
    x = 2 * mb * mt / (mb + mt) * r.interaction_radius**2 * over / HBARC_MEV_FM**2
    return (-1 + math.sqrt(1 + 4 * x)) / 2


# ---------------------------------------------------------------------------------------------------------------
# The channel table: the seam between Stage 1 and the fusion and evaporation codes


@dataclass
class ResidueChannel:
    """One exit channel as the event generator uses it.

    Filled from the setup (:meth:`ChannelTable.from_experiment`) or by another code
    (:meth:`ChannelTable.from_cross_sections`): only ``particles`` and ``cross_section_mb`` are needed, the rest
    falls back to defaults."""

    #: The particles, one entry each, in emission order: ("n", "n", "n", "n").
    particles: tuple
    cross_section_mb: float
    residue: str = ""
    #: The channel as the setup writes it ("4n"); :func:`particles_label` if not given.
    label: str = ""
    temperature_mev: Optional[float] = None
    entry_energy_mev: Optional[float] = None
    entry_spin: Optional[float] = None
    #: Where the cross section comes from ("setup", "PACE4", "physim fusion + evaporation", ...).
    source: str = ""


@dataclass
class ChannelTable:
    """The exit channels of a fusion-evaporation reaction with their cross sections: the input Stage 1 needs.

    The setup's ``[[reaction.channels]]`` are one way to fill it; :mod:`physim.nuclear.fusion` and
    :mod:`physim.nuclear.evaporation` fill the same object (:meth:`from_cross_sections`, :meth:`scaled_to`), and
    :meth:`to_setup` writes any table back as ``[[reaction.channels]]``, so the setup file keeps its shape whichever
    code made the numbers."""

    beam: str
    target: str
    channels: list = field(default_factory=list)
    #: Defaults for the channels that give none (None: from E*).
    temperature_mev: Optional[float] = None
    alignment: float = ALIGNMENT
    side_feeding: float = SIDE_FEEDING

    @property
    def compound(self) -> str:
        b, t = data.nuclide(self.beam), data.nuclide(self.target)
        return data.nuclide((b.Z + t.Z, b.A + t.A)).name

    @property
    def total_mb(self) -> float:
        return float(sum(ch.cross_section_mb for ch in self.channels))

    @property
    def labels(self) -> list:
        return [ch.label for ch in self.channels]

    def shares(self) -> np.ndarray:
        """σ_i / Σσ, the probability that an event is in each channel."""
        s = np.array([ch.cross_section_mb for ch in self.channels], dtype=float)
        return s / s.sum()

    def compound_at(self, beam_energy) -> CompoundState:
        return compound(self.beam, self.target, beam_energy)

    def index(self, label: str) -> int:
        """The channel by its label ("4n") or its residue ("158Er")."""
        for i, ch in enumerate(self.channels):
            if label in (ch.label, ch.residue, particles_label(ch.particles)):
                return i
        raise KeyError(f"no channel {label!r}; the channels are {', '.join(self.labels)}")

    @classmethod
    def from_experiment(cls, experiment) -> ChannelTable:
        """The setup's ``[[reaction.channels]]``."""
        from .experiment import _q, main_nuclide
        from .rates import beam_ion

        fe = experiment.fusion
        if fe is None:
            raise ValueError("the setup is not a fusion-evaporation reaction ([reaction] type = "
                             "\"fusion-evaporation\")")
        table = cls(beam_ion(experiment), main_nuclide(experiment.target),
                    temperature_mev=_q(fe.temperature).to("MeV") if fe.temperature is not None else None,
                    alignment=fe.alignment if fe.alignment is not None else ALIGNMENT,
                    side_feeding=_q(fe.side_feeding).to("%") / 100 if fe.side_feeding is not None else SIDE_FEEDING)
        for ch in fe.channels:
            parts = parse_particles(ch.particles)
            table.channels.append(ResidueChannel(
                parts, ch.cross_section_mb, ch.residue_name, ch.particles,
                _q(ch.temperature).to("MeV") if ch.temperature is not None else None,
                _q(ch.entry_energy).to("MeV") if ch.entry_energy is not None else None,
                ch.entry_spin, ch.source or "setup"))
        return table

    @classmethod
    def from_cross_sections(cls, beam: str, target: str, cross_sections: dict, source: str,
                            **defaults) -> ChannelTable:
        """A table from another code's results: ``cross_sections`` maps a channel's particles ("4n", or a tuple
        such as ("n", "n", "n", "n")) to its cross section in mb. Channels of zero cross section are left out;
        ``defaults`` are :class:`ChannelTable` fields (temperature_mev, alignment, side_feeding)."""
        table = cls(data.nuclide(beam).name, data.nuclide(target).name, **defaults)
        cn = data.nuclide(table.compound)
        for key, sigma in cross_sections.items():
            if sigma <= 0:
                continue
            parts = tuple(key) if not isinstance(key, str) else parse_particles(key)
            residue = data.nuclide(residue_of((cn.Z, cn.A), parts)).name
            table.channels.append(ResidueChannel(parts, float(sigma), residue,
                                                 key if isinstance(key, str) else particles_label(parts),
                                                 source=source))
        return table

    def scaled_to(self, total_mb: float, source: Optional[str] = None) -> ChannelTable:
        """The same shares, scaled to a total (the fusion cross section of Stage 2)."""
        k = total_mb / self.total_mb
        return replace(self, channels=[replace(ch, cross_section_mb=ch.cross_section_mb * k,
                                               source=source or ch.source) for ch in self.channels])

    def to_setup(self) -> list:
        """The table as ``[[reaction.channels]]`` sections, quantities with their units."""
        out = []
        for ch in self.channels:
            d = {"residue": ch.residue, "particles": ch.label or particles_label(ch.particles),
                 "cross_section": f"{ch.cross_section_mb:.10g} mb"}
            if ch.temperature_mev is not None:
                d["temperature"] = f"{ch.temperature_mev:.6g} MeV"
            if ch.entry_energy_mev is not None:
                d["entry_energy"] = f"{ch.entry_energy_mev:.6g} MeV"
            if ch.entry_spin is not None:
                d["entry_spin"] = ch.entry_spin
            if ch.source and ch.source != "setup":
                d["source"] = ch.source
            out.append(d)
        return out

    def entry_spin(self, i: int, beam_energy: float) -> float:
        """Where channel i's cascade starts, ħ: the table's, or (2/3) ℓ_max of the sharp cut-off (the mean spin of
        a triangular distribution) less what the particles take away (:data:`SPIN_REMOVED`), rounded down to the
        residue's kind of spin (whole for even A, half-whole for odd A)."""
        ch = self.channels[i]
        if ch.entry_spin is not None:
            return float(ch.entry_spin)
        j = 2.0 / 3.0 * grazing_l(self.beam, self.target, beam_energy) - sum(SPIN_REMOVED[p] for p in ch.particles)
        odd = data.nuclide(ch.residue).A % 2
        j = max(j, 0.5 * odd)
        return math.floor(j - 0.5 * odd) + 0.5 * odd

    def threshold(self, i: int) -> float:
        """The least E* (MeV) channel i needs: the separation energies of its particles and, for charged ones,
        their Coulomb barriers."""
        cn = data.nuclide(self.compound)
        za, need = (cn.Z, cn.A), 0.0
        for p in self.channels[i].particles:
            need += separation_energy(za, p) + emission_barrier(za, p)
            pz, pa = parse_nuclide(LIGHT[p])
            za = (za[0] - pz, za[1] - pa)
        return need


# ---------------------------------------------------------------------------------------------------------------
# Events


def _isotropic(rng, n: int) -> np.ndarray:
    cos = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2 * math.pi, n)
    sin = np.sqrt(1 - cos * cos)
    return np.stack([sin * np.cos(phi), sin * np.sin(phi), cos], axis=-1)


def _boost(energy: np.ndarray, momentum: np.ndarray, beta: np.ndarray) -> tuple:
    """Four-momenta (E, p) from a frame moving with velocity ``beta`` (units of c) to the laboratory."""
    b2 = np.einsum("ij,ij->i", beta, beta)
    g = 1.0 / np.sqrt(np.maximum(1.0 - b2, 1e-300))
    bp = np.einsum("ij,ij->i", beta, momentum)
    k = np.where(b2 > 0, (g - 1.0) * bp / np.where(b2 > 0, b2, 1.0), 0.0)
    p = momentum + (k + g * energy)[:, None] * beta
    return g * (energy + bp), p


def _two_body_momentum(m: np.ndarray, m1: np.ndarray, m2: np.ndarray) -> np.ndarray:
    """Momentum of either product of the decay of mass m into m1 + m2 at rest."""
    x = np.maximum((m * m - (m1 + m2) ** 2) * (m * m - (m1 - m2) ** 2), 0.0)
    return np.sqrt(x) / (2 * m)


#: The per-event columns of :class:`ResidueEvents`, with units.
COLUMNS = {
    "event": "event number",
    "channel": "index into ResidueEvents.channels",
    "depth": "depth of the reaction below the target's front face, along its normal, mg/cm2",
    "beam_energy": "beam energy at the reaction, MeV",
    "excitation": "excitation energy of the compound nucleus, MeV",
    "beta_cn": "speed of the compound nucleus, units of c",
    "beta": "speed of the residue, units of c",
    "theta": "lab polar angle of the residue, degrees",
    "phi": "lab azimuth of the residue, degrees",
    "energy": "kinetic energy of the residue, MeV",
    "entry_excitation": "excitation energy the particles leave the residue, MeV",
    "entry_spin": "spin the residue's cascade starts from, ħ",
    "x": "where on the target the reaction happened, across the beam, mm",
    "y": "the same, vertically, mm",
    "weight": "rate this event stands for, per second",
}
#: The per-particle columns (:attr:`ResidueEvents.particles`).
PARTICLE_COLUMNS = {
    "event": "event number",
    "kind": "which light particle: index into LIGHT ('n', 'p', 'd', 't', 'a')",
    "step": "1 for the first particle the event emits, 2 for the second, ...",
    "energy": "kinetic energy in the laboratory, MeV",
    "energy_cm": "kinetic energy in the frame of the compound nucleus's first emitter (relative motion), MeV",
    "theta": "lab polar angle, degrees",
    "phi": "lab azimuth, degrees",
}


@dataclass
class ResidueEvents:
    """Simulated fusion-evaporation events: one row per event in :attr:`columns` (see :data:`COLUMNS`), one per
    emitted particle in :attr:`particles` (see :data:`PARTICLE_COLUMNS`). Each event stands for the same rate, so
    the residues of a channel per second are its events × that rate."""

    columns: dict
    particles: dict
    table: ChannelTable
    n_events: int
    seed: int
    beam_time_s: float
    #: Beam particles per second × target nuclei per cm² (the luminosity, cm⁻² s⁻¹).
    luminosity: float
    #: β_CN for the nominal beam energy at the middle of the target, as the Doppler correction assumes it.
    beta_nominal: float
    notes: list = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.columns["event"])

    def __getitem__(self, name: str) -> np.ndarray:
        return self.columns[name]

    @property
    def channels(self) -> list:
        return self.table.labels

    def select(self, channel: Optional[str] = None) -> np.ndarray:
        m = np.ones(len(self), dtype=bool)
        if channel is not None:
            m &= self.columns["channel"] == self.table.index(channel)
        return m

    def rate(self, channel: Optional[str] = None) -> tuple:
        """(residues per second, statistical standard error)."""
        w = self.columns["weight"][self.select(channel)]
        return float(w.sum()), float(math.sqrt(np.sum(w**2)))

    def counts(self, channel: Optional[str] = None) -> float:
        """Residues made in the planned beam time."""
        return self.rate(channel)[0] * self.beam_time_s

    def per_channel(self) -> dict:
        """{channel label: (residues per second, statistical error)}."""
        return {label: self.rate(label) for label in self.channels}

    def expected(self, channel: Optional[str] = None) -> float:
        """What the table says: σ × luminosity, residues per second."""
        sigma = (self.table.total_mb if channel is None
                 else self.table.channels[self.table.index(channel)].cross_section_mb)
        return sigma * 1e-27 * self.luminosity


def _rest_energy(kind: str) -> float:
    return data.nuclide(LIGHT[kind]).nuclear_mass_mev


#: How often a channel's chain of emissions is drawn again for an event whose residue keeps enough excitation
#: energy to emit one more particle.
CHAIN_ATTEMPTS = 20


def _needs(compound_za: tuple, particles: tuple) -> tuple:
    """What each emission of a channel needs at least (its separation energy and barrier), and the residue."""
    need, zz = [], compound_za
    for kind in particles:
        need.append(separation_energy(zz, kind) + emission_barrier(zz, kind))
        pz, pa = parse_nuclide(LIGHT[kind])
        zz = (zz[0] - pz, zz[1] - pa)
    return need, zz


def _chain(rng, ch: ResidueChannel, table: ChannelTable, compound_za: tuple, need: list, velocity: np.ndarray,
           mass: np.ndarray, e_star: np.ndarray, limit: Optional[float] = None) -> dict:
    """One chain of emissions of channel ``ch`` for each compound nucleus given (velocity in c, mass with E*, E*):
    the residue's velocity, mass and E*, whether it ran short of energy, and per particle (kind, lab kinetic
    energy, relative kinetic energy ε, lab direction). With ``limit``, the last particle takes at least what
    leaves the residue below ``limit`` (drawn from the spectrum's tail above it; ``forced`` says where)."""
    velocity, mass, e_star = velocity.copy(), mass.copy(), e_star.copy()
    n = len(e_star)
    za, particles = compound_za, []
    short = np.zeros(n, dtype=bool)
    forced = np.zeros(n, dtype=bool)
    for step, kind in enumerate(ch.particles, start=1):
        pz, pa = parse_nuclide(LIGHT[kind])
        d_za = (za[0] - pz, za[1] - pa)
        m_x, m_d = _rest_energy(kind), data.nuclide(d_za).nuclear_mass_mev
        available = e_star - separation_energy(za, kind)
        # Each particle leaves behind what the ones after it need, so a channel the nucleus can afford completes.
        cap = available - sum(need[step:])
        barrier = emission_barrier(za, kind)
        temp = ch.temperature_mev if ch.temperature_mev is not None else table.temperature_mev
        t = np.full(n, temp) if temp is not None else default_temperature(e_star, za[1])
        # ε − V from x e^(−x/T): the sum of two exponential draws; redrawn while it takes more than it may.
        above = rng.exponential(1.0, (n, 2)).sum(axis=1) * t
        for _ in range(30):
            over = barrier + above > cap
            if not over.any():
                break
            above[over] = rng.exponential(1.0, (int(over.sum()), 2)).sum(axis=1) * t[over]
        # Still above the cap after 30 draws, or a cap below the barrier (tunnelling): the most it may take. Where
        # the nucleus cannot afford the channel at all (cap < 0) the particle leaves at rest and the residue is left
        # below its ground state by what is missing; energy and momentum stay exact.
        eps = np.maximum(np.minimum(barrier + above, np.maximum(cap, 0.0)), 0.0)
        if limit is not None and step == len(ch.particles):
            # x e^(−x/T) above c = floor − V: x = c + y, y from (c + y) e^(−y/T), which is an exponential with
            # probability c / (c + T) and the sum of two otherwise. Exact, and never below the floor.
            floor = available - limit * (1 - 1e-9)
            forced = floor > eps
            c = np.maximum(floor[forced] - barrier, 0.0)
            tf = t[forced]
            one = rng.uniform(0.0, 1.0, len(c)) < c / (c + tf)
            y = np.where(one, rng.exponential(1.0, len(c)), rng.exponential(1.0, (len(c), 2)).sum(axis=1)) * tf
            eps[forced] = np.minimum(barrier + c + y, np.maximum(cap[forced], 0.0))
        left = available - eps
        if step == len(ch.particles):
            short = cap < 0
        m1 = np.full(n, m_x)
        m2 = m_d + left
        p = _two_body_momentum(mass, m1, m2)
        u = _isotropic(rng, n)
        e1, p1 = _boost(np.sqrt(p * p + m1 * m1), p[:, None] * u, velocity)
        e2, p2 = _boost(np.sqrt(p * p + m2 * m2), -p[:, None] * u, velocity)
        velocity = p2 / e2[:, None]
        mass, e_star = m2, left
        norm = np.linalg.norm(p1, axis=1)
        particles.append((kind, e1 - m1, eps, p1 / np.where(norm > 0, norm, 1.0)[:, None]))
        za = d_za
    return {"velocity": velocity, "mass": mass, "e_star": e_star, "short": short, "particles": particles,
            "forced": forced if limit is not None else np.zeros(n, dtype=bool)}


def simulate_residues(experiment, events: int = 100_000, seed: int = 1) -> ResidueEvents:
    """Generate ``events`` fusion-evaporation events for a setup with ``[reaction] type = "fusion-evaporation"``.
    The same setup and seed give the same events, bit for bit."""
    from .rates import _q, beam_ion, spot_sigma_mm, stack

    if events < 1:
        raise ValueError("events must be at least 1")
    table = ChannelTable.from_experiment(experiment)
    rng = np.random.default_rng(int(seed))
    layers = stack(experiment)
    target_za = parse_nuclide(table.target)
    atoms = layers[0].nuclides.get(target_za, 0.0)
    pps = experiment.beam.particles_per_second
    luminosity = pps * atoms
    profile = compound_profile(experiment)
    n = int(events)

    # -- where and at which beam energy ------------------------------------------------------------------------
    channel = rng.choice(len(table.channels), size=n, p=table.shares())
    depth = rng.uniform(0.0, layers[0].thickness, n)
    e_mean = np.interp(depth, profile["depth"], profile["beam_energy"])
    e_sig = np.interp(depth, profile["depth"], profile["beam_sigma"])
    beam_energy = np.maximum(e_mean + e_sig * rng.standard_normal(n), 1e-6)
    spot = spot_sigma_mm(experiment)
    x, y = spot * rng.standard_normal(n), spot * rng.standard_normal(n)
    cn = compound(beam_ion(experiment), table.target, beam_energy)
    cn_za = parse_nuclide(cn.nuclide)
    excitation = np.asarray(cn.excitation, dtype=float)

    # -- the evaporation, channel by channel, particle by particle --------------------------------------------------
    velocity = np.zeros((n, 3))
    velocity[:, 2] = cn.beta
    mass = data.nuclide(cn.nuclide).nuclear_mass_mev + excitation    # the emitter's mass with its E*
    e_star = excitation.copy()
    entry_spin = np.zeros(n)
    p_event, p_kind, p_step, p_energy, p_ecm, p_dir = [], [], [], [], [], []
    kinds = list(LIGHT)
    short, hot = [0] * len(table.channels), [0] * len(table.channels)
    mid = float(np.interp(layers[0].thickness / 2, profile["depth"], profile["beam_energy"]))
    # One stream of numbers per channel, drawn in channel order: the same seed gives the same events.
    for i, ch in enumerate(table.channels):
        rows = np.flatnonzero(channel == i)
        entry_spin[rows] = table.entry_spin(i, mid)
        need, zz = _needs(cn_za, ch.particles)
        # The residue must not be able to emit another particle: its E* ends below its lowest threshold.
        limit = min(separation_energy(zz, k) + emission_barrier(zz, k) for k in ("n", "p", "a"))
        pending = rows
        for attempt in range(CHAIN_ATTEMPTS):
            if not len(pending):
                break
            last = attempt == CHAIN_ATTEMPTS - 1
            out = _chain(rng, ch, table, cn_za, need, velocity[pending], mass[pending], e_star[pending],
                         limit if last else None)
            ok = out["e_star"] < limit if not last else np.ones(len(pending), dtype=bool)
            if last:
                hot[i] += int(np.count_nonzero(out["forced"]))
            done = pending[ok]
            velocity[done], mass[done], e_star[done] = out["velocity"][ok], out["mass"][ok], out["e_star"][ok]
            short[i] += int(np.count_nonzero(out["short"][ok]))
            for step, (kind, energy, eps, dirs) in enumerate(out["particles"], start=1):
                p_event.append(done)
                p_kind.append(np.full(len(done), kinds.index(kind)))
                p_step.append(np.full(len(done), step))
                p_energy.append(energy[ok])
                p_ecm.append(eps[ok])
                p_dir.append(dirs[ok])
            pending = pending[~ok]

    beta = np.linalg.norm(velocity, axis=1)
    gam = 1.0 / np.sqrt(1.0 - beta * beta)
    residue_mass = mass                                   # with what E* is left
    weight = np.full(n, table.total_mb * 1e-27 * luminosity / n)
    cols = {"event": np.arange(n), "channel": channel, "depth": depth, "beam_energy": beam_energy,
            "excitation": excitation, "beta_cn": np.asarray(cn.beta, dtype=float), "beta": beta,
            "theta": np.degrees(np.arccos(np.clip(velocity[:, 2] / np.where(beta > 0, beta, 1.0), -1, 1))),
            "phi": np.degrees(np.arctan2(velocity[:, 1], velocity[:, 0])), "energy": residue_mass * (gam - 1.0),
            "entry_excitation": e_star, "entry_spin": entry_spin, "x": x, "y": y, "weight": weight,
            "_velocity": velocity}
    if p_event:
        dirs = np.concatenate(p_dir)
        parts = {"event": np.concatenate(p_event), "kind": np.concatenate(p_kind), "step": np.concatenate(p_step),
                 "energy": np.concatenate(p_energy), "energy_cm": np.concatenate(p_ecm),
                 "theta": np.degrees(np.arccos(np.clip(dirs[:, 2], -1, 1))),
                 "phi": np.degrees(np.arctan2(dirs[:, 1], dirs[:, 0]))}
        order = np.lexsort((parts["step"], parts["event"]))
        parts = {k: v[order] for k, v in parts.items()}
    else:
        parts = {k: np.zeros(0) for k in PARTICLE_COLUMNS}
    mid = float(np.interp(layers[0].thickness / 2, profile["depth"], profile["beam_energy"]))
    nominal = compound(beam_ion(experiment), table.target, mid)
    notes = [f"{n} events in {len(table.channels)} channels of {cn.nuclide}* (E* {excitation.mean():.4g} MeV on "
             f"average, Q = {cn.q_value:.4g} MeV); they stand for {table.total_mb * 1e-27 * luminosity:.4g} "
             "residues per second."]
    for i, k in enumerate(short):
        if k:
            notes.append(f"{k} events of channel {table.channels[i].label} had less excitation energy than the "
                         "channel needs; their residue is left below its ground state by what was missing.")
    for i, k in enumerate(hot):
        if k > 0.01 * np.count_nonzero(channel == i):
            notes.append(f"{k} of the {int(np.count_nonzero(channel == i))} events of channel "
                         f"{table.channels[i].label} kept enough excitation energy to emit another particle after "
                         f"{CHAIN_ATTEMPTS} draws of their chain; their last particle was drawn from the spectrum's "
                         "tail to carry the excess. The compound nucleus has more energy than this channel "
                         "usually carries away.")
    return ResidueEvents(cols, parts, table, n, int(seed), _q(experiment.run.beam_time).to("s"), luminosity,
                         float(nominal.beta), notes)


# ---------------------------------------------------------------------------------------------------------------
# The residues' level schemes and cascades


def residue_scheme(experiment, residue: str, max_energy_kev: float = 6000.0):
    """The level scheme of a residue: the setup's ``[levels.<residue>]``, else the local ENSDF copy (levels up to
    ``max_energy_kev``), else None."""
    if residue in experiment.levels:
        return experiment.levels[residue]
    return _ensdf_scheme(residue, max_energy_kev)


_ENSDF_SCHEMES: dict = {}


def _ensdf_scheme(residue: str, max_energy_kev: float):
    key = (residue, max_energy_kev)
    if key not in _ENSDF_SCHEMES:
        from . import ensdf
        from .levels import LevelScheme

        try:
            scheme = LevelScheme.from_adopted(ensdf.adopted(residue))
            scheme = scheme.select(max_energy_kev, connected_only=False) if scheme.levels else None
        except (ensdf.EnsdfMissing, ValueError, OSError):
            scheme = None
        _ENSDF_SCHEMES[key] = scheme
    return _ENSDF_SCHEMES[key]


def _multipoles(t, ji: Optional[float], jf: Optional[float]) -> tuple:
    """(L1, L2, δ) of a transition: from its multipolarity, else the lowest multipole the spins allow (dipole for
    ΔJ = 0)."""
    from .levels import components

    parts = components(t.multipolarity or "")
    if parts:
        l1 = int(parts[0][1])
        l2 = int(parts[1][1]) if len(parts) > 1 else 0
        delta = float(t.mixing.value) if (len(parts) > 1 and t.mixing is not None) else 0.0
        return l1, l2, delta
    if ji is None or jf is None:
        return 0, 0, 0.0
    return max(int(round(abs(ji - jf))), 1), 0, 0.0


def _gaussian_tensors(j: float, sigma_over_j: float) -> dict:
    """B_k = t_k0 / t_00 (k = 2, 4) of a level of spin J whose magnetic substates about the beam are populated as
    exp(−m²/2σ²), σ = ``sigma_over_j`` × J: the alignment of a spin perpendicular to the beam."""
    from .orientation import tensors_of

    if j < 1:
        return {2: 0.0, 4: 0.0}
    ms = np.array([-j + i for i in range(int(round(2 * j)) + 1)])
    sigma = max(sigma_over_j * j, 1e-3)
    p = np.exp(-ms**2 / (2 * sigma**2))
    rho = np.diag(p / p.sum()).astype(complex)
    t = tensors_of(rho, j, 4)
    return {k: float(t[(k, 0)].real / t[(0, 0)].real) if (k, 0) in t else 0.0 for k in (2, 4)}


@dataclass
class ResidueCascade:
    """How a residue's γ-ray cascade runs through its level scheme: where it enters (by spin) and the branching at
    each level, with the alignment carried down."""

    scheme: object
    alignment: float = ALIGNMENT
    side_feeding: float = SIDE_FEEDING

    def __post_init__(self):
        s = self.scheme
        energies = [lv.energy.value for lv in s.levels]
        self.ground = int(np.argmin(energies))
        #: The yrast levels: the lowest level of each known spin (not the ground state's own spin and below).
        best: dict = {}
        for i, lv in enumerate(s.levels):
            if lv.spin is None:
                continue
            if lv.spin not in best or energies[i] < energies[best[lv.spin]]:
                best[lv.spin] = i
        g_spin = s.levels[self.ground].spin or 0.0
        self.yrast = sorted(((j, i) for j, i in best.items() if i != self.ground and j > g_spin - 1e-9),
                            key=lambda x: x[0])
        #: Per level: [(final, total share, γ share, L1, L2, δ, energy MeV)], the shares summing to 1.
        self.branches: dict = {}
        for i in range(len(s.levels)):
            out = [t for t in s.transitions if t.initial == i and t.energy is not None]
            if not out:
                continue
            raw = [float(t.intensity.value) if t.intensity is not None else 0.0 for t in out]
            if not any(raw):
                raw = [1.0 if k == 0 else 0.0 for k in range(len(out))]
            totals = []
            for t, r in zip(out, raw):
                alpha = float(t.conversion.value) if t.conversion is not None else 0.0
                totals.append((r * (1 + alpha), 1.0 / (1 + alpha)))
            norm = sum(x for x, _ in totals)
            if norm <= 0:
                continue
            self.branches[i] = [(t.final, x / norm, g, *_multipoles(t, s.levels[i].spin, s.levels[t.final].spin),
                                 1e-3 * float(t.energy.value)) for t, (x, g) in zip(out, totals)]

    def feeding(self, entry_spin: float) -> dict:
        """{level: share} of the cascades a residue entering at ``entry_spin`` starts from: the top yrast level not
        above it, and the lower yrast levels sharing :attr:`side_feeding` evenly."""
        below = [i for j, i in self.yrast if j <= entry_spin + 1e-9]
        if not below:
            return {self.ground: 1.0}
        top, lower = below[-1], below[:-1]
        if not lower or self.side_feeding <= 0:
            return {top: 1.0}
        out = {top: 1.0 - self.side_feeding}
        for i in lower:
            out[i] = out.get(i, 0.0) + self.side_feeding / len(lower)
        return out

    def populations(self, entry_spin: float) -> dict:
        """{level: probability that a cascade passes through it} for one entry spin."""
        pop = dict(self.feeding(entry_spin))
        order = sorted(range(len(self.scheme.levels)), key=lambda i: -self.scheme.levels[i].energy.value)
        for i in order:
            p = pop.get(i, 0.0)
            if p <= 0:
                continue
            for final, share, *_ in self.branches.get(i, ()):
                pop[final] = pop.get(final, 0.0) + p * share
        return pop

    def line_yields(self, entry_spin: float) -> list:
        """The γ rays per residue of each transition: [(energy MeV, γ rays per residue, initial, final)],
        strongest first."""
        pop = self.populations(entry_spin)
        out = []
        for i, br in self.branches.items():
            for final, share, g, *_rest, e in br:
                y = pop.get(i, 0.0) * share * g
                if y > 0:
                    out.append((e, y, i, final))
        return sorted(out, key=lambda x: -x[1])

    def run(self, entry_spin: np.ndarray, rng) -> dict:
        """One cascade per entry: the γ rays, as arrays ``source`` (which entry), ``initial``, ``final``,
        ``energy`` (MeV) and the distribution coefficients ``a2`` and ``a4`` about the beam (W = 1 + a2 P2 +
        a4 P4)."""
        n = len(entry_spin)
        level = np.zeros(n, dtype=int)
        for j in np.unique(entry_spin):
            rows = np.flatnonzero(entry_spin == j)
            feed = self.feeding(float(j))
            keys = np.array(list(feed))
            level[rows] = keys[rng.choice(len(keys), size=len(rows), p=np.array(list(feed.values())))]
        spins = [lv.spin for lv in self.scheme.levels]
        b2, b4 = np.zeros(n), np.zeros(n)
        for lv in np.unique(level):
            m = level == lv
            t = _gaussian_tensors(spins[lv] or 0.0, self.alignment)
            b2[m], b4[m] = t[2], t[4]
        alive = np.ones(n, dtype=bool)
        out = {k: [] for k in ("source", "initial", "final", "energy", "a2", "a4")}
        for _ in range(4 * len(self.scheme.levels) + 4):
            alive &= np.isin(level, list(self.branches))
            alive &= level != self.ground
            idx = np.flatnonzero(alive)
            if not len(idx):
                break
            u = rng.uniform(0.0, 1.0, len(idx))
            g_draw = rng.uniform(0.0, 1.0, len(idx))
            now = level[idx].copy()  # each cascade takes one step per pass
            for lv in np.unique(now):
                here = now == lv
                sel, us, gs = idx[here], u[here], g_draw[here]
                br = self.branches[int(lv)]
                cum = np.cumsum([b[1] for b in br])
                pick = np.minimum(np.searchsorted(cum, us * cum[-1], side="right"), len(br) - 1)
                ji = spins[lv]
                for k, (final, _, g, l1, l2, delta, e) in enumerate(br):
                    rows = sel[pick == k]
                    if not len(rows):
                        continue
                    jf = spins[final]
                    known = ji is not None and jf is not None and l1 > 0
                    emitted = rows[gs[pick == k] < g]
                    if len(emitted):
                        a2 = angular.distribution_coefficient(2, ji, jf, l1, l2, delta) if known else 0.0
                        a4 = angular.distribution_coefficient(4, ji, jf, l1, l2, delta) if known else 0.0
                        out["source"].append(emitted)
                        out["initial"].append(np.full(len(emitted), int(lv)))
                        out["final"].append(np.full(len(emitted), final))
                        out["energy"].append(np.full(len(emitted), e))
                        out["a2"].append(b2[emitted] * a2)
                        out["a4"].append(b4[emitted] * a4)
                    if known:
                        u2, u4 = (_u(k_, ji, jf, l1, l2, delta) for k_ in (2, 4))
                        b2[rows] *= u2
                        b4[rows] *= u4
                    else:
                        b2[rows] = 0.0
                        b4[rows] = 0.0
                    level[rows] = final
        if not out["source"]:
            return {k: np.zeros(0, dtype=int if k in ("source", "initial", "final") else float) for k in out}
        res = {k: np.concatenate(v) for k, v in out.items()}
        order = np.argsort(res["source"], kind="stable")
        return {k: v[order] for k, v in res.items()}


def _u(k: int, ji: float, jf: float, l1: int, l2: int, delta: float) -> float:
    """U_k of an unobserved transition, mixed multipoles weighted by their intensities."""
    u = angular.u_coefficient(k, ji, jf, l1)
    if l2 and delta:
        u = (u + delta**2 * angular.u_coefficient(k, ji, jf, l2)) / (1 + delta**2)
    return u


# ---------------------------------------------------------------------------------------------------------------
# The γ rays in the crystals


@dataclass
class ResidueGammas:
    """γ rays of the residues that reached a crystal (see :data:`GAMMA_COLUMNS`). Singles: no particle is asked
    for; each row stands for a rate."""

    columns: dict
    events: ResidueEvents
    crystals: list
    seed: int
    #: β the Doppler correction assumes (the compound nucleus's, nominal beam energy, middle of the target).
    beta_nominal: float
    notes: list = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.columns["event"])

    def __getitem__(self, name: str) -> np.ndarray:
        return self.columns[name]

    @property
    def crystal_names(self) -> list:
        return [c.name for c in self.crystals]

    def select(self, crystal: Optional[str] = None, channel: Optional[str] = None, counted: bool = True,
               full_energy: bool = False) -> np.ndarray:
        """Boolean mask over the γ rays: of a crystal or whole γ-ray detector, of a channel, above threshold,
        leaving their full energy."""
        m = np.ones(len(self), dtype=bool)
        if crystal is not None:
            idx = [i for i, c in enumerate(self.crystals)
                   if c.name == crystal or (c.element is not None and c.name.rsplit(" ", 1)[0] == crystal)]
            if not idx:
                raise KeyError(f"no γ-ray detector or crystal named {crystal!r}; the crystals are "
                               f"{', '.join(self.crystal_names)}")
            m &= np.isin(self.columns["crystal"], idx)
        if channel is not None:
            m &= self.columns["channel"] == self.events.table.index(channel)
        if counted:
            m &= self.columns["counted"]
        if full_energy:
            m &= np.abs(self.columns["deposited"] - self.columns["energy_lab"]) < 1e-9
        return m

    def rate(self, crystal: Optional[str] = None, channel: Optional[str] = None) -> tuple:
        """(γ rays counted per second, statistical error)."""
        w = self.columns["weight"][self.select(crystal, channel)]
        return float(w.sum()), float(math.sqrt(np.sum(w**2)))

    def counts(self, crystal: Optional[str] = None, channel: Optional[str] = None) -> float:
        return self.rate(crystal, channel)[0] * self.events.beam_time_s

    def spectrum(self, crystal: Optional[str] = None, corrected: bool = True, channel: Optional[str] = None,
                 bins=400, range: Optional[tuple] = None) -> tuple:
        """Counts per bin in the planned beam time: the Doppler-corrected energies (or the measured ones). ``bins``
        is a number of bins or their edges (MeV)."""
        m = self.select(crystal, channel)
        x = self.columns["corrected" if corrected else "measured"][m]
        if range is None and np.ndim(bins) == 0:
            range = (0.0, max(float(x.max()) * 1.05 if len(x) else 1.0, 0.1))
        return np.histogram(x, bins=bins, range=range, weights=self.columns["weight"][m] * self.events.beam_time_s)


def doppler_correct(beta: float, crystal_dirs: np.ndarray, measured: np.ndarray) -> np.ndarray:
    """E₀ = E γ (1 − β cos θ), θ the crystal centre's angle to the beam: the emitter assumed along the beam."""
    g = 1.0 / math.sqrt(1.0 - beta * beta)
    return measured * g * (1.0 - beta * crystal_dirs[:, 2])


def simulate_residue_gammas(experiment, events: int = 100_000, seed: int = 1,
                            residue_events: Optional[ResidueEvents] = None, gammas_per_event: int = 4,
                            plain: bool = False) -> ResidueGammas:
    """The residues' γ rays in the setup's crystals. Each cascade is emitted ``gammas_per_event`` times over, each
    with a share of the event's rate (the rates are unchanged, the spectra smoother). The same setup and seed
    give the same γ rays, bit for bit."""
    from .gamma_events import _deposit, _first_crystal, _own_shares, _unit, crystals_of
    from .response import FWHM_PER_SIGMA

    if not experiment.gamma_detectors:
        raise ValueError("the setup has no γ-ray detectors")
    ev = residue_events if residue_events is not None else simulate_residues(experiment, events, seed)
    table = ev.table
    c = ev.columns
    rng = np.random.default_rng([int(seed), 75])
    crystals, bare = crystals_of(experiment, plain=True)
    centres = np.array([cr.centre for cr in crystals])
    normals = _unit(centres)
    radii = np.array([cr.radius for cr in crystals])
    k_rep = max(int(gammas_per_event), 1)

    parts, notes, skipped = [], list(ev.notes), []
    first = 0
    for i, ch in enumerate(table.channels):
        rows = np.flatnonzero(c["channel"] == i)
        scheme = residue_scheme(experiment, ch.residue)
        if scheme is None or not scheme.transitions:
            skipped.append(f"{ch.label} ({ch.residue})")
            continue
        cas = _cascade_of(scheme, table)
        rows = np.repeat(rows, k_rep)
        g = cas.run(c["entry_spin"][rows], rng)
        src = rows[g["source"]]
        cascade = g["source"] + first
        first += len(rows)
        m = len(src)
        if not m:
            continue
        # Directions in the residue's frame from W(θ) about the beam, by rejection.
        cos = np.empty(m)
        todo = np.arange(m)
        top = 1.0 + np.abs(g["a2"]) + np.abs(g["a4"])
        for _ in range(200):
            x = rng.uniform(-1.0, 1.0, len(todo))
            p2, p4 = (3 * x * x - 1) / 2, (35 * x**4 - 30 * x * x + 3) / 8
            w = 1.0 + g["a2"][todo] * p2 + g["a4"][todo] * p4
            ok = rng.uniform(0.0, 1.0, len(todo)) * top[todo] <= w
            cos[todo[ok]] = x[ok]
            todo = todo[~ok]
            if not len(todo):
                break
        cos[todo] = rng.uniform(-1.0, 1.0, len(todo))
        phi = rng.uniform(0.0, 2 * math.pi, m)
        sin = np.sqrt(1 - cos * cos)
        rest = np.stack([sin * np.cos(phi), sin * np.sin(phi), cos], axis=-1)
        vel = c["_velocity"][src]
        e_lab, p_lab = _boost(g["energy"], g["energy"][:, None] * rest, vel)
        lab = p_lab / e_lab[:, None]
        origin = np.zeros((m, 3))
        origin[:, 0], origin[:, 1] = c["x"][src], c["y"][src]
        which, reached = _first_crystal(origin, lab, centres, normals, radii)
        keep = np.flatnonzero(reached)
        parts.append({"event": src[keep], "cascade": cascade[keep], "channel": np.full(len(keep), i),
                      "initial": g["initial"][keep], "final": g["final"][keep], "crystal": which[keep],
                      "energy0": g["energy"][keep], "energy_lab": e_lab[keep], "lab": lab[keep],
                      "beta": c["beta"][src[keep]], "weight": c["weight"][src[keep]] / k_rep})
    if skipped:
        notes.append("No level scheme (in the setup or the local ENSDF copy) for " + ", ".join(skipped)
                     + ": their residues emit no γ rays here.")
    if not parts:
        cols = {k: np.zeros(0) for k in GAMMA_COLUMNS}
        return ResidueGammas(cols, ev, crystals, int(seed), ev.beta_nominal, notes)
    cols = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    lab = cols.pop("lab")
    n = len(cols["event"])
    which, energy_lab = cols["crystal"], cols["energy_lab"]
    deposited, measured = np.zeros(n), np.zeros(n)
    counted, keep = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    own = _own_shares(crystals)
    for k, cr in enumerate(crystals):
        m = which == k
        if not m.any():
            continue
        e = energy_lab[m]
        r = bare[k]
        cover = (sum(r.coverage) if cr.element is None else r.coverage[cr.element]) * own[k]
        p_int = r.transmission(e) * np.minimum(r.total_efficiency(e, cr.element) / cover / r.transmission(e), 1.0)
        interact = rng.uniform(0.0, 1.0, m.sum()) <= p_int
        dep, _ = _deposit(r, e, rng)
        meas = dep + rng.normal(0.0, 1.0, len(dep)) * r.fwhm(np.maximum(dep, 1e-3)) / FWHM_PER_SIGMA
        idx = np.flatnonzero(m)
        deposited[idx], measured[idx] = dep, meas
        counted[idx] = interact & (meas >= cr.threshold)
        keep[idx] = interact
    cols.update(deposited=deposited, measured=measured, counted=counted,
                theta=np.degrees(np.arccos(np.clip(lab[:, 2], -1, 1))), phi=np.degrees(np.arctan2(lab[:, 1], lab[:, 0])))
    cols = {k: v[keep] for k, v in cols.items()}
    cols["corrected"] = doppler_correct(ev.beta_nominal, normals[cols["crystal"]], cols["measured"])
    notes.append(f"{len(cols['event'])} γ rays interacted in the crystals ({k_rep} emissions of each cascade). "
                 f"Doppler correction for β = {ev.beta_nominal:.5f} along the beam. Each level is aligned with "
                 f"σ/J = {table.alignment:g}; {100 * table.side_feeding:.0f} % side feeding. No summing, add-back "
                 "or lifetimes.")
    return ResidueGammas(cols, ev, crystals, int(seed), ev.beta_nominal, notes)


def _cascade_of(scheme, table: ChannelTable) -> ResidueCascade:
    key = (id(scheme), table.alignment, table.side_feeding)
    cached = _CASCADES.get(key)
    if cached is None or cached.scheme is not scheme:
        cached = _CASCADES[key] = ResidueCascade(scheme, table.alignment, table.side_feeding)
    return cached


_CASCADES: dict = {}


# ---------------------------------------------------------------------------------------------------------------
# Rates and the Plan table


def residue_table(experiment, lines: int = 3) -> dict:
    """The Plan tab's residues: per channel its cross section, share, residues per second, residues in the run,
    the beam time for the counts wanted in its strongest line (detected in the γ array's full-energy peaks), and
    its strongest lines (energy, γ rays per residue, full-energy counts per second in the array). Also the compound
    nucleus at the middle of the target, and notes (channels closed, no level scheme, placeholder values).

    The line rates use each detector's full-energy-peak efficiency at the line's energy, without the alignment's
    angular distribution or the Doppler shift; the events (:func:`simulate_residue_gammas`) have both."""
    from .rates import _q, beam_ion, stack

    table = ChannelTable.from_experiment(experiment)
    layers = stack(experiment)
    atoms = layers[0].nuclides.get(parse_nuclide(table.target), 0.0)
    pps = experiment.beam.particles_per_second
    lum = pps * atoms
    profile = compound_profile(experiment)
    mid = float(np.interp(layers[0].thickness / 2, profile["depth"], profile["beam_energy"]))
    cn = compound(beam_ion(experiment), table.target, mid)
    beam_time = _q(experiment.run.beam_time).to("s")
    wanted = experiment.run.counts_wanted
    effs: dict = {}

    def efficiency(e_mev: float) -> float:
        key = round(e_mev, 6)
        if key not in effs:
            effs[key] = float(sum(gd.peak_efficiency(e_mev, experiment) for gd in experiment.gamma_detectors))
        return effs[key]

    rows, notes = [], []
    for i, ch in enumerate(table.channels):
        rate = ch.cross_section_mb * 1e-27 * lum
        spin = table.entry_spin(i, mid)
        scheme = residue_scheme(experiment, ch.residue)
        strongest = []
        if scheme is not None and scheme.transitions:
            for e, y, a, b in _cascade_of(scheme, table).line_yields(spin)[:lines]:
                det = rate * y * efficiency(e) if experiment.gamma_detectors else None
                strongest.append({"energy_kev": 1e3 * e, "per_residue": y, "detected_per_s": det,
                                  "label": f"{_jpi(scheme, a)} → {_jpi(scheme, b)}"})
        need = table.threshold(i)
        top = strongest[0]["detected_per_s"] if strongest else None
        rows.append({"channel": ch.label, "residue": ch.residue, "cross_section_mb": ch.cross_section_mb,
                     "share": ch.cross_section_mb / table.total_mb, "rate_per_s": rate,
                     "residues_in_run": rate * beam_time, "entry_spin": spin, "threshold_mev": need,
                     "open": need < cn.excitation, "source": ch.source, "lines": strongest,
                     "scheme": None if scheme is None else ("setup" if ch.residue in experiment.levels else "ENSDF"),
                     "beam_time_s": (wanted / top if (wanted and top) else (math.inf if wanted else None))})
        if need >= cn.excitation:
            notes.append(f"{ch.label} ({ch.residue}) needs E* of at least {need:.3g} MeV (separation energies and "
                         f"barriers); the compound nucleus has {cn.excitation:.3g} MeV at the middle of the target.")
        if scheme is None:
            notes.append(f"No level scheme for {ch.residue}: add [levels.{ch.residue}] or a local ENSDF copy to see "
                         "its γ rays.")
        if (ch.source or "").lower() == "placeholder":
            notes.append(f"The cross section of {ch.label} is a placeholder: enter one from PACE4 or a measurement.")
    return {"compound": cn.nuclide, "beam_energy_mev": mid, "e_cm_mev": cn.e_cm, "q_value_mev": cn.q_value,
            "excitation_mev": cn.excitation, "excitation_range_mev": (float(profile["excitation"].min()),
                                                                      float(profile["excitation"].max())),
            "beta": cn.beta, "beta_range": (float(profile["beta"].min()), float(profile["beta"].max())),
            "total_mb": table.total_mb, "luminosity": lum, "rows": rows, "notes": notes,
            "counts_wanted": wanted}


def _jpi(scheme, i: int) -> str:
    lv = scheme.levels[i]
    return lv.label if (lv.spin is not None or lv.jpi) else f"{lv.energy.value:g} keV"


def warnings(experiment, table: Optional[dict] = None) -> list:
    """Plain-word warnings about a fusion-evaporation setup: below the Coulomb barrier, channels the excitation
    energy cannot open, placeholder cross sections. ``table`` is the setup's :func:`residue_table`, if at hand."""
    from .rates import beam_ion
    from .rutherford import Rutherford
    from .experiment import main_nuclide

    out = []
    r = Rutherford(beam_ion(experiment), main_nuclide(experiment.target), experiment.beam.energy_mev)
    if r.e_cm <= r.coulomb_barrier_cm:
        out.append(f"The beam is below the Coulomb barrier ({r.e_cm:.3g} MeV in the CM against "
                   f"{r.coulomb_barrier_cm:.3g} MeV): fusion is orders of magnitude weaker there than any table "
                   "above the barrier says.")
    out += (table if table is not None else residue_table(experiment))["notes"]
    return out


__all__ = ["COLUMNS", "GAMMA_COLUMNS", "LIGHT", "PARTICLE_COLUMNS", "ChannelTable", "CompoundState", "ResidueCascade",
           "ResidueChannel", "ResidueEvents", "ResidueGammas", "compound", "compound_profile", "default_temperature",
           "doppler_correct", "emission_barrier", "grazing_l", "parse_particles", "particles_label", "residue_of",
           "residue_scheme", "residue_table", "separation_energy", "simulate_residue_gammas", "simulate_residues",
           "warnings"]
