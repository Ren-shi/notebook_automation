"""Count rates, beam time and expected peaks: the numbers a beam-time proposal needs.

The rate into a detector (or one strip) from one kind of nucleus in the target is

    counts/s = (beam particles/s) × (nuclei/cm²) × ∫ dσ/dΩ dΩ,

integrated over the detector face, for the scattered beam (the ejectile) and for the recoiling target nucleus, and
averaged through the target's thickness (the beam slows down, and Rutherford's cross section grows as 1/E²).
A backing contributes in the same way, so a carbon backing shows up as its own peak::

    from physim.nuclear import Experiment
    from physim.nuclear.rates import Rates

    r = Rates(Experiment.example("alpha_on_gold"))
    r.per_detector()             # {detector: counts/s}
    r.counts_in_run("A45")       # counts in the planned beam time
    r.beam_time_for("A45")       # seconds to collect the counts the setup asks for
    r.peaks("A45")               # where the peaks sit and how wide they are
    r.warnings()                 # detectors counting too fast, Rutherford's limits, layout problems

Angles in degrees, energies in MeV, rates per second. Only elastic (Rutherford) scattering is modelled so far.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, fields, replace
from functools import lru_cache
from typing import Optional

import numpy as np

from . import data
from .detectors import Array, angles
from .kinematics import TwoBody, lab_points
from .quantity import Quantity
from .rutherford import Rutherford
from .stopping import Stopping

#: Detector count rate above which :meth:`Rates.warnings` warns, per second: silicon detectors and their
#: electronics start to suffer from pile-up and dead time at a few thousand counts per second.
MAX_RATE = 5000.0
#: Smallest lab angle counted, degrees: Rutherford's cross section diverges at 0°.
THETA_FLOOR = 0.5
#: Quadrature order over the face of a detector partly hidden behind another. The shadow's edge makes the integrand
#: jump, so the rate converges only slowly with the order: 48 is within 0.5% of the exact visible solid angle of a
#: disc half hidden by a smaller one (12 can be off by a few percent).
SHADOW_ORDER = 48
#: Lowest energy (MeV) a particle must leave the reaction with to be followed, unless a detector threshold is lower
#: still. Rutherford's cross section diverges for distant collisions, which send recoils out at almost 90° with
#: almost no energy; they never leave the target or reach a measurable signal.
ENERGY_FLOOR = 0.01
#: FWHM = 2√(2 ln 2) σ.
FWHM_PER_SIGMA = 2.3548200450309493
#: mb → cm².
MB_CM2 = 1e-27
#: Longest path through the target, as a multiple of its thickness, for tracks almost in its plane.
MAX_PATH_FACTOR = 1e3

PARTICLES = ("ejectile", "recoil")


def _q(x) -> Optional[Quantity]:
    if x is None:
        return None
    return Quantity.parse(x) if isinstance(x, str) else x


# ---------------------------------------------------------------------------------------------------------------
# The target stack and the reaction channels


@dataclass
class StackLayer:
    """One layer the beam crosses: the target, then its backing (downstream)."""

    name: str
    material: data.Material
    #: Thickness along the layer's normal, mg/cm².
    thickness: float
    #: Nuclei per cm² of each (Z, A).
    nuclides: dict


@dataclass
class Channel:
    """Scattering on one kind of nucleus in one layer: elastic, or (with ``excitation``) Coulomb excitation of a
    state of the target nucleus or of the projectile."""

    layer: int
    layer_name: str
    nuclide: data.Nuclide
    atoms_per_cm2: float
    #: The excited state (:class:`~physim.nuclear.experiment.Excitation`) for a Coulomb-excitation channel.
    excitation: object = None

    @property
    def label(self) -> str:
        base = f"{self.nuclide.name} ({self.layer_name})"
        if self.excitation is None:
            return base
        e = self.excitation
        return f"{base}, {e.excite} excited to {e.energy_mev * 1e3:g} keV"

    @property
    def kinematics_args(self) -> dict:
        """Keyword arguments for :class:`~physim.nuclear.kinematics.TwoBody`: the excitation energy and where."""
        if self.excitation is None:
            return {}
        return {"excitation_mev": self.excitation.energy_mev,
                "excite": "recoil" if self.excitation.excite == "target" else "ejectile"}


def stack(experiment) -> list:
    """The target and its backing as :class:`StackLayer` objects, in the order the beam crosses them."""
    t = experiment.target
    layers = [(t.material_data(), _q(t.thickness), "target")]
    if t.backing is not None:
        layers.append((t.backing.material_data(), _q(t.backing.thickness), "backing"))
    # Contaminants (backlog item 74) on the downstream face: thin layers of their own, each scattering elastically.
    for name, th in t.contaminants or []:
        layers.append((data.material(name), _q(th), f"{name} contaminant"))
    return [StackLayer(name, mat, mat.areal_density_mg_cm2(th), mat.atoms_per_cm2(th)) for mat, th, name in layers]


def recoil_stack(experiment, distance_um: Optional[float] = None) -> list:
    """The recoil's path for a lifetime measurement, as :class:`~physim.nuclear.transport.Segment` objects from the
    front of the target: the layers of :func:`stack` (target, backing, contaminants), then for a plunger (``method``
    "rdds") a vacuum gap of ``distance_um`` plus the plunger's offset and the foil after it (named "stopper" or
    "degrader" by its kind), and for a differential plunger a second gap and the second foil ("degrader2").

    ``distance_um`` is a target-to-stopper distance as the setup gives it (``Lifetime.distances``, without the
    offset); a plunger needs one. Without a lifetime section, or for DSAM, the stack alone (the backing slows the
    recoil). The beam's path is :func:`stack`, unchanged: the stopper is in the recoil's way only."""
    from .transport import Segment

    out = [Segment(lay.material, lay.thickness, name=lay.name) for lay in stack(experiment)]
    lt = getattr(experiment, "lifetime", None)
    if lt is None or lt.method != "rdds":
        return out
    if distance_um is None:
        raise ValueError("a plunger's recoil path needs a target-to-stopper distance (distance_um)")
    if distance_um < 0:
        raise ValueError(f"the distance must not be negative, not {distance_um} um")
    out.append(Segment(None, gap_um=float(distance_um) + lt.offset_um, name="gap"))
    if lt.stopper is not None:
        mat = lt.stopper.material_data()
        name = "stopper" if lt.stopper.kind == "stopper" else "degrader"
        out.append(Segment(mat, mat.areal_density_mg_cm2(_q(lt.stopper.thickness)), name=name))
    if lt.degrader2 is not None and lt.stopper is not None and lt.stopper.kind == "differential":
        mat = lt.degrader2.material_data()
        out.append(Segment(None, gap_um=lt.degrader2.gap_um, name="gap2"))
        out.append(Segment(mat, mat.areal_density_mg_cm2(_q(lt.degrader2.thickness)), name="degrader2"))
    return out


def channels(layers: list, experiment=None) -> list:
    """One elastic :class:`Channel` per nuclide per layer and, for a Coulomb-excitation setup, the excitation
    channels: of the main nuclide of the target (target excitation), or of the beam by every nuclide (projectile
    excitation)."""
    out = [Channel(i, lay.name, data.nuclide(key), n) for i, lay in enumerate(layers)
           for key, n in lay.nuclides.items() if n > 0]
    exc = getattr(experiment, "excitation", None) if experiment is not None else None
    if exc is not None:
        if exc.excite == "target":
            key, n = max(layers[0].nuclides.items(), key=lambda kv: kv[1])
            out.append(Channel(0, layers[0].name, data.nuclide(key), n, exc))
        else:
            out += [Channel(i, lay.name, data.nuclide(key), n, exc) for i, lay in enumerate(layers)
                    for key, n in lay.nuclides.items() if n > 0]
    return out


def rutherford_factor(e_cm: float, e_star: float) -> float:
    """Alder and Winther's Rutherford factor of an inelastic cross section over the elastic Rutherford cross section
    at the initial CM energy ``e_cm`` (MeV), at the same CM angle: v_i/v_f = √(E_cm / (E_cm − E*)) (see
    :meth:`physim.nuclear.coulex.Coulex.rutherford_cm`). 1 at or below the threshold, where nothing is excited."""
    return math.sqrt(e_cm / (e_cm - e_star)) if e_cm > e_star > 0 else 1.0


class _CoulexReaction:
    """Coulomb excitation at one beam energy, shaped like :class:`~physim.nuclear.rutherford.Rutherford`:
    ``kinematics`` (with Q = −E*) and ``cross_section_cm``. The excitation probability P(θ) is that of
    ``coulex``, at this energy (:meth:`DepthProbability.at`). The cross section is P times Alder and Winther's
    inelastic Rutherford factor, the elastic Rutherford cross section at this energy times v_i/v_f
    (:func:`rutherford_factor`), as in :meth:`physim.nuclear.coulex.Coulex.cross_section_cm` and in GOSIA; with a
    level scheme, v_f is that of the ``[reaction]`` state, whose kinematics the channel follows."""

    def __init__(self, beam: str, channel: Channel, energy: float, coulex):
        self.kinematics = TwoBody(beam, channel.nuclide.name, energy, **channel.kinematics_args)
        self._ruth = Rutherford(beam, channel.nuclide.name, energy)
        self._factor = rutherford_factor(self._ruth.e_cm, channel.excitation.energy_mev)
        #: Anything with ``probability(theta_cm)``: see :func:`excitation_probability`.
        self._coulex = coulex

    def cross_section_cm(self, theta_cm):
        return (self._factor * np.asarray(self._ruth.cross_section_cm(theta_cm))
                * np.asarray(self._coulex.probability(theta_cm)))


@lru_cache(maxsize=64)
def _coulex_at(beam: str, nuclide: str, energy: float, excite: str, e_star: float, multipolarity: str,
               b_up: float):
    from .coulex import Coulex

    return Coulex(beam, nuclide, energy, excite=excite, energy=e_star, multipolarity=multipolarity, b_up=b_up)


def layer_energy(experiment, layer: int, layers: Optional[list] = None) -> float:
    """The mean beam energy (MeV) in the middle of the part of layer ``layer`` the beam crosses: the layer's
    middle, or halfway to where the beam stops when it stops inside the layer (a thick backing); 0 when the beam
    stops before it."""
    layers = layers or stack(experiment)
    lay = layers[layer]
    e = float(beam_energy_at(experiment, layer, [lay.thickness / 2], layers)[0])
    if e > 0:
        return e
    entry = float(beam_energy_at(experiment, layer, [0.0], layers)[0])
    if entry <= 0:
        return 0.0
    c = math.cos(math.radians(tilt_deg(experiment)))
    reached = min(float(stopping(beam_ion(experiment), lay.material).range(entry)) * c, lay.thickness)
    return float(beam_energy_at(experiment, layer, [reached / 2], layers)[0])


def coulex_for(experiment, channel: Channel, layers: Optional[list] = None, energy: Optional[float] = None):
    """The :class:`~physim.nuclear.coulex.Coulex` of a channel, at beam energy ``energy`` (MeV), by default the
    energy in the middle of its layer (of the part of it the beam crosses, when the beam stops inside:
    :func:`layer_energy`). Raises ValueError when the beam stops before the layer, or (from
    :class:`~physim.nuclear.coulex.Coulex`) when the energy is below the excitation threshold."""
    layers = layers or stack(experiment)
    if energy is None:
        energy = layer_energy(experiment, channel.layer, layers)
        if not energy > 0:
            raise ValueError(f"the beam stops before {layers[channel.layer].name}: no Coulomb excitation there")
    e = channel.excitation
    return _coulex_at(beam_ion(experiment), channel.nuclide.name, round(float(energy), 9), e.excite, e.energy_mev,
                      e.multipolarity, e.b_up_e2fm)


class SchemeProbability:
    """The probability that a collision excites the nucleus at all, against the CM angle, from a solution of its
    level scheme (:func:`physim.nuclear.gamma.cascade_excitation`): the populations of all its excited levels
    after the collision, summed, between the angles of the solution's table (and falling to zero at 0°)."""

    def __init__(self, solution):
        grid, states = solution.table()
        self.solution = solution
        self._grid = np.r_[0.0, grid]
        self._total = np.r_[0.0, [sum(s.direct.values()) for s in states]]

    def probability(self, theta_cm):
        return np.interp(np.asarray(theta_cm, dtype=float), self._grid, self._total)


def excitation_probability(experiment, channel: Channel, layers: Optional[list] = None,
                           energy: Optional[float] = None):
    """What gives the excitation probability P(θ) of a Coulomb-excitation channel at beam energy ``energy`` (MeV;
    by default the middle of its layer), an object with ``probability(theta_cm)``: with a level scheme that the
    events follow (:func:`physim.nuclear.gamma.cascade_excitation`), the summed populations of its solution for
    this channel (:class:`SchemeProbability`: to all orders, or first order if the setup says so); otherwise the
    first-order :class:`~physim.nuclear.coulex.Coulex` of the ``[reaction]`` state.

    This is P at one energy. Through a thick layer the beam slows down and P changes (by orders of magnitude for
    inverse kinematics on a carbon target): the rates, the events and the analysis take P through the depth from
    :func:`depth_probability`, which solves it here at a few energies."""
    from .gamma import cascade_excitation

    solution = cascade_excitation(experiment, channel=channel, energy=energy)
    if solution is None:
        return coulex_for(experiment, channel, layers, energy)
    kept = solution.__dict__.get("_scheme_probability")
    if kept is None:
        kept = solution.__dict__["_scheme_probability"] = SchemeProbability(solution)
    return kept


#: Relative spread of the beam energy through a layer (its energy loss over its entry energy) below which the
#: excitation probability is taken at the middle of the layer alone, as for a thin target.
THIN_SPAN = 2e-3
#: (relative energy spread through the layer, energies at which P is solved): P is interpolated between them in
#: ln P, a polynomial through Chebyshev points. Six points follow a 38 % energy loss (62Ni on 1 mg/cm² of carbon,
#: where P changes 150-fold) to 3e-4 at the worst angle; fewer do for thinner layers.
DEPTH_NODES = ((0.03, 3), (0.12, 4), (0.25, 5), (0.4, 6), (0.6, 7), (math.inf, 8))
#: Depth points (Gauss–Legendre) of :meth:`DepthProbability.probability`, and per slice of :meth:`DepthProbability.
#: slices`.
DEPTH_AVERAGE_POINTS = 24
SLICE_POINTS = 4
#: Largest change of ln P across one slice of a layer in the event generator (:meth:`DepthProbability.slice_count`),
#: and the most slices. The rates of the events are right with any number of slices (P is averaged over each);
#: the slices give the excitations their depth, and so the inelastic peaks their energies: with n slices the
#: generator keeps about 1 − 1/n² of the shift of the peak (coulex_ni58, where P doubles through the target: 2
#: slices keep 70 % of its 0.24 MeV, the 5 it gets now 96 %).
SLICE_LOG_STEP = 0.15
MAX_SLICES = 16
#: A layer through which ln P changes by less than this stays one slice: the depth the excitations move to then
#: shifts their peaks by less than 2.5 % of the layer's own energy spread, and the events of such a setup (thin,
#: or nearly so: coulex_ni62_pt194, a 2.7 % energy loss) stay those of before.
SLICE_MIN_LOG_CHANGE = 0.3
_DEPTH_KEPT: dict = {}
_DEPTH_KEPT_MOST = 32


class _NoExcitation:
    """P = 0 at every angle: a beam energy at or just above the excitation threshold."""

    def probability(self, theta_cm):
        return 0.0 if np.ndim(theta_cm) == 0 else np.zeros(np.shape(theta_cm))


#: CM angles (degrees) on which an interpolated P is tabulated, then interpolated (cubic) between them.
INTERPOLATION_GRID = np.linspace(0.0, 180.0, 721)


class _InterpolatedProbability:
    """P(θ) at one beam energy, interpolated in ln P between solutions at nearby energies (``nodes``, with the
    Lagrange ``weights`` of this energy: :class:`DepthProbability`), tabulated every 0.25° and interpolated
    between (four-point cubic, as :meth:`physim.nuclear.coulex.Coulex.probability`)."""

    def __init__(self, nodes: list, weights: np.ndarray, energy: float):
        self.energy = energy
        grid = INTERPOLATION_GRID
        ps = np.array([np.asarray(n.probability(grid), dtype=float) for n in nodes])
        w = np.asarray(weights, dtype=float)[:, None]
        positive = np.all(ps > 0, axis=0)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            log = np.exp((w * np.log(np.where(ps > 0, ps, 1.0))).sum(axis=0))
        self._values = np.where(positive, log, np.maximum((w * ps).sum(axis=0), 0.0))

    def probability(self, theta_cm):
        th = np.atleast_1d(np.asarray(theta_cm, dtype=float))
        grid, values = INTERPOLATION_GRID, self._values
        h = grid[1] - grid[0]
        t = np.clip(np.nan_to_num(th), 0.0, 180.0)
        i = np.clip(np.floor(t / h).astype(int) - 1, 0, len(grid) - 4)
        x = (t - grid[i]) / h
        y0, y1, y2, y3 = (values[i + k] for k in range(4))
        out = (-y0 * (x - 1) * (x - 2) * (x - 3) / 6 + y1 * x * (x - 2) * (x - 3) / 2
               - y2 * x * (x - 1) * (x - 3) / 2 + y3 * x * (x - 1) * (x - 2) / 6)
        out = np.where((th <= 0) | np.isnan(th), 0.0, np.maximum(out, 0.0))
        return float(out[0]) if np.ndim(theta_cm) == 0 else out.reshape(np.shape(theta_cm))


class DepthProbability:
    """The excitation probability of a Coulomb-excitation channel through the depth of its layer.

    The beam slows down in the layer, and the excitation probability P(θ) changes with its energy. P is solved
    (:func:`excitation_probability`) at a few beam energies between where the beam enters the layer and where it
    leaves (or stops: :data:`DEPTH_NODES`) and interpolated between them in ln P; a thin layer keeps the one
    solution at its middle (:data:`THIN_SPAN`), so its results are those of the mid-layer P.

    - :meth:`at` gives P at one beam energy (the rates take it at each of their depth points).
    - :meth:`probability` gives P averaged over the depth at given CM angles, each depth weighted by the elastic
      Rutherford cross section there (∝ 1/E_cm² at a fixed CM angle): ⟨P⟩ = ∫ σ_R(E(z), θ) P(E(z), θ) dz /
      ∫ σ_R(E(z), θ) dz.
    - :meth:`ratio` gives σ_Coulex/σ_Rutherford(elastic) integrated over the layer, ⟨P · v_i/v_f⟩: the
      Coulomb-excitation cross section carries Alder and Winther's inelastic Rutherford factor
      (:func:`rutherford_factor`). An analysis that divides the excited counts by the elastic ones compares them
      with this.
    - :meth:`slices` cuts the layer into slices, each with P (or the ratio) averaged over it (the event
      generator).

    Get one with :func:`depth_probability`."""

    def __init__(self, experiment, channel: Channel, layers: Optional[list] = None):
        layers = layers or stack(experiment)
        self.channel = channel
        self.beam = beam_ion(experiment)
        lay = layers[channel.layer]
        self.thickness = lay.thickness
        entry = float(beam_energy_at(experiment, channel.layer, [0.0], layers)[0])
        if not entry > 0:
            raise ValueError(f"the beam stops before {lay.name}: no Coulomb excitation there")
        c = math.cos(math.radians(tilt_deg(experiment)))
        #: Depth (mg/cm², along the normal) down to which the beam crosses the layer, and its energy there.
        self.reached = min(float(stopping(self.beam, lay.material).range(entry)) * c, lay.thickness)
        self._experiment, self._layers = experiment, layers
        self.e_front = entry
        self.e_back = float(beam_energy_at(experiment, channel.layer, [self.reached], layers)[0])
        e_star = channel.excitation.energy_mev
        m1, m2 = data.nuclide(self.beam).atomic_mass_u, channel.nuclide.atomic_mass_u
        # Below the threshold (E_cm = E*) nothing is excited; just above it P is vanishingly small.
        threshold = e_star * (m1 + m2) / m2
        self.e_low = max(self.e_back, min(1.05 * threshold, entry))
        span = (self.e_front - self.e_low) / self.e_front
        if span < THIN_SPAN:
            self.energies = np.array([layer_energy(experiment, channel.layer, layers)])
        else:
            k = next(n for limit, n in DEPTH_NODES if span < limit)
            x = np.cos(np.pi * (np.arange(k) + 0.5) / k)
            self.energies = (self.e_front + self.e_low) / 2 + (self.e_front - self.e_low) / 2 * x
        self.nodes = [excitation_probability(experiment, channel, layers, float(e)) for e in self.energies]
        self._at: dict = {}

    @property
    def thin(self) -> bool:
        """Whether P is the one at the middle of the layer (a thin layer)."""
        return len(self.nodes) == 1

    def at(self, energy: float):
        """P at beam energy ``energy`` (MeV): an object with ``probability(theta_cm)``. Clamped to the energies
        of the layer; zero below the excitation threshold."""
        if self.thin:
            return self.nodes[0]
        e = min(max(float(energy), self.e_low), self.e_front)
        if self.closed(energy):
            return _NoExcitation()
        key = round(e, 9)
        kept = self._at.get(key)
        if kept is None:
            if len(self._at) > 256:
                self._at.clear()
            kept = self._at[key] = self.between(self.nodes, e)
        return kept

    def closed(self, energy: float) -> bool:
        """Whether beam energy ``energy`` (MeV) is near or below the excitation threshold, where :meth:`at` is 0."""
        return not self.thin and float(energy) < self.e_low and self.e_low > self.e_back

    def lagrange(self, energies) -> np.ndarray:
        """Weights of the solutions at :attr:`energies` for the beam energies ``energies`` (MeV, clamped to the
        layer's), one row per solution: a quantity solved at :attr:`energies` is interpolated between them by
        their weighted sum (the polynomial through the Chebyshev points)."""
        e = np.atleast_1d(np.asarray(energies, dtype=float))
        if self.thin:
            return np.ones((1, len(e)))
        e = np.clip(e, self.e_low, self.e_front)
        es = self.energies
        return np.array([np.prod([(e - es[m]) / (es[j] - es[m]) for m in range(len(es)) if m != j], axis=0)
                         for j in range(len(es))])

    def between(self, nodes: list, energy: float):
        """What ``nodes`` (objects with ``probability(theta_cm)``, one per energy of :attr:`energies`) give at beam
        energy ``energy`` (MeV, inside the layer's), interpolated in ln P between them: an object with
        ``probability(theta_cm)``. :meth:`at` is this for the whole excitation's probability."""
        if self.thin:
            return nodes[0]
        return _InterpolatedProbability(nodes, self.lagrange(energy)[:, 0], float(energy))

    def _rutherford_weight(self, energies) -> tuple:
        """The elastic Rutherford cross section at a fixed CM angle at each beam energy, relative (∝ 1/E_cm²), and
        the inelastic Rutherford factor v_i/v_f there (:func:`rutherford_factor`)."""
        weight, factor = [], []
        e_star = self.channel.excitation.energy_mev
        for e in energies:
            if e > 0:
                r = Rutherford(self.beam, self.channel.nuclide.name, float(e))
                weight.append(float(r.cross_section_cm(90.0)))
                factor.append(rutherford_factor(r.e_cm, e_star))
            else:
                weight.append(0.0)
                factor.append(1.0)
        return np.array(weight), np.array(factor)

    def _average(self, z_lo: float, z_hi: float, theta_cm, points: int, ratio: bool = False):
        x, w = np.polynomial.legendre.leggauss(points)
        zs = z_lo + (x + 1) / 2 * (z_hi - z_lo)
        es = beam_energy_at(self._experiment, self.channel.layer, zs, self._layers)
        weight, factor = self._rutherford_weight(es)
        wr = w * weight
        if wr.sum() <= 0:
            return np.zeros(np.shape(theta_cm)) if np.ndim(theta_cm) else 0.0
        ps = np.array([np.asarray(self.at(e).probability(theta_cm), dtype=float) for e in es])
        out = np.tensordot(wr * factor if ratio else wr, ps, axes=1) / wr.sum()
        return float(out) if np.ndim(theta_cm) == 0 else out

    def _mid_factor(self) -> float:
        """v_i/v_f at the energy of a thin layer's one solution."""
        r = Rutherford(self.beam, self.channel.nuclide.name, float(self.energies[0]))
        return rutherford_factor(r.e_cm, self.channel.excitation.energy_mev)

    def probability(self, theta_cm, points: int = DEPTH_AVERAGE_POINTS):
        """⟨P⟩: P at CM angles ``theta_cm`` (degrees, of the scattered beam) averaged over the depth of the layer
        the beam crosses, weighted by the elastic Rutherford cross section at each depth. What elastic scattering
        loses to the excitation (1 − ⟨P⟩); the ratio of the Coulomb-excitation cross section to the elastic one is
        :meth:`ratio`."""
        if self.thin:
            return self.nodes[0].probability(theta_cm)
        return self._average(0.0, self.reached, theta_cm, points)

    def ratio(self, theta_cm, points: int = DEPTH_AVERAGE_POINTS):
        """σ_Coulex / σ_Rutherford(**elastic**) at CM angles ``theta_cm`` (degrees, of the scattered beam), each
        integrated over the depth of the layer: ⟨P · v_i/v_f⟩ with the elastic Rutherford weight of each depth.
        The Coulomb-excitation cross section carries Alder and Winther's inelastic Rutherford factor, v_i/v_f
        above the elastic one (:meth:`physim.nuclear.coulex.Coulex.rutherford_cm`); this is what an analysis that
        divides the excited counts by the elastic ones (at the same CM angle) compares them with. GOSIA's ratios
        divide by its own, inelastic, Rutherford cross section instead (``"R_G"``), which is larger."""
        if self.thin:
            p = self.nodes[0].probability(theta_cm)
            return self._mid_factor() * (float(p) if np.ndim(theta_cm) == 0 else np.asarray(p, dtype=float))
        return self._average(0.0, self.reached, theta_cm, points, ratio=True)

    def slice_count(self, theta_cm: float = 170.0) -> int:
        """How many slices :meth:`slices` should cut the layer into: enough that ln P changes by at most
        :data:`SLICE_LOG_STEP` across each (at CM angle ``theta_cm``), at most :data:`MAX_SLICES`; one when it
        changes by less than :data:`SLICE_MIN_LOG_CHANGE` through the whole layer."""
        if self.thin:
            return 1
        p = [float(self.at(e).probability(theta_cm)) for e in (self.e_front, self.e_low)]
        if min(p) <= 0:
            return MAX_SLICES
        change = abs(math.log(p[0] / p[1]))
        if change < SLICE_MIN_LOG_CHANGE:
            return 1
        return int(min(max(math.ceil(change / SLICE_LOG_STEP), 1), MAX_SLICES))

    def slices(self, n: int, theta_cm, ratio: bool = False) -> list:
        """The layer cut into ``n`` slices of equal thickness: [(depth of the slice's front, its thickness, P at
        ``theta_cm`` averaged over the slice as :meth:`probability` averages it)]. With ``ratio``, the average is
        :meth:`ratio`'s instead (P × v_i/v_f, the weight of the excitation events on the elastic Rutherford cross
        section). The part of the layer beyond where the beam stops has P = 0."""
        edges = np.linspace(0.0, self.thickness, n + 1)
        mid = self._mid_factor() if self.thin and ratio else 1.0
        out = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            if self.thin:
                p = mid * np.asarray(self.nodes[0].probability(theta_cm), dtype=float)
            elif lo >= self.reached:
                p = np.zeros(np.shape(theta_cm))
            else:
                p = np.asarray(self._average(lo, min(hi, self.reached), theta_cm, SLICE_POINTS, ratio),
                               dtype=float)
            out.append((float(lo), float(hi - lo), p))
        return out


def _depth_key(experiment, channel: Channel, layers: list) -> tuple:
    from .gamma import excitation_order

    exc = channel.excitation
    scheme = experiment.levels.get("target" if exc.excite == "target" else "beam")
    beam = experiment.beam
    return (beam_ion(experiment), beam.energy_mev, tilt_deg(experiment), channel.layer,
            tuple((lay.material, lay.thickness) for lay in layers[:channel.layer + 1]), channel.nuclide.name,
            exc.excite, exc.energy_mev, exc.multipolarity, exc.b_up_e2fm, getattr(exc, "emission", None),
            excitation_order(experiment),
            None if scheme is None else json.dumps(scheme.to_dict(), sort_keys=True, default=str))


def depth_probability(experiment, channel: Channel, layers: Optional[list] = None) -> DepthProbability:
    """The :class:`DepthProbability` of a Coulomb-excitation channel: its excitation probability through the
    depth of its layer. Kept by the setup's physics, so the rates, the events and the analysis of one setup share
    it. Raises ValueError when the beam stops before the layer.

    For an analysis (σ_Coulex/σ_Rutherford(elastic) per detector): ``depth_probability(experiment, channel).ratio(
    theta_cm)`` is that ratio at those CM angles, integrated over the layer (P × v_i/v_f averaged with the
    elastic Rutherford weight of each depth); ``.probability(theta_cm)`` is ⟨P⟩ alone."""
    layers = layers or stack(experiment)
    key = _depth_key(experiment, channel, layers)
    kept = _DEPTH_KEPT.get(key)
    if kept is None:
        kept = DepthProbability(experiment, channel, layers)
        if len(_DEPTH_KEPT) >= _DEPTH_KEPT_MOST:
            _DEPTH_KEPT.pop(next(iter(_DEPTH_KEPT)))
        _DEPTH_KEPT[key] = kept
    return kept


class _ElasticReaction:
    """Elastic scattering beside Coulomb excitation, shaped like :class:`~physim.nuclear.rutherford.Rutherford`:
    a collision either excites the nucleus or leaves it as it was, so elastic scattering is the Rutherford cross
    section times 1 − P(θ), P being the excitation probability of the channel that excites it (``coulex``, at this
    energy). Elastic and excited together are Rutherford scattering, one particle per collision (to (v_i/v_f − 1) P:
    the excitation carries Alder and Winther's inelastic Rutherford factor, :func:`rutherford_factor`)."""

    def __init__(self, beam: str, channel: Channel, energy: float, coulex):
        self._ruth = Rutherford(beam, channel.nuclide.name, energy)
        self.kinematics = self._ruth.kinematics
        self._coulex = coulex

    def cross_section_cm(self, theta_cm):
        p = np.clip(np.asarray(self._coulex.probability(theta_cm), dtype=float), 0.0, 1.0)
        return np.asarray(self._ruth.cross_section_cm(theta_cm)) * (1.0 - p)


def excited_partner(experiment, channel: Channel, layers: Optional[list] = None) -> Optional[Channel]:
    """The Coulomb-excitation channel on the same nuclide in the same layer as an elastic ``channel``: the one
    whose excitation takes its collisions away from elastic scattering. None when the reaction excites nothing
    there."""
    if channel.excitation is not None or getattr(experiment, "excitation", None) is None:
        return None
    for ch in channels(layers or stack(experiment), experiment):
        if ch.excitation is not None and ch.layer == channel.layer and ch.nuclide.name == channel.nuclide.name:
            return ch
    return None


def reaction_at(experiment, channel: Channel, energy: float, layers: Optional[list] = None):
    """The reaction of a channel at one beam energy: an object with ``kinematics`` and ``cross_section_cm``."""
    beam = beam_ion(experiment)
    if channel.excitation is None:
        partner = excited_partner(experiment, channel, layers)
        if partner is None:
            return Rutherford(beam, channel.nuclide.name, energy)
        return _ElasticReaction(beam, channel, energy, depth_probability(experiment, partner, layers).at(energy))
    return _CoulexReaction(beam, channel, energy, depth_probability(experiment, channel, layers).at(energy))


class _YieldReaction:
    """Scattering on a nuclide weighted by the γ rays of one transition per collision, shaped like
    :class:`~physim.nuclear.rutherford.Rutherford`: the Rutherford cross section at this energy times the yield
    at each CM angle (``grid``, ``per_collision``: from a level scheme's solution at the middle of the layer).
    The particle keeps the kinematics of elastic scattering."""

    def __init__(self, beam: str, nuclide: str, energy: float, grid, per_collision):
        self._ruth = Rutherford(beam, nuclide, energy)
        self.kinematics = self._ruth.kinematics
        self._grid, self._y = np.asarray(grid, dtype=float), np.asarray(per_collision, dtype=float)

    def cross_section_cm(self, theta_cm):
        return np.asarray(self._ruth.cross_section_cm(theta_cm)) * np.interp(np.asarray(theta_cm, dtype=float),
                                                                           self._grid, self._y)


def beam_ion(experiment) -> str:
    return data.nuclide((experiment.beam.Z, experiment.beam.A)).name


def tilt_deg(experiment) -> float:
    t = experiment.target.tilt
    return _q(t).to("deg") if t is not None else 0.0


def energy_sigma(experiment) -> float:
    """Standard deviation of the beam energy, MeV (the setup gives the FWHM, as an energy or a percentage)."""
    spread = _q(experiment.beam.energy_spread)
    if spread is None:
        return 0.0
    fwhm = spread.to("MeV") if spread.kind == "energy" else spread.canonical / 100 * experiment.beam.energy_mev
    return fwhm / FWHM_PER_SIGMA


def spot_sigma_mm(experiment) -> float:
    """Standard deviation of the beam spot in x and in y, mm (the setup gives the FWHM)."""
    spot = _q(experiment.beam.spot_size)
    return 0.0 if spot is None else spot.to("mm") / FWHM_PER_SIGMA


@lru_cache(maxsize=256)
def stopping(ion: str, material: data.Material) -> Stopping:
    """Cached :class:`~physim.nuclear.stopping.Stopping`."""
    return Stopping(ion, material)


@lru_cache(maxsize=256)
def _w_table(st: Stopping) -> np.ndarray:
    return np.asarray(st.transport_table()["w"])


def _after(st: Stopping, e, path) -> np.ndarray:
    """Mean energy after ``path`` (mg/cm²), element by element; 0 where the ion stops."""
    e, path = np.broadcast_arrays(np.asarray(e, dtype=float), np.asarray(path, dtype=float))
    left = np.asarray(st.range(e)) - path
    out = np.exp(np.interp(np.log(np.maximum(left, 1e-300)), np.log(st._range), np.log(st._e)))
    return np.where((left > 0) & (e > 0), np.where(path > 0, out, e), 0.0)


def _carry(st: Stopping, e: float, var: float, path: float, straggle: bool = True) -> tuple:
    """Mean energy and variance after ``path``: an incoming spread is carried by (S_out/S_in)², and straggling
    adds S_out² (W(E_in) − W(E_out)) (the same formula the event generator uses)."""
    if e <= 0 or path <= 0:
        return max(e, 0.0), var if e > 0 else 0.0
    out = float(_after(st, e, path))
    if out <= 0:
        return 0.0, 0.0
    s_in, s_out = st.stopping_power(e), st.stopping_power(out)
    var = var * (s_out / s_in) ** 2
    if straggle:
        w = _w_table(st)
        ln = np.log(st._e)
        var += s_out**2 * max(float(np.interp(math.log(e), ln, w) - np.interp(math.log(out), ln, w)), 0.0)
    return out, var


def beam_energy_at(experiment, layer: int, z, layers: Optional[list] = None) -> np.ndarray:
    """Mean beam energy (MeV) at depth ``z`` (mg/cm² along the normal) into layer ``layer``."""
    layers = layers or stack(experiment)
    ion = beam_ion(experiment)
    c = math.cos(math.radians(tilt_deg(experiment)))
    e = experiment.beam.energy_mev
    for lay in layers[:layer]:
        e = float(_after(stopping(ion, lay.material), e, lay.thickness / c))
    z = np.asarray(z, dtype=float)
    return _after(stopping(ion, layers[layer].material), np.full(z.shape, e), z / c)


def _exit_paths(experiment, layers: list, layer: int, z, dirs: np.ndarray) -> tuple:
    """For tracks leaving from depth ``z`` of ``layer`` along ``dirs``: (forward mask, [(layer index, path in
    mg/cm², forward pass)]). Forward tracks meet layers ``layer``, ``layer + 1``, ... and leave by the back face;
    backward ones meet ``layer``, ``layer − 1``, ... and leave by the front. Paths grow as 1/cos up to a cap."""
    cos_n = np.asarray(dirs) @ _normal(experiment)
    total = sum(lay.thickness for lay in layers)
    fwd = cos_n > 0
    with np.errstate(divide="ignore"):
        inv = 1.0 / np.abs(cos_n)
    z = np.asarray(z, dtype=float)
    steps = []
    for forward, js in ((True, range(layer, len(layers))), (False, range(layer, -1, -1))):
        for j in js:
            th = layers[j].thickness
            dz = ((th - z) if forward else z) if j == layer else np.full_like(cos_n, th)
            dz = np.broadcast_to(dz, np.broadcast(dz, cos_n).shape)
            with np.errstate(invalid="ignore"):
                path = np.where(dz > 0, np.minimum(dz * inv, MAX_PATH_FACTOR * total), 0.0)
            steps.append((j, path, forward))
    return fwd, steps


def exit_energy(experiment, layers: list, layer: int, z, ion: str, e, dirs: np.ndarray) -> np.ndarray:
    """Mean energy of particles of kind ``ion`` leaving the stack from depth ``z`` of ``layer`` along ``dirs``."""
    e = np.asarray(e, dtype=float)
    fwd, steps = _exit_paths(experiment, layers, layer, z, dirs)
    out = {True: e.copy(), False: e.copy()}
    for j, path, forward in steps:
        out[forward] = _after(stopping(ion, layers[j].material), out[forward], path)
    return np.where(fwd, out[True], out[False])


def _normal(experiment) -> np.ndarray:
    t = math.radians(tilt_deg(experiment))
    return np.array([math.sin(t), 0.0, math.cos(t)])


def energy_cut(experiment) -> float:
    """Lowest energy (MeV) a particle must leave the reaction with to be followed: the lowest detector threshold,
    but at least :data:`ENERGY_FLOOR`."""
    thresholds = [_q(d.threshold).to("MeV") if d.threshold is not None else 0.0 for d in experiment.detectors]
    return max(min(thresholds, default=0.0), ENERGY_FLOOR)


def cm_ranges(experiment, layers: list, channel: Channel, array: Optional[Array] = None,
              theta_floor: float = THETA_FLOOR) -> list:
    """The separate ranges of the ejectile's CM angle, [(smallest, largest), ...] in degrees and in rising order,
    for which the ejectile or the recoil can reach a detector with at least :func:`energy_cut`, over the energies
    the beam has in the channel's layer; empty if neither can.

    Each detector covers its own band of lab angles, so detectors on both sides of 90° leave a gap between their
    ranges. The recoils of distant collisions fill the lab angles just below 90°, where the cross section is
    largest; a range that held the gap would spend nearly every event there, on recoils that meet no detector."""
    array = array or Array.from_experiment(experiment)
    bands = [g.theta_range() for g in array]
    if not bands:
        return []
    spot = spot_sigma_mm(experiment)
    halo = experiment.beam
    # The beam reaches 3σ of its spot from the axis, or the halo's radius if it has one (backlog item 74).
    reach = max(3 * spot, _q(halo.halo_radius).to("mm") if halo.halo_fraction is not None and halo.halo_radius
                is not None and _q(halo.halo_fraction).value > 0 else 0.0)
    d_min = min(float(np.linalg.norm(g.centre)) for g in array)
    widen = math.degrees(reach / d_min) + 0.2
    bands = [(max(lo - widen, theta_floor), min(hi + widen, 180.0)) for lo, hi in bands]
    lay = layers[channel.layer]
    e_front, e_back = beam_energy_at(experiment, channel.layer, [0.0, lay.thickness], layers)
    de = 3 * energy_sigma(experiment)
    cut = energy_cut(experiment)
    grid = np.linspace(0.05, 180.0, 3600)
    mask = np.zeros(grid.shape, dtype=bool)
    for e in (e_front + de, max(e_back - de, 0.5 * e_back)):
        if e <= 0:
            continue
        try:
            tb = TwoBody(beam_ion(experiment), channel.nuclide.name, float(e), **channel.kinematics_args)
        except ValueError:
            continue
        for pt in (tb.at_cm(grid), tb.recoil_for(grid)):
            seen = np.zeros(grid.shape, dtype=bool)
            for lo, hi in bands:
                seen |= (pt.theta_lab >= lo) & (pt.theta_lab <= hi)
            mask |= seen & (pt.energy >= cut)
    if not mask.any():
        return []
    step = grid[1] - grid[0]
    edges = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(np.int8)))
    out: list = []
    for first, after in zip(edges[::2], edges[1::2]):
        lo, hi = max(float(grid[first] - step), 0.01), min(float(grid[after - 1] + step), 180.0)
        if out and lo <= out[-1][1] + step:  # one grid point apart: the same range
            out[-1] = (out[-1][0], hi)
        else:
            out.append((lo, hi))
    return out


def cm_acceptance(experiment, layers: list, channel: Channel, array: Optional[Array] = None,
                  theta_floor: float = THETA_FLOOR) -> Optional[tuple]:
    """(smallest, largest) CM angle of the ejectile, degrees, for which the ejectile or the recoil can reach a
    detector with at least :func:`energy_cut`, over the energies the beam has in the channel's layer; ``None`` if
    neither can. The detectors may cover only parts of it: see :func:`cm_ranges`."""
    ranges = cm_ranges(experiment, layers, channel, array, theta_floor)
    return (ranges[0][0], ranges[-1][1]) if ranges else None


# ---------------------------------------------------------------------------------------------------------------
# Rates


@dataclass
class RateRow:
    """Rate of one particle kind from one channel into one detector segment."""

    detector: str
    segment: tuple
    channel: str
    particle: str
    #: Counts per second.
    rate: float


@dataclass
class Peak:
    """An expected peak in a detector's measured-energy spectrum."""

    detector: str
    segment: Optional[tuple]
    channel: str
    particle: str
    #: 1, or 2 for the second (lower-energy) kinematic solution where the kinematics are double-valued.
    branch: int
    #: Counts per second in the peak.
    rate: float
    #: Mean measured energy, MeV, and its standard deviation.
    mean: float
    sigma: float
    #: Contributions to the width (standard deviations, MeV): "target thickness", "kinematic" (spread of angles
    #: over the face), "geometry" (the two together, as they combine), "beam energy spread", "straggling",
    #: "resolution". ``sigma`` adds geometry, energy spread, straggling and resolution in quadrature.
    components: dict = field(default_factory=dict)
    #: Mean energy reaching the detector face, MeV.
    energy_at_face: float = 0.0
    threshold: float = 0.0
    #: Fraction of the particles whose mean measured energy (before resolution smearing) is above the threshold.
    counted_fraction: float = 1.0

    @property
    def fwhm(self) -> float:
        return FWHM_PER_SIGMA * self.sigma

    @property
    def visible(self) -> bool:
        """The peak lies above the detector threshold (and the particles reach the detector)."""
        return self.mean > 0 and self.mean >= self.threshold


class Rates:
    """Analytic count rates and expected peaks for every detector of an experiment.

    ``depth_points`` Gauss–Legendre points average over the target's thickness and ``order`` × ``order`` points
    integrate over each detector segment; the defaults are accurate to well below 0.1% for the examples.

    ``only`` computes one detector alone (the others still cast their shadows): with a low ``order`` this is the
    quick estimate shown while a detector is dragged. ``previous`` is an earlier :class:`Rates`; the rows of every
    detector that has not changed, and that no other detector shadows before or after, are taken from it.
    """

    def __init__(self, experiment, depth_points: int = 4, order: int = 12, theta_floor: float = THETA_FLOOR,
                 min_energy: Optional[float] = None, only: Optional[str] = None, previous: Optional[Rates] = None):
        self.experiment = experiment
        self.array = Array.from_experiment(experiment)
        self.layers = stack(experiment)
        self.channels = channels(self.layers, experiment)
        self.beam = beam_ion(experiment)
        self.particles_per_second = experiment.beam.particles_per_second
        self.beam_time_s = _q(experiment.run.beam_time).to("s")
        self.counts_wanted = experiment.run.counts_wanted
        self.theta_floor = theta_floor
        #: Particles leaving the reaction with less energy than this (MeV) are not counted (see :func:`energy_cut`).
        self.min_energy = energy_cut(experiment) if min_energy is None else min_energy
        self._depth_points, self._order = depth_points, order
        #: Detectors partly hidden behind others: integrated on a finer grid, as the edge of the shadow cuts across
        #: the face (fully hidden ones simply count nothing).
        shadows = self.array.shadowing()
        self._partly_hidden = {behind for (_, behind), frac in shadows.items() if frac < 1 - 1e-6}
        self._shadowed = {behind for _, behind in shadows}
        self._only = only
        #: Everything but the detectors that decides a detector's rows, up to the beam's intensity: the rows are
        #: proportional to the particles per second, and the run conditions (beam time, counts wanted, dead time,
        #: coincidence window) do not enter them, so an edit of either reuses every row.
        beam = tuple((f.name, getattr(experiment.beam, f.name)) for f in fields(experiment.beam)
                     if f.name != "current")
        exc = experiment.excitation
        # The level scheme of the excited nucleus sets the excitation probability (:func:`excitation_probability`).
        scheme = experiment.levels.get("target" if exc.excite == "target" else "beam") if exc is not None else None
        self._physics = repr((beam, experiment.target, experiment.reaction, experiment.excitation, depth_points,
                              order, theta_floor, self.min_energy, scheme.to_dict() if scheme is not None else None))
        self._peak_cache: dict = {}
        #: Detectors whose rows were taken from ``previous``.
        self.reused = self._unchanged(previous)
        self.rows = self._compute(previous)

    def _unchanged(self, previous: Optional[Rates]) -> set:
        if previous is None or previous._only is not None or previous._physics != self._physics:
            return set()
        old = {g.name: d for g, d in zip(previous.array.geometries, previous.array._setup)}
        out = set()
        for g, d in zip(self.array.geometries, self.array._setup):
            if (d is not None and g.name in old and old[g.name] == d and g.name not in self._shadowed
                    and g.name not in previous._shadowed):
                out.add(g.name)
        return out

    def _depth_nodes(self, layer: int) -> tuple:
        x, w = np.polynomial.legendre.leggauss(self._depth_points)
        th = self.layers[layer].thickness
        return (x + 1) / 2 * th, w / 2

    def _order_for(self, g) -> int:
        """Quadrature order per segment: finer for a partly hidden detector, so the face as a whole is sampled about
        as finely as ``SHADOW_ORDER`` over one segment (its many segments already make the grid fine)."""
        if g.name not in self._partly_hidden:
            return self._order
        return max(self._order, math.ceil(SHADOW_ORDER / math.sqrt(len(g.segments))))

    def _compute(self, previous: Optional[Rates] = None) -> list:
        rows = []
        for g in self.array:
            if self._only is not None and g.name != self._only:
                continue
            if g.name in self.reused:
                scale = self.particles_per_second / previous.particles_per_second
                if scale == 1.0:
                    rows += [r for r in previous.rows if r.detector == g.name]
                    self._peak_cache.update({k: v for k, v in previous._peak_cache.items() if k[0] == g.name})
                else:  # another beam current: the same rows and peaks, in proportion
                    rows += [replace(r, rate=r.rate * scale) for r in previous.rows if r.detector == g.name]
                    self._peak_cache.update({k: [replace(pk, rate=pk.rate * scale) for pk in v]
                                             for k, v in previous._peak_cache.items() if k[0] == g.name})
                continue
            segs = g.segments
            parts = [g.directions(seg, self._order_for(g)) for seg in segs]
            dirs = np.concatenate([p[0] for p in parts])
            dom = np.concatenate([p[1] for p in parts])
            index = np.repeat(np.arange(len(segs)), [len(p[1]) for p in parts])
            theta = angles(dirs)[0]
            # Directions that meet another detector first never reach this one.
            dom = np.where((theta >= self.theta_floor) & self.array.visible(g.name, dirs), dom, 0.0)
            for ch in self.channels:
                zs, wz = self._depth_nodes(ch.layer)
                es = beam_energy_at(self.experiment, ch.layer, zs, self.layers)
                sums = {p: np.zeros(len(dirs)) for p in PARTICLES}
                reacs, ws = [], []
                for e, w in zip(es, wz):
                    if e <= 0:
                        continue
                    try:
                        reacs.append(reaction_at(self.experiment, ch, float(e), self.layers))
                    except ValueError:  # below the excitation threshold, or no energy left
                        continue
                    ws.append(w)
                if reacs:
                    ws = np.array(ws)[:, None]
                    for p in PARTICLES:
                        for sigma, _ in self._branches_many(reacs, theta, p):
                            sums[p] += (ws * sigma).sum(axis=0)
                for p in PARTICLES:
                    per_seg = np.bincount(index, sums[p] * dom, minlength=len(segs))
                    rate = self.particles_per_second * ch.atoms_per_cm2 * per_seg * MB_CM2
                    rows += [RateRow(g.name, seg, ch.label, p, float(r)) for seg, r in zip(segs, rate) if r > 0]
        return rows

    def _branches_many(self, reacs: list, theta: np.ndarray, particle: str) -> list:
        """[(lab cross section mb/sr, lab energy MeV)] for the first and the second kinematic branch at lab angles
        ``theta``, for several reactions at once (one per beam energy; :func:`physim.nuclear.kinematics.lab_points`):
        arrays of shape (len(reacs), len(theta)), zero (energy NaN) where the particle cannot go, leaves with less
        than ``min_energy``, or the ejectile's CM angle is below the floor."""
        out = []
        for pt in lab_points([r.kinematics for r in reacs], theta, particle):
            th_ej = pt.theta_cm if particle == "ejectile" else 180.0 - pt.theta_cm
            ok = (~np.isnan(th_ej)) & (th_ej >= self.theta_floor) & (np.nan_to_num(pt.energy) >= self.min_energy)
            if not ok.any():  # a branch the particle never reaches (the second one, mostly)
                out.append((np.zeros(ok.shape), np.full(ok.shape, np.nan)))
                continue
            safe = np.where(ok, th_ej, 90.0)
            with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
                sigma = np.array([np.asarray(r.cross_section_cm(safe[i])) if ok[i].any() else np.zeros(ok.shape[1])
                                  for i, r in enumerate(reacs)])
                sigma = sigma * pt.jacobian
            out.append((np.where(ok, sigma, 0.0), np.where(ok, pt.energy, np.nan)))
        return out

    # -- totals ---------------------------------------------------------------------------------------------------

    def _counted_fraction(self, detector: str, channel: str, particle: str) -> float:
        peaks = [pk for pk in self._all_peaks(detector, None) if pk.channel == channel and pk.particle == particle]
        total = sum(pk.rate for pk in peaks)
        return sum(pk.rate * pk.counted_fraction for pk in peaks) / total if total > 0 else 0.0

    def _select(self, detector: str, segment: Optional[tuple], counted: bool) -> list:
        """The rows of a detector (or one segment), scaled by the counted fraction with ``counted``. Kept per
        choice: the rows do not change once computed, and every total, table and warning asks for them again."""
        key = (detector, None if segment is None else tuple(segment), bool(counted))
        kept = self.__dict__.setdefault("_selected", {})
        if key not in kept:
            kept[key] = self._select_rows(*key)
        return list(kept[key])

    def _select_rows(self, detector: str, segment: Optional[tuple], counted: bool) -> list:
        if detector not in [g.name for g in self.array]:
            raise KeyError(detector)
        rows = [r for r in self.rows if r.detector == detector and (segment is None or r.segment == segment)]
        if counted:
            frac = {(c, p): self._counted_fraction(detector, c, p) for c, p in {(r.channel, r.particle) for r in rows}}
            rows = [RateRow(r.detector, r.segment, r.channel, r.particle, r.rate * frac[(r.channel, r.particle)])
                    for r in rows if frac[(r.channel, r.particle)] > 0]
        return rows

    @property
    def measured(self) -> str:
        """What a Coulomb-excitation measurement counts, and so what ``counts_wanted`` and the beam time refer to:
        ``"all"`` particles for elastic scattering; for Coulomb excitation, the excitation events seen in a particle
        detector together with a γ ray (``"coincidences"``), or without γ detectors the excitation events alone
        (``"excitations"``)."""
        if not any(c.excitation is not None for c in self.channels):
            return "all"
        return "coincidences" if self.experiment.gamma_detectors else "excitations"

    def gamma_efficiency(self, detector: Optional[str] = None) -> tuple:
        """(full-energy-peak efficiency of all γ detectors together at the energy of the excited state's γ ray,
        whether any of it comes from the typical response model rather than from the setup's ``efficiency`` or
        ``efficiency_curve``). See :mod:`physim.nuclear.response`.

        With ``detector`` (a particle detector), each γ detector's efficiency is multiplied by its factor from
        :func:`physim.nuclear.gamma.correlation_table`: the γ rays seen in coincidence with that detector are
        not emitted evenly in all directions."""
        dets = self.experiment.gamma_detectors
        exc = self.experiment.excitation
        energy = exc.energy_mev if exc is not None else 1.332492
        if getattr(self, "_gamma_effs", None) is None:
            self._gamma_effs = [g.peak_efficiency(energy, self.experiment) for g in dets]
        effs = self._gamma_effs
        if detector is not None and exc is not None:
            factors = self.correlation().get(detector, {})
            effs = [e * factors.get(g.name or f"G{i + 1}", 1.0) for i, (g, e) in enumerate(zip(dets, effs))]
        return min(sum(effs), 1.0), any(g.efficiency is None and not g.efficiency_curve for g in dets)

    def correlation(self) -> dict:
        """{particle detector: {γ detector: factor}}: the γ rays each γ detector sees in coincidence with each
        particle detector, relative to isotropic emission. Computed once."""
        if getattr(self, "_correlation", None) is None:
            from .gamma import correlation_table

            out: dict = {}
            if self.experiment.excitation is not None and self.experiment.gamma_detectors:
                for row in correlation_table(self.experiment):
                    out.setdefault(row["particle_detector"], {})[row["gamma_detector"]] = row["detector_factor"]
            self._correlation = out
        return self._correlation

    def other_peak(self) -> Optional[dict]:
        """The peak of the other nucleus of the collision that the events hold when the setup has its level
        scheme (:func:`physim.nuclear.gamma.other_excitation`): its first transition to the ground state.
        {"label", "nuclide", "energy_kev", "emitted_per_s": {particle detector: those γ rays emitted per second
        with a counted particle in it}, "efficiency": {particle detector: the full-energy efficiency of all γ-ray
        detectors together at that energy, with the angular correlation}, "coincidence_per_s": {particle
        detector: counts per second in the full-energy peak}}. None without that scheme or without γ-ray
        detectors. Computed once: the Rutherford cross section times the transition's yield per collision, over
        each detector's face and the depth of the target."""
        if "_other_peak" in self.__dict__:
            return self._other_peak
        from .gamma import correlation_table, other_excitation, transition_lines

        self._other_peak = None
        exp = self.experiment
        sol = other_excitation(exp) if exp.excitation is not None and exp.gamma_detectors else None
        if sol is None:
            return None
        levels = sol.scheme.levels
        lev_f = sol.ground
        lev_i = min((n for n in range(len(levels)) if n != lev_f), key=lambda n: levels[n].energy.value)
        line = next((ln for ln in transition_lines(exp) if ln["nucleus"] == 1 and (ln["initial"], ln["final"])
                     == (lev_i, lev_f)), None)
        if line is None:
            return None
        grid, states = sol.table()
        grid, per_collision = np.r_[0.0, grid], np.r_[0.0, [st.gamma_yield(lev_i, lev_f) for st in states]]
        ch = next(c for c in self.channels if c.excitation is None and c.layer == 0
                  and c.nuclide.name == sol.target)
        zs, wz = self._depth_nodes(0)
        es = beam_energy_at(exp, 0, zs, self.layers)
        reacs = [_YieldReaction(self.beam, ch.nuclide.name, float(e), grid, per_collision) for e in es if e > 0]
        ws = np.array([w for e, w in zip(es, wz) if e > 0])[:, None]
        factors: dict = {}
        for row in correlation_table(exp, nucleus="other"):
            factors.setdefault(row["particle_detector"], {})[row["gamma_detector"]] = row["detector_factor"]
        e_line = 1e-3 * line["energy_kev"]
        effs = {(gd.name or f"G{i + 1}"): gd.peak_efficiency(e_line, exp) for i, gd in enumerate(exp.gamma_detectors)}
        emitted, efficiency, seen = {}, {}, {}
        for g in self.array:
            parts = [g.directions(seg, self._order_for(g)) for seg in g.segments]
            dirs = np.concatenate([p[0] for p in parts])
            dom = np.concatenate([p[1] for p in parts])
            theta = angles(dirs)[0]
            dom = np.where((theta >= self.theta_floor) & self.array.visible(g.name, dirs), dom, 0.0)
            total = 0.0
            for p in PARTICLES:
                per_dir = np.zeros(len(dirs))
                if reacs:
                    for sigma, _ in self._branches_many(reacs, theta, p):
                        per_dir += (ws * sigma).sum(axis=0)
                rate = self.particles_per_second * ch.atoms_per_cm2 * float(np.sum(per_dir * dom)) * MB_CM2
                total += rate * self._counted_fraction(g.name, ch.label, p)
            emitted[g.name] = total
            efficiency[g.name] = min(sum(e * factors.get(g.name, {}).get(n, 1.0) for n, e in effs.items()), 1.0)
            seen[g.name] = total * efficiency[g.name]
        self._other_peak = {"label": line["label"], "nuclide": line["nuclide"], "energy_kev": line["energy_kev"],
                            "initial": lev_i, "final": lev_f, "emitted_per_s": emitted, "efficiency": efficiency,
                            "coincidence_per_s": seen}
        return self._other_peak

    def rate(self, detector: str, segment: Optional[tuple] = None, counted: bool = True, what: str = "all") -> float:
        """Counts per second in a detector (or one segment). With ``counted``, each channel's rate is scaled by the
        fraction of its particles measured above the threshold (estimated for the whole detector from mean
        energies, see :meth:`peaks`; the event generator counts each particle exactly). Without it, every particle
        reaching the face is counted.

        ``what`` picks the events: ``"all"``, ``"excitations"`` (Coulomb-excitation channels only) or
        ``"coincidences"`` (excitations with the γ ray in a γ detector, see :meth:`gamma_efficiency`); ``"measured"``
        is :attr:`measured`."""
        if what == "measured":
            what = self.measured
        if what not in ("all", "excitations", "coincidences"):
            raise ValueError(f"what must be 'all', 'excitations', 'coincidences' or 'measured', not {what!r}")
        rows = self._select(detector, segment, counted)
        if what != "all":
            excited = {c.label for c in self.channels if c.excitation is not None}
            rows = [r for r in rows if r.channel in excited]
        total = float(sum(r.rate for r in rows))
        return total * self.gamma_efficiency(detector)[0] if what == "coincidences" else total

    def per_detector(self, counted: bool = True) -> dict:
        """{detector: counts per second}."""
        return {g.name: self.rate(g.name, counted=counted) for g in self.array}

    def per_segment(self, detector: str, counted: bool = True) -> dict:
        """{segment: counts per second} for one detector."""
        out = {seg: 0.0 for seg in self.array[detector].segments}
        for r in self._select(detector, None, counted):
            out[r.segment] += r.rate
        return out

    def by_channel(self, detector: str, counted: bool = False) -> dict:
        """{(channel, particle): counts per second} for one detector."""
        out: dict = {}
        for r in self._select(detector, None, counted):
            out[(r.channel, r.particle)] = out.get((r.channel, r.particle), 0.0) + r.rate
        return out

    def counts_in_run(self, detector: str, segment: Optional[tuple] = None, counted: bool = True,
                      what: str = "all") -> float:
        """Counts expected in the planned beam time (``what`` as in :meth:`rate`)."""
        return self.rate(detector, segment, counted, what) * self.beam_time_s

    def beam_time_for(self, detector: str, counts: Optional[int] = None, segment: Optional[tuple] = None,
                      counted: bool = True, what: str = "all") -> float:
        """Seconds of beam needed to collect ``counts`` (default: the setup's ``counts_wanted``) of the events picked
        by ``what`` (as in :meth:`rate`); the relative statistical error is then 1/√counts. Infinite if the detector
        sees none."""
        counts = self.counts_wanted if counts is None else counts
        if counts is None:
            raise ValueError("give counts, or set run.counts_wanted in the setup")
        r = self.rate(detector, segment, counted, what)
        return counts / r if r > 0 else math.inf

    def relative_error(self, detector: str, segment: Optional[tuple] = None, counted: bool = True,
                       what: str = "all") -> float:
        """Relative statistical error 1/√N of the counts collected in the planned beam time."""
        n = self.counts_in_run(detector, segment, counted, what)
        return 1 / math.sqrt(n) if n > 0 else math.inf

    # -- peaks ----------------------------------------------------------------------------------------------------

    def peaks(self, detector: str, segment: Optional[tuple] = None, min_share: float = 1e-3) -> list:
        """Expected peaks in a detector (or one segment): one per channel, particle and kinematic branch carrying
        at least ``min_share`` of the detector's rate (Coulomb-excitation peaks always), strongest first."""
        peaks = self._all_peaks(detector, segment)
        total = sum(pk.rate for pk in peaks)
        excited = {ch.label for ch in self.channels if ch.excitation is not None}
        # Coulomb-excitation peaks are what the experiment is for: always listed, however weak.
        return [pk for pk in peaks if total > 0 and (pk.rate >= min_share * total or (pk.channel in excited
                                                                                        and pk.rate > 0))]

    def _all_peaks(self, detector: str, segment) -> list:
        key = (detector, None if segment is None else tuple(segment))
        if key not in self._peak_cache:
            peaks = []
            for ch in self.channels:
                for p in PARTICLES:
                    peaks += self._peaks_for(detector, key[1], ch, p)
            self._peak_cache[key] = sorted(peaks, key=lambda pk: -pk.rate)
        return self._peak_cache[key]

    def _peaks_for(self, detector: str, segment, ch: Channel, particle: str, n_depth: int = 24,
                   order: int = 16) -> list:
        exp, layers = self.experiment, self.layers
        g = self.array[detector]
        setup = exp.detectors[[x.name for x in self.array].index(detector)]
        if detector in self._partly_hidden and segment is None:
            order = max(order, SHADOW_ORDER)
        dirs, dom = g.directions(segment, order)
        theta = angles(dirs)[0]
        dom = np.where((theta >= self.theta_floor) & self.array.visible(detector, dirs), dom, 0.0)
        lay = layers[ch.layer]
        zs = (np.arange(n_depth) + 0.5) / n_depth * lay.thickness
        es = beam_energy_at(exp, ch.layer, zs, layers)
        ion = self.beam if particle == "ejectile" else ch.nuclide.name
        det_mat = setup.material_data()
        det_st = stopping(ion, det_mat)
        dead = det_mat.areal_density_mg_cm2(_q(setup.dead_layer)) if setup.dead_layer is not None else 0.0
        active = det_mat.areal_density_mg_cm2(_q(setup.thickness))
        fwhm = _q(setup.resolution).to("MeV") if setup.resolution is not None else 0.0
        threshold = _q(setup.threshold).to("MeV") if setup.threshold is not None else 0.0
        cos_inc = np.clip(-(dirs @ g.n), 1e-9, 1.0)

        def response(e_face, ci):
            e_a = _after(det_st, e_face, dead / ci)
            return e_a - _after(det_st, e_a, active / ci)

        reacs, js = [], []
        for j, e in enumerate(es):
            if e <= 0:
                continue
            try:
                reacs.append(reaction_at(self.experiment, ch, float(e), self.layers))
            except ValueError:
                continue
            js.append(j)
        branches = self._branches_many(reacs, theta, particle) if reacs else []
        out = []
        for branch in (0, 1):
            energy = np.full((n_depth, len(dirs)), np.nan)
            sigma = np.zeros((n_depth, len(dirs)))
            if branch < len(branches):
                sigma[js], energy[js] = branches[branch]
            weight = np.where(np.isnan(energy), 0.0, sigma * dom[None, :])
            if weight.sum() <= 0:
                continue
            rate = self.particles_per_second * ch.atoms_per_cm2 * weight.sum() / n_depth * MB_CM2
            e_lab = np.nan_to_num(energy)
            z2 = np.broadcast_to(zs[:, None], e_lab.shape)
            d2 = np.broadcast_to(dirs[None, :, :], e_lab.shape + (3,))
            e_face = exit_energy(exp, layers, ch.layer, z2.ravel(), ion, e_lab.ravel(), d2.reshape(-1, 3))
            e_face = e_face.reshape(e_lab.shape)
            measured = response(e_face, cos_inc[None, :])
            wsum = weight.sum()
            mean = float((weight * measured).sum() / wsum)
            geometry = math.sqrt(max(float((weight * (measured - mean) ** 2).sum() / wsum), 0.0))
            per_node = (weight * measured).sum(axis=0) / np.maximum(weight.sum(axis=0), 1e-300)
            node_w = weight.sum(axis=0)
            kinematic = math.sqrt(max(float((node_w * (per_node - mean) ** 2).sum() / node_w.sum()), 0.0))
            thick_var = (weight * (measured - per_node[None, :]) ** 2).sum() / wsum
            spread, straggling = self._spread_and_straggling(ch, particle, branch, dirs, node_w, cos_inc, det_st,
                                                             dead, active)
            res = fwhm / FWHM_PER_SIGMA
            total = math.sqrt(geometry**2 + spread**2 + straggling**2 + res**2)
            out.append(Peak(detector, segment, ch.label, particle, branch + 1, rate, mean, total, {
                "target thickness": math.sqrt(max(float(thick_var), 0.0)), "kinematic": kinematic,
                "geometry": geometry, "beam energy spread": spread, "straggling": straggling, "resolution": res,
            }, float((weight * e_face).sum() / wsum), threshold,
                float((weight * ((measured >= threshold) & (measured > 0))).sum() / wsum)))
        return out

    def _spread_and_straggling(self, ch, particle, branch, dirs, node_w, cos_inc, det_st, dead, active) -> tuple:
        """Widths (σ, MeV) from the beam energy spread and from straggling, carried to the measured energy along
        the detector's mean direction from the middle of the layer."""
        exp, layers = self.experiment, self.layers
        d = (node_w[:, None] * dirs).sum(axis=0)
        d /= np.linalg.norm(d)
        ci = float((node_w * cos_inc).sum() / node_w.sum())
        theta = float(angles(d)[0])
        c = math.cos(math.radians(tilt_deg(exp)))
        z_mid = layers[ch.layer].thickness / 2
        ion = self.beam if particle == "ejectile" else ch.nuclide.name
        results = []
        for var0, straggle in ((energy_sigma(exp) ** 2, False), (0.0, True)):
            e, var = exp.beam.energy_mev, var0
            for j, lay in enumerate(layers[: ch.layer + 1]):
                dz = z_mid if j == ch.layer else lay.thickness
                e, var = _carry(stopping(self.beam, lay.material), e, var, dz / c, straggle)
            if e <= 0:
                results.append(0.0)
                continue

            def lab(energy):
                try:
                    tb = TwoBody(self.beam, ch.nuclide.name, energy, **ch.kinematics_args)
                except ValueError:  # below the excitation threshold at the middle of the layer
                    return math.nan
                pt = tb.at_lab(theta, particle)[branch]
                return float(pt.energy) if pt is not None else math.nan

            h = 1e-4 * e
            e_out = lab(e)
            if not (np.isfinite(e_out) and np.isfinite(lab(e - h))):
                results.append(0.0)
                continue
            var *= ((lab(e + h) - lab(e - h)) / (2 * h)) ** 2
            fwd, steps = _exit_paths(exp, layers, ch.layer, z_mid, d[None, :])
            for j, path, forward in steps:
                if forward == bool(fwd[0]):
                    e_out, var = _carry(stopping(ion, layers[j].material), e_out, var, float(np.ravel(path)[0]),
                                        straggle)
            e_a, var = _carry(det_st, e_out, var, dead / ci, straggle)
            e_b = float(_after(det_st, e_a, active / ci))
            if e_b > 0:  # punches through: the deposit is E_a − E_b
                _, var_b = _carry(det_st, e_a, 0.0, active / ci, straggle)
                factor = 1 - det_st.stopping_power(e_b) / det_st.stopping_power(e_a)
                var = factor**2 * var + var_b
            results.append(math.sqrt(max(var, 0.0)))
        return results[0], results[1]

    # -- warnings -------------------------------------------------------------------------------------------------

    def _all_orders(self) -> bool:
        """Whether the excitation follows a level scheme to all orders (then the events and these rates do)."""
        from .gamma import cascade_excitation, excitation_order

        try:
            return cascade_excitation(self.experiment) is not None and excitation_order(self.experiment) == "all"
        except Exception:  # noqa: BLE001 - a warning's wording must not stop the warnings
            return False

    def warnings(self, max_rate: float = MAX_RATE) -> list:
        """Plain-language warnings: detectors counting faster than ``max_rate`` per second, too few counts in the
        planned run, Rutherford's formula not applying, the beam stopping in the target, and layout problems."""
        out = list(self.array.warnings())
        last = len(self.layers) - 1
        if float(beam_energy_at(self.experiment, last, [self.layers[last].thickness], self.layers)[0]) <= 0:
            out.append("The beam stops inside the target: reactions happen only down to the depth it reaches.")
        seen = set()
        for ch in self.channels:
            if ch.excitation is not None and not layer_energy(self.experiment, ch.layer, self.layers) > 0:
                continue  # the beam stops before this layer
            acc = cm_acceptance(self.experiment, self.layers, ch, self.array, self.theta_floor)
            if acc is None:
                continue
            if ch.excitation is not None:
                found = coulex_for(self.experiment, ch, self.layers).warnings(
                    theta_cm_max=acc[1], all_orders=self._all_orders())
            else:
                found = Rutherford(self.beam, ch.nuclide.name, self.experiment.beam.energy_mev).warnings(
                    theta_cm_max=acc[1])
            for w in found:
                text = f"{ch.label}: {w}"
                if text not in seen:
                    seen.add(text)
                    out.append(text)
        hours = self.beam_time_s / 3600
        shadows = self.array.shadowing()
        for name, r in self.per_detector().items():
            if r > max_rate:
                out.append(f"{name} counts {r:.3g} per second, above {max_rate:.0f}/s: expect pile-up and dead time "
                           f"(and radiation damage at forward angles). Lower the beam current or move it back.")
            hidden_by = [front for (front, behind), frac in shadows.items() if behind == name and frac >= 1 - 1e-6]
            if r == 0 and hidden_by:
                out.append(f"{name} records nothing: {hidden_by[0]} stops every particle before it gets there. Move "
                           f"one of them, or remove {name}.")
            elif r == 0:
                out.append(f"{name} records nothing: no scattered particle reaches it above its threshold.")
            elif self.counts_wanted is not None and r * self.beam_time_s < self.counts_wanted:
                need = self.counts_wanted / r / 3600
                out.append(f"{name} collects {r * self.beam_time_s:.0f} counts in {hours:g} h; "
                           f"{self.counts_wanted} counts need {need:.3g} h.")
        return out
