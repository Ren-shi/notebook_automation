"""γ rays from a nucleus excited in the reaction: Doppler shift and Doppler broadening.

The excited nucleus decays in flight. A γ ray emitted at angle α to the nucleus's velocity β has the lab energy

    E_γ = E₀ √(1 − β²) / (1 − β cos α).

For each pair of a particle detector (which fixes the direction of the scattered particle, and so of the excited
nucleus) and a γ detector, physim samples the particle detector's face, the γ detector's face and the depth of the
reaction in the target, and reports the mean γ energy and its spread::

    from physim.nuclear import Experiment
    from physim.nuclear.gamma import doppler_table

    for row in doppler_table(Experiment.example("coulex_ni58")):
        print(row["particle_detector"], row["gamma_detector"], row["mean_kev"], row["fwhm_kev"])

Assumptions: the state decays after the nucleus leaves the target (lifetimes of picoseconds and more, thin
targets), isotropically; the nucleus moves with its velocity just after the reaction, slowed to the target exit.
"""

from __future__ import annotations

import json
import math
import threading
from collections import OrderedDict
from functools import lru_cache
from typing import Optional

import numpy as np

from .detectors import Array, _gauss, angles
from .rates import _q, beam_energy_at, beam_ion, exit_energy, stack, tilt_deg

FWHM_PER_SIGMA = 2.3548200450309493


def doppler_energy(e0: float, beta, cos_alpha):
    """E₀ √(1 − β²) / (1 − β cos α), elementwise."""
    beta = np.asarray(beta, dtype=float)
    return e0 * np.sqrt(1 - beta**2) / (1 - beta * np.asarray(cos_alpha, dtype=float))


def _disc_directions(direction, half_angle_deg: float, n: int = 12) -> np.ndarray:
    """Unit vectors spread uniformly over a cone (a disc detector seen from the target): n² points."""
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    e1 = np.cross(d, [0.0, 1.0, 0.0] if abs(d[1]) < 0.9 else [1.0, 0.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    t = math.tan(math.radians(half_angle_deg))
    # Uniform in area on the disc: radius ∝ √u.
    u = (np.arange(n) + 0.5) / n
    a = 2 * math.pi * (np.arange(n) + 0.5) / n
    r = t * np.sqrt(u)
    pts = d[None, :] + (r[:, None, None] * (np.cos(a)[None, :, None] * e1 + np.sin(a)[None, :, None] * e2)).reshape(
        -1, 3)
    return pts / np.linalg.norm(pts, axis=1, keepdims=True)


def doppler_table(experiment, depth_points: int = 8, order: int = 10) -> list:
    """For every particle detector × γ detector: the Doppler-shifted γ energy (mean and spread) of the excited
    nucleus, for the scattered particle reaching that particle detector (first kinematic solution).

    Each sample (a point of the particle detector, a direction into the γ detector, a depth in the target) is
    weighted by the Coulomb-excitation cross section there; the kinematics and the excitation probability are those
    at the beam energy of each depth (:func:`~physim.nuclear.rates.depth_probability`)."""
    from . import data
    from .kinematics import TwoBody
    from .rates import channels, depth_probability, rutherford_factor
    from .rutherford import Rutherford

    exc = experiment.excitation
    if exc is None:
        raise ValueError("the setup has no excited state ([reaction] type = \"coulex\")")
    layers = stack(experiment)
    lay = layers[0]
    target = data.nuclide(max(lay.nuclides.items(), key=lambda kv: kv[1])[0]).name
    beam = beam_ion(experiment)
    zs = (np.arange(depth_points) + 0.5) / depth_points * lay.thickness
    es = beam_energy_at(experiment, 0, zs, layers)
    e0 = exc.energy_mev
    excite = "recoil" if exc.excite == "target" else "ejectile"
    emitter = target if exc.excite == "target" else beam
    m_emitter = data.nuclide(emitter).nuclear_mass_mev + e0
    channel = next(ch for ch in channels(layers, experiment)
                   if ch.excitation is not None and ch.layer == 0 and ch.nuclide.name == target)
    depth = depth_probability(experiment, channel, layers)
    # One entry per crystal: a clover's four crystals are corrected for the Doppler shift separately.
    gammas = []
    for i, gd in enumerate(experiment.gamma_detectors):
        res = _q(gd.resolution).to("MeV") / FWHM_PER_SIGMA if gd.resolution is not None else 0.0
        for label, centre, radius in gd.elements():
            d = math.sqrt(sum(c * c for c in centre))
            gammas.append(((gd.name or f"G{i + 1}") + (f" {label}" if label else ""),
                           _disc_directions(tuple(c / d for c in centre), math.degrees(math.atan2(radius, d))),
                           res))
    rows = []
    for g in Array.from_experiment(experiment):
        dirs, dom = g.directions(None, order)
        theta = angles(dirs)[0]
        phi = np.radians(angles(dirs)[1])
        vel, wts, betas = [], [], []
        for z, e in zip(zs, es):
            if e <= e0:
                continue
            try:
                tb = TwoBody(beam, target, float(e), excitation_mev=e0, excite=excite)
            except ValueError:  # below the excitation threshold at this depth
                continue
            first, _ = tb.at_lab(theta, "ejectile")
            th_cm = np.asarray(first.theta_cm)
            ok = ~np.isnan(th_cm)
            if not ok.any():
                continue
            safe = np.where(ok, th_cm, 90.0)
            # P times Alder and Winther's inelastic Rutherford factor: the elastic one times v_i/v_f.
            ruth = Rutherford(beam, target, float(e))
            sigma = np.where(ok, np.asarray(depth.at(float(e)).probability(safe))
                             * rutherford_factor(ruth.e_cm, e0) * np.asarray(ruth.cross_section_cm(safe))
                             * np.nan_to_num(np.asarray(first.jacobian)), 0.0)
            if excite == "recoil":
                rec = tb.recoil_for(np.where(ok, th_cm, 90.0))
                th_e, ph_e, t_e = np.radians(rec.theta_lab), phi + math.pi, np.asarray(rec.energy)
            else:
                th_e, ph_e, t_e = np.radians(first.theta_lab), phi, np.asarray(first.energy)
            v = np.stack([np.sin(th_e) * np.cos(ph_e), np.sin(th_e) * np.sin(ph_e), np.cos(th_e)], axis=-1)
            t_e = exit_energy(experiment, layers, 0, np.full(theta.shape, z), emitter, np.nan_to_num(t_e), v)
            gam = 1 + t_e / m_emitter
            beta = np.sqrt(np.maximum(1 - 1 / gam**2, 0.0))
            keep = ok & (t_e > 0) & (sigma > 0)
            vel.append(v[keep])
            wts.append((sigma * dom)[keep])
            betas.append(beta[keep])
        if not vel:
            continue
        v, w, beta = np.concatenate(vel), np.concatenate(wts), np.concatenate(betas)
        if w.sum() <= 0:
            continue
        for name, gdirs, res in gammas:
            eg = doppler_energy(e0, beta[:, None], v @ gdirs.T)
            ww = np.repeat(w[:, None], gdirs.shape[0], axis=1)
            mean = float(np.average(eg, weights=ww))
            sd = math.sqrt(float(np.average((eg - mean) ** 2, weights=ww)))
            rows.append({
                "particle_detector": g.name, "gamma_detector": name, "emitter": emitter, "e0_kev": e0 * 1e3,
                "mean_kev": mean * 1e3, "shift_kev": (mean - e0) * 1e3,
                "doppler_fwhm_kev": FWHM_PER_SIGMA * sd * 1e3,
                "fwhm_kev": FWHM_PER_SIGMA * math.sqrt(sd**2 + res**2) * 1e3,
                "beta_mean": float(np.average(beta, weights=w)),
            })
    return rows


def excitation_of(experiment, isotropic=None, energy=None):
    """The :class:`~physim.nuclear.orientation.Excitation` of the setup's excited state (the ``[reaction]``
    section: one state reached from a 0⁺ ground state), at beam energy ``energy`` (MeV), by default the energy in
    the middle of the target. ``isotropic`` overrides the setup's ``emission`` switch."""
    from . import data

    exc = experiment.excitation
    if exc is None:
        raise ValueError("the setup has no excited state ([reaction] type = \"coulex\")")
    layers = stack(experiment)
    lay = layers[0]
    target = data.nuclide(max(lay.nuclides.items(), key=lambda kv: kv[1])[0]).name
    beam = beam_ion(experiment)
    # The same energy, rounded the same way, as the rates use: the orbit integrals are then shared.
    if energy is None:
        energy = beam_energy_at(experiment, 0, [lay.thickness / 2], layers)[0]
    e_mid = round(float(energy), 9)
    flat = (exc.emission == "isotropic") if isotropic is None else isotropic
    return _excitation(beam, target, e_mid, exc.excite, exc.energy_mev, exc.multipolarity, exc.b_up_e2fm, flat,
                       e1_polarisation(experiment))


#: Accuracy of the coupled equations' probabilities in the events (relative, roughly; see
#: :class:`~physim.nuclear.coupled.CoupledChannels`): far below what any run's statistics show.
EVENT_TOLERANCE = 1e-4
#: The solutions kept (one per beam, partner nucleus, energy, level scheme, order and E1 polarisation), and how many.
_SOLVED: OrderedDict = OrderedDict()
_SOLVED_MOST = 48
_SOLVING = threading.Lock()


def excitation_order(experiment) -> str:
    """ "all" or "first": how the events populate the levels of a level scheme (``[reaction] order``)."""
    exc = experiment.excitation
    return (getattr(exc, "order", None) or "all") if exc is not None else "all"


def e1_polarisation(experiment) -> float:
    """GOSIA's parameter d of the E1 polarisation correction in the setup's collisions: its ``[reaction]
    e1_polarisation``, or GOSIA's default, 0.0038 (:data:`~physim.nuclear.coulex.GOSIA_DIPOLE_POLARISATION`), when the
    setup does not set it (a setup that is not a Coulomb excitation cannot); 0 when it leaves the correction out."""
    from .coulex import GOSIA_DIPOLE_POLARISATION

    exc = experiment.excitation
    d = getattr(exc, "e1_polarisation", None) if exc is not None else None
    return GOSIA_DIPOLE_POLARISATION if d is None else float(d)


def deorientation_of(experiment) -> str:
    """ "vacuum" or "none": whether the nuclei of the setup that decay in vacuum are deoriented there
    (:mod:`physim.nuclear.deorientation`): its ``[reaction] deorientation``, "vacuum" when the setup does not say. A
    lifetime measurement (a ``[lifetime]`` section) is "none": a plunger has its own ``lifetime.deorientation``, or none."""
    if getattr(experiment, "lifetime", None) is not None:
        return "none"
    exc = experiment.excitation
    return (getattr(exc, "deorientation", None) or "vacuum") if exc is not None else "vacuum"


def deorientation_parameters(experiment, nucleus: str = "reaction"):
    """The :class:`~physim.nuclear.deorientation.Parameters` (GOSIA's, with the setup's g factor) of the deorientation in
    vacuum of the ``[reaction]`` nucleus, or of the ``"other"`` one of the collision, which has g = Z/A; None when there
    is none (:func:`deorientation_of`)."""
    from .deorientation import Parameters

    if deorientation_of(experiment) == "none":
        return None
    g = getattr(experiment.excitation, "deorientation_g", None) if nucleus == "reaction" else None
    return Parameters(g_factor=g)


def other_role(experiment) -> str:
    """ "projectile" or "target": the nucleus of the collision that the ``[reaction]`` does not name."""
    return "projectile" if experiment.excitation.excite == "target" else "target"


def other_excitation(experiment, isotropic=None, order=None):
    """The excitation of the other nucleus of the collision (the target when the ``[reaction]`` excites the beam,
    the beam when it excites the target), from its level scheme in the setup's ``[levels]``: the solution of its
    coupled equations (or first order, as ``[reaction] order`` says) for scattering on the target's main nuclide,
    at the beam energy in the middle of the target. None when the setup has no level scheme for it, or one
    without a matrix element from the ground state. With it, the simulated events hold that nucleus's γ rays too
    (:func:`physim.nuclear.gamma_events.simulate_gammas`)."""
    if experiment.excitation is None:
        return None
    return cascade_excitation(experiment, isotropic, order, nucleus="other")


def emitters(experiment) -> list:
    """The nuclei whose γ rays the events hold, the ``[reaction]`` nucleus first: dicts with "nucleus" (0 or 1,
    as the γ-ray events' ``nucleus`` column), "nuclide", "role" ("projectile" or "target"), "correction" (the
    Doppler correction that is right for it, "projectile" or "recoil"), "scheme" (the level scheme its γ rays
    follow, or None for the one state of the ``[reaction]``) and "solution" (the excitation of that scheme, or
    None)."""
    from . import data

    exc = experiment.excitation
    if exc is None:
        return []
    lay = stack(experiment)[0]
    names = {"projectile": beam_ion(experiment),
             "target": data.nuclide(max(lay.nuclides.items(), key=lambda kv: kv[1])[0]).name}
    out = []
    for n, role in enumerate((exc.excite, other_role(experiment))):
        solution = cascade_excitation(experiment) if n == 0 else other_excitation(experiment)
        if n == 1 and solution is None:
            break
        out.append({"nucleus": n, "nuclide": names[role], "role": role,
                    "correction": "recoil" if role == "target" else "projectile",
                    "scheme": solution.scheme if solution is not None else None, "solution": solution})
    return out


def transition_lines(experiment) -> list:
    """The γ-ray lines the events can hold, the ``[reaction]`` nucleus first: dicts with "nucleus", "nuclide",
    "initial", "final" (level indices in that nucleus's scheme), "energy_kev", "label" (such as
    "62Ni 2+ → 0+") and "correction" (the Doppler correction that sharpens it). For the ``[reaction]`` state
    without a scheme, its one line."""
    exc = experiment.excitation
    out = []
    for em in emitters(experiment):
        if em["solution"] is None:
            out.append({"nucleus": em["nucleus"], "nuclide": em["nuclide"], "initial": 1, "final": 0,
                        "energy_kev": 1e3 * exc.energy_mev,
                        "label": f"{em['nuclide']} {1e3 * exc.energy_mev:.1f} keV", "correction": em["correction"]})
            continue
        levels = em["solution"].scheme.levels
        names = level_names(levels)
        for i, f in em["solution"].transitions():
            if em["solution"].decay(i, f).gamma <= 0:
                continue
            out.append({"nucleus": em["nucleus"], "nuclide": em["nuclide"], "initial": int(i), "final": int(f),
                        "energy_kev": float(levels[i].energy.value - levels[f].energy.value),
                        "label": f"{em['nuclide']} {names[i]} → {names[f]}", "correction": em["correction"]})
    return out


def level_names(levels: list) -> list:
    """A name for each level of a scheme: its spin and parity, with its place among the levels of that spin and
    parity as a subscript when there are several ("2+₁", "2+₂"; the only 4⁺ level is "4+")."""
    def jpi(level) -> str:
        text = getattr(level, "jpi", None)
        if text:
            return str(text)
        spin = getattr(level, "spin", None)
        return f"{spin:g}" if spin is not None else "?"

    texts = [jpi(lv) for lv in levels]
    order = sorted(range(len(levels)), key=lambda n: levels[n].energy.value)
    seen: dict = {}
    out = list(texts)
    for n in order:
        seen[texts[n]] = seen.get(texts[n], 0) + 1
        if texts.count(texts[n]) > 1:
            out[n] = texts[n] + str(seen[texts[n]]).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
    return out


def has_static_moment(scheme) -> bool:
    """Whether a level scheme sets a diagonal matrix element (a static moment, such as the E2 one that gives the
    spectroscopic quadrupole moment Q_s): then the excitation needs the coupled equations even for one
    transition (the reorientation effect)."""
    return any(me.a == me.b for me in getattr(scheme, "matrix_elements", []) or [])


def cascade_excitation(experiment, isotropic=None, order=None, channel=None, nucleus: str = "reaction",
                       energy=None):
    """The excitation of the excited nucleus's whole level scheme (the setup's ``[levels]`` of the target or the
    beam, whichever the reaction excites), when its decays follow more than one transition: then each excitation
    starts in a level drawn from this solution's populations and decays by a cascade (backlog item 71). None
    when the setup has no such scheme, so a single state keeps the one-γ-ray chain of item 55; but a single
    transition whose scheme sets a diagonal matrix element (:func:`has_static_moment`, such as Q_s of the 2⁺
    state) is solved to all orders too, unless the order is "first": the reorientation effect of Q_s needs the
    coupled equations.

    It is the coupled equations' solution to all orders (:class:`~physim.nuclear.coupled.CoupledChannels`: the
    levels reached in two steps, the reorientation, and the orientation of every level), or with ``order``
    "first" (the setup's ``[reaction] order`` if None) first-order theory
    (:class:`~physim.nuclear.orientation.Excitation`: each level in one step from the ground state). Either has the
    E1 polarisation correction of the setup (:func:`e1_polarisation`).

    ``channel`` (a :class:`~physim.nuclear.rates.Channel` with an excitation) picks the nucleus the excited one
    scatters on and the layer; by default the target's main nuclide. The solution is at beam energy ``energy``
    (MeV), by default the energy in the middle of the layer, and is kept with its table of angles: the events,
    the rates and the γ rays of one setup share it. Through a thick layer the rates and the events take it at
    several energies (:func:`physim.nuclear.rates.depth_probability`).

    ``nucleus`` "other" gives the same for the other nucleus of the collision (:func:`other_excitation`), whose
    scheme may hold a single transition."""
    from . import data
    from .coupled import CoupledChannels
    from .orientation import Excitation

    exc = experiment.excitation
    if exc is None:
        return None
    if nucleus not in ("reaction", "other"):
        raise ValueError("nucleus must be 'reaction' or 'other'")
    excite = exc.excite if nucleus == "reaction" else other_role(experiment)
    role = "target" if excite == "target" else "beam"
    scheme = experiment.levels.get(role)
    if scheme is None:
        return None
    layers = stack(experiment)
    layer = 0 if channel is None else channel.layer
    lay = layers[layer]
    target = (data.nuclide(max(lay.nuclides.items(), key=lambda kv: kv[1])[0]).name if channel is None
              else channel.nuclide.name)
    if energy is None:
        energy = beam_energy_at(experiment, layer, [lay.thickness / 2], layers)[0]
    e_mid = round(float(energy), 9)
    flat = (exc.emission == "isotropic") if isotropic is None else isotropic
    order = order or excitation_order(experiment)
    if order not in ("all", "first"):
        raise ValueError("order must be 'all' or 'first'")
    beam = beam_ion(experiment)
    polarisation = e1_polarisation(experiment)
    key = (beam, target, e_mid, excite, bool(flat), order, polarisation,
           json.dumps(scheme.to_dict(), sort_keys=True, default=str))
    with _SOLVING:  # one thread solves; the others then find the solution
        if key in _SOLVED:
            _SOLVED.move_to_end(key)
            return _SOLVED[key]
        try:
            if order == "first":
                ex = Excitation(beam, target, e_mid, scheme, excite=excite, isotropic=flat,
                                e1_polarisation=polarisation)
            else:
                ex = CoupledChannels(beam, target, e_mid, scheme, excite=excite, isotropic=flat,
                                     tolerance=EVENT_TOLERANCE, e1_polarisation=polarisation)
        except ValueError:
            ex = None
        # The [reaction] nucleus keeps its one-γ-ray chain for a single transition, unless a static moment asks
        # for the coupled equations (to all orders); the other nucleus has no such chain, so one transition is
        # enough for it.
        least = 1 if nucleus == "other" or (order == "all" and has_static_moment(scheme)) else 2
        if ex is not None and (not ex.paths or len(ex.transitions()) < least):
            ex = None
        if ex is not None:
            ex.table()
        _SOLVED[key] = ex
        while len(_SOLVED) > _SOLVED_MOST:
            _SOLVED.popitem(last=False)
        return ex


def state_transition(experiment) -> tuple:
    """(initial, final) level indices of the γ ray of the ``[reaction]`` state to the ground state, as the γ-ray
    events number them: in the level scheme the events follow (the level nearest the state's energy), or (1, 0)
    for a single state."""
    ex = cascade_excitation(experiment)
    if ex is None:
        return 1, 0
    e_kev = 1e3 * experiment.excitation.energy_mev
    levels = ex.scheme.levels
    excited = [n for n in range(len(levels)) if n != ex.ground]
    return min(excited, key=lambda n: abs(levels[n].energy.value - levels[ex.ground].energy.value - e_kev)), ex.ground


def start_levels(ex) -> list:
    """The levels an excitation can start in, by the table of ``ex`` (a solution of
    :func:`cascade_excitation`): those populated after the collision at some angle. In first order they are the
    levels a matrix element joins to the ground state; to all orders every level a chain of matrix elements
    reaches."""
    _, states = ex.table()
    return sorted({n for s in states for n, v in s.direct.items() if v > 0})


def unreached_levels(experiment) -> list:
    """The levels of the excited nucleus's scheme that the simulated events never populate, as indices into the
    scheme: no chain of matrix elements joins them to the ground state (or, in first order, no matrix element
    does), and no populated level decays to them. Empty without a scheme that the events follow."""
    ex = cascade_excitation(experiment)
    if ex is None:
        return []
    _, states = ex.table()
    seen = {n for s in states for n in s.tensors if s.population(n) > 0}
    return [n for n in range(len(ex.scheme.levels)) if n != ex.ground and n not in seen]


@lru_cache(maxsize=64)
def _excitation(beam: str, target: str, e_mid: float, excite: str, e_star: float, multipolarity: str,
                b_up: float, isotropic: bool, polarisation: float):
    from .orientation import Excitation, simple_scheme

    scheme = simple_scheme(target if excite == "target" else beam, e_star, multipolarity, b_up)
    return Excitation(beam, target, e_mid, scheme, excite=excite, isotropic=isotropic, e1_polarisation=polarisation)


class _LevelPopulation:
    """The population of one level (with the feeding from above) after a collision, against the CM angle, from a
    solution's table; shaped like :class:`~physim.nuclear.coulex.Coulex` (``probability``)."""

    def __init__(self, solution, level: int):
        self._grid, states = solution.table()
        self._values = np.array([s.population(level) for s in states])

    def probability(self, theta_cm):
        return np.interp(np.asarray(theta_cm, dtype=float), self._grid, self._values)


#: A density matrix interpolated between beam energies may have eigenvalues down to this below zero (what the
#: polynomial between the solutions is off by); one with a lower eigenvalue is made physical again.
NEGATIVE_EIGENVALUE = 1e-4
#: Density matrices interpolated at once.
DENSITY_BLOCK = 20_000


class DepthOrientation:
    """The population and orientation of the levels of the ``[reaction]`` nucleus at the beam energies the beam has
    through the target: those of ``solutions``, solved at the energies of ``depth`` (the
    :class:`~physim.nuclear.rates.DepthProbability` of the excitation), and interpolated between them the same way.
    A thin target has one solution, at its middle. Get one with :func:`depth_orientation`."""

    def __init__(self, depth, solutions: list):
        self.depth, self.solutions = depth, solutions
        self._tables: dict = {}
        self._densities: dict = {}

    @property
    def thin(self) -> bool:
        """Whether there is one solution, at the middle of the target."""
        return len(self.solutions) == 1

    def population(self, level: int, energy: float):
        """The population of ``level`` at beam energy ``energy`` (MeV): an object with ``probability(theta_cm)``,
        interpolated in ln P. None near or below the excitation threshold."""
        if not self.thin and self.depth.closed(energy):
            return None
        nodes = [_LevelPopulation(s, level) for s in self.solutions]
        return nodes[0] if self.thin else self.depth.between(nodes, energy)

    def coefficients(self, initial: int, final: int, theta_cm, energy, fade: Optional[dict] = None) -> dict:
        """{(k, q): a_kq} of a transition (:meth:`~physim.nuclear.orientation.Populated.coefficients`) at each CM
        angle of ``theta_cm`` (degrees, an array) and the beam energy there (MeV, an array like it or a number).

        With ``fade`` (:meth:`Faded.at`: {k: {level: an array like the angles}}), deoriented in vacuum: the tensors that
        the level ``initial`` decays with are the sums of what the collision leaves each level with, times those."""
        theta = np.asarray(theta_cm, dtype=float)
        out: dict = {}
        for (grid, tensors, populated), w, s in zip(self._table(initial, final), self._weights(theta, energy),
                                                    self.solutions):
            p = np.interp(theta, grid, populated)
            safe = np.where(p > 0, p, 1.0)
            for key, v in tensors.items():
                if key == (0, 0):
                    a = np.ones(len(theta), dtype=complex)
                elif fade is not None and key[0] in fade:
                    a = np.zeros(len(theta), dtype=complex)
                    for n, share in fade[key[0]].items():
                        t = s.tensor_table(n).get(key)
                        if t is not None:
                            a += share * (np.interp(theta, grid, t.real) + 1j * np.interp(theta, grid, t.imag))
                    a *= s.distribution_coefficient(initial, final, key[0]) / safe
                else:
                    a = (np.interp(theta, grid, v.real) + 1j * np.interp(theta, grid, v.imag)) / safe
                out[key] = out.get(key, 0.0) + w * a
        return out

    def density(self, level: int, theta_cm, energy) -> np.ndarray:
        """The density matrices (n × (2I + 1)², trace 1, hermitian, M ascending from −I, in the orbit's frame) that the
        collision leaves ``level`` in, before any decay, at each CM angle of ``theta_cm`` (degrees, an array) and the
        beam energy there (MeV, an array like it or a number).

        Each solution's matrix is normalised at the angle, as
        :meth:`~physim.nuclear.orientation.Excitation.density_table` interpolated in it (unoriented where the level is
        not excited), and the matrices are weighted as :meth:`coefficients` weighs the coefficients: a transition's
        distribution from the one is that from the other. Some weights are negative, so an eigenvalue can be too. One
        below −:data:`NEGATIVE_EIGENVALUE` (close to 0°, or so deep that the level is hardly excited) has the negative
        eigenvalues of its matrix set to zero and the trace to 1 again; smaller ones are left as they are."""
        grid, rho, trace = self._density_table(level)
        theta = np.asarray(theta_cm, dtype=float)
        d = int(round(math.sqrt(rho.shape[-1] / 2)))
        out = np.empty((len(theta), d, d), dtype=complex)
        energy = np.broadcast_to(np.asarray(energy, dtype=float), theta.shape)
        unit = np.eye(d).ravel()
        for block in (slice(a, a + DENSITY_BLOCK) for a in range(0, len(theta), DENSITY_BLOCK)):
            t = np.clip(theta[block], grid[0], grid[-1])
            i = np.clip(np.searchsorted(grid, t) - 1, 0, len(grid) - 2)
            f = (t - grid[i]) / (grid[i + 1] - grid[i])
            weights = np.broadcast_to(self._weights(t, energy[block]), (len(self.solutions), len(t)))
            total = (1 - f) * trace[:, i] + f * trace[:, i + 1]
            share = np.where(total > 0, weights / np.where(total > 0, total, 1.0), 0.0)
            low, high = share * (1 - f), share * f
            flat = np.empty((len(t), 2 * d * d))
            # The events between the same two angles of the grid take the same matrices.
            order = np.argsort(i, kind="stable")
            cuts = np.r_[0, np.flatnonzero(np.diff(i[order])) + 1, len(t)]
            for a, b in zip(cuts[:-1], cuts[1:]):
                here, j = order[a:b], i[order[a]]
                flat[here] = np.concatenate([low[:, here], high[:, here]]).T @ np.concatenate([rho[:, j], rho[:, j + 1]])
            flat[:, ::2] += (weights * (total <= 0)).sum(axis=0)[:, None] / d * unit
            out[block] = _physical(flat.view(complex).reshape(len(t), d, d))
        return out

    def _weights(self, theta: np.ndarray, energy) -> np.ndarray:
        """The weight of each solution (rows) at the beam energies ``energy`` of the CM angles ``theta``."""
        return (np.ones((1, 1)) if self.thin
                else self.depth.lagrange(np.broadcast_to(np.asarray(energy, dtype=float), theta.shape)))

    def _density_table(self, level: int) -> tuple:
        """(grid of CM angles, density matrices of the level as the collision leaves it: solutions × grid ×
        2 (2I + 1)² real numbers, real and imaginary part in turn, their traces: solutions × grid)."""
        if level not in self._densities:
            each = [s.density_table(level) for s in self.solutions]
            rho = np.array([r for _, r in each])
            trace = np.einsum("sgaa->sg", rho).real
            flat = np.ascontiguousarray(rho.reshape(*rho.shape[:2], -1)).view(float)
            self._densities[level] = (each[0][0], flat, trace)
        return self._densities[level]

    def largest(self, initial: int, final: int) -> float:
        """The largest ``Σ |a_kq|`` of a transition's distribution, which 4π W nowhere exceeds: at the CM angles of the
        solutions' grid, and at enough energies to see between them."""
        grid = self._table(initial, final)[0][0]
        energies = [0.0] if self.thin else np.linspace(self.depth.e_low, self.depth.e_front, 65)
        return max(float(np.max(sum(np.abs(v) for v in self.coefficients(initial, final, grid, e).values())))
                   for e in energies)

    def _table(self, initial: int, final: int) -> list:
        """[(grid of CM angles, {(k, q): a_kq × population}, population)] of each solution."""
        if (initial, final) not in self._tables:
            out = []
            for s in self.solutions:
                grid, coeff, _ = s.coefficient_table(initial, final)
                populated = np.array([st.population(initial) for st in s.table()[1]])
                out.append((grid, {key: np.asarray(v) * populated for key, v in coeff.items()}, populated))
            self._tables[(initial, final)] = out
        return self._tables[(initial, final)]


class Faded:
    """What deorientation in vacuum leaves of the orientation with which the nuclei of a level decay, averaged over their
    decays, for a solution ``ex`` with the mean life (ps; 0: decays at once) of each level of its scheme in
    ``lives_ps`` and the model's ``parameters``. It is what the simulated events average to.

    Feeding by radiation that is not observed leaves each rank k by itself: a level hands on F = (branch) × U_k of its
    tensors. With D the tensors the collision leaves, T the time at which the nucleus leaves the target, τ the mean
    lives and G_k of each level (:meth:`~physim.nuclear.deorientation.TwoState.average`), the levels hold at T

        N = e^{AT} D,    A = −(1 − F)/τ,

    and the tensors the nuclei decay with, summed over the decays before T (with the orientation the collision or the
    feeding gave) and after (each level deoriented while it lives, and what it feeds with it), are

        (1 − F)⁻¹ (1 − e^{AT}) D + (1 − G F)⁻¹ G e^{AT} D,

    which is (1 − F)⁻¹ D, the feeding of :class:`~physim.nuclear.orientation.Populated`, for G = 1 or a nucleus that never
    leaves the target. Levels that decay at once hand on what they get and are left out of N."""

    def __init__(self, ex, level: int, lives_ps, parameters):
        from . import angular
        from .deorientation import model

        levels = ex.scheme.levels
        above = [level]
        for n in above:
            above += [i for i, f in ex.transitions() if f == n and i not in above and levels[i].spin is not None]
        spins = [float(levels[n].spin) for n in above]
        self.ranks = [k for k in (2, 4, 6) if k <= round(2 * spins[0])]
        self.above = above
        tau = np.array([lives_ps[n] for n in above], dtype=float)
        self._timed, self._prompt = np.flatnonzero(tau > 0), np.flatnonzero(tau <= 0)
        self._tau = tau[self._timed]
        self._models = [model(ex.scheme.nuclide, spins[a], parameters) for a in self._timed]
        self._feeds = {}
        for k in self.ranks:
            feed = np.zeros((len(above), len(above)))
            for (i, f) in ex.transitions():
                if i not in above or f not in above:
                    continue
                d, a, b = ex.decay(i, f), above.index(i), above.index(f)
                parts = [(d.l1, 1.0)] if not d.l2 else [(d.l1, 1 / (1 + d.mixing**2)),
                                                        (d.l2, d.mixing**2 / (1 + d.mixing**2))]
                feed[b, a] = d.branch * sum(w * angular.u_coefficient(k, spins[a], spins[b], lam) for lam, w in parts
                                            if lam and angular.triangle(spins[a], lam, spins[b])
                                            and k <= round(2 * min(spins[a], spins[b])))
            self._feeds[k] = feed

    @property
    def off(self) -> bool:
        """Whether nothing is deoriented: no level on the way has a mean life, or the model is off."""
        return not len(self._timed) or not self.ranks or all(m.off for m in self._models)

    def at(self, left_ps: np.ndarray, beta: np.ndarray) -> dict:
        """{k: {level: factor of the tensors of rank k that the collision leaves that level with}} for nuclei that leave
        the target ``left_ps`` (ps; inf: never) after the collision at the speed ``beta``: arrays alike."""
        left, beta = np.asarray(left_ps, dtype=float), np.asarray(beta, dtype=float)
        n, timed, prompt = len(left), self._timed, self._prompt
        one = np.eye(len(timed))
        out_of = np.isfinite(left)
        g = [m.average(np.full(n, t), np.where(out_of, beta, 0.0)) for m, t in zip(self._models, self._tau)]
        out: dict = {}
        for k in self.ranks:
            feed = self._feeds[k]
            at_once = np.linalg.inv(np.eye(len(prompt)) - feed[np.ix_(prompt, prompt)])
            phi = (feed[np.ix_(timed, timed)]
                   + feed[np.ix_(timed, prompt)] @ at_once @ feed[np.ix_(prompt, timed)])
            # What the collision (or a level that decays at once) gives each level that lives a while.
            start = np.zeros((len(timed), len(self.above)))
            start[np.arange(len(timed)), timed] = 1.0
            start[:, prompt] = feed[np.ix_(timed, prompt)] @ at_once
            rate = -(one - phi) / self._tau
            still = np.zeros((n, len(timed), len(timed)))
            still[out_of] = _exp(rate * left[out_of][:, None, None])
            keep = np.stack([gk.get(k, np.ones(n)) for gk in g], axis=1)[:, :, None]
            decays = (np.linalg.inv(one - phi) @ (one - still) + np.linalg.solve(one - keep * phi, keep * still)) @ start
            if 0 in timed:
                row = decays[:, int(np.flatnonzero(timed == 0)[0]), :]
            else:
                mine = np.zeros((len(prompt), len(self.above)))
                mine[np.arange(len(prompt)), prompt] = 1.0
                row = (at_once @ (mine + feed[np.ix_(prompt, timed)] @ decays))[:, int(np.flatnonzero(prompt == 0)[0]), :]
            out[k] = {level: row[:, a] for a, level in enumerate(self.above)}
        return out


def _exp(a: np.ndarray) -> np.ndarray:
    """The exponentials of square matrices (n × m × m), by scaling and squaring."""
    norm = float(np.abs(a).sum(axis=-1).max(initial=0.0))
    halvings = max(0, math.ceil(math.log2(norm / 0.25))) if norm > 0 else 0
    x = a / 2.0**halvings
    out = np.eye(a.shape[-1]) + x
    term = x
    for n in range(2, 12):
        term = term @ x / n
        out = out + term
    for _ in range(halvings):
        out = out @ out
    return out


def _way_out(experiment, layers: list, depth: float, dirs: np.ndarray, emitter: str, energy: np.ndarray,
             mass_mev: float) -> tuple:
    """(the time (ps) at which nuclei made at ``depth`` with ``energy`` (MeV) along ``dirs`` leave the target and what
    is behind it, inf for one that stops in there; their speed outside, units of c)."""
    n = len(energy)
    at = np.full(n, depth)
    if any(not lay.material.density_g_cm3 for lay in layers):  # no clock without a density
        out = exit_energy(experiment, layers, 0, at, emitter, energy, dirs)
        return np.where(out > 0, 0.0, np.inf), _speed(out, mass_mev)
    from .gamma_events import _slowed

    s = _slowed(experiment, layers, at, dirs, emitter, energy, np.full(n, np.inf), mass_mev)
    return s.left_ps, s.beta


def _physical(rho: np.ndarray) -> np.ndarray:
    """Density matrices (n × d × d, trace 1) without eigenvalues below −:data:`NEGATIVE_EIGENVALUE`: those that have one
    get their negative eigenvalues set to zero and trace 1; the others are unchanged."""
    d = rho.shape[-1]
    try:
        np.linalg.cholesky(rho + NEGATIVE_EIGENVALUE * np.eye(d))  # much faster than the eigenvalues of every one
        return rho
    except np.linalg.LinAlgError:
        pass
    values, vectors = np.linalg.eigh(rho)
    bad = values[:, 0] < -NEGATIVE_EIGENVALUE
    if bad.any():
        kept = np.maximum(values[bad], 0.0)
        kept /= kept.sum(axis=1)[:, None]
        rho[bad] = (vectors[bad] * kept[:, None, :]) @ np.swapaxes(vectors[bad].conj(), 1, 2)
    return rho


def depth_orientation(experiment, solution=None) -> DepthOrientation:
    """The :class:`DepthOrientation` of the ``[reaction]`` nucleus through the target, for collisions on its main
    nuclide: with the solutions of its level scheme when the events follow one (``solution``, that of
    :func:`cascade_excitation`), else of the one state (:func:`excitation_of`), at the beam energies at which
    :func:`~physim.nuclear.rates.depth_probability` solves the excitation probability. A thin target
    (:data:`~physim.nuclear.rates.THIN_SPAN`), a lifetime measurement (its γ rays leave from wherever the plunger or
    the backing has the emitter) and a scheme that cannot be solved at one of those energies have the one solution
    at the middle of the target."""
    from . import data
    from .rates import channels, depth_probability

    mid = solution if solution is not None else excitation_of(experiment)
    layers = stack(experiment)
    target = data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
    channel = next((ch for ch in channels(layers, experiment)
                    if ch.excitation is not None and ch.layer == 0 and ch.nuclide.name == target), None)
    if channel is None or getattr(experiment, "lifetime", None) is not None:
        return DepthOrientation(None, [mid])
    depth = depth_probability(experiment, channel, layers)
    if depth.thin:
        return DepthOrientation(None, [mid])
    if solution is None:
        each = [excitation_of(experiment, energy=e) for e in depth.energies]
    else:
        each = [cascade_excitation(experiment, channel=channel, energy=e) for e in depth.energies]
    if any(s is None for s in each):
        return DepthOrientation(None, [mid])
    return DepthOrientation(depth, each)


#: Equal shares of the decays (by time) over which the emitter's speed is averaged when its level has a half-life.
DECAY_SHARES = 4
#: Depths (Gauss–Legendre) of a target whose middle the table takes over which the deorientation is averaged: how long
#: a nucleus takes to get out.
EXIT_DEPTHS = 3


def correlation_table(experiment, order: int = 8, crystal_points: int = 6, nucleus: str = "reaction",
                      transition: Optional[tuple] = None) -> list:
    """For every particle detector × γ-ray crystal: how many γ rays the crystal sees when the particle is seen in
    that detector, relative to isotropic emission (1 = isotropic).

    The factor is 4π W / yield in the laboratory, averaged over the crystal's face (``crystal_points``² points)
    and over the particle detector (``order`` points each way, :func:`_table_points`: equally spaced round an annulus,
    and cut and crowded where the cone of a particle's largest lab angle crosses the face, beyond which it does not go
    and up to which its cross section per lab solid angle goes as 1/√ of the distance), each point of which is weighted
    by the excitation cross section there. Both the scattered beam and the recoil reaching the detector are counted, each with the
    orbit it belongs to. It holds the particle–γ angular correlation and the forward throw of γ rays from a
    moving nucleus; with ``emission = "isotropic"`` only the second.

    The γ ray is that of the ``[reaction]`` state to the ground state. With a level scheme whose cascades the
    events follow (:func:`cascade_excitation`) the factor is from that scheme's solution, as the events are: the
    state with the orientation the other levels leave it, and the feeding from above. Without one it is from the
    state alone, in first order (:func:`excitation_of`).

    With deorientation in vacuum (:func:`deorientation_of`) the orientation is what the hyperfine field leaves of it,
    averaged over the decays as the events have them: each level with a half-life is deoriented (:class:`Faded`) while it
    lives outside the target, from where in the target the nucleus starts and with its own speed outside, and the levels
    above feed it with what they have left. A nucleus that decays in the target or stops there is not deoriented, nor is
    a level without a half-life.

    Each row also gives "detector_factor" for the γ-ray detector as a whole (its crystals averaged by the solid
    angle each covers). Shadowing is not applied.

    **Through the depth of the target** the factor is averaged as the rates and the events are
    (:func:`~physim.nuclear.rates.depth_probability`): at each depth the kinematics, the excitation cross section
    (the weight) and the orientation are those of the beam energy there, and the emitting nucleus slows down on its
    way out (to where it decays, if its level has a half-life), as the γ-ray events have it, of a single state and of
    a level scheme. (The feeding of the state is that of each depth too, where the events draw the level an excitation
    starts in from the populations in the middle of the target.) A thin target (the beam loses less than
    :data:`~physim.nuclear.rates.THIN_SPAN` of its energy) and a lifetime measurement (a ``[lifetime]`` section) keep
    the kinematics of the middle of the target, with the emitter at its speed just after the reaction; so does the
    table of the other nucleus, whose excitation the events take there too.

    ``nucleus`` "other" gives the table for the other nucleus of the collision
    (:func:`other_excitation`; empty without its level scheme): for the transition from its first excited level
    to the ground state, each point weighted by that level's population with the feeding from above.

    ``transition`` (initial, final), of a nucleus with a level scheme: the table for that one instead."""
    ex, solution, transition = _table_solution(experiment, nucleus, transition)
    if solution is None:
        return []
    # Kept with the solution, by the detectors, by the target the beam slows down in (the solution is that of its
    # middle, with the setup's E1 polarisation; the nuclei leave it, and its layers' densities say when) and by the
    # deorientation, which the solution does not have: the plan asks again at every change of anything else.
    scheme = experiment.levels.get("target" if experiment.excitation.excite == "target" else "beam")
    key = repr((order, crystal_points, transition, experiment.detectors, experiment.gamma_detectors,
                experiment.chamber if hasattr(experiment, "chamber") else None, experiment.beam.energy_mev,
                tilt_deg(experiment), [(lay.material, lay.thickness) for lay in stack(experiment)],
                getattr(experiment, "lifetime", None), deorientation_parameters(experiment, nucleus),
                None if scheme is None else json.dumps(scheme.to_dict(), sort_keys=True, default=str)))
    kept = solution.__dict__.setdefault("_correlation_tables", {})
    if key not in kept:
        if len(kept) > 8:
            kept.clear()
        kept[key] = _correlation_rows(experiment, ex, order, crystal_points, transition, nucleus)
    return [dict(row) for row in kept[key]]


def _orbit_axes(theta_cm: np.ndarray, phi_deg: np.ndarray, flip: bool) -> np.ndarray:
    """The orbit's frame in the laboratory of many collisions (n × 3 × 3, the axes as columns):
    :meth:`physim.nuclear.orientation.Populated.axes` for arrays. ``flip`` for projectile excitation."""
    th, ph = np.radians(theta_cm), np.radians(phi_deg)
    beam = np.array([0.0, 0.0, 1.0])
    out = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=-1)
    s, c = np.sin(th / 2)[:, None], np.cos(th / 2)[:, None]
    side = np.stack([np.cos(ph), np.sin(ph), np.zeros(len(th))], axis=-1)
    x = np.where(s > 1e-9, (out - beam) / (2 * np.maximum(s, 1e-12)), side)
    y = np.where(c > 1e-9, (out + beam) / (2 * np.maximum(c, 1e-12)), side)
    z = np.cross(x, y)
    if flip:
        x, y = -x, -y
    return np.stack([x, y, z], axis=-1)


def _speed(energy, mass_mev: float) -> np.ndarray:
    """Speed (units of c) of nuclei of kinetic energy ``energy`` (MeV) and mass ``mass_mev``."""
    return np.sqrt(np.maximum(1 - 1 / (1 + np.maximum(energy, 0.0) / mass_mev) ** 2, 0.0))


def _decay_speed(experiment, layers: list, depth: float, dirs: np.ndarray, emitter: str, energy: np.ndarray,
                 mass_mev: float, tau_s: float) -> np.ndarray:
    """The speed (units of c) at which nuclei made at ``depth`` with ``energy`` (MeV) along ``dirs`` decay, as the
    γ-ray events have it: after leaving the target, or for a level of mean life ``tau_s`` (s) averaged over its
    decays, each wherever on the way out it happens."""
    n = len(energy)
    at = np.full(n, depth)
    if not tau_s > 0:
        return _speed(exit_energy(experiment, layers, 0, at, emitter, energy, dirs), mass_mev)
    from .gamma_events import _slowed

    times = [-tau_s * math.log(1 - (j + 0.5) / DECAY_SHARES) for j in range(DECAY_SHARES)]
    return np.mean([_speed(_slowed(experiment, layers, at, dirs, emitter, energy, np.full(n, t), mass_mev).energy,
                           mass_mev) for t in times], axis=0)


def _faded(experiment, ex, solution, level: int, nucleus: str) -> Optional[Faded]:
    """The :class:`Faded` of a level of ``ex`` (``solution``, or the ``[reaction]`` state alone), or None when the setup
    deorients nothing of it."""
    from .gamma_events import _mean_lives

    parameters = deorientation_parameters(experiment, nucleus)
    if parameters is None or ex.isotropic or ex.scheme.levels[level].spin is None:
        return None
    lives = _mean_lives(experiment, solution, np.arange(len(ex.scheme.levels)))
    if lives is None:
        return None
    faded = Faded(ex, level, 1e12 * lives, parameters)
    return None if faded.off else faded


def _crossings(start: np.ndarray, along: np.ndarray, lo: float, hi: float, limits) -> list:
    """Where the line ``start`` + s ``along`` (mm) crosses the cones of polar angle ``limits`` (degrees) about the
    beam axis, for s strictly between ``lo`` and ``hi``: sorted."""
    out = []
    for limit in limits:
        c = math.cos(math.radians(limit))
        a = along[2] ** 2 - c * c * float(along @ along)
        b = 2 * (start[2] * along[2] - c * c * float(start @ along))
        k = start[2] ** 2 - c * c * float(start @ start)
        if abs(a) < 1e-12 * max(float(along @ along), 1.0):
            roots = [] if abs(b) < 1e-300 else [-k / b]
        else:
            d = b * b - 4 * a * k
            roots = [] if d < 0 else [(-b - math.sqrt(d)) / (2 * a), (-b + math.sqrt(d)) / (2 * a)]
        span = hi - lo
        for s in roots:
            z = start[2] + s * along[2]
            # On the cone itself, not on its mirror image behind the target.
            if lo + 1e-9 * span < s < hi - 1e-9 * span and z * c >= -1e-12 * abs(z):
                out.append(float(s))
    return sorted(out)


#: The share of a line's points that each of its pieces takes where a particle's largest lab angle cuts it.
CUT_SHARE = 5 / 8


def _line_nodes(lo: float, hi: float, cuts: list, order: int) -> tuple:
    """Gauss points and weights along [lo, hi], in pieces between the ``cuts``. A piece that ends at a cut has its
    points crowded quadratically to both ends, so a share that goes as 1/√ of the distance to a cut is integrated as a
    smooth one. A line without cuts has ``order`` points; one with cuts :data:`CUT_SHARE` of them in each piece (the
    pieces are smooth: at order 8 the table of the 98Ru example is then within 0.015 % of a fine one)."""
    if not cuts:
        x, w = _gauss(order)
        return (hi - lo) / 2 * x + (lo + hi) / 2, (hi - lo) / 2 * w
    x, w = _gauss(max(3, math.ceil(CUT_SHARE * order)))
    ends = [lo, *cuts, hi]
    s, ws = [], []
    for a, b in zip(ends[:-1], ends[1:]):
        t = math.pi * (x + 1) / 4
        s.append(a + (b - a) * np.sin(t) ** 2)
        ws.append((b - a) * math.pi / 4 * np.sin(2 * t) * w)
    return np.concatenate(s), np.concatenate(ws)


def _table_points(g, order: int, limits=(), segment: Optional[tuple] = None) -> tuple:
    """The points over a particle detector's face (or one ``segment`` of it) for :func:`correlation_table`: unit
    directions and the solid angle each stands for. An annulus (or disc) takes ``order`` points along each radius, on
    ``order`` radii: equally spaced round a whole ring (exact for the azimuthal harmonics of the angular correlation up
    to ``order`` − 1), at Gauss points across a sector; a rectangle ``order`` × ``order`` Gauss points. Where the cone of
    one of the ``limits`` (the largest lab angle of a particle whose kinematics are double-valued, degrees) crosses the
    face, each line of points is cut there and crowded to the cut (:func:`_line_nodes`): the particle's cross section
    per lab solid angle goes as 1/√ of the distance to it."""
    pts, area = [], []
    x, w = _gauss(order)
    if g.shape == "rectangle":
        u0, u1, v0, v1 = -g.width / 2, g.width / 2, -g.height / 2, g.height / 2
        if segment is not None:
            du, dv = g.width / g.strips_x, g.height / g.strips_y
            u0, v0 = u0 + segment[0] * du, v0 + segment[1] * dv
            u1, v1 = u0 + du, v0 + dv
        for xu, wu in zip(x, w):
            start = g.centre + ((u1 - u0) / 2 * xu + (u1 + u0) / 2) * g.u
            s, ws = _line_nodes(v0, v1, _crossings(start, g.v, v0, v1, limits), order)
            pts.append(start + s[:, None] * g.v)
            area.append((u1 - u0) / 2 * wu * ws)
    else:
        r0, r1, a0, a1 = g.inner_radius, g.outer_radius, 0.0, 2 * math.pi
        if segment is not None:
            dr, da = (r1 - r0) / g.rings, 2 * math.pi / g.sectors
            r0, a0 = r0 + segment[0] * dr, segment[1] * da
            r1, a1 = r0 + dr, a0 + da
        if a1 - a0 > 2 * math.pi - 1e-9:
            angles_, steps = 2 * math.pi * (np.arange(order) + 0.5) / order, np.full(order, 2 * math.pi / order)
        else:
            angles_, steps = (a1 - a0) / 2 * x + (a1 + a0) / 2, (a1 - a0) / 2 * w
        for phi, step in zip(angles_, steps):
            along = math.cos(phi) * g.u + math.sin(phi) * g.v
            s, ws = _line_nodes(r0, r1, _crossings(g.centre, along, r0, r1, limits), order)
            pts.append(g.centre + s[:, None] * along)
            area.append(step * s * ws)
    p, w = np.concatenate(pts), np.concatenate(area)
    r = np.linalg.norm(p, axis=1)
    return p / r[:, None], w * np.clip(-(p @ g.n) / r, 0, None) / r**2


def _limits_on(g, tb) -> tuple:
    """The largest lab angles (degrees) of the particles of ``tb`` (:class:`~physim.nuclear.kinematics.TwoBody`) whose
    kinematics are double-valued, where they lie on the face of ``g``."""
    lo, hi = g.theta_range()
    return tuple(sorted(a for a in (tb.max_angle(p) for p in ("ejectile", "recoil") if tb.double_valued(p))
                        if lo < a < hi))


def _mean_fade(parts: list) -> dict:
    """The mean of [(share, :meth:`Faded.at`)]."""
    return {k: {n: sum(w * f[k][n] for w, f in parts) for n in parts[0][1][k]} for k in parts[0][1]}


class _Integrand:
    """What :func:`correlation_table` integrates over a particle detector, for the ``transition`` (initial, final) of
    ``solution``, the solution of a level scheme (by default from its first excited level to the ground state), or
    without a solution for the ``[reaction]`` state alone: :attr:`steps`, the depths of the target with their share
    and kinematics, and :meth:`seen`, the weight of the particles along given directions and what each crystal sees
    of their γ rays."""

    def __init__(self, experiment, solution, crystal_points: int, transition=None, nucleus: str = "reaction"):
        from . import data
        from .coupled import CoupledChannels
        from .gamma_events import _mean_lives
        from .kinematics import TwoBody
        from .rates import rutherford_factor
        from .rutherford import Rutherford

        exc = experiment.excitation
        self.experiment = experiment
        if solution is not None:
            ex = solution
            levels = ex.scheme.levels
            if transition is None:
                lev_f = ex.ground
                lev_i = min((n for n in range(len(levels)) if n != lev_f), key=lambda n: levels[n].energy.value)
            else:
                lev_i, lev_f = transition
            # The collision leaves the nucleus in the initial level: its excitation energy sets the kinematics (not
            # the γ-ray energy, which differs for a line that does not go to the ground state).
            e_star = 1e-3 * (levels[lev_i].energy.value - levels[ex.ground].energy.value)
        else:
            ex = excitation_of(experiment)
            e_star = exc.energy_mev
            lev_i, lev_f = 1, 0
        self.ex, self.lev_i, self.lev_f = ex, lev_i, lev_f
        self.excite = "recoil" if ex.excite == "target" else "ejectile"
        self.flip = ex.excite == "projectile"
        self.m_emitter = data.nuclide(ex.scheme.nuclide).nuclear_mass_mev + e_star
        crystals = []
        for i, gd in enumerate(experiment.gamma_detectors):
            for label, centre, radius in gd.elements():
                d = math.sqrt(sum(c * c for c in centre))
                cover = (1 - math.cos(math.atan2(radius, d))) / 2
                crystals.append((gd.name or f"G{i + 1}", label,
                                 _disc_directions(tuple(c / d for c in centre), math.degrees(math.atan2(radius, d)),
                                                  crystal_points), cover))
        self.crystals = crystals
        self.all_dirs = np.concatenate([c[2] for c in crystals]) if crystals else np.zeros((0, 3))
        layers = self.layers = stack(experiment)
        orient = self.orient = (DepthOrientation(None, [ex]) if nucleus == "other"
                                else depth_orientation(experiment, solution))
        lives = None if orient.thin else _mean_lives(experiment, solution, np.array([lev_i]))
        self.tau = 0.0 if lives is None else float(lives[0])
        self.faded = _faded(experiment, ex, solution, lev_i, nucleus)
        x, wq = np.polynomial.legendre.leggauss(EXIT_DEPTHS)
        self.exits = list(zip((x + 1) / 2 * layers[0].thickness, wq / 2))
        # The depths (None: the middle of a thin target, the emitter at its speed just after the reaction), the share
        # of each and the beam energy there.
        points = [(None, 1.0, ex.beam_energy)]
        if not orient.thin:
            depth = orient.depth
            # As deep as the beam still has the energy to excite the nucleus.
            zs = np.linspace(0.0, depth.reached, 129)
            es = np.asarray(beam_energy_at(experiment, 0, zs, layers), dtype=float)
            deepest = float(zs[es >= depth.e_low].max(initial=zs[1]))
            x, wq = np.polynomial.legendre.leggauss(len(depth.energies))
            zs = (x + 1) / 2 * deepest
            points = list(zip(zs, wq, (float(e) for e in beam_energy_at(experiment, 0, zs, layers))))
        # (In first order only the levels reached in one step are excited; the coupled equations reach the others too.)
        reached = lev_i in ex.paths or isinstance(ex, CoupledChannels)
        if not reached or ex.decay(lev_i, lev_f).gamma <= 0:
            points = []
        self.steps = []
        for z, share, e_beam in points:
            populated = orient.population(lev_i, e_beam)
            if populated is None:
                continue
            try:
                tb = TwoBody(ex.beam, ex.target, e_beam, excitation_mev=e_star, excite=self.excite)
            except ValueError:  # below the threshold at this depth
                continue
            ruth = Rutherford(ex.beam, ex.target, e_beam)
            # P times Alder and Winther's inelastic Rutherford factor: the elastic one times v_i/v_f.
            self.steps.append((z, share * rutherford_factor(ruth.e_cm, e_star), e_beam, populated, tb, ruth))

    def seen(self, step: tuple, theta: np.ndarray, phi: np.ndarray, dom: np.ndarray, particle: str):
        """For ``particle`` ("ejectile" or "recoil") along the directions (θ, φ in degrees) of the face, each standing
        for the solid angle ``dom``, at one of :attr:`steps`: for each orbit that reaches them, (the directions it
        reaches, as indices, the weight of each, and 4π W in the laboratory averaged over each crystal's face: an
        (n, crystals) array)."""
        from .gamma_events import _directions, _distribution

        experiment, ex, layers, crystals, all_dirs = self.experiment, self.ex, self.layers, self.crystals, self.all_dirs
        z, share, e_beam, populated, tb, ruth = step
        for point in tb.at_lab(theta, particle):
            if point is None:
                continue
            # The CM angle and azimuth of the scattered beam on each orbit.
            th_p = np.asarray(point.theta_cm, dtype=float)
            jac = np.nan_to_num(np.asarray(point.jacobian, dtype=float))
            th_ej = th_p if particle == "ejectile" else 180.0 - th_p
            ph_ej = phi if particle == "ejectile" else phi + 180.0
            ok = ~np.isnan(th_p) & (jac > 0) & (th_ej > 0) & (th_ej <= 180)
            safe = np.where(ok, th_ej, 90.0)
            w = np.where(ok, share * np.asarray(populated.probability(safe), dtype=float)
                         * np.asarray(ruth.cross_section_cm(safe)) * jac * dom, 0.0)
            ok = w > 0
            if not ok.any():
                continue
            th_ej, ph_ej, w = th_ej[ok], ph_ej[ok], w[ok]
            em = tb.recoil_for(th_ej) if self.excite == "recoil" else tb.at_cm(th_ej)
            u = _directions(np.asarray(em.theta_lab), ph_ej + (180.0 if self.excite == "recoil" else 0.0))
            energy = np.nan_to_num(np.asarray(em.energy, dtype=float))
            beta = (_speed(energy, self.m_emitter) if z is None else
                    _decay_speed(experiment, layers, z, u, ex.scheme.nuclide, energy, self.m_emitter, self.tau))
            # 4π W in the laboratory along every direction into the crystals, from every orbit: at the
            # direction in the emitter's rest frame, times how much the motion gathers the γ rays along it.
            b = beta[:, None]
            cos_a = u @ all_dirs.T
            par = (cos_a - b) / (1 - b * cos_a)
            scale = np.sqrt(1 - b**2) / (1 - b * cos_a)
            rest = ((all_dirs[None, :, :] - cos_a[:, :, None] * u[:, None, :]) * scale[:, :, None]
                    + par[:, :, None] * u[:, None, :])
            local = np.einsum("nmi,nij->nmj", rest, _orbit_axes(th_ej, ph_ej, self.flip))
            fade = None
            if self.faded is not None:
                fade = _mean_fade([(share_z, self.faded.at(*_way_out(experiment, layers, at, u, ex.scheme.nuclide,
                                                                     energy, self.m_emitter)))
                                   for at, share_z in (self.exits if z is None else [(z, 1.0)])])
            a_kq = {key: v[:, None] for key, v in
                    self.orient.coefficients(self.lev_i, self.lev_f, th_ej, e_beam, fade).items()}
            seen = (_distribution(a_kq, np.arccos(np.clip(local[:, :, 2], -1.0, 1.0)),
                                  np.arctan2(local[:, :, 1], local[:, :, 0]))
                    * (1 - b**2) / (1 - b * cos_a) ** 2)
            yield np.flatnonzero(ok), w, seen.reshape(len(w), len(crystals), -1).mean(axis=2)

    def by_detector(self, factors: np.ndarray) -> tuple:
        """The γ-ray detectors' names, and each one's factor from those of its crystals (the last axis of
        ``factors``), averaged by the solid angle each crystal covers."""
        names = list(dict.fromkeys(c[0] for c in self.crystals))
        out = []
        for name in names:
            same = [j for j, x in enumerate(self.crystals) if x[0] == name]
            cover = np.array([self.crystals[j][3] for j in same])
            out.append(factors[..., same] @ cover / cover.sum())
        return names, np.stack(out, axis=-1) if out else np.zeros(factors.shape[:-1] + (0,))


def _correlation_rows(experiment, solution, order: int, crystal_points: int, transition=None,
                      nucleus: str = "reaction") -> list:
    """The rows of :func:`correlation_table`: for the ``transition`` (initial, final) of ``solution``, the
    solution of a level scheme (by default from its first excited level to the ground state), or without a
    solution for the ``[reaction]`` state alone."""
    it = _Integrand(experiment, solution, crystal_points, transition, nucleus)
    crystals = it.crystals
    rows = []
    for g in Array.from_experiment(experiment):
        faces = {}
        sums = np.zeros(len(crystals))
        weight = 0.0
        for step in it.steps:
            # The face's points, cut where a particle's largest lab angle at this beam energy lies on it.
            limits = _limits_on(g, step[4])
            if limits not in faces:
                dirs, dom = _table_points(g, order, limits)
                faces[limits] = (dom, *angles(dirs))
            dom, theta, phi = faces[limits]
            for particle in ("ejectile", "recoil"):
                for _, w, seen in it.seen(step, theta, phi, dom, particle):
                    sums += w @ seen
                    weight += float(w.sum())
        factors = sums / weight if weight > 0 else np.ones(len(crystals))
        for c, (name, label, _, cover) in enumerate(crystals):
            same = [j for j, x in enumerate(crystals) if x[0] == name]
            whole = sum(factors[j] * crystals[j][3] for j in same) / sum(crystals[j][3] for j in same)
            rows.append({"particle_detector": g.name, "gamma_detector": name, "crystal": label,
                         "factor": float(factors[c]), "detector_factor": float(whole)})
    return rows


#: Points along each side of a segment for :func:`segment_correlation`, at the least.
SEGMENT_TABLE_ORDER = 2
#: Points of the particle detectors taken at once by :func:`segment_correlation`.
SEGMENT_CHUNK = 512


def segment_correlation(experiment, parts: list, order: int = 8, crystal_points: int = 6, nucleus: str = "reaction",
                        transition: Optional[tuple] = None) -> tuple:
    """:func:`correlation_table`'s factor of each γ-ray detector for a part of a particle detector: ``parts`` are
    (particle detector, segment (i, j), whether the particle is the recoil). The factor is that of the particles of
    that kind in that segment alone, so a gate that keeps some segments, or one kind of particle, or the particles of
    some energy (which in one segment tells the kinds apart) has the factor of what it keeps. Each segment takes
    points enough that the face as a whole has about as many as the table (``order`` × ``order``), and
    :data:`SEGMENT_TABLE_ORDER` along each side at the least; they are cut and crowded where a particle's largest lab
    angle crosses the segment (:func:`_table_points`).

    Returns (the γ-ray detectors' names, an array of the factors, one row per part); a row is NaN where the model
    excites nothing with a particle of that kind in that segment."""
    ex, solution, transition = _table_solution(experiment, nucleus, transition)
    if ex is None and nucleus == "other":
        return [], np.full((len(parts), 0), np.nan)
    it = _Integrand(experiment, ex, crystal_points, transition, nucleus)
    sums, weight = np.zeros((len(parts), len(it.crystals))), np.zeros(len(parts))
    array = Array.from_experiment(experiment)
    for name in dict.fromkeys(p[0] for p in parts):
        g = array[name]
        n = max(SEGMENT_TABLE_ORDER, math.ceil(order / math.sqrt(len(g.segments))))
        for recoil in (False, True):
            mine = [k for k, p in enumerate(parts) if p[0] == name and bool(p[2]) == recoil]
            if not mine:
                continue
            particle = "recoil" if recoil else "ejectile"
            for step in it.steps:
                tb = step[4]
                limits = tuple(sorted(tb.max_angle(p) for p in ("ejectile", "recoil") if tb.double_valued(p)))
                pts = [_table_points(g, n, limits, parts[k][1]) for k in mine]
                dom = np.concatenate([p[1] for p in pts])
                owner = np.repeat(mine, [len(p[1]) for p in pts])
                theta, phi = angles(np.concatenate([p[0] for p in pts]))
                # In chunks: each point holds an array over every direction into the crystals.
                for a in range(0, len(dom), SEGMENT_CHUNK):
                    part = slice(a, a + SEGMENT_CHUNK)
                    for idx, w, seen in it.seen(step, theta[part], phi[part], dom[part], particle):
                        np.add.at(sums, owner[part][idx], w[:, None] * seen)
                        np.add.at(weight, owner[part][idx], w)
    with np.errstate(invalid="ignore", divide="ignore"):
        factors = sums / weight[:, None]
    return it.by_detector(factors)


def _table_solution(experiment, nucleus: str, asked) -> tuple:
    """For :func:`correlation_table`: the level scheme's solution its rows are of (None for the ``[reaction]`` state
    alone, or for the other nucleus without its scheme), the solution it is kept with, and the transition."""
    transition = None
    if nucleus == "other":
        ex = other_excitation(experiment)
        if ex is None:
            return None, None, None
    else:
        ex = cascade_excitation(experiment)
        if ex is not None:
            transition = state_transition(experiment)
    if asked is not None and ex is not None:
        transition = tuple(asked)
    return ex, ex if ex is not None else excitation_of(experiment), transition


__all__ = ["DepthOrientation", "Faded", "cascade_excitation", "correlation_table", "deorientation_of",
           "deorientation_parameters", "depth_orientation", "doppler_energy", "doppler_table", "e1_polarisation", "emitters",
           "excitation_of", "excitation_order", "has_static_moment", "level_names", "other_excitation", "other_role",
           "segment_correlation", "start_levels", "state_transition", "transition_lines", "unreached_levels"]
