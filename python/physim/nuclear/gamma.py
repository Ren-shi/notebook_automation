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

import numpy as np

from .detectors import Array, angles
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
    return _excitation(beam, target, e_mid, exc.excite, exc.energy_mev, exc.multipolarity, exc.b_up_e2fm, flat)


#: Accuracy of the coupled equations' probabilities in the events (relative, roughly; see
#: :class:`~physim.nuclear.coupled.CoupledChannels`): far below what any run's statistics show.
EVENT_TOLERANCE = 1e-4
#: The solutions kept (one per beam, partner nucleus, energy, level scheme and order), and how many.
_SOLVED: OrderedDict = OrderedDict()
_SOLVED_MOST = 48
_SOLVING = threading.Lock()


def excitation_order(experiment) -> str:
    """ "all" or "first": how the events populate the levels of a level scheme (``[reaction] order``)."""
    exc = experiment.excitation
    return (getattr(exc, "order", None) or "all") if exc is not None else "all"


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
    (:class:`~physim.nuclear.orientation.Excitation`: each level in one step from the ground state).

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
    key = (beam, target, e_mid, excite, bool(flat), order, json.dumps(scheme.to_dict(), sort_keys=True,
                                                                     default=str))
    with _SOLVING:  # one thread solves; the others then find the solution
        if key in _SOLVED:
            _SOLVED.move_to_end(key)
            return _SOLVED[key]
        try:
            if order == "first":
                ex = Excitation(beam, target, e_mid, scheme, excite=excite, isotropic=flat)
            else:
                ex = CoupledChannels(beam, target, e_mid, scheme, excite=excite, isotropic=flat,
                                     tolerance=EVENT_TOLERANCE)
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
                b_up: float, isotropic: bool):
    from .orientation import Excitation, simple_scheme

    scheme = simple_scheme(target if excite == "target" else beam, e_star, multipolarity, b_up)
    return Excitation(beam, target, e_mid, scheme, excite=excite, isotropic=isotropic)


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

    def coefficients(self, initial: int, final: int, theta_cm, energy) -> dict:
        """{(k, q): a_kq} of a transition (:meth:`~physim.nuclear.orientation.Populated.coefficients`) at each CM
        angle of ``theta_cm`` (degrees, an array) and the beam energy there (MeV, an array like it or a number)."""
        theta = np.asarray(theta_cm, dtype=float)
        out: dict = {}
        for (grid, tensors, populated), w in zip(self._table(initial, final), self._weights(theta, energy)):
            p = np.interp(theta, grid, populated)
            safe = np.where(p > 0, p, 1.0)
            for key, v in tensors.items():
                if key == (0, 0):
                    a = np.ones(len(theta), dtype=complex)
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


def correlation_table(experiment, order: int = 8, crystal_points: int = 6, nucleus: str = "reaction") -> list:
    """For every particle detector × γ-ray crystal: how many γ rays the crystal sees when the particle is seen in
    that detector, relative to isotropic emission (1 = isotropic).

    The factor is 4π W / yield in the laboratory, averaged over the crystal's face (``crystal_points``² points)
    and over the particle detector (``order`` points each way), each point of which is weighted by the excitation
    cross section there. Both the scattered beam and the recoil reaching the detector are counted, each with the
    orbit it belongs to. It holds the particle–γ angular correlation and the forward throw of γ rays from a
    moving nucleus; with ``emission = "isotropic"`` only the second.

    The γ ray is that of the ``[reaction]`` state to the ground state. With a level scheme whose cascades the
    events follow (:func:`cascade_excitation`) the factor is from that scheme's solution, as the events are: the
    state with the orientation the other levels leave it, and the feeding from above. Without one it is from the
    state alone, in first order (:func:`excitation_of`).

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
    to the ground state, each point weighted by that level's population with the feeding from above."""
    transition = None
    if nucleus == "other":
        ex = other_excitation(experiment)
        if ex is None:
            return []
    else:
        ex = cascade_excitation(experiment)
        if ex is not None:
            transition = state_transition(experiment)
    solution = ex if ex is not None else excitation_of(experiment)
    # Kept with the solution, by the detectors and by the target the beam slows down in (the solution is that of
    # its middle): the plan asks again at every change of anything else.
    scheme = experiment.levels.get("target" if experiment.excitation.excite == "target" else "beam")
    key = repr((order, crystal_points, transition, experiment.detectors, experiment.gamma_detectors,
                experiment.chamber if hasattr(experiment, "chamber") else None, experiment.beam.energy_mev,
                tilt_deg(experiment), [(lay.material, lay.thickness) for lay in stack(experiment)],
                getattr(experiment, "lifetime", None),
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


def _correlation_rows(experiment, solution, order: int, crystal_points: int, transition=None,
                      nucleus: str = "reaction") -> list:
    """The rows of :func:`correlation_table`: for the ``transition`` (initial, final) of ``solution``, the
    solution of a level scheme (by default from its first excited level to the ground state), or without a
    solution for the ``[reaction]`` state alone."""
    from . import data
    from .gamma_events import _directions, _distribution, _mean_lives
    from .kinematics import TwoBody
    from .rates import rutherford_factor
    from .rutherford import Rutherford

    exc = experiment.excitation
    if solution is not None:
        ex = solution
        levels = ex.scheme.levels
        if transition is None:
            lev_f = ex.ground
            lev_i = min((n for n in range(len(levels)) if n != lev_f), key=lambda n: levels[n].energy.value)
        else:
            lev_i, lev_f = transition
        e0 = 1e-3 * (levels[lev_i].energy.value - levels[lev_f].energy.value)
    else:
        ex = excitation_of(experiment)
        e0 = exc.energy_mev
        lev_i, lev_f = 1, 0
    excite = "recoil" if ex.excite == "target" else "ejectile"
    flip = ex.excite == "projectile"
    e_star = 1e-3 * (ex.scheme.levels[lev_i].energy.value - ex.scheme.levels[ex.ground].energy.value)
    m_emitter = data.nuclide(ex.scheme.nuclide).nuclear_mass_mev + e0
    crystals = []
    for i, gd in enumerate(experiment.gamma_detectors):
        for label, centre, radius in gd.elements():
            d = math.sqrt(sum(c * c for c in centre))
            cover = (1 - math.cos(math.atan2(radius, d))) / 2
            crystals.append((gd.name or f"G{i + 1}", label,
                             _disc_directions(tuple(c / d for c in centre), math.degrees(math.atan2(radius, d)),
                                              crystal_points), cover))
    layers = stack(experiment)
    orient = DepthOrientation(None, [ex]) if nucleus == "other" else depth_orientation(experiment, solution)
    lives = None if orient.thin else _mean_lives(experiment, solution, np.array([lev_i]))
    tau = 0.0 if lives is None else float(lives[0])
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
    if lev_i not in ex.paths or ex.decay(lev_i, lev_f).gamma <= 0:
        points = []
    steps = []
    for z, share, e_beam in points:
        populated = orient.population(lev_i, e_beam)
        if populated is None:
            continue
        try:
            tb = TwoBody(ex.beam, ex.target, e_beam, excitation_mev=e0, excite=excite)
        except ValueError:  # below the threshold at this depth
            continue
        ruth = Rutherford(ex.beam, ex.target, e_beam)
        # P times Alder and Winther's inelastic Rutherford factor: the elastic one times v_i/v_f.
        steps.append((z, share * rutherford_factor(ruth.e_cm, e_star), e_beam, populated, tb, ruth))
    rows = []
    all_dirs = np.concatenate([c[2] for c in crystals]) if crystals else np.zeros((0, 3))
    for g in Array.from_experiment(experiment):
        dirs, dom = g.directions(None, order)
        theta, phi = angles(dirs)
        sums = np.zeros(len(crystals))
        weight = 0.0
        for z, share, e_beam, populated, tb, ruth in steps:
            for particle in ("ejectile", "recoil"):
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
                    em = tb.recoil_for(th_ej) if excite == "recoil" else tb.at_cm(th_ej)
                    u = _directions(np.asarray(em.theta_lab), ph_ej + (180.0 if excite == "recoil" else 0.0))
                    energy = np.nan_to_num(np.asarray(em.energy, dtype=float))
                    beta = (_speed(energy, m_emitter) if z is None else
                            _decay_speed(experiment, layers, z, u, ex.scheme.nuclide, energy, m_emitter, tau))
                    # 4π W in the laboratory along every direction into the crystals, from every orbit: at the
                    # direction in the emitter's rest frame, times how much the motion gathers the γ rays along it.
                    b = beta[:, None]
                    cos_a = u @ all_dirs.T
                    par = (cos_a - b) / (1 - b * cos_a)
                    scale = np.sqrt(1 - b**2) / (1 - b * cos_a)
                    rest = ((all_dirs[None, :, :] - cos_a[:, :, None] * u[:, None, :]) * scale[:, :, None]
                            + par[:, :, None] * u[:, None, :])
                    local = np.einsum("nmi,nij->nmj", rest, _orbit_axes(th_ej, ph_ej, flip))
                    a_kq = {key: v[:, None] for key, v in
                            orient.coefficients(lev_i, lev_f, th_ej, e_beam).items()}
                    seen = (_distribution(a_kq, np.arccos(np.clip(local[:, :, 2], -1.0, 1.0)),
                                          np.arctan2(local[:, :, 1], local[:, :, 0]))
                            * (1 - b**2) / (1 - b * cos_a) ** 2)
                    sums += w @ seen.reshape(len(w), len(crystals), -1).mean(axis=2)
                    weight += float(w.sum())
        factors = sums / weight if weight > 0 else np.ones(len(crystals))
        for c, (name, label, _, cover) in enumerate(crystals):
            same = [j for j, x in enumerate(crystals) if x[0] == name]
            whole = sum(factors[j] * crystals[j][3] for j in same) / sum(crystals[j][3] for j in same)
            rows.append({"particle_detector": g.name, "gamma_detector": name, "crystal": label,
                         "factor": float(factors[c]), "detector_factor": float(whole)})
    return rows


__all__ = ["DepthOrientation", "cascade_excitation", "correlation_table", "depth_orientation",
           "doppler_energy", "doppler_table", "emitters",
           "excitation_of", "excitation_order", "has_static_moment", "level_names", "other_excitation", "other_role",
           "start_levels", "state_transition", "transition_lines", "unreached_levels"]
