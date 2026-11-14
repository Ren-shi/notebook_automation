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
from .rates import _q, beam_energy_at, beam_ion, exit_energy, stack

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
    weighted by the Coulomb-excitation cross section there; the kinematics are exact at each depth, the excitation
    probability is that at mid-target energy."""
    from . import data
    from .coulex import Coulex
    from .kinematics import TwoBody

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
    e_mid = float(beam_energy_at(experiment, 0, [lay.thickness / 2], layers)[0])
    cx = Coulex(beam, target, e_mid, excite=exc.excite, energy=e0, multipolarity=exc.multipolarity,
                b_up=exc.b_up_e2fm)
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
            tb = TwoBody(beam, target, float(e), excitation_mev=e0, excite=excite)
            first, _ = tb.at_lab(theta, "ejectile")
            th_cm = np.asarray(first.theta_cm)
            ok = ~np.isnan(th_cm)
            if not ok.any():
                continue
            sigma = np.where(ok, np.asarray(cx.cross_section_cm(np.where(ok, th_cm, 90.0)))
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


def excitation_of(experiment, isotropic=None):
    """The :class:`~physim.nuclear.orientation.Excitation` of the setup's excited state (the ``[reaction]``
    section: one state reached from a 0⁺ ground state), at the beam energy in the middle of the target.
    ``isotropic`` overrides the setup's ``emission`` switch."""
    from . import data

    exc = experiment.excitation
    if exc is None:
        raise ValueError("the setup has no excited state ([reaction] type = \"coulex\")")
    layers = stack(experiment)
    lay = layers[0]
    target = data.nuclide(max(lay.nuclides.items(), key=lambda kv: kv[1])[0]).name
    beam = beam_ion(experiment)
    # The same energy, rounded the same way, as the rates use: the orbit integrals are then shared.
    e_mid = round(float(beam_energy_at(experiment, 0, [lay.thickness / 2], layers)[0]), 9)
    flat = (exc.emission == "isotropic") if isotropic is None else isotropic
    return _excitation(beam, target, e_mid, exc.excite, exc.energy_mev, exc.multipolarity, exc.b_up_e2fm, flat)


#: Accuracy of the coupled equations' probabilities in the events (relative, roughly; see
#: :class:`~physim.nuclear.coupled.CoupledChannels`): far below what any run's statistics show.
EVENT_TOLERANCE = 1e-4
#: The solutions kept (one per beam, partner nucleus, energy, level scheme and order), and how many.
_SOLVED: OrderedDict = OrderedDict()
_SOLVED_MOST = 16
_SOLVING = threading.Lock()


def excitation_order(experiment) -> str:
    """ "all" or "first": how the events populate the levels of a level scheme (``[reaction] order``)."""
    exc = experiment.excitation
    return (getattr(exc, "order", None) or "all") if exc is not None else "all"


def deorientation_of(experiment) -> str:
    """ "none" or "vacuum": whether a nucleus that decays after leaving the target is deoriented
    (``[reaction] deorientation``, :mod:`physim.nuclear.deorientation`)."""
    exc = experiment.excitation
    return (getattr(exc, "deorientation", None) or "none") if exc is not None else "none"


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


def cascade_excitation(experiment, isotropic=None, order=None, channel=None, nucleus: str = "reaction"):
    """The excitation of the excited nucleus's whole level scheme (the setup's ``[levels]`` of the target or the
    beam, whichever the reaction excites), when its decays follow more than one transition: then each excitation
    starts in a level drawn from this solution's populations and decays by a cascade (backlog item 71). None
    when the setup has no such scheme, so a single state keeps the one-γ-ray chain of item 55.

    It is the coupled equations' solution to all orders (:class:`~physim.nuclear.coupled.CoupledChannels`: the
    levels reached in two steps, the reorientation, and the orientation of every level), or with ``order``
    "first" (the setup's ``[reaction] order`` if None) first-order theory
    (:class:`~physim.nuclear.orientation.Excitation`: each level in one step from the ground state).

    ``channel`` (a :class:`~physim.nuclear.rates.Channel` with an excitation) picks the nucleus the excited one
    scatters on and the layer; by default the target's main nuclide. The solution is at the beam energy in the
    middle of the layer, and is kept with its table of angles: the events, the rates and the γ rays of one setup
    share it.

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
    e_mid = round(float(beam_energy_at(experiment, layer, [lay.thickness / 2], layers)[0]), 9)
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
        # The [reaction] nucleus keeps its one-γ-ray chain for a single transition; the other nucleus has no
        # such chain, so one transition is enough for it.
        if ex is not None and (not ex.paths or len(ex.transitions()) < (2 if nucleus == "reaction" else 1)):
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


@lru_cache(maxsize=16)
def _excitation(beam: str, target: str, e_mid: float, excite: str, e_star: float, multipolarity: str,
                b_up: float, isotropic: bool):
    from .orientation import Excitation, simple_scheme

    scheme = simple_scheme(target if excite == "target" else beam, e_star, multipolarity, b_up)
    return Excitation(beam, target, e_mid, scheme, excite=excite, isotropic=isotropic)


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
    angle each covers). The kinematics are those at the middle of the target; shadowing is not applied.

    With ``[reaction] deorientation = "vacuum"`` every level with a half-life is deoriented before it decays, as
    in the events: its tensors are multiplied by G_k averaged over its decay times (GOSIA's time integral,
    :meth:`~physim.nuclear.deorientation.TwoState.integral`, at the speed the nucleus leaves the target with),
    times the chance that it decays after leaving the target, to which the chance that it decays inside adds 1;
    a nucleus that stops in the target is not deoriented. The levels above feed each level after their own
    attenuation (:class:`~physim.nuclear.deorientation.Faded`). Where the nucleus leaves is worked out at each
    CM angle of the solution's grid, from the middle of the target, along the azimuth 0.

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
    # Kept with the solution, by the detectors: the plan asks again at every change of anything else.
    fade = deorientation_of(experiment) == "vacuum"
    key = repr((order, crystal_points, transition, experiment.detectors, experiment.gamma_detectors,
                experiment.chamber if hasattr(experiment, "chamber") else None,
                (experiment.target, experiment.levels) if fade else None))
    kept = solution.__dict__.setdefault("_correlation_tables", {})
    if key not in kept:
        if len(kept) > 8:
            kept.clear()
        kept[key] = _correlation_rows(experiment, ex, order, crystal_points, transition, fade)
    return [dict(row) for row in kept[key]]


def level_lives(experiment, ex, single: bool) -> list:
    """The mean life (s) of each level of ``ex.scheme``, 0 without a half-life: from the scheme itself, or for
    the ``[reaction]`` state alone (``single``) from the setup's scheme level at the state's energy, as the
    events take it."""
    from .gamma_events import LN2

    if not single:
        return [lv.half_life.value / LN2 if lv.half_life is not None and lv.half_life.value > 0 else 0.0
                for lv in ex.scheme.levels]
    exc = experiment.excitation
    lives = [0.0] * len(ex.scheme.levels)
    scheme = experiment.levels.get("target" if exc.excite == "target" else "beam")
    for lv in (scheme.levels if scheme is not None else []):
        if abs(lv.energy.value - 1e3 * exc.energy_mev) < 1.0 and lv.half_life is not None and lv.half_life.value > 0:
            lives[1] = lv.half_life.value / LN2
            break
    return lives


def deorientation_factors(experiment, ex, excite: str, tb, m_emitter: float, lives: list, step: float = 2.0) -> dict:
    """{level: G_k on the CM-angle grid of ``ex.table(step)`` (angles × ranks)} for the levels with a mean life
    in ``lives``: the attenuation averaged over each level's decay, for a nucleus that starts in the middle of the
    target and leaves along the two-body direction at azimuth 0 (see :func:`correlation_table`)."""
    from .deorientation import model
    from .gamma_events import _slowing

    grid, _ = ex.table(step)
    layers = stack(experiment)
    point = tb.recoil_for(grid) if excite == "recoil" else tb.at_cm(grid, "ejectile")
    th = np.radians(np.nan_to_num(np.asarray(point.theta_lab, dtype=float)))
    energy = np.nan_to_num(np.asarray(point.energy, dtype=float))
    dirs = np.stack([np.sin(th), np.zeros_like(th), np.cos(th)], axis=-1)
    depth = np.full(len(grid), layers[0].thickness / 2)
    e_out, _, t_exit = _slowing(experiment, layers, depth, dirs, ex.scheme.nuclide, energy, np.full(len(grid), np.inf),
                                m_emitter, exits=True)
    gam = 1 + np.maximum(e_out, 0.0) / m_emitter
    beta = np.sqrt(np.maximum(1 - 1 / gam**2, 0.0))
    left = np.isfinite(t_exit) & (energy > 0)
    out = {}
    for n, tau in enumerate(lives):
        spin = ex.scheme.levels[n].spin
        if tau <= 0 or spin is None or spin <= 0:
            continue
        g = model(ex.scheme.nuclide, float(spin)).integral(np.full(len(grid), tau * 1e12), beta)
        after = np.where(left, np.exp(-np.where(left, t_exit, 0.0) / tau), 0.0)[:, None]
        out[n] = (1 - after) + after * g
    return out


def _correlation_rows(experiment, solution, order: int, crystal_points: int, transition=None,
                      fade: bool = False) -> list:
    """The rows of :func:`correlation_table`: for the ``transition`` (initial, final) of ``solution``, the
    solution of a level scheme (by default from its first excited level to the ground state), or without a
    solution for the ``[reaction]`` state alone; with ``fade``, deoriented in vacuum."""
    from . import data
    from .kinematics import TwoBody

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
        excite = "recoil" if ex.excite == "target" else "ejectile"

        def weigh(state):
            return state.population(lev_i)
    else:
        ex = excitation_of(experiment)
        e0 = exc.energy_mev
        excite = "recoil" if exc.excite == "target" else "ejectile"
        lev_i, lev_f = 1, 0

        def weigh(state):
            return state.direct.get(1, 0.0)
    tb = TwoBody(ex.beam, ex.target, ex.beam_energy, excitation_mev=e0, excite=excite)
    m_emitter = data.nuclide(ex.scheme.nuclide).nuclear_mass_mev + e0
    cx = ex.paths[lev_i][0][2] if lev_i in ex.paths else None
    states_of = ex.near
    if fade and not ex.isotropic:
        from .deorientation import Faded

        factors = deorientation_factors(experiment, ex, excite, tb, m_emitter, level_lives(experiment, ex,
                                                                                         solution is None))
        if factors:
            states_of = Faded(ex, factors).near
    crystals = []
    for i, gd in enumerate(experiment.gamma_detectors):
        for label, centre, radius in gd.elements():
            d = math.sqrt(sum(c * c for c in centre))
            cover = (1 - math.cos(math.atan2(radius, d))) / 2
            crystals.append((gd.name or f"G{i + 1}", label,
                             _disc_directions(tuple(c / d for c in centre), math.degrees(math.atan2(radius, d)),
                                              crystal_points), cover))
    rows = []
    all_dirs = np.concatenate([c[2] for c in crystals]) if crystals else np.zeros((0, 3))
    for g in Array.from_experiment(experiment):
        dirs, dom = g.directions(None, order)
        theta, phi = angles(dirs)
        sums = np.zeros(len(crystals))
        weight = 0.0
        if cx is not None:
            for particle in ("ejectile", "recoil"):
                for point in tb.at_lab(theta, particle):
                    if point is None:
                        continue
                    th_p = np.asarray(point.theta_cm, dtype=float)
                    jac = np.nan_to_num(np.asarray(point.jacobian, dtype=float))
                    for k in range(len(theta)):
                        if np.isnan(th_p[k]) or jac[k] <= 0:
                            continue
                        # The CM angle and azimuth of the scattered beam on this orbit.
                        th_ej = th_p[k] if particle == "ejectile" else 180.0 - th_p[k]
                        ph_ej = phi[k] if particle == "ejectile" else phi[k] + 180.0
                        if not 0 < th_ej <= 180:
                            continue
                        state = states_of(th_ej)
                        w = weigh(state) * float(cx.rutherford_cm(th_ej)) * jac[k] * dom[k]
                        if w <= 0:
                            continue
                        em = tb.recoil_for(th_ej) if excite == "recoil" else tb.at_cm(th_ej)
                        ph_em = math.radians(ph_ej + (180.0 if excite == "recoil" else 0.0))
                        th_em = math.radians(float(em.theta_lab))
                        gam = 1 + float(em.energy) / m_emitter
                        beta = math.sqrt(max(1 - 1 / gam**2, 0.0))
                        vel = beta * np.array([math.sin(th_em) * math.cos(ph_em), math.sin(th_em) * math.sin(ph_em),
                                               math.cos(th_em)])
                        y = state.gamma_yield(lev_i, lev_f)
                        if y <= 0:
                            continue
                        seen = state.w_lab(lev_i, lev_f, all_dirs, ph_ej, vel).reshape(len(crystals), -1).mean(axis=1)
                        sums += w * seen * 4 * math.pi / y
                        weight += w
        factors = sums / weight if weight > 0 else np.ones(len(crystals))
        for c, (name, label, _, cover) in enumerate(crystals):
            same = [j for j, x in enumerate(crystals) if x[0] == name]
            whole = sum(factors[j] * crystals[j][3] for j in same) / sum(crystals[j][3] for j in same)
            rows.append({"particle_detector": g.name, "gamma_detector": name, "crystal": label,
                         "factor": float(factors[c]), "detector_factor": float(whole)})
    return rows


__all__ = ["cascade_excitation", "correlation_table", "doppler_energy", "doppler_table", "emitters",
           "excitation_of", "excitation_order", "level_names", "other_excitation", "other_role", "start_levels",
           "state_transition", "transition_lines", "unreached_levels"]
