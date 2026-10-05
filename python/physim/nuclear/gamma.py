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

import math
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


@lru_cache(maxsize=16)
def _excitation(beam: str, target: str, e_mid: float, excite: str, e_star: float, multipolarity: str,
                b_up: float, isotropic: bool):
    from .orientation import Excitation, simple_scheme

    scheme = simple_scheme(target if excite == "target" else beam, e_star, multipolarity, b_up)
    return Excitation(beam, target, e_mid, scheme, excite=excite, isotropic=isotropic)


def correlation_table(experiment, order: int = 4, crystal_points: int = 3) -> list:
    """For every particle detector × γ-ray crystal: how many γ rays the crystal sees when the particle is seen in
    that detector, relative to isotropic emission (1 = isotropic).

    The factor is 4π W / yield in the laboratory, averaged over the crystal's face and over the particle
    detector, each point of which is weighted by the excitation cross section there. Both the scattered beam and
    the recoil reaching the detector are counted, each with the orbit it belongs to. It holds the particle–γ
    angular correlation and the forward throw of γ rays from a moving nucleus; with ``emission = "isotropic"``
    only the second.

    Each row also gives "detector_factor" for the γ-ray detector as a whole (its crystals averaged by the solid
    angle each covers). The kinematics are those at the middle of the target; shadowing is not applied."""
    from . import data
    from .kinematics import TwoBody

    exc = experiment.excitation
    ex = excitation_of(experiment)
    e0 = exc.energy_mev
    excite = "recoil" if exc.excite == "target" else "ejectile"
    tb = TwoBody(ex.beam, ex.target, ex.beam_energy, excitation_mev=e0, excite=excite)
    m_emitter = data.nuclide(ex.scheme.nuclide).nuclear_mass_mev + e0
    cx = ex.paths[1][0][2] if 1 in ex.paths else None
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
                        state = ex.near(th_ej)
                        w = state.direct.get(1, 0.0) * float(cx.rutherford_cm(th_ej)) * jac[k] * dom[k]
                        if w <= 0:
                            continue
                        em = tb.recoil_for(th_ej) if excite == "recoil" else tb.at_cm(th_ej)
                        ph_em = math.radians(ph_ej + (180.0 if excite == "recoil" else 0.0))
                        th_em = math.radians(float(em.theta_lab))
                        gam = 1 + float(em.energy) / m_emitter
                        beta = math.sqrt(max(1 - 1 / gam**2, 0.0))
                        vel = beta * np.array([math.sin(th_em) * math.cos(ph_em), math.sin(th_em) * math.sin(ph_em),
                                               math.cos(th_em)])
                        y = state.gamma_yield(1, 0)
                        if y <= 0:
                            continue
                        seen = state.w_lab(1, 0, all_dirs, ph_ej, vel).reshape(len(crystals), -1).mean(axis=1)
                        sums += w * seen * 4 * math.pi / y
                        weight += w
        factors = sums / weight if weight > 0 else np.ones(len(crystals))
        for c, (name, label, _, cover) in enumerate(crystals):
            same = [j for j, x in enumerate(crystals) if x[0] == name]
            whole = sum(factors[j] * crystals[j][3] for j in same) / sum(crystals[j][3] for j in same)
            rows.append({"particle_detector": g.name, "gamma_detector": name, "crystal": label,
                         "factor": float(factors[c]), "detector_factor": float(whole)})
    return rows


__all__ = ["correlation_table", "doppler_energy", "doppler_table", "excitation_of"]
