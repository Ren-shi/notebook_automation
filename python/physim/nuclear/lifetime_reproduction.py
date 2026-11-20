"""Published lifetime measurements reproduced with physim, and the simulated runs of the lifetime validation checks.

Each reproduction builds the experiment from physim's example after the paper (the setup, the level scheme, the
lifetimes the paper found), simulates it with physim's own event generators and transport, analyses the simulated
spectra the way the paper analysed its data with :mod:`physim.nuclear.lifetime_analysis`, and compares the τ it gets
with the paper's (``tests/reference/nuclear/lifetimes/published/``; the yardstick is the git-ignored
``planning/lifetimes-validation.md``)::

    from physim.nuclear import lifetime_reproduction as lr

    lr.ru98()     # Radeck et al., PRC 85, 014301 (2012): RDDS after Coulex, τ(2⁺₁) with the deorientation correction
    lr.cr50()     # Beckers et al., NIM A 1042, 167418 (2022): three-foil RDDS after fusion, γγ DDCM (their Eq. 13)
    lr.sm143()    # Rajbanshi et al., PLB 782, 143 (2018): DSAM after fusion on a Au backing, line-shape fits

The simulation starts from the paper's lifetimes (a reproduction cannot know τ otherwise), so what is checked is
the chain: the paper's setup, physim's simulation of it (kinematics, slowing, the cascade's timing, deorientation)
and physim's analysis done as the paper did it return the paper's τ within the paper's uncertainty and the
simulation's statistics, with the paper's own corrections (the deorientation of ⁹⁸Ru's 2⁺₁, the gates of the
three-foil DDCM, SRIM stopping and side feeding for ¹⁴³Sm). Where a paper does not give a number the simulation
needs, the value used is marked illustrative here and in the examples.

The run sizes default to what the validation suite can afford (seconds per reproduction); ``scripts/
reproduce_lifetimes.py`` runs them larger. Results are dicts of plain numbers (and lists) with the comparison:
``passed``, ``worst`` (the largest deviation in units of the allowed one) and ``detail`` (a sentence).

The NaPyTau comparison and the one-level RDDS closure read ``tests/reference/nuclear/lifetimes/napytau/
ddcm_simulated.csv``: peak areas from a simulated plunger run of the ⁹⁸Ru example without deorientation
(:func:`napytau_input`, written by ``scripts/lifetime_reference.py napytau-simulated`` together with NaPyTau's
answers on them).
"""

from __future__ import annotations

import copy
import dataclasses
import math
from contextlib import contextmanager
from pathlib import Path
from typing import Optional

import numpy as np

from . import data
from . import lifetime_analysis as la
from . import lifetimes as lt

#: Where the published values are, inside the reference folder.
PUBLISHED = "lifetimes/published"
#: The simulated plunger series that NaPyTau analysed, inside the reference folder.
NAPYTAU_SIMULATED = "lifetimes/napytau/ddcm_simulated.csv"
#: The NaPyTau polynomial degrees compared.
NAPYTAU_DEGREES = (2, 3, 4)

_CACHE: dict = {}


def _reference(ref) -> Optional[Path]:
    from .validation import reference_dir

    return Path(ref) if ref is not None else reference_dir()


def published(ref=None) -> Optional[dict]:
    """The published lifetimes (``lifetimes.csv``: rows as dicts with floats) and setups (``setups.toml``) of the
    reproduced experiments; None if the files are not there."""
    from ._toml import loads
    from .validation import read_reference

    base = _reference(ref)
    if base is None or not (base / PUBLISHED / "lifetimes.csv").exists():
        return None
    rows = []
    for r in read_reference(base / PUBLISHED / "lifetimes.csv")[1]:
        row = dict(r)
        for k in ("tau_ps", "plus_ps", "minus_ps", "e_gamma_kev"):
            row[k] = float(r[k]) if r[k] not in ("", None) else None
        rows.append(row)
    setups = loads((base / PUBLISHED / "setups.toml").read_text(encoding="utf-8"))
    return {"lifetimes": rows, "setups": setups}


def _paper(pub: dict, case: str, state: str, role: str = "result") -> dict:
    found = [r for r in pub["lifetimes"] if r["case"] == case and r["state"] == state and r["role"] == role]
    if not found:
        raise KeyError(f"no published {role} for {case} {state}")
    return found[0]


def _sigma(row: dict, deviation: float) -> float:
    """The paper's uncertainty on the side of the deviation (plus when physim is above)."""
    return float(row["plus_ps"] if deviation > 0 else row["minus_ps"])


def _quantity(text) -> float:
    return float(str(text).split()[0])


# ---------------------------------------------------------------------------------------------------------------
# ⁹⁸Ru: RDDS after Coulex (Radeck et al. 2012)

#: The detector rings analysed (the example's clovers at 41.5° and 138.5°; the paper used its rings 1 and 6 for the
#: weaker states and all six for the 2⁺₁).
RU98_RINGS = (41.5, 138.5)
#: The paper's offset of the zero distance (electrical contact), µm.
RU98_OFFSET_UM = 8.0


def coulex_alignment(experiment, particle: str = "recoil", order: int = 8) -> dict:
    """The unattenuated angular distribution {2: a₂, 4: a₄} of the γ rays of a one-state Coulex setup in
    coincidence with its particle detectors, from physim's Coulex solution, the way
    :func:`physim.nuclear.gamma.correlation_table` works it out but counting only one of the two particles in the
    detectors: ``particle`` "recoil" (the target nucleus, all that reaches a detector behind a plunger's stopper) or
    "ejectile" ("both" for both, as correlation_table). The particle–γ correlation W at each crystal (4π W/yield,
    averaged over the crystal's face and the particle detector, weighted by the excitation cross section) over its
    value for isotropic emission (only the forward throw of the moving nucleus) is W(θ); averaged per ring,
    W = 1 + a₂P₂ + a₄P₄ at the forward and backward rings (averaged) and at 90° gives a₂ and a₄. What the paper took
    from its Coulex calculation for the deorientation correction. Also ``w``, W per ring."""
    from . import gamma
    from .detectors import Array, angles
    from .kinematics import TwoBody

    ex = gamma.excitation_of(experiment)
    exc = experiment.excitation
    e0 = exc.energy_mev
    excite = "recoil" if exc.excite == "target" else "ejectile"
    tb = TwoBody(ex.beam, ex.target, ex.beam_energy, excitation_mev=e0, excite=excite)
    m_emitter = data.nuclide(ex.scheme.nuclide).nuclear_mass_mev + e0
    cx = ex.paths[1][0][2]
    crystals = []
    for gd in experiment.gamma_detectors:
        for _label, centre, radius in gd.elements():
            d = math.sqrt(sum(c * c for c in centre))
            crystals.append((_quantity(gd.theta), gamma._disc_directions(tuple(c / d for c in centre),
                                                                          math.degrees(math.atan2(radius, d)), 6)))
    dirs_all = np.concatenate([c[1] for c in crystals])
    particles = ("ejectile", "recoil") if particle == "both" else (particle,)
    corr, flat = np.zeros(len(crystals)), np.zeros(len(crystals))
    for g in Array.from_experiment(experiment):
        dirs, dom = g.directions(None, order)
        theta, phi = angles(dirs)
        for which in particles:
            for point in tb.at_lab(theta, which):
                if point is None:
                    continue
                th_p = np.asarray(point.theta_cm, dtype=float)
                jac = np.nan_to_num(np.asarray(point.jacobian, dtype=float))
                for k in range(len(theta)):
                    if np.isnan(th_p[k]) or jac[k] <= 0:
                        continue
                    th_ej = th_p[k] if which == "ejectile" else 180.0 - th_p[k]
                    ph_ej = phi[k] if which == "ejectile" else phi[k] + 180.0
                    if not 0 < th_ej <= 180:
                        continue
                    state = ex.near(th_ej)
                    w = state.direct.get(1, 0.0) * float(cx.rutherford_cm(th_ej)) * jac[k] * dom[k]
                    y = state.gamma_yield(1, 0)
                    if w <= 0 or y <= 0:
                        continue
                    em = tb.recoil_for(th_ej) if excite == "recoil" else tb.at_cm(th_ej)
                    ph_em = math.radians(ph_ej + (180.0 if excite == "recoil" else 0.0))
                    th_em = math.radians(float(em.theta_lab))
                    gam = 1 + float(em.energy) / m_emitter
                    beta = math.sqrt(max(1 - 1 / gam**2, 0.0))
                    vel = beta * np.array([math.sin(th_em) * math.cos(ph_em), math.sin(th_em) * math.sin(ph_em),
                                           math.cos(th_em)])
                    seen = state.w_lab(1, 0, dirs_all, ph_ej, vel).reshape(len(crystals), -1).mean(axis=1)
                    # isotropic emission in the nucleus's frame: the forward throw alone
                    g2 = 1 / (1 - beta**2)
                    cosv = dirs_all @ (vel / beta) if beta > 0 else np.zeros(len(dirs_all))
                    boost = (1 / (g2 * (1 - beta * cosv) ** 2)).reshape(len(crystals), -1).mean(axis=1)
                    corr += w * seen * 4 * math.pi / y
                    flat += w * boost
    rings: dict = {}
    for (th, _d), c, f in zip(crystals, corr, flat):
        rings.setdefault(round(th, 3), []).append(c / f)
    w = {th: float(np.mean(v)) for th, v in rings.items()}
    off = [th for th in w if abs(th - 90.0) > 1.0]
    if not off or not any(abs(th - 90.0) <= 1.0 for th in w):
        raise ValueError("the alignment needs crystals at 90° and at one other angle")
    c = math.cos(math.radians(min(off)))
    w_off = float(np.mean([w[th] for th in off]))
    w_90 = float(np.mean([w[th] for th in w if abs(th - 90.0) <= 1.0]))
    a = np.array([[_p2(c), _p4(c)], [_p2(0.0), _p4(0.0)]])
    a2, a4 = np.linalg.solve(a, np.array([w_off - 1.0, w_90 - 1.0]))
    return {2: float(a2), 4: float(a4), "w": w}


def _p2(c):
    return 0.5 * (3 * c * c - 1)


def _p4(c):
    return (35 * c**4 - 30 * c * c + 3) / 8


def _ring_spectra(experiment, g, rings, edges, mask) -> dict:
    """Per ring: the weighted spectrum of the γ rays in ``mask`` (full detector response), scaled to effective
    counts (Σw/Σw² per γ ray, so that its errors are those of counts)."""
    theta = {gd.name: _quantity(gd.theta) for gd in experiment.gamma_detectors}
    detector = np.array([name.split(" ")[0] for name in g.crystal_names])[g["crystal"]]
    out = {}
    for r in rings:
        names = [n for n, th in theta.items() if abs(th - r) < 0.05]
        m = mask & np.isin(detector, names)
        w = g["weight"][m]
        scale = float(w.sum() / (w * w).sum()) if w.size else 0.0
        out[r] = np.histogram(1e3 * g["measured"][m], edges, weights=w)[0] * scale
    return out


def coulex_plunger(experiment, distances_um, *, particles: int = 2_000_000, seed: int = 1, rings=RU98_RINGS,
                   gammas_per_event: int = 1, bin_kev: float = 0.5, gate: str = "recoil") -> dict:
    """A simulated plunger run of a Coulex setup at each distance: the particle events once
    (:func:`physim.nuclear.events.simulate`, ``particles`` reactions), kept when the detected particle is the
    target recoil (``gate="recoil"``: what reaches a detector behind the stopper, which takes the scattered beam), and
    the γ rays at each distance from the same particles (:func:`physim.nuclear.gamma_events.simulate_gammas`, its own
    seed per distance), as spectra per detector ring in effective counts. Also the shifted peak's centroid at each
    ring at the largest distance and, for reference, the true unshifted share (decays in the stopper)."""
    from . import events as particle_events
    from . import gamma_events

    ev = particle_events.simulate(experiment, int(particles), int(seed))
    if gate == "recoil":
        keep = np.asarray(ev["recoil"], dtype=bool)
        ev = dataclasses.replace(ev, columns={k: v[keep] for k, v in ev.columns.items()})
    e0 = 1e3 * experiment.excitation.energy_mev
    edges = np.arange(e0 - 45.0, e0 + 45.0 + 1e-9, bin_kev)
    spectra = {r: [] for r in rings}
    truth, n_gamma = [], []
    last = None
    for i, d in enumerate(distances_um):
        g = gamma_events.simulate_gammas(experiment, particle_events=ev, seed=int(seed) + 101 * (i + 1),
                                         distance_um=float(d), gammas_per_event=int(gammas_per_event))
        ok = np.asarray(g["counted"], dtype=bool)
        for r, h in _ring_spectra(experiment, g, rings, edges, ok).items():
            spectra[r].append(h)
        stopper = g.segments.index("stopper") if "stopper" in g.segments else len(g.segments) - 1
        w = g["weight"]
        truth.append(float(np.sum(w[g["segment"] >= stopper]) / np.sum(w)))
        n_gamma.append(int(ok.sum()))
        last = g
    spectra = {r: np.array(v) for r, v in spectra.items()}
    # The shifted line's centroid at the largest distance, where it holds nearly all decays: the counts from
    # 2.5 FWHM beyond E₀ out to its kinematic end (the recoil's speed with 30% room) less a flat background read
    # below the line on the far side.
    beta_guess = _recoil_beta(experiment)
    centres = 0.5 * (edges[1:] + edges[:-1])
    fwhm = 2.1
    centroid = {}
    for r in rings:
        c = math.cos(math.radians(r))
        sign = 1.0 if c >= 0 else -1.0
        off = sign * (centres - e0)
        win = (off > 2.5 * fwhm) & (off <= 1.3 * beta_guess * abs(c) * e0 + 3 * fwhm)
        bg = (-off >= 2.5 * fwhm + 10.0) & (-off <= 2.5 * fwhm + 40.0)
        y = spectra[r][-1]
        level = float(y[bg].mean()) if bg.any() else 0.0
        centroid[r] = float(np.sum((y[win] - level) * centres[win]) / np.sum(y[win] - level))
    return {"distances_um": [float(d) for d in distances_um], "edges_kev": edges, "spectra": spectra, "e0_kev": e0,
            "centroid_kev": centroid, "cos_ring": ring_cos(experiment, rings), "true_unshifted": truth, "gammas": n_gamma, "particles": len(ev),
            "segments": list(last.segments) if last is not None else []}


def ring_cos(experiment, rings) -> dict:
    """The mean cos θ over the faces of each ring's crystals (θ from the beam): what a Doppler shift at the ring
    measures, a little less in magnitude than the cosine of the ring's nominal angle for crystals seen over ±13°."""
    from . import gamma

    out = {}
    for r in rings:
        vals = []
        for gd in experiment.gamma_detectors:
            if abs(_quantity(gd.theta) - r) > 0.05:
                continue
            for _label, centre, radius in gd.elements():
                d = math.sqrt(sum(c * c for c in centre))
                dirs = gamma._disc_directions(tuple(c / d for c in centre), math.degrees(math.atan2(radius, d)), 6)
                vals.append(float(np.mean(dirs[:, 2])))
        out[r] = float(np.mean(vals))
    return out


def _recoil_beta(experiment) -> float:
    nuclide, _made, out = experiment.lifetime_recoil()
    mass = data.nuclide(nuclide).nuclear_mass_mev
    return math.sqrt(1 - (mass / (mass + out)) ** 2)


def beta_from_centroids(run: dict) -> float:
    """The recoil's speed along the beam from the shifted line's centroids at a forward and a backward ring (first
    order, as the paper measured its v/c = 4.50(6)%): β = (E_f − E_b) / (E₀ (cos θ_f − cos θ_b)), with each ring's
    mean cos θ over its crystals' faces (:func:`ring_cos`)."""
    rings = sorted(run["centroid_kev"])
    f, b = rings[0], rings[-1]
    cos = run.get("cos_ring") or {r: math.cos(math.radians(r)) for r in rings}
    cf, cb = cos[f], cos[b]
    return (run["centroid_kev"][f] - run["centroid_kev"][b]) / (run["e0_kev"] * (cf - cb))


def plunger_areas(run: dict, ring: float, beta: float, fwhm_kev: float = 2.1, method: str = "window") -> dict:
    """The unshifted and shifted intensities at one ring at every distance.

    ``method="window"`` (default): the counts in the stopped peak's window, E₀ ± 2.5 FWHM, and in the rest of the
    line out to its kinematic end, E₀ (1 + β |cos θ|) + 3 FWHM on the shifted side (β the speed from the centroids,
    with 25% room for the spread of speeds), each less a flat background read from 10–40 keV below E₀ − 2.5 FWHM
    on the far side. The decays in the stopper while the recoil still slows (a tail between the peaks) fall in the
    shifted window: a constant share of the unshifted decays, which a free normalisation or the efficiency ratio
    (:func:`efficiency_ratio`) takes up. ``method="fit"``: :func:`lifetime_analysis.peak_areas` (Gaussians on a
    straight line), which misplaces counts when the shifted peak is broad and flat-topped, as with clovers seen
    over ±13° around 41.5° (the reason the window is the default here)."""
    edges = run["edges_kev"]
    e0 = run["e0_kev"]
    if method == "fit":
        areas = la.peak_areas(list(run["spectra"][ring]), e0, beta, ring, fwhm_kev, edges_kev=edges,
                              fwhm_shifted_kev=4.0)
        s, u, err = (np.array(x, dtype=float) for x in zip(*areas))
        return {"unshifted": u, "shifted": s, "err_unshifted": err[:, 1], "err_shifted": err[:, 0]}
    if method != "window":
        raise ValueError("method is 'window' or 'fit'")
    centres = 0.5 * (edges[1:] + edges[:-1])
    width = np.diff(edges)
    sign = 1.0 if math.cos(math.radians(ring)) >= 0 else -1.0
    half = 2.5 * fwhm_kev
    far = e0 + sign * 1.25 * beta * abs(math.cos(math.radians(ring))) * e0 + sign * 3 * fwhm_kev
    un = np.abs(centres - e0) <= half
    sh = (sign * (centres - e0) > half) & (sign * (centres - far) <= 0)
    bg = (-sign * (centres - e0) >= half + 10.0) & (-sign * (centres - e0) <= half + 40.0)
    spectra = np.asarray(run["spectra"][ring], dtype=float)
    level = spectra[:, bg].sum(axis=1) / width[bg].sum()  # counts per keV
    n_u, n_s = spectra[:, un].sum(axis=1), spectra[:, sh].sum(axis=1)
    b_u, b_s = level * width[un].sum(), level * width[sh].sum()
    var_level = spectra[:, bg].sum(axis=1) / width[bg].sum() ** 2
    return {"unshifted": n_u - b_u, "shifted": n_s - b_s,
            "err_unshifted": np.sqrt(np.maximum(n_u, 1.0) + var_level * width[un].sum() ** 2),
            "err_shifted": np.sqrt(np.maximum(n_s, 1.0) + var_level * width[sh].sum() ** 2)}


def efficiency_ratio(unshifted, shifted, err_unshifted, err_shifted) -> float:
    """The detection efficiency of the shifted component over the unshifted one, k, from runs with the same number
    of reactions at every distance: I_u + I_s/k is then the same at every distance (least squares, weights from the
    errors). It holds the Lorentz boost of the solid angle, the efficiency at the two energies and the stopper's
    slowing tail counted as shifted: the paper's "efficiency of shifted vs unshifted" correction."""
    u, s = np.asarray(unshifted, dtype=float), np.asarray(shifted, dtype=float)
    w = 1.0 / (np.asarray(err_unshifted, dtype=float) ** 2 + np.asarray(err_shifted, dtype=float) ** 2)
    # u = c − q s: weighted linear fit of u against s
    a = np.c_[np.ones_like(s), -s] * np.sqrt(w)[:, None]
    c, q = np.linalg.lstsq(a, u * np.sqrt(w), rcond=None)[0]
    return float(1.0 / q)


def ru98(ref=None, *, particles: int = 2_000_000, seed: int = 1, gammas_per_event: int = 1) -> dict:
    """Radeck et al.'s ⁹⁸Ru 2⁺₁ lifetime from a simulated plunger run of ``coulex_ru98_plunger``: the paper's
    distances (3–450 µm; those between the ends illustrative), offset 8 µm, deorientation from its Table I, the
    particle–γ coincidences with the ²⁴Mg recoils (the stopper takes the scattered ⁹⁸Ru), the rings at 41.5° and
    138.5°. Analysed as the paper did: the speed from the shifted peak's centroids, the peak areas per ring, and a
    decay-curve fit per ring of the unshifted intensity per particle (a free normalisation, as their
    P(d) = a(1 − e^(−(d−d₀)/vτ)) of the shifted one; it takes up the efficiency and the Lorentz boost of each
    component), corrected for the deorientation: divided by its angular distribution W_U(d) = 1 + Σ a_k G_k(t)
    P_k(cos θ), G_k from Table I at the flight time and a_k from physim's Coulex correlation for the particle
    detector (:func:`coulex_alignment`). τ is the rings' weighted mean, compared with 8.36(29) ps within 2σ of the
    combined uncertainty. The 4⁺₁ and 2⁺₂ need two-step Coulex with matrix elements the paper does not give: not
    reproduced."""
    from .experiment import Experiment

    pub = published(ref)
    if pub is None:
        return {"available": False, "reason": f"{PUBLISHED} not found"}
    key = ("ru98", particles, seed, gammas_per_event, str(_reference(ref)))
    if key in _CACHE:
        return _CACHE[key]
    paper = _paper(pub, "ru98_rdds_coulex", "2+1")
    exp = Experiment.example("coulex_ru98_plunger")
    distances = [_quantity(d) for d in exp.lifetime.distances]
    offset = _quantity(exp.lifetime.offset) if exp.lifetime.offset is not None else RU98_OFFSET_UM
    dor = exp.lifetime.deorientation
    alpha, gamma_ps = dor.alpha, dor.gamma_per_ps
    a_k = coulex_alignment(exp)
    run = coulex_plunger(exp, distances, particles=particles, seed=seed, gammas_per_event=gammas_per_event)
    beta = beta_from_centroids(run)
    t = lt.flight_times(np.array(distances) + offset, beta, proper=True)
    g_u = lt.deorientation(t, alpha, gamma_ps)
    rings = {}
    for r in RU98_RINGS:
        a = plunger_areas(run, r, beta)
        u, su = a["unshifted"], a["err_unshifted"]
        plain = la.decay_curve(distances, u, su, beta, offset_um=offset, free_norm=True)
        c = math.cos(math.radians(r))
        w_u = 1 + a_k[2] * np.asarray(g_u[2]) * _p2(c) + a_k[4] * np.asarray(g_u[4]) * _p4(c)
        corrected = la.decay_curve(distances, u / w_u, su / w_u, beta, offset_um=offset, free_norm=True)
        rings[r] = {"tau_ps": corrected.tau_ps, "stat_ps": corrected.stat_ps, "uncorrected_ps": plain.tau_ps,
                    "uncorrected_stat_ps": plain.stat_ps, "w_unshifted": w_u.tolist(), "unshifted": u.tolist(),
                    "err_unshifted": su.tolist(), "shifted": a["shifted"].tolist(),
                    "model": [p["model"] for p in corrected.per_distance] if corrected.per_distance else None}
    w = np.array([1 / rings[r]["stat_ps"] ** 2 for r in RU98_RINGS])
    tau = float(np.sum(w * [rings[r]["tau_ps"] for r in RU98_RINGS]) / w.sum())
    stat = float(1 / math.sqrt(w.sum()))
    wp = np.array([1 / rings[r]["uncorrected_stat_ps"] ** 2 for r in RU98_RINGS])
    tau_plain = float(np.sum(wp * [rings[r]["uncorrected_ps"] for r in RU98_RINGS]) / wp.sum())
    dev = tau - paper["tau_ps"]
    allowed = 2 * math.hypot(_sigma(paper, dev), stat)
    out = {"available": True, "tau_ps": tau, "stat_ps": stat, "uncorrected_ps": tau_plain, "beta": beta,
           "paper_tau_ps": paper["tau_ps"], "paper_sigma_ps": _sigma(paper, dev), "deviation_ps": dev,
           "allowed_ps": allowed, "worst": abs(dev) / allowed, "passed": abs(dev) <= allowed, "rings": rings,
           "a_k": {2: a_k[2], 4: a_k[4]}, "distances_um": distances, "offset_um": offset,
           "true_unshifted": run["true_unshifted"], "gammas": run["gammas"], "particles": run["particles"],
           "simulated_tau_ps": float(list(lt.level_taus(_excited_scheme(exp)).values())[0][0])}
    out["detail"] = (f"τ(2⁺₁) = {tau:.2f} ± {stat:.2f} ps with the deorientation correction (rings "
                     + ", ".join(f"{r:g}°: {rings[r]['tau_ps']:.2f} ± {rings[r]['stat_ps']:.2f}" for r in RU98_RINGS)
                     + f"), against the paper's {paper['tau_ps']:.2f}({100 * paper['plus_ps']:.0f}) ps: "
                     f"{dev:+.2f} ps, {abs(dev) / allowed:.0%} of the 2σ allowed; without the correction "
                     f"{tau_plain:.2f} ps. β from the shifted centroids {100 * beta:.2f}% (paper 4.50(6)%); "
                     f"{sum(run['gammas'])} γ rays with ²⁴Mg recoils over {len(distances)} distances; a₂ = "
                     f"{a_k[2]:.3f}, a₄ = {a_k[4]:.3f} from physim's Coulex correlation. The 4⁺₁ (2.31 ps) and 2⁺₂ "
                     "(1.7 ps) are not reproduced: their two-step excitation needs matrix elements the paper does "
                     "not give")
    _CACHE[key] = out
    return out


def _excited_scheme(experiment):
    from . import gamma

    return gamma.excitation_of(experiment).scheme


# -- the simulated series NaPyTau analysed (and the one-level closure) ---------------------------------------------


def napytau_input(*, particles: int = 4_000_000, seed: int = 11, gammas_per_event: int = 4,
                  ring: float = 41.5) -> dict:
    """The simulated plunger series for the NaPyTau comparison and the one-level closure: the ⁹⁸Ru example without
    its deorientation (one exponential), the forward ring, peak areas at every distance, the speed from the
    shifted centroids. Written to ``napytau/ddcm_simulated.csv`` by ``scripts/lifetime_reference.py``; it takes about
    a minute at the default size."""
    from .experiment import Experiment

    exp = Experiment.example("coulex_ru98_plunger")
    exp.lifetime.deorientation = None
    distances = [_quantity(d) for d in exp.lifetime.distances]
    offset = _quantity(exp.lifetime.offset) if exp.lifetime.offset is not None else RU98_OFFSET_UM
    run = coulex_plunger(exp, distances, particles=particles, seed=seed, gammas_per_event=gammas_per_event)
    beta = beta_from_centroids(run)
    a = plunger_areas(run, ring, beta)
    tau = float(list(lt.level_taus(_excited_scheme(exp)).values())[0][0])
    return {"distances_um": distances, "offset_um": offset, "beta": beta, "tau_ps": tau, "ring_deg": ring,
            "particles": particles, "seed": seed, "gammas_per_event": gammas_per_event, **a}


def read_napytau_simulated(ref=None) -> Optional[dict]:
    """``napytau/ddcm_simulated.csv``: the simulated areas, the speed, offset and true τ, and NaPyTau's τ(x) for each
    polynomial degree; None if the file is not there."""
    from .validation import read_reference

    base = _reference(ref)
    if base is None or not (base / NAPYTAU_SIMULATED).exists():
        return None
    header, rows = read_reference(base / NAPYTAU_SIMULATED)
    col = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    out = {"distances_um": col["distance_um"], "unshifted": col["unshifted"], "shifted": col["shifted"],
           "err_unshifted": col["err_unshifted"], "err_shifted": col["err_shifted"],
           "napytau": {d: col[f"napytau_tau_deg{d}_ps"] for d in NAPYTAU_DEGREES},
           "beta": float(header["beta"].split()[0]), "offset_um": float(header["offset_um"].split()[0]),
           "tau_ps": float(header["tau_ps"].split()[0])}
    return out


def ddcm_vs_napytau(ref=None) -> Optional[dict]:
    """physim's DDCM in NaPyTau's procedure (one polynomial in t through the shifted points, unweighted, τ(x) =
    unshifted / P′(t), no time dilation; :func:`lifetime_analysis.ddcm` with ``fit="shifted"``) on the simulated
    areas, against NaPyTau's τ(x) for degrees 2–4 at every distance (both on the intensities normalised by
    I_u + I_s, NaPyTau's distances with the offset added). Also physim's own DDCM (joint fit) for comparison."""
    sim = read_napytau_simulated(ref)
    if sim is None:
        return None
    d, u, s = sim["distances_um"], sim["unshifted"], sim["shifted"]
    errs = (sim["err_unshifted"], sim["err_shifted"])
    per = {}
    worst = 0.0
    for deg in NAPYTAU_DEGREES:
        r = la.ddcm(d, u, s, errs, sim["beta"], polynomial=deg, fit="shifted", weighted=False, proper=False,
                    offset_um=sim["offset_um"])
        ours = np.array([p["tau_ps"] for p in r.per_distance])
        theirs = sim["napytau"][deg]
        ok = np.isfinite(ours) & np.isfinite(theirs) & (theirs > 0)
        dev = np.abs(ours[ok] / theirs[ok] - 1)
        worst = max(worst, float(dev.max()))
        per[deg] = {"physim": ours.tolist(), "napytau": theirs.tolist(), "worst": float(dev.max()),
                    "points": int(ok.sum())}
    joint = la.ddcm(d, u, s, errs, sim["beta"], offset_um=sim["offset_um"])
    return {"per_degree": per, "worst": worst, "joint_tau_ps": joint.tau_ps, "joint_stat_ps": joint.stat_ps,
            "tau_ps": sim["tau_ps"], "distances_um": d.tolist()}


def closure_one_level(ref=None) -> Optional[dict]:
    """τ back from the simulated series of :func:`napytau_input` (one level, known τ), with time dilation: the
    decay-curve fit of the unshifted intensity with a free normalisation (the same number of reactions at every
    distance), and the DDCM (joint fit) with the shifted intensities divided by the efficiency ratio
    (:func:`efficiency_ratio`), against the simulation's τ. (A fit of the unshifted share would need the decays in
    the stopper while the recoil slows, counted in the shifted window, put back; these two do not.)"""
    sim = read_napytau_simulated(ref)
    if sim is None:
        return None
    d, u, s = sim["distances_um"], sim["unshifted"], sim["shifted"]
    su, ss = sim["err_unshifted"], sim["err_shifted"]
    k = efficiency_ratio(u, s, su, ss)
    dc = la.decay_curve(d, u, su, sim["beta"], offset_um=sim["offset_um"], free_norm=True)
    dd = la.ddcm(d, u, s / k, (su, ss / k), sim["beta"], offset_um=sim["offset_um"])
    return {"tau_ps": sim["tau_ps"], "decay_curve": (dc.tau_ps, dc.stat_ps), "ddcm": (dd.tau_ps, dd.stat_ps),
            "efficiency_ratio": k, "distances_um": d.tolist(), "unshifted": u.tolist(),
            "model": [p["model"] for p in dc.per_distance]}


# ---------------------------------------------------------------------------------------------------------------
# ⁵⁰Cr: three-foil RDDS after fusion-evaporation (Beckers et al. 2022)

#: The pairs (target–degrader x₁, degrader–stopper x₂) in µm. The paper names 9 sets with x₁ from 17 to 750 µm and
#: x₂ from 40 to 335 µm and gives (42, 40), (62, 67) and (122, 127); the others are illustrative, spread over those
#: ranges as its Figs. 8 and 9 show.
CR50_PAIRS = ((17, 40), (42, 40), (62, 67), (90, 80), (122, 127), (200, 150), (350, 200), (500, 250), (750, 335))


def fusion_decays(experiment, residues: int, seed: int) -> dict:
    """The studied residue's decays: its starting states (:func:`physim.nuclear.residues.simulate_residues`, the
    residues of ``residues`` reactions) and its cascade with each γ ray's proper time after formation
    (:meth:`physim.nuclear.residues.ResidueCascade.run`: side-feeding band and level lifetimes)."""
    from . import residues as res

    cas = res.lifetime_cascade(experiment)
    start = res._sample(experiment, cas["index"], int(residues), int(seed))
    n = len(start["energy"])
    out = cas["cascade"].run(np.full(n, float(cas["entry_spin"])), np.random.default_rng([int(seed), 3]),
                             timing_rng=np.random.default_rng([int(seed), 4]))
    return {"start": start, "gammas": {k: np.asarray(v) for k, v in out.items()}, "residue": cas["residue"],
            "names": cas["names"], "cascade": cas["cascade"], "n": n}


def _pair_stack(experiment, x1: float, x2: float) -> list:
    from .rates import recoil_stack

    segs = recoil_stack(experiment, float(x1))
    return [dataclasses.replace(s, gap_um=float(x2)) if s.name == "gap2" else s for s in segs]


def three_foil_components(experiment, decays: dict, pair) -> dict:
    """Where each γ ray of ``decays`` is emitted for one pair of distances: 0 fully shifted (target and first gap),
    1 reduced shifted (in the degrader and the second gap: after the recoil reached the degrader, before the
    stopper), 2 unshifted (in the stopper); and its speed then (straight paths, physim's stopping)."""
    from . import transport

    segs = _pair_stack(experiment, *pair)
    st, g = decays["start"], decays["gammas"]
    src = g["source"].astype(int)
    species = decays["residue"]
    mass = data.nuclide(species).nuclear_mass_mev
    lab = transport.lab_time(segs, species, mass, st["energy"][src], st["direction"][src], st["depth"][src],
                             g["time_ps"])
    s = transport.slow(segs, species, mass, st["energy"][src], st["direction"][src], st["depth"][src], lab)
    names = [x.name for x in segs]
    deg, sto = names.index("degrader"), names.index("degrader2")
    seg = np.asarray(s.segment)
    comp = np.where(seg >= sto, 2, np.where(seg >= deg, 1, 0))
    return {"component": comp, "beta": np.asarray(s.beta), "segment": seg, "gaps": (names.index("gap"),
                                                                                 names.index("gap2"))}


def _counts(mask_a, comp_a) -> np.ndarray:
    return np.array([np.sum(mask_a & (comp_a == k)) for k in range(3)], dtype=float)


def cr50(ref=None, *, residues: int = 400_000, seed: int = 1) -> dict:
    """Beckers et al.'s ⁵⁰Cr 2⁺₁ and 4⁺₁ lifetimes from a simulated three-foil plunger run of
    ``fusion_cr50_plunger``: ⁵⁰Cr residues from physim's fusion events, their cascade timed by its lifetimes and
    side feeding, flown and slowed through target, gap, degrader, gap and stopper at each pair of
    :data:`CR50_PAIRS`. The γγ coincidences are counted where the γ rays are emitted (ideal detectors and complete
    gates; as many as the residues simulated, about the paper's 10⁴–10⁵ per distance, its Fig. 8). Analysed as the
    paper did: the speeds before and after the degrader from the shifted components, γγ DDCM with direct gates
    (:func:`lifetime_analysis.ddcm_three_foil`, its Eq. 13), for the 2⁺₁ the degrader ("t₁") part gated on the 4⁺
    feeder's fully shifted component and the stopper ("t₂") part on its fully and reduced shifted ones, each
    analysed apart and averaged; for the 4⁺₁ the t₁ part gated on the 6⁺ feeder. The speeds are those of the
    decays in the two gaps (the shifted and reduced-shifted peaks' centroids, as the paper measured them).
    Compared with the paper's values and the other methods' it quotes (Table 1): within the 2σ spread of the
    published values. Also returns the singles analysis with the 4⁺ feeder's correction, each part apart and
    averaged (the closure "fed from above"). The parts are kept apart because the paper's x_eff = x₁ + x₂ v/v′
    leaves out the time the recoil takes through the degrader, which shifts the stopper points against the
    degrader points (the paper's remark on absolute distances); each part on its own is not affected (a constant
    shift of the time axis does not change a DDCM)."""
    from .experiment import Experiment

    pub = published(ref)
    if pub is None:
        return {"available": False, "reason": f"{PUBLISHED} not found"}
    key = ("cr50", residues, seed, str(_reference(ref)))
    if key in _CACHE:
        return _CACHE[key]
    exp = Experiment.example("fusion_cr50_plunger")
    dec = fusion_decays(exp, residues, seed)
    names = dec["names"]
    g = dec["gammas"]
    lev = {nm: i for i, nm in enumerate(names)}
    i2, i4, i6 = lev["2+"], lev["4+"], lev["6+"]
    src = g["source"].astype(int)
    n = dec["n"]

    def per_residue(initial, final, values):
        out = np.full(n, -1)
        m = (g["initial"] == initial) & (g["final"] == final)
        out[src[m]] = values[m]
        return out

    pairs = np.array(CR50_PAIRS, dtype=float)
    g1_2, g2_2, g1_4, sing_a, sing_b = [], [], [], [], []
    beta_fs, beta_rs = [], []
    for pair in pairs:
        c = three_foil_components(exp, dec, pair)
        comp = c["component"]
        a = per_residue(i2, 0, comp)  # 2⁺ → 0⁺
        b = per_residue(i4, i2, comp)  # 4⁺ → 2⁺
        cc = per_residue(i6, i4, comp)  # 6⁺ → 4⁺
        g1_2.append(_counts(b == 0, a))
        g2_2.append(_counts((b == 0) | (b == 1), a))
        g1_4.append(_counts(cc == 0, b))
        sing_a.append(_counts(a >= 0, a))
        sing_b.append(_counts(b >= 0, b))
        m_a = (g["initial"] == i2) & (g["final"] == 0)
        gap1, gap2 = c["gaps"]
        beta_fs.append(np.mean(c["beta"][m_a & (c["segment"] == gap1)]))
        beta_rs.append(np.mean(c["beta"][m_a & (c["segment"] == gap2)]))
    g1_2, g2_2, g1_4 = np.array(g1_2), np.array(g2_2), np.array(g1_4)
    sing_a, sing_b = np.array(sing_a), np.array(sing_b)
    beta, beta_red = float(np.mean(beta_fs)), float(np.mean(beta_rs))

    def errs(x):
        return tuple(np.sqrt(np.maximum(x[:, k], 1.0)) for k in range(3))

    def cols(x):
        return x[:, 0], x[:, 1], x[:, 2]

    t1 = la.ddcm_three_foil(pairs, *cols(g1_2), errs(g1_2), beta, beta_red, gated=True, points="degrader")
    t2 = la.ddcm_three_foil(pairs, *cols(g2_2), errs(g2_2), beta, beta_red, gated=True, points="stopper")
    w = np.array([1 / t1.stat_ps**2, 1 / t2.stat_ps**2])
    tau2 = float((w[0] * t1.tau_ps + w[1] * t2.tau_ps) / w.sum())
    stat2 = float(1 / math.sqrt(w.sum()))
    t4 = la.ddcm_three_foil(pairs, *cols(g1_4), errs(g1_4), beta, beta_red, gated=True, points="degrader")
    # Singles with the 4⁺ feeder's correction (all its decays go to the 2⁺): the closure fed from above.
    single_parts = {part: la.ddcm_three_foil(pairs, *cols(sing_a), errs(sing_a) + errs(sing_b), beta, beta_red,
                                             *cols(sing_b), alpha=1.0, points=part)
                    for part in ("degrader", "stopper")}
    ws = np.array([1 / r.stat_ps**2 for r in single_parts.values()])
    singles_tau = float(np.sum(ws * [r.tau_ps for r in single_parts.values()]) / ws.sum())
    taus = dec["cascade"].taus_ps
    states = {}
    worst = 0.0
    passed = True
    for state, (tau, stat) in (("2+1", (tau2, stat2)), ("4+1", (t4.tau_ps, t4.stat_ps))):
        rows = [r for r in pub["lifetimes"] if r["case"] == "cr50_rdds_fusion" and r["state"] == state]
        lo = min(r["tau_ps"] - 2 * r["minus_ps"] for r in rows)
        hi = max(r["tau_ps"] + 2 * r["plus_ps"] for r in rows)
        paper = _paper(pub, "cr50_rdds_fusion", state)
        inside = lo <= tau <= hi
        half = (hi - lo) / 2
        worst = max(worst, abs(tau - (lo + hi) / 2) / half)
        passed &= inside
        dev = tau - paper["tau_ps"]
        states[state] = {"tau_ps": tau, "stat_ps": stat, "paper_tau_ps": paper["tau_ps"],
                         "paper_sigma_ps": _sigma(paper, dev), "band_ps": [lo, hi], "inside": inside,
                         "deviation_sigma": dev / math.hypot(_sigma(paper, dev), stat),
                         "published": [(r["method"], r["tau_ps"], r["plus_ps"], r["source"]) for r in rows]}
    states["2+1"]["parts"] = {"t1": (t1.tau_ps, t1.stat_ps), "t2": (t2.tau_ps, t2.stat_ps)}
    out = {"available": True, "states": states, "worst": worst, "passed": passed, "beta": beta,
           "beta_reduced": beta_red, "pairs_um": pairs.tolist(), "residues": n,
           "simulated_tau_ps": {"2+1": float(taus[i2]), "4+1": float(taus[i4]), "6+1": float(taus[i6])},
           "side_feeding_ps": {"2+1": float(dec["cascade"].side_feeding_ps[i2]),
                               "4+1": float(dec["cascade"].side_feeding_ps[i4])},
           "singles": (singles_tau, float(1 / math.sqrt(ws.sum()))),
           "singles_parts": {k: (r.tau_ps, r.stat_ps) for k, r in single_parts.items()},
           "gated_2": {"t1": g1_2.tolist(), "t2": g2_2.tolist()}, "gated_4": g1_4.tolist(),
           "per_distance": {"t1": t1.per_distance, "t2": t2.per_distance, "4+1": t4.per_distance}}
    s2, s4 = states["2+1"], states["4+1"]
    out["detail"] = (f"τ(2⁺₁) = {s2['tau_ps']:.2f} ± {s2['stat_ps']:.2f} ps (t₁ part {t1.tau_ps:.2f} ± "
                     f"{t1.stat_ps:.2f}, t₂ part {t2.tau_ps:.2f} ± {t2.stat_ps:.2f}) against the paper's "
                     f"{s2['paper_tau_ps']:.1f}(4) ps, the published values' 2σ band {s2['band_ps'][0]:.1f}–"
                     f"{s2['band_ps'][1]:.1f} ps: {'inside' if s2['inside'] else 'outside'}; τ(4⁺₁) = "
                     f"{s4['tau_ps']:.2f} ± {s4['stat_ps']:.2f} ps (t₁ part, gated on the 6⁺ feeder) against "
                     f"{s4['paper_tau_ps']:.1f}(3) ps, band {s4['band_ps'][0]:.1f}–{s4['band_ps'][1]:.1f} ps: "
                     f"{'inside' if s4['inside'] else 'outside'}. β {100 * beta:.2f}% then {100 * beta_red:.2f}% "
                     f"(paper 3.9% and 1.5%); {n} ⁵⁰Cr residues, γγ counted at emission (ideal detectors, complete "
                     "gates)")
    _CACHE[key] = out
    return out


# ---------------------------------------------------------------------------------------------------------------
# ¹⁴³Sm: DSAM after fusion-evaporation on a Au backing (Rajbanshi et al. 2018)

#: Relative γ-ray intensities of the band members' in-band transitions (all of each level's decays, %), which set
#: the side feeding of each level by intensity balance. The paper does not give them: illustrative, falling with
#: spin as in its spectra, DB II weaker than DB I.
SM143_INTENSITY = {"45/2- (DB I)": 100.0, "47/2- (DB I)": 80.0, "49/2- (DB I)": 65.0, "51/2- (DB I)": 30.0,
                   "53/2- (DB I)": 12.0, "51/2- (DB II)": 25.0, "53/2- (DB II)": 15.0, "55/2(-) (DB II)": 8.0,
                   "57/2(-) (DB II)": 4.0}
#: The angles the paper fitted with a Doppler shift (it fitted 90° too, which has none).
SM143_ANGLES = (148.0, 123.0)
#: Counts in the line of the strongest transition (239.5 keV) at each angle; the others scale with
#: :data:`SM143_INTENSITY`. Illustrative: about what the paper's gated spectra (its Fig. 2) hold.
SM143_COUNTS = 6000.0
#: The steps of each side-feeding cascade, each with its share of the paper's side-feeding time. LINESHAPE, which the
#: paper used, takes five; one level with the paper's (effective) time is used here, in the simulation and the fit
#: alike, because physim's Bateman solution (a matrix exponential at every kernel time) for the 56 levels of five
#: steps takes 0.5 s per line shape against 0.03 s for 20 levels, too slow for the validation suite.
SM143_SIDE_STEPS = 1
#: The lines of the closure (the paper's Fig. 2: 444.8 and 527.5 keV).
SM143_CLOSURE = ("49/2- (DB I)", "51/2- (DB II)")
#: SRIM's stopping powers are good to about this (the paper's systematic, not in its errors).
SM143_STOPPING = 0.05
#: The detectors' FWHM (keV) for the line shapes (the example's clovers).
SM143_FWHM_KEV = 2.1


def _sm143_labels(scheme) -> list:
    """The paper's names of the band levels (``"51/2- (DB I)"``, as in lifetimes.csv) for the scheme's levels, None
    for the others."""
    out = []
    seen: dict = {}
    for lv in scheme.levels:
        jpi = lv.jpi if hasattr(lv, "jpi") else None
        out.append(jpi)
    labels = []
    for i, jpi in enumerate(out):
        if jpi in ("51/2-", "53/2-"):
            k = seen.get(jpi, 0)
            seen[jpi] = k + 1
            # the lower of the two is DB II's (Fig. 1: 10165 and 10750 keV against DB I's 10216 and 10944 keV)
            labels.append(f"{jpi} ({'DB II' if k == 0 else 'DB I'})")
        elif jpi in ("55/2-", "57/2-"):
            labels.append(jpi.replace("-", "(-)") + " (DB II)")
        elif jpi in ("45/2-", "47/2-", "49/2-"):
            labels.append(jpi + " (DB I)")
        else:
            labels.append(None)
    return labels


def sm143_cascade(scheme, side_feeding_ps: dict, steps: int = None) -> dict:
    """The cascade for the line shapes, as LINESHAPE models it: the band levels with their mean lives and
    branchings from the setup's level scheme, each fed by a side-feeding cascade of ``steps`` (default
    :data:`SM143_SIDE_STEPS`) steps with the paper's side-feeding time (Table I) shared evenly, the side feeding of each level its
    intensity (:data:`SM143_INTENSITY`) less what the band brings it. Returns ``taus_ps``, ``feeding``,
    ``populations`` (scheme levels, then the side steps), ``labels``, ``index`` {label: level}, ``e_gamma``
    {label: in-band γ-ray energy keV}, ``branch`` (each level's share of decays per final level)."""
    from .gamma import level_names

    labels = _sm143_labels(scheme)
    index = {lab: i for i, lab in enumerate(labels) if lab is not None}
    taus = list(lt_level_taus_array(scheme))
    n = len(scheme.levels)
    branch = np.zeros((n, n))
    e_gamma = {}
    for i in range(n):
        out = [t for t in scheme.transitions if t.initial == i]
        tot = sum(float(t.intensity.value) for t in out if t.intensity is not None)
        for t in out:
            branch[i, t.final] += float(t.intensity.value) / tot
        if labels[i] is not None and out:
            # the in-band (strongest) transition
            e_gamma[labels[i]] = float(max(out, key=lambda t: float(t.intensity.value)).energy.value)
    intensity = np.zeros(n)
    for lab, i in index.items():
        intensity[i] = SM143_INTENSITY[lab]
    side = intensity - branch.T @ intensity  # out less what comes in from the band
    if np.any(side[list(index.values())] < -1e-9):
        raise ValueError("the band's intensities leave a level with negative side feeding")
    side = np.clip(side, 0.0, None)
    steps = int(SM143_SIDE_STEPS if steps is None else steps)
    m = n + steps * len(index)
    feeding = np.zeros((m, m))
    feeding[:n, :n] = branch
    pops = np.zeros(m)
    taus_all = taus + [0.0] * (steps * len(index))
    k = n
    for lab, i in index.items():
        tau_sf = float(side_feeding_ps[lab])
        below = i
        for _ in range(steps):
            taus_all[k] = tau_sf / steps
            feeding[k, below] = 1.0
            below = k
            k += 1
        pops[below] = side[i]
    return {"taus_ps": np.array(taus_all), "feeding": feeding, "populations": pops / pops.sum(), "labels": labels,
            "index": index, "e_gamma": e_gamma, "names": level_names(scheme.levels), "intensity": intensity}


def _feeders(feeding: np.ndarray, level: int) -> np.ndarray:
    """The levels whose decays reach ``level`` (and the level itself), in order: the sub-cascade a line shape of
    that level needs."""
    reach = np.linalg.inv(np.eye(len(feeding)) - feeding)  # (I − F)⁻¹: paths of any length, row = from
    return np.flatnonzero((reach[:, level] > 1e-12) | (np.arange(len(feeding)) == level))


def lt_level_taus_array(scheme) -> np.ndarray:
    """Each level's mean life (ps) from its half-life or B values (:func:`lifetimes.level_taus`); inf for a level
    that does not decay within the scheme."""
    from .gamma import level_names

    found = lt.level_taus(scheme)
    return np.array([found[nm][0] if nm in found else np.inf for nm in level_names(scheme.levels)])


@contextmanager
def _stopping_in(species: str, material, st):
    """Within the block, the transport slows ``species`` in ``material`` (only) with the Stopping ``st``."""
    from . import transport

    saved, clocks = transport.stopping, dict(transport._CLOCKS)

    def pick(ion, mat):
        if ion == species and getattr(mat, "name", None) == getattr(material, "name", object()):
            return st
        return saved(ion, mat)

    transport.stopping = pick
    transport._CLOCKS.clear()
    try:
        yield
    finally:
        transport.stopping = saved
        transport._CLOCKS.clear()
        transport._CLOCKS.update(clocks)


def dsam_kernel(experiment, *, residues: int = 20_000, n: int = 3000, times: int = 80, seed: int = 1,
                kernel_seed: int = 2, stopping=None):
    """The DSAM kernel of the setup's studied residue: its starting states from physim's fusion events, slowed
    and scattered through target and backing (:func:`physim.nuclear.transport.kernel`, ``n`` samples at ``times``
    lab times), with physim's stopping or, with ``stopping`` = (material, Stopping), that Stopping for the residue in
    that material."""
    from . import residues as res
    from . import transport
    from .rates import recoil_stack

    table = res.ChannelTable.from_experiment(experiment)
    idx = res._studied(experiment, table, None)
    species = table.channels[idx].residue
    mass = data.nuclide(species).nuclear_mass_mev
    start = res._sample(experiment, idx, int(residues), int(seed))
    segs = recoil_stack(experiment)

    def build():
        t_stop = transport.time_to_stop(segs, species, mass, start["energy"], start["depth"])
        t_max = 1.5 * float(np.max(t_stop[np.isfinite(t_stop)]))
        grid = np.r_[0.0, np.geomspace(1e-3 * t_max, t_max, int(times) - 1)]
        return transport.kernel(segs, species, mass, start["energy"], start["direction"], start["depth"], grid,
                                scatter=True, rng=np.random.default_rng([int(seed), int(kernel_seed)]), n=int(n))

    if stopping is None:
        return build()
    material, st = stopping
    with _stopping_in(species, material, st):
        return build()


def _fold_uniform(edges: np.ndarray, shape: np.ndarray, fwhm_kev: float) -> np.ndarray:
    from .validation import _fold

    return _fold(edges, shape, fwhm_kev)


def _rebin(fine_edges: np.ndarray, contents: np.ndarray, edges: np.ndarray) -> np.ndarray:
    from .validation import _rebin as rebin

    return rebin(fine_edges, contents, edges)


def sm143(ref=None, *, residues: int = 20_000, n: int = 3000, times: int = 60, seed: int = 1,
          side_steps: int = None) -> dict:
    """Rajbanshi et al.'s ¹⁴³Sm lifetimes from simulated DSAM line shapes of ``fusion_sm143_dsam``: the residues
    from physim's fusion events slow and scatter through the 0.8 mg/cm² ¹²⁴Sn target and the 13 mg/cm² Au backing
    with physim's stopping (the "true" kernel), and the band's cascade (:func:`sm143_cascade`: the paper's
    lifetimes, branchings and side-feeding times) gives each line's shape at 148° and 123°, folded with the
    detectors' resolution, with Poisson noise (:data:`SM143_COUNTS` in the strongest line; illustrative) on a flat
    background. Analysed as the paper did: LINESHAPE's way, a kernel from a separate set of recoils slowed with
    SRIM's stopping in the Au (lifetimes/trim/srim_stopping.csv), fits from the top of each band down
    (:func:`lifetime_analysis.fit_lineshape`; the feeders' τ fixed at the values just fitted, the side feeding at
    the paper's times: gates from below), each angle apart and their weighted mean. Compared with Table I: within
    2σ of the paper's uncertainty, the fit's and SRIM's ±5% stopping added in quadrature; for the two upper limits,
    the fitted τ less 2σ below the limit. Also the closure: the same spectra fitted with a kernel from a separate
    set of recoils slowed with physim's own stopping."""
    from .experiment import Experiment
    from .validation import _lifetime_file, srim_stopping

    pub = published(ref)
    base = _reference(ref)
    srim = None if base is None else _lifetime_file(base, "trim/srim_stopping.csv")
    if pub is None or srim is None:
        return {"available": False, "reason": f"{PUBLISHED} or lifetimes/trim/srim_stopping.csv not found"}
    key = ("sm143", residues, n, times, seed, side_steps, str(base))
    if key in _CACHE:
        return _CACHE[key]
    exp = Experiment.example("fusion_sm143_dsam")
    scheme = exp.levels["143Sm"]
    sf = pub["setups"]["sm143_dsam_fusion"]["side_feeding_ps"]
    side = {}
    for lab in SM143_INTENSITY:
        jpi = lab.split(" ")[0]
        side[lab] = float(sf.get(lab, sf.get(jpi)))
    cas = sm143_cascade(scheme, side, side_steps)
    gold = data.material("Au")
    k_true = dsam_kernel(exp, residues=residues, n=n, times=times, seed=seed, kernel_seed=2)
    k_srim = dsam_kernel(exp, residues=residues, n=n, times=times, seed=seed, kernel_seed=3,
                         stopping=(gold, srim_stopping(base, "143Sm", gold)))
    k_own = dsam_kernel(exp, residues=residues, n=n, times=times, seed=seed, kernel_seed=4)
    taus, feed, pops = cas["taus_ps"], cas["feeding"], cas["populations"]
    rng = np.random.default_rng([int(seed), 143])
    beta0 = float(np.max(k_true.beta0))
    spectra = {}
    for lab, i in cas["index"].items():
        e0 = cas["e_gamma"][lab]
        half = e0 * 1.3 * beta0 + 4 * SM143_FWHM_KEV
        edges = np.arange(e0 - half, e0 + half + 1e-9, 0.5)
        spectra[lab] = {"edges": edges, "counts": {}}
        # The shape on a fine grid that has no edge at E₀ (the stopped peak is a spike there), folded with the
        # resolution, then put in the 0.5 keV bins.
        fine = np.arange(edges[0] - 10.0, edges[-1] + 10.0, 0.02) + 0.0071
        for th in SM143_ANGLES:
            row = lt.lineshape(e0, taus, feed, pops, k_true, th, fine)[i]
            shape = _rebin(fine, _fold_uniform(fine, np.nan_to_num(row), SM143_FWHM_KEV), edges)
            amp = SM143_COUNTS * cas["intensity"][i] / 100.0
            mu = amp * shape / shape.sum() + 0.02 * amp * shape.max() / shape.sum()
            spectra[lab]["counts"][th] = rng.poisson(mu).astype(float)
    order = sorted(cas["index"], key=lambda lab: -scheme.levels[cas["index"][lab]].energy.value)

    def fit_all(kernel, labels=None) -> dict:
        """Fits from the top down, each level's feeders at the values just fitted; with ``labels``, only those
        levels, their feeders at the simulated values (the closure)."""
        fitted = dict(enumerate(taus))
        out = {}
        for lab in (order if labels is None else labels):
            i = cas["index"][lab]
            per = {}
            sub = _feeders(feed, i)
            for th in SM143_ANGLES:
                f = la.Feeding([fitted[j] for j in sub], feed[np.ix_(sub, sub)], pops[sub],
                               int(np.flatnonzero(sub == i)[0]))
                r = la.fit_lineshape(spectra[lab]["counts"][th], spectra[lab]["edges"], cas["e_gamma"][lab], kernel,
                                     th, f, "linear", fwhm_kev=SM143_FWHM_KEV, stopping_uncertainty=0.0)
                per[th] = (r.tau_ps, r.stat_ps, r.params.get("plus_ps"), r.params.get("minus_ps"))
            w = np.array([1 / per[th][1] ** 2 for th in SM143_ANGLES])
            tau = float(np.sum(w * [per[th][0] for th in SM143_ANGLES]) / w.sum())
            out[lab] = {"tau_ps": tau, "stat_ps": float(1 / math.sqrt(w.sum())), "angles": per}
            if labels is None:
                fitted[i] = tau
        return out

    srim_fits = fit_all(k_srim)
    own_fits = fit_all(k_own, list(SM143_CLOSURE))
    levels = {}
    worst = 0.0
    passed = True
    for lab in order:
        row = _paper(pub, "sm143_dsam_fusion", lab)
        fit = srim_fits[lab]
        tau, stat = fit["tau_ps"], fit["stat_ps"]
        if row["kind"] == "upper":
            allowed = 2 * stat
            dev = (tau - 2 * stat) - row["tau_ps"]
            ratio = max(dev, 0.0) / (row["tau_ps"]) if dev > 0 else 0.0
            ok = dev <= 0
            levels[lab] = {"tau_ps": tau, "stat_ps": stat, "paper_tau_ps": row["tau_ps"], "kind": "upper",
                           "passed": ok, "worst": 0.0 if ok else 1.0 + ratio}
        else:
            dev = tau - row["tau_ps"]
            allowed = 2 * math.sqrt(_sigma(row, dev) ** 2 + stat**2 + (SM143_STOPPING * row["tau_ps"]) ** 2)
            ok = abs(dev) <= allowed
            levels[lab] = {"tau_ps": tau, "stat_ps": stat, "paper_tau_ps": row["tau_ps"], "kind": "value",
                           "paper_sigma_ps": _sigma(row, dev), "allowed_ps": allowed, "passed": ok,
                           "worst": abs(dev) / allowed}
        levels[lab].update(e_gamma_kev=cas["e_gamma"][lab], side_feeding_ps=side[lab],
                           simulated_tau_ps=float(taus[cas["index"][lab]]), angles=fit["angles"])
        if lab in own_fits:
            levels[lab].update(own_tau_ps=own_fits[lab]["tau_ps"], own_stat_ps=own_fits[lab]["stat_ps"])
        worst = max(worst, levels[lab]["worst"])
        passed &= ok
    out = {"available": True, "levels": levels, "worst": worst, "passed": passed, "order": order,
           "beta0": float(np.mean(k_true.beta0)), "spectra": spectra, "angles": list(SM143_ANGLES)}
    vals = [lab for lab in order if levels[lab]["kind"] == "value"]
    ratios = np.array([levels[lab]["tau_ps"] / levels[lab]["simulated_tau_ps"] for lab in order])
    out["srim_over_simulated"] = [float(ratios.min()), float(ratios.max()), float(np.mean(ratios))]
    out["detail"] = ("; ".join(f"{lab} ({levels[lab]['e_gamma_kev']:.1f} keV) {levels[lab]['tau_ps']:.2f} ± "
                               f"{levels[lab]['stat_ps']:.2f} against {levels[lab]['paper_tau_ps']:.2f} ps"
                               + ("" if levels[lab]["kind"] == "value" else " (limit)") for lab in order)
                     + f". {sum(levels[lab]['passed'] for lab in order)} of {len(order)} within the tolerance "
                     f"(worst {worst:.0%} of it). Fitted with SRIM's stopping in the Au as the paper did, the spectra "
                     f"made with physim's (β₀ {100 * out['beta0']:.2f}%, paper about 1.6%): the fits are "
                     f"{100 * (ratios.min() - 1):.0f}–{100 * (ratios.max() - 1):.0f}% above the simulated lifetimes "
                     "(the paper's), because physim's stopping slows ¹⁴³Sm in Au faster than SRIM's (its range is 19% "
                     "shorter); the two upper limits, simulated at the limit, then come out above it")
    _CACHE[key] = out
    return out


def dsam_closure(ref=None, **kw) -> dict:
    """τ back from the simulated ¹⁴³Sm line shapes of :func:`sm143` (side feeding and the band feeding them) with a
    kernel from separate recoils slowed with the same (physim's) stopping, for the two lines of the paper's Fig. 2
    (444.8 and 527.5 keV), their feeders at the simulated lifetimes: within 2σ of the fit."""
    rep = sm143(ref, **kw)
    if not rep.get("available"):
        return rep
    out = {}
    for lab in SM143_CLOSURE:
        lv = rep["levels"][lab]
        dev = lv["own_tau_ps"] - lv["simulated_tau_ps"]
        out[lab] = {"tau_ps": lv["own_tau_ps"], "stat_ps": lv["own_stat_ps"], "true_ps": lv["simulated_tau_ps"],
                    "pull": dev / lv["own_stat_ps"], "e_gamma_kev": lv["e_gamma_kev"]}
    return {"available": True, "levels": out, "worst": max(abs(v["pull"]) / 2 for v in out.values())}


__all__ = ["CR50_PAIRS", "RU98_RINGS", "SM143_INTENSITY", "beta_from_centroids", "closure_one_level",
           "coulex_alignment", "coulex_plunger", "cr50", "ddcm_vs_napytau", "dsam_closure", "dsam_kernel",
           "fusion_decays", "napytau_input", "plunger_areas", "published", "read_napytau_simulated", "ru98",
           "sm143", "sm143_cascade", "three_foil_components"]
