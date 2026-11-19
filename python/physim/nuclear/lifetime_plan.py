"""Planning, predicting and analysing a lifetime measurement for the planner (lifetime wave 3A).

The planner's three lifetime views (:meth:`physim.nuclear.planner.Planner.lifetimes`,
:meth:`~physim.nuclear.planner.Planner.lifetime_spectra`, :meth:`~physim.nuclear.planner.Planner.analyse_lifetimes`)
are thin calls into this module, which reads the setup and calls the pure physics: :mod:`physim.nuclear.lifetimes`
(τ, flight times, sensitive range, peak separation, F(τ), line shapes), :mod:`physim.nuclear.transport` (the
slowing-down kernel), :mod:`physim.nuclear.residues` (fusion residues: their cascade, plunger spectra and DSAM line
shapes) and :mod:`physim.nuclear.lifetime_analysis` (peak areas, DDCM, decay curve, line-shape fit, B(Eλ)).

**The nucleus and its levels.** Coulomb excitation: the excited nucleus (the projectile or the target) with the
setup's level scheme of it (a cascade), or the single state of ``[reaction]`` (its B(E2) gives τ). Fusion-
evaporation: the studied residue (``[run] channel``, else the largest cross section) with its cascade and
side-feeding bands (:func:`physim.nuclear.residues.lifetime_cascade`). τ of each level is
:func:`physim.nuclear.lifetimes.level_taus` (its half-life, else its B(Eλ↓)), with its source.

**The recoil's speed** is :meth:`physim.nuclear.experiment.Experiment.lifetime_recoil`: one typical recoil made in
the middle of the target along the beam (Coulex: a head-on collision; fusion: the compound nucleus's speed). β
"made" is its speed where it is made, β "after the target" its speed leaving the target and its backing (the
plunger's flight). The events and runs have the full spread of speeds and angles; the planning numbers are those of
this one recoil.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from . import lifetimes as lt

#: The lifetimes each method measures (ps): a plunger (RDDS) about 1 ps to 1 ns, DSAM about 10 fs to 1 ps.
RDDS_RANGE_PS = (1.0, 1000.0)
DSAM_RANGE_PS = (0.01, 1.0)
#: Target relative precision on τ when none is given.
DEFAULT_PRECISION = 0.1
#: The shifted and unshifted peaks count as resolved when they are at least this many FWHM apart (two Gaussians of
#: one width are fully apart at about 1.7 FWHM; the shifted peak is broader, from the spread of speeds and angles).
RESOLVED_FWHM = 2.0
#: A Coulex plunger without deorientation is warned about for levels longer-lived than this (ps).
DEORIENTATION_TAU_PS = 1.0
#: A τ given to a level without one, for the cascade's arithmetic (ps): it decays at once.
PROMPT_PS = 1e-6
#: Lines weaker than this share of the strongest are left out of the level list.
MIN_YIELD = 0.01
#: The most levels listed.
MAX_LEVELS = 12

#: Units of B(λ↓) by multipole.
B_UNITS = {"E1": "e2fm2", "E2": "e2fm4", "E3": "e2fm6", "E4": "e2fm8", "M1": "mu_N2", "M2": "mu_N2fm2",
           "M3": "mu_N2fm4"}


def method_for(tau_ps: Optional[float]) -> str:
    """"rdds", "dsam" or "outside" for a mean life (ps): see :data:`RDDS_RANGE_PS` and :data:`DSAM_RANGE_PS`. A
    level without a known τ is "outside"."""
    if tau_ps is None or not math.isfinite(tau_ps) or tau_ps <= 0:
        return "outside"
    if RDDS_RANGE_PS[0] <= tau_ps <= RDDS_RANGE_PS[1]:
        return "rdds"
    if DSAM_RANGE_PS[0] <= tau_ps < DSAM_RANGE_PS[1]:
        return "dsam"
    return "outside"


def relative_precision(x: float) -> float:
    """σ_τ/τ × √N from one distance's unshifted share R = e^(−x), x = t/τ, for N counts of the line (one level, the
    binomial error of R): √(e^x − 1)/x (as :func:`physim.nuclear.lifetimes.sensitive_range` uses it)."""
    if x <= 0:
        return math.inf
    if x > 700:
        return math.inf
    return math.sqrt(math.expm1(x)) / x


def _beta(energy_mev: float, mass_mev: float) -> float:
    if energy_mev is None or energy_mev <= 0:
        return 0.0
    g = 1.0 + energy_mev / mass_mev
    return math.sqrt(1.0 - 1.0 / (g * g))


def _q(x):
    from .quantity import Quantity

    return x if isinstance(x, Quantity) else Quantity.parse(str(x))


# ---------------------------------------------------------------------------------------------------------------
# The nucleus studied: scheme, levels, lines, rates


class Subject:
    """What a lifetime measurement of the setup is about (see the module's notes): the nucleus, its level scheme,
    the levels worth listing (the studied one first), each level's strongest line and its rate in the γ-ray
    detectors, the cascade (τ, feeding, populations) and the recoil."""

    def __init__(self, planner):
        from . import data
        from .gamma import level_names

        self.planner = planner
        exp = planner.experiment
        self.experiment = exp
        self.reaction = exp.reaction
        self.problems: list = []
        self.notes: list = []
        self.scheme = None
        self.nucleus = None
        self.levels: list = []         # indices into the scheme, the studied level first
        self.lines: dict = {}          # level -> (final, γ-ray energy keV)
        self.rates: dict = {}          # level -> full-energy counts per second of its line, or None
        self.yields: dict = {}         # (initial, final) -> relative strength
        self.taus: dict = {}           # level name -> (τ ps, source)
        self.cascade = None            # (taus_ps, feeding, populations) over the scheme's levels (and side bands)
        self.residue_index = None
        self.residue_rate = None
        if exp.excitation is not None:
            self._coulex()
        elif exp.fusion is not None:
            self._fusion()
        else:
            self.problems.append("A lifetime measurement goes with Coulomb excitation or fusion-evaporation, not "
                                 f"with {exp.reaction} scattering: set the reaction type first.")
        if self.scheme is not None:
            self.names = level_names(self.scheme.levels)
            self.taus = lt.level_taus(self.scheme)
        else:
            self.names = []
        recoil = exp.lifetime_recoil()
        self.recoil = recoil
        if recoil is not None:
            self.species = recoil[0]
            self.mass_mev = data.nuclide(recoil[0]).nuclear_mass_mev
            self.beta_made = _beta(recoil[1], self.mass_mev)
            self.beta_after = _beta(recoil[2], self.mass_mev)
        else:
            self.species = self.nucleus
            self.mass_mev = data.nuclide(self.nucleus).nuclear_mass_mev if self.nucleus else None
            self.beta_made = self.beta_after = 0.0
        lt_section = exp.lifetime
        self.method = lt_section.method if lt_section is not None else None
        if self.method == "dsam" or (self.method is None and recoil is not None and recoil[2] <= 0):
            self.beta = self.beta_made
        else:
            self.beta = self.beta_after if self.beta_after > 0 else self.beta_made
        self._studied_first()

    def _studied_first(self) -> None:
        """Put the studied level first: Coulex's excited state; for fusion the level of the strongest line among
        those whose τ suits the setup's method (else any method, else the strongest line)."""
        if self.experiment.fusion is None or len(self.levels) < 2:
            return
        def method(i):
            found = self.taus.get(self.names[i])
            return method_for(found[0] if found else None)

        for ok in ((lambda m: m == self.method), (lambda m: m != "outside")):
            fits = [i for i in self.levels if ok(method(i))]
            if fits:
                first = fits[0]  # the levels are in order of their lines' strength
                self.levels = [first] + [i for i in self.levels if i != first]
                return

    # -- Coulomb excitation ---------------------------------------------------------------------------------------

    def _coulex(self) -> None:
        from .gamma import cascade_excitation, excitation_of

        exp = self.experiment
        exc = exp.excitation
        try:
            cex = cascade_excitation(exp)
        except Exception as e:  # noqa: BLE001 - the single state is still there
            cex = None
            self.notes.append(f"The level scheme's excitation could not be worked out ({e}); the single state is "
                              "used.")
        try:
            self.scheme = cex.scheme if cex is not None else excitation_of(exp).scheme
        except Exception as e:  # noqa: BLE001 - said as a problem
            self.problems.append(f"The excited state could not be worked out: {e}")
            return
        scheme = self.scheme
        self.nucleus = scheme.nuclide
        e_state = 1e3 * exc.energy_mev
        studied = min(range(1, len(scheme.levels)), key=lambda i: abs(scheme.levels[i].energy.value - e_state),
                      default=None)
        if studied is None:
            self.problems.append("The level scheme has no excited level.")
            return
        feeding = _feeding_matrix(scheme)
        n = len(scheme.levels)
        populations = np.zeros(n)
        levels = [studied]
        share = {studied: 1.0}
        if cex is not None:
            pops = self.planner.populations()
            if pops.get("available"):
                for row in pops["levels"]:
                    populations[row["index"]] = row["direct_mb"]
                total_gamma = {(g["initial"], g["final"]): g["gamma_mb"] for g in pops["gammas"]}
                self.yields = dict(total_gamma)
                ref = max((v for (i, _f), v in total_gamma.items() if i == studied), default=0.0)
                for (i, f), v in sorted(total_gamma.items(), key=lambda kv: -kv[1]):
                    if i not in levels and len(levels) < MAX_LEVELS and v >= MIN_YIELD * max(total_gamma.values()):
                        levels.append(i)
                for i in levels:
                    mine = max((v for (a, _f), v in total_gamma.items() if a == i), default=0.0)
                    share[i] = mine / ref if ref > 0 else 0.0
        if not populations.any():
            populations[studied] = 1.0
        self.levels = levels
        for i in levels:
            branches = [t for t in scheme.transitions if t.initial == i and t.energy is not None
                        and t.energy.value > 0]
            if not branches:
                continue
            t = max(branches, key=lambda t: self.yields.get((i, t.final), 0.0) or
                    (t.intensity.value if t.intensity is not None else 0.0))
            self.lines[i] = (t.final, float(t.energy.value))
        taus = lt.level_taus(scheme)
        from .gamma import level_names

        names = level_names(scheme.levels)
        tau_arr = np.array([taus.get(nm, (PROMPT_PS, ""))[0] for nm in names], dtype=float)
        tau_arr[0] = np.inf
        self.cascade = (tau_arr, feeding, populations)
        # The line's rate: the planner's particle–γ coincidences of the excited state, scaled to each line by its
        # γ-ray cross section and the full-energy efficiency at its energy.
        try:
            plan = self.planner.plan()
            coinc = sum(g["coincidence_per_s"] for g in plan["gamma_detectors"])
        except Exception:  # noqa: BLE001 - no rate, said in the warnings
            coinc = None
        eff_ref = _efficiency(exp, e_state)
        for i in levels:
            if i not in self.lines:
                continue
            e = self.lines[i][1]
            if coinc is None or not exp.gamma_detectors or eff_ref <= 0:
                self.rates[i] = None
            else:
                self.rates[i] = coinc * share.get(i, 0.0) * _efficiency(exp, e) / eff_ref

    # -- fusion-evaporation ---------------------------------------------------------------------------------------

    def _fusion(self) -> None:
        from .residues import lifetime_cascade

        exp = self.experiment
        try:
            cas = lifetime_cascade(exp)
        except Exception as e:  # noqa: BLE001 - no scheme: the view says why
            self.problems.append(str(e)[:1].upper() + str(e)[1:])
            return
        self.cascade_info = cas
        scheme = cas["cascade"].scheme
        self.scheme = scheme
        self.nucleus = cas["residue"]
        self.residue_index = cas["index"]
        n = len(scheme.levels)
        ys = {ab: y for ab, y in cas["yields"].items() if ab[0] < n and ab[1] < n}
        self.yields = ys
        if not ys:
            self.problems.append(f"No γ-ray line of {self.nucleus} is fed in the cascade.")
            return
        top = max(ys.values())
        levels = []
        for (a, b), y in sorted(ys.items(), key=lambda kv: -kv[1]):
            if a not in levels and y >= MIN_YIELD * top and len(levels) < MAX_LEVELS:
                levels.append(a)
                self.lines[a] = (b, next(e for i, f, e, _g in cas["lines"] if i == a and f == b))
        self.levels = levels
        self.cascade = (np.asarray(cas["taus_ps"], dtype=float), cas["feeding"], cas["populations"])
        try:
            row = self.planner.residues()["rows"][self.residue_index]
            self.residue_rate = float(row["rate_per_s"])
        except Exception:  # noqa: BLE001 - no rate, said in the warnings
            self.residue_rate = None
        for a in levels:
            b, e = self.lines[a]
            if self.residue_rate is None or not exp.gamma_detectors:
                self.rates[a] = None
            else:
                self.rates[a] = self.residue_rate * ys[(a, b)] * _efficiency(exp, e)

    # -- helpers ----------------------------------------------------------------------------------------------------

    def pick(self, level) -> int:
        """A level of the scheme from its index, its name ("2+₁", "2+1"), or its energy ("783 keV", 783.3)."""
        if self.scheme is None:
            raise ValueError("there is no level scheme to pick a level from")
        n = len(self.scheme.levels)
        if isinstance(level, (int, np.integer)) and not isinstance(level, bool):
            if not 0 < int(level) < n:
                raise ValueError(f"no excited level {level}: the scheme has levels 1 to {n - 1}")
            return int(level)
        text = str(level).strip()
        plain = {nm.translate(_SUBSCRIPTS): i for i, nm in enumerate(self.names)}
        if text in self.names:
            return self.names.index(text)
        if text in plain:
            return plain[text]
        try:
            e = float(_q(text).to("keV")) if any(c.isalpha() for c in text) else float(text)
        except Exception:  # noqa: BLE001
            raise ValueError(f"no level {text!r}: give a level name ({', '.join(self.names[1:8])}, ...), its "
                             "number in the scheme or its energy (\"652 keV\")") from None
        i = min(range(1, n), key=lambda k: abs(self.scheme.levels[k].energy.value - e))
        if abs(self.scheme.levels[i].energy.value - e) <= 2.0:
            return i
        # A γ-ray energy names the level its line leaves (the listed levels' lines first, then any transition).
        lines = [(a, line[1]) for a, line in self.lines.items()] + [
            (t.initial, t.energy.value) for t in self.scheme.transitions if t.energy is not None]
        found = [a for a, eg in lines if abs(eg - e) <= 2.0]
        if found:
            return found[0]
        raise ValueError(f"no level or γ-ray line at {e:g} keV in the scheme of {self.nucleus}")

    def line(self, i: int) -> Optional[tuple]:
        """(final level, γ-ray energy keV) of level ``i``'s strongest line."""
        if i in self.lines:
            return self.lines[i]
        branches = [t for t in self.scheme.transitions if t.initial == i and t.energy is not None
                    and t.energy.value > 0]
        if not branches:
            return None
        t = max(branches, key=lambda t: t.intensity.value if t.intensity is not None else 0.0)
        return (t.final, float(t.energy.value))

    def tau(self, i: int) -> tuple:
        """(τ ps or None, source or None) of level ``i``."""
        found = self.taus.get(self.names[i])
        return (float(found[0]), found[1]) if found else (None, None)


_SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def _feeding_matrix(scheme) -> np.ndarray:
    """(levels × levels) share of each level's decays to each lower level, from the γ-ray intensities and
    conversion (equal shares where the intensities are missing)."""
    n = len(scheme.levels)
    f = np.zeros((n, n))
    for i in range(n):
        branches = [t for t in scheme.transitions if t.initial == i and t.energy is not None and t.energy.value > 0]
        if not branches:
            continue
        w = []
        for t in branches:
            inten = t.intensity.value if t.intensity is not None else 1.0
            conv = t.conversion.value if t.conversion is not None else 0.0
            w.append(max(inten, 0.0) * (1.0 + conv))
        total = sum(w) or float(len(w))
        for t, x in zip(branches, w):
            f[i, t.final] += (x if sum(w) else 1.0) / total
    return f


def _efficiency(experiment, energy_kev: float, theta: Optional[float] = None) -> float:
    """The full-energy-peak efficiency at ``energy_kev`` of all γ-ray detectors, or those at polar angle
    ``theta`` (degrees, to 0.1°)."""
    total = 0.0
    for gd in experiment.gamma_detectors:
        if theta is not None and abs(round(float(_q(gd.theta).to("deg")), 1) - theta) > 0.05:
            continue
        total += float(gd.peak_efficiency(1e-3 * energy_kev, experiment))
    return total


def detector_resolution(experiment) -> list:
    """[(name, θ in degrees, FWHM function of E in keV → keV)] for each γ-ray detector."""
    from .gamma_events import crystals_of

    crystals, bare = crystals_of(experiment, plain=True)
    out = []
    for k, gd in enumerate(experiment.gamma_detectors):
        first = next(j for j, cr in enumerate(crystals) if cr.detector == k)
        out.append((gd.name or f"G{k + 1}", round(float(_q(gd.theta).to("deg")), 1), _kev_fwhm(bare[first].fwhm)))
    return out


def _kev_fwhm(fwhm_mev):
    """A crystal's FWHM (a function of MeV, in MeV) as a function of keV in keV, for numbers or arrays."""
    def fwhm(e_kev):
        v = 1e3 * np.asarray(fwhm_mev(1e-3 * np.asarray(e_kev, dtype=float)), dtype=float)
        return float(v) if v.ndim == 0 else v

    return fwhm


def rings_of(experiment) -> list:
    """[(θ in degrees, FWHM function keV → keV, crystal opening half-angle in rad)] per distinct polar angle."""
    from .gamma_events import crystals_of

    crystals, bare = crystals_of(experiment, plain=True)
    seen = {}
    for k, gd in enumerate(experiment.gamma_detectors):
        th = round(float(_q(gd.theta).to("deg")), 1)
        if th in seen:
            continue
        first = next(j for j, cr in enumerate(crystals) if cr.detector == k)
        cr = crystals[first]
        dist = float(np.linalg.norm(cr.centre)) or 1.0
        seen[th] = (_kev_fwhm(bare[first].fwhm), float(cr.radius) / dist)
    return [(th, *seen[th]) for th in sorted(seen)]


# ---------------------------------------------------------------------------------------------------------------
# DSAM kernel for planning


def planning_kernel(subject: Subject, n: int = 1000, times: int = 120, seed: int = 1):
    """The slowing-down kernel of the typical recoil (made in the middle of the target along the beam at its
    energy) through the target and its backing, with ZBL scattering (:func:`physim.nuclear.transport.kernel`, ``n``
    sample recoils), and the time it takes to stop (ps). (None, None) when the recoil does not stop there."""
    from . import transport
    from .rates import stack

    exp = subject.experiment
    if subject.recoil is None or subject.recoil[1] <= 0:
        return None, None
    layers = stack(exp)
    segs = [transport.Segment(lay.material, lay.thickness, name=lay.name) for lay in layers]
    depth = 0.5 * float(layers[0].thickness)
    e = np.array([subject.recoil[1]])
    t_stop = float(transport.time_to_stop(segs, subject.species, subject.mass_mev, e, np.array([depth]))[0])
    if not math.isfinite(t_stop) or t_stop <= 0:
        return None, None
    grid = np.r_[0.0, np.geomspace(1e-4 * 1.5 * t_stop, 1.5 * t_stop, max(int(times), 3) - 1)]
    k = transport.kernel(segs, subject.species, subject.mass_mev, np.full(n, subject.recoil[1]),
                         np.tile([0.0, 0.0, 1.0], (n, 1)), np.full(n, depth), grid, scatter=True,
                         rng=np.random.default_rng([int(seed), 79]), n=n)
    return k, t_stop


def stopper_time_ps(subject: Subject) -> Optional[float]:
    """The time (ps) the plunger's recoil takes to stop in its stopper (the last foil), from its speed after the
    target (and the degrader, for a differential plunger): decays during it give a DSAM-like tail under the
    unshifted peak."""
    from . import transport
    from .rates import recoil_stack

    exp = subject.experiment
    if exp.lifetime is None or exp.lifetime.method != "rdds" or subject.recoil is None or subject.recoil[2] <= 0:
        return None
    segs = recoil_stack(exp, 0.0)
    names = [s.name for s in segs]
    g1 = names.index("gap")
    foils = [s for s in segs[g1 + 1:] if not s.vacuum]
    if not foils:
        return None
    try:
        e = np.array([subject.recoil[2]])
        if len(foils) > 1:  # through the degrader first
            from .rates import _after, stopping

            first = foils[0]
            e = np.atleast_1d(_after(stopping(subject.species, first.material), e, first.thickness_mg_cm2))
            if e[0] <= 0:
                return float(transport.time_to_stop([first], subject.species, subject.mass_mev,
                                                    np.array([subject.recoil[2]]), np.zeros(1))[0])
        t = float(transport.time_to_stop([foils[-1]], subject.species, subject.mass_mev, e, np.zeros(1))[0])
    except Exception:  # noqa: BLE001 - a planning number only
        return None
    return t if math.isfinite(t) else None


def degraded_beta(subject: Subject) -> Optional[float]:
    """The recoil's speed behind the plunger's first foil, for a degrader or differential plunger (None for a
    stopper, or when it stops in the foil)."""
    exp = subject.experiment
    lt_ = exp.lifetime
    if lt_ is None or lt_.method != "rdds" or lt_.stopper is None or lt_.stopper.kind == "stopper":
        return None
    if subject.recoil is None or subject.recoil[2] <= 0:
        return None
    from .experiment import _slowing

    try:
        _, _, after = _slowing(subject.species, lt_.stopper, subject.recoil[2])
    except ValueError:
        return None
    return _beta(after, subject.mass_mev) if after > 0 else None


# ---------------------------------------------------------------------------------------------------------------
# Planner.lifetimes


def _empty(subject: Optional[Subject], experiment, precision: float) -> dict:
    lt_ = experiment.lifetime
    return {"available": False, "reason": None, "reaction": experiment.reaction,
            "nucleus": subject.nucleus if subject is not None else None,
            "method": lt_.method if lt_ is not None else None,
            "distances_um": list(lt_.distances_um) if lt_ is not None and lt_.method == "rdds" else [],
            "offset_um": float(lt_.offset_um) if lt_ is not None else 0.0,
            "precision": float(precision), "levels": [], "time_to_stop_ps": None, "problems": [], "warnings": []}


def plan(planner, level=None, precision: float = DEFAULT_PRECISION) -> dict:
    """See :meth:`physim.nuclear.planner.Planner.lifetimes`."""
    from .planner import LIFETIME_REACTIONS

    exp = planner.experiment
    if not 0 < precision < 1:
        raise ValueError("the precision is a relative uncertainty on τ, between 0 and 1 (0.1 for 10%)")
    if exp.reaction not in LIFETIME_REACTIONS:
        out = _empty(None, exp, precision)
        out["reason"] = (f"A lifetime measurement goes with Coulomb excitation or fusion-evaporation, not with "
                         f"{exp.reaction} scattering: set the reaction type first.")
        out["problems"] = list(planner.problems)
        return out
    subject = Subject(planner)
    out = _empty(subject, exp, precision)
    out["problems"] = list(planner.problems) + list(subject.problems)
    if not exp.gamma_detectors:
        out["problems"].append("The setup has no γ-ray detectors: a lifetime measurement sees the γ rays.")
    if subject.scheme is None or not subject.levels:
        out["reason"] = (out["problems"][-1] if out["problems"] else
                         "No level of the studied nucleus has a γ-ray line.")
        return out
    if level is not None:
        wanted = [subject.pick(level)]
    else:
        wanted = list(subject.levels)
    lt_section = exp.lifetime
    out.update(available=True, beta_made=subject.beta_made, beta_after_target=subject.beta_after,
               recoil=None if subject.recoil is None else
               {"nuclide": subject.recoil[0], "energy_made_mev": float(subject.recoil[1]),
                "energy_after_target_mev": float(subject.recoil[2])},
               studied=subject.names[subject.levels[0]], notes=list(subject.notes))
    warnings = out["warnings"]
    dsam_possible = (subject.method == "dsam" or (subject.method is None and subject.recoil is not None
                                                  and subject.recoil[2] <= 0))
    kernel = None
    if dsam_possible:
        try:
            kernel, t_stop = planning_kernel(subject)
            out["time_to_stop_ps"] = t_stop
        except Exception as e:  # noqa: BLE001 - F is left out; the warnings say why
            warnings.append(f"The slowing-down kernel could not be worked out ({e}): no F(τ).")
    elif subject.method == "rdds":
        out["time_to_stop_ps"] = stopper_time_ps(subject)
    resolution = detector_resolution(exp) if exp.gamma_detectors else []
    flights = lt_section.flight_distances_um if (lt_section is not None and lt_section.method == "rdds") else []
    no_tau, long_lived = [], []
    for i in wanted:
        name = subject.names[i]
        line = subject.line(i)
        tau, source = subject.tau(i)
        if tau is not None and not math.isfinite(tau):
            tau = None
        meth = method_for(tau)
        e_gamma = line[1] if line else None
        beta = subject.beta
        row = {"level": name, "index": i, "energy_kev": float(subject.scheme.levels[i].energy.value),
               "gamma_kev": e_gamma, "final": subject.names[line[0]] if line else None,
               "tau_ps": tau, "tau_source": source, "beta": beta, "method": meth,
               "recommended_distances_um": [], "sensitive_um": None, "separation": [], "attenuation": None,
               "hours_per_distance": None, "rate_per_s": subject.rates.get(i), "distances_used_um": [],
               "counts_per_distance": None, "note": None}
        if tau is None:
            no_tau.append(name)
            row["note"] = "No known lifetime: give the level a half-life, or its transition a B(E2), in the scheme."
        if meth == "rdds" and beta > 0:
            row["recommended_distances_um"] = lt.recommend_distances(tau, beta)
            row["sensitive_um"] = [float(x) for x in lt.sensitive_range(tau, beta)]
        if e_gamma is not None:
            b_sep = subject.beta if meth != "dsam" else subject.beta_made
            for det, theta, fwhm in resolution:
                sep = float(lt.peak_separation(e_gamma, b_sep, theta)) if b_sep > 0 else 0.0
                fw = fwhm(e_gamma)
                row["separation"].append({"detector": det, "theta_deg": theta, "separation_kev": sep,
                                          "fwhm_kev": fw, "resolved": bool(abs(sep) >= RESOLVED_FWHM * fw)})
        if kernel is not None and tau is not None:
            row["attenuation"] = float(lt.dsam_attenuation(tau, kernel))
        if meth == "rdds" and beta > 0:
            used = flights if flights else [d for d in row["recommended_distances_um"]]
            row["distances_used_um"] = [float(d) for d in used]
            xs = [float(lt.flight_times(d, beta, proper=True)) / tau for d in used if d > 0]
            inv = sum(relative_precision(x) ** -2 for x in xs if math.isfinite(relative_precision(x)))
            if inv > 0:
                n_counts = 1.0 / (precision ** 2 * inv)
                row["counts_per_distance"] = n_counts
                rate = row["rate_per_s"]
                if rate is not None and rate > 0:
                    row["hours_per_distance"] = n_counts / rate / 3600.0
        if tau is not None and tau > DEORIENTATION_TAU_PS:
            long_lived.append(name)
        out["levels"].append(row)
    _plan_warnings(subject, out, no_tau, long_lived)
    return out


def _deorientation_text(long_lived: list) -> str:
    return (f"No deorientation set for this Coulex plunger, and {', '.join(long_lived)} "
            f"live{'s' if len(long_lived) == 1 else ''} longer than {DEORIENTATION_TAU_PS:g} ps: the recoil's "
            "alignment fades in vacuum (hyperfine interaction), so the shifted and unshifted γ rays have different "
            "angular distributions and the shifted/unshifted ratio is biased at each detector angle. Set "
            "lifetime.deorientation (α_k and Γ_k, or g with C_k; Radeck et al., PRC 85, 014301 for 98Ru), or plan "
            "detectors near the angles where the effect is small.")


def deorientation_warning(experiment) -> Optional[str]:
    """The warning for a Coulex plunger without deorientation whose excited levels (the ``[reaction]`` state, or
    the levels of the excited nucleus's scheme) live longer than :data:`DEORIENTATION_TAU_PS`; None otherwise.
    Cheap: it reads only the level scheme, for the planner's warnings banner."""
    lt_ = experiment.lifetime
    if experiment.excitation is None or lt_ is None or lt_.method != "rdds" or lt_.deorientation is not None:
        return None
    from .gamma import cascade_excitation, excitation_of

    try:
        cex = cascade_excitation(experiment)
        scheme = cex.scheme if cex is not None else excitation_of(experiment).scheme
    except Exception:  # noqa: BLE001 - no scheme, no warning: the lifetime plan says what is wrong
        return None
    long_lived = [name for name, (tau, _src) in lt.level_taus(scheme).items() if tau > DEORIENTATION_TAU_PS
                  and math.isfinite(tau) and tau < 1e6]
    return _deorientation_text(long_lived) if long_lived else None


def _plan_warnings(subject: Subject, out: dict, no_tau: list, long_lived: list) -> None:
    exp = subject.experiment
    lt_ = exp.lifetime
    w = out["warnings"]
    studied = out["levels"][0] if out["levels"] else None
    if exp.excitation is not None and lt_ is not None and lt_.method == "rdds" and lt_.deorientation is None \
            and long_lived:
        w.append(_deorientation_text(long_lived))
    if no_tau:
        w.append(f"No known lifetime for {', '.join(no_tau)}: give each a half-life, or its transition a B(E2), "
                 "in the level scheme ([levels]); the plan has no distances or beam time for them.")
    if studied is None:
        return
    if subject.beta <= 0:
        w.append("The recoil does not leave the target with a speed (it stops, or flies backwards): no Doppler "
                 "shift to measure a lifetime from in flight.")
    elif subject.beta < 0.005:
        w.append(f"The recoil is slow (β = {100 * subject.beta:.2f}% of c): the Doppler shifts are small; check "
                 "the peak separations.")
    m = studied["method"]
    if lt_ is not None and m in ("rdds", "dsam") and lt_.method != m:
        word = {"rdds": "a plunger (RDDS)", "dsam": "DSAM"}
        w.append(f"{studied['level']} (τ = {studied['tau_ps']:.3g} ps) suits {word[m]}, but the setup plans "
                 f"{word[lt_.method]}: switch the method with set_lifetime(method=\"{m}\").")
    if studied["tau_ps"] is not None and m == "outside":
        tau = studied["tau_ps"]
        w.append(f"{studied['level']} (τ = {tau:.3g} ps) is outside the plunger (1 ps–1 ns) and DSAM (10 fs–1 ps) "
                 + ("ranges: longer lives need electronic timing (LaBr₃)." if tau > RDDS_RANGE_PS[1] else
                    "ranges: shorter lives need other methods (e.g. Coulex or NRF)."))
    if lt_ is not None and lt_.method == "rdds" and m == "rdds" and studied["sensitive_um"]:
        lo, hi = studied["sensitive_um"]
        inside = [d for d in lt_.flight_distances_um if lo <= d <= hi]
        if len(inside) < 3:
            w.append(f"Only {len(inside)} of the plunger's distances (with the offset) fall in the sensitive range "
                     f"of {studied['level']} ({lo:.3g}–{hi:.3g} µm); recommended: "
                     f"{', '.join(f'{d:g}' for d in studied['recommended_distances_um'])} µm.")
    if lt_ is not None and lt_.method == "rdds" and studied["separation"] and \
            not any(s["resolved"] for s in studied["separation"]):
        w.append(f"The shifted and unshifted peaks of {studied['level']} are not resolved at any γ-ray detector "
                 f"(separation under {RESOLVED_FWHM:g} FWHM): move detectors forward or backward, away from 90°.")
    if studied["rate_per_s"] is None:
        w.append("No line rate: the beam time per distance needs the γ-ray detectors and the reaction's rate.")
    if lt_ is not None and lt_.method == "dsam":
        w.append("DSAM line shapes depend on the stopping powers (about 5–10%), and two of physim's DSAM checks are "
                 "known failures (the report's lifetime section lists them): treat F(τ) as a guide.")


# ---------------------------------------------------------------------------------------------------------------
# Planner.lifetime_spectra


def _gauss_fold(shape: np.ndarray, centres: np.ndarray, fwhm_kev: float) -> np.ndarray:
    """A spectrum on evenly spaced bins folded with a Gaussian of the given FWHM, the same length out."""
    if fwhm_kev <= 0 or len(centres) < 2:
        return shape
    width = centres[1] - centres[0]
    sigma = fwhm_kev / 2.354820045 / width
    half = int(min(math.ceil(5 * sigma), len(shape) - 1))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) / max(sigma, 1e-9)) ** 2)
    k /= k.sum()
    full = np.convolve(shape, k, mode="full")
    return full[half:half + len(shape)]


def _gauss_bins(edges: np.ndarray, centre: float, sigma: float) -> np.ndarray:
    from math import erf

    s = max(sigma, 1e-9) * math.sqrt(2.0)
    c = np.array([erf((x - centre) / s) for x in edges])
    return 0.5 * np.diff(c)


def spectra(planner, distance_um=None, level=None, bins: int = 200, events: int = 10_000) -> dict:
    """See :meth:`physim.nuclear.planner.Planner.lifetime_spectra`."""
    exp = planner.experiment
    out = {"available": False, "reason": None, "method": None, "edges_kev": [], "spectra": [], "level": None,
           "gamma_kev": None, "model": None, "seconds_per_distance": None, "notes": []}
    lt_ = exp.lifetime
    if lt_ is None:
        out["reason"] = ("The setup has no lifetime measurement: add a plunger (set_lifetime) or DSAM "
                         "(set_lifetime(method=\"dsam\") with a backing).")
        return out
    out["method"] = lt_.method
    if not exp.gamma_detectors:
        out["reason"] = "The setup has no γ-ray detectors."
        return out
    subject = Subject(planner)
    if subject.scheme is None or not subject.levels:
        out["reason"] = subject.problems[-1] if subject.problems else "No level with a γ-ray line."
        return out
    i = subject.pick(level) if level is not None else subject.levels[0]
    line = subject.line(i)
    if line is None:
        out["reason"] = f"{subject.names[i]} has no γ-ray line."
        return out
    final, e0 = line
    out.update(level=subject.names[i], gamma_kev=e0)
    beam_s = float(_q(exp.run.beam_time).to("s"))
    if lt_.method == "rdds":
        dists = [float(distance_um)] if distance_um is not None else list(lt_.distances_um)
        if not dists:
            out["reason"] = "The plunger has no distances: give some with set_lifetime(distances=[...])."
            return out
        per = beam_s / max(len(lt_.distances_um), 1)
    else:
        if distance_um is not None:
            raise ValueError("a DSAM measurement has no plunger distances: leave distance_um out")
        dists = [None]
        per = beam_s
    out["seconds_per_distance"] = per
    rings = rings_of(exp)
    try:
        if exp.fusion is not None:
            _fusion_spectra(subject, out, i, final, e0, dists, rings, per, bins, events)
        else:
            _coulex_spectra(subject, out, i, final, e0, dists, rings, per, bins)
    except ValueError as e:
        out["reason"] = str(e)
        out["spectra"] = []
        return out
    out["available"] = True
    return out


def _fusion_spectra(subject, out, i, final, e0, dists, rings, per, bins, events) -> None:
    from .residues import dsam_lineshapes, plunger_spectra

    exp = subject.experiment
    thetas = [th for th, _fw, _o in rings]
    rate = subject.residue_rate or 0.0
    if subject.residue_rate is None:
        out["notes"].append("No residue rate: the spectra are per residue made, not counts.")
        rate, per = 1.0, 1.0
    if exp.lifetime.method == "rdds":
        res = plunger_spectra(exp, distances_um=dists, theta_deg=thetas, bins=bins, events=events,
                              resolution=False, min_yield=0.0)
        found = next((x for x in res["lines"] if (x["initial"], x["final"]) == (i, final)), None)
        if found is None:
            raise ValueError(f"the line {subject.names[i]} → {subject.names[final]} is not fed in the cascade")
        edges = np.asarray(found["edges_kev"], dtype=float)
        centres = 0.5 * (edges[1:] + edges[:-1])
        comp = res["components"]
        out["model"] = ("residues: plunger_spectra (residues drawn by simulate_residues, the cascade with "
                        "side feeding, the shares at each residue's own flight time)")
        for k, d in enumerate(dists):
            for r, (theta, fwhm, _open) in enumerate(rings):
                eff = _efficiency(exp, e0, theta)
                counts = _gauss_fold(found["total"][k, r], centres, fwhm(e0)) * rate * eff * per
                share = float(np.nan_to_num(comp["unshifted"][i][k]))
                b1 = comp["beta_shifted"]
                out["spectra"].append({"distance_um": float(d), "theta_deg": theta, "counts": counts.tolist(),
                                       "shifted_kev": float(lt.doppler_energy(e0, b1, math.cos(math.radians(theta))))
                                       if b1 > 0 else e0,
                                       "unshifted_kev": e0, "unshifted_share": share})
        if comp.get("beta_degraded"):
            out["notes"].append(f"Behind the degrader the residues fly at β = {100 * comp['beta_degraded']:.2f}% "
                                "(the reduced-shifted peak).")
    else:
        res = dsam_lineshapes(exp, theta_deg=thetas, bins=bins, events=min(events, 4000), n=800, times=60,
                              resolution=False, min_yield=0.0)
        found = next((x for x in res["lines"] if (x["initial"], x["final"]) == (i, final)), None)
        if found is None:
            raise ValueError(f"the line {subject.names[i]} → {subject.names[final]} is not fed in the cascade")
        edges = np.asarray(found["edges_kev"], dtype=float)
        centres = 0.5 * (edges[1:] + edges[:-1])
        out["model"] = ("residues: dsam_lineshapes (residues drawn by simulate_residues, slowed with ZBL "
                        "scattering through the target and backing, the cascade with side feeding)")
        for r, (theta, fwhm, _open) in enumerate(rings):
            eff = _efficiency(exp, e0, theta)
            counts = _gauss_fold(found["shapes"][r], centres, fwhm(e0)) * rate * eff * per
            out["spectra"].append({"distance_um": None, "theta_deg": theta, "counts": counts.tolist(),
                                   "shifted_kev": float(found["centroid_kev"][r]), "unshifted_kev": e0,
                                   "attenuation": float(found["attenuation"][r])})
        out["notes"] += list(res.get("notes", []))
    out["edges_kev"] = edges.tolist()


def _share_after(taus, feeding, pops, i: int, t_ps: float) -> float:
    """The share of level ``i``'s decays after proper time ``t_ps`` (the ground state, index 0, left out)."""
    n = len(taus) - 1
    u = np.asarray(lt.unshifted(taus[1:], feeding[1:, 1:], pops[1:], [t_ps]), dtype=float).reshape(n, -1)
    return float(np.nan_to_num(u[i - 1, 0]))


def _coulex_spectra(subject, out, i, final, e0, dists, rings, per, bins) -> None:
    exp = subject.experiment
    lt_ = exp.lifetime
    taus, feeding, pops = subject.cascade
    taus = np.where(np.isfinite(taus) & (taus > 0), taus, np.where(np.isinf(taus), np.inf, PROMPT_PS))
    rate = subject.rates.get(i)
    eff_all = _efficiency(exp, e0)
    if rate is None or eff_all <= 0:
        out["notes"].append("No line rate: the spectra are shares of the line's γ rays, not counts.")
        rate, per, eff_all = 1.0, 1.0, 1.0
    beta = subject.beta
    bmax = max(subject.beta_made, beta, 1e-4)
    edges = e0 * np.linspace(1 - 1.3 * bmax - 0.002, 1 + 1.3 * bmax + 0.002, bins + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    live = [k for k in range(len(taus)) if k > 0]
    lam_taus = np.where(np.isinf(taus), 1e30, taus)
    if lt_.method == "rdds":
        out["model"] = ("analytic (Coulex): lifetimes.unshifted for the level's cascade at the flight time of the "
                        "typical recoil (Experiment.lifetime_recoil: head-on, one speed); the shifted peak a "
                        "Gaussian at its Doppler energy, as wide as the detector's resolution and the spread of "
                        "angles over a crystal's face; decays while the recoil slows in the target or stopper are "
                        "put in the unshifted peak. The spread of speeds and angles of real recoils broadens and "
                        "moves the shifted peak: take a run for that.")
        if beta <= 0:
            raise ValueError("the recoil does not leave the target: no plunger spectra")
        b2 = degraded_beta(subject)
        gap2 = lt_.degrader2.gap_um if (lt_.degrader2 is not None and lt_.stopper is not None
                                        and lt_.stopper.kind == "differential") else None
        for d in dists:
            t1 = float(lt.flight_times(d + lt_.offset_um, beta, proper=True))
            u1 = _share_after(lam_taus, feeding, pops, i, t1)
            if lt_.stopper is not None and lt_.stopper.kind == "degrader":
                u2 = 0.0
            elif gap2 is not None and b2:
                t2 = t1 + float(lt.flight_times(gap2, b2, proper=True))
                u2 = _share_after(lam_taus, feeding, pops, i, t2)
            else:
                u2 = u1
            for theta, fwhm, opening in rings:
                cos = math.cos(math.radians(theta))
                sin = abs(math.sin(math.radians(theta)))
                ring_rate = rate * _efficiency(exp, e0, theta) / eff_all
                fw = fwhm(e0)
                sig0 = fw / 2.354820045
                comps = np.zeros(len(centres))
                e_s = float(lt.doppler_energy(e0, beta, cos))
                sig_s = math.hypot(sig0, e0 * beta * sin * opening / 2.0)
                comps += (1.0 - u1) * _gauss_bins(edges, e_s, sig_s)
                if b2 and u1 - u2 > 0:
                    e_d = float(lt.doppler_energy(e0, b2, cos))
                    comps += (u1 - u2) * _gauss_bins(edges, e_d, math.hypot(sig0, e0 * b2 * sin * opening / 2.0))
                comps += u2 * _gauss_bins(edges, e0, sig0)
                out["spectra"].append({"distance_um": float(d), "theta_deg": theta,
                                       "counts": (comps * ring_rate * per).tolist(), "shifted_kev": e_s,
                                       "unshifted_kev": e0, "unshifted_share": float(u2)})
    else:
        out["model"] = ("Coulex DSAM: lifetimes.lineshape of the level's cascade on the slowing-down kernel of the "
                        "typical recoil (made head-on in the middle of the target, slowed with ZBL scattering "
                        "through the target and backing), folded with the detector's resolution.")
        kernel, _t = planning_kernel(subject, n=2000)
        if kernel is None:
            raise ValueError("the recoil does not stop in the target and its backing: no DSAM line shape")
        for theta, fwhm, _open in rings:
            shape = lt.lineshape(e0, lam_taus[1:], feeding[1:, 1:], pops[1:], kernel, theta, edges)
            row = np.nan_to_num(np.asarray(shape, dtype=float).reshape(len(live), -1)[i - 1])
            ring_rate = rate * _efficiency(exp, e0, theta) / eff_all
            total = row.sum()
            centroid = float(row @ centres / total) if total > 0 else e0
            counts = _gauss_fold(row, centres, fwhm(e0)) * ring_rate * per
            out["spectra"].append({"distance_um": None, "theta_deg": theta, "counts": counts.tolist(),
                                   "shifted_kev": centroid, "unshifted_kev": e0})
    out["edges_kev"] = edges.tolist()


# ---------------------------------------------------------------------------------------------------------------
# Planner.analyse_lifetimes

#: The analysis methods Planner.analyse_lifetimes takes.
ANALYSES = ("ddcm", "decay_curve", "ddcm_three_foil", "lineshape")
#: γ-ray detector rings closer to 90° than this |cos θ| show too small a Doppler shift to analyse.
MIN_COS = 0.25


def _beta_from_shift(counts, edges, e0: float, theta: float, fwhm: float, beta_guess: float) -> Optional[float]:
    """The recoil's speed from where the shifted peak is: the centroid of the counts above the background between
    the unshifted peak and twice the planned shift, E_s/E₀ − 1 ≈ β cos θ (first order). None when nothing is
    there."""
    c = math.cos(math.radians(theta))
    centres = 0.5 * (edges[1:] + edges[:-1])
    shift = e0 * 2.0 * beta_guess * c
    if c > 0:
        lo, hi = e0 + 1.5 * fwhm, e0 + shift + 3 * fwhm
        side = (centres > hi) & (centres <= hi + 10 * fwhm)
    else:
        lo, hi = e0 + shift - 3 * fwhm, e0 - 1.5 * fwhm
        side = (centres < lo) & (centres >= lo - 10 * fwhm)
    m = (centres >= lo) & (centres <= hi)
    if not m.any():
        return None
    bg = float(np.median(counts[side])) if side.any() else 0.0
    y = np.clip(counts[m] - bg, 0.0, None)
    if y.sum() < 10:
        return None
    e_s = float(y @ centres[m] / y.sum())
    beta = (e_s / e0 - 1.0) / c
    if not 0.2 * beta_guess <= beta <= 3.0 * beta_guess:
        return None
    return beta


def _peaks_between(counts, edges, lo: float, hi: float, fwhm: float) -> list:
    """The centroids (keV) of the local maxima of the counts between ``lo`` and ``hi``, smoothed over a FWHM, that
    reach a third of the highest, from low to high energy."""
    centres = 0.5 * (edges[1:] + edges[:-1])
    m = (centres >= min(lo, hi)) & (centres <= max(lo, hi))
    if m.sum() < 5:
        return []
    width = float(np.mean(np.diff(edges)))
    k = max(int(round(fwhm / width)), 1)
    smooth = np.convolve(counts, np.ones(k) / k, mode="same")[m]
    x = centres[m]
    y = smooth - float(np.min(smooth))
    if y.max() <= 0:
        return []
    peaks = [j for j in range(1, len(y) - 1) if y[j] >= y[j - 1] and y[j] > y[j + 1] and y[j] >= y.max() / 3]
    out = []
    for j in peaks:
        near = np.abs(x - x[j]) <= fwhm
        w = np.clip(counts[m][near] - float(np.min(smooth)), 0.0, None)
        out.append(float(w @ x[near] / w.sum()) if w.sum() > 0 else float(x[j]))
    merged = []
    for e in out:  # one peak per FWHM
        if not merged or abs(e - merged[-1]) > fwhm:
            merged.append(e)
    return merged


def _beta_of(e_s: float, e0: float, theta: float) -> float:
    """β from a peak at ``e_s`` of a line ``e0`` at ``theta`` (the relativistic Doppler formula solved for β)."""
    c = math.cos(math.radians(theta))
    r = e_s / e0
    # r = √(1 − β²)/(1 − βc): (r²c² + 1) β² − 2r²c β + (r² − 1) = 0, the root near (r − 1)/c.
    a, b, cc = r * r * c * c + 1.0, -2.0 * r * r * c, r * r - 1.0
    disc = max(b * b - 4 * a * cc, 0.0)
    roots = [(-b + s * math.sqrt(disc)) / (2 * a) for s in (1, -1)]
    guess = (r - 1.0) / c
    return min(roots, key=lambda x: abs(x - guess))


def integrated_areas(counts, edges, e0: float, theta: float, beta: float, fwhm: float) -> tuple:
    """(shifted, unshifted, (error shifted, error unshifted)): the line's counts split at E₀ + half the planned
    shift (at least 3σ of the resolution), each side summed over the background, which is a straight line through
    the counts beyond both ends. The unshifted side runs from E₀ − 3σ and holds the decays while the recoil slows
    in the stopper (the tail towards the shifted side), as :func:`physim.nuclear.lifetimes.unshifted` counts them;
    the shifted side reaches twice the planned shift, so a broad spread of speeds (Coulomb excitation, where the
    particle detector sees beam and recoils over a range of angles) is counted whole, which a Gaussian fit is not."""
    centres = 0.5 * (edges[1:] + edges[:-1])
    sigma = fwhm / 2.354820045
    shift = float(lt.peak_separation(e0, beta, theta))
    sign = 1.0 if shift >= 0 else -1.0
    split = e0 + sign * max(3 * sigma, 0.5 * abs(shift))
    far = e0 + 2.0 * shift + sign * 3 * fwhm
    near = e0 - sign * 3 * sigma
    lo, hi = min(near, far), max(near, far)
    side_lo = (centres >= lo - 8 * fwhm) & (centres < lo - fwhm)
    side_hi = (centres > hi + fwhm) & (centres <= hi + 8 * fwhm)
    side = side_lo | side_hi
    if side.sum() >= 2 and np.ptp(centres[side]) > 0:
        a, b = np.polyfit(centres[side], counts[side], 1)
        bg = np.clip(a * centres + b, 0.0, None)
        n_side = int(side.sum())
    else:
        bg, n_side = np.zeros(len(centres)), 0
    in_u = (centres >= min(near, split)) & (centres < max(near, split))
    in_s = (centres >= min(split, far)) & (centres < max(split, far))
    out = []
    for m in (in_s, in_u):
        raw = float(counts[m].sum())
        b = float(bg[m].sum())
        var_b = (b * b / max(float(counts[side].sum()), 1.0)) if n_side else 0.0
        out.append((max(raw - b, 0.0), math.sqrt(max(raw, 1.0) + var_b)))
    (s, es), (u, eu) = out
    return s, u, (es, eu)


def _areas(counts, edges, energies: list, widths: list) -> list:
    """[(area, error)] of single peaks at ``energies`` (keV) fitted one by one (a Gaussian on a line); ``widths``
    are their starting FWHM (keV)."""
    from .analysis import FWHM_PER_SIGMA, fit_peak

    gaps = [abs(a - b) for k, a in enumerate(energies) for b in energies[k + 1:]]
    half = max(1.5, min(4.0, 0.5 * min(gaps) / max(widths))) if gaps else 4.0
    out = []
    for e, w in zip(energies, widths):
        f = fit_peak(counts, edges, e, w / FWHM_PER_SIGMA, half, max_sigma=6 * w / FWHM_PER_SIGMA)
        out.append((max(float(f.area), 0.0), float(f.area_unc)))
    return out


def _default_series(planner, experiment):
    if planner.folder is None:
        return None
    found = planner.folder.distance_series(experiment)
    return found[0] if found else None


def _dsam_runs(planner, experiment) -> list:
    """The stored runs of this setup (its physics, labels aside) that a DSAM line shape is read from."""
    from . import runs as _runs

    if planner.folder is None:
        return []
    want = experiment.to_dict()
    out = []
    for p in planner.folder.run_folders():
        try:
            run = _runs.RunData.load(p)
        except Exception:  # noqa: BLE001 - a run that cannot be read is left out
            continue
        if run.kind in ("source", "alignment") or _runs.physics_changes(run.setup, want):
            continue
        out.append(run)
    return out


def analyse(planner, level=None, series=None, method: Optional[str] = None) -> dict:
    """See :meth:`physim.nuclear.planner.Planner.analyse_lifetimes`."""
    from . import runs as _runs

    out = {"available": False, "reason": None, "level": None, "method": None, "tau_ps": None, "stat_ps": None,
           "total_ps": None, "budget": {}, "points": [], "curve": {"x": [], "y": [], "kind": None},
           "b_down": None, "b_unit": None, "compare": [], "notes": [], "beta": None, "rings_deg": []}
    if method is not None and method not in ANALYSES:
        raise ValueError(f"unknown analysis {method!r}; the analyses are {', '.join(ANALYSES)}")
    exp = planner.experiment if series is None else series.experiment
    lt_ = exp.lifetime
    if lt_ is None:
        out["reason"] = "The setup has no lifetime measurement ([lifetime])."
        return out
    if series is not None and _runs.physics_changes(planner.experiment.to_dict(), exp.to_dict()):
        from .planner import Planner

        planner = Planner(exp)  # the runs' own setup, for its scheme, rates and recoil
    subject = Subject(planner)
    if subject.scheme is None or not subject.levels:
        out["reason"] = subject.problems[-1] if subject.problems else "No level with a γ-ray line."
        return out
    i = subject.pick(level) if level is not None else subject.levels[0]
    line = subject.line(i)
    if line is None:
        out["reason"] = f"{subject.names[i]} has no γ-ray line."
        return out
    out["level"] = subject.names[i]
    final, e0 = line
    try:
        if lt_.method == "rdds":
            if method == "lineshape":
                raise ValueError("a line-shape fit is for DSAM; a plunger is analysed by ddcm, decay_curve or "
                                 "ddcm_three_foil")
            if series is None:
                series = _default_series(planner, exp)
            if series is None or not len(series):
                out["reason"] = ("No plunger runs of this setup yet: take one run per distance "
                                 "(start_run(distance_um=...)), then analyse.")
                return out
            result = _analyse_rdds(subject, out, i, final, e0, series, method)
        else:
            if method not in (None, "lineshape"):
                raise ValueError("a DSAM measurement is analysed by its line shape (method \"lineshape\")")
            runs = _dsam_runs(planner, exp)
            if not runs:
                out["reason"] = "No runs of this DSAM setup yet: take a run, then analyse its line shapes."
                return out
            result = _analyse_dsam(subject, out, i, final, e0, runs)
    except ValueError as e:
        out["reason"] = str(e)[:1].upper() + str(e)[1:]
        out["points"], out["curve"] = [], {"x": [], "y": [], "kind": None}
        return out
    out.update(available=True, tau_ps=float(result.tau_ps), stat_ps=float(result.stat_ps),
               total_ps=float(result.total_ps), budget={k: float(v) for k, v in result.budget.items()})
    out["notes"] += list(result.notes)
    from .lifetime_analysis import b_from_result

    lam = _multipole(subject.scheme, i, final)
    out["b_multipole"] = lam
    try:
        b, unc = b_from_result(result, subject.scheme, i, final=final, lam=lam)
        out["b_down"] = [float(b), float(unc)]
        out["b_unit"] = B_UNITS.get(lam, lam)
    except ValueError as e:
        out["notes"].append(f"No B({lam}↓): {e}.")
    out["compare"] = compare(subject, i, final, lam)
    return out


def _multipole(scheme, i: int, final: int) -> str:
    """The multipole B(λ↓) is given for: E2 for an E2 or mixed M1+E2 line (or one without a multipolarity), else
    the line's own."""
    tr = next((t for t in scheme.transitions if t.initial == i and t.final == final), None)
    mult = ((tr.multipolarity if tr is not None else "") or "").replace("(", "").replace(")", "").upper()
    if mult in ("", "E2") or "+" in mult:
        return "E2"
    return mult if mult in B_UNITS else "E2"


def _feeder(subject, i: int, final: int) -> Optional[tuple]:
    """(level, γ-ray keV) of the strongest line into level ``i``, when it is at least a tenth of ``i``'s own."""
    mine = subject.yields.get((i, final), 0.0)
    into = [(a, y) for (a, b), y in subject.yields.items() if b == i and a != i]
    if not into or mine <= 0:
        return None
    a, y = max(into, key=lambda kv: kv[1])
    line = subject.line(a)
    if y < 0.1 * mine or line is None or line[0] != i:
        return None
    return a, line[1]


def _rings_used(subject, rings, summed, edges, e0) -> tuple:
    """The rings whose shifted and unshifted peaks stand apart, and the recoil's speed seen in each."""
    used, betas, weights = [], [], []
    beta_plan = subject.beta
    for r, (theta, fwhm, _open) in enumerate(rings):
        fw = fwhm(e0)
        if abs(math.cos(math.radians(theta))) < MIN_COS or abs(float(lt.peak_separation(e0, beta_plan, theta))) < fw:
            continue
        b = _beta_from_shift(summed[r], edges, e0, theta, fw, beta_plan)
        used.append(r)
        betas.append(b if b is not None else beta_plan)
        weights.append(float(summed[r].sum()))
    return used, betas, weights


def _analyse_rdds(subject, out, i, final, e0, series, method):
    from . import dataviews

    exp = series.experiment
    lt_ = exp.lifetime
    kind = lt_.stopper.kind if lt_.stopper is not None else "stopper"
    if method is None:
        method = "ddcm_three_foil" if kind == "differential" else "ddcm"
    if method == "ddcm_three_foil" and kind != "differential":
        raise ValueError("ddcm_three_foil is for a differential (two-foil) plunger")
    if kind == "differential" and method != "ddcm_three_foil":
        out["notes"].append("A differential plunger analysed as a single gap: the reduced-shifted peak is fitted "
                            "as part of the shifted one's surroundings, not as a component.")
    spec = dataviews.distance_spectra(series, doppler=True)
    edges = np.asarray(spec["edges_kev"], dtype=float)
    dists = sorted(spec["by_distance"])
    if len(dists) < 3:
        raise ValueError(f"runs at {len(dists)} distance(s): a lifetime fit needs at least three")
    if subject.beta <= 0:
        raise ValueError("the recoil does not leave the target: no plunger analysis")
    rings = rings_of(exp)
    total = {d: spec["by_distance"][d]["counts"] for d in dists}
    summed = sum(total.values())
    used, betas, weights = _rings_used(subject, rings, summed, edges, e0)
    if not used:
        raise ValueError(f"no γ-ray detector ring sees the shifted and unshifted peaks apart (|cos θ| ≥ {MIN_COS:g} "
                         "and a separation of at least one FWHM)")
    beta = float(np.average(betas, weights=np.maximum(weights, 1.0)))
    beta_err = float(np.std(betas)) if len(betas) > 1 else 0.0
    out["beta"] = beta
    out["rings_deg"] = [rings[r][0] for r in used]
    if method != "ddcm_three_foil":
        out["notes"].append(f"Recoil speed from the shifted peaks' positions: β = {100 * beta:.3f}% (planned "
                            f"{100 * subject.beta:.3f}%); rings at {', '.join(f'{rings[r][0]:g}°' for r in used)} "
                            "added together.")
    if lt_.deorientation is not None:
        out["notes"].append("The setup has deorientation in vacuum: the shifted and unshifted γ rays have different "
                            "angular distributions, and no correction is applied here "
                            "(lifetime_analysis.deorientation_correction, with the angular distribution's a_k): τ "
                            "can be biased (for 98Ru's example by about +10%).")
    feeder = _feeder(subject, i, final)
    common = dict(subject=subject, out=out, i=i, final=final, e0=e0, dists=dists, total=total, edges=edges,
                  rings=rings, used=used, betas=betas, beta=beta, beta_err=beta_err, feeder=feeder, lt_=lt_)
    if method == "ddcm_three_foil":
        return _three_foil(**common)
    return _two_components(method=method, **common)


def _two_components(subject, out, i, final, e0, dists, total, edges, rings, used, betas, beta, beta_err, feeder,
                    lt_, method):
    from . import lifetime_analysis as la

    exp = subject.experiment
    n = len(dists)

    integrate = exp.excitation is not None
    if integrate:
        out["notes"].append("Coulomb excitation: the shifted component is summed over its region (up to twice the "
                            "planned shift) and the unshifted one from E₀ − 3σ to half the shift, both over a straight "
                            "background, as the recoils' speeds spread too widely for one Gaussian.")
    else:
        out["notes"].append("Peak areas from Gaussian fits (lifetime_analysis.peak_areas), the shifted peak's width "
                            "free.")

    def areas(e_line, r, b):
        theta, fwhm, _o = rings[r]
        fw = fwhm(e_line)
        if integrate:
            return [integrated_areas(total[d][r], edges, e_line, theta, b, fw) for d in dists]
        return la.peak_areas([(total[d][r], edges) for d in dists], e_line, b, theta, fw, fwhm_shifted_kev=1.5 * fw)

    s_, u_, vs, vu = (np.zeros(n) for _ in range(4))
    fs_, fu_, fvs, fvu = (np.zeros(n) for _ in range(4))
    feeder_ok = feeder is not None and method == "ddcm"
    for r, b in zip(used, betas):
        for k, (s, u, (es, eu)) in enumerate(areas(e0, r, b)):
            s_[k] += s
            u_[k] += u
            vs[k] += es * es
            vu[k] += eu * eu
        if feeder_ok:
            try:
                for k, (s, u, (es, eu)) in enumerate(areas(feeder[1], r, b)):
                    fs_[k] += s
                    fu_[k] += u
                    fvs[k] += es * es
                    fvu[k] += eu * eu
            except ValueError:
                feeder_ok = False
    tot = s_ + u_
    if np.any(tot <= 0):
        raise ValueError("a distance has no counts in the line: take longer runs")
    s_r = np.sqrt(s_ * s_ * np.maximum(vu, 1.0) + u_ * u_ * np.maximum(vs, 1.0)) / tot ** 2
    if method == "decay_curve":
        res = la.decay_curve(dists, u_ / tot, s_r, beta, offset_um=lt_.offset_um, beta_error=beta_err)
        out["notes"].append("Decay curve of the level alone (one exponential): its feeding is not modelled; the "
                            "DDCM corrects for it.")
        for k, d in enumerate(dists):
            out["points"].append({"distance_um": float(d), "unshifted": float(u_[k] / tot[k]),
                                  "shifted": float(s_[k] / tot[k]), "error": float(s_r[k]), "tau_ps": None,
                                  "tau_err_ps": None})
        out["curve"] = {"x": [row["distance_um"] for row in res.per_distance],
                        "y": [row["model"] for row in res.per_distance],
                        "kind": "unshifted share (decay-curve fit)"}
        out["method"] = "decay_curve"
        return res
    kw = {}
    errors = (np.sqrt(np.maximum(vu, 1.0)), np.sqrt(np.maximum(vs, 1.0)))
    if feeder_ok:
        scale = _feeder_scale(exp, used, rings, e0, feeder[1])
        errors = errors + (np.sqrt(np.maximum(fvu, 1.0)) * scale, np.sqrt(np.maximum(fvs, 1.0)) * scale)
        kw = dict(feeder_unshifted=fu_ * scale, feeder_shifted=fs_ * scale, alpha=_alpha(subject, i, final))
        out["notes"].append(f"Feeding corrected with {subject.names[feeder[0]]} → {subject.names[i]} "
                            f"({feeder[1]:.1f} keV), its intensities corrected for the efficiency.")
    elif feeder is not None:
        out["notes"].append(f"The feeder {subject.names[feeder[0]]} → {subject.names[i]} could not be fitted: no "
                            "feeding correction.")
    res = la.ddcm(dists, u_, s_, errors, beta, offset_um=lt_.offset_um, beta_error=beta_err, **kw)
    for k, d in enumerate(dists):
        row = res.per_distance[k]
        out["points"].append({"distance_um": float(d), "unshifted": float(u_[k] / tot[k]),
                              "shifted": float(s_[k] / tot[k]), "error": float(s_r[k]),
                              "tau_ps": row["tau_ps"], "tau_err_ps": row["error_ps"], "used": row["used"]})
    out["curve"] = {"x": [row["distance_um"] for row in res.per_distance],
                    "y": [row["curve"] for row in res.per_distance], "kind": "shifted intensity (DDCM fit)"}
    out["method"] = "ddcm"
    return res


def _three_foil(subject, out, i, final, e0, dists, total, edges, rings, used, betas, beta, beta_err, feeder, lt_):
    from . import lifetime_analysis as la

    exp = subject.experiment
    n = len(dists)
    b2_plan = degraded_beta(subject)
    if not b2_plan:
        raise ValueError("the recoil's speed behind the degrader is unknown (it stops in the degrader?)")
    b2 = b2_plan * beta / subject.beta
    # Two shifted peaks: the speed in front of the degrader from the outermost, behind it from the one between.
    summed = sum(total.values())
    fast, slow = [], []
    for r in used:
        theta, fwhm, _o = rings[r]
        fw = fwhm(e0)
        c = math.cos(math.radians(theta))
        reach = e0 + 2.0 * float(lt.peak_separation(e0, subject.beta, theta)) + math.copysign(3 * fw, c)
        found = _peaks_between(summed[r], edges, e0 + math.copysign(2 * fw, c), reach, fw)
        if c < 0:
            found = found[::-1]
        if found:
            fast.append(_beta_of(found[-1], e0, theta))
        if len(found) >= 2:
            slow.append(_beta_of(found[-2], e0, theta))
    if fast:
        beta = float(np.mean(fast))
        beta_err = float(np.std(fast)) if len(fast) > 1 else 0.0
        betas = list(fast) if len(fast) == len(used) else [beta] * len(used)
        out["beta"] = beta
        b2 = float(np.mean(slow)) if slow else b2_plan * beta / subject.beta
        out["notes"].append(f"Three-foil: β = {100 * beta:.3f}% from the fully shifted peaks, β′ = {100 * b2:.3f}% "
                            + ("from the reduced-shifted ones." if slow else "scaled from the planned speed behind "
                               "the degrader."))
    comps = {"line": np.zeros((3, n)), "line_var": np.zeros((3, n)), "feed": np.zeros((3, n)),
             "feed_var": np.zeros((3, n))}
    feeder_ok = feeder is not None
    for r, b in zip(used, betas):
        theta, fwhm, _o = rings[r]
        c = math.cos(math.radians(theta))
        lines = [("line", e0)] + ([("feed", feeder[1])] if feeder_ok else [])
        for key, e_line in lines:
            fw = fwhm(e_line)
            energies = [float(lt.doppler_energy(e_line, b, c)), float(lt.doppler_energy(e_line, b2 * b / beta, c)),
                        e_line]
            if min(abs(energies[0] - energies[1]), abs(energies[1] - energies[2])) < 2.5 * fw:
                if key == "line":
                    raise ValueError("the fully shifted, reduced-shifted and unshifted peaks are not apart (under "
                                     "2.5 FWHM) at the rings used")
                feeder_ok = False
                continue
            try:
                for k, d in enumerate(dists):
                    for j, (a, s) in enumerate(_areas(total[d][r], edges, energies, [1.5 * fw, 1.3 * fw, fw])):
                        comps[key][j, k] += a
                        comps[key + "_var"][j, k] += s * s
            except ValueError:
                if key == "line":
                    raise
                feeder_ok = False
    fs, rs, uu = comps["line"]
    gap2 = lt_.degrader2.gap_um if lt_.degrader2 is not None else 0.0
    pairs = np.array([(d, gap2) for d in dists])
    errors = tuple(np.sqrt(np.maximum(comps["line_var"], 1.0)))
    kw = {}
    if feeder_ok:
        scale = _feeder_scale(exp, used, rings, e0, feeder[1])
        ffs, frs, fu = comps["feed"] * scale
        errors = errors + tuple(np.sqrt(np.maximum(comps["feed_var"], 1.0)) * scale)
        kw = dict(feeder_fully_shifted=ffs, feeder_reduced_shifted=frs, feeder_unshifted=fu,
                  alpha=_alpha(subject, i, final))
        out["notes"].append(f"Feeding corrected with {subject.names[feeder[0]]} → {subject.names[i]} "
                            f"({feeder[1]:.1f} keV).")
    elif feeder is not None:
        out["notes"].append(f"The feeder {subject.names[feeder[0]]} → {subject.names[i]} could not be fitted: no "
                            "feeding correction.")
    res = la.ddcm_three_foil(pairs, fs, rs, uu, errors, beta, b2, offset_um=lt_.offset_um, beta_error=beta_err,
                             **kw)
    tot = fs + rs + uu
    if np.any(tot <= 0):
        raise ValueError("a distance has no counts in the line: take longer runs")
    stopper_rows = res.per_distance[n:]
    su = np.sqrt(np.maximum(comps["line_var"][2], 1.0))
    for k, d in enumerate(dists):
        row = stopper_rows[k]
        out["points"].append({"distance_um": float(d), "unshifted": float(uu[k] / tot[k]),
                              "shifted": float((fs[k] + rs[k]) / tot[k]), "error": float(su[k] / tot[k]),
                              "tau_ps": row["tau_ps"], "tau_err_ps": row["error_ps"],
                              "fully_shifted": float(fs[k] / tot[k]), "reduced_shifted": float(rs[k] / tot[k])})
    out["curve"] = {"x": [row["distance_um"] for row in res.per_distance],
                    "y": [row["curve"] for row in res.per_distance],
                    "kind": "decays before each foil (DDCM fit), against x₁ + x₂ β/β′"}
    out["beta_reduced"] = b2
    out["method"] = "ddcm_three_foil"
    return res


def _feeder_scale(experiment, used, rings, e0, e_feeder) -> float:
    """Efficiency correction of the feeder's areas into the level's line: ε(E₀)/ε(E_feeder) over the rings used."""
    a = sum(_efficiency(experiment, e0, rings[r][0]) for r in used)
    b = sum(_efficiency(experiment, e_feeder, rings[r][0]) for r in used)
    return a / b if b > 0 else 1.0


def _alpha(subject, i, final) -> float:
    """The share of level ``i``'s decays that go through the measured transition (the DDCM's α)."""
    try:
        share = float(subject.cascade[1][i, final])
        return share if share > 0 else 1.0
    except Exception:  # noqa: BLE001 - no cascade: one branch
        return 1.0


def _analyse_dsam(subject, out, i, final, e0, runs):
    from . import dataviews
    from . import lifetime_analysis as la
    from .runs import DistanceSeries

    exp = subject.experiment
    spec = dataviews.distance_spectra(DistanceSeries(exp.to_dict(), {0.0: runs}), doppler=True)
    edges = np.asarray(spec["edges_kev"], dtype=float)
    counts = spec["by_distance"][0.0]["counts"]
    rings = rings_of(exp)
    taus, feeding, pops = subject.cascade
    if exp.fusion is not None:
        from .residues import _reached, dsam_lineshapes

        k = dsam_lineshapes(exp, theta_deg=[rings[0][0]], bins=50, events=8000, n=2000, times=120,
                            resolution=False, min_yield=1.0)["kernel"]
        idx, t_r, f_r, p_r = _reached(subject.cascade_info)
        where = np.flatnonzero(idx == i)
        if not len(where):
            raise ValueError(f"{subject.names[i]} is not reached by the cascade")
        t_r = np.where(np.isfinite(t_r) & (t_r > 0), t_r, PROMPT_PS)
        feed = la.Feeding(list(t_r), f_r, list(p_r), int(where[0]))
        beta0 = float(np.mean(k.beta0))
        out["notes"].append("Kernel: residues drawn by simulate_residues, slowed with ZBL scattering through the "
                            "target and backing; the feeders' lifetimes fixed at the scheme's (side feeding as "
                            "the setup's band model).")
    else:
        k, _t = planning_kernel(subject, n=3000)
        if k is None:
            raise ValueError("the recoil does not stop in the target and its backing")
        t = np.where(np.isfinite(taus) & (taus > 0), taus, PROMPT_PS)[1:]
        feed = la.Feeding(list(t), feeding[1:, 1:], list(pops[1:]), i - 1)
        beta0 = subject.beta_made
        out["notes"].append("Kernel: the typical recoil (made head-on in the middle of the target), slowed with ZBL "
                            "scattering through the target and backing.")
    out["beta"] = beta0
    results = []
    for r, (theta, fwhm, _o) in enumerate(rings):
        if abs(math.cos(math.radians(theta))) < MIN_COS + 0.05:
            continue
        span = e0 * (1.3 * beta0 + 0.004) + 5 * fwhm(e0)
        sel = np.flatnonzero((edges[:-1] >= e0 - span) & (edges[1:] <= e0 + span))
        if len(sel) < 20 or counts[r][sel].sum() < 50:
            continue
        lo, hi = int(sel[0]), int(sel[-1]) + 1
        try:
            res = la.fit_lineshape(counts[r][lo:hi], edges[lo:hi + 1], e0, k, theta, feed, fwhm_kev=fwhm)
        except ValueError as e:
            out["notes"].append(f"The line shape at {theta:g}° could not be fitted: {e}.")
            continue
        if math.isfinite(res.stat_ps) and res.stat_ps > 0:
            results.append((theta, res, lo, hi, r))
    if not results:
        raise ValueError("no ring has a line shape to fit (too few counts, or detectors near 90°): take a longer "
                         "run")
    w = np.array([1.0 / x[1].stat_ps ** 2 for x in results])
    tau = float(np.sum(w * np.array([x[1].tau_ps for x in results])) / w.sum())
    stat = float(1.0 / math.sqrt(w.sum()))
    keys = set().union(*(x[1].budget for x in results))
    budget = {kk: float(np.sum(w * np.array([x[1].budget.get(kk, 0.0) for x in results])) / w.sum()) for kk in keys}
    budget["statistics"] = stat
    theta, best, lo, hi, r = max(results, key=lambda x: 1.0 / x[1].stat_ps)
    centres = 0.5 * (edges[lo:hi] + edges[lo + 1:hi + 1])
    out["curve"] = {"x": centres.tolist(), "y": list(best.params.get("model", [])),
                    "kind": f"line shape at {theta:g}° (fit)", "data": counts[r][lo:hi].tolist()}
    out["rings_deg"] = [x[0] for x in results]
    out["method"] = "lineshape"
    notes = [f"Line shapes fitted at {', '.join(f'{x[0]:g}°' for x in results)}; τ is their weighted mean."]
    notes += list(best.notes)
    return la.TauResult(tau, stat, budget, "line shape", None, notes,
                        {"per_ring": [{"theta_deg": x[0], "tau_ps": x[1].tau_ps, "stat_ps": x[1].stat_ps}
                                      for x in results]})


def compare(subject, i: int, final: int, lam: str) -> list:
    """The level's τ and B(λ↓) from elsewhere: the setup's level scheme (its half-life or matrix element, with its
    source: ENSDF, typed, derived), the local ENSDF copy, the setup's Coulex B(E2) (the ``[reaction]`` state), and
    shell-model matrix elements in the scheme."""
    from . import ensdf as _ensdf

    out = []
    scheme = subject.scheme
    e_gamma = subject.line(i)[1]

    def b_of(tau):
        from .lifetime_analysis import TauResult, b_from_result

        try:
            return float(b_from_result(TauResult(float(tau), 0.0, {}, "given"), scheme, i, final=final, lam=lam)[0])
        except (ValueError, ZeroDivisionError):
            return None

    def add(source, tau, b=None):
        if tau is None or not math.isfinite(tau) or tau <= 0:
            return
        if any(x["source"] == source for x in out):
            return
        out.append({"source": source, "tau_ps": float(tau), "b_down": b if b is not None else b_of(tau)})

    tau, source = subject.tau(i)
    if tau is not None:
        add(f"level scheme: {source}", tau)
    exc = subject.experiment.excitation
    if exc is not None and exc.multipolarity == "E2" and \
            abs(scheme.levels[i].energy.value - 1e3 * exc.energy_mev) < 1.0:
        j0 = scheme.levels[0].spin if scheme.levels[0].spin is not None else 0.0
        j = scheme.levels[i].spin if scheme.levels[i].spin is not None else 2.0
        b_down = float(exc.b_up_e2fm) * (2 * j0 + 1) / (2 * j + 1)
        add("setup's Coulex B(E2)", float(lt.tau_from_b(e_gamma, "E2", b_down, branching=_alpha(subject, i, final))),
            b_down)
    for t in scheme.transitions:
        if t.initial != i or t.energy is None:
            continue
        me = scheme.matrix_element(i, t.final, lam)
        if me is None or "shell" not in str(me.source).lower():
            continue
        b = scheme.b(i, t.final, lam)
        if b:
            share = _alpha(subject, i, t.final)
            add(f"shell model ({me.source})", float(lt.tau_from_b(t.energy.value, lam, b, branching=share)),
                float(b) if t.final == final else None)
    if _ensdf.available() and subject.nucleus:
        try:
            from .gamma import level_names
            from .levels import LevelScheme

            e_level = scheme.levels[i].energy.value
            ens = LevelScheme.from_ensdf(subject.nucleus, max_energy_kev=e_level + 50.0)
            taus = lt.level_taus(ens)
            names = level_names(ens.levels)
            k = min(range(len(ens.levels)), key=lambda kk: abs(ens.levels[kk].energy.value - e_level))
            if abs(ens.levels[k].energy.value - e_level) <= 2.0 and names[k] in taus:
                add("ENSDF (local copy)", taus[names[k]][0])
        except Exception:  # noqa: BLE001 - no ENSDF entry: nothing to compare with
            pass
    return out
