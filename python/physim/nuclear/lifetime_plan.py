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
# DSAM backing advice (lifetime wave 4A)

#: The usual DSAM backings, densest first (elements; their densities are the material table's).
BACKING_CANDIDATES = ("Pt", "Au", "W", "Ta", "Pb", "Ag", "Cu", "Ni", "Fe", "C")
#: A backing counts as good for a level when its sensitivity |dF/d ln τ| is at least this share of the best
#: candidate's: σ_τ at most 25% larger for the same precision on F. The setup's backing below it is warned about.
BACKING_GOOD = 0.8
#: The thickness advised is this many times the thickness that stops the typical recoil: room for the recoils'
#: spread of energies and directions and for the stopping's 5–10%.
BACKING_MARGIN = 1.5
#: Sample recoils per candidate's kernel (the same random numbers for every candidate, so they compare fairly).
BACKING_SAMPLES = 250


def _round_up(x: float, digits: int = 2) -> float:
    """``x`` rounded up to ``digits`` significant figures."""
    if not x > 0:
        return 0.0
    scale = 10.0 ** (math.floor(math.log10(x)) - digits + 1)
    return math.ceil(x / scale - 1e-9) * scale


def _um_of(mg_cm2: float, density: Optional[float]) -> Optional[float]:
    return None if not density else 10.0 * mg_cm2 / density


def backing_candidates(subject: Subject, seed: int = 1) -> list:
    """The usual DSAM backings (:data:`BACKING_CANDIDATES`, plus the setup's own backing material when it is not
    one of them) behind the setup's target, each as thick as :data:`BACKING_MARGIN` times what stops the typical
    recoil (made in the middle of the target along the beam): per backing its material, density, the thickness to
    stop the recoil and the thickness advised (mg/cm² and µm), the recoil's time to stop and the slowing-down time
    α (:func:`physim.nuclear.lifetimes.slowing_time`), whether it is the setup's backing material, and its
    slowing-down ``kernel`` (target and backing, ZBL scattering, :data:`BACKING_SAMPLES` recoils). These do not
    depend on τ: :func:`backing_advice` ranks them for a level. Empty when the recoil does not reach the backing."""
    from . import data, transport
    from .rates import _after, stack, stopping

    exp = subject.experiment
    if subject.recoil is None or subject.recoil[1] <= 0:
        return []
    target = stack(exp)[0]
    depth = 0.5 * float(target.thickness)
    e_made = float(subject.recoil[1])
    e_in = float(_after(stopping(subject.species, target.material), e_made, depth))
    if e_in <= 0:
        return []
    setup_mat = exp.target.backing.material_data() if exp.target.backing is not None else None
    materials = [data.material(name) for name in BACKING_CANDIDATES]
    if setup_mat is not None and setup_mat.name not in {m.name for m in materials}:
        materials.append(setup_mat)
    out = []
    for mat in materials:
        stop_mg = float(stopping(subject.species, mat).range(e_in))
        if not stop_mg > 0:
            continue
        thick = BACKING_MARGIN * stop_mg
        segs = [transport.Segment(target.material, float(target.thickness), name="target"),
                transport.Segment(mat, thick, name="backing")]
        t_stop = float(transport.time_to_stop(segs, subject.species, subject.mass_mev, np.array([e_made]),
                                              np.array([depth]))[0])
        if not math.isfinite(t_stop) or t_stop <= 0:
            continue
        n = BACKING_SAMPLES
        grid = np.r_[0.0, np.geomspace(1e-4 * 1.5 * t_stop, 1.5 * t_stop, 79)]
        kernel = transport.kernel(segs, subject.species, subject.mass_mev, np.full(n, e_made),
                                  np.tile([0.0, 0.0, 1.0], (n, 1)), np.full(n, depth), grid, scatter=True,
                                  rng=np.random.default_rng([int(seed), 79]), n=n)
        advised = _round_up(thick)
        rho = mat.density_g_cm3
        out.append({"material": mat.name, "density_g_cm3": rho, "thickness_to_stop_mg_cm2": stop_mg,
                    "thickness_to_stop_um": _um_of(stop_mg, rho), "thickness_mg_cm2": advised,
                    "thickness_um": _um_of(advised, rho), "time_to_stop_ps": t_stop,
                    "slowing_time_ps": float(lt.slowing_time(kernel)),
                    "setup": setup_mat is not None and mat.name == setup_mat.name, "kernel": kernel})
    return out


def backing_advice(tau_ps: float, candidates: list) -> list:
    """The candidates of :func:`backing_candidates` ranked for a level of mean life ``tau_ps``, best first: each
    with ``attenuation`` F(τ) and ``sensitivity`` |dF/d ln τ| (:func:`physim.nuclear.lifetimes.dsam_sensitivity`;
    σ_τ/τ = σ_F / sensitivity, so the larger the better: at most about 1/4, where the slowing-down time matches τ
    and F is near 1/2), ``relative_sensitivity`` (its share of the best's) and ``good`` (at least
    :data:`BACKING_GOOD` of the best's), without the kernels."""
    rows = []
    for c in candidates:
        row = {k: v for k, v in c.items() if k != "kernel"}
        row["attenuation"] = float(lt.dsam_attenuation(tau_ps, c["kernel"]))
        row["sensitivity"] = float(lt.dsam_sensitivity(tau_ps, c["kernel"]))
        rows.append(row)
    rows.sort(key=lambda r: -r["sensitivity"])
    best = rows[0]["sensitivity"] if rows else 0.0
    for r in rows:
        r["relative_sensitivity"] = r["sensitivity"] / best if best > 0 else 0.0
        r["good"] = bool(r["relative_sensitivity"] >= BACKING_GOOD)
    return rows


def _kept(planner, key, make):
    """``make()`` kept in the planner's results for this setup (shared by the plan and the warnings banner)."""
    cache = getattr(planner, "_cache", None)
    if cache is None:
        return make()
    if key not in cache:
        cache[key] = make()
    return cache[key]


# ---------------------------------------------------------------------------------------------------------------
# The beam in the plunger's foils: heating and the foils' own Coulomb excitation (lifetime wave 4A)

#: The beam power (W) a plunger's foil takes: a stretched foil of a few mg/cm², cooled only by radiation, sheds
#: about half a watt from a 1 cm² beam spot at about 900 K (both faces, emissivity about 0.1); above that it heats
#: further, sags and the distance drifts, and thin Au or Cu foils can melt.
STOPPER_POWER_W = 0.5
#: A foil's γ-ray line within this many FWHM of a studied line's unshifted or shifted peak is warned about.
NEAR_FWHM = 3.0
#: Isotopes of a foil's material below this share of its atoms are left out of its Coulomb excitation.
STOPPER_ISOTOPE_SHARE = 0.1
#: A foil's lines weaker than this share of its strongest are left out.
STOPPER_MIN_LINE = 0.01
#: The beam energies (shares of its energy entering the foil, down to where it leaves) at which the foil's Coulomb
#: excitation is worked out and integrated over the foil; below half the energy it is left out (much weaker).
STOPPER_ENERGY_STEPS = (1.0, 0.75, 0.5)
#: The foil's nuclei's levels up to this energy (keV).
STOPPER_MAX_KEV = 3000.0
#: Angle step (degrees) of the foil's Coulomb-excitation cross sections.
STOPPER_ANGLE_STEP = 2.0
_MEV_J = 1.602176634e-13


def _stopper_scheme(experiment, nuclide: str):
    """(level scheme of ``nuclide``, or None, and where it came from or why there is none): a scheme of the setup
    for that nuclide, else the local ENSDF copy's adopted levels up to :data:`STOPPER_MAX_KEV`."""
    for scheme in (experiment.levels or {}).values():
        if getattr(scheme, "nuclide", None) == nuclide and scheme.transitions:
            return scheme, "the setup's level scheme"
    from . import ensdf
    from .levels import LevelScheme

    try:
        return LevelScheme.from_ensdf(nuclide, STOPPER_MAX_KEV), "ENSDF"
    except ensdf.EnsdfMissing as e:
        return None, str(e)
    except Exception as e:  # noqa: BLE001 - no scheme: said, not guessed
        return None, f"its ENSDF levels could not be read ({e})"


def _gamma_cross_sections(ex, step: float) -> dict:
    """{(initial, final): σ of the γ ray, mb} of an :class:`~physim.nuclear.orientation.Excitation`, as its
    ``cross_sections`` works them out (each level's direct excitation probability times the Rutherford cross section
    of its orbit, integrated over CM angles ``step`` apart; the decays share them out), from the populations alone:
    the substates' orientation, which costs ten times more and does not change a cross section, is left out."""
    grid = np.arange(step, 180.0 + step / 2, step)
    th = np.radians(grid)
    probs = {n: np.zeros(len(grid)) for n in ex.paths}
    for k, t in enumerate(grid):
        for n, m in ex._direct(float(t)).items():
            probs[n][k] = float(np.trace(m).real)
    levels = ex.scheme.levels
    populated = {n: float(2 * math.pi * np.sum(p * np.asarray(ex.paths[n][0][2].rutherford_cm(grid)) * np.sin(th))
                          * math.radians(step)) for n, p in probs.items()}
    for n in sorted(set(populated) | {d.initial for d in ex._decays.values()}, key=lambda k: -levels[k].energy.value):
        for (i, f), d in ex._decays.items():
            if i == n and populated.get(n):
                populated[f] = populated.get(f, 0.0) + d.branch * populated[n]
    return {(i, f): populated.get(i, 0.0) * d.gamma for (i, f), d in ex._decays.items()}


def _foil_coulex(experiment, beam: str, mat, e_in: float, e_out: float) -> dict:
    """The foil's own Coulomb excitation by the beam entering it at ``e_in`` (MeV) and leaving at ``e_out`` (0 if
    it stops): per isotope (share ≥ :data:`STOPPER_ISOTOPE_SHARE`) its γ-ray lines with their yield per beam
    particle, the cross sections integrated over the foil from ``e_in`` down to ``e_out`` or half of ``e_in``
    (trapezoids at :data:`STOPPER_ENERGY_STEPS`). Isotopes without level data are listed in ``missing``."""
    from . import data
    from .gamma import level_names
    from .orientation import Excitation
    from .rates import stopping

    st = stopping(beam, mat)
    energies = [e_in * f for f in STOPPER_ENERGY_STEPS if e_in * f > e_out]
    if e_out > 0 and (not energies or energies[-1] > e_out):
        energies.append(e_out)
    per_mg = {}
    for z, a, share in mat.atoms:
        per_mg[(z, a)] = per_mg.get((z, a), 0.0) + share * 1e-3 * 6.02214076e23 / mat.molar_mass
    total = sum(per_mg.values()) or 1.0
    lines, missing, sources, notes = [], [], {}, []
    for key, n_mg in sorted(per_mg.items(), key=lambda kv: -kv[1]):
        if n_mg / total < STOPPER_ISOTOPE_SHARE:
            continue
        nuc = data.nuclide(key).name
        scheme, source = _stopper_scheme(experiment, nuc)
        if scheme is None:
            missing.append({"nuclide": nuc, "why": source})
            continue
        sources[nuc] = source
        names = level_names(scheme.levels)
        sigma = []
        try:
            for e in energies:
                sigma.append(_gamma_cross_sections(Excitation(beam, nuc, e, scheme, excite="target"),
                                                   STOPPER_ANGLE_STEP))
        except Exception as ex:  # noqa: BLE001 - said, not guessed
            missing.append({"nuclide": nuc, "why": f"its Coulomb excitation could not be worked out ({ex})"})
            continue
        ranges = [float(st.range(e)) for e in energies]
        for (i, f) in sorted(set().union(*sigma)) if sigma else []:
            y = 0.0
            for k in range(len(energies) - 1):
                path = ranges[k] - ranges[k + 1]  # mg/cm² between the two energies
                y += 0.5 * (sigma[k].get((i, f), 0.0) + sigma[k + 1].get((i, f), 0.0)) * 1e-27 * n_mg * path
            e_line = scheme.levels[i].energy.value - scheme.levels[f].energy.value
            if y > 0 and e_line > 0:
                lines.append({"nuclide": nuc, "gamma_kev": float(e_line), "transition": f"{names[i]} → {names[f]}",
                              "yield_per_beam": float(y)})
    if lines:
        top = max(x["yield_per_beam"] for x in lines)
        lines = sorted((x for x in lines if x["yield_per_beam"] >= STOPPER_MIN_LINE * top),
                       key=lambda x: -x["yield_per_beam"])
    if energies and energies[-1] > 0 and e_out <= 0 and len(energies) > 1:
        notes.append(f"Coulomb excitation below {energies[-1]:.3g} MeV, half the energy entering the foil, is "
                     "left out (much weaker).")
    return {"lines": lines, "missing": missing, "sources": sources, "energies_mev": energies, "notes": notes}


def _foils(experiment) -> list:
    """The plunger's foils in the beam's way: [(name, material, thickness mg/cm²)] from the recoil's path."""
    from .rates import recoil_stack

    lt_ = experiment.lifetime
    if lt_ is None or lt_.method != "rdds" or lt_.stopper is None:
        return []
    segs = recoil_stack(experiment, 0.0)
    g1 = [s.name for s in segs].index("gap")
    return [(s.name, s.material, float(s.thickness_mg_cm2)) for s in segs[g1 + 1:] if not s.vacuum]


def stopper_beam(planner) -> list:
    """What the beam does in the plunger's foils (the stopper; for a differential plunger the degrader and the
    second foil), kept for the setup: the beam passes the target and flies on into them, in normal and inverse
    kinematics alike. Per foil: ``foil``, ``material``, ``thickness_mg_cm2``, the beam's energy entering and
    leaving (``beam_in_mev``, ``beam_out_mev``; 0 when it stops there), ``deposited_mev``, ``power_w`` (particle
    current × deposited energy) against ``power_limit_w`` (:data:`STOPPER_POWER_W`), and the foil's own Coulomb
    excitation by the beam (:func:`_foil_coulex`): ``lines`` with ``yield_per_beam`` and ``rate_per_s`` (full-energy
    counts per second in all γ-ray detectors, singles), ``missing`` [{``nuclide``, ``why``}] for isotopes without level
    data, ``sources``. Empty without a plunger, or when the beam stops before it."""
    def make():
        from .rates import _after, beam_energy_at, beam_ion, stack, stopping

        exp = planner.experiment
        foils = _foils(exp)
        if not foils:
            return []
        layers = stack(exp)
        e = float(beam_energy_at(exp, len(layers) - 1, [layers[-1].thickness], layers)[0])
        beam = beam_ion(exp)
        pps = float(exp.beam.particles_per_second)
        out = []
        for name, mat, thick in foils:
            if e <= 0:
                break
            e_out = float(_after(stopping(beam, mat), e, thick))
            dep = e - e_out
            row = {"foil": name, "material": mat.name, "thickness_mg_cm2": thick, "beam_in_mev": e,
                   "beam_out_mev": e_out, "deposited_mev": dep, "power_w": pps * dep * _MEV_J,
                   "power_limit_w": STOPPER_POWER_W, "particles_per_s": pps}
            cx = _foil_coulex(exp, beam, mat, e, e_out)
            for line in cx["lines"]:
                line["rate_per_s"] = (line["yield_per_beam"] * pps * _efficiency(exp, line["gamma_kev"])
                                      if exp.gamma_detectors else None)
            row.update(lines=cx["lines"], missing=cx["missing"], sources=cx["sources"], notes=cx["notes"])
            out.append(row)
            e = e_out
        return out

    return _kept(planner, ("lifetimes", "stopper"), make)


def nucleus_lines(experiment) -> tuple:
    """([(level name, γ-ray keV)], β in the plunger's gap) of the studied nucleus, for the foils' lines to be held
    against: Coulomb excitation's every transition of the excited nucleus's scheme; fusion-evaporation's studied
    residue's lines fed with at least :data:`MIN_YIELD` of the strongest. Cheap: no rates, for the banner."""
    from . import data

    lines = []
    if experiment.excitation is not None:
        from .gamma import cascade_excitation, excitation_of, level_names

        try:
            cex = cascade_excitation(experiment)
            scheme = cex.scheme if cex is not None else excitation_of(experiment).scheme
        except Exception:  # noqa: BLE001 - no lines: nothing to hold the foils' lines against
            scheme = None
        if scheme is not None:
            names = level_names(scheme.levels)
            lines = [(names[t.initial], float(t.energy.value)) for t in scheme.transitions
                     if t.energy is not None and t.energy.value > 0]
    elif experiment.fusion is not None:
        from .gamma import level_names
        from .residues import lifetime_cascade

        try:
            cas = lifetime_cascade(experiment)
        except Exception:  # noqa: BLE001
            cas = None
        if cas is not None and cas["yields"]:
            names = level_names(cas["cascade"].scheme.levels)
            top = max(cas["yields"].values())
            for i, f, e, _g in cas["lines"]:
                if cas["yields"].get((i, f), 0.0) >= MIN_YIELD * top and i < len(names):
                    lines.append((names[i], float(e)))
    recoil = experiment.lifetime_recoil()
    beta = _beta(recoil[2], data.nuclide(recoil[0]).nuclear_mass_mev) if recoil is not None else 0.0
    return lines, beta


def stopper_warnings(planner) -> list:
    """The warnings about the beam in the plunger's foils (:func:`stopper_beam`): a foil taking more beam power
    than :data:`STOPPER_POWER_W`, a foil's own Coulomb-excitation line within :data:`NEAR_FWHM` FWHM of a studied
    line's unshifted or shifted peak (:func:`nucleus_lines`), and a foil whose nuclei have no level data (said,
    not guessed). The planner's banner and :meth:`~physim.nuclear.planner.Planner.lifetimes` show the same."""
    exp = planner.experiment
    try:
        foils = stopper_beam(planner)
    except Exception as e:  # noqa: BLE001 - a planning aid: say it could not be worked out
        return [f"The beam in the plunger's foils could not be worked out ({e})."]
    if not foils:
        return []
    out = []
    current_pna = float(exp.beam.particle_current_pna)
    for f in foils:
        where = f"the plunger's {f['thickness_mg_cm2']:.3g} mg/cm² {f['material']} {f['foil']}"
        if f["power_w"] > f["power_limit_w"]:
            how = "it stops there" if f["beam_out_mev"] <= 0 else f"it leaves at {f['beam_out_mev']:.3g} MeV"
            out.append(f"The beam deposits {f['power_w']:.2g} W in {where} ({current_pna:.3g} pnA × "
                       f"{f['deposited_mev']:.3g} MeV of its {f['beam_in_mev']:.3g} MeV; {how}): more than the "
                       f"{f['power_limit_w']:g} W a thin stretched foil cooled only by radiation takes (about half a "
                       "watt from a 1 cm² beam spot at 900 K), so it heats, sags and the distance drifts. Keep the "
                       f"current under {current_pna * f['power_limit_w'] / f['power_w']:.2g} pnA, spread the beam "
                       "spot, or use a cooled stopper.")
    lines, beta = nucleus_lines(exp)
    resolution = detector_resolution(exp) if exp.gamma_detectors else []
    for f in foils:
        where = f"the plunger's {f['material']} {f['foil']}"
        near = []
        for x in f["lines"]:
            hit = _near(x["gamma_kev"], lines, beta, resolution)
            if hit is not None:
                near.append((x, hit))
        if near:
            parts = []
            for x, (level, e_line, what, d, fw) in near[:4]:
                rate = f", {x['rate_per_s']:.3g}/s in the γ-ray detectors (singles)" if x.get("rate_per_s") else ""
                parts.append(f"{x['nuclide']} {x['gamma_kev']:.1f} keV ({x['transition']}{rate}) is {abs(d):.1f} keV "
                             f"from the {e_line:.1f} keV line of {level} {what} (FWHM {fw:.2g} keV)")
            out.append(f"The beam ({f['beam_in_mev']:.3g} MeV) Coulomb-excites {where}, and its lines sit near the "
                       "studied ones: " + "; ".join(parts) + ". Gate on the particles or on the feeding γ rays, or "
                       "choose a stopper material without lines there.")
        if f["missing"]:
            nucs = ", ".join(m["nuclide"] for m in f["missing"])
            out.append(f"No level data for {where} ({nucs}): its own Coulomb excitation by the beam "
                       f"({f['beam_in_mev']:.3g} MeV entering it) could not be checked for lines near the studied "
                       "ones. With a local ENSDF copy physim checks it.")
    return out


def _near(e_foil: float, lines: list, beta: float, resolution: list) -> Optional[tuple]:
    """The closest studied line's peak within :data:`NEAR_FWHM` FWHM of a foil's line at ``e_foil`` keV:
    (level, its γ-ray keV, "unshifted" or "shifted at θ°", distance keV, FWHM keV), or None. The unshifted peak and
    the shifted one at each γ-ray detector's angle (speed ``beta``); without detectors, 2 keV FWHM."""
    best = None
    dets = resolution or [("", None, lambda e: 2.0)]
    for level, e0 in lines:
        for _name, theta, fwhm in dets:
            fw = float(fwhm(e0))
            peaks = [("unshifted", e0)]
            if theta is not None and beta > 0:
                peaks.append((f"shifted at {theta:g}°", e0 + float(lt.peak_separation(e0, beta, theta))))
            for what, e in peaks:
                d = e_foil - e
                if abs(d) <= NEAR_FWHM * fw and (best is None or abs(d) / fw < abs(best[3]) / best[4]):
                    best = (level, e0, what, d, fw)
    return best


# ---------------------------------------------------------------------------------------------------------------
# Planner.lifetimes


def _empty(subject: Optional[Subject], experiment, precision: float) -> dict:
    lt_ = experiment.lifetime
    return {"available": False, "reason": None, "reaction": experiment.reaction,
            "nucleus": subject.nucleus if subject is not None else None,
            "method": lt_.method if lt_ is not None else None,
            "distances_um": list(lt_.distances_um) if lt_ is not None and lt_.method == "rdds" else [],
            "offset_um": float(lt_.offset_um) if lt_ is not None else 0.0,
            "precision": float(precision), "levels": [], "time_to_stop_ps": None, "problems": [], "warnings": [],
            "stopper": []}


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
               "counts_per_distance": None, "note": None, "backing_advice": None}
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
    # DSAM backing advice for the levels DSAM would measure (every level with a τ in a DSAM setup).
    advised = [r for r in out["levels"] if r["tau_ps"] is not None
               and (r["method"] == "dsam" or (subject.method == "dsam" and r["method"] != "outside"))]
    if advised:
        try:
            candidates = _kept(planner, ("lifetimes", "backings"), lambda: backing_candidates(subject))
        except Exception as e:  # noqa: BLE001 - no advice; the warnings say why
            candidates = []
            warnings.append(f"The DSAM backing advice could not be worked out ({e}).")
        for r in advised:
            r["backing_advice"] = backing_advice(r["tau_ps"], candidates) if candidates else None
    out["stopper"] = [{k: v for k, v in f.items()} for f in stopper_beam(planner)] \
        if subject.method == "rdds" else []
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
        w.extend(_backing_warnings(subject, studied))
    if lt_ is not None and lt_.method == "rdds":
        w.extend(stopper_warnings(subject.planner))


def _backing_warnings(subject: Subject, studied: dict) -> list:
    """A DSAM setup's backing against the advice for the studied level: warned about when its sensitivity is under
    :data:`BACKING_GOOD` of the best candidate's. (A backing that does not stop the recoil, or none when the recoil
    leaves the target, is already a problem of the setup.)"""
    advice = studied.get("backing_advice")
    if not advice or subject.experiment.target.backing is None:
        return []
    best = advice[0]
    tau = studied["tau_ps"]
    best_text = (f"{best['material']} (F = {best['attenuation']:.2f}, |dF/d ln τ| = {best['sensitivity']:.3f}; "
                 f"{best['thickness_mg_cm2']:.2g} mg/cm²" + (f" = {best['thickness_um']:.2g} µm" if best.get(
                     "thickness_um") else "") + f", slowing-down time {best['slowing_time_ps']:.2g} ps)")
    out = []
    mine = next((r for r in advice if r["setup"]), None)
    if mine is not None and not mine["good"]:
        out.append(f"The {mine['material']} backing gives {studied['level']} (τ = {tau:.3g} ps) F = "
                   f"{mine['attenuation']:.2f} with |dF/d ln τ| = {mine['sensitivity']:.3f}, "
                   f"{100 * mine['relative_sensitivity']:.0f}% of the best of the usual backings, {best_text}: σ_τ "
                   f"would be {1 / max(mine['relative_sensitivity'], 1e-9):.2g} times larger for the same line-shape "
                   f"precision. Its slowing-down time is {mine['slowing_time_ps']:.2g} ps against τ.")
    return out


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
#: The efficiency ratio k = ε_shifted/ε_unshifted per ring of a plunger analysis: from the setup ("model": the Lorentz
#: boost of the solid angle and the full-energy efficiency at the two energies), fitted from the runs ("fit": I_u +
#: I_s/k the same at every distance, as Radeck et al.; biased upwards at low counts, ~+4% on τ for ⁹⁸Ru's runs of
#: 200 000 reactions), or 1 ("off").
EFFICIENCY_RATIO = "model"


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
        out["notes"].append(f"Rings at {', '.join(f'{rings[r][0]:g}°' for r in used)} added together; the "
                            f"shifted peaks are where β = {100 * beta:.3f}% puts them (planned "
                            f"{100 * subject.beta:.3f}%).")
    feeder = _feeder(subject, i, final)
    common = dict(subject=subject, out=out, i=i, final=final, e0=e0, dists=dists, total=total, edges=edges,
                  rings=rings, used=used, betas=betas, beta=beta, beta_err=beta_err, feeder=feeder, lt_=lt_)
    if method == "ddcm_three_foil":
        return _three_foil(**common)
    real = np.array([sum(float(getattr(run, "real_s", 0.0) or 0.0) for run in series.runs[d]) for d in dists])
    return _two_components(method=method, real_s=real, **common)


def _ring_cos(experiment, thetas) -> dict:
    """The mean cos θ over the faces of each ring's crystals: what a Doppler shift at the ring measures (a little
    less in magnitude than the cosine of the ring's angle, for crystals seen over ±13°)."""
    from .lifetime_reproduction import ring_cos

    try:
        return ring_cos(experiment, thetas)
    except Exception:  # noqa: BLE001 - a ring without crystals of its own: the nominal angle
        return {th: math.cos(math.radians(th)) for th in thetas}


def _alignment(subject, i: int, final: int) -> tuple:
    """({2: a₂, 4: a₄}, note) of the measured γ rays' unattenuated angular distribution for the deorientation
    correction, or (None, why not). Coulex's first excited state only: physim's Coulex particle–γ correlation for
    the particles the setup's detectors see (:func:`lifetime_reproduction.coulex_alignment`, both nuclei: physim's
    stopper does not take the scattered beam, so the detectors see the beam and the recoils, as the runs have
    them). Kept with the planner until the setup changes."""
    exp = subject.experiment
    if exp.excitation is None:
        return None, "the deorientation correction needs the alignment of a Coulomb-excited state"
    if (i, final) != (1, 0):
        return None, "the deorientation correction is made for the first excited state's line only"
    cache = getattr(subject.planner, "_cache", {})
    key = ("lifetime", "alignment")
    if key not in cache:
        from .lifetime_reproduction import coulex_alignment

        try:
            a = coulex_alignment(exp, particle="both")
            cache[key] = ({2: float(a[2]), 4: float(a[4])}, None)
        except Exception as e:  # noqa: BLE001 - no correlation (no crystals at 90° and elsewhere, say)
            cache[key] = (None, f"no alignment for the deorientation correction ({e})")
    return cache[key]


def _component_areas(counts_by_distance, edges, e_line, theta, beta, fwhm, cos_eff) -> dict:
    """Per distance (arrays): ``unshifted``, ``shifted``, ``shifted_beta`` (shifted counts × their speed) and their
    errors at one ring, from :func:`lifetime_analysis.window_areas` when the peaks are apart, else from
    :func:`lifetime_analysis.peak_areas`' fits (the speed of every shifted count then ``beta``). Also ``windows``."""
    from . import lifetime_analysis as la

    keys = ("unshifted", "shifted", "shifted_beta", "err_unshifted", "err_shifted", "err_shifted_beta")
    out = {k: [] for k in keys}
    windows = abs(float(lt.peak_separation(e_line, beta, theta))) >= la.WINDOW_MIN_SEPARATION * fwhm
    for counts in counts_by_distance:
        if windows:
            w = la.window_areas(counts, edges, e_line, beta, theta, fwhm, cos_theta=cos_eff)
        else:
            (s, u, (es, eu)), = la.peak_areas([(counts, edges)], e_line, beta, theta, fwhm, method="fit",
                                              fwhm_shifted_kev=1.5 * fwhm)
            w = {"unshifted": u, "shifted": s, "err_unshifted": eu, "err_shifted": es, "shifted_beta": s * beta,
                 "err_shifted_beta": es * beta}
        for k in keys:
            out[k].append(float(w[k]))
    out = {k: np.asarray(v) for k, v in out.items()}
    out["windows"] = windows
    # The spread of speeds from the largest distance, where nearly every decay is in flight.
    out["speeds"] = w.get("speeds") if windows else (np.array([beta]), np.array([1.0]))
    return out


def _efficiency_ratio(u, s, eu, es) -> Optional[float]:
    """k = ε_shifted/ε_unshifted from I_u + I_s/k being the same at every distance (the runs' counts per real
    second), or None when the distances do not pin it (the shares hardly change, or k comes out at the edge of 0.5–2).
    Both intensities are noisy, so k minimises χ²(k) = Σ (I_u + I_s/k − c)² / (σ_u² + σ_s²/k²), c profiled: a
    regression of I_u on I_s alone would take I_s as exact and pull k high at low counts."""
    from .lifetime_analysis import _minimise_1d

    u, s, eu, es = (np.asarray(x, dtype=float) for x in (u, s, eu, es))
    if np.ptp(s / np.maximum(u + s, 1e-300)) < 0.3:
        return None

    def chi2(ln_k):
        k = math.exp(ln_k)
        w = 1.0 / np.maximum(eu * eu + (es / k) ** 2, 1e-300)
        total = u + s / k
        c = float(np.sum(w * total) / np.sum(w))
        return float(np.sum(w * (total - c) ** 2))  # the effective-variance χ² (Orear, NIM 203, 157 (1982))

    lo, hi = math.log(0.5), math.log(2.0)
    try:
        ln_k = _minimise_1d(chi2, lo, hi)[0]
    except Exception:  # noqa: BLE001 - a fit that does not settle
        return None
    if not math.isfinite(ln_k) or min(ln_k - lo, hi - ln_k) < 1e-3:
        return None
    return math.exp(ln_k)


def _two_components(subject, out, i, final, e0, dists, total, edges, rings, used, betas, beta, beta_err, feeder,
                    lt_, method, real_s=None):
    """DDCM or decay curve from the shifted and unshifted components at the rings ``used``.

    Per ring: the components in windows (:func:`lifetime_analysis.window_areas`, the stopper's slowing tail with
    the unshifted counts), per real second of the runs at each distance; with deorientation in vacuum, each divided
    by its angular distribution, W_U = 1 + Σ a_k G_k(t) P_k and W_S with G̃_k (Radeck et al.), a_k from physim's
    Coulex correlation and τ iterated; the shifted ones divided by the efficiency ratio k = ε_s/ε_u
    (:data:`EFFICIENCY_RATIO`). The rings are then added. The DDCM takes the shifted counts weighted by
    their own speed (τ(x) = I_u / (c d(Σ β I_s)/dx), exact for any spread of speeds, which Coulex gives: the
    particle detector sees the beam and the recoils over a range of angles), normalised by I_u + I_s/k."""
    from . import lifetime_analysis as la

    exp = subject.experiment
    n = len(dists)
    d_arr = np.asarray(dists, dtype=float)
    real = np.asarray(real_s if real_s is not None else np.ones(n), dtype=float)
    per_s = np.where(real > 0, float(np.mean(real[real > 0])) / np.where(real > 0, real, 1.0), 1.0) \
        if np.any(real > 0) else np.ones(n)
    thetas = [rings[r][0] for r in used]
    cos_eff = _ring_cos(exp, thetas)
    lines = [("line", e0)] + ([("feed", feeder[1])] if feeder is not None and method == "ddcm" else [])
    comps = {}
    feeder_ok = len(lines) == 2
    for key, e_line in lines:
        for r, b in zip(used, betas):
            theta, fwhm, _o = rings[r]
            try:
                a = _component_areas([total[d][r] for d in dists], edges, e_line, theta, b, fwhm(e_line),
                                     cos_eff[theta])
            except ValueError:
                if key == "line":
                    raise
                feeder_ok = False
                break
            for k in ("unshifted", "shifted", "shifted_beta", "err_unshifted", "err_shifted", "err_shifted_beta"):
                a[k] = a[k] * per_s
            comps[(key, r)] = a
    windows = all(comps[("line", r)]["windows"] for r in used)
    if windows:
        out["notes"].append("Peak areas in windows (lifetime_analysis.window_areas): the unshifted one from E₀ − 3σ "
                            "to half the shift, with the decays while the recoil slows in the stopper; the shifted "
                            "one up to twice the planned shift; both over a straight background. Per real second of "
                            "the runs at each distance.")
    else:
        out["notes"].append("Peak areas from Gaussian fits (lifetime_analysis.peak_areas, method \"fit\"): the peaks "
                            "are too close for windows at some ring; the shifted counts all take the planned speed.")

    # Deorientation in vacuum: the alignment the unshifted and shifted γ rays keep differs.
    a_k = None
    dor = lt_.deorientation
    if dor is not None:
        a_k, why = _alignment(subject, i, final)
        if a_k is None:
            out["notes"].append(f"The setup has deorientation in vacuum, but {why}: no correction, τ can be biased.")

    def attenuation(tau_ps, b):
        if a_k is None:
            return None
        t = lt.flight_times(d_arr + lt_.offset_um, b, proper=True)
        g_u = lt.deorientation(t, dor.alpha, dor.gamma_per_ps)
        g_s = lt.deorientation_shifted(t, max(tau_ps, 1e-6), dor.alpha, dor.gamma_per_ps)
        return {r: (la._w(rings[r][0], a_k, g_u), la._w(rings[r][0], a_k, g_s)) for r in used}

    def summed(w_rings):
        """The rings added after the deorientation and efficiency corrections: arrays and the k per ring."""
        acc = {k: np.zeros(n) for k in ("u", "s", "sb", "vu", "vs", "vsb")}
        ks = {}
        for r in used:
            a = comps[("line", r)]
            w_u, w_s = w_rings[r] if w_rings is not None else (1.0, 1.0)
            u, eu = a["unshifted"] / w_u, a["err_unshifted"] / w_u
            s, es = a["shifted"] / w_s, a["err_shifted"] / w_s
            sb, esb = a["shifted_beta"] / w_s, a["err_shifted_beta"] / w_s
            k = (_efficiency_ratio(u, s, eu, es) if EFFICIENCY_RATIO == "fit"
                 else k_model.get(r) if EFFICIENCY_RATIO == "model" else None)
            ks[r] = k
            k = k or 1.0
            acc["u"] += u
            acc["s"] += s / k
            acc["sb"] += sb / k
            acc["vu"] += eu * eu
            acc["vs"] += (es / k) ** 2
            acc["vsb"] += (esb / k) ** 2
            values, weights = a["speeds"]
            acc.setdefault("speeds", []).append((values, weights / (k * np.atleast_1d(w_s)[-1])))
        return acc, ks

    def fit(acc):
        u, s, sb = acc["u"], acc["s"], acc["sb"]
        tot = u + s
        if np.any(tot <= 0):
            raise ValueError("a distance has no counts in the line: take longer runs")
        b_bar = float(np.sum(sb) / np.sum(s)) if np.sum(s) > 0 else beta
        if not 0 < b_bar < 1:
            b_bar = beta
        if method == "decay_curve":
            s_r = np.sqrt(s * s * np.maximum(acc["vu"], 1.0) + u * u * np.maximum(acc["vs"], 1.0)) / tot ** 2
            speeds = (np.concatenate([v for v, _w in acc["speeds"]]), np.concatenate([w for _v, w in acc["speeds"]]))
            # A free normalisation (Radeck et al.'s a(1 − e^(−t/τ))): the share stays below 1 at the shortest
            # distances (decays in the target and over the offset, the efficiency ratio's own uncertainty).
            return la.decay_curve(dists, u / tot, s_r, b_bar, offset_um=lt_.offset_um, beta_error=beta_err,
                                  speeds=speeds, free_norm=True), b_bar
        errors = (np.sqrt(np.maximum(acc["vu"], 1.0)), np.sqrt(np.maximum(acc["vsb"], 1e-300)) / b_bar)
        kw = {}
        if feeder_ok:
            scale = _feeder_scale(exp, used, rings, e0, feeder[1])
            fu = sum(comps[("feed", r)]["unshifted"] for r in used) * scale
            fs = sum(comps[("feed", r)]["shifted"] for r in used) * scale
            fvu = sum(comps[("feed", r)]["err_unshifted"] ** 2 for r in used) * scale ** 2
            fvs = sum(comps[("feed", r)]["err_shifted"] ** 2 for r in used) * scale ** 2
            errors = errors + (np.sqrt(np.maximum(fvu, 1.0)), np.sqrt(np.maximum(fvs, 1.0)))
            kw = dict(feeder_unshifted=fu, feeder_shifted=fs, alpha=_alpha(subject, i, final))
        res = la.ddcm(dists, u, sb / b_bar, errors, b_bar, norm=tot, offset_um=lt_.offset_um, beta_error=beta_err,
                      **kw)
        return res, b_bar

    ring_betas = [float(np.sum(comps[("line", r)]["shifted_beta"]) / np.sum(comps[("line", r)]["shifted"]))
                  for r in used if np.sum(comps[("line", r)]["shifted"]) > 0]
    if len(ring_betas) > 1:  # the rings' mean speeds differ by the Doppler shift's systematics (cos θ, background)
        beta_err = float(np.std(ring_betas))
    # The efficiency ratio the setup gives: the Lorentz boost of the solid angle, (1 − β²)/(1 − β cos θ)², times the
    # full-energy efficiency at the shifted energy over that at E₀, at each ring's mean speed.
    k_model = {}
    for r in used:
        a = comps[("line", r)]
        if np.sum(a["shifted"]) <= 0:
            continue
        b_r = float(np.sum(a["shifted_beta"]) / np.sum(a["shifted"]))
        c_r = cos_eff[rings[r][0]]
        e_s = float(lt.doppler_energy(e0, b_r, c_r))
        eff0 = _efficiency(exp, e0, rings[r][0])
        if eff0 > 0 and 0 < b_r < 1:
            k_model[r] = (1 - b_r * b_r) / (1 - b_r * c_r) ** 2 * _efficiency(exp, e_s, rings[r][0]) / eff0
    acc, ks = summed(None)
    res, b_bar = fit(acc)
    if a_k is not None:
        for _ in range(6):  # G̃_k of the shifted component depends on τ: iterate
            acc, ks = summed(attenuation(res.tau_ps, b_bar))
            new, b_bar = fit(acc)
            done = abs(new.tau_ps - res.tau_ps) <= 1e-3 * res.tau_ps
            res = new
            if done:
                break
        out["notes"].append(f"Deorientation in vacuum corrected: each component divided by its angular distribution "
                            f"(a₂ = {a_k[2]:.3f}, a₄ = {a_k[4]:.3f} from physim's Coulex correlation for the particles "
                            f"the detectors see; G_k from the setup's α_k and Γ_k at the flight time, G̃_k averaged "
                            f"over the decays in flight for the fitted τ). The correction is good to about 2% of τ "
                            f"(⁹⁸Ru's simulated runs: ring-averaged a_k, one flight time per distance); that is not "
                            f"in the budget.")
    kept = [f"{rings[r][0]:g}°: {k:.3f}" for r, k in ks.items() if k is not None]
    how = {"model": "the Lorentz boost of the solid angle times the full-energy efficiency at the shifted energy "
                    "over that at E₀, at each ring's mean speed",
           "fit": "fitted from I_u + I_s/k being the same at every distance"}.get(EFFICIENCY_RATIO)
    if kept:
        out["notes"].append(f"Efficiency of the shifted over the unshifted component (k, {how}): "
                            + ", ".join(kept) + ".")
    if how and any(k is None for k in ks.values()):
        out["notes"].append("No efficiency ratio at " + ", ".join(f"{rings[r][0]:g}°" for r, k in ks.items()
                                                                 if k is None) + ": k = 1 there.")
    out["beta"] = b_bar
    u, s = acc["u"], acc["s"]
    tot = u + s
    s_r = np.sqrt(s * s * np.maximum(acc["vu"], 1.0) + u * u * np.maximum(acc["vs"], 1.0)) / tot ** 2
    if method == "decay_curve":
        out["notes"].append(f"Decay curve of the level alone (one exponential per speed, over the spread of speeds "
                            f"the shifted line shows at the largest distance; mean β = {100 * b_bar:.3f}%; a free "
                            "normalisation): its feeding is not modelled; the DDCM corrects for it.")
        for k, d in enumerate(dists):
            out["points"].append({"distance_um": float(d), "unshifted": float(u[k] / tot[k]),
                                  "shifted": float(s[k] / tot[k]), "error": float(s_r[k]), "tau_ps": None,
                                  "tau_err_ps": None})
        out["curve"] = {"x": [row["distance_um"] for row in res.per_distance],
                        "y": [row["model"] for row in res.per_distance],
                        "kind": "unshifted share (decay-curve fit)"}
        out["method"] = "decay_curve"
        return res
    out["notes"].append(f"DDCM with the shifted counts weighted by their own speed (mean β = {100 * b_bar:.3f}%), "
                        "so a spread of the recoils' speeds does not bias τ.")
    if feeder_ok:
        out["notes"].append(f"Feeding corrected with {subject.names[feeder[0]]} → {subject.names[i]} "
                            f"({feeder[1]:.1f} keV), its intensities corrected for the efficiency.")
    elif feeder is not None:
        out["notes"].append(f"The feeder {subject.names[feeder[0]]} → {subject.names[i]} could not be measured: no "
                            "feeding correction.")
    for k, d in enumerate(dists):
        row = res.per_distance[k]
        out["points"].append({"distance_um": float(d), "unshifted": float(u[k] / tot[k]),
                              "shifted": float(s[k] / tot[k]), "error": float(s_r[k]),
                              "tau_ps": row["tau_ps"], "tau_err_ps": row["error_ps"], "used": row["used"]})
    out["curve"] = {"x": [row["distance_um"] for row in res.per_distance],
                    "y": [row["curve"] for row in res.per_distance],
                    "kind": "speed-weighted shifted intensity (DDCM fit)"}
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
