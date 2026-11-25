"""Validation of the nuclear planner against established tools and published data.

Every capability in the physics register (``docs/physics-register/``) has checks here: at least one against a tool
(LISE++, SRIM, ...) and one against literature, data or an exact analytic result. Each check compares physim with a
reference and reports whether the deviation is inside the tolerance the register states::

    from physim.nuclear import validation

    for r in validation.run():                    # every check; reference files are found automatically
        print(r.capability, r.kind, r.status, r.worst, r.tolerance)
    print(validation.report())                    # the same as a Markdown table

Tool references are files in ``tests/reference/nuclear/`` (see its README for the format and provenance header).
They come from closed programs run by hand, so a check whose file does not exist yet reports ``"pending"`` and the
register shows 🟡. The automated tests (``tests/python/test_nuclear_validation.py``) require every check the register
marks ✅ to pass, and deliberately broken formulas to fail.
"""

from __future__ import annotations

import copy
import csv
import math
import os
import re
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from . import data

#: The register's capability names (rows of the nuclear planner table).
CAPABILITIES = (
    "Atomic masses and material data",
    "Two-body reaction kinematics",
    "Stopping power and range",
    "Energy and angular straggling",
    "Rutherford cross section",
    "Distance of closest approach and validity checks",
    "Coulomb trajectories",
    "Detector solid angles and response",
    "Count rates and beam time",
    "Monte Carlo spectra",
    "Coulomb excitation, first-order semiclassical",
    "γ-ray Doppler shift and broadening",
    "Fusion-evaporation residues from a channel table",
    "Fusion cross section (Bass barrier, Wong)",
    "Evaporation Monte Carlo (Hauser–Feshbach)",
    "Neutron detectors and channel selection",
    "Shell-model calculations (KSHELL)",
    "Lifetimes: recoil distance (RDDS)",
    "Lifetimes: Doppler-shift attenuation (DSAM)",
)


@dataclass
class Result:
    """The outcome of one check."""

    capability: str
    #: "tool" (LISE++, SRIM, ...) or "literature" (published data or exact analytic results).
    kind: str
    name: str
    #: What physim is compared with, including where the numbers come from.
    reference: str
    #: The stated tolerance, in words.
    tolerance: str
    #: "pass", "fail" or "pending" (the reference has not been produced yet).
    status: str
    #: Worst deviation found, in the units the tolerance is stated in (often a relative deviation).
    worst: float = math.nan
    points: int = 0
    detail: str = ""
    #: Data for a comparison plot: {"x": ..., "physim": ..., "reference": ..., "xlabel": ..., "ylabel": ...}.
    plot: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.status == "pass"


@dataclass
class Check:
    capability: str
    kind: str
    name: str
    reference: str
    tolerance: str
    #: Called with the reference folder (or None); returns (passed or None if pending, worst, points, detail, plot).
    run: Callable


CHECKS: list = []


def check(capability: str, kind: str, reference: str, tolerance: str):
    """Register a check function (decorator)."""
    if capability not in CAPABILITIES or kind not in ("tool", "literature"):
        raise ValueError(f"unknown capability or kind: {capability!r}, {kind!r}")

    def wrap(fn):
        CHECKS.append(Check(capability, kind, fn.__name__, reference, tolerance, fn))
        return fn

    return wrap


# ---------------------------------------------------------------------------------------------------------------
# Reference files


def reference_dir() -> Optional[Path]:
    """The folder of tool reference files: ``$PHYSIM_NUCLEAR_REFERENCE``, or ``tests/reference/nuclear`` in the
    source tree this package lives in or the current directory is inside; ``None`` if there is none (an installed
    wheel)."""
    env = os.environ.get("PHYSIM_NUCLEAR_REFERENCE")
    if env:
        return Path(env)
    for start in (Path(__file__).resolve(), Path.cwd().resolve()):
        for parent in [start, *start.parents]:
            cand = parent / "tests" / "reference" / "nuclear"
            if cand.is_dir():
                return cand
    return None


def read_reference(path) -> tuple:
    """(header, rows) of a reference file: the ``# key: value`` lines (continuation lines appended) and the CSV
    rows as dicts."""
    header: dict = {}
    last = None
    rows_text = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                body = line[1:].strip()
                key, sep, value = body.partition(":")
                if sep and key and " " not in key.strip().replace(" by", "") and not line.startswith("#  "):
                    last = key.strip()
                    header[last] = value.strip()
                elif last is not None:
                    header[last] += " " + body
            else:
                rows_text.append(line)
    return header, list(csv.DictReader(rows_text))


def _file(ref: Optional[Path], name: str) -> Optional[Path]:
    if ref is None:
        return None
    p = Path(ref) / name
    return p if p.exists() else None


def _rel(a, b) -> np.ndarray:
    return np.abs(np.asarray(a, dtype=float) / np.asarray(b, dtype=float) - 1)


# ---------------------------------------------------------------------------------------------------------------
# Masses and material data

#: Q-values as tabulated in nuclear data compilations, each quoted to the precision given (MeV).
PUBLISHED_Q = [
    (["2H", "3H"], ["4He", "n"], 17.589),  # D-T fusion
    (["2H", "2H"], ["3He", "n"], 3.269),
    (["2H", "2H"], ["3H", "p"], 4.033),
    (["2H", "3He"], ["4He", "p"], 18.353),
    (["n", "6Li"], ["3H", "4He"], 4.783),
    (["n", "p"], ["2H"], 2.224566),  # deuteron binding energy
    (["238U"], ["234Th", "4He"], 4.270),  # alpha decay
    (["12C", "4He"], ["16O"], 7.162),  # 12C(alpha, gamma)16O
]


@check("Atomic masses and material data", "tool", "LISE++ isotope_mass, 47 nuclides from 1H to 238U "
       "(masses_lise.csv); masses AME2020 only estimates (78Ni) are skipped", "2 keV, or the AME2020 uncertainty "
       "if larger (100Sn: 240 keV)")
def masses_vs_lise(ref):
    path = _file(ref, "masses_lise.csv")
    if path is None:
        return None, math.nan, 0, "masses_lise.csv not found", {}
    _, rows = read_reference(path)
    a, diff_kev, allowed, skipped = [], [], [], []
    for r in rows:
        n = data.nuclide((int(r["Z"]), int(r["A"])))
        if n.estimated:
            skipped.append(n.name)
            continue
        a.append(n.A)
        diff_kev.append((n.nuclear_mass_mev - float(r["nuclear_mass_u"]) * data.U_KEV * 1e-3) * 1e3)
        allowed.append(max(2.0, n.mass_excess_unc_kev))
    diff_kev, allowed = np.array(diff_kev), np.array(allowed)
    ok = bool(np.all(np.abs(diff_kev) <= allowed))
    worst = float(np.max(np.abs(diff_kev) / allowed))  # as a fraction of the allowed deviation
    plot = {"x": a, "physim": diff_kev, "reference": np.zeros_like(diff_kev), "xlabel": "mass number A",
            "ylabel": "physim − LISE++ (keV)", "residual": True}
    detail = (f"{len(a)} measured masses, worst {worst:.0%} of the allowed deviation; skipped estimates: "
              f"{', '.join(skipped) or 'none'}")
    return ok, worst, len(a), detail, plot


@check("Atomic masses and material data", "literature", "Published Q-values (D-T fusion, D-D, ³He, ⁶Li(n,t), "
       "deuteron binding, ²³⁸U α decay, ¹²C(α,γ))", "1 keV")
def q_values(ref):
    ours = [data.q_value(i, o) for i, o, _ in PUBLISHED_Q]
    pub = [q for _, _, q in PUBLISHED_Q]
    diff = np.abs(np.array(ours) - np.array(pub)) * 1e3
    worst = float(diff.max())
    plot = {"x": pub, "physim": ours, "reference": pub, "xlabel": "published Q (MeV)", "ylabel": "Q (MeV)"}
    return worst <= 1.0, worst, len(pub), f"worst {worst:.3f} keV", plot


# ---------------------------------------------------------------------------------------------------------------
# Kinematics


@check("Two-body reaction kinematics", "tool", "LISE++ kinematic calculator (kinematics_*.csv; more cases requested "
       "in pending/kinematics_lisepp.md)", "1 keV, or the precision LISE++ displays")
def kinematics_vs_lise(ref):
    from .kinematics import TwoBody

    files = sorted(Path(ref).glob("kinematics_*.csv")) if ref is not None else []
    if not files:
        return None, math.nan, 0, "no kinematics reference files", {}
    worst, n, ok = 0.0, 0, True
    xs, ours, theirs = [], [], []
    for path in files:
        for row in read_reference(path)[1]:
            r = TwoBody(row["beam"], row["target"], float(row["beam_energy_mev"]), ejectile=row["ejectile"] or None,
                        excitation_mev=float(row["excitation_mev"] or 0))
            first, second = r.at_lab(float(row["theta_lab_deg"]), row["particle"])
            point = first if row["branch"] in ("", "1") else second
            tol = max(1e-3, float(row.get("energy_tol_mev") or 0))
            dev = abs(point.energy - float(row["energy_mev"]))
            ok &= dev <= tol
            if row.get("theta_cm_deg"):
                ok &= abs(point.theta_cm - float(row["theta_cm_deg"])) <= 0.01
            worst = max(worst, dev * 1e3)
            n += 1
            xs.append(float(row["theta_lab_deg"]))
            ours.append(point.energy)
            theirs.append(float(row["energy_mev"]))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "lab angle (deg)", "ylabel": "energy (MeV)"}
    return bool(ok), worst, n, f"worst {worst:.2f} keV over {n} points", plot


@check("Two-body reaction kinematics", "literature", "Non-relativistic kinematic factor K(θ) (textbook formula) at "
       "1 keV, four mass ratios", "1e-6 in K")
def kinematics_classical_limit(ref):
    from .kinematics import TwoBody, kinematic_factor

    worst, n = 0.0, 0
    plot = {}
    for beam, target in (("4He", "197Au"), ("1H", "12C"), ("12C", "1H"), ("16O", "4He")):
        r = TwoBody(beam, target, 0.001)
        th = np.array([1.0, 10.0, 30.0, 60.0, 90.0, 120.0, 150.0, 179.0])
        th = th[th < r.max_angle() - 1e-6]
        first, _ = r.at_lab(th)
        k = kinematic_factor(r.m1, r.m2, th)
        worst = max(worst, float(np.max(np.abs(first.energy / 0.001 - k))))
        n += len(th)
        if beam == "4He":
            plot = {"x": th, "physim": first.energy / 0.001, "reference": k, "xlabel": "lab angle (deg)",
                    "ylabel": "E / E_beam (α on Au)"}
    return worst <= 1e-6, worst, n, f"worst {worst:.1e}", plot


# ---------------------------------------------------------------------------------------------------------------
# Stopping and straggling

LISE_SYMBOLS = {6: "C", 13: "Al", 14: "Si", 28: "Ni", 29: "Cu", 47: "Ag", 73: "Ta", 79: "Au"}

#: Worst deviation from LISE++ (ATIMA 1.4) allowed, by quantity, ion class and energy band (MeV/u). Each is the
#: measured worst case over 8 targets (2 of them interpolated) rounded up; see docs/physics-register/stopping.md.
STOPPING_TOLERANCE = {
    ("stopping", "p,a", "1-10"): 0.07, ("stopping", "p,a", ">=10"): 0.04, ("stopping", "p,a", "0.1-1"): 0.13,
    ("stopping", "C-Ar", ">=10"): 0.04, ("stopping", "C-Ar", "1-10"): 0.11, ("stopping", "C-Ar", "0.1-1"): 0.25,
    ("stopping", "Kr,Xe", ">=10"): 0.10, ("stopping", "Kr,Xe", "1-10"): 0.12, ("stopping", "Kr,Xe", "0.1-1"): 0.55,
    ("range", "p,a", ">=10"): 0.035, ("range", "p,a", "1-10"): 0.06, ("range", "p,a", "0.1-1"): 0.21,
    ("range", "C-Ar", ">=10"): 0.05, ("range", "C-Ar", "1-10"): 0.15, ("range", "C-Ar", "0.1-1"): 0.39,
    ("range", "Kr,Xe", ">=10"): 0.09, ("range", "Kr,Xe", "1-10"): 0.53, ("range", "Kr,Xe", "0.1-1"): 1.0,
    ("energy_after", "p,a", ">=10"): 0.005, ("energy_after", "p,a", "1-10"): 0.012,
    ("energy_after", "p,a", "0.1-1"): 0.053, ("energy_after", "C-Ar", ">=10"): 0.008,
    ("energy_after", "C-Ar", "1-10"): 0.048, ("energy_after", "C-Ar", "0.1-1"): 0.115,
    ("energy_after", "Kr,Xe", ">=10"): 0.022, ("energy_after", "Kr,Xe", "1-10"): 0.06,
    ("energy_after", "Kr,Xe", "0.1-1"): 0.19,
}


def lise_stopping_rows(ref) -> list:
    """The ATIMA 1.4 rows of ``stopping_lise.csv`` (empty if the file is missing)."""
    path = _file(ref, "stopping_lise.csv")
    if path is None:
        return []
    return [r for r in read_reference(path)[1] if r["model"] == "ATIMA 1.4"]


def _lise_ion(z: int, a: int) -> str:
    return data.nuclide((z, a)).name if z > 1 else "1H"


@check("Stopping power and range", "tool", "LISE++ ATIMA 1.4: 7 ions × 8 elemental targets × 0.1–300 MeV/u "
       "(stopping_lise.csv)", "by ion class and energy band, 0.5%–100% (STOPPING_TOLERANCE)")
def stopping_vs_lise(ref):
    from .stopping import Stopping

    rows = lise_stopping_rows(ref)
    if not rows:
        return None, math.nan, 0, "stopping_lise.csv not found", {}
    worst: dict = {}
    cache: dict = {}
    xs, ratio = [], []
    for r in rows:
        z, a, zt, e_u = int(r["Z"]), int(r["A"]), int(r["Zt"]), float(r["e_mev_u"])
        if e_u < 0.1:
            continue  # below 0.1 MeV/u the codes themselves differ by tens of percent; not tested
        ion = _lise_ion(z, a)
        s = cache.setdefault((ion, zt), Stopping(ion, LISE_SYMBOLS[zt]))
        cls = "p,a" if z <= 2 else ("C-Ar" if z <= 18 else "Kr,Xe")
        band = "0.1-1" if e_u < 1 else ("1-10" if e_u < 10 else ">=10")
        e = e_u * a
        got = {
            "stopping": s.stopping_power(e) / float(r["stopping_mev_mg_cm2"]),
            "range": s.range(e) / float(r["range_mg_cm2"]),
            "energy_after": s.energy_after(e, float(r["thickness_mg_cm2"])) / a / float(r["e_after_mev_u"]),
        }
        for q, rr in got.items():
            key = (q, cls, band)
            worst[key] = max(worst.get(key, 0.0), abs(rr - 1))
        if z <= 2:
            xs.append(e_u)
            ratio.append(got["stopping"])
    excess = {k: v / STOPPING_TOLERANCE[k] for k, v in worst.items()}
    failures = {k: round(v, 4) for k, v in worst.items() if v > STOPPING_TOLERANCE[k]}
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(ratio)), "xlabel": "energy (MeV/u), p and α",
            "ylabel": "stopping power, physim / LISE++", "logx": True}
    worst_frac = max(excess.values())
    detail = "all bands inside their tolerance" if not failures else f"outside: {failures}"
    return not failures, worst_frac, len(rows), detail + f" (worst {worst_frac:.0%} of a band's tolerance)", plot


def _ion_class(z: int) -> str:
    return "p,a" if z <= 2 else ("C-Ar" if z <= 18 else "Kr,Xe")


def _energy_band(e_u: float) -> str:
    return "0.1-1" if e_u < 1 else ("1-10" if e_u < 10 else ">=10")


#: The LISE++ bands, except Kr and Xe at 10 MeV/u and above: 11% against SRIM, which is 7% above measured stopping of
#: I and Xe in light targets at 5–10 MeV/u (stopping_measured_heavy.csv), where physim is 5% below it.
SRIM_STOPPING_TOLERANCE = {**STOPPING_TOLERANCE, ("stopping", "Kr,Xe", ">=10"): 0.11}


def _stopping_vs_srim(ref, targets: tuple) -> tuple:
    """Total stopping power (electronic + nuclear) against SR Module, judged by the LISE++ bands
    (SRIM_STOPPING_TOLERANCE). SRIM's ranges are projected ranges, which are shorter than physim's path-length
    (CSDA) ranges where ions scatter (light ions in heavy targets, slow ions), so they are reported, not judged."""
    from .stopping import Stopping

    tolerance = SRIM_STOPPING_TOLERANCE

    path = _file(ref, "stopping_srim.csv")
    if path is None:
        return None, math.nan, 0, "stopping_srim.csv not yet produced (SR Module runs requested)", {}
    rows = [r for r in read_reference(path)[1] if r["target"] in targets]
    worst: dict = {}
    per_target: dict = {}
    range_10: dict = {}
    xs, ratio = [], []
    for r in rows:
        nuc = data.nuclide(r["ion"])
        e_u, rho = float(r["e_mev_u"]), float(r["density_g_cm3"])
        s = Stopping(r["ion"], r["target"])
        e = e_u * nuc.A
        d = s.stopping_power(e) / (float(r["electronic_mev_mg_cm2"]) + float(r["nuclear_mev_mg_cm2"])) - 1
        key = ("stopping", _ion_class(nuc.Z), _energy_band(e_u))
        worst[key] = max(worst.get(key, 0.0), abs(d))
        t = r["target"]
        if abs(d) / tolerance[key] > per_target.get(t, (0.0,))[0]:
            per_target[t] = (abs(d) / tolerance[key], f"{r['ion']} {e_u:g} MeV/u {d:+.1%}")
        if e_u >= 10:
            dr = float(s.range(e)) / rho * 10 / float(r["projected_range_um"]) - 1
            range_10[t] = max(range_10.get(t, 0.0), abs(dr))
        xs.append(e_u)
        ratio.append(1 + d)
    failures = {k: round(v, 3) for k, v in worst.items() if v > tolerance[k]}
    worst_frac = max(v / tolerance[k] for k, v in worst.items())
    detail = ("all bands inside their tolerance" if not failures else f"outside: {failures}") + (
        f" (worst {worst_frac:.0%} of a band's tolerance); worst per target: "
        + ", ".join(f"{t} {w[1]}" for t, w in per_target.items())
        + "; path-length vs projected range at 10 MeV/u (not judged): "
        + ", ".join(f"{t} {v:.1%}" for t, v in range_10.items()))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(ratio)), "xlabel": "energy (MeV/u)",
            "ylabel": "stopping power, physim / SRIM", "logx": True}
    return not failures, worst_frac, len(rows), detail, plot


@check("Stopping power and range", "tool", "SRIM 2013 SR Module: 5 ions (⁴He to ¹³²Xe) × Au, Si × 0.1–10 MeV/u "
       "(stopping_srim.csv)", "the LISE++ bands for stopping power (STOPPING_TOLERANCE; Kr, Xe ≥ 10 MeV/u: 11%, SRIM "
       "being above measured data there); ranges reported, not judged")
def stopping_vs_srim(ref):
    return _stopping_vs_srim(ref, ("Au", "Si"))


@check("Stopping power and range", "tool", "SRIM 2013 SR Module, compounds: 5 ions × Mylar, CD₂, CsI × "
       "0.1–10 MeV/u, compound correction 1 (stopping_srim.csv)", "the same bands as for elements "
       "(STOPPING_TOLERANCE); ranges reported, not judged")
def stopping_compounds_vs_srim(ref):
    return _stopping_vs_srim(ref, ("Mylar", "CD2", "CsI"))


#: physim material for each target of stopping_measured_hydrogen.csv, and the group it is judged in.
MEASURED_HYDROGEN_TARGETS = {"H2": ("H", "H₂/D₂ gas"), "D2": ("D2", "H₂/D₂ gas"), "Mylar": ("Mylar", "Mylar"),
                             "Polypropylene": ("Polyethylene", "(CH₂)ₙ"), "Polyethylene": ("Polyethylene", "(CH₂)ₙ"),
                             "CH2": ("Polyethylene", "(CH₂)ₙ")}


@check("Stopping power and range", "literature", "Measured electronic stopping of C, O, Ar, Kr, Xe in H₂, D₂, "
       "(CH₂)ₙ polymers and Mylar, 0.1–15 MeV/u, from the IAEA Stopping Power Database "
       "(stopping_measured_hydrogen.csv; one outlier data set excluded, see its header)",
       "median deviation per target group, ion and energy band inside the LISE++ band for that ion class and energy "
       "band (STOPPING_TOLERANCE)")
def stopping_hydrogen_vs_data(ref):
    from .stopping import Stopping

    path = _file(ref, "stopping_measured_hydrogen.csv")
    if path is None:
        return None, math.nan, 0, "stopping_measured_hydrogen.csv not found", {}
    rows = [r for r in read_reference(path)[1] if r["used"] == "yes" and float(r["e_mev_u"]) >= 0.1]
    cache: dict = {}
    cells: dict = {}
    srim_cells: dict = {}
    xs, ratio = [], []
    for r in rows:
        material, group = MEASURED_HYDROGEN_TARGETS[r["target"]]
        ion = f"{r['isotope']}{r['ion']}"
        z, e_u = data.nuclide(ion).Z, float(r["e_mev_u"])
        s = cache.setdefault((ion, material), Stopping(ion, material))
        meas = float(r["electronic_mev_mg_cm2"])
        d = float(s.stopping_power(e_u * int(r["isotope"]), part="electronic")) / meas - 1
        key = (group, r["ion"], _ion_class(z), _energy_band(e_u))
        cells.setdefault(key, []).append(d)
        srim_cells.setdefault(key, []).append(float(r["srim_electronic_mev_mg_cm2"]) / meas - 1)
        xs.append(e_u)
        ratio.append(1 + d)
    frac, fails, srim_fails = {}, {}, []
    for key, devs in cells.items():
        tol = STOPPING_TOLERANCE[("stopping", key[2], key[3])]
        med = float(np.median(devs))
        frac[key] = abs(med) / tol
        if abs(med) > tol:
            fails[key] = med
        if abs(float(np.median(srim_cells[key]))) > tol:
            srim_fails.append(key)
    worst = max(frac.values())
    label = lambda k: f"{k[1]} in {k[0]} {k[3]} MeV/u"  # noqa: E731
    detail = (f"{len(cells)} cells; " + ("all inside their band" if not fails else "outside: " + ", ".join(
        f"{label(k)} median {v:+.1%} (band {STOPPING_TOLERANCE[('stopping', k[2], k[3])]:.0%}, "
        f"{len(cells[k])} points)" for k, v in sorted(fails.items(), key=lambda kv: kv[1])))
        + f"; worst {worst:.0%} of a band. SRIM 2013 by the same rule (context): "
        + ("all inside" if not srim_fails else "outside in " + ", ".join(label(k) for k in srim_fails)))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(ratio)), "xlabel": "energy (MeV/u)",
            "ylabel": "electronic stopping, physim / measured", "logx": True}
    return not fails, worst, len(rows), detail, plot


#: Bands of the heavy-ion check: median deviation from measured stopping allowed per ion class and energy band.
HEAVY_DATA_TOLERANCE = {"<0.1": 0.25, "0.1-1": 0.12, "1-10": 0.08}


def _heavy_class(z: int) -> str:
    return "Ar-Ni" if z <= 28 else ("Cu-Ag" if z <= 47 else "I-U")


def _heavy_band(e_u: float) -> str:
    return "<0.1" if e_u < 0.1 else ("0.1-1" if e_u < 1 else "1-10")


@check("Stopping power and range", "literature", "Measured electronic stopping of ions from Ar to U (21 elements) in 12 "
       "elemental solids from C to Pb, 5 keV/u–10 MeV/u, from the IAEA Stopping Power Database "
       "(stopping_measured_heavy.csv): slow heavy recoils in targets, backings and plunger stoppers",
       "median deviation per ion class (Ar–Ni, Cu–Ag, I–U) and energy band within 25% below 0.1 MeV/u, 12% at "
       "0.1–1 MeV/u, 8% at 1–10 MeV/u (HEAVY_DATA_TOLERANCE)")
def stopping_heavy_vs_data(ref):
    from .stopping import Stopping

    path = _file(ref, "stopping_measured_heavy.csv")
    if path is None:
        return None, math.nan, 0, "stopping_measured_heavy.csv not found", {}
    rows = read_reference(path)[1]
    cache: dict = {}
    cells: dict = {}
    srim_cells: dict = {}
    by_target: dict = {}
    xs, ratio = [], []
    for r in rows:
        z, a, e_u = int(r["z1"]), int(r["a1"]), float(r["e_mev_u"])
        s = cache.setdefault((z, a, r["target"]), Stopping(data.nuclide((z, a)).name, r["target"]))
        meas = float(r["stopping_mev_mg_cm2"])
        d = float(s.stopping_power(e_u * s.mass_u, part="electronic")) / meas - 1
        key = (_heavy_class(z), _heavy_band(e_u))
        cells.setdefault(key, []).append(d)
        if r["srim_mev_mg_cm2"]:
            srim_cells.setdefault(key, []).append(float(r["srim_mev_mg_cm2"]) / meas - 1)
        group = "light" if int(r["z2"]) <= 14 else ("mid" if int(r["z2"]) <= 32 else "heavy")
        by_target.setdefault((group, key[1]), []).append(d)
        xs.append(e_u)
        ratio.append(1 + d)
    med = {k: float(np.median(v)) for k, v in cells.items()}
    frac = {k: abs(m) / HEAVY_DATA_TOLERANCE[k[1]] for k, m in med.items()}
    fails = {k: m for k, m in med.items() if abs(m) > HEAVY_DATA_TOLERANCE[k[1]]}
    cell = lambda k: f"{k[0]} {k[1]} MeV/u {med[k]:+.1%} ({len(cells[k])})"  # noqa: E731
    detail = (("all cells inside their band" if not fails else "outside: " + ", ".join(cell(k) for k in fails))
              + f" (worst {max(frac.values()):.0%} of a band); medians: " + ", ".join(cell(k) for k in sorted(med))
              + "; by target (C–Si light, Ti–Ge mid, Ag–Pb heavy): " + ", ".join(
                  f"{g} {b} {float(np.median(v)):+.1%}" for (g, b), v in sorted(by_target.items()))
              + "; SRIM 2013 medians (context): " + ", ".join(
                  f"{k[0]} {k[1]} {float(np.median(v)):+.1%}" for k, v in sorted(srim_cells.items())))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(ratio)), "xlabel": "energy (MeV/u)",
            "ylabel": "electronic stopping, physim / measured", "logx": True}
    return not fails, max(frac.values()), len(rows), detail, plot


@check("Stopping power and range", "literature", "NIST PSTAR/ASTAR CSDA ranges (ICRU 49), p and α in C, Si, Au, "
       "water, 2 keV–100 MeV", "1%")
def ranges_vs_nist(ref):
    from . import stopping as st

    worst, n = 0.0, 0
    plot = {}
    energies = np.array([0.002, 0.01, 0.1, 1.0, 10.0, 100.0])
    for ion, prog in (("1H", "pstar"), ("4He", "astar")):
        for material in ("C", "Si", "Au", "Water"):
            key = st._star_compound(material) or st._elements_in_tables()[data.element(material).Z]
            e, _, _, r = st._table(prog, key)
            s = st.Stopping(ion, material)
            nist = np.exp(np.interp(np.log(energies), np.log(e), np.log(r))) * 1e3
            ours = s.range(energies)
            worst = max(worst, float(np.max(_rel(ours, nist))))
            n += len(energies)
            if ion == "4He" and material == "Au":
                plot = {"x": energies, "physim": ours, "reference": nist, "xlabel": "energy (MeV)",
                        "ylabel": "CSDA range of α in Au (mg/cm²)", "logx": True, "logy": True}
    return worst <= 0.01, worst, n, f"worst {worst:.2%}", plot


@check("Energy and angular straggling", "tool", "LISE++ straggling_energy (FWHM) ≥ 30 MeV/u in C and Au "
       "(stopping_lise.csv)", "25%")
def straggling_vs_lise(ref):
    from .stopping import Stopping

    rows = [r for r in lise_stopping_rows(ref) if float(r["e_mev_u"]) >= 30 and int(r["Zt"]) in (6, 79)]
    if not rows:
        return None, math.nan, 0, "stopping_lise.csv not found", {}
    xs, ours, theirs = [], [], []
    for r in rows:
        z, a, zt, e_u = int(r["Z"]), int(r["A"]), int(r["Zt"]), float(r["e_mev_u"])
        fwhm = 2.3548 * Stopping(_lise_ion(z, a), LISE_SYMBOLS[zt]).straggling(
            e_u * a, float(r["thickness_mg_cm2"]), steps=40) / a
        xs.append(e_u)
        ours.append(fwhm)
        theirs.append(float(r["straggling_mev_u"]))
    worst = float(np.max(_rel(ours, theirs)))
    plot = {"x": xs, "physim": np.array(ours) / np.array(theirs), "reference": np.ones(len(xs)),
            "xlabel": "energy (MeV/u)", "ylabel": "straggling FWHM, physim / LISE++", "logx": True}
    return worst <= 0.25, worst, len(rows), f"worst {worst:.1%}", plot


@check("Energy and angular straggling", "literature", "Bohr's formula for a thin layer, 100 MeV protons in Si, and "
       "the textbook constant 4π r_e² (m_e c²)² N_A = 0.1569 MeV² cm²/g", "0.2%")
def straggling_bohr(ref):
    from . import stopping as st

    s = st.Stopping("1H", "Si")
    t = 1.0
    worst, xs, ours, theirs = 0.0, [], [], []
    for e in (50.0, 100.0, 300.0):
        beta2 = 1 - (1 / (1 + e / 938.272)) ** 2
        z_over_a = 14 / data.element("Si").molar_mass
        bohr = math.sqrt(0.1569 * z_over_a * t * 1e-3 * (1 - beta2 / 2) / (1 - beta2))
        got = s.straggling(e, t)
        worst = max(worst, abs(got / bohr - 1))
        xs.append(e)
        ours.append(got)
        theirs.append(bohr)
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "proton energy (MeV)",
            "ylabel": "σ after 1 mg/cm² Si (MeV)"}
    return worst <= 2e-3, worst, len(xs), f"worst {worst:.2%}", plot


# ---------------------------------------------------------------------------------------------------------------
# Rutherford scattering

#: Geiger and Marsden, Phil. Mag. 25, 604 (1913), Table II, gold, first series: lab angle (deg) and number of
#: scintillations N in equal times. Transcribed from https://www.chemteam.info/Chem-History/GeigerMarsden-1913/
#: GeigerMarsden-1913.html; every row's N / (1/sin⁴(φ/2)) reproduces the paper's last column.
GEIGER_MARSDEN_GOLD = [(150, 33.1), (135, 43.0), (120, 51.9), (105, 69.5), (75, 211), (60, 477), (45, 1435),
                       (37.5, 3300), (30, 7800), (22.5, 27300), (15, 132000)]


def _rutherford_rows(ref) -> list:
    path = _file(ref, "rutherford_lisepp.csv")
    return read_reference(path)[1] if path is not None else []


def _rutherford_value(row) -> float:
    from .rutherford import Rutherford

    r = Rutherford(row["beam"], row["target"], float(row["beam_energy_mev"]))
    q = row["quantity"]
    th = float(row["theta_cm_deg"]) if row["theta_cm_deg"] else None
    if q == "cross_section_cm":
        return r.cross_section_cm(th)
    if q == "cross_section_lab":
        return r.cross_section_lab(r.kinematics.at_cm(th).theta_lab)[0]
    if q == "grazing_angle":
        return r.grazing_angle()
    if q == "head_on_distance":
        return r.d0
    if q == "impact_parameter":
        return r.impact_parameter(th)
    raise ValueError(f"unknown quantity {q!r}")


def _rutherford_vs_lise(ref, quantities) -> tuple:
    rows = [r for r in _rutherford_rows(ref) if r["quantity"] in quantities]
    if not rows:
        return None, math.nan, 0, "rutherford_lisepp.csv not found", {}
    ok, worst = True, 0.0
    names, ours, theirs = [], [], []
    for row in rows:
        got, want = _rutherford_value(row), float(row["value"])
        dev = abs(got / want - 1)
        ok &= dev <= float(row["rel_tol"]) + 1e-12
        worst = max(worst, dev)
        names.append(row["quantity"])
        ours.append(got)
        theirs.append(want)
    plot = {"x": names, "physim": ours, "reference": theirs, "xlabel": "", "ylabel": "value", "categorical": True}
    return bool(ok), worst, len(rows), f"worst {worst:.2%}", plot


@check("Rutherford cross section", "tool", "LISE++ kinematic calculator, ⁴He + ¹⁹⁷Au at 20 MeV, 50° CM, CM and lab "
       "(rutherford_lisepp.csv)", "0.1%")
def rutherford_vs_lise(ref):
    return _rutherford_vs_lise(ref, ("cross_section_cm", "cross_section_lab"))


@check("Rutherford cross section", "literature", "Geiger and Marsden 1913, Table II (gold, 15°–150°)",
       "counts/theory within ±25% of its mean (their own data scatter ±18%)")
def rutherford_geiger_marsden(ref):
    from .rutherford import Rutherford

    r = Rutherford("4He", "197Au", 7.7)  # RaC' α particles
    angles = np.array([a for a, _ in GEIGER_MARSDEN_GOLD], dtype=float)
    counts = np.array([n for _, n in GEIGER_MARSDEN_GOLD])
    theory = np.array([r.cross_section_lab(a)[0] for a in angles])
    ratio = counts / theory
    dev = np.abs(ratio / ratio.mean() - 1)
    worst = float(dev.max())
    plot = {"x": angles, "physim": theory * ratio.mean(), "reference": counts, "xlabel": "lab angle (deg)",
            "ylabel": "scintillations (theory scaled to the mean)", "logy": True}
    return worst <= 0.25, worst, len(angles), f"worst {worst:.0%}", plot


@check("Distance of closest approach and validity checks", "tool", "LISE++ head-on distance d₀, b(60°) (read off "
       "its plots) and grazing angle, ⁴He + ¹⁹⁷Au at 20 MeV (rutherford_lisepp.csv)", "2% (plot reading)")
def closest_approach_vs_lise(ref):
    return _rutherford_vs_lise(ref, ("grazing_angle", "head_on_distance", "impact_parameter"))


@check("Distance of closest approach and validity checks", "literature", "Onset of non-Rutherford scattering, α + Pb "
       "(Farwell and Wegner 1954; requested in pending/rutherford_above_barrier.md)", "onset energy within 10%")
def closest_approach_above_barrier(ref):
    path = _file(ref, "rutherford_above_barrier.csv")
    if path is None:
        return None, math.nan, 0, "rutherford_above_barrier.csv not yet produced (data requested)", {}
    from .rutherford import Rutherford

    rows = read_reference(path)[1]
    worst = 0.0
    for row in rows:
        lo, hi = 1.0, 200.0
        th = float(row["theta_cm_deg"])
        for _ in range(60):  # the beam energy at which the grazing angle reaches θ
            mid = 0.5 * (lo + hi)
            if Rutherford(row["beam"], row["target"], mid).grazing_angle() > th:
                lo = mid
            else:
                hi = mid
        worst = max(worst, abs(0.5 * (lo + hi) / float(row["onset_mev"]) - 1))
    return worst <= 0.10, worst, len(rows), f"worst {worst:.1%}", {}


@check("Coulomb trajectories", "literature", "Exact deflection of a Coulomb orbit, 2 arctan(d₀/2b), for orbits "
       "integrated by the engine (α + Au, 4 impact parameters)", "1e-6°")
def trajectories_analytic(ref):
    from .rutherford import Rutherford

    r = Rutherford("4He", "197Au", 20.0)
    bs = r.d0 * np.array([0.05, 0.3, 1.0, 3.0])
    got, want = [], []
    for orbit in r.trajectories(bs, distance=500 * r.d0):
        got.append(orbit.deflection())
        want.append(r.angle_for_impact(orbit.b))
    worst = float(np.max(np.abs(np.array(got) - np.array(want))))
    plot = {"x": bs, "physim": got, "reference": want, "xlabel": "impact parameter (fm)",
            "ylabel": "deflection (deg)"}
    return worst <= 1e-6, worst, len(bs), f"worst {worst:.1e}°", plot


# ---------------------------------------------------------------------------------------------------------------
# Detectors, rates, events


@check("Detector solid angles and response", "literature", "Closed forms: disc on axis 2π(1 − cos α), centred "
       "rectangle 4 arctan(ab/(d√(a²+b²+d²))), annulus on axis", "1e-10")
def solid_angles_closed_form(ref):
    from .detectors import Geometry, direction

    worst, xs, ours, theirs = 0.0, [], [], []
    for d, radius in ((100.0, 2.5), (40.0, 30.0), (10.0, 50.0)):
        g = Geometry("disc", "circle", [0, 0, d], [0, 0, -1], radius=radius)
        exact = 2 * math.pi * (1 - d / math.hypot(d, radius)) * 1e3
        xs.append(math.degrees(math.atan(radius / d)))
        ours.append(g.solid_angle())
        theirs.append(exact)
    for d, a, b in ((150.0, 25.0, 25.0), (50.0, 40.0, 10.0)):
        g = Geometry("r", "rectangle", direction(37.0, 20.0) * d, -direction(37.0, 20.0), width=2 * a, height=2 * b)
        exact = 4 * math.atan(a * b / (d * math.sqrt(a * a + b * b + d * d))) * 1e3
        xs.append(math.degrees(math.atan(math.hypot(a, b) / d)))
        ours.append(g.solid_angle())
        theirs.append(exact)
    g = Geometry("cd", "annular", [0, 0, -40], [0, 0, 1], inner_radius=9.0, outer_radius=41.0)
    xs.append(math.degrees(math.atan(41 / 40)))
    ours.append(g.solid_angle())
    theirs.append(2 * math.pi * (40 / math.hypot(40, 9) - 40 / math.hypot(40, 41)) * 1e3)
    worst = float(np.max(_rel(ours, theirs)))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "half-angle of the face (deg)",
            "ylabel": "solid angle (msr)", "logy": True}
    return worst <= 1e-10, worst, len(xs), f"worst {worst:.1e}", plot


@check("Count rates and beam time", "literature", "I n (dσ/dΩ) Ω for the small discs of the α + Au example "
       "(30°, 60°, 135°), cross section at mid-target", "1%")
def rates_small_detector(ref):
    from . import Experiment
    from .rates import MB_CM2, Rates, stack
    from .rutherford import Rutherford
    from .stopping import Stopping

    exp = Experiment.example("alpha_on_gold")
    r = Rates(exp)
    layer = stack(exp)[0]
    n = layer.nuclides[(79, 197)]
    pps = exp.beam.particles_per_second
    e_mid = Stopping("4He", "Au").energy_after(5.5, layer.thickness / 2)
    ruth = Rutherford("4He", "197Au", e_mid)
    xs, ours, theirs = [], [], []
    for name in ("A30", "A60", "A135"):
        g = r.array[name]
        sigma, _ = ruth.cross_section_lab(g.mean_theta())
        theirs.append(pps * n * sigma * g.solid_angle() * 1e-3 * MB_CM2)
        ours.append(r.by_channel(name)[("197Au (target)", "ejectile")])
        xs.append(g.mean_theta())
    worst = float(np.max(_rel(ours, theirs)))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "detector angle (deg)",
            "ylabel": "counts per second", "logy": True}
    return worst <= 0.01, worst, len(xs), f"worst {worst:.2%}", plot


@check("Count rates and beam time", "tool", "LISE++ Rutherford cross sections at the detector angles of the "
       "examples (requested in pending/rates_lisepp.md)", "0.5% per detector")
def rates_vs_lise(ref):
    path = _file(ref, "rates_lisepp.csv")
    if path is None:
        return None, math.nan, 0, "rates_lisepp.csv not yet produced (requested)", {}
    from .rutherford import Rutherford

    rows = read_reference(path)[1]
    ours = [Rutherford(r["beam"], r["target"], float(r["beam_energy_mev"])).cross_section_lab(
        float(r["theta_lab_deg"]))[0] for r in rows]
    theirs = [float(r["cross_section_lab_mb_sr"]) for r in rows]
    worst = float(np.max(_rel(ours, theirs)))
    plot = {"x": [float(r["theta_lab_deg"]) for r in rows], "physim": ours, "reference": theirs,
            "xlabel": "lab angle (deg)", "ylabel": "dσ/dΩ lab (mb/sr)", "logy": True}
    return worst <= 0.005, worst, len(rows), f"worst {worst:.2%}", plot


@check("Monte Carlo spectra", "literature", "Peak mean and width of the α + Au example at 30°, against kinematics plus "
       "mean energy loss and the quadrature sum of the widths (10⁶ events)", "mean 4 standard errors, width 5%")
def spectra_vs_analytic(ref):
    from . import Experiment
    from .events import simulate
    from .rates import Rates

    exp = Experiment.example("alpha_on_gold")
    pk = Rates(exp).peaks("A30")[0]
    ev = simulate(exp, 1_000_000, seed=1)
    m = ev.select("A30", particle="ejectile")
    x, w = ev["measured"][m], ev["weight"][m]
    mean = float(np.average(x, weights=w))
    sd = math.sqrt(float(np.average((x - mean) ** 2, weights=w)))
    se = sd / math.sqrt(m.sum())
    ok = abs(mean - pk.mean) <= max(4 * se, 1e-3) and abs(sd / pk.sigma - 1) <= 0.05
    worst = max(abs(mean - pk.mean) / max(4 * se, 1e-3), abs(sd / pk.sigma - 1) / 0.05)
    h, edges = ev.spectrum("A30", bins=60, particle="ejectile")
    centres = 0.5 * (edges[1:] + edges[:-1])
    gauss = h.sum() * (edges[1] - edges[0]) / (pk.sigma * math.sqrt(2 * math.pi)) * np.exp(
        -0.5 * ((centres - pk.mean) / pk.sigma) ** 2)
    plot = {"x": centres, "physim": h, "reference": gauss, "xlabel": "measured energy (MeV)",
            "ylabel": "counts per bin (Monte Carlo against the analytic peak)"}
    return bool(ok), worst, int(m.sum()), (f"mean {mean:.4f} vs {pk.mean:.4f} MeV, σ {sd * 1e3:.1f} vs "
                                           f"{pk.sigma * 1e3:.1f} keV"), plot


@check("Monte Carlo spectra", "tool", "TRIM (SRIM 2013) transmitted-energy spectra, 4 layers, from physim's own "
       "TRIM.IN exports (straggling_trim.csv)", "mean 1%, width 20%")
def spectra_vs_trim(ref):
    path = _file(ref, "straggling_trim.csv")
    if path is None:
        return None, math.nan, 0, "straggling_trim.csv not yet produced (TRIM runs requested)", {}
    from .stopping import Stopping

    rows = read_reference(path)[1]
    ok, worst, notes = True, 0.0, []
    for r in rows:
        s = Stopping(r["ion"], r["material"], float(r["density"]) if r.get("density") else None)  # g/cm³
        mean = s.energy_after(float(r["energy_mev"]), r["thickness"])
        sigma = s.straggling(float(r["energy_mev"]), r["thickness"])
        d_mean = mean / float(r["mean_mev"]) - 1
        d_sigma = sigma / float(r["sigma_mev"]) - 1
        ok &= abs(d_mean) <= 0.01 and abs(d_sigma) <= 0.20 and r.get("opened_unchanged", "yes") == "yes"
        worst = max(worst, abs(d_mean) / 0.01, abs(d_sigma) / 0.20)
        notes.append(f"{r['ion']} {r['energy_mev']} MeV in {r['thickness']} {r['material']}: mean {d_mean:+.2%}, "
                     f"σ {d_sigma:+.0%}")
    return bool(ok), worst, len(rows), "; ".join(notes), {}


# ---------------------------------------------------------------------------------------------------------------
# Coulomb excitation and γ rays


@check("Coulomb excitation, first-order semiclassical", "literature", "Exact results of the first-order theory: the "
       "orbit integrals at θ = 180°, ξ = 0 (2, 2/3, 4/15 for E1, E2, E3) and the closed-form P(180°) for "
       "¹⁶O + ⁵⁸Ni (2⁺)", "1e-6")
def coulex_closed_form(ref):
    from .coulex import HBARC_MEV_FM, Coulex, orbit_integrals
    from .rutherford import E2_MEV_FM

    worst = 0.0
    for lam, exact in ((1, 2.0), (2, 2 / 3), (3, 4 / 15)):
        worst = max(worst, float(np.max(np.abs(orbit_integrals(lam, 1.0, 0.0) / exact - 1))))
    c = Coulex("16O", "58Ni", 30.0, energy=1e-9, b_up=695.0)
    pref = (4 * math.pi * 8 * E2_MEV_FM / (HBARC_MEV_FM * 5)) ** 2 * c.b_up / 5
    closed = pref * (2 / 3) ** 2 * (5 / (4 * math.pi)) / (c.beta * c.a**2) ** 2
    got = c.probability(180.0, exact=True)
    worst = max(worst, abs(got / closed - 1))
    plot = {"x": ["I(E1)", "I(E2)", "I(E3)", "P(180°)"], "physim": [
        float(orbit_integrals(1, 1.0, 0.0)[0].real), float(orbit_integrals(2, 1.0, 0.0)[0].real),
        float(orbit_integrals(3, 1.0, 0.0)[0].real), got], "reference": [2.0, 2 / 3, 4 / 15, closed],
        "categorical": True}
    return worst <= 1e-6, worst, 4, f"worst {worst:.1e}", plot


@check("Coulomb excitation, first-order semiclassical", "tool", "GOSIA 20110524.14, first order (matrix element ÷10, "
       "weak coupling), E1 polarisation at its default (d = 0.0038) in both, ¹⁶O + ⁵⁸Ni at 30 MeV, five CM angles "
       "(coulex_gosia.csv)", "2% in P(θ)")
def coulex_vs_gosia(ref):
    path = _file(ref, "coulex_gosia.csv")
    if path is None:
        return None, math.nan, 0, "coulex_gosia.csv not found", {}
    from .coulex import GOSIA_DIPOLE_POLARISATION, Coulex

    rows = read_reference(path)[1]
    ours, theirs, xs = [], [], []
    for r in rows:
        c = Coulex(r["beam"], r["target"], float(r["beam_energy_mev"]), energy=float(r["state_mev"]),
                   b_up=f"{r['b_e2_e2b2']} e2b2", e1_polarisation=GOSIA_DIPOLE_POLARISATION)
        ours.append(c.probability(float(r["theta_cm_deg"]), exact=True))
        theirs.append(float(r["probability"]))
        xs.append(float(r["theta_cm_deg"]))
    worst = float(np.max(_rel(ours, theirs)))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "CM angle (deg)",
            "ylabel": "excitation probability", "logy": True}
    return worst <= 0.02, worst, len(rows), f"worst {worst:.2%}", plot


@check("Coulomb excitation, first-order semiclassical", "literature", "Reproduced experiment: Allmond et al., Phys. "
       "Rev. C 90, 034309 (2014), ⁶²Ni at 111.6 MeV on 1 mg/cm² carbon, carbon recoils in two rings (14–28°, "
       "28–44°): the example coulex_ni62_c_allmond against the paper's energy losses, GOSIA 20110524.14 (nuclear "
       "masses and the E1 polarisation off in both; and as run, integer masses, DIP 3.8 and d = 0.0038), its own simulated runs (closure over several "
       "seeds, and one as a stored run) and the paper's per-ring <0⁺‖E2‖2⁺> rebuilt into ratios through GOSIA's "
       "ratio(M) (the paper gives no yields: this tests that physim extracts the same B(E2) from the same data as "
       "GOSIA did); all orders judged, first order shown beside it (physim.nuclear.coulex_reproduction; "
       "coulex/allmond2014.toml)",
       "⁵⁸Ni, ⁶⁰Ni losses within their 0.8 MeV, ⁶²Ni 42.4(8) MeV; σ_C/σ_R per ring within 1% of GOSIA as run "
       "(missed by +1.3/+1.2% before physim had the E1 polarisation) and 0.5% of GOSIA with nuclear masses and DIP off; "
       "closure and the "
       "stored run within a quarter of the paper's statistical error + 3 σ; B(E2) per ring within the paper's 4.1%; "
       "Q_s = ±0.1 e b and ENSDF feeding as in GOSIA within 0.5 percentage points and 2% of its change")
def reproduction_allmond2014(ref):
    """``worst`` is the largest deviation in units of its band. (Until physim had the E1 polarisation the check failed,
    for a known and explained miss of one band, GOSIA as run.)"""
    from . import coulex_reproduction as cr

    r = cr.run(ref)
    if r is None:
        return None, math.nan, 0, f"{cr.REFERENCE_FILE} not found", {}
    c = r["checks"]
    xs, ours, theirs = [], [], []
    for x in c["gosia"]["rows"]:
        xs.append(f"{x['ring']} R vs GOSIA {'matched' if x['setting'] == 'matched' else 'as run'}")
        ours.append(x["physim"])
        theirs.append(x["gosia"])
    for x in c["paper"]["rows"]:
        xs.append(f"{x['ring']} B(E2) vs paper")
        ours.append(x["b_physim_e2b2"])
        theirs.append(x["b_paper_e2b2"])
    for x in c["closure"]["rows"]:
        xs.append(f"closure {x['rings']}")
        ours.append(x["b_e2fm4"])
        theirs.append(x["b_in_e2fm4"])
    points = sum(len(v["rows"]) for v in c.values())
    plot = {"x": xs, "physim": ours, "reference": theirs, "categorical": True}
    for x in c["stored"]["rows"]:
        xs.append(f"stored run {x['what']} {x['rings']}")
        ours.append(x["physim"])
        theirs.append(x["reference"])
    return r["passed"], r["worst"], points, r["detail"], plot


@check("γ-ray Doppler shift and broadening", "literature", "Lorentz boost of the photon four-vector (relativistic "
       "Doppler formula with aberration), 50 random velocities and angles", "1e-12")
def doppler_lorentz(ref):
    from .gamma import doppler_energy

    rng = np.random.default_rng(7)
    beta = rng.uniform(0.0, 0.5, 50)
    cos_rest = rng.uniform(-1.0, 1.0, 50)
    gam = 1 / np.sqrt(1 - beta**2)
    e_lab = gam * (1 + beta * cos_rest)
    cos_lab = (cos_rest + beta) / (1 + beta * cos_rest)
    got = doppler_energy(1.0, beta, cos_lab)
    worst = float(np.max(np.abs(got / e_lab - 1)))
    order = np.argsort(cos_lab)
    plot = {"x": cos_lab[order], "physim": got[order], "reference": e_lab[order], "xlabel": "cos α (lab)",
            "ylabel": "E_γ / E₀"}
    return worst <= 1e-12, worst, 50, f"worst {worst:.1e}", plot


# ---------------------------------------------------------------------------------------------------------------
# Fusion-evaporation: the compound nucleus, fusion, evaporation, neutron detectors

#: Bass 1980 barrier heights V_B (MeV, CM) and positions R_B (fm): I. Dutt and R. K. Puri, Phys. Rev. C 81, 064609
#: (2010), Table II, column "Bass 80" (arXiv:1005.5213).
BASS80_DUTT_PURI = [
    ("12C", "17O", 7.79, 8.13), ("4He", "209Bi", 20.30, 11.00), ("48Ca", "48Ca", 51.39, 10.40),
    ("16O", "144Sm", 61.34, 10.82), ("28Si", "92Zr", 71.21, 10.51), ("16O", "208Pb", 75.92, 11.60),
    ("40Ca", "90Zr", 99.28, 10.71), ("86Kr", "208Pb", 298.65, 13.15),
]

#: σ_fus is compared with PACE4 where physim's CM energy is at least this many times physim's Bass barrier. Nearer
#: the barrier the two transmissions differ by construction (PACE4's sharp cut-off against physim's parabolic
#: tunnelling), and those points are reported, not judged.
FUSION_ABOVE_BARRIER = 1.1
#: Events per system in the evaporation comparisons with PACE4 (PACE4 itself ran 20 000 cascades).
EVAPORATION_EVENTS = 20000


def _pace4_files(ref, kind: str) -> list:
    """PACE4 reference files: "evaporation" (``pace4_<beam>_<target>_<E>MeV.csv``) or "fusion"
    (``pace4_fusion_<beam>_<target>.csv``), with (beam, target) from the name."""
    if ref is None:
        return []
    out = []
    for path in sorted(Path(ref).glob("pace4_*.csv")):
        parts = path.stem.split("_")
        if kind == "fusion" and parts[1] == "fusion":
            out.append((path, parts[2], parts[3]))
        elif kind == "evaporation" and parts[1] != "fusion" and parts[-1].endswith("MeV"):
            out.append((path, parts[1], parts[2]))
    order = {b: i for i, b in enumerate(("4He", "12C", "16O", "32S", "48Ca"))}
    return sorted(out, key=lambda t: (order.get(t[1], 99), t[1]))


_EVAPORATION_CACHE: dict = {}


def physim_evaporation(beam: str, target: str, excitation_mev: float, sigma_j_mb: tuple,
                       events: int = EVAPORATION_EVENTS, seed: int = 1):
    """physim's evaporation of the compound nucleus of ``beam`` + ``target`` at PACE4's E* and σ_J (cached: the
    suite runs these comparisons many times)."""
    from . import evaporation

    key = (beam, target, float(excitation_mev), tuple(sigma_j_mb), int(events), seed)
    if key not in _EVAPORATION_CACHE:
        b, t = data.nuclide(beam), data.nuclide(target)
        _EVAPORATION_CACHE[key] = evaporation.evaporate((b.Z + t.Z, b.A + t.A), float(excitation_mev),
                                                        sigma_l_mb=list(sigma_j_mb), events=events, seed=seed)
    return _EVAPORATION_CACHE[key]


def _pace4_evaporation(path) -> tuple:
    """(header numbers, PACE4 shares by (Z, A)) of one evaporation reference file."""
    header, rows = read_reference(path)
    total = sum(int(r["events"]) for r in rows)
    shares = {(int(r["Z"]), int(r["A"])): int(r["events"]) / total for r in rows}
    info = {"excitation_mev": float(header["excitation_mev"]),
            "sigma_j_mb": tuple(float(x) for x in header["sigma_j_mb"].split(",")),
            "mean": {k: float(v) for k, v in (kv.split("=") for kv in header["mean_cm_energy_mev"].split(","))},
            "multiplicity": {k: float(v) for k, v in (kv.split("=") for kv in header["multiplicity"].split(","))},
            "events": total}
    return info, shares


def _nuclide_label(za) -> str:
    return data.nuclide(za).name


@check("Fusion-evaporation residues from a channel table", "literature", "Hand calculation from AME2020 nuclear "
       "masses: E* = √s − M_CN and β = p_b / (E_b + m_t), ¹²C + ¹⁵⁰Nd at 60 MeV, ¹⁶O + ¹⁴⁸Nd at 85 MeV, "
       "⁴⁸Ca + ¹²⁴Sn at 200 MeV", "1 keV in E*, 1e-6 in β")
def compound_nucleus_hand_calculation(ref):
    from . import residues

    worst, names, ours, theirs = 0.0, [], [], []
    ok = True
    for beam, target, e in (("12C", "150Nd", 60.0), ("16O", "148Nd", 85.0), ("48Ca", "124Sn", 200.0)):
        b, t = data.nuclide(beam), data.nuclide(target)
        mb, mt = b.nuclear_mass_mev, t.nuclear_mass_mev
        mc = data.nuclide((b.Z + t.Z, b.A + t.A)).nuclear_mass_mev
        p = math.sqrt(e**2 + 2 * e * mb)
        total = e + mb + mt
        e_star = math.sqrt(total**2 - p**2) - mc
        beta = p / total
        cn = residues.compound(beam, target, e)
        d_e, d_b = abs(cn.excitation - e_star), abs(cn.beta - beta)
        ok &= d_e <= 1e-3 and d_b <= 1e-6
        worst = max(worst, d_e / 1e-3, d_b / 1e-6)
        names += [f"E* {beam}", f"β {beam}"]
        ours += [cn.excitation, cn.beta]
        theirs += [e_star, beta]
    plot = {"x": names, "physim": ours, "reference": theirs, "categorical": True}
    return bool(ok), worst, 6, f"worst {worst:.1e} of the tolerance", plot


@check("Fusion-evaporation residues from a channel table", "literature", "Residues per channel against the channel "
       "table × luminosity, the fusion_dy158 example (2 × 10⁵ events)", "every channel within 4 standard errors; "
       "the total within 1e-9")
def residue_counts_vs_luminosity(ref):
    from . import residues
    from .experiment import Experiment

    ev = residues.simulate_residues(Experiment.example("fusion_dy158"), events=200_000, seed=11)
    total_dev = abs(ev.rate()[0] / ev.expected() - 1)
    worst, xs, ours, theirs = 0.0, [], [], []
    for label in ev.channels:
        r, err = ev.rate(label)
        worst = max(worst, abs(r - ev.expected(label)) / (4 * err))
        xs.append(label)
        ours.append(r)
        theirs.append(ev.expected(label))
    ok = worst <= 1.0 and total_dev <= 1e-9
    plot = {"x": xs, "physim": ours, "reference": theirs, "categorical": True}
    return bool(ok), worst, len(xs), (f"worst channel {worst * 4:.1f} standard errors; total off by "
                                      f"{total_dev:.1e}"), plot


@check("Fusion cross section (Bass barrier, Wong)", "tool", "PACE4 4.36.0 (Bass model, L-diffuseness 2): excitation "
       "functions of the five systems of theory/evaporation, 8–10 energies each from the barrier to 2 V_B "
       "(pace4_fusion_*.csv)", f"15% where E_cm ≥ {FUSION_ABOVE_BARRIER} V_B (physim's Bass barrier); nearer the "
       "barrier reported, not judged")
def fusion_vs_pace4(ref):
    from .fusion import Fusion

    files = _pace4_files(ref, "fusion")
    if not files:
        return None, math.nan, 0, "pace4_fusion_*.csv not found", {}
    worst, n, notes, near = 0.0, 0, [], []
    xs, ratio, l_notes = [], [], []
    for path, beam, target in files:
        f = Fusion(beam, target)
        vb = f.barrier().height_mev
        sys_worst, sys_l = 0.0, 0.0
        for r in read_reference(path)[1]:
            e = float(r["beam_energy_mev"])
            x = float(f.e_cm(e)) / vb
            theirs = float(r["sigma_fus_mb"])
            ours = float(f.cross_section(e))
            if x < FUSION_ABOVE_BARRIER:
                near.append(f"{beam} {e:g} MeV ({x:.2f} V_B): {ours:.0f} vs {theirs:.0f} mb")
                continue
            d = ours / theirs - 1
            sys_worst = max(sys_worst, abs(d))
            n += 1
            xs.append(x)
            ratio.append(1 + d)
            sj = np.array([float(s) for s in r["sigma_j_mb"].split(";")])
            l_pace = float((np.arange(sj.size) * sj).sum() / sj.sum())
            l_ours = f.compound(e).mean_l
            sys_l = max(sys_l, abs(l_ours / l_pace - 1))
        worst = max(worst, sys_worst)
        notes.append(f"{beam} + {target} {sys_worst:.1%}")
        l_notes.append(f"{beam} {sys_l:.0%}")
    detail = (f"worst above {FUSION_ABOVE_BARRIER} V_B per system: " + ", ".join(notes)
              + "; near the barrier (not judged): " + "; ".join(near)
              + "; mean ℓ, worst per system (not judged): " + ", ".join(l_notes))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(xs)), "xlabel": "E_cm / V_B (physim's Bass barrier)",
            "ylabel": "σ_fus, physim / PACE4"}
    return worst <= 0.15, worst, n, detail, plot


@check("Fusion cross section (Bass barrier, Wong)", "literature", "Bass 1980 barriers V_B and R_B for 8 systems, "
       "Dutt and Puri, Phys. Rev. C 81, 064609 (2010), Table II", "0.03 MeV and 0.03 fm (the table's precision)")
def fusion_barrier_dutt_puri(ref):
    from .fusion import Fusion

    worst, names, ours, theirs = 0.0, [], [], []
    for beam, target, vb, rb in BASS80_DUTT_PURI:
        b = Fusion(beam, target).barrier()
        worst = max(worst, abs(b.height_mev - vb) / 0.03, abs(b.radius_fm - rb) / 0.03)
        names.append(f"{beam}+{target}")
        ours.append(b.height_mev)
        theirs.append(vb)
    plot = {"x": theirs, "physim": ours, "reference": theirs, "xlabel": "V_B, Dutt and Puri (MeV)",
            "ylabel": "V_B (MeV)"}
    return worst <= 1.0, worst, 2 * len(names), f"worst {worst:.0%} of the tolerance", plot


@check("Fusion cross section (Bass barrier, Wong)", "literature", "Exact results of the model: Wong's closed form "
       "against the fixed-barrier partial-wave sum (3 systems, V_B − 6 to V_B + 40 MeV); the sharp cut-off limit "
       "πR_B²(1 − V_B/E) at 1.5–3 V_B; the thick-target average against the integral over the energy loss "
       "(¹⁶O + ²⁰⁸Pb, 90 MeV, 2 mg/cm²)", "1%, 1e-6 and 1%")
def fusion_closed_forms(ref):
    from .fusion import Fusion
    from .stopping import Stopping

    worst = 0.0
    xs, ours, theirs = [], [], []
    for system in (("16O", "208Pb"), ("48Ca", "48Ca"), ("28Si", "92Zr")):
        f = Fusion(*system)
        b = f.barrier()
        for de in (-6.0, -3.0, 0.0, 2.0, 5.0, 15.0, 40.0):
            e = b.height_mev + de
            total, wong = f.partial_waves_cm(e, l_dependent=False).sum(), float(f.wong_cm(e))
            worst = max(worst, abs(total / wong - 1) / 0.01)
            if system == ("16O", "208Pb"):
                xs.append(de)
                ours.append(total)
                theirs.append(wong)
    f = Fusion("16O", "208Pb")
    b = f.barrier()
    for factor in (1.5, 2.0, 3.0):
        e = factor * b.height_mev
        geometric = math.pi * b.radius_fm**2 * (1 - b.height_mev / e) * 10.0
        worst = max(worst, abs(float(f.wong_cm(e)) / geometric - 1) / 1e-6)
    e_in, t = 90.0, 2.0
    tt = f.thick_target(e_in, t)
    stop = Stopping("16O", "208Pb")
    e = np.linspace(tt.e_out_mev, e_in, 2001)
    integrand = f.cross_section(e) / stop.stopping_power(e)
    expected = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(e)) / t
    worst = max(worst, abs(tt.sigma_mb / expected - 1) / 0.01)
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "E_cm − V_B (MeV), ¹⁶O + ²⁰⁸Pb",
            "ylabel": "σ_fus (mb): partial-wave sum against Wong", "logy": True}
    return worst <= 1.0, worst, 25, f"worst {worst:.0%} of its tolerance", plot


#: Measured fusion excitation functions are judged from this many times physim's Bass barrier: just above the
#: barrier for a spherical target; higher for a deformed one, whose rotational states spread the barrier over a
#: distribution several MeV wide (a one-dimensional barrier cannot follow it there).
FUSION_DATA_FROM = {"142Nd": 1.05, "150Nd": 1.15, "154Sm": 1.15}
#: ... up to this many times the barrier. Above it Visakh et al. (2021) find fusion hindered in both ¹⁶O + Nd
#: systems (they need a larger diffuseness than any global potential has, and attribute it to dynamical effects),
#: which a one-dimensional Bass barrier does not describe. Points above are reported, not judged.
FUSION_DATA_TO = 1.3
#: Allowed ``|physim / measured − 1|``: backlog 75's 15% for Stage 2 combined in quadrature with the measurements'
#: total errors (15% for Visakh et al., 10% for Stokstad et al.), rounded up.
FUSION_DATA_TOLERANCE = 0.20


@check("Fusion cross section (Bass barrier, Wong)", "literature", "Measured fusion excitation functions: ¹⁶O + "
       "¹⁴²Nd, ¹⁵⁰Nd (Visakh et al., Phys. Rev. C 104, 054602 (2021), Table I) and ¹⁶O + ¹⁵⁴Sm (Stokstad et al., "
       "Phys. Rev. C 21, 2427 (1980)), from EXFOR (fusion_*_literature.csv)",
       f"{FUSION_DATA_TOLERANCE:.0%} from 1.05 V_B (spherical ¹⁴²Nd) or 1.15 V_B (deformed ¹⁵⁰Nd, ¹⁵⁴Sm) to "
       f"{FUSION_DATA_TO} V_B (physim's Bass barrier); below (coupled channels) and above (fusion hindrance) "
       "reported, not judged")
def fusion_vs_data(ref):
    from .fusion import Fusion

    files = sorted(Path(ref).glob("fusion_*_literature.csv")) if ref is not None else []
    if not files:
        return None, math.nan, 0, "fusion_*_literature.csv not found", {}
    systems: dict = {}
    for path in files:
        for r in read_reference(path)[1]:
            systems.setdefault((r["beam"], r["target"]), []).append(r)
    worst, n, judged, below, above = 0.0, 0, [], [], []
    xs, ratio = [], []
    for (beam, target), rows in systems.items():
        f = Fusion(beam, target)
        vb = f.barrier().height_mev
        lo = FUSION_DATA_FROM.get(target, 1.15)
        sys_worst, sys_n, sys_below, sys_above = 0.0, 0, [], []
        for r in rows:
            e_cm = float(r["e_cm_mev"]) if r.get("e_cm_mev") else float(f.e_cm(float(r["e_lab_mev"])))
            x = e_cm / vb
            ours = float(f.cross_section_cm(e_cm))
            q = ours / float(r["sigma_fus_mb"])
            if x < lo:
                sys_below.append((x, q))
            elif x > FUSION_DATA_TO:
                sys_above.append((x, q))
            else:
                sys_worst = max(sys_worst, abs(q - 1))
                sys_n += 1
                xs.append(x)
                ratio.append(q)
        n += sys_n
        worst = max(worst, sys_worst)
        judged.append(f"{target} {sys_worst:.1%} ({sys_n} points)" if sys_n else
                      f"{target} no point in the window (data end at {max(x for x, _ in sys_below):.2f} V_B)")
        if sys_below:
            qs = [q for _, q in sys_below]
            below.append(f"{target} {min(qs):.3g}–{max(qs):.2f} at {min(x for x, _ in sys_below):.2f}–"
                         f"{max(x for x, _ in sys_below):.2f} V_B")
        if sys_above:
            qs = [q for _, q in sys_above]
            above.append(f"{target} {min(qs):.2f}–{max(qs):.2f} at {min(x for x, _ in sys_above):.2f}–"
                         f"{max(x for x, _ in sys_above):.2f} V_B")
    detail = ("worst in the window: " + ", ".join(judged) + "; physim/measured below (not judged): "
              + "; ".join(below) + "; above (not judged): " + ("; ".join(above) or "none"))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(xs)), "xlabel": "E_cm / V_B (physim's Bass barrier)",
            "ylabel": "σ_fus, physim / measured"}
    return worst <= FUSION_DATA_TOLERANCE, worst, n, detail, plot


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "tool", "PACE4 4.36.0, 20 000 cascades, five systems from ⁴He + "
       "¹⁶⁵Ho to ⁴⁸Ca + ¹²⁴Sn (pace4_*MeV.csv); physim runs 20 000 events from PACE4's E* and σ_J",
       "the strongest three residues the same as PACE4's, in the same order, each share within 30%")
def evaporation_channels_vs_pace4(ref):
    files = _pace4_files(ref, "evaporation")
    if not files:
        return None, math.nan, 0, "pace4_*MeV.csv not found", {}
    ok, worst, notes = True, 0.0, []
    xs, ratio = [], []
    for path, beam, target in files:
        info, pace = _pace4_evaporation(path)
        r = physim_evaporation(beam, target, info["excitation_mev"], info["sigma_j_mb"])
        mine: dict = {}
        for row in r.channel_table():
            if row["Z"] is not None:
                mine[(row["Z"], row["A"])] = mine.get((row["Z"], row["A"]), 0.0) + row["share"]
        top_p = sorted(pace, key=lambda k: -pace[k])[:3]
        top_m = sorted(mine, key=lambda k: -mine[k])[:3]
        devs = [mine.get(k, 0.0) / pace[k] - 1 for k in top_p]
        same = top_m == top_p
        sys_ok = same and all(abs(d) <= 0.30 for d in devs)
        ok &= sys_ok
        worst = max(worst, *(abs(d) for d in devs))
        xs += [f"{_nuclide_label(k)} ({beam})" for k in top_p]
        ratio += [1 + d for d in devs]
        notes.append(f"{beam} + {target}: " + ", ".join(
            f"{_nuclide_label(k)} {pace[k]:.3f}/{mine.get(k, 0.0):.3f}" for k in top_p)
            + ("" if same else " (physim's order: " + ", ".join(_nuclide_label(k) for k in top_m) + ")")
            + (" ✓" if sys_ok else " ✗"))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(xs)), "categorical": True}
    return bool(ok), worst, 3 * len(files), "shares PACE4/physim: " + "; ".join(notes), plot


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "tool", "PACE4 4.36.0, mean CM energies of the emitted neutrons "
       "and protons, the same five systems and runs (pace4_*MeV.csv)", "⟨ε_n⟩ and ⟨ε_p⟩ within 10%; α reported, "
       "not judged")
def evaporation_spectra_vs_pace4(ref):
    files = _pace4_files(ref, "evaporation")
    if not files:
        return None, math.nan, 0, "pace4_*MeV.csv not found", {}
    worst, notes, xs, ratio = 0.0, [], [], []
    for path, beam, target in files:
        info, _ = _pace4_evaporation(path)
        r = physim_evaporation(beam, target, info["excitation_mev"], info["sigma_j_mb"])
        d = {k: r.mean_energy_mev(k) / info["mean"][k] - 1 for k in ("n", "p", "a")}
        worst = max(worst, abs(d["n"]), abs(d["p"]))
        xs += [f"n ({beam})", f"p ({beam})"]
        ratio += [1 + d["n"], 1 + d["p"]]
        notes.append(f"{beam}: n {info['mean']['n']:.2f}/{r.mean_energy_mev('n'):.2f} ({d['n']:+.0%}), "
                     f"p {info['mean']['p']:.2f}/{r.mean_energy_mev('p'):.2f} ({d['p']:+.0%}), "
                     f"α {d['a']:+.0%}")
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(xs)), "categorical": True}
    return worst <= 0.10, worst, 2 * len(files), "mean CM energy PACE4/physim (MeV): " + "; ".join(notes), plot


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "literature", "Weisskopf's evaporation spectrum: with only "
       "neutrons, the first neutron from ¹⁶⁸Yb at 30 MeV and spin 0 follows ε σ_inv(ε) ρ(U − ε), 0.5–8 MeV "
       "(4 × 10⁴ events)", "χ²/dof < 2")
def evaporation_weisskopf(ref):
    from . import evaporation as ev

    z, a, e = 70, 168, 30.0
    sep = (data.nuclide((z, a - 1)).nuclear_mass_mev + data.nuclide((0, 1)).nuclear_mass_mev
           - data.nuclide((z, a)).nuclear_mass_mev)
    r = ev.evaporate((z, a), e, spin=0, events=40000, seed=3, decays=("n",))
    first = r.energy_mev[r.offsets[:-1]]
    edges = np.arange(0.5, 8.0 + 1e-9, 0.25)
    counts, _ = np.histogram(first, edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    t = ev.neutron_transmission(mid, a - 1, 20)
    g = ((2 * np.arange(21) + 1)[:, None] * t).sum(axis=0)  # ε σ_inv(ε) up to a constant
    rho = ev.level_density(z, a - 1, e - sep - mid, 0)
    shape = g * rho
    expected = shape / shape.sum() * counts.sum()
    chi2 = float(((counts - expected) ** 2 / expected).sum() / (mid.size - 1))
    plot = {"x": mid, "physim": counts, "reference": expected, "xlabel": "first neutron's energy (MeV)",
            "ylabel": "events per 0.25 MeV", "logy": True}
    return chi2 < 2.0, chi2 / 2.0, mid.size, f"χ²/dof = {chi2:.2f}", plot


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "literature", "Energy and momentum conservation in every event: "
       "¹⁶²Dy at 59.3 MeV with a fusion-like spin distribution, 4000 events", "1e-4 MeV of 1.5 × 10⁵ MeV; 1e-6 MeV/c")
def evaporation_conservation(ref):
    from . import evaporation as ev

    sigma_l = np.array([(2 * l + 1) / (1 + math.exp((l - 30) / 2.0)) for l in range(50)])
    run_ = ev.evaporate("162Dy", 59.3, sigma_l_mb=sigma_l, events=4000, seed=7)
    z, a = run_.compound
    total = data.nuclide((z, a)).nuclear_mass_mev + run_.excitation_mev
    m_res = np.array([data.nuclide((int(zz), int(aa))).nuclear_mass_mev
                      for zz, aa in zip(run_.residue_z, run_.residue_a)])
    e_res = np.sqrt((run_.recoil_momentum_mev_c**2).sum(axis=1) + (m_res + run_.entry_energy_mev) ** 2)
    e_emit = np.bincount(run_.event, weights=run_.energy_mev + run_.mass_mev, minlength=run_.events)
    d_e = float(np.abs(e_emit + e_res - total).max())
    p_sum = np.stack([np.bincount(run_.event, weights=run_.momentum_mev_c[:, i], minlength=run_.events)
                      for i in range(3)], axis=1)
    d_p = float(np.abs(p_sum + run_.recoil_momentum_mev_c).max())
    worst = max(d_e / 1e-4, d_p / 1e-6)
    return worst <= 1.0, worst, run_.events, f"worst ΔE {d_e:.1e} MeV, Δp {d_p:.1e} MeV/c", {}


def _gemini_files(ref) -> list:
    """GEMINI++ reference files ``gemini_<beam>_<target>_<E>MeV.csv``, with (beam, target) from the name."""
    if ref is None:
        return []
    out = [(p, p.stem.split("_")[1], p.stem.split("_")[2]) for p in sorted(Path(ref).glob("gemini_*MeV.csv"))]
    order = {b: i for i, b in enumerate(("4He", "12C", "16O", "32S", "48Ca"))}
    return sorted(out, key=lambda t: (order.get(t[1], 99), t[1]))


def _kv(text: str) -> dict:
    return {k.strip(): float(v) for k, v in (kv.split("=") for kv in text.split(","))}


def _gemini_evaporation(path) -> tuple:
    """(header numbers, GEMINI++ residue events by (Z, A)) of one GEMINI++ reference file."""
    header, rows = read_reference(path)
    counts = {(int(r["Z"]), int(r["A"])): int(r["events"]) for r in rows}
    spin = _kv(header["spin"])
    info = {"excitation_mev": float(header["excitation_mev"]), "sigma_fus_mb": float(header["sigma_fus_mb"]),
            "l0": spin["l0_hbar"], "diffuseness": spin["diffuseness_hbar"],
            "mean": _kv(header["mean_cm_energy_mev"]), "multiplicity": _kv(header["multiplicity"]),
            "events": sum(counts.values())}
    return info, counts


def gemini_spins(info: dict) -> tuple:
    """GEMINI++'s spin distribution, σ_ℓ ∝ (2ℓ+1) / (1 + exp((ℓ − L0)/Δ)), normalised to its σ_fus (mb)."""
    l0, dif = info["l0"], info["diffuseness"]
    ell = np.arange(int(l0 + 12 * dif) + 1)
    w = (2 * ell + 1) / (1 + np.exp((ell - l0) / dif))
    return tuple(w / w.sum() * info["sigma_fus_mb"])


def _residue_shares(r) -> dict:
    out: dict = {}
    for row in r.channel_table():
        if row["Z"] is not None:
            out[(row["Z"], row["A"])] = out.get((row["Z"], row["A"]), 0.0) + row["share"]
    return out


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "tool", "GEMINI++ 3.3.0 (in LISE++), Hauser–Feshbach mode, "
       "⁴He + ¹⁶⁵Ho, ¹²C + ¹⁵⁰Nd (20 000 events) and ¹⁶O + ¹⁴⁴Sm (1500; GEMINI++ hangs on more, and on ³²S and ⁴⁸Ca) "
       "(gemini_*MeV.csv); physim runs 20 000 events from GEMINI++'s E* and spin distribution",
       "the strongest three residues the same as GEMINI++'s, in the same order (residues whose GEMINI++ shares agree "
       "within 2 standard errors count as tied), each share within 30%")
def evaporation_channels_vs_gemini(ref):
    files = _gemini_files(ref)
    if not files:
        return None, math.nan, 0, "gemini_*MeV.csv not found", {}
    ok, worst, notes = True, 0.0, []
    xs, ratio = [], []
    for path, beam, target in files:
        info, counts = _gemini_evaporation(path)
        n = info["events"]
        gem = {k: c / n for k, c in counts.items()}
        err = {k: math.sqrt(c) / n for k, c in counts.items()}
        mine = _residue_shares(physim_evaporation(beam, target, info["excitation_mev"], gemini_spins(info)))
        top_g = sorted(gem, key=lambda k: -gem[k])[:3]
        top_m = sorted(mine, key=lambda k: -mine[k])[:3]
        same = all(abs(gem.get(km, 0.0) - gem[kg]) <= 2 * math.hypot(err.get(km, 1 / n), err[kg])
                   for km, kg in zip(top_m, top_g))
        devs = [mine.get(k, 0.0) / gem[k] - 1 for k in top_g]
        sys_ok = same and all(abs(d) <= 0.30 for d in devs)
        ok &= sys_ok
        worst = max(worst, *(abs(d) for d in devs))
        xs += [f"{_nuclide_label(k)} ({beam})" for k in top_g]
        ratio += [1 + d for d in devs]
        notes.append(f"{beam} + {target}: " + ", ".join(
            f"{_nuclide_label(k)} {gem[k]:.3f}/{mine.get(k, 0.0):.3f}" for k in top_g)
            + ("" if top_m == top_g else " (physim's order: " + ", ".join(_nuclide_label(k) for k in top_m)
               + (", tied within statistics" if same else "") + ")")
            + (" ✓" if sys_ok else " ✗"))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(xs)), "categorical": True}
    return bool(ok), worst, 3 * len(files), "shares GEMINI++/physim: " + "; ".join(notes), plot


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "tool", "GEMINI++ 3.3.0, mean CM energies of the emitted "
       "neutrons and protons, the same three systems and runs (gemini_*MeV.csv)", "⟨ε_n⟩ and ⟨ε_p⟩ within 10%; α "
       "reported, not judged")
def evaporation_spectra_vs_gemini(ref):
    files = _gemini_files(ref)
    if not files:
        return None, math.nan, 0, "gemini_*MeV.csv not found", {}
    worst, notes, xs, ratio = 0.0, [], [], []
    for path, beam, target in files:
        info, _ = _gemini_evaporation(path)
        r = physim_evaporation(beam, target, info["excitation_mev"], gemini_spins(info))
        d = {k: r.mean_energy_mev(k) / info["mean"][k] - 1 for k in ("n", "p", "a")}
        worst = max(worst, abs(d["n"]), abs(d["p"]))
        xs += [f"n ({beam})", f"p ({beam})"]
        ratio += [1 + d["n"], 1 + d["p"]]
        notes.append(f"{beam}: n {info['mean']['n']:.2f}/{r.mean_energy_mev('n'):.2f} ({d['n']:+.0%}), "
                     f"p {info['mean']['p']:.2f}/{r.mean_energy_mev('p'):.2f} ({d['p']:+.0%}), α {d['a']:+.0%}")
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(xs)), "categorical": True}
    return worst <= 0.10, worst, 2 * len(files), "mean CM energy GEMINI++/physim (MeV): " + "; ".join(notes), plot


#: How far outside the band spanned by PACE4 and GEMINI++ physim's ⟨ε_n⟩ may lie (relative to the nearer edge).
NEUTRON_BAND_MARGIN = 0.10


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "tool", "PACE4 4.36.0 and GEMINI++ 3.3.0 together: mean CM "
       "energy of the evaporated neutrons for the three systems both codes ran (pace4_*MeV.csv, gemini_*MeV.csv); "
       "physim from PACE4's E* and σ_J as in the PACE4 check", "⟨ε_n⟩ inside the band the two codes span, widened by "
       "10% at each edge: the two codes disagree by 25–35%, so their spread is the size of the model uncertainty, "
       "and 10% is backlog 75's tolerance against one code")
def evaporation_neutron_energy_band(ref):
    gem = {(b, t): p for p, b, t in _gemini_files(ref)}
    pace = {(b, t): p for p, b, t in _pace4_files(ref, "evaporation")}
    both = [k for k in pace if k in gem]
    if not both:
        return None, math.nan, 0, "no system with both pace4_*MeV.csv and gemini_*MeV.csv", {}
    worst, ok, notes, xs, ours, mids = 0.0, True, [], [], [], []
    for beam, target in both:
        p_info, _ = _pace4_evaporation(pace[(beam, target)])
        g_info, _ = _gemini_evaporation(gem[(beam, target)])
        e = physim_evaporation(beam, target, p_info["excitation_mev"], p_info["sigma_j_mb"]).mean_energy_mev("n")
        lo, hi = sorted((p_info["mean"]["n"], g_info["mean"]["n"]))
        d = e / hi - 1 if e > hi else (e / lo - 1 if e < lo else 0.0)
        ok &= abs(d) <= NEUTRON_BAND_MARGIN
        worst = max(worst, abs(d))
        xs.append(f"{beam} + {target}")
        ours.append(e)
        mids.append(0.5 * (lo + hi))
        notes.append(f"{beam}: physim {e:.2f} MeV, band {lo:.2f}–{hi:.2f}" + (f" ({d:+.0%} outside)" if d else ""))
    plot = {"x": xs, "physim": ours, "reference": mids, "categorical": True}
    return bool(ok), worst, len(both), "; ".join(notes) + f" (worst {worst:.0%} outside the band)", plot


#: Evaporation events per beam energy in the comparison with measured residue cross sections.
RESIDUE_EVENTS = 4000
#: Allowed deviation of physim's xn cross section from the band of the three measurements (each ± its error):
#: backlog 75's 30% for the strongest channels.
RESIDUE_TOLERANCE = 0.30
#: The channels judged where they are the strongest measured channel (1n is the strongest only below the barrier).
RESIDUE_JUDGED = ("2n", "3n", "4n")


def _measured_curves(rows) -> dict:
    """{(ref, channel): (energies, σ, relative error)} sorted by energy."""
    out: dict = {}
    for r in rows:
        out.setdefault((r["ref"], r["channel"]), []).append(
            (float(r["e_lab_mev"]), float(r["sigma_mb"]), float(r["uncertainty_mb"]) / float(r["sigma_mb"])))
    return {k: tuple(np.array(c) for c in zip(*sorted(v))) for k, v in out.items()}


def _measured_band(curves: dict, channel: str, e: float):
    """(low, high) of the measurements of one channel at energy e: each set interpolated (log σ, linear in E) where
    it covers e, ± its error; None where no set does."""
    lo, hi = math.inf, 0.0
    for (_, ch), (en, sig, rel) in curves.items():
        if ch != channel or not en[0] <= e <= en[-1]:
            continue
        s = float(np.exp(np.interp(e, en, np.log(sig))))
        u = float(np.interp(e, en, rel))
        lo, hi = min(lo, s * (1 - u)), max(hi, s * (1 + u))
    return (lo, hi) if hi > 0 else None


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "literature", "Measured residue cross sections ¹⁶⁵Ho(α,xn), "
       "1n–4n and 6n, 13–66 MeV: Usman et al. 2020, Kormazeva et al. 2022, Gadkari et al. 1997, from EXFOR "
       "(evaporation_4He_165Ho_literature.csv); physim as the planner computes it: σ_fus × the residue's share of "
       f"{RESIDUE_EVENTS} evaporation events at each measured energy (thin target)",
       f"2n, 3n, 4n where each is the strongest measured channel (4n only up to its measured peak: (α,5n) was not "
       f"measured) and E_cm ≥ {FUSION_ABOVE_BARRIER} V_B: within "
       f"{RESIDUE_TOLERANCE:.0%} of the band the three measurements span (each ± its error); the tails above each "
       "peak (pre-equilibrium emission, which physim leaves out), 1n and 6n reported, not judged")
def residues_vs_data(ref):
    from .fusion import Fusion

    path = _file(ref, "evaporation_4He_165Ho_literature.csv")
    if path is None:
        return None, math.nan, 0, "evaporation_4He_165Ho_literature.csv not found", {}
    rows = read_reference(path)[1]
    curves = _measured_curves(rows)
    f = Fusion("4He", "165Ho")
    vb = f.barrier().height_mev
    zc, ac = f.compound_nuclide.Z, f.compound_nuclide.A
    physim: dict = {}

    def sigma(e: float, channel: str) -> float:
        if e not in physim:
            cn = f.compound(e)
            r = physim_evaporation("4He", "165Ho", cn.excitation_mev, tuple(cn.sigma_l_mb), events=RESIDUE_EVENTS)
            physim[e] = (cn.sigma_mb, _residue_shares(r))
        s, shares = physim[e]
        return s * shares.get((zc, ac - int(channel[:-1])), 0.0)

    at: dict = {}  # (ref, energy) -> {channel: σ}
    for r in rows:
        at.setdefault((r["ref"], float(r["e_lab_mev"])), {})[r["channel"]] = float(r["sigma_mb"])

    def strongest(ref_name: str, e: float, ch: str) -> bool:
        """Whether ch is the strongest channel at e in that data set. Above a channel's measured peak the next
        channel takes over; where that one was not measured ((α,5n)), the strongest channel is unknown."""
        here = at[(ref_name, e)]
        if max(here, key=here.get) != ch:
            return False
        if (ref_name, f"{int(ch[:-1]) + 1}n") in curves:
            return True
        en, sig, _ = curves[(ref_name, ch)]
        return e <= en[int(np.argmax(sig))]

    worst, judged, xs, ratio, labels = 0.0, set(), [], [], {}
    others: dict = {}
    for r in rows:
        e, ch = float(r["e_lab_mev"]), r["channel"]
        if float(f.e_cm(e)) / vb < FUSION_ABOVE_BARRIER:
            continue
        ours = sigma(e, ch)
        if ch in RESIDUE_JUDGED and strongest(r["ref"], e, ch):
            if (ch, e) in judged:
                continue
            judged.add((ch, e))
            lo, hi = _measured_band(curves, ch, e)
            d = ours / hi - 1 if ours > hi else (ours / lo - 1 if ours < lo else 0.0)
            if abs(d) > worst:
                worst = abs(d)
                labels["worst"] = f"{ch} at {e:g} MeV: physim {ours:.0f} mb, measured {lo:.0f}–{hi:.0f} mb"
            xs.append(e)
            ratio.append(ours / float(r["sigma_mb"]))
        else:
            others.setdefault(ch, []).append((e, ours / float(r["sigma_mb"])))
    per = {}
    for ch in RESIDUE_JUDGED:
        es = sorted(e for c, e in judged if c == ch)
        per[ch] = f"{ch} {es[0]:g}–{es[-1]:g} MeV" if es else f"{ch} never strongest"
    rest = "; ".join(f"{ch} {min(q for _, q in v):.2g}–{max(q for _, q in v):.2g} ({len(v)} points)"
                     for ch, v in sorted(others.items()))
    near50 = [f"{q:.2f} ({e:g} MeV)" for e, q in sorted(others.get("3n", [])) if 49 <= e <= 51]
    detail = (f"{len(judged)} judged points ({', '.join(per.values())}); worst {worst:.0%} outside the band"
              + (f" ({labels['worst']})" if labels else "") + "; physim/measured elsewhere (not judged): " + rest
              + ("; 3n at 49–51 MeV (pre-equilibrium tail): " + ", ".join(near50) if near50 else ""))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(xs)), "xlabel": "α energy (MeV, lab)",
            "ylabel": "σ(α,xn), physim / measured, where judged"}
    return worst <= RESIDUE_TOLERANCE, worst, len(judged), detail, plot


@check("Evaporation Monte Carlo (Hauser–Feshbach)", "literature", "Reproduced experiment: Usman et al., Nucl. "
       "Instrum. Methods B 469, 42 (2020), ⁴He + ¹⁶⁵Ho stacked-foil activation at RIKEN (EXFOR E2650), end to end "
       "through the planner: each of the 26 Ho foils as the fusion_tm165 setup at the foil's energy, the channel "
       "table filled from physim's fusion and evaporation (Fill from physim, 20 000 chains), the residues simulated "
       "and σ inferred from the residues made (physim.nuclear.reproduction; saved as reproduction_usman2020.csv); "
       "Kormazeva 2022 and Gadkari 1997 as cross-checks",
       "Above 1.1 V_B: the peaks of 2n, 3n, 4n within the energy spread + 1 MeV; the strongest channel the same at "
       "every foil where the measurement resolves it; the strongest three channels above 10% of the measured total "
       "within 30% of the band of the three measurements (± error, ± energy uncertainty), up to each channel's "
       "peak. Reported, not judged: the tails above the peaks (pre-equilibrium emission, which physim leaves out), "
       "1n, σ(3n)/σ(4n) at 45–50 MeV")
def reproduction_usman2020(ref):
    from . import reproduction as rp

    path = _file(ref, rp.REFERENCE_FILE)
    if path is None:
        return None, math.nan, 0, f"{rp.REFERENCE_FILE} not found", {}
    rows = read_reference(path)[1]
    rep = rp.run(energies=rp.foil_energies(rows))
    c = rp.compare(rep.rows, rows)
    judged = [x for k in ("peaks", "strongest", "magnitudes") for x in c[k] if x["judged"]]
    worst = max(c["worst"]["peaks"], c["worst"]["magnitudes"])
    wrong = [f"{x['e_lab_mev']:g} MeV (physim {x['physim']}, measured {x['measured']})" for x in c["strongest"]
             if x["judged"] and not x["passed"]]
    detail = (rp.summary(c) + (f"; strongest channel wrong at {', '.join(wrong)}" if wrong else "") + ". "
              + rp.reported(c) + f"; worst {worst:.2f} of the peak and magnitude tolerances")
    mags = [x for x in c["magnitudes"] if x["judged"]]
    plot = {"x": [x["e_lab_mev"] for x in mags], "physim": [x["ratio"] for x in mags],
            "reference": np.ones(len(mags)), "xlabel": "α energy in the foil (MeV, lab)",
            "ylabel": "σ(α,xn), physim (planner) / Usman 2020, where judged"}
    return c["passed"], worst, len(judged), detail, plot


@check("Neutron detectors and channel selection", "literature", "Knoll, Radiation Detection and Measurement (4th ed., "
       "ch. 15), recoil-proton scintillator: ε = n_Hσ_H/Σ · (1 − e^(−ΣD)) · (1 − E_B/E_n), EJ-301 cell, 1–8 MeV",
       "1e-12")
def neutron_efficiency_knoll(ref):
    from . import neutrons as N

    cell = N.NeutronDetector()
    s = N.material(cell.material)
    d_cm = cell.depth_mm / 10
    e_b = cell.proton_threshold()
    xs, ours, theirs = [], [], []
    for energy in (1.0, 1.5, 2.5, 4.0, 8.0):
        sh = s["hydrogen"] * float(N.hydrogen_cross_section(energy)) * 1e-24
        sc = s["carbon"] * float(N.carbon_cross_section(energy)) * 1e-24
        xs.append(energy)
        theirs.append(sh / (sh + sc) * (1 - math.exp(-(sh + sc) * d_cm)) * (1 - e_b / energy))
        ours.append(float(cell.intrinsic_efficiency(energy)))
    worst = float(np.max(_rel(ours, theirs)))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "neutron energy (MeV)",
            "ylabel": "intrinsic efficiency"}
    return worst <= 1e-12, worst, len(xs), f"worst {worst:.1e}", plot


@check("Neutron detectors and channel selection", "literature", "The binomial expectation: closed forms for 3n at "
       "ε = 0.3, and the event generator's selection matrix (7 channels, cross-talk, γ rays, particle "
       "misidentification; 6 × 10⁴ events) against the analytic one", "closed forms 1e-12; every Monte Carlo cell "
       "within 4.5 binomial standard errors")
def channel_selection_binomial(ref):
    from . import neutrons as N

    eps = 0.3
    m3 = N.selection_matrix(N.channels_from(["n", "2n", "3n"]), N.Detection(neutron_efficiency=eps))
    closed = [(m3.share("3n", "3n"), eps**3), (m3.share("3n", "2n"), 3 * eps**2 * (1 - eps)),
              (m3.share("3n", "n"), 3 * eps * (1 - eps) ** 2), (m3.share("3n", N.UNASSIGNED), (1 - eps) ** 3)]
    d_closed = max(abs(a / b - 1) for a, b in closed)
    channels = N.channels_from(["n", "2n", "3n", "pn", "p2n", "2pn", "αn"], gammas=8)
    full = N.Detection(neutron_efficiency=0.25, proton_efficiency=0.6, alpha_efficiency=0.4, crosstalk=0.06,
                       gamma_efficiency=0.1, gamma_as_neutron=0.02, neutron_as_gamma=0.03, proton_as_alpha=0.02,
                       alpha_as_proton=0.05)
    events = N.isotropic_events(channels, events=60_000, seed=21)
    mc = N.simulate_selection(events, channels, full, seed=22)
    m = N.selection_matrix(channels, full)
    n = mc.events[:, None]
    err = np.sqrt(m.matrix * (1 - m.matrix) / n) + 1.0 / n
    z = float((np.abs(mc.matrix - m.matrix) / err).max())
    worst = max(z / 4.5, d_closed / 1e-12)
    plot = {"x": m.matrix.ravel(), "physim": mc.matrix.ravel(), "reference": m.matrix.ravel(),
            "xlabel": "share, binomial expectation", "ylabel": "share, Monte Carlo"}
    return worst <= 1.0, worst, int(m.matrix.size) + 4, (f"worst Monte Carlo cell {z:.1f} standard errors; closed "
                                                         f"forms {d_closed:.1e}"), plot


@check("Neutron detectors and channel selection", "literature", "ENDF/B-VIII.0 total neutron cross section of "
       "natural carbon (¹²C and ¹³C, NNDC), 0.5–20 MeV, averaged over ±2.5% and ±10% of the energy "
       "(neutron_carbon_endf.csv). An evaluated library is judged data, like NIST's PSTAR: measurements combined "
       "and fitted by evaluators, published in Nucl. Data Sheets 148, 1 (2018), not a code run on physim's inputs",
       "physim's table (bins of ±2.9%) within 0.2% of both averages where the cross section is smooth (below "
       "1.8 MeV, ±10% below 1.7 MeV); above, within 7% of the ±2.5% average (resonances narrower than a bin) and "
       "5% of the ±10% average")
def carbon_vs_endf(ref):
    from . import neutrons as N

    path = _file(ref, "neutron_carbon_endf.csv")
    if path is None:
        return None, math.nan, 0, "neutron_carbon_endf.csv not found", {}
    rows = [r for r in read_reference(path)[1] if r["quantity"] == "natC_total"]
    worst, xs, ours, theirs = 0.0, [], [], []
    for r in rows:
        e, mean, wide = float(r["energy_mev"]), float(r["mean_b"]), float(r["wide_b"])
        point = float(N.carbon_cross_section(e))
        x = np.linspace(0.9 * e, 1.1 * e, 4001)
        y = np.asarray(N.carbon_cross_section(x), dtype=float)
        averaged = float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x))) / (0.2 * e)
        worst = max(worst, abs(point / mean - 1) / (0.002 if e < 1.8 else 0.07),
                    abs(averaged / wide - 1) / (0.002 if e < 1.7 else 0.05))
        xs.append(e)
        ours.append(point)
        theirs.append(mean)
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "neutron energy (MeV)",
            "ylabel": "σ_tot, natural carbon (b)", "logx": True}
    return worst <= 1.0, worst, 2 * len(rows), f"worst {worst:.0%} of its tolerance", plot


# ---------------------------------------------------------------------------------------------------------------
# Shell model (KSHELL)

#: Published calculations with the same Hamiltonians as the KSHELL fixtures (see the file's header).
SHELLMODEL_BENCHMARKS = "shellmodel/benchmarks.csv"

#: ENSDF's ²⁸Si levels and B(E2) below 10 MeV as physim reads them, kept so the context runs without ENSDF.
SHELLMODEL_ENSDF_SI28 = "shellmodel/si28_ensdf.csv"

#: Published dimensions of full-space calculations: (nucleus, space, Z_val, N_val, 2J or None for the M-scheme
#: dimension at M = 0, value, source). Positive parity. Exact counts.
PUBLISHED_SHELL_DIMENSIONS = [
    ("²⁸Si", "sd", 6, 6, None, 93_710, "Yoshida and Shimizu, PTEP 2022 (arXiv:2105.08256v2), Sec. 2.2"),
    ("⁴⁸Cr", "fp", 4, 4, None, 1_963_461, "Pittel, Thakur and Sandulescu, IJMPE 17, 122 (2009) (arXiv:0808.1303), "
     "Sec. 3.1"),
    ("⁴⁸Cr", "fp", 4, 4, 0, 41_355, "same, Sec. 3.1 and Tables 1, 2"),
    ("⁴⁸Cr", "fp", 4, 4, 4, 182_421, "same, Sec. 3.1 and Table 2"),
    ("⁴⁸Cr", "fp", 4, 4, 8, 246_979, "same, Table 2"),
    ("⁴⁸Cr", "fp", 4, 4, 12, 226_259, "same, Table 2"),
    ("⁵⁶Ni", "fp", 8, 8, None, 1_087_455_228, "same, Sec. 3.2"),
    ("⁵⁶Ni", "fp", 8, 8, 0, 15_443_684, "same, Sec. 3.2"),
]


def _printed_step(text: str) -> float:
    """One unit in the last digit of a number as printed: ``"78.3"`` → 0.1, ``"110"`` → 1."""
    text = text.strip()
    return 10.0 ** -len(text.split(".")[1]) if "." in text else 1.0


def _shellmodel_value(result, quantity: str, state: str) -> Optional[float]:
    """physim's value of one benchmark row from a :class:`~physim.nuclear.shellmodel.ShellModelResult`."""
    from . import shellmodel as sm

    if quantity == "dimension":  # "M=0+": 2M = 0, positive parity
        m2, parity = state.removeprefix("M=")[:-1], 1 if state.endswith("+") else -1
        return result.dimensions.get((2 * int(m2), parity))
    if quantity == "B(E2)":
        a, b = (sm.parse_label(s) for s in state.split(">"))
        t = result.transition("E2", a, b)
        if t is None:
            return None
        return t.b_down if t.initial == a else t.b_up
    s = result.state(*sm.parse_label(state))
    if s is None:
        return None
    return {"energy": s.energy_mev, "Q_s": s.q_s}.get(quantity)


@check("Shell-model calculations (KSHELL)", "tool", "Published calculations with the same Hamiltonian, against "
       "physim's reading of the bundled KSHELL's logs (shellmodel/benchmarks.csv; logs in shellmodel/kshell/): "
       "²⁸Si/USDB E(0⁺₁), B(E2; 2⁺₁→0⁺₁), B(E2; 4⁺₁→2⁺₁) and Q_s(2⁺₁) with e_p 1.31, e_n 0.46 from ANTOINE (Frycz et "
       "al., Phys. Rev. C 110, 054326 (2024), arXiv:2404.14506, Tables 2, 3, 5); ⁴⁸Cr/KB3 exact energies of 0⁺₁, "
       "2⁺₁, 4⁺₁, 6⁺₁, 0⁺₂ (Pittel et al., Int. J. Mod. Phys. E 17, 122 (2009), arXiv:0808.1303, Tables 1, 2); the "
       "M-scheme dimensions KSHELL printed", "one unit in the last digit the source prints (its rounding, and "
       "constants such as the oscillator length that it does not state); dimensions exact")
def shellmodel_vs_published(ref):
    from . import shellmodel as sm

    path = _file(ref, SHELLMODEL_BENCHMARKS)
    if path is None:
        return None, math.nan, 0, f"{SHELLMODEL_BENCHMARKS} not found", {}
    rows = read_reference(path)[1]
    results: dict = {}
    worst, worst_text, labels, ours, theirs, outside = 0.0, "", [], [], [], []
    for r in rows:
        folder = path.parent / r["folder"]
        if folder not in results:
            results[folder] = sm.read_folder(folder)
        value, published = _shellmodel_value(results[folder], r["quantity"], r["state"]), float(r["value"])
        name = f"{r['case']} {r['quantity']} {r['state']}"
        if value is None:
            outside.append(f"{name}: not in the logs")
            worst = math.inf
            continue
        allowed = 0.0 if r["quantity"] == "dimension" else _printed_step(r["value"])
        dev = abs(value - published)
        frac = (0.0 if dev == 0 else math.inf) if allowed == 0 else dev / allowed
        if frac > worst:
            worst, worst_text = frac, f"{name}, {value:.6g} against {r['value']}"
        if frac > 1:
            outside.append(f"{name}: physim {value:.6g} against {r['value']} ({r['source']})")
        labels.append(name)
        ours.append(value)
        theirs.append(published)
    plot = {"x": labels, "physim": ours, "reference": theirs, "categorical": True,
            "xlabel": "benchmark", "ylabel": "physim / published"}
    detail = (f"{len(rows)} numbers from {len(results)} sets of logs, "
              + ("all within one unit of the printed digit" if not outside else "outside: " + "; ".join(outside))
              + f"; worst {worst:.0%} of the tolerance ({worst_text})")
    return not outside, worst, len(rows), detail, plot


def _si28_ensdf_scheme(path):
    """A :class:`~physim.nuclear.levels.LevelScheme` from :data:`SHELLMODEL_ENSDF_SI28`."""
    from .ensdf import spin_parity
    from .levels import Level, LevelScheme, Value

    header, rows = read_reference(path)
    scheme = LevelScheme("28Si", reference=header.get("tool", "ENSDF"))
    for r in rows:
        if r["kind"] == "level":
            j, parity, _ = spin_parity(r["jpi"])
            unc = float(r["unc"]) if r["unc"] else None
            scheme.levels.append(Level(Value(float(r["value"]), unc, "ensdf"), j, parity, r["jpi"], None))
    for r in rows:
        if r["kind"] == "BE2":
            i, f = int(r["initial"]), int(r["final"])
            me = math.sqrt(float(r["value"]) * (2 * scheme.levels[i].spin + 1))
            scheme.set_matrix_element(i, f, "E2", me, source="ensdf")
    return scheme


def shellmodel_si28_vs_ensdf(ref) -> Optional[dict]:
    """²⁸Si/USDB from the fixtures against ENSDF, as context (shell-model papers report it; it is not a test of
    physim): every shell-model state matched to a known level by J^π and order, the deviations of the excitation
    energies and their rms, and B(E2) (e_p 1.31, e_n 0.46) over ENSDF's where ENSDF has one. ``None`` when the
    fixtures or the ENSDF snapshot are missing."""
    from . import shellmodel as sm

    base = _file(ref, "shellmodel/kshell/si28_usdb")
    snap = _file(ref, SHELLMODEL_ENSDF_SI28)
    if base is None or snap is None:
        return None
    energies = sm.read_folder(base / "logs")
    charged = sm.read_folder(base / "logs_eff_1.31_0.46")
    scheme = _si28_ensdf_scheme(snap)
    matches = [m for m in sm.match(energies, scheme) if m.level is not None]
    excited = [m for m in matches if m.level != 0]
    level_of = {m.state.key: m.level for m in matches}
    be2 = []
    for t in charged.transitions:
        if t.kind != "E2" or t.diagonal or t.initial not in level_of or t.final not in level_of:
            continue
        known = scheme.b(level_of[t.initial], level_of[t.final], "E2")
        if known:
            be2.append((sm.label(*t.initial), sm.label(*t.final), t.b_down, known))
    rms = math.sqrt(sum(m.delta_kev**2 for m in excited) / len(excited)) if excited else math.nan
    return {"matches": matches, "rms_kev": rms, "excited": len(excited), "be2": be2}


def _shellmodel_context_text(ctx: Optional[dict]) -> str:
    if ctx is None:
        return "Context against ENSDF not available (no shell-model fixtures or ENSDF snapshot)"
    devs = ", ".join(f"{m.state.label} {m.delta_kev:+.0f}" for m in ctx["matches"] if m.level != 0)
    ratios = ", ".join(f"{a}→{b} {ours / known:.2f}" for a, b, ours, known in ctx["be2"])
    return (f"Context, not judged: ²⁸Si/USDB against ENSDF (shellmodel/si28_ensdf.csv), excitation energies rms "
            f"{ctx['rms_kev']:.0f} keV over {ctx['excited']} excited states matched by J^π and order (keV: {devs}); "
            f"B(E2), e_p 1.31, e_n 0.46, over ENSDF's: {ratios}")


@check("Shell-model calculations (KSHELL)", "literature", "Published dimensions of full-space calculations, exact "
       "counts: ²⁸Si sd shell at M = 0 (Yoshida and Shimizu, PTEP 2022, arXiv:2105.08256, Sec. 2.2); ⁴⁸Cr fp shell "
       "at M = 0 and for J = 0, 2, 4, 6, ⁵⁶Ni at M = 0 and J = 0 (Pittel et al., arXiv:0808.1303, Sec. 3 and "
       "Table 2), against physim's own count (shellmodel_spaces), which sizes every run. Reported as context, not "
       "judged: ²⁸Si/USDB excitation energies and B(E2) against ENSDF (rms deviation, as shell-model papers give it)",
       "exact; the ENSDF comparison is reported, not judged")
def shellmodel_dimensions(ref):
    from . import shellmodel_spaces as sp

    spaces = {"sd": sp.SD_SPACE, "fp": sp.FP_SPACE}
    labels, ours, theirs, wrong = [], [], [], []
    for nucleus, space, z, n, j2, published, _ in PUBLISHED_SHELL_DIMENSIONS:
        s = spaces[space]
        mine = sp.dimensions(s, z, n, 1).get(0, 0) if j2 is None else sp.j_dimension(s, z, n, 1, j2)
        name = f"{nucleus} " + ("M = 0" if j2 is None else f"J = {j2 // 2}")
        if mine != published:
            wrong.append(f"{name}: physim {mine}, published {published}")
        labels.append(name)
        ours.append(mine)
        theirs.append(published)
    plot = {"x": labels, "physim": ours, "reference": theirs, "categorical": True, "xlabel": "dimension",
            "ylabel": "physim / published"}
    detail = (f"{len(labels)} dimensions, " + ("all equal" if not wrong else "different: " + "; ".join(wrong))
              + ". " + _shellmodel_context_text(shellmodel_si28_vs_ensdf(ref)))
    return not wrong, 0.0 if not wrong else math.inf, len(labels), detail, plot


# ---------------------------------------------------------------------------------------------------------------
# Lifetime measurements: the plunger (RDDS) and the Doppler-shift attenuation method (DSAM)
#
# The references are in tests/reference/nuclear/lifetimes/ (its README gives their provenance): closed forms
# (analytic/), TRIM for the recoil transport (trim/), APCAD's and AhKin's DSAM examples (apcad/, ahkin/), NaPyTau's
# DDCM (napytau/) and the published experiments (published/). The recoils' stopping against measured data is the
# check stopping_heavy_vs_data (Stopping power and range), not repeated here. The DDCM against NaPyTau, τ recovered
# from simulated runs (closure) and the three published experiments simulate and analyse with
# physim.nuclear.lifetime_reproduction (wave 3C), at run sizes of seconds each.

RDDS = "Lifetimes: recoil distance (RDDS)"
DSAM = "Lifetimes: Doppler-shift attenuation (DSAM)"
#: The lifetime references' folder, inside the reference folder.
LIFETIME_REFERENCE = "lifetimes"


def _lifetime_file(ref, name: str) -> Optional[Path]:
    return _file(ref, f"{LIFETIME_REFERENCE}/{name}")


def _missing(name: str) -> tuple:
    return None, math.nan, 0, f"{LIFETIME_REFERENCE}/{name} not found", {}


def _share_deviation(got, want) -> np.ndarray:
    """|got/want − 1|, and for a share that underflowed to 0 in the reference, 0 if physim's is below 1e-300 too."""
    got, want = np.asarray(got, dtype=float), np.asarray(want, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(got / np.where(want > 0, want, 1.0) - 1)
    return np.where(want > 0, rel, np.where(np.abs(got) < 1e-300, 0.0, np.inf))


@check(RDDS, "literature", "Closed form for one level without feeding, the recoil along the beam at one speed: "
       "unshifted share exp(−d/(γβcτ)), and exp(−d/(βcτ)) without time dilation; τ = 1, 10, 100 ps, β = 2, 5, 10%, "
       "d = 1 µm to 5 mm (lifetimes/analytic/rdds_one_level.csv)", "1e-9 relative")
def rdds_one_level_analytic(ref):
    from . import lifetimes

    path = _lifetime_file(ref, "analytic/rdds_one_level.csv")
    if path is None:
        return _missing("analytic/rdds_one_level.csv")
    rows = read_reference(path)[1]
    worst, xs, ours, theirs = 0.0, [], [], []
    for r in rows:
        tau, beta, d = float(r["tau_ps"]), float(r["beta"]), float(r["distance_um"])
        dilated = lifetimes.unshifted([tau], None, [1.0], lifetimes.flight_times(d, beta, proper=True))[0]
        plain = lifetimes.unshifted([tau], None, [1.0], lifetimes.flight_times(d, beta))[0]
        worst = max(worst, float(_share_deviation(dilated, float(r["unshifted"]))),
                    float(_share_deviation(plain, float(r["unshifted_no_dilation"]))))
        if tau == 10.0 and beta == 0.05 and float(r["unshifted"]) > 1e-12:
            xs.append(d)
            ours.append(float(dilated))
            theirs.append(float(r["unshifted"]))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "target-to-stopper distance (µm), τ = 10 ps, "
            "β = 5%", "ylabel": "unshifted share", "logx": True, "logy": True}
    return worst <= 1e-9, worst, 2 * len(rows), f"{len(rows)} distances with and without time dilation, worst " \
        f"{worst:.1e}", plot


@check(RDDS, "literature", "Bateman's equations for a level fed from one level above: the unshifted shares of both "
       "transitions for three lifetime pairs (one with τ₁ ≈ τ₂) and two population splits "
       "(lifetimes/analytic/rdds_two_levels.csv)", "1e-9 relative")
def rdds_two_levels_bateman(ref):
    from . import lifetimes

    path = _lifetime_file(ref, "analytic/rdds_two_levels.csv")
    if path is None:
        return _missing("analytic/rdds_two_levels.csv")
    rows = read_reference(path)[1]
    feeding = np.array([[0.0, 1.0], [0.0, 0.0]])  # the upper level decays to the lower one
    worst, xs, ours, theirs = 0.0, [], [], []
    for r in rows:
        tau1, tau2, p1, p2 = (float(r[k]) for k in ("tau1_ps", "tau2_ps", "p1", "p2"))
        t = lifetimes.flight_times(float(r["distance_um"]), float(r["beta"]))
        upper, lower = lifetimes.unshifted([tau1, tau2], feeding, [p1, p2], t)
        worst = max(worst, float(_share_deviation(upper, float(r["unshifted_upper"]))),
                    float(_share_deviation(lower, float(r["unshifted_lower"]))))
        if (tau1, tau2, p1) == (2.0, 10.0, 1.0) and float(r["unshifted_lower"]) > 1e-12:
            xs.append(float(r["distance_um"]))
            ours.append(float(lower))
            theirs.append(float(r["unshifted_lower"]))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "distance (µm), τ = 2 ps feeding τ = 10 ps",
            "ylabel": "unshifted share, lower transition", "logx": True, "logy": True}
    return worst <= 1e-9, worst, 2 * len(rows), f"{len(rows)} distances × 2 transitions, worst {worst:.1e}", plot


def ru98_deorientation(ref) -> Optional[dict]:
    """Radeck et al.'s Table I deorientation parameters (α_k, Γ_k in 1/ps), the recoil speed (µm/ps) and τ(2⁺₁) (ps),
    from lifetimes/published/setups.toml and lifetimes.csv; None if the files are not there."""
    from ._toml import loads

    setups, taus = _lifetime_file(ref, "published/setups.toml"), _lifetime_file(ref, "published/lifetimes.csv")
    if setups is None or taus is None:
        return None
    s = loads(setups.read_text(encoding="utf-8"))["ru98_rdds_coulex"]

    def number(key: str, text: str) -> float:
        return float(re.search(rf"{key} ([0-9.]+)\(", text).group(1))

    tau = [float(r["tau_ps"]) for r in read_reference(taus)[1]
           if r["case"] == "ru98_rdds_coulex" and r["state"] == "2+1" and r["role"] == "result"]
    return {"alpha": {k: number(f"alpha_{k}", s["analysis"]) for k in (2, 4)},
            "gamma": {k: number(f"Gamma_{k}", s["analysis"]) for k in (2, 4)},
            "v_um_ps": number("v =", s["recoil_velocity"]), "tau_ps": tau[0]}


@check(RDDS, "literature", "Deorientation in vacuum, ⁹⁸Ru 2⁺₁ (Radeck et al., PRC 85, 014301 (2012), Table I: α₂ 0.18, "
       "Γ₂ 0.142/ps, α₄ 0.24, Γ₄ 0.314/ps at v = 13.49 µm/ps, τ = 8.36 ps; lifetimes/published/): the unshifted "
       "attenuation G_k(d) = α_k + (1 − α_k) e^(−Γ_k d/v) and the shifted one, Radeck's Eq. 7 by Gauss–Legendre "
       "quadrature, for d = 0–450 µm", "1e-8")
def deorientation_ru98(ref):
    from . import lifetimes

    p = ru98_deorientation(ref)
    if p is None:
        return _missing("published/setups.toml")
    d = np.linspace(0.0, 450.0, 91)
    t = d / p["v_um_ps"]
    tau, alpha, gamma = p["tau_ps"], p["alpha"], p["gamma"]
    unshifted = lifetimes.deorientation(t, alpha, gamma)
    shifted = lifetimes.deorientation_shifted(t, tau, alpha, gamma)
    x, w = np.polynomial.legendre.leggauss(64)
    worst, curves = 0.0, {}
    for k in (2, 4):
        g_closed = alpha[k] + (1 - alpha[k]) * np.exp(-gamma[k] * t)
        g_eq7 = np.ones_like(t)
        for i, tt in enumerate(t):
            if tt > 0:  # ∫₀ᵀ G_k e^(−t/τ) dt / ∫₀ᵀ e^(−t/τ) dt
                nodes, weights = 0.5 * tt * (x + 1), 0.5 * tt * w
                decay = np.exp(-nodes / tau)
                g_eq7[i] = np.sum(weights * (alpha[k] + (1 - alpha[k]) * np.exp(-gamma[k] * nodes)) * decay) / np.sum(
                    weights * decay)
        worst = max(worst, float(np.max(np.abs(unshifted[k] - g_closed))), float(np.max(np.abs(shifted[k] - g_eq7))))
        curves[k] = {"unshifted": np.asarray(unshifted[k]), "shifted": np.asarray(shifted[k]), "eq7": g_eq7}
    plot = {"x": d, "physim": curves[2]["shifted"], "reference": curves[2]["eq7"], "xlabel": "distance (µm)",
            "ylabel": "G̃₂ (shifted γ rays)", "curves": curves}
    detail = (f"G₂, G₄ and G̃₂, G̃₄ at 91 distances, worst {worst:.1e}; at 450 µm "
              f"G₂ = {curves[2]['unshifted'][-1]:.3f} (unshifted), G̃₂ = {curves[2]['shifted'][-1]:.3f} (shifted), G₄ = {curves[4]['unshifted'][-1]:.3f}, "
              f"G̃₄ = {curves[4]['shifted'][-1]:.3f}")
    return worst <= 1e-8, worst, 4 * len(d), detail, plot


@check(DSAM, "literature", "Stopping proportional to speed, v(t) = v₀ e^(−t/α): the attenuation F(τ) = α/(α + τ) "
       "exactly, and the distribution of v/v₀ at emission in 20 bins as the line shape at 0°; α = 0.2, 0.5, 1 ps, "
       "τ = 0.01–5 ps (lifetimes/analytic/dsam_exponential.csv)", "F 1e-6 relative; bin shares 2e-5")
def dsam_exponential_analytic(ref):
    from . import lifetimes
    from .gamma import doppler_energy

    path = _lifetime_file(ref, "analytic/dsam_exponential.csv")
    if path is None:
        return _missing("analytic/dsam_exponential.csv")
    rows = read_reference(path)[1]
    beta0, e0 = 0.04, 1000.0
    edges = doppler_energy(e0, beta0 * np.linspace(0.0, 1.0, 21), 1.0)
    kernels: dict = {}
    worst_f, worst_share, xs, ours, theirs = 0.0, 0.0, [], [], []
    for r in rows:
        alpha, tau = float(r["alpha_ps"]), float(r["tau_ps"])
        kernel = kernels.setdefault(alpha, lifetimes.exponential_kernel(beta0, alpha))
        f = lifetimes.dsam_attenuation(tau, kernel)
        worst_f = max(worst_f, abs(f / float(r["attenuation_F"]) - 1))
        shares = lifetimes.lineshape(e0, [tau], None, [1.0], kernel, 0.0, edges)[0]
        expected = np.array([float(v) for v in r["share_v_over_v0_in_20_bins"].split(";")])
        worst_share = max(worst_share, float(np.max(np.abs(shares - expected))))
        xs.append(tau / alpha)
        ours.append(f)
        theirs.append(float(r["attenuation_F"]))
    order = np.argsort(xs)
    worst = max(worst_f / 1e-6, worst_share / 2e-5)
    plot = {"x": np.array(xs)[order], "physim": np.array(ours)[order], "reference": np.array(theirs)[order],
            "xlabel": "τ / α", "ylabel": "attenuation F(τ)", "logx": True}
    detail = (f"{len(rows)} (α, τ): worst F {worst_f:.1e}, worst bin share {worst_share:.1e} ({worst:.0%} of the "
              "tolerance)")
    return worst <= 1.0, worst, len(rows) * 21, detail, plot


# -- the recoil transport against TRIM ----------------------------------------------------------------------------

#: Tolerances against TRIM with TRIM's own (SRIM 2013) stopping, so that only the scattering differs (as in
#: tests/python/test_nuclear_transport.py): TRIM adds electronic energy-loss straggling, which physim's Monte Carlo
#: leaves out (range straggling a little larger), and drops the distant collisions physim adds as a Gaussian turn
#: (its median angle a few per cent smaller). The transmitted fraction is absolute.
TRIM_TRANSPORT_TOLERANCE = {"projected_range_um": 0.03, "range_straggling_um": 0.15, "lateral_sigma_um": 0.06,
                            "transmitted": 0.02, "mean_energy_mev": 0.02, "rms_angle_deg": 0.06,
                            "median_angle_deg": 0.10, "angle_90_deg": 0.06}
#: Recoils followed per TRIM run (TRIM ran 5000).
TRIM_RECOILS = 10000


def _trim_runs(ref) -> Optional[dict]:
    path = _lifetime_file(ref, "trim/transport_trim.csv")
    if path is None:
        return None
    out: dict = {}
    for r in read_reference(path)[1]:
        out.setdefault((r["case"], r["run"]), {"row": r})[r["quantity"]] = float(r["value"])
    return out


def srim_stopping(ref, species: str, material):
    """physim's Stopping for ``species`` in ``material`` with SRIM 2013's electronic and nuclear stopping
    (lifetimes/trim/srim_stopping.csv) in place of its own, on its own grid and with its range integrated the same
    way: what TRIM slowed the recoils with. A copy: the cached stopping is not changed."""
    from .rates import stopping

    path = _lifetime_file(ref, "trim/srim_stopping.csv")
    rows = [r for r in read_reference(path)[1] if r["ion"] == species]
    e = np.array([float(r["energy_mev"]) for r in rows])
    se = np.array([float(r["electronic_mev_mg_cm2"]) for r in rows])
    sn = np.array([float(r["nuclear_mev_mg_cm2"]) for r in rows])
    st = copy.copy(stopping(species, material))
    ln = np.log(st._e)
    se_g = np.exp(np.interp(ln, np.log(e), np.log(se)))
    sn_g = np.exp(np.interp(ln, np.log(e), np.log(sn)))
    return _with_stopping(st, se_g, sn_g)


def _with_stopping(st, se_mev_mg_cm2: np.ndarray, sn_mev_mg_cm2: np.ndarray):
    """``st`` (a copy of a Stopping) with these electronic and nuclear stopping powers on its grid, and its range
    integrated from them as Stopping does."""
    st._se, st._sn = se_mev_mg_cm2 * 1e3, sn_mev_mg_cm2 * 1e3  # MeV cm²/g
    total = se_mev_mg_cm2 + sn_mev_mg_cm2
    f = st._e / total
    ln = np.log(st._e)
    st._range = np.r_[0.0, np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(ln))] + 2 * st._e[0] / total[0]
    return st


@contextmanager
def recoil_stopping(species: str, st):
    """Within the block, the recoil transport slows ``species`` with the Stopping ``st`` (its time-to-stop clocks
    recomputed for it); everything is put back afterwards."""
    from . import transport

    saved, clocks = transport.stopping, dict(transport._CLOCKS)
    transport.stopping = lambda ion, material: st if ion == species else saved(ion, material)
    transport._CLOCKS.clear()
    try:
        yield
    finally:
        transport.stopping = saved
        transport._CLOCKS.clear()
        transport._CLOCKS.update(clocks)


def _transport_vs_trim(ref, case: str) -> tuple:
    from . import transport
    from .rates import stopping

    runs = _trim_runs(ref)
    if runs is None or _lifetime_file(ref, "trim/srim_stopping.csv") is None:
        return _missing("trim/transport_trim.csv")
    first = runs[(case, "range")]["row"]
    species, material, energy = first["ion"], data.material(first["material"]), float(first["energy_mev"])
    mass = data.nuclide(species).nuclear_mass_mev
    um = float(material.density_g_cm3) / float(first["density_g_cm3"])  # physim's µm into TRIM's
    st = srim_stopping(ref, species, material)
    labels, ours, theirs, notes = [], [], [], []
    worst = 0.0
    z = np.array([0.0, 0.0, 1.0])
    for seed, run in enumerate(("range", "quarter", "half"), start=11):
        ref_run = runs[(case, run)]
        thickness = float(ref_run["row"]["thickness_mg_cm2"])
        with recoil_stopping(species, st):
            s = transport.slow([transport.Segment(material, thickness)], species, mass, np.full(TRIM_RECOILS, energy),
                               z, 0.0, np.full(TRIM_RECOILS, np.inf), scatter=True, rng=seed)
        if run == "range":
            pos = s.position_um[s.stopped] * um
            got = {"projected_range_um": pos[:, 2].mean(), "range_straggling_um": pos[:, 2].std(),
                   "lateral_sigma_um": math.sqrt(0.5 * float((pos[:, 0] ** 2 + pos[:, 1] ** 2).mean()))}
        else:
            out = s.segment == 1
            angle = np.degrees(np.arccos(np.clip(s.direction[out, 2], -1.0, 1.0)))
            got = {"transmitted": out.mean(), "mean_energy_mev": s.energy[out].mean(),
                   "rms_angle_deg": math.sqrt(float((angle**2).mean())), "median_angle_deg": float(np.median(angle)),
                   "angle_90_deg": float(np.percentile(angle, 90))}
        for key, value in got.items():
            want = ref_run[key]
            dev = abs(value - want) if key == "transmitted" else abs(value / want - 1)
            worst = max(worst, dev / TRIM_TRANSPORT_TOLERANCE[key])
            labels.append(f"{run}: {key.replace('_um', '').replace('_mev', '').replace('_deg', '')}")
            ours.append(float(value))
            theirs.append(want)
            if key != "transmitted":
                notes.append(f"{run} {key.rsplit('_', 1)[0].replace('_', ' ')} {value / want - 1:+.1%}")
    own = float(stopping(species, material).range(energy) / st.range(energy))
    detail = (f"{species} at β = {float(first['beta']):.3f} ({energy:.2f} MeV) into {first['material']}, "
              f"{TRIM_RECOILS} recoils per run, SRIM's stopping: " + ", ".join(notes)
              + f" (worst {worst:.0%} of its tolerance). With physim's own stopping the range is {own - 1:+.0%} "
              "(stopping_heavy_vs_data checks that stopping against measured data)")
    plot = {"x": labels, "physim": ours, "reference": theirs, "categorical": True}
    return worst <= 1.0, worst, len(labels), detail, plot


@check(RDDS, "tool", "TRIM (SRIM 2013), 5000 ions per run: ⁹⁸Ru at β = 4.5% (Radeck et al.'s recoil) into Cu, the "
       "plunger's stopper: projected range, range straggling and lateral spread, and the energy and angles after "
       "about a quarter and half of the range (lifetimes/trim/transport_trim.csv); physim's Monte Carlo slowing with "
       "SRIM's stopping (srim_stopping.csv), so that only the scattering differs",
       "range 3%, range straggling 15%, lateral spread 6%, transmitted fraction 0.02, mean energy 2%, rms and 90% "
       "angle 6%, median angle 10% (TRIM_TRANSPORT_TOLERANCE)")
def rdds_transport_vs_trim(ref):
    return _transport_vs_trim(ref, "ru98_cu")


@check(DSAM, "tool", "TRIM (SRIM 2013), 5000 ions per run: ¹⁴³Sm at β = 1.6% (Rajbanshi et al.'s residue) into Au, "
       "the DSAM backing: projected range, range straggling and lateral spread, and the energy and angles after "
       "about a quarter and half of the range (lifetimes/trim/transport_trim.csv); physim's Monte Carlo slowing with "
       "SRIM's stopping (srim_stopping.csv), so that only the scattering differs",
       "range 3%, range straggling 15%, lateral spread 6%, transmitted fraction 0.02, mean energy 2%, rms and 90% "
       "angle 6%, median angle 10% (TRIM_TRANSPORT_TOLERANCE)")
def dsam_transport_vs_trim(ref):
    return _transport_vs_trim(ref, "sm143_au")


# -- the lifetime analysis on simulated runs, and the published experiments (wave 3C) ----------------------------

#: NaPyTau against physim's DDCM in NaPyTau's procedure: τ(x) at each distance within this (relative).
NAPYTAU_TOLERANCE = 0.01


@check(RDDS, "tool", "NaPyTau (Cologne's DDCM, reimplemented; commit in lifetimes/napytau/ddcm_simulated.csv) on the "
       "peak areas of a simulated plunger run of the ⁹⁸Ru example (no deorientation, ring 41.5°, 11 distances; "
       "lifetime_reproduction.napytau_input): physim's DDCM in NaPyTau's procedure (one polynomial of degree 2, 3, 4 "
       "through the shifted points, unweighted, no time dilation) on the same areas, τ(x) at every distance",
       "τ within 1% of NaPyTau's (its central values; its uncertainties are not usable)")
def rdds_ddcm_vs_napytau(ref):
    from . import lifetime_reproduction as lr

    res = lr.ddcm_vs_napytau(ref)
    if res is None:
        return _missing("napytau/ddcm_simulated.csv")
    worst = res["worst"]
    deg = max(lr.NAPYTAU_DEGREES)
    per = res["per_degree"]
    points = sum(v["points"] for v in per.values())
    detail = (f"τ(x) at {points} distance×degree points, worst {worst:.1e} relative ("
              + ", ".join(f"degree {d}: {v['worst']:.1e}" for d, v in per.items())
              + f"); physim's own DDCM (joint fit, quadratic pieces) on the same areas {res['joint_tau_ps']:.2f} ± "
              f"{res['joint_stat_ps']:.2f} ps for the simulated {res['tau_ps']:.2f} ps")
    plot = {"x": res["distances_um"], "physim": per[deg]["physim"], "reference": per[deg]["napytau"],
            "xlabel": f"distance (µm), polynomial degree {deg}", "ylabel": "τ(x) (ps)", "logx": True}
    return worst <= NAPYTAU_TOLERANCE, worst / NAPYTAU_TOLERANCE, points, detail, plot


@check(RDDS, "literature", "Closure: τ recovered from simulated plunger runs of known τ: one level (the ⁹⁸Ru example's "
       "2⁺₁, 8.36 ps, no deorientation; the areas of lifetimes/napytau/ddcm_simulated.csv) by the decay-curve fit and "
       "the DDCM, and fed from above (the ⁵⁰Cr example's 2⁺₁, 12.2 ps, fed by its 4⁺₁; a seeded three-foil fusion "
       "run, lifetime_reproduction.cr50) by the three-foil DDCM in singles with the feeder's correction",
       "within the fit's uncertainty (2σ)")
def rdds_closure(ref):
    from . import lifetime_reproduction as lr

    one = lr.closure_one_level(ref)
    if one is None:
        return _missing("napytau/ddcm_simulated.csv")
    fed = lr.cr50(ref)
    if not fed.get("available"):
        return _missing("published/lifetimes.csv")
    rows = [("one level, decay curve", *one["decay_curve"], one["tau_ps"]),
            ("one level, DDCM", *one["ddcm"], one["tau_ps"]),
            ("fed from above, three-foil DDCM (singles)", *fed["singles"], fed["simulated_tau_ps"]["2+1"])]
    pulls = [(tau - true) / stat for _n, tau, stat, true in rows]
    worst = max(abs(x) for x in pulls) / 2
    detail = "; ".join(f"{name}: {tau:.2f} ± {stat:.2f} ps for {true:.2f} ({(tau - true) / stat:+.1f}σ)"
                       for name, tau, stat, true in rows)
    plot = {"x": one["distances_um"], "physim": one["model"], "reference": one["unshifted"],
            "xlabel": "distance (µm), one level", "ylabel": "unshifted counts", "logx": True}
    return worst <= 1.0, worst, len(rows), detail, plot


@check(RDDS, "literature", "Reproduced experiment: ⁹⁸Ru, RDDS after Coulex below the barrier (Radeck et al., PRC 85, "
       "014301 (2012)): τ(2⁺₁) = 8.36(29) ps with the deorientation correction (lifetimes/published/). The example "
       "coulex_ru98_plunger simulated at its 11 distances (particle–γ with the ²⁴Mg recoils), analysed as the paper: "
       "β from the shifted centroids, decay-curve fits per ring corrected with Table I's deorientation "
       "(lifetime_reproduction.ru98). The 4⁺₁ (2.31(16) ps) and 2⁺₂ (1.7(2) ps) are not reproduced: their two-step "
       "excitation needs matrix elements the paper does not give", "within the combined uncertainty (2σ)")
def rdds_ru98_reproduction(ref):
    from . import lifetime_reproduction as lr

    res = lr.ru98(ref)
    if not res.get("available"):
        return _missing("published/lifetimes.csv")
    rings = sorted(res["rings"])
    r0 = res["rings"][rings[0]]
    plot = {"x": res["distances_um"], "physim": r0["model"], "reference": list(np.array(r0["unshifted"])
                                                                            / np.array(r0["w_unshifted"])),
            "xlabel": f"distance (µm), ring {rings[0]:g}°", "ylabel": "unshifted counts / W_U (deorientation)",
            "logx": True}
    return res["passed"], res["worst"], len(res["distances_um"]) * len(rings), res["detail"], plot


@check(RDDS, "literature", "Reproduced experiment: ⁵⁰Cr, RDDS after fusion-evaporation with a three-foil plunger "
       "(Beckers et al., NIM A 1042, 167418 (2022)): τ(2⁺₁) = 12.2(4) ps, τ(4⁺₁) = 2.4(3) ps by γγ DDCM, against "
       "the Coulex, RDDS and DSAM values for the same states (lifetimes/published/). The example fusion_cr50_plunger "
       "simulated at nine pairs of distances, analysed as the paper: direct gates as its Eq. 13, the degrader and "
       "stopper parts apart (lifetime_reproduction.cr50)", "within the 2σ spread of the published values")
def rdds_cr50_reproduction(ref):
    from . import lifetime_reproduction as lr

    res = lr.cr50(ref)
    if not res.get("available"):
        return _missing("published/lifetimes.csv")
    t1 = res["per_distance"]["t1"]
    plot = {"x": [p["distance_um"] for p in t1], "physim": [p["tau_ps"] for p in t1],
            "reference": [res["states"]["2+1"]["paper_tau_ps"]] * len(t1),
            "xlabel": "x₁ (µm), 2⁺₁, t₁ part", "ylabel": "τ(x) (ps)", "logx": True}
    return res["passed"], res["worst"], 2, res["detail"], plot


#: The papers' plunger distances (µm) for their studied level (lifetimes/published/setups.toml): (label, example,
#: shortest, longest). ⁹⁸Ru: 3–450 µm for τ(2⁺₁) = 8.36 ps (Radeck et al. 2012); ⁵⁰Cr: 17–750 µm to the degrader for
#: τ(2⁺₁) = 12.2 ps (Beckers et al. 2022).
PAPER_DISTANCES = (("⁹⁸Ru 2⁺₁", "coulex_ru98_plunger", 3.0, 450.0), ("⁵⁰Cr 2⁺₁", "fusion_cr50_plunger", 17.0, 750.0))
#: How far (a factor) the paper's range may fall short of the planner's recommended distances at either end.
PLANNER_RANGE_FACTOR = 1.5


@check(RDDS, "literature", "The planner's choices against the papers': the plunger distances Planner.lifetimes "
       "recommends for the studied level of the examples made after the papers (coulex_ru98_plunger, "
       "fusion_cr50_plunger: log-spaced over the sensitive range, 0.18 to 5.1 decay lengths at physim's recoil speed) "
       "against the distances the papers chose: ⁹⁸Ru 3–450 µm for τ(2⁺₁) = 8.36 ps (Radeck et al. 2012), ⁵⁰Cr "
       "17–750 µm to the degrader for τ(2⁺₁) = 12.2 ps (Beckers et al. 2022) (lifetimes/published/setups.toml)",
       "the paper's range covers the planner's recommended distances to within a factor 1.5 at each end: the paper's "
       "shortest at most 1.5 × the planner's shortest, its longest at least the planner's longest / 1.5 (worst: the "
       "larger of the two ratios paper/planner at the short end and planner/paper at the long end). The factor allows "
       "physim's recoil speeds differing from the papers' (by up to 25%) and the factor-2 definition of the sensitive range; shorter "
       "distances a paper adds for its shorter-lived levels are not held against the planner")
def rdds_planner_distances(ref):
    from .planner import Planner

    worst, rows, plot_x, plot_p, plot_r = 0.0, [], [], [], []
    for label, example, lo, hi in PAPER_DISTANCES:
        out = Planner.example(example).lifetimes()
        if not out.get("available") or not out["levels"] or not out["levels"][0]["recommended_distances_um"]:
            return None, math.nan, 0, f"{example}: no recommended distances ({out.get('reason')})", {}
        row = out["levels"][0]
        rec = row["recommended_distances_um"]
        short, long_ = lo / min(rec), max(rec) / hi
        worst = max(worst, short, long_)
        rows.append(f"{label} (planner's {row['level']}, τ = {row['tau_ps']:.3g} ps, β = {100 * row['beta']:.2f}%): "
                    f"planner {min(rec):.3g}–{max(rec):.3g} µm, paper {lo:g}–{hi:g} µm (short end {short:.2f}, long "
                    f"end {long_:.2f})")
        plot_x += [f"{label} shortest", f"{label} longest"]
        plot_p += [min(rec), max(rec)]
        plot_r += [lo, hi]
    plot = {"x": plot_x, "physim": plot_p, "reference": plot_r, "categorical": True}
    return worst <= PLANNER_RANGE_FACTOR, worst, len(rows), "; ".join(rows), plot


@check(DSAM, "literature", "Closure: τ recovered by the line-shape fit from simulated DSAM spectra of known τ (with "
       "side feeding): the ¹⁴³Sm 444.8 and 527.5 keV lines (the paper's Fig. 2) of lifetime_reproduction.sm143, "
       "fitted at 148° and 123° with a kernel of separate recoils slowed with the same stopping",
       "within the fit's uncertainty (2σ)")
def dsam_closure(ref):
    from . import lifetime_reproduction as lr

    res = lr.dsam_closure(ref)
    if not res.get("available"):
        return _missing("published/lifetimes.csv")
    lv = res["levels"]
    detail = "; ".join(f"{lab} ({v['e_gamma_kev']:.1f} keV): {v['tau_ps']:.3f} ± {v['stat_ps']:.3f} ps for "
                       f"{v['true_ps']:.3f} ({v['pull']:+.1f}σ)" for lab, v in lv.items())
    plot = {"x": list(lv), "physim": [v["tau_ps"] for v in lv.values()],
            "reference": [v["true_ps"] for v in lv.values()], "categorical": True}
    return res["worst"] <= 1.0, res["worst"], len(lv), detail, plot


@check(DSAM, "literature", "Reproduced experiment: ¹⁴³Sm, DSAM after fusion-evaporation on a 13 mg/cm² Au backing "
       "(Rajbanshi et al., PLB 782, 143 (2018)): seven lifetimes 0.41–1.02 ps and two upper limits in two dipole "
       "bands (lifetimes/published/). The example fusion_sm143_dsam (the paper's level scheme, branchings and "
       "side-feeding times) simulated with physim's stopping, analysed as the paper: line-shape fits from the top "
       "of each band down with SRIM's stopping, at 148° and 123° (lifetime_reproduction.sm143)",
       "within the published uncertainty plus SRIM's ±5% stopping (2σ)")
def dsam_sm143_reproduction(ref):
    from . import lifetime_reproduction as lr

    res = lr.sm143(ref)
    if not res.get("available"):
        return _missing("published/lifetimes.csv")
    lv = res["levels"]
    order = res["order"]
    plot = {"x": [f"{lv[k]['e_gamma_kev']:.0f}" for k in order], "physim": [lv[k]["tau_ps"] for k in order],
            "reference": [lv[k]["paper_tau_ps"] for k in order], "categorical": True}
    return res["passed"], res["worst"], len(order), res["detail"], plot


#: The paper's DSAM backing (lifetimes/published/setups.toml): ¹⁴³Sm on 13 mg/cm² Au (Rajbanshi et al. 2018).
PAPER_BACKING = ("fusion_sm143_dsam", "Au", 13.0)
#: The share of the best candidate's sensitivity the paper's backing must reach (the planner's BACKING_GOOD).
BACKING_SHARE = 0.8


@check(DSAM, "literature", "The planner's choices against the papers': the DSAM backing advice of Planner.lifetimes "
       "(the usual backings behind the target, ranked per level by the sensitivity |dF/d ln τ| of F(τ) on each "
       "one's slowing-down kernel) for the levels of the example made after the paper (fusion_sm143_dsam: τ "
       "0.41–1.02 ps), against the backing the paper chose: 13 mg/cm² Au (Rajbanshi et al. 2018; "
       "lifetimes/published/setups.toml)",
       "for every level advised on, Au's sensitivity at least 80% of the best candidate's (worst: the largest "
       "1 − S_Au/S_best, at most 0.2), and the paper's 13 mg/cm² at least the planner's thickness to stop ¹⁴³Sm in "
       "Au")
def dsam_planner_backing(ref):
    from .planner import Planner

    example, material, thickness = PAPER_BACKING
    out = Planner.example(example).lifetimes()
    rows = [r for r in out.get("levels", []) if r.get("backing_advice")]
    if not rows:
        return None, math.nan, 0, f"{example}: no backing advice ({out.get('reason')})", {}
    worst, parts, stops, x, got, best = 0.0, [], True, [], [], []
    for r in rows:
        adv = r["backing_advice"]
        mine = next((c for c in adv if c["material"] == material), None)
        if mine is None:
            return False, math.inf, len(rows), f"{material} is not among the planner's candidates", {}
        worst = max(worst, 1.0 - mine["relative_sensitivity"])
        stops = stops and thickness >= mine["thickness_to_stop_mg_cm2"]
        rank = adv.index(mine) + 1
        parts.append(f"{r['level']} ({r['tau_ps']:.2f} ps): {material} {100 * mine['relative_sensitivity']:.0f}% of "
                     f"the best ({adv[0]['material']}), rank {rank} of {len(adv)}")
        x.append(f"{r['gamma_kev']:.0f} keV")
        got.append(mine["sensitivity"])
        best.append(adv[0]["sensitivity"])
    to_stop = next(c for c in rows[0]["backing_advice"] if c["material"] == material)["thickness_to_stop_mg_cm2"]
    detail = (f"thickness to stop ¹⁴³Sm in {material}: {to_stop:.2f} mg/cm² (paper {thickness:g}); "
              + "; ".join(parts))
    plot = {"x": x, "physim": got, "reference": best, "categorical": True}
    return worst <= 1.0 - BACKING_SHARE and stops, worst, len(rows), detail, plot


# -- DSAM line shapes against APCAD and AhKin ---------------------------------------------------------------------


def _fold(edges: np.ndarray, shape: np.ndarray, fwhm_kev: float) -> np.ndarray:
    """A line shape on uniform bins folded with a Gaussian of this FWHM (the detector's response)."""
    step = edges[1] - edges[0]
    sigma = fwhm_kev / math.sqrt(8 * math.log(2))
    half = int(math.ceil(6 * sigma / step))
    x = np.arange(-half, half + 1) * step  # symmetric, odd length: no shift in "same" mode
    g = np.exp(-0.5 * (x / sigma) ** 2)
    return np.convolve(shape, g / g.sum(), mode="same")


def _rebin(fine_edges: np.ndarray, contents: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Contents on ``fine_edges`` put into the bins ``edges`` (linear within a fine bin)."""
    return np.diff(np.interp(edges, fine_edges, np.r_[0.0, np.cumsum(contents)]))


def _moments(x: np.ndarray, y: np.ndarray) -> tuple:
    """(centroid, rms width) of the counts y at x."""
    c = float(np.sum(y * x) / np.sum(y))
    return c, math.sqrt(abs(float(np.sum(y * (x - c) ** 2) / np.sum(y))))


#: APCAD's example (lifetimes/apcad/: Macro_StopSim_example.mac, README.txt; StopSim's source settles the rest): a
#: 400 MeV ¹²⁰Sn beam on 6.637168 µm C (1.5 mg/cm² at StopSim's 2.26 g/cm³) on 11.210762 µm Cu (10.0 mg/cm² at
#: 8.92 g/cm³). StopSim tracks the beam-like ion (TrackBeamLike), excited to 1000 keV in the carbon only
#: (KillInLayer2), with the CM angle of the projectile drawn from XSec_example.dat (each entry dσ integrated over its
#: 1° bin, linear in the beam energy); a nucleus leaving the target keeps its exit velocity. The spectra: τ = 1 ps,
#: FWHM 3 keV, background 400 − 0.252 E + 5.6·10⁻⁵ E² (E in keV), 29 angles 20°–160°, 1 keV bins.
APCAD_SETUP = {"beam": "120Sn", "target": "12C", "beam_energy_mev": 400.0, "target_mg_cm2": 1.5,
               "backing_mg_cm2": 10.0, "e_gamma_kev": 1000.0, "tau_ps": 1.0, "fwhm_kev": 3.0,
               "background": (400.0, -0.252, 5.6e-5)}
#: The first-order centroid shift and the rms width, physim over APCAD, within this (about what a ±10% change of
#: the stopping powers does to the shift).
APCAD_TOLERANCE = 0.10
#: Sample nuclei in physim's kernel for the APCAD and AhKin comparisons.
LINESHAPE_RECOILS = 4000


def _apcad_kernel(ref):
    """physim's DSAM kernel for APCAD's example: the excited ¹²⁰Sn from APCAD's cross-section table and physim's
    kinematics, slowed through the rest of the carbon and the copper by physim's transport (with scattering)."""
    from . import transport
    from .kinematics import TwoBody
    from .rates import stopping

    s = APCAD_SETUP
    xs = np.loadtxt(_lifetime_file(ref, "apcad/XSec_example.dat"))
    energies, angles = np.unique(xs[:, 0]), np.unique(xs[:, 1])
    table = np.zeros((len(energies), len(angles)))
    table[np.searchsorted(energies, xs[:, 0]), np.searchsorted(angles, xs[:, 1])] = xs[:, 2]
    rng = np.random.default_rng(5)
    carbon, copper = data.material("C"), data.material("Cu")
    beam_stopping = stopping(s["beam"], carbon)
    grid = np.linspace(0.0, s["target_mg_cm2"], 61)
    beam_e = np.array([s["beam_energy_mev"]] + [beam_stopping.energy_after(s["beam_energy_mev"], f"{d} mg/cm2")
                                                for d in grid[1:]])
    depth = rng.uniform(0.0, s["target_mg_cm2"], 8 * LINESHAPE_RECOILS)
    e_at = np.interp(depth, grid, beam_e)
    total = np.interp(e_at, energies, table.sum(axis=1))
    keep = rng.uniform(0.0, total.max(), len(total)) < total  # where the reactions happen: ∝ σ(E)
    depth, e_at = depth[keep][:2 * LINESHAPE_RECOILS], e_at[keep][:2 * LINESHAPE_RECOILS]
    theta = np.empty(len(depth))
    for i, e in enumerate(e_at):
        j = int(np.clip(np.searchsorted(energies, e) - 1, 0, len(energies) - 2))
        w = (e - energies[j]) / (energies[j + 1] - energies[j])
        row = (1 - w) * table[j] + w * table[j + 1]
        k = rng.choice(len(angles), p=row / row.sum())
        theta[i] = angles[k] - rng.uniform(0.0, 1.0)  # the entry at a covers (a − 1°, a]
    theta = np.clip(theta, 1e-3, 180.0)
    energy, theta_lab = np.empty(len(theta)), np.empty(len(theta))
    rounded = np.round(e_at)
    for e in np.unique(rounded):
        m = rounded == e
        p = TwoBody(s["beam"], s["target"], float(e), excitation_mev=s["e_gamma_kev"] / 1e3,
                    excite="ejectile").at_cm(theta[m], "ejectile")
        energy[m], theta_lab[m] = p.energy, p.theta_lab
    phi = rng.uniform(0.0, 2 * math.pi, len(theta))
    sin = np.sin(np.radians(theta_lab))
    direction = np.c_[sin * np.cos(phi), sin * np.sin(phi), np.cos(np.radians(theta_lab))]
    segments = [transport.Segment(carbon, s["target_mg_cm2"], name="target"),
                transport.Segment(copper, s["backing_mg_cm2"], name="backing")]
    mass = data.nuclide(s["beam"]).nuclear_mass_mev + s["e_gamma_kev"] / 1e3
    kern = transport.kernel(segments, s["beam"], mass, energy, direction, depth, np.linspace(0.0, 5.0, 201), rng=6,
                            n=LINESHAPE_RECOILS)
    return kern, float(stopping(s["beam"], copper).range(float(np.median(energy))))


@check(DSAM, "tool", "APCAD (Stahl, Leske, Lettmann, Pietralla; StopSim with Geant4), its OneShape example: a 1000 keV "
       "line of τ = 1 ps in the excited ¹²⁰Sn after 400 MeV ¹²⁰Sn on 1.5 mg/cm² C on 10 mg/cm² Cu, 29 detector angles "
       "20°–160° (lifetimes/apcad/). physim's kernel (its kinematics with APCAD's cross-section table, its stopping "
       "and scattering) through lifetimes.lineshape, folded with the 3 keV FWHM, against the spectra less their stated "
       "background: shape only, not a lifetime fit",
       "first-order centroid shift (half the difference between θ and 180° − θ, summed over 20°–85°) and the mean rms "
       "width, physim/APCAD within 10% (APCAD_TOLERANCE); χ²/bin reported")
def dsam_lineshape_vs_apcad(ref):
    from . import lifetimes

    path = _lifetime_file(ref, "apcad/OneShape_example.csv")
    if path is None or _lifetime_file(ref, "apcad/XSec_example.dat") is None:
        return _missing("apcad/OneShape_example.csv")
    s = APCAD_SETUP
    rows = read_reference(path)[1]
    angle = np.array([float(r["angle_deg"]) for r in rows])
    energy = np.array([float(r["energy_kev"]) for r in rows])
    counts = np.array([float(r["counts"]) for r in rows])
    b0, b1, b2 = s["background"]
    background = b0 + b1 * energy + b2 * energy**2
    kern, range_mg_cm2 = _apcad_kernel(ref)
    left = float(np.mean(kern.segment[-1] == 2))
    fine = np.arange(600.0, 1400.0001, 0.5)
    # The window at each angle, the same for both, from the kinematics alone: the excited ¹²⁰Sn moves at most at the
    # beam's speed and within asin(m(¹²C)/m(¹²⁰Sn)) = 5.7° of the beam, so its γ rays lie between E₀ (at rest) and
    # E₀ √(1 − β²)/(1 − β cos ψ) for ψ within 5.7° of the detector's angle; three FWHM more on each side.
    mass = data.nuclide(s["beam"]).nuclear_mass_mev
    beta = math.sqrt(1 - (mass / (mass + s["beam_energy_mev"])) ** 2)
    cone = math.degrees(math.asin(data.nuclide(s["target"]).nuclear_mass_mev / mass))
    shift_p, shift_d, width_p, width_d, chi, dof, off_chi, off_n = {}, {}, {}, {}, 0.0, 0, 0.0, 0
    for a in np.unique(angle):
        m = angle == a
        e, y, bg = energy[m], counts[m], background[m]
        edges = np.r_[e - 0.5, e[-1] + 0.5]
        shape = lifetimes.lineshape(s["e_gamma_kev"], [s["tau_ps"]], None, [1.0], kern, float(a), fine)[0]
        model = _rebin(fine, _fold(fine, shape, s["fwhm_kev"]), edges)
        psi = np.radians(np.clip(np.linspace(a - cone, a + cone, 25), 0.0, 180.0))
        ends = np.r_[1.0, math.sqrt(1 - beta**2) / (1 - beta * np.cos(psi))] * s["e_gamma_kev"]
        win = (e >= ends.min() - 3 * s["fwhm_kev"]) & (e <= ends.max() + 3 * s["fwhm_kev"])
        shift_p[a], width_p[a] = _moments(e[win], model[win])
        shift_d[a], width_d[a] = _moments(e[win], (y - bg)[win])
        shift_p[a] -= s["e_gamma_kev"]
        shift_d[a] -= s["e_gamma_kev"]
        w = 1.0 / np.maximum(y, 1.0)
        amp = float(np.sum(((y - bg) * model * w)[win]) / np.sum((model**2 * w)[win]))
        chi += float(np.sum(((y - bg - amp * model) ** 2 * w)[win]))
        dof += int(win.sum()) - 1
        off = ~win & (y > 0)  # (the spectra's last bin, at 1500 keV, is empty)
        off_chi += float(np.sum(((y - bg) ** 2 * w)[off]))
        off_n += int(off.sum())
    pairs = [(a, 180.0 - a) for a in shift_p if a < 90 and (180.0 - a) in shift_p]
    first_p = sum(shift_p[a] - shift_p[b] for a, b in pairs) / 2
    first_d = sum(shift_d[a] - shift_d[b] for a, b in pairs) / 2
    ratio_shift = first_p / first_d
    ratios_w = np.array([width_p[a] / width_d[a] for a in width_p])
    ratio_width = float(ratios_w.mean())
    worst = max(abs(ratio_shift - 1), abs(ratio_width - 1)) / APCAD_TOLERANCE
    angles = sorted(shift_p)
    detail = (f"first-order shift physim/APCAD {ratio_shift:.3f} (at 20° {shift_p[20.0]:+.1f} against "
              f"{shift_d[20.0]:+.1f} keV, at 160° {shift_p[160.0]:+.1f} against {shift_d[160.0]:+.1f} keV); rms width "
              f"physim/APCAD mean {ratio_width:.3f} ({ratios_w.min():.2f}–{ratios_w.max():.2f} by angle), both within "
              f"each angle's kinematic window; χ²/bin {chi / dof:.1f} there (amplitude free per angle, the stated "
              f"background; the spectra's χ²/bin outside {off_chi / off_n:.2f}). In physim {left:.0%} of the ¹²⁰Sn leave the "
              f"10 mg/cm² Cu (range at the median energy {range_mg_cm2:.1f} mg/cm²) and fly on at their exit speed; "
              "APCAD's line has more of the stopped peak and its flying component lower in energy: StopSim slows the "
              "ions more, with Geant4's ion stopping (G4IonParametrisedLossModel, scaled above Z = 26)")
    plot = {"x": angles, "physim": [shift_p[a] for a in angles], "reference": [shift_d[a] for a in angles],
            "xlabel": "detector angle (deg)", "ylabel": "centroid shift (keV)"}
    return worst <= 1.0, worst, len(angles), detail, plot


def _numbers(path, labels: bool = False):
    """The values of an AhKin input file: the lines' first fields (comments after # dropped), or {label: value}."""
    out = [] if not labels else {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        body = line.split("#", 1)[0].split()
        if not body:
            continue
        if labels:
            if len(body) >= 2:
                out[body[0]] = body[1]
        else:
            out.append(body)
    return out


def ahkin_case(ref) -> Optional[dict]:
    """AhKin's ⁸³Y example (lifetimes/ahkin/y83_in_ta/): the cascade, detector, resolution, background, the measured
    line and AhKin's best line and its 1σ edges; None without the files."""
    folder = _lifetime_file(ref, "ahkin/y83_in_ta")
    names = ("in_AhKin.dat", "in_casc.dat", "in_stp.dat", "pcos.dat", "EL1192GFAA53.txt", "ech0.131308.dat",
             "ech0.0986328.dat", "ech0.165235.dat", "Y83_in_Ta.txt")
    if folder is None or any(not (folder / n).exists() for n in names):
        return None
    general = _numbers(folder / "in_AhKin.dat")
    e0, (cal0, cal1), reso, background = float(general[0][0]), map(float, general[1][:2]), general[3], general[5]
    casc = [line.split() for line in (folder / "in_casc.dat").read_text(encoding="utf-8").splitlines() if line.split()]
    n = int(casc[0][0])
    taus = [float(v) for v in casc[1]]
    feeding = np.zeros((n, n))  # AhKin's B[i][j]: from j to i; physim's rows are "from"
    for i in range(n):
        for j, v in enumerate(casc[2 + i][:i]):
            feeding[j, i] = float(v)
    stp = _numbers(folder / "in_stp.dat", labels=True)
    measured = np.loadtxt(folder / "EL1192GFAA53.txt")
    cos = np.loadtxt(folder / "pcos.dat", comments="#")
    return {"e0": e0, "cal": (cal0, cal1), "fwhm": float(reso[0]) + float(reso[1]) * e0 + float(reso[2]) * e0**2,
            "background": float(background[0]), "taus": taus, "feeding": feeding,
            "populations": np.array([float(v) for v in casc[2 + n]]), "e_p": float(stp["E_p"]) / 1e3,
            "fsig": float(stp["fsig"]), "residue": data.nuclide((int(stp["Zr"]), int(stp["Ar"]))).name,
            "a_beam": int(stp["Ap"]), "a_target": int(stp["At"]), "cos": cos[cos[:, 1] > 0],
            "channel": measured[:, 0], "measured": measured[:, 1],
            "best": np.loadtxt(folder / "ech0.131308.dat")[:, 1], "tau": 0.131308,
            "edges": {0.0986328: np.loadtxt(folder / "ech0.0986328.dat")[:, 1],
                      0.165235: np.loadtxt(folder / "ech0.165235.dat")[:, 1]},
            "stopping_table": np.loadtxt(folder / "Y83_in_Ta.txt")[1:]}


def _ahkin_line(case: dict, ahkin_stopping: bool) -> np.ndarray:
    """physim's line for AhKin's example (counts per channel, with AhKin's flat background, scaled to the measured
    counts): the residues start at the compound nucleus's speed with AhKin's Gaussian spread, along the beam, and
    slow in Ta by physim's transport (with scattering) and either AhKin's stopping table or physim's own; the
    cascade's times through lifetimes.lineshape at five detector angles weighted by AhKin's cosine distribution,
    folded with its resolution."""
    from . import lifetimes, transport
    from .rates import stopping

    beam, target = "32S", "58Ni"  # the README's ³²S + ⁵⁸Ni at 135 MeV; in_stp.dat gives the mass numbers
    assert (data.nuclide(beam).A, data.nuclide(target).A) == (case["a_beam"], case["a_target"])
    species = case["residue"]
    mass = data.nuclide(species).nuclear_mass_mev
    m_b, m_t = data.nuclide(beam).nuclear_mass_mev, data.nuclide(target).nuclear_mass_mev
    p = math.sqrt(case["e_p"] * (case["e_p"] + 2 * m_b))
    beta_cn = p / (case["e_p"] + m_b + m_t)
    e_r = mass * (1 / math.sqrt(1 - beta_cn**2) - 1)
    rng = np.random.default_rng(1)
    energies = e_r * (1 + case["fsig"] * rng.standard_normal(2 * LINESHAPE_RECOILS))
    energies = energies[energies > 0.05 * e_r]
    tantalum = data.material("Ta")
    st = copy.copy(stopping(species, tantalum))
    if ahkin_stopping:  # AhKin's total stopping (keV/µm), physim's nuclear part inside it
        t = case["stopping_table"]
        total = np.exp(np.interp(np.log(st._e), np.log(t[:, 0] / 1e3), np.log(t[:, 1] / (100 * float(
            tantalum.density_g_cm3)))))
        nuclear = np.minimum(np.asarray(st._sn) * 1e-3, 0.95 * total)
        st = _with_stopping(st, total - nuclear, nuclear)
    with recoil_stopping(species, st):
        kern = transport.kernel([transport.Segment(tantalum, 10.0, name="backing")], species, mass, energies,
                                np.array([0.0, 0.0, 1.0]), 0.0, np.linspace(0.0, 1.6, 321), rng=2,
                                n=LINESHAPE_RECOILS)
    cos, w = case["cos"][:, 0], case["cos"][:, 1]
    group = np.minimum((np.cumsum(w) / w.sum() * 5 - 1e-9).astype(int), 4)  # five groups of equal weight
    cal0, cal1 = case["cal"]
    edges = cal0 + cal1 * np.r_[case["channel"] - 0.5, case["channel"][-1] + 0.5]
    fine = np.arange(edges[0] - 10.0, edges[-1] + 10.0, 0.2)
    taus = list(case["taus"])
    taus[-1] = case["tau"]  # the level of interest is the last; AhKin's τ
    shape = np.zeros(len(fine) - 1)
    for g in range(5):
        m = group == g
        theta = math.degrees(math.acos(float(np.average(cos[m], weights=w[m]))))
        shape += w[m].sum() * lifetimes.lineshape(case["e0"], taus, case["feeding"], case["populations"], kern, theta,
                                                  fine)[-1]
    line = _rebin(fine, _fold(fine, shape / w.sum(), case["fwhm"]), edges)
    signal = case["measured"] - case["background"]
    return case["background"] + line * signal.sum() / line.sum()


@check(DSAM, "tool", "AhKin (Garzón, Rodríguez, Cristancho, Tao, CPC 2020), its example: the 1193 keV line of ⁸³Y "
       "(³²S + ⁵⁸Ni at 135 MeV, residues slowing in Ta) at 52.8° and its best τ = 0.1313 ps "
       "(lifetimes/ahkin/y83_in_ta/). physim with AhKin's inputs (its stopping table, ten-level cascade with side "
       "feeding, detector cosine distribution, 10% energy spread, resolution, background) and its own transport and "
       "line shape, against AhKin's best line, its lines at the ends of its 1σ interval (0.0986, 0.1652 ps) and the "
       "measured line: shape only, not a lifetime fit",
       "centroid shift inside AhKin's 1σ band (the shifts of its two edge lines); rms width within 10% of AhKin's; "
       "χ²/bin against the measured line no larger than that of AhKin's edge lines (its own χ²_min + 1)")
def dsam_lineshape_vs_ahkin(ref):
    case = ahkin_case(ref)
    if case is None:
        return _missing("ahkin/y83_in_ta/")
    cal0, cal1 = case["cal"]
    x = cal0 + cal1 * case["channel"]
    bg, measured = case["background"], case["measured"]
    err2 = np.maximum(measured, 1.0)

    def chi(line):
        return float(np.sum((measured - line) ** 2 / err2) / (len(measured) - 1))

    def shift_width(line):
        c, w = _moments(x, line - bg)
        return c - case["e0"], w

    ours = _ahkin_line(case, ahkin_stopping=True)
    own = _ahkin_line(case, ahkin_stopping=False)
    s_p, w_p = shift_width(ours)
    s_a, w_a = shift_width(case["best"])
    band = sorted(shift_width(line)[0] for line in case["edges"].values())
    chi_band = max(chi(line) for line in case["edges"].values())
    in_band = band[0] <= s_p <= band[1]
    half = (band[1] - band[0]) / 2
    worst = max(abs(s_p - (band[0] + band[1]) / 2) / half, abs(w_p / w_a - 1) / 0.10, chi(ours) / chi_band)
    ok = in_band and abs(w_p / w_a - 1) <= 0.10 and chi(ours) <= chi_band
    detail = (f"centroid shift {s_p:+.2f} keV against AhKin's {s_a:+.2f} (1σ band {band[0]:+.2f} to {band[1]:+.2f}): "
              f"{'inside' if in_band else 'outside'}; rms width {w_p:.2f} against {w_a:.2f} keV "
              f"({w_p / w_a - 1:+.0%}); χ²/bin against the measured line {chi(ours):.1f} against AhKin's "
              f"{chi(case['best']):.1f} (edges {chi_band:.1f}). physim's scattering (checked against TRIM) turns the "
              "residues, so their sideways motion broadens the line; AhKin follows the velocity along the beam only, "
              "reduced by Blaugrund's mean cosine. With physim's own stopping (above AhKin's table at these "
              f"energies): shift {shift_width(own)[0]:+.2f} keV, width {shift_width(own)[1]:.2f} keV, χ²/bin "
              f"{chi(own):.1f}")
    plot = {"x": x, "physim": ours, "reference": measured, "ahkin": case["best"], "xlabel": "energy (keV)",
            "ylabel": "counts per channel"}
    return ok, worst, len(measured), detail, plot


# ---------------------------------------------------------------------------------------------------------------
# Running and reporting


def run(capability: Optional[str] = None, kind: Optional[str] = None, reference: Optional[Path] = None,
        names=None) -> list:
    """Run the checks (all, or those of one capability and/or kind, or the checks named in ``names``) and return
    their :class:`Result`."""
    ref = reference_dir() if reference is None else Path(reference)
    out = []
    for c in CHECKS:
        if (capability is not None and c.capability != capability) or (kind is not None and c.kind != kind):
            continue
        if names is not None and c.name not in names:
            continue
        passed, worst, points, detail, plot = c.run(ref)
        status = "pending" if passed is None else ("pass" if passed else "fail")
        out.append(Result(c.capability, c.kind, c.name, c.reference, c.tolerance, status, float(worst), points,
                          detail, plot))
    return out


def report(results: Optional[list] = None) -> str:
    """A Markdown table of check results."""
    results = run() if results is None else results
    mark = {"pass": "✅", "fail": "❌", "pending": "🟡"}
    lines = ["| Capability | Kind | Check | Tolerance | Result |", "|---|---|---|---|---|"]
    for r in results:
        lines.append(f"| {r.capability} | {r.kind} | {r.reference} | {r.tolerance} | {mark[r.status]} {r.detail} |")
    return "\n".join(lines)


def plot(result: Result, ax=None):
    """Draw a check's comparison: physim against the reference, or the residual. Returns the axes."""
    from ..plot import _plt, tidy

    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 3.2))
    p = result.plot
    if not p:
        ax.text(0.5, 0.5, f"{result.name}: {result.status}", transform=ax.transAxes, ha="center")
        return ax
    x = np.arange(len(p["x"])) if p.get("categorical") else np.asarray(p["x"], dtype=float)
    if p.get("categorical"):  # different quantities and units: compare as ratios
        ax.axhline(1.0, color="#888888", lw=0.8)
        ax.plot(x, np.asarray(p["physim"], dtype=float) / np.asarray(p["reference"], dtype=float), "o", ms=4,
                color="#d62728")
        p = dict(p, ylabel="physim / reference", residual=True)
    elif p.get("residual"):
        ax.axhline(0.0, color="#888888", lw=0.8)
        ax.plot(x, p["physim"], "o", ms=4, color="#1f77b4")
    else:
        ax.plot(x, p["reference"], "s", ms=5, mfc="none", color="#444444", label="reference")
        ax.plot(x, p["physim"], "o", ms=3, color="#d62728", label="physim")
    if p.get("categorical"):
        ax.set_xticks(x)
        ax.set_xticklabels(p["x"], rotation=20, fontsize=7)
    if p.get("logx"):
        ax.set_xscale("log")
    if p.get("logy"):
        ax.set_yscale("log")
    ax.set_xlabel(p.get("xlabel", ""))
    ax.set_ylabel(p.get("ylabel", ""))
    ax.set_title(f"{result.capability}: {result.kind} ({result.status})", fontsize=9, loc="left")
    tidy(ax, legend=not p.get("residual"))
    return ax
