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

import csv
import math
import os
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


def _stopping_vs_srim(ref, targets: tuple) -> tuple:
    """Total stopping power (electronic + nuclear) against SR Module, judged by the LISE++ bands of
    STOPPING_TOLERANCE. SRIM's ranges are projected ranges, which are shorter than physim's path-length (CSDA)
    ranges where ions scatter (light ions in heavy targets, slow ions), so they are reported, not judged."""
    from .stopping import Stopping

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
        if abs(d) / STOPPING_TOLERANCE[key] > per_target.get(t, (0.0,))[0]:
            per_target[t] = (abs(d) / STOPPING_TOLERANCE[key], f"{r['ion']} {e_u:g} MeV/u {d:+.1%}")
        if e_u >= 10:
            dr = float(s.range(e)) / rho * 10 / float(r["projected_range_um"]) - 1
            range_10[t] = max(range_10.get(t, 0.0), abs(dr))
        xs.append(e_u)
        ratio.append(1 + d)
    failures = {k: round(v, 3) for k, v in worst.items() if v > STOPPING_TOLERANCE[k]}
    worst_frac = max(v / STOPPING_TOLERANCE[k] for k, v in worst.items())
    detail = ("all bands inside their tolerance" if not failures else f"outside: {failures}") + (
        f" (worst {worst_frac:.0%} of a band's tolerance); worst per target: "
        + ", ".join(f"{t} {w[1]}" for t, w in per_target.items())
        + "; path-length vs projected range at 10 MeV/u (not judged): "
        + ", ".join(f"{t} {v:.1%}" for t, v in range_10.items()))
    plot = {"x": xs, "physim": ratio, "reference": np.ones(len(ratio)), "xlabel": "energy (MeV/u)",
            "ylabel": "stopping power, physim / SRIM", "logx": True}
    return not failures, worst_frac, len(rows), detail, plot


@check("Stopping power and range", "tool", "SRIM 2013 SR Module: 5 ions (⁴He to ¹³²Xe) × Au, Si × 0.1–10 MeV/u "
       "(stopping_srim.csv)", "the LISE++ bands for stopping power (STOPPING_TOLERANCE); ranges reported, not judged")
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
       "weak coupling), ¹⁶O + ⁵⁸Ni at 30 MeV, five CM angles (coulex_gosia.csv)", "2% in P(θ)")
def coulex_vs_gosia(ref):
    path = _file(ref, "coulex_gosia.csv")
    if path is None:
        return None, math.nan, 0, "coulex_gosia.csv not found", {}
    from .coulex import Coulex

    rows = read_reference(path)[1]
    ours, theirs, xs = [], [], []
    for r in rows:
        c = Coulex(r["beam"], r["target"], float(r["beam_energy_mev"]), energy=float(r["state_mev"]),
                   b_up=f"{r['b_e2_e2b2']} e2b2")
        ours.append(c.probability(float(r["theta_cm_deg"]), exact=True))
        theirs.append(float(r["probability"]))
        xs.append(float(r["theta_cm_deg"]))
    worst = float(np.max(_rel(ours, theirs)))
    plot = {"x": xs, "physim": ours, "reference": theirs, "xlabel": "CM angle (deg)",
            "ylabel": "excitation probability", "logy": True}
    return worst <= 0.02, worst, len(rows), f"worst {worst:.2%}", plot


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
