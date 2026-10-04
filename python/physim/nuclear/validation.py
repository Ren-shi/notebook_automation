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


@check("Monte Carlo spectra", "tool", "TRIM transmitted-energy spectra (requested in pending/stopping_srim.md, part "
       "2)", "mean 1%, width 20%")
def spectra_vs_trim(ref):
    path = _file(ref, "straggling_trim.csv")
    if path is None:
        return None, math.nan, 0, "straggling_trim.csv not yet produced (TRIM runs requested)", {}
    from .stopping import Stopping

    rows = read_reference(path)[1]
    ok, worst = True, 0.0
    for r in rows:
        s = Stopping(r["ion"], r["material"], r.get("density") or None)
        mean = s.energy_after(float(r["energy_mev"]), r["thickness"])
        sigma = s.straggling(float(r["energy_mev"]), r["thickness"])
        d_mean = abs(mean / float(r["mean_mev"]) - 1)
        d_sigma = abs(sigma / float(r["sigma_mev"]) - 1)
        ok &= d_mean <= 0.01 and d_sigma <= 0.20
        worst = max(worst, d_mean / 0.01, d_sigma / 0.20)
    return bool(ok), worst, len(rows), "", {}


# ---------------------------------------------------------------------------------------------------------------
# Running and reporting


def run(capability: Optional[str] = None, kind: Optional[str] = None, reference: Optional[Path] = None) -> list:
    """Run the checks (all, or those of one capability and/or kind) and return their :class:`Result`."""
    ref = reference_dir() if reference is None else Path(reference)
    out = []
    for c in CHECKS:
        if (capability is not None and c.capability != capability) or (kind is not None and c.kind != kind):
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
