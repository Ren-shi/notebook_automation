"""How a γ-ray detector answers a γ ray: efficiency against energy, attenuation on the way, the shape of the
spectrum, and runs with a calibration source.

No photon is tracked inside the crystal. The response is built from a few numbers per crystal type
(:data:`CRYSTALS`), each either from a data sheet or marked typical::

    from physim.nuclear import Experiment, response

    exp = Experiment.example("coulex_ni58")
    r = response.Response(exp, exp.gamma_detectors[0])
    r.peak_efficiency(1.332)                 # full-energy-peak efficiency at 1332 keV, as a fraction of all γ rays
    r.fwhm(1.332)                            # resolution, MeV
    run = response.source_run(exp, "152Eu", activity="37 kBq", time="1 h")
    run.efficiency_points("Ge90")            # what the simulated run says the efficiency is, line by line

**The model.** For a γ ray of energy E emitted at the target:

- *Solid angle:* the fraction of all directions the crystal faces cover.
- *Attenuation:* transmission exp(−Σ μᵢ xᵢ) through the chamber wall, the absorbers of the setup and the housing
  window, with NIST mass attenuation coefficients (:func:`mass_attenuation`).
- *Interaction:* a γ ray entering a crystal of length L interacts with probability k (1 − exp(−μ L)). The factor k
  allows for rays near the edge, the central hole of a coaxial crystal and its dead layers.
- *Peak-to-total:* of the γ rays that interact, the fraction P/T(E) = min(cap, P/T(1332 keV) (E / 1332 keV)^−q)
  leaves its full energy in the crystal.
- *The rest:* a Compton continuum from the Klein–Nishina formula for one scattering, a smaller flat part between
  the Compton edge and the peak for several scatterings, and above 1.022 MeV the single- and double-escape peaks.
- *Resolution:* FWHM(E)² = noise² + (FWHM(1332 keV)² − noise²) E / 1332 keV.

For germanium, k is fixed so that a 50 mm × 70 mm crystal has the relative efficiency of the clover's
specification sheet (21.5% of a 3 inch × 3 inch NaI crystal at 1332 keV and 25 cm, that is, 1.2 × 10⁻³).

These are typical responses, not the calibration of any real detector. A measured curve in the setup
(``efficiency_curve``) replaces the efficiency model for that detector.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np

from . import data
from .quantity import Quantity

DATA = Path(__file__).resolve().parent / "data"
FWHM_PER_SIGMA = 2.3548200450309493
ELECTRON_MEV = 0.51099895
#: Reference energy of resolutions and peak-to-total ratios, MeV (the upper ⁶⁰Co line).
REFERENCE_MEV = 1.332492
#: Absolute full-energy-peak efficiency of a 3 inch × 3 inch NaI(Tl) crystal at 1332 keV and 25 cm: the standard
#: against which the relative efficiency of a germanium detector is quoted (IEEE Std 325).
NAI_STANDARD = 1.2e-3


@dataclass(frozen=True)
class Crystal:
    """The numbers that describe one crystal material's response. ``typical`` lists those that are typical values
    rather than taken from a document."""

    material: str
    density_g_cm3: float
    #: Peak-to-total ratio at 1332 keV, its power-law slope q, and its largest value (at low energy).
    peak_to_total: float
    slope: float
    cap: float
    #: The factor k on the interaction probability.
    collection: float
    #: Resolution (FWHM) at 1332 keV, MeV, used where the setup gives none; and the electronic noise, MeV.
    fwhm: float
    noise: float
    #: The window in front of the crystal: (material, thickness in mm).
    window: tuple
    #: Of the γ rays that interact, the share that ends in the double-escape peak 1 MeV above the pair threshold
    #: (it grows as (E − 1.022 MeV)^1.5 up to ``escape_cap``); the single-escape peak has ``single_escape`` times
    #: that. And the share of the continuum that is flat between the Compton edge and the peak.
    double_escape: float
    escape_cap: float
    single_escape: float
    multiple: float
    sources: tuple
    typical: tuple


#: Crystal materials the response is defined for.
CRYSTALS = {
    "Ge": Crystal(
        "Ge", 5.323, peak_to_total=0.18, slope=0.68, cap=0.95, collection=0.0, fwhm=2.1e-3, noise=0.9e-3,
        window=("Al", 1.5), double_escape=0.03, escape_cap=0.12, single_escape=0.6, multiple=0.15,
        sources=("Mirion Technologies (Canberra), 'Clover Detectors', specification sheet C39840 (2017): relative "
                 "efficiency 21 to 22% per crystal, resolution 2.1 keV at 1332 keV and 1.05 keV at 122 keV",
                 "IEEE Std 325: relative efficiency is against 1.2e-3 for a 3 x 3 inch NaI(Tl) at 25 cm",
                 "G. F. Knoll, Radiation Detection and Measurement (4th ed.), ch. 12: efficiency and peak-to-total "
                 "of coaxial germanium detectors"),
        typical=("peak_to_total", "slope", "cap", "window", "double_escape", "single_escape", "multiple")),
    "LaBr3": Crystal(
        "LaBr3", 5.08, peak_to_total=0.20, slope=0.75, cap=0.95, collection=1.0, fwhm=28e-3, noise=0.0,
        window=("Al", 0.5), double_escape=0.03, escape_cap=0.12, single_escape=0.6, multiple=0.15,
        sources=("Saint-Gobain Crystals, 'Lanthanum Bromide Scintillators Performance Summary' (2021): density "
                 "5.08 g/cm3, resolution 2.9% at 662 keV, 2.1% at 1332 keV, 1.6% at 2615 keV",),
        typical=("peak_to_total", "slope", "cap", "collection", "window", "double_escape", "single_escape",
                 "multiple")),
}

#: The germanium crystal whose relative efficiency fixes the factor k: diameter and length in mm, and the relative
#: efficiency (the middle of the specification's 21 to 22%).
GE_REFERENCE = (50.0, 70.0, 0.215)


def _q(x) -> Optional[Quantity]:
    if x is None:
        return None
    return x if isinstance(x, Quantity) else Quantity.parse(x)


# -- attenuation ----------------------------------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _attenuation_tables() -> dict:
    out: dict = {}
    with open(DATA / "nist_attenuation.csv", encoding="utf-8") as f:
        for row in csv.DictReader(line for line in f if not line.startswith("#")):
            out.setdefault(int(row["Z"]), []).append((float(row["energy_mev"]), float(row["mu_rho_cm2_g"])))
    tables = {}
    for z, rows in out.items():
        e = np.array([r[0] for r in rows])
        # At an absorption edge the energy is listed twice; nudge the lower row down so interpolation is defined.
        for i in range(1, len(e)):
            if e[i] <= e[i - 1]:
                e[i - 1] = e[i] * (1 - 1e-9)
        tables[z] = (np.log(e), np.log(np.array([r[1] for r in rows])))
    return tables


def attenuation_elements() -> list:
    """Atomic numbers of the elements with attenuation data."""
    return sorted(_attenuation_tables())


def mass_attenuation(material: str, energy_mev) -> np.ndarray:
    """Mass attenuation coefficient μ/ρ of a material, cm²/g, at γ-ray energies in MeV (5 keV to 20 MeV).

    Interpolated on log–log axes in the NIST tables of the elements; a compound is the sum over its elements
    weighted by mass. Raises ``ValueError`` for a material with an element that has no table."""
    m = data.material(material)
    tables = _attenuation_tables()
    e = np.asarray(energy_mev, dtype=float)
    if np.any(e < 5e-3) or np.any(e > 20.0):
        raise ValueError("attenuation coefficients are tabulated from 5 keV to 20 MeV")
    masses = {}
    for z, a, n in m.atoms:
        masses[z] = masses.get(z, 0.0) + n * data.nuclide((z, a)).atomic_mass_u
    total = sum(masses.values())
    out = np.zeros_like(e)
    for z, mass in masses.items():
        if z not in tables:
            have = ", ".join(data.element(k).symbol for k in sorted(tables))
            raise ValueError(f"no attenuation data for {data.element(z).symbol} (in '{material}'); the elements "
                             f"with data are {have}")
        le, lmu = tables[z]
        out = out + mass / total * np.exp(np.interp(np.log(e), le, lmu))
    return out


def areal_density_g_cm2(material: str, thickness) -> float:
    """A thickness (a length, or an areal density) of a material as g/cm²."""
    return data.material(material).areal_density_mg_cm2(_q(thickness)) * 1e-3


def transmission(material: str, thickness, energy_mev) -> np.ndarray:
    """The fraction of γ rays that pass through ``thickness`` of a material without interacting,
    exp(−(μ/ρ) ρ x), at normal incidence."""
    return np.exp(-mass_attenuation(material, energy_mev) * areal_density_g_cm2(material, thickness))


# -- one detector's response ----------------------------------------------------------------------------------------


def _interaction(crystal: Crystal, length_mm: float, energy_mev) -> np.ndarray:
    mu = mass_attenuation(crystal.material, energy_mev) * crystal.density_g_cm3
    return 1.0 - np.exp(-mu * length_mm / 10.0)


def _peak_to_total(crystal: Crystal, energy_mev) -> np.ndarray:
    e = np.asarray(energy_mev, dtype=float)
    return np.minimum(crystal.cap, crystal.peak_to_total * (e / REFERENCE_MEV) ** (-crystal.slope))


@lru_cache(maxsize=None)
def _collection(material: str) -> float:
    c = CRYSTALS[material]
    if c.collection > 0:
        return c.collection
    # Germanium: the factor that gives the reference crystal its data-sheet relative efficiency.
    diameter, length, relative = GE_REFERENCE
    coverage = (1.0 - math.cos(math.atan2(diameter / 2, 250.0))) / 2.0
    bare = float(_interaction(c, length, REFERENCE_MEV) * _peak_to_total(c, REFERENCE_MEV))
    return relative * NAI_STANDARD / (coverage * bare)


@dataclass
class Absorber:
    """Material between the target and the crystals. ``origin`` says where it comes from: "chamber wall",
    "absorber" (the setup's list) or "window"."""

    material: str
    g_cm2: float
    origin: str


class Response:
    """The response of one γ-ray detector of an experiment (all its crystals together, or one of them)."""

    def __init__(self, experiment, detector):
        #: The experiment the detector stands in (``None``: no chamber wall).
        self.experiment = experiment
        self.detector = detector
        self.material = detector.material or "Ge"
        if self.material not in CRYSTALS:
            raise ValueError(f"no response for a crystal of '{self.material}'; the crystal materials are "
                             f"{', '.join(CRYSTALS)}")
        self.crystal = CRYSTALS[self.material]
        self.elements = detector.elements()
        radius = detector.crystal_radius_mm()
        #: Crystal length, mm; a detector given only a radius is taken as long as it is wide.
        self.length_mm = (_q(detector.crystal_length).to("mm") if detector.crystal_length is not None
                          else 2.0 * radius)
        self.assumed_length = detector.crystal_length is None
        self.coverage = [(1.0 - math.cos(math.atan2(r, math.hypot(*c)))) / 2.0 for _, c, r in self.elements]
        self.absorbers = self._absorbers()
        self._curve = self._measured_curve()

    def _absorbers(self) -> list:
        out = []
        ch = self.experiment.chamber if self.experiment is not None else None
        if ch is not None and ch.wall_thickness is not None and ch.wall_material is not None:
            if _q(self.detector.distance).to("mm") >= _q(ch.radius).to("mm"):
                out.append(Absorber(ch.wall_material, areal_density_g_cm2(ch.wall_material, ch.wall_thickness),
                                    "chamber wall"))
        for material, thickness in self.detector.absorbers or ():
            out.append(Absorber(material, areal_density_g_cm2(material, thickness), "absorber"))
        material, mm = self.crystal.window
        out.append(Absorber(material, areal_density_g_cm2(material, Quantity(mm, "mm")), "window"))
        return out

    def _measured_curve(self) -> Optional[tuple]:
        curve = self.detector.efficiency_curve
        if not curve:
            return None
        pts = sorted((_q(e).to("MeV"), _q(eff).to("%") / 100.0) for e, eff in curve)
        return np.log([p[0] for p in pts]), np.log([p[1] for p in pts])

    @property
    def source(self) -> str:
        """Where the efficiency comes from: "curve" (measured points of the setup), "fixed" (the setup's single
        ``efficiency``) or "model"."""
        if self._curve is not None:
            return "curve"
        return "fixed" if self.detector.efficiency is not None else "model"

    def transmission(self, energy_mev) -> np.ndarray:
        """Fraction of γ rays that reach the crystals through the wall, the absorbers and the window."""
        e = np.asarray(energy_mev, dtype=float)
        t = np.ones_like(e)
        for a in self.absorbers:
            t = t * np.exp(-mass_attenuation(a.material, e) * a.g_cm2)
        return t

    def peak_to_total(self, energy_mev) -> np.ndarray:
        """Of the γ rays that interact in a crystal, the fraction that leaves its full energy there."""
        return _peak_to_total(self.crystal, energy_mev)

    def _model_peak(self, energy_mev, element: Optional[int] = None) -> np.ndarray:
        e = np.asarray(energy_mev, dtype=float)
        cover = sum(self.coverage) if element is None else self.coverage[element]
        intrinsic = _collection(self.material) * _interaction(self.crystal, self.length_mm, e) * self.peak_to_total(e)
        return cover * self.transmission(e) * np.minimum(intrinsic, 1.0)

    def peak_efficiency(self, energy_mev, element: Optional[int] = None) -> np.ndarray:
        """Full-energy-peak efficiency at γ-ray energies in MeV: the fraction of all γ rays emitted at the target
        that leave their full energy in the detector (``element`` picks one crystal of a clover).

        A measured ``efficiency_curve`` is interpolated on log–log axes, and beyond its ends follows the model's
        shape. A single ``efficiency`` in the setup is used at every energy."""
        e = np.asarray(energy_mev, dtype=float)
        share = 1.0 if element is None else self.coverage[element] / sum(self.coverage)
        if self._curve is not None:
            le, leff = self._curve
            inside = np.exp(np.interp(np.log(e), le, leff))
            lo, hi = math.exp(le[0]), math.exp(le[-1])
            model = self._model_peak(e)
            below = math.exp(leff[0]) * model / self._model_peak(lo)
            above = math.exp(leff[-1]) * model / self._model_peak(hi)
            return share * np.where(e < lo, below, np.where(e > hi, above, inside))
        if self.detector.efficiency is not None:
            return share * np.full_like(e, _q(self.detector.efficiency).to("%") / 100.0)
        return self._model_peak(e, element)

    def total_efficiency(self, energy_mev, element: Optional[int] = None) -> np.ndarray:
        """Fraction of all γ rays emitted that leave any energy in the detector: peak efficiency / (P/T)."""
        return self.peak_efficiency(energy_mev, element) / self.peak_to_total(energy_mev)

    def fwhm(self, energy_mev) -> np.ndarray:
        """Energy resolution (FWHM), MeV, at γ-ray energies in MeV. The setup's ``resolution`` is the value at
        1332 keV."""
        e = np.asarray(energy_mev, dtype=float)
        ref = _q(self.detector.resolution).to("MeV") if self.detector.resolution is not None else self.crystal.fwhm
        noise = min(self.crystal.noise, ref)
        return np.sqrt(noise**2 + (ref**2 - noise**2) * e / REFERENCE_MEV)

    def escape_fractions(self, energy_mev: float) -> tuple:
        """(single escape, double escape): the shares of the interacting γ rays that end in each escape peak."""
        x = energy_mev - 2 * ELECTRON_MEV
        if x <= 0:
            return 0.0, 0.0
        de = min(self.crystal.escape_cap, self.crystal.double_escape * x**1.5)
        return self.crystal.single_escape * de, de

    def shape(self, energy_mev: float, edges: np.ndarray) -> np.ndarray:
        """The spectrum of one interacting γ ray of ``energy_mev``, as probabilities per bin of ``edges`` (MeV):
        the full-energy peak, the escape peaks and the Compton continuum, each with the detector's resolution.
        Sums to 1, apart from what falls outside the bins."""
        e0 = float(energy_mev)
        pt = float(self.peak_to_total(e0))
        se, de = self.escape_fractions(e0)
        se, de = min(se, (1 - pt) / 2), min(de, (1 - pt) / 2)
        out = pt * _gauss_bins(e0, float(self.fwhm(e0)), edges)
        for share, energy in ((se, e0 - ELECTRON_MEV), (de, e0 - 2 * ELECTRON_MEV)):
            if share > 0:
                out = out + share * _gauss_bins(energy, float(self.fwhm(energy)), edges)
        rest = 1 - pt - se - de
        if rest > 0:
            out = out + rest * _smeared(compton_continuum(e0, edges, self.crystal.multiple), edges,
                                        float(self.fwhm(max(compton_edge(e0), 0.02))))
        return out

    def describe(self) -> dict:
        """The response in words and numbers, for the app and the report."""
        return {"material": self.material, "source": self.source, "crystals": len(self.elements),
                "length_mm": self.length_mm, "assumed_length": self.assumed_length,
                "coverage": sum(self.coverage),
                "absorbers": [{"material": a.material, "g_cm2": a.g_cm2, "origin": a.origin} for a in self.absorbers],
                "typical": list(self.crystal.typical), "sources": list(self.crystal.sources)}


def compton_edge(energy_mev: float) -> float:
    """Largest energy a γ ray gives an electron in one Compton scattering, MeV."""
    return energy_mev * 2 * energy_mev / (ELECTRON_MEV + 2 * energy_mev)


def compton_continuum(energy_mev: float, edges: np.ndarray, multiple: float = 0.0) -> np.ndarray:
    """Probability per bin of the energy left by a γ ray that does not end in a peak.

    A share 1 − ``multiple`` follows the Klein–Nishina spectrum of the electron from one scattering, from 0 to
    the Compton edge; the share ``multiple`` is flat between the edge and the full energy (several scatterings)."""
    k = energy_mev / ELECTRON_MEV
    edge = compton_edge(energy_mev)
    fine = np.linspace(0.0, edge, 2001)
    t = (fine[:-1] + fine[1:]) / 2
    s = t / energy_mev
    # Klein–Nishina, as dσ/dT against the electron's energy T.
    kn = 2 + s**2 / (k**2 * (1 - s) ** 2) + s / (1 - s) * (s - 2 / k)
    cdf = np.concatenate([[0.0], np.cumsum(kn)])
    cdf /= cdf[-1]
    single = np.diff(np.interp(np.clip(edges, 0.0, edge), fine, cdf))
    flat = np.diff(np.clip((edges - edge) / (energy_mev - edge), 0.0, 1.0))
    return (1 - multiple) * single + multiple * flat


def _gauss_bins(mean: float, fwhm: float, edges: np.ndarray) -> np.ndarray:
    from math import erf

    sigma = max(fwhm / FWHM_PER_SIGMA, 1e-9)
    cdf = np.array([0.5 * (1 + erf((x - mean) / (sigma * math.sqrt(2)))) for x in edges])
    return np.diff(cdf)


def _smeared(hist: np.ndarray, edges: np.ndarray, fwhm: float) -> np.ndarray:
    """A histogram folded with a Gaussian of constant width (uniform bins)."""
    width = edges[1] - edges[0]
    sigma = fwhm / FWHM_PER_SIGMA / width
    if sigma < 0.3:
        return hist
    n = int(math.ceil(4 * sigma))
    x = np.arange(-n, n + 1)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    return np.convolve(hist, kernel / kernel.sum(), mode="same")


# -- calibration sources --------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Line:
    energy_mev: float
    #: Photons per decay, and its uncertainty.
    intensity: float
    intensity_unc: float
    kind: str


@dataclass(frozen=True)
class Source:
    """A calibration source: its photon lines, from the DDEP evaluation named in ``evaluation``."""

    nuclide: str
    half_life_s: float
    lines: tuple
    evaluation: str

    def strong(self, least: float = 0.02) -> list:
        """The γ lines of at least ``least`` photons per decay."""
        return [ln for ln in self.lines if ln.kind == "gamma" and ln.intensity >= least]


@lru_cache(maxsize=None)
def _sources() -> dict:
    rows: dict = {}
    with open(DATA / "calibration_sources.csv", encoding="utf-8") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            rows.setdefault(r["nuclide"], []).append(r)
    return {n: Source(n, float(rs[0]["half_life_s"]),
                      tuple(Line(float(r["energy_kev"]) * 1e-3, float(r["intensity_percent"]) / 100,
                                 float(r["intensity_unc_percent"]) / 100, r["kind"]) for r in rs),
                      rs[0]["evaluation"]) for n, rs in rows.items()}


def source_names() -> list:
    """The calibration sources with decay data."""
    return sorted(_sources(), key=lambda n: (int("".join(c for c in n if c.isdigit())), n))


def source(nuclide: str) -> Source:
    """A calibration source by nuclide ("152Eu", "60Co", ...)."""
    s = _sources().get(nuclide)
    if s is None:
        raise ValueError(f"no decay data for a '{nuclide}' source; the sources are {', '.join(source_names())}")
    return s


@dataclass
class SourceRun:
    """A simulated run with a calibration source at the target position: one spectrum per crystal."""

    source: Source
    #: Decays during the run.
    decays: float
    time_s: float
    edges: np.ndarray
    #: {crystal name: counts per bin}, and the same before the counting statistics were drawn.
    spectra: dict
    expected: dict
    #: {crystal name: (Response, element index or None)}.
    responses: dict = field(default_factory=dict)

    def names(self) -> list:
        return list(self.spectra)

    def detector(self, name: str) -> np.ndarray:
        """Counts per bin of a detector: one crystal by its name ("Clo A"), or all crystals of a detector summed
        by the detector's name."""
        if name in self.spectra:
            return self.spectra[name]
        parts = [h for n, h in self.spectra.items() if n.rsplit(" ", 1)[0] == name]
        if not parts:
            raise KeyError(f"no γ-ray detector named {name!r}; the names are {', '.join(self.spectra)}")
        return np.sum(parts, axis=0)

    def _response(self, name: str) -> tuple:
        if name in self.responses:
            return self.responses[name]
        r = next(r for n, (r, _) in self.responses.items() if n.rsplit(" ", 1)[0] == name)
        return r, None

    def peak_area(self, name: str, energy_mev: float) -> tuple:
        """(area, uncertainty) of the peak at ``energy_mev``: the counts within ±3σ, less the mean level of the
        bands from 3σ to 6σ on each side (leaving out bins near the source's other lines), corrected for the part
        of the peak outside the window."""
        h = self.detector(name)
        r, _ = self._response(name)
        sigma = float(r.fwhm(energy_mev)) / FWHM_PER_SIGMA
        centres = (self.edges[:-1] + self.edges[1:]) / 2
        d = np.abs(centres - energy_mev)
        peak = d <= 3 * sigma
        side = (d > 3 * sigma) & (d <= 6 * sigma)
        # Bands are taken clear of the source's other lines, which are known.
        for ln in self.source.lines:
            if abs(ln.energy_mev - energy_mev) > 1e-9:
                side &= np.abs(centres - ln.energy_mev) > 3.5 * float(r.fwhm(ln.energy_mev)) / FWHM_PER_SIGMA
        if peak.sum() < 3 or side.sum() < 2:
            raise ValueError("the bins are too wide for this peak; use more bins")
        background = h[side].mean() * peak.sum()
        # The bins taken as the peak cover a little less or more than ±3σ: correct for the part of the Gaussian
        # they hold.
        inside = np.flatnonzero(peak)
        held = float(_gauss_bins(energy_mev, sigma * FWHM_PER_SIGMA,
                                 np.array([self.edges[inside[0]], self.edges[inside[-1] + 1]]))[0])
        # The bands hold a little of the peak's tails, which the background subtraction takes away again.
        tails = _gauss_bins(energy_mev, sigma * FWHM_PER_SIGMA, self.edges)[side].sum()
        held -= float(tails) * peak.sum() / side.sum()
        area = (h[peak].sum() - background) / held
        unc = math.sqrt(h[peak].sum() + h[side].sum() * (peak.sum() / side.sum()) ** 2) / held
        return float(area), float(unc)

    def blended(self, name: str, energy_mev: float) -> bool:
        """Whether another line of the source, at least a hundredth as strong, lies within 6σ of this one (the
        two peaks then overlap, and the simple peak area is not this line's alone)."""
        r, _ = self._response(name)
        sigma = float(r.fwhm(energy_mev)) / FWHM_PER_SIGMA
        this = next(ln for ln in self.source.lines if abs(ln.energy_mev - energy_mev) < 1e-9)
        return any(ln is not this and abs(ln.energy_mev - energy_mev) < 6 * sigma
                   and ln.intensity >= 0.01 * this.intensity for ln in self.source.lines)

    def efficiency_points(self, name: str, least: float = 0.02) -> list:
        """The full-energy-peak efficiency the run gives at each strong line that stands alone: a list of
        {"energy_mev", "efficiency", "uncertainty", "true"}, "true" being the efficiency that was put in."""
        r, element = self._response(name)
        out = []
        for ln in self.source.strong(least):
            if self.blended(name, ln.energy_mev):
                continue
            area, unc = self.peak_area(name, ln.energy_mev)
            emitted = self.decays * ln.intensity
            out.append({"energy_mev": ln.energy_mev, "efficiency": area / emitted,
                        "uncertainty": unc / emitted,
                        "true": float(r.peak_efficiency(ln.energy_mev, element))})
        return out

    def peak_to_total(self, name: str, least: float = 0.02) -> dict:
        """The peak-to-total ratio of the run: the summed areas of the strong γ peaks over all counts in the
        spectrum. Returns {"value", "uncertainty", "true"}; "true" is what the response put in, Σ I ε_peak over
        Σ I ε_total for the same lines and all lines.

        This is how the ratio is measured with a source of one or two lines (¹³⁷Cs, ⁶⁰Co)."""
        r, element = self._response(name)
        areas = [self.peak_area(name, ln.energy_mev) for ln in self.source.strong(least)]
        total = float(self.detector(name).sum())
        peak = sum(a for a, _ in areas)
        unc = math.sqrt(sum(u * u for _, u in areas)) / total
        lines = [ln for ln in self.source.lines if ln.energy_mev >= 5e-3]
        true = (sum(ln.intensity * float(r.peak_efficiency(ln.energy_mev, element))
                    for ln in self.source.strong(least))
                / sum(ln.intensity * float(r.total_efficiency(ln.energy_mev, element)) for ln in lines))
        return {"value": peak / total, "uncertainty": unc, "true": true}


def source_run(experiment, nuclide: str = "152Eu", activity="37 kBq", time="1 h", seed: int = 1,
               bin_kev: Optional[float] = None, max_mev: Optional[float] = None) -> SourceRun:
    """Simulate a run with a calibration source at the target position, without beam.

    Each line of the source gives, in each crystal, activity × time × intensity × total efficiency interacting
    γ rays on average, spread over the spectrum by :meth:`Response.shape`; the counts in each bin are then drawn
    from a Poisson distribution, which is exact for independent γ rays. The activity is that at the start of the
    run and decays with the source's half-life.

    Not included: two γ rays of one decay summing in a crystal, the room background, dead time. ``bin_kev``
    defaults to a third of the narrowest resolution."""
    src = source(nuclide)
    a0 = _q(activity).to("Bq")
    t = _q(time).to("s")
    lam = math.log(2) / src.half_life_s
    decays = a0 * (1 - math.exp(-lam * t)) / lam
    responses = {}
    for i, gd in enumerate(experiment.gamma_detectors):
        r = Response(experiment, gd)
        for k, (label, _, _) in enumerate(r.elements):
            name = (gd.name or f"G{i + 1}") + (f" {label}" if label else "")
            responses[name] = (r, k if len(r.elements) > 1 else None)
    if not responses:
        raise ValueError("the setup has no γ-ray detectors")
    top = max_mev or 1.08 * max(ln.energy_mev for ln in src.lines) + 0.05
    if bin_kev is None:
        bin_kev = max(0.25, min(1e3 * float(r.fwhm(0.1)) for r, _ in responses.values()) / 3)
    edges = np.arange(0.0, top + bin_kev * 1e-3, bin_kev * 1e-3)
    rng = np.random.default_rng(seed)
    spectra, expected = {}, {}
    for name, (r, k) in responses.items():
        mu = np.zeros(len(edges) - 1)
        for ln in src.lines:
            if ln.energy_mev < 5e-3:
                continue
            interacting = decays * ln.intensity * float(r.total_efficiency(ln.energy_mev, k))
            mu += interacting * r.shape(ln.energy_mev, edges)
        expected[name] = mu
        spectra[name] = rng.poisson(mu).astype(float)
    return SourceRun(src, decays, t, edges, spectra, expected, responses)


__all__ = ["CRYSTALS", "Absorber", "Crystal", "Line", "NAI_STANDARD", "REFERENCE_MEV", "Response", "Source",
           "SourceRun", "areal_density_g_cm2", "attenuation_elements", "compton_continuum", "compton_edge",
           "mass_attenuation", "source", "source_names", "source_run", "transmission"]
