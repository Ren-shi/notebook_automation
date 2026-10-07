"""The run record: every number of the plan and of the analysis with its explanation, and the record as one page
for the app and for printing.

Each :class:`Explanation` gives a quantity in four parts: the formula, the same formula with this run's numbers
substituted, what the quantity means physically, and the assumptions and a reference. The numbers come from the
calculation that produced the result (they are passed in, never computed again here), and each explanation
carries a rule that recomputes its value from the numbers shown, so a test can check that the two agree::

    from physim.nuclear.planner import Planner
    from physim.nuclear.record import explanations, record_html

    p = Planner.example("coulex_ni58")
    for x in explanations(p):                 # the plan: rates, excitation, Doppler shift, beam time, ...
        print(x.key, x.value, x.unit)
        assert abs(x.recompute() - x.value) <= 1e-6 * abs(x.value)
    html = record_html(p)                     # the record, with the analysis chain if one has been run

Written for a reader who knows nuclear physics and is new to Coulomb excitation.
"""

from __future__ import annotations

import html
import math
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from . import data
from .coulex import HBARC_MEV_FM
from .levels import quadrupole_factor, rate_per_b, weisskopf_unit
from .rutherford import E2_MEV_FM

#: Where the theory pages of the documentation are, for the links in each explanation.
DOCS = "https://ren-shi.github.io/notebook_automation/"


@dataclass
class Explanation:
    """One calculated quantity and its explanation."""

    key: str
    title: str
    value: float
    unit: str
    formula: str
    #: The formula with this run's numbers in it.
    substituted: str
    meaning: str
    assumptions: str
    reference: str
    #: The numbers shown in ``substituted``, and the rule that gives ``value`` from them.
    inputs: dict = field(default_factory=dict)
    rule: Optional[Callable] = field(default=None, repr=False)

    def recompute(self) -> float:
        """The value computed again from :attr:`inputs` by :attr:`rule` (``value`` if there is no rule)."""
        return float(self.rule(**self.inputs)) if self.rule is not None else self.value

    def as_dict(self) -> dict:
        return {"key": self.key, "title": self.title, "value": self.value, "unit": self.unit,
                "formula": self.formula, "substituted": self.substituted, "meaning": self.meaning,
                "assumptions": self.assumptions, "reference": self.reference, "inputs": dict(self.inputs)}


def _g(x: float, digits: int = 4) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "—"
    return f"{x:.{digits}g}"


# -- the plan -------------------------------------------------------------------------------------------------------


def explanations(planner, detector: Optional[str] = None, gamma_detector: Optional[str] = None) -> list:
    """The explanations of the plan's numbers for one particle detector (the first if none is named), one γ-ray
    detector (the first) and the whole: solid angle, elastic rate, counts in the run, beam time, the excitation
    probability, the γ-ray efficiency, the coincidence rate and the Doppler shift."""
    exp = planner.experiment
    r = planner._rates()
    out = []
    g = r.array[detector] if detector else r.array.geometries[0]
    name = g.name
    omega = g.solid_angle()
    theta_lo, theta_hi = g.theta_range()
    out.append(Explanation(
        "solid_angle", f"Solid angle of {name}", omega, "msr",
        "Ω = ∫ (r̂ · n̂) / r² dA over the face",
        f"Ω = {_g(omega)} msr, the face seen from the target between θ = {theta_lo:.1f}° and {theta_hi:.1f}°",
        "The fraction of all directions the detector covers: 4π sr is everything, so this detector catches "
        f"{omega / (4e3 * math.pi):.2%} of particles sent out evenly.",
        "The face is flat and the target a point (or a spot of the given size); no gaps between strips.",
        f"{DOCS}nuclear-setup.html#detector-models-and-the-chamber",
        {"omega_msr": omega}, lambda omega_msr: omega_msr))

    elastic_channels = [ch for ch in r.channels if ch.excitation is None]
    ch = elastic_channels[0] if elastic_channels else r.channels[0]
    rate_all = r.rate(name, counted=False)
    pps = r.particles_per_second
    sigma_eff = rate_all / (pps * ch.atoms_per_cm2 * 1e-27) if pps * ch.atoms_per_cm2 > 0 else 0.0
    out.append(Explanation(
        "rate", f"Particles per second in {name}", rate_all, "1/s",
        "rate = I_beam × N_target × ∫ dσ/dΩ dΩ, with dσ/dΩ = (Z₁Z₂e²/4E)² / sin⁴(θ/2) (Rutherford, in the CM frame)",
        f"rate = {_g(pps, 3)} particles/s × {_g(ch.atoms_per_cm2, 3)} nuclei/cm² × {_g(sigma_eff, 3)} mb × "
        f"10⁻²⁷ cm²/mb = {_g(rate_all)} /s",
        "Every beam particle sees the target's nuclei side by side; the cross section is the area each nucleus "
        "presents for scattering into the detector. Rutherford's formula rises steeply towards small angles, so "
        "forward detectors count fastest.",
        "Pure Coulomb scattering below the barrier; the target's thickness is averaged over; particles that "
        "another detector or dead material stops first are left out.",
        f"{DOCS}physics-register/rutherford.html",
        {"pps": pps, "atoms_per_cm2": ch.atoms_per_cm2, "sigma_mb": sigma_eff},
        lambda pps, atoms_per_cm2, sigma_mb: pps * atoms_per_cm2 * sigma_mb * 1e-27))

    what = r.measured
    counts = r.counts_in_run(name, what=what)
    rate_counted = r.rate(name, what=what)
    kind = {"all": "particles above threshold", "excitations": "excitation events",
            "coincidences": "particle–γ coincidences"}[what]
    out.append(Explanation(
        "counts", f"Counts in {name} in the run", counts, "counts",
        "N = rate × T",
        f"N = {_g(rate_counted)} /s × {_g(r.beam_time_s)} s = {_g(counts)} {kind}",
        f"What the measurement collects in {name} in the planned beam time: {kind}, which is what the "
        "measurement uses. The relative statistical uncertainty of a count is 1/√N.",
        "A steady beam; no dead time.",
        f"{DOCS}nuclear-setup.html#rates-beam-time-and-simulated-spectra",
        {"rate": rate_counted, "time_s": r.beam_time_s}, lambda rate, time_s: rate * time_s))

    if r.counts_wanted:
        rate_measured = r.rate(name, what=what)
        t = r.beam_time_for(name, what=what)
        out.append(Explanation(
            "beam_time", f"Beam time for {r.counts_wanted} {what} in {name}", t / 3600, "h",
            "T = N_wanted / rate",
            f"T = {r.counts_wanted} / {_g(rate_measured)} /s = {_g(t)} s = {_g(t / 3600)} h",
            f"How long the beam must run for {name} to collect the counts asked for, which fixes the statistical "
            f"precision 1/√N = {1 / math.sqrt(r.counts_wanted):.1%}.",
            "The rate stays what the plan says; for Coulomb excitation the counts are particle–γ coincidences.",
            f"{DOCS}nuclear-setup.html#rates-beam-time-and-simulated-spectra",
            {"wanted": r.counts_wanted, "rate": rate_measured}, lambda wanted, rate: wanted / rate / 3600))

    exc = exp.excitation
    if exc is not None:
        from .rates import coulex_for

        chx = next(c for c in r.channels if c.excitation is not None)
        cx = coulex_for(exp, chx, r.layers)
        theta = 150.0
        p_val = float(cx.probability(theta, exact=True))
        lam = cx.lam
        pref = (4 * math.pi * cx.z_exciting * E2_MEV_FM / (HBARC_MEV_FM * (2 * lam + 1))) ** 2
        amp2 = p_val / (pref * cx.b_up / (2 * lam + 1)) if p_val > 0 else 0.0
        out.append(Explanation(
            "excitation_probability", f"Excitation probability at θ_CM = {theta:g}°", p_val, "",
            "P(θ) = (4π Z e² / (ħc (2λ+1)))² × B(Eλ↑)/(2λ+1) × Σ_μ |Y_λμ(π/2, 0) I_λμ(θ, ξ) / (β aλ)|²",
            f"P = {_g(pref, 4)} × {_g(cx.b_up, 4)} e²fm^{2 * lam}/{2 * lam + 1} × {_g(amp2, 4)} = {_g(p_val, 4)}",
            "The chance that one collision, scattering to this angle, lifts the nucleus to the state: the electric "
            "field of the passing partner acts on it for a moment, and the state's B(Eλ) says how readily it "
            f"responds. At ξ = {cx.xi:.2f} the collision is {'fast enough' if cx.xi < 1 else 'slow'} compared with "
            "the state's period; the closer the approach (larger θ), the larger P.",
            "First-order perturbation theory on the classical Rutherford orbit (Alder and Winther); one state "
            "from a 0⁺ ground state; no second step and no reorientation.",
            f"{DOCS}physics-register/coulomb-excitation.html",
            {"pref": pref, "b_up": cx.b_up, "two_lam_plus_one": 2 * lam + 1, "amp2": amp2},
            lambda pref, b_up, two_lam_plus_one, amp2: pref * b_up / two_lam_plus_one * amp2))

        if exp.gamma_detectors:
            from .response import Response

            gd = next((x for x in exp.gamma_detectors if x.name == gamma_detector), exp.gamma_detectors[0])
            resp = Response(exp, gd)
            e0 = exc.energy_mev
            cover = sum(resp.coverage)
            trans = float(resp.transmission(e0))
            intrinsic = float(resp.peak_efficiency(e0)) / (cover * trans) if cover * trans > 0 else 0.0
            eff = float(resp.peak_efficiency(e0))
            out.append(Explanation(
                "gamma_efficiency", f"Full-energy-peak efficiency of {gd.name} at {1e3 * e0:.0f} keV", eff, "",
                "ε = Ω/4π × T(E) × ε_int(E)",
                f"ε = {_g(cover, 4)} × {_g(trans, 4)} × {_g(intrinsic, 4)} = {_g(eff, 4)}",
                "Of all γ rays sent out from the target, the fraction that leaves its full energy in this detector: "
                "the solid angle it covers, times the share that passes the material in front, times the share "
                "that interacts in the crystal and ends in the peak rather than the Compton continuum.",
                "A typical response of such a crystal unless the setup gives a measured efficiency or curve; "
                "normal incidence through the absorbers.",
                f"{DOCS}nuclear-setup.html#γ-ray-response",
                {"coverage": cover, "transmission": trans, "intrinsic": intrinsic},
                lambda coverage, transmission, intrinsic: coverage * transmission * intrinsic))

            exc_rate = r.rate(name, what="excitations")
            eff_all = r.gamma_efficiency(name)[0]
            coinc = r.rate(name, what="coincidences")
            out.append(Explanation(
                "coincidences", f"Particle–γ coincidences per second with {name}", coinc, "1/s",
                "rate_pγ = rate_excitations × Σ_g ε_g F_g",
                f"rate_pγ = {_g(exc_rate)} /s × {_g(eff_all, 4)} = {_g(coinc)} /s",
                "An excitation seen in the particle detector counts for the measurement only when its γ ray is "
                "caught too. Σ ε_g F_g is the γ-ray detectors' efficiency for γ rays in coincidence with this "
                "particle detector, F_g being the angular-correlation factor (γ rays are not sent out evenly).",
                "The angular correlation of first-order excitation; every excitation decays by this γ ray.",
                f"{DOCS}nuclear-setup.html#orientation-and-the-particle–γ-correlation",
                {"excitations": exc_rate, "efficiency": eff_all}, lambda excitations, efficiency: excitations * efficiency))

            rows = [d for d in planner.gamma()["doppler"] if d["particle_detector"] == name]
            if rows:
                d = rows[0]
                beta = d["beta_mean"]
                cos_alpha = (1 - e0 * math.sqrt(1 - beta**2) / (d["mean_kev"] * 1e-3)) / beta if beta > 0 else 0.0
                out.append(Explanation(
                    "doppler", f"Doppler-shifted γ-ray energy in {d['gamma_detector']} with {name}",
                    d["mean_kev"], "keV",
                    "E_γ = E₀ √(1 − β²) / (1 − β cos α)",
                    f"E_γ = {1e3 * e0:.1f} keV × √(1 − {beta:.4f}²) / (1 − {beta:.4f} × {cos_alpha:.3f}) = "
                    f"{d['mean_kev']:.1f} keV",
                    "The nucleus emits in flight, so the γ ray is shifted up when sent forward and down when sent "
                    "back, by about β cos α. The spread of angles over the two detectors gives the peak a width of "
                    f"{d['doppler_fwhm_kev']:.1f} keV before the detector's resolution, which the Doppler correction "
                    "removes as far as the segments allow.",
                    "The nucleus decays after leaving the target, moving with its speed at the exit; the mean is "
                    "weighted by the excitation cross section.",
                    f"{DOCS}nuclear-setup.html#particle–γ-events",
                    {"e0_kev": 1e3 * e0, "beta": beta, "cos_alpha": cos_alpha},
                    lambda e0_kev, beta, cos_alpha: e0_kev * math.sqrt(1 - beta**2) / (1 - beta * cos_alpha)))
    return out


# -- the analysis chain ---------------------------------------------------------------------------------------------


def analysis_explanations(planner, result) -> list:
    """The chain from the peak area to B(E2) and the shape, from a :class:`~physim.nuclear.analysis.Result`."""
    exp = planner.experiment
    exc = exp.excitation
    out = []
    gate = next((s["text"] for s in result.steps if s["step"] == "yield"), "")
    gate_eff = 1.0
    if "gate efficiency" in gate:
        gate_eff = float(gate.split("gate efficiency")[1].split(")")[0])
    out.append(Explanation(
        "area", "Peak area", result.area, "counts",
        "A = Σ counts within ±3σ of the peak − background under them − random coincidences",
        f"A = {_g(result.fit.area)} (counts over the fitted line) "
        f"− {_g(result.fit.randoms)} (randoms) = {_g(result.area)} ± {_g(math.hypot(result.fit.area_unc, math.sqrt(max(result.fit.randoms, 0))))}",
        "The γ rays of the transition that the detectors recorded in full, after the Doppler correction put them "
        "in one peak, over the continuum of other γ rays and the chance coincidences.",
        "A Gaussian peak on a straight-line background; the random coincidences expected from the singles rates "
        "and the coincidence window.",
        f"{DOCS}nuclear-setup.html#automatic-analysis",
        {"fit_area": result.fit.area, "randoms": result.fit.randoms}, lambda fit_area, randoms: fit_area - randoms))
    yield_ = result.area / (result.efficiency * result.correlation * result.branch * gate_eff)
    out.append(Explanation(
        "yield", "Excitations in the gated rings", yield_, "",
        "Y = A / (ε × F × b_γ × g)",
        f"Y = {_g(result.area)} / ({_g(result.efficiency, 4)} × {_g(result.correlation, 4)} × {_g(result.branch)} "
        f"× {_g(gate_eff, 4)}) = {_g(yield_)}",
        "How many excitations the gated particles stand for: the peak area divided by the share of γ rays the "
        "detectors catch (efficiency ε and angular-correlation factor F), the share of decays that give this γ "
        "ray (b_γ, less than 1 with internal conversion) and the share of excitations the particle gate keeps (g).",
        "ε from the response, F from the first-order orientation, g from the simulation.",
        f"{DOCS}nuclear-setup.html#automatic-analysis",
        {"area": result.area, "efficiency": result.efficiency, "correlation": result.correlation,
         "branch": result.branch, "gate": gate_eff},
        lambda area, efficiency, correlation, branch, gate: area / (efficiency * correlation * branch * gate)))
    p_mean = yield_ / result.norm_counts if result.norm_counts else float("nan")
    norm_text = next((s["text"] for s in result.steps if s["step"] == "normalisation"), "")
    if result.normalisation == "rutherford":
        p_per_b = p_mean / result.b_e2fm4 if result.b_e2fm4 else float("nan")
        out.append(Explanation(
            "mean_probability", "Mean excitation probability in the gated rings", p_mean, "",
            "⟨P⟩ = Y / N_elastic",
            f"⟨P⟩ = {_g(yield_)} / {_g(result.norm_counts)} = {_g(p_mean, 4)}",
            "The excitations per elastically scattered particle in the same rings. Normalising to the elastic "
            "particles cancels the beam current and the target thickness, which are hard to know precisely.",
            "The same rings see both; the elastic particles are counted above threshold.",
            f"{DOCS}nuclear-setup.html#automatic-analysis",
            {"yield_": yield_, "n_elastic": result.norm_counts}, lambda yield_, n_elastic: yield_ / n_elastic))
        out.append(Explanation(
            "b_e2", "B(E2↑) from the mean probability", result.b_e2fm4, "e²fm⁴",
            "B(E2↑) = ⟨P⟩ / (⟨P⟩/B)_theory, since in first order P ∝ B(E2↑)",
            f"B(E2↑) = {_g(p_mean, 4)} / {_g(p_per_b, 4)} per e²fm⁴ = {_g(result.b_e2fm4)} e²fm⁴ "
            f"= {_g(result.b_e2b2)} e²b²",
            "The excitation probability measures the reduced transition probability directly: (⟨P⟩/B)_theory is "
            "the first-order probability per unit B(E2), averaged over the same rings with the elastic "
            "distribution. B(E2↑) is the square of the matrix element ⟨2‖M(E2)‖0⟩: how strongly the electric "
            "quadrupole operator connects the ground state to the 2⁺ state.",
            "First-order theory; the matrix element enters only through its size, so no sign is obtained.",
            f"{DOCS}physics-register/coulomb-excitation.html",
            {"p_mean": p_mean, "p_per_b": p_per_b}, lambda p_mean, p_per_b: p_mean / p_per_b))
    else:
        out.append(Explanation(
            "b_e2", "B(E2↑) from the reference transition", result.b_e2fm4, "e²fm⁴",
            "B(E2↑) = B_ref × (Y / Y_ref) × (⟨P⟩/B)_ref / (⟨P⟩/B)_theory",
            norm_text + f"; B(E2↑) = {_g(result.b_e2fm4)} e²fm⁴",
            "The yield of the unknown transition over that of a transition whose B(E2) is known, in the same "
            "collision: the beam, the target and most of the efficiency cancel in the ratio.",
            "The reference's angular correlation is taken as even; first order for both.",
            f"{DOCS}nuclear-setup.html#automatic-analysis",
            {"b": result.b_e2fm4}, lambda b: b))
    a_mass = data.nuclide(_excited(exp)).A
    wu = weisskopf_unit(2, a_mass)
    out.append(Explanation(
        "weisskopf", "B(E2↑) in Weisskopf units", result.b_wu, "W.u.",
        "B_W(E2) = (1/4π) (3/5)² (1.2 A^⅓ fm)⁴ e²;  B / B_W",
        f"B_W = {_g(wu, 4)} e²fm⁴ for A = {a_mass};  {_g(result.b_e2fm4)} / {_g(wu, 4)} = {_g(result.b_wu, 3)} W.u.",
        "The strength one proton moving alone would give. Tens of Weisskopf units mean many nucleons move "
        "together: a collective, deformed nucleus. A few W.u. or less is a single-particle transition.",
        "The Weisskopf estimate with R = 1.2 A^⅓ fm.",
        f"{DOCS}nuclear-setup.html#level-schemes",
        {"b": result.b_e2fm4, "b_w": wu}, lambda b, b_w: b / b_w))
    sh = result.shape
    if sh:
        z = data.nuclide(_excited(exp)).Z
        r0 = 1.2 * a_mass ** (1 / 3)
        out.append(Explanation(
            "beta2", "Deformation β₂ from B(E2↑)", sh["beta2"], "",
            "β₂ = (4π / (3 Z R₀²)) √(B(E2↑)/e²), R₀ = 1.2 A^⅓ fm",
            f"β₂ = 4π / (3 × {z} × {_g(r0, 4)}² fm²) × √{_g(result.b_e2fm4)} fm² = {_g(sh['beta2'], 3)}",
            "The size of the quadrupole deformation: how far the nucleus's charge is from a sphere, as a fraction "
            "of its radius. B(E2) gives only the size; whether the nucleus is prolate or oblate is the sign of the "
            "quadrupole moment, which first-order excitation does not measure.",
            "A uniformly charged, axially symmetric rotor.",
            f"{DOCS}nuclear-setup.html#automatic-analysis",
            {"z": z, "r0": r0, "b": result.b_e2fm4}, lambda z, r0, b: 4 * math.pi / (3 * z * r0**2) * math.sqrt(b)))
        out.append(Explanation(
            "q0", "Intrinsic quadrupole moment of a rotor", sh["q0_efm2"], "e fm²",
            "B(E2; 0 → 2) = (5/16π) e² Q₀²  →  Q₀ = √(16π B(E2↑) / 5) / e",
            f"Q₀ = √(16π × {_g(result.b_e2fm4)} / 5) = {_g(sh['q0_efm2'], 4)} e fm²;  Q_s(2⁺) = −(2/7) Q₀ = "
            f"{_g(sh['qs_2plus_efm2'], 4)} e fm²",
            "The quadrupole moment of the nucleus in its own frame if it were a rigid rotor, and the spectroscopic "
            "moment of its 2⁺ state in the laboratory. A measured Q_s(2⁺) of this size and sign (negative) means "
            "a prolate rotor; the opposite sign, oblate. The sign needs the reorientation effect (item 58).",
            "The rotational model: Q_s(2⁺) = −(2/7) Q₀ for K = 0.",
            f"{DOCS}nuclear-setup.html#level-schemes",
            {"b": result.b_e2fm4}, lambda b: math.sqrt(16 * math.pi * b / 5)))
        b_down = result.b_e2fm4 / 5
        tau = 1 / (rate_per_b(2, 1e3 * exc.energy_mev) * b_down)
        out.append(Explanation(
            "lifetime", "Mean lifetime of the 2⁺ state from B(E2)", tau * 1e12, "ps",
            "τ = 1 / T(E2), T(E2) = (8π(λ+1) / (λ[(2λ+1)!!]²)) (e²/ħ) (E_γ/ħc)^5 B(E2↓), B(E2↓) = B(E2↑)/(2J_f+1)",
            f"B(E2↓) = {_g(result.b_e2fm4)} / 5 = {_g(b_down)} e²fm⁴;  T = {_g(rate_per_b(2, 1e3 * exc.energy_mev), 4)} "
            f"/s per e²fm⁴ × {_g(b_down)} = {_g(1 / tau, 4)} /s;  τ = {_g(tau * 1e12, 3)} ps",
            "A strong B(E2) means a short-lived state: the same matrix element that lets the field excite the "
            "state lets the state decay. The lifetime decides whether the nucleus decays in flight (the Doppler "
            "shift) or after stopping.",
            "No internal conversion, no other branch; the 2⁺ → 0⁺ transition only.",
            f"{DOCS}nuclear-setup.html#level-schemes",
            {"rate_per_b": rate_per_b(2, 1e3 * exc.energy_mev), "b_down": b_down},
            lambda rate_per_b, b_down: 1e12 / (rate_per_b * b_down)))
    out.append(Explanation(
        "uncertainty", "Total relative uncertainty of B(E2↑)", result.total_unc, "",
        "δB/B = √(δ_stat² + Σ δ_sys,i²)",
        f"δB/B = √({_g(result.statistical, 3)}² + " + " + ".join(f"{_g(v, 3)}²" for v in result.budget.values())
        + f") = {_g(result.total_unc, 3)}",
        "The statistical part comes from the counts in the peak and in the normalisation; each systematic part is "
        "what B(E2) moves by when one input is changed by its own uncertainty. Added in quadrature as independent.",
        "The systematic inputs are independent and their effects linear.",
        f"{DOCS}nuclear-setup.html#automatic-analysis",
        {"stat": result.statistical, "sys": list(result.budget.values())},
        lambda stat, sys: math.sqrt(stat**2 + sum(v**2 for v in sys))))
    return out


def _excited(experiment) -> str:
    from .analysis import ex_nuclide

    return ex_nuclide(experiment)


# -- the record -----------------------------------------------------------------------------------------------------


def provenance(experiment) -> list:
    """Every piece of nuclear data the plan uses and where it comes from: (quantity, value, source)."""
    rows = [("Masses", "AME2020", "data"), ("Stopping powers", "NIST PSTAR/ASTAR and tables", "data"),
            ("γ-ray attenuation", "NIST XCOM", "data")]
    exc = experiment.excitation
    if exc is not None:
        rows.append((f"State at {1e3 * exc.energy_mev:g} keV, B({exc.multipolarity}↑)", exc.b_up, "user"))
    for role, scheme in experiment.levels.items():
        for i, lev in enumerate(scheme.levels):
            rows.append((f"{role}: level {i} {lev.label} energy", f"{lev.energy.value:g} keV", lev.energy.source))
        for me in scheme.matrix_elements:
            rows.append((f"{role}: ⟨{me.b}‖M({me.multipolarity})‖{me.a}⟩", f"{me.value.value:.4g} e fm^λ",
                         me.value.source + (f" ({me.value.note})" if me.value.note else "")))
        if scheme.reference:
            rows.append((f"{role}: level scheme", scheme.reference, "ensdf"))
    for gd in experiment.gamma_detectors:
        src = "user" if (gd.efficiency is not None or gd.efficiency_curve) else "assumed (typical response)"
        rows.append((f"Efficiency of {gd.name}", "measured" if src == "user" else "typical response", src))
    return rows


LEFT_OUT = (
    "Excitation in more than one step and the reorientation effect (first-order theory only).",
    "Deorientation; and lifetimes only where the level scheme gives a half-life (the nucleus then slows down along "
    "its path and decays on the way, with its orientation unchanged; no angular straggling).",
    "Cascades: followed through the level scheme when the setup has one (one γ ray per excitation otherwise); the "
    "angular correlation between successive γ rays of a cascade is not included.",
    "Pile-up, summing, and reactions on contaminants of the target.",
    "Photon transport inside the crystals: the response is parametrised and typical unless measured.",
)


def record_html(planner, result=None, title: Optional[str] = None, scene_png: Optional[str] = None) -> str:
    """The run record as one self-contained HTML page: the setup (with the scene if ``scene_png`` is given as
    base64), the data and its provenance, the method, every number with its explanation, the analysis chain, the
    uncertainty budget and what is left out. ``result`` is a :class:`~physim.nuclear.analysis.Result`."""
    exp = planner.experiment
    e = html.escape
    parts = []
    parts.append(f"<h1>{e(title or exp.title)}: run record</h1>")
    b, t = exp.beam, exp.target
    parts.append("<h2>1. The setup</h2><ul>"
                 f"<li>Beam: {e(b.nuclide)} at {e(str(b.energy))}, {e(str(b.current))}.</li>"
                 f"<li>Target: {e(t.material)}, {e(str(t.thickness))}.</li>"
                 f"<li>Particle detectors: {e(', '.join(d.name or f'D{i + 1}' for i, d in enumerate(exp.detectors)))}.</li>"
                 f"<li>γ-ray detectors: {e(', '.join(g.name or f'G{i + 1}' for i, g in enumerate(exp.gamma_detectors))) or 'none'}.</li>"
                 + (f"<li>Reaction: Coulomb excitation of the {e(exp.excitation.excite)} to "
                    f"{1e3 * exp.excitation.energy_mev:g} keV ({e(exp.excitation.multipolarity)}).</li>"
                    if exp.excitation else "<li>Reaction: elastic scattering.</li>") + "</ul>")
    if scene_png:
        parts.append(f'<p><img src="data:image/png;base64,{scene_png}" alt="the scene" style="max-width:100%"></p>')
    parts.append("<h2>2. The nuclear data and where it comes from</h2>")
    parts.append("<table><thead><tr><th>Quantity</th><th>Value</th><th>Source</th></tr></thead><tbody>"
                 + "".join(f"<tr><td>{e(q)}</td><td>{e(str(v))}</td><td>{e(s)}</td></tr>"
                           for q, v, s in provenance(exp)) + "</tbody></table>"
                 "<p>Sources: <em>ensdf</em> is the evaluated data as read; <em>derived</em> is computed from "
                 "ENSDF's values by physim; <em>assumed</em> is a typical value, not a measurement; <em>user</em> "
                 "was given in the setup.</p>")
    parts.append("<h2>3. The method, step by step</h2><ol>"
                 "<li><b>Scattering.</b> The beam particle follows a Rutherford orbit past the target nucleus; "
                 "the scattering angle fixes how close it comes.</li>"
                 "<li><b>Excitation.</b> The passing charge's field lifts the nucleus to the state with "
                 "probability P(θ), proportional to B(Eλ↑), and leaves it oriented.</li>"
                 "<li><b>Detection of the particle.</b> The scattered beam or the recoil reaches a ring and sector "
                 "of a silicon detector, losing energy on the way out of the target and in the dead layer.</li>"
                 "<li><b>Decay.</b> The nucleus leaves the target and emits its γ ray in flight, with the angular "
                 "correlation its orientation gives, Doppler-shifted by its motion.</li>"
                 "<li><b>Detection of the γ ray.</b> A crystal records the full energy, a Compton deposit or an "
                 "escape peak, with its resolution; chance coincidences and the room background join the "
                 "spectrum.</li>"
                 "<li><b>Doppler correction.</b> From the segment and crystal that fired, and two-body kinematics, "
                 "the energy is corrected back to the nucleus's frame.</li></ol>")
    parts.append("<h2>4. Every number, explained</h2>")
    for x in explanations(planner):
        parts.append(_explanation_html(x))
    from .logbook import run_explanations

    taken = run_explanations(planner) if getattr(planner, "run", None) is not None else []
    if taken:
        parts.append(f"<h3>Run {planner.run.number}: its counters, its gates, and the Plan against what it "
                     "measured</h3>")
        for x in taken:
            parts.append(_explanation_html(x))
    if result is not None:
        parts.append("<h2>5. From the peak area to B(E2) and the shape</h2>")
        for x in analysis_explanations(planner, result):
            parts.append(_explanation_html(x))
        parts.append("<h2>6. The uncertainty budget</h2>"
                     "<table><thead><tr><th>Source</th><th>Relative (%)</th></tr></thead><tbody>"
                     f"<tr><td>statistical</td><td>{100 * result.statistical:.2f}</td></tr>"
                     + "".join(f"<tr><td>{e(k)}</td><td>{100 * v:.2f}</td></tr>" for k, v in result.budget.items())
                     + f"<tr><td><b>total</b></td><td><b>{100 * result.total_unc:.2f}</b></td></tr></tbody></table>"
                     + (f"<p>The Monte Carlo sample itself adds {100 * result.monte_carlo:.1f}%, which more "
                        "simulated events reduce. "
                        if getattr(result, "whole_run_statistical", None) is None else
                        f"<p>The statistics come from {html.escape(result.statistics_from)}; the whole run would give "
                        f"{100 * result.whole_run_statistical:.1f}%. ")
                     + f"Against the value put in, the result pulls by {result.pull:+.2f} standard deviations.</p>")
        parts.append("<h2>7. What the simulation leaves out</h2>")
    else:
        parts.append("<h2>5. What the simulation leaves out</h2>")
    parts.append("<ul>" + "".join(f"<li>{e(x)}</li>" for x in LEFT_OUT) + "</ul>")
    return ("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>"
            f"{e(title or exp.title)}: run record</title><style>{record_css()}</style></head><body>"
            + "\n".join(parts) + "</body></html>")


def record_css(scope: str = "body") -> str:
    """The record's style, for the page itself (``scope`` "body") or for a box of the app that holds its body."""
    s = scope
    return (f"{s}{{font-family:Georgia,serif;max-width:900px;margin:0 auto;line-height:1.45}}"
            f"{s} h1{{font-size:1.6em}}{s} h2{{font-size:1.25em;margin-top:1.2em}}{s} h3{{font-size:1.05em}}"
            f"{s} table{{border-collapse:collapse;margin:0.5em 0}}{s} td,{s} th{{border:1px solid #ccc;"
            f"padding:3px 8px;text-align:left}}{s} .x{{border:1px solid #D5D9D6;border-radius:8px;"
            f"padding:8px 14px;margin:10px 0}}{s} .x h3{{margin:0 0 4px}}{s} .x .v{{font-family:monospace;"
            f"font-size:1.1em}}{s} .x .f{{font-family:monospace;background:rgba(127,127,127,0.12);padding:4px 8px;"
            f"margin:4px 0}}{s} .x .m{{margin:4px 0}}{s} .x .a{{opacity:0.75;font-size:0.92em}}"
            + ("@media print{body{max-width:none}}" if scope == "body" else ""))



def _explanation_html(x: Explanation) -> str:
    e = html.escape
    unit = f" {e(x.unit)}" if x.unit else ""
    return (f'<div class="x" id="{e(x.key)}"><h3>{e(x.title)}: <span class="v">{_g(x.value)}{unit}</span></h3>'
            f'<div class="f">{e(x.formula)}</div><div class="f">{e(x.substituted)}</div>'
            f'<div class="m">{e(x.meaning)}</div>'
            f'<div class="a">Assumptions: {e(x.assumptions)} <a href="{e(x.reference)}">Theory</a></div></div>')


__all__ = ["DOCS", "Explanation", "LEFT_OUT", "analysis_explanations", "explanations", "provenance", "record_css",
           "record_html"]
