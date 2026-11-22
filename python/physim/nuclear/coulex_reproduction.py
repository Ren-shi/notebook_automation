"""A published Coulomb-excitation experiment reproduced with physim: Allmond et al., PRC 90, 034309 (2014), ⁶²Ni.

**The experiment.** J. M. Allmond et al., "High-precision B(E2) measurements of semi-magic ⁵⁸,⁶⁰,⁶²,⁶⁴Ni by Coulomb
excitation", Phys. Rev. C 90, 034309 (2014), at the HRIBF (Oak Ridge). A ⁶²Ni beam at 1.8 MeV/u (111.6 MeV, ~5 pA)
on ~1 mg/cm² of natural carbon, in inverse kinematics: the carbon recoils are counted in two full rings of the
BareBall CsI array (ring 2 at 14–28°, ring 3 at 28–44° lab) and the γ rays in nine CLARION clovers at 21.75 cm.
GOSIA turns σ_Coulex/σ_Rutherford of each ring into <0⁺‖E2‖2⁺>: 0.3038(37) e b in ring 2 and 0.3002(24) e b in
ring 3 (statistical only), 0.301(6) e b in all, B(E2; 0⁺ → 2⁺) = 0.0906(37) e²b² (Table I). The setup is the
example ``coulex_ni62_c_allmond`` (what is the paper's and what is assumed is in its description); the paper's
numbers, the GOSIA references and the bands are ``tests/reference/nuclear/coulex/allmond2014.toml``.

**The paper gives no yields.** Its σ_C/σ_R per ring are not printed, only the matrix elements GOSIA made of them.
So check 4 rebuilds each ring's ratio from the paper's matrix element through GOSIA's own ratio(<0‖E2‖2>) curve
(the outside tools' fit, R = a M² + b M⁴, to ten GOSIA runs of this setup) and turns it back into a B(E2) with
physim's curve. This tests whether physim extracts the same B(E2) from the same data as GOSIA did; it is not an
independent measurement of B(E2).

**The checks** (bands fixed before the run, with their reasons, in the reference file; :func:`run`)::

    from physim.nuclear import coulex_reproduction as cr

    r = cr.run()          # the five checks; r["passed"], r["checks"][name]["passed"], r["detail"]
    print(cr.table(r))    # one line per number: physim, reference, band, pass/FAIL

1. *Thickness*: physim's stopping through the example's carbon gives the measured ⁵⁸Ni and ⁶⁰Ni energy losses
   (42.7(8) and 42.1(8) MeV at 1.8 MeV/u) within their errors, and 42.4(8) MeV for ⁶²Ni (SRIM and LISE++ ATIMA,
   each tuned the same way, give 42.39 and 42.42 MeV).
2. *physim against GOSIA*: σ_C/σ_R(elastic) per ring, integrated over the target, at <0⁺‖E2‖2⁺> = 0.29734 e b,
   Q = 0, to all orders. Within 0.5 % of GOSIA run with the nuclear masses and without E1 polarisation (the two
   simplifications of GOSIA's that physim does not make), within 2 % of GOSIA as the outside tools ran it (integer
   mass numbers, DIP 3.8; these move the ratio by about −0.6 % and −0.5 %).
3. *Closure*: the example simulated and analysed (gate "all"): B(E2) back per ring and for both rings within 1 %
   + 3 Monte Carlo standard deviations (the 1 % is the known residual of the analysis: the γ-ray correlation table
   takes the kinematics of the middle of the target).
4. *The paper*: per ring, the paper's matrix element → GOSIA's ratio → physim's B(E2), within the paper's total
   uncertainty of 4.1 % of its per-ring B(E2); the pulls against the per-ring statistical errors are reported.
5. *Sensitivities*: Q_s(2⁺) = ±0.1 e b and the feeding from ENSDF's higher states (0⁺₂, 2⁺₂, 4⁺₁) move the ratio
   as in GOSIA, within 0.5 percentage points.

Results are dicts of plain numbers: ``checks`` (one dict per check with its values, bands and ``passed``),
``passed``, ``worst`` (the largest deviation in units of its band) and ``detail`` (a sentence).
"""

from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Optional

import numpy as np

#: The example setup made after the paper.
EXAMPLE = "coulex_ni62_c_allmond"
#: The paper's numbers, the GOSIA references and the bands, inside the reference folder.
REFERENCE_FILE = "coulex/allmond2014.toml"
RINGS = ("Ring2", "Ring3")
#: The ring keys of the reference file.
REF_RING = {"Ring2": "ring2", "Ring3": "ring3"}
#: The example the ENSDF level scheme of ⁶²Ni (0⁺₂, 2⁺₂, 4⁺₁ and their matrix elements) is taken from.
ENSDF_SCHEME_EXAMPLE = "coulex_ni62_pt194"
#: Default run size of the closure: about 25 s; the Monte Carlo uncertainty of B(E2) is about 1.2 % (ring 2),
#: 0.7 % (ring 3) and 0.6 % (both).
EVENTS = 500_000
SEED = 11

_CACHE: dict = {}


def reference(ref=None) -> Optional[dict]:
    """The reference file (a dict), or None if the reference folder or the file is not there."""
    from ._toml import loads
    from .validation import reference_dir

    base = Path(ref) if ref is not None else reference_dir()
    if base is None or not (base / REFERENCE_FILE).is_file():
        return None
    return loads((base / REFERENCE_FILE).read_text(encoding="utf-8"))


# -- setups ------------------------------------------------------------------------------------------------------


def _two_levels(me_eb: float, q_eb: float) -> dict:
    """[levels.beam] of ⁶²Ni with its 0⁺ and 2⁺₁ only: <0⁺‖E2‖2⁺> = ``me_eb`` and Q_s(2⁺) = ``q_eb`` (the diagonal
    element is set even for Q = 0, so that the coupled equations solve it to all orders)."""
    from .levels import quadrupole_factor

    return {
        "nuclide": "62Ni",
        "level": [{"energy": "0 keV", "jpi": "0+"},
                  {"energy": "1172.91 keV", "jpi": "2+", "half_life": "1.442e-12 s"}],
        "transition": [{"from": 1, "to": 0, "energy": "1172.89 keV", "multipolarity": "E2"}],
        "matrix_element": [
            {"from": 0, "to": 1, "multipolarity": "E2", "value": f"{100 * me_eb:.10g} efm2"},
            {"from": 1, "to": 1, "multipolarity": "E2", "value": f"{100 * q_eb / quadrupole_factor(2.0):.10g} efm2"},
        ],
    }


def _ensdf_levels(me_eb: float) -> dict:
    """[levels.beam] of ⁶²Ni up to the 4⁺₁ at 2336 keV from ENSDF (the example ``coulex_ni62_pt194``; the 0⁺₃ at
    2891 keV, not in the GOSIA run it is compared with, left out), <0⁺‖E2‖2⁺> = ``me_eb`` and Q_s(2⁺) = 0, as in
    the outside tools' GOSIA run with the higher states."""
    from .experiment import Experiment

    beam = copy.deepcopy(Experiment.example(ENSDF_SCHEME_EXAMPLE).to_dict()["levels"]["beam"])
    keep = 5
    beam["level"] = beam["level"][:keep]
    beam["transition"] = [t for t in beam["transition"] if t["from"] < keep and t["to"] < keep]
    beam["matrix_element"] = [m for m in beam["matrix_element"] if m["from"] < keep and m["to"] < keep]
    for m in beam["matrix_element"]:
        pair = (m["from"], m["to"])
        if pair in ((0, 1), (1, 1)):
            m["value"] = f"{100 * me_eb:.10g} efm2" if pair == (0, 1) else "0 efm2"
            m.pop("value_unc", None)
            m["source"] = "user"
            m.pop("note", None)
    return beam


def setup(me_eb: Optional[float] = None, q_eb: Optional[float] = None, levels: Optional[str] = None):
    """The example, optionally with <0⁺‖E2‖2⁺> = ``me_eb`` (e b; B(E2↑) = M² in the [reaction]) and a level scheme
    for the beam: ``levels`` "two" (0⁺, 2⁺ with Q_s = ``q_eb``, 0 if not given; all orders) or "ensdf" (up to the
    4⁺₁, Q_s = 0)."""
    from .experiment import Experiment

    d = Experiment.example(EXAMPLE).to_dict()
    if me_eb is None:
        me_eb = math.sqrt(Experiment.example(EXAMPLE).excitation.b_up_e2fm) / 100.0
    d["reaction"]["b_up"] = f"{(100 * me_eb) ** 2:.10g} e2fm4"
    if levels == "two":
        d["levels"] = {"beam": _two_levels(me_eb, q_eb or 0.0)}
    elif levels == "ensdf":
        d["levels"] = {"beam": _ensdf_levels(me_eb)}
    elif levels is not None:
        raise ValueError(f"levels must be None, 'two' or 'ensdf', not {levels!r}")
    return Experiment.from_dict(d)


def ratios(exp) -> dict:
    """σ_Coulex/σ_Rutherford(elastic) per ring, integrated over the target (the rates, every particle counted)."""
    from .rates import Rates

    rt = Rates(exp)
    return {r: rt.rate(r, counted=False, what="excitations") / rt.rate(r, counted=False) for r in RINGS}


def ratio_fit(points: dict) -> tuple:
    """(a, b) of R = a M² + b M⁴ through {M: R} by least squares (M in e b), the outside tools' form."""
    m = np.array(sorted(points), dtype=float)
    r = np.array([points[k] for k in sorted(points)], dtype=float)
    a, b = np.linalg.lstsq(np.c_[m**2, m**4], r, rcond=None)[0]
    return float(a), float(b)


def ratio_of(a: float, b: float, me_eb: float) -> float:
    return a * me_eb**2 + b * me_eb**4


def me_of(a: float, b: float, ratio: float) -> float:
    """The inverse of :func:`ratio_of`: M = sqrt[(−a + sqrt(a² + 4bR)) / (2b)]."""
    if b == 0:
        return math.sqrt(ratio / a)
    return math.sqrt((-a + math.sqrt(a * a + 4 * b * ratio)) / (2 * b))


# -- the checks --------------------------------------------------------------------------------------------------


def _within(value: float, want: float, band: float) -> tuple:
    dev = value - want
    return abs(dev) <= band, abs(dev) / band


def check_thickness(refd: dict) -> dict:
    from .experiment import Experiment
    from .stopping import Stopping

    from .rates import stack

    exp = Experiment.example(EXAMPLE)
    thick = float(stack(exp)[0].thickness)  # mg/cm²
    material = exp.target.material
    t = refd["paper"]["energy_loss"]
    rows, worst, ok = [], 0.0, True
    for ion in ("58Ni", "60Ni", "62Ni"):
        e_in = 1.8 * int(ion[:2])
        loss = e_in - float(Stopping(ion, material).energy_after(e_in, thick))
        want, unc = t[ion]["value_mev"], t[ion]["unc_mev"]
        passed, w = _within(loss, want, unc)
        rows.append({"ion": ion, "entry_mev": e_in, "loss_mev": loss, "reference_mev": want, "band_mev": unc,
                     "source": t[ion]["source"], "passed": passed})
        worst, ok = max(worst, w), ok and passed
    return {"thickness_mg_cm2": float(thick), "material": material, "rows": rows, "passed": ok, "worst": worst}


def check_gosia(refd: dict, base: dict) -> dict:
    g = refd["gosia"]
    bands = refd["bands"]["gosia"]
    rows, worst, ok = [], 0.0, True
    for setting, band in (("matched", bands["matched_rel"]), ("as_run", bands["as_run_rel"])):
        for r in RINGS:
            want = g[setting][REF_RING[r]]["ratio_elastic"]
            dev = base[r] / want - 1
            passed = abs(dev) <= band
            rows.append({"setting": setting, "ring": r, "physim": base[r], "gosia": want, "deviation": dev,
                         "band": band, "passed": passed})
            worst, ok = max(worst, abs(dev) / band), ok and passed
    return {"rows": rows, "passed": ok, "worst": worst}


def check_closure(refd: dict, events: int, seed: int) -> dict:
    from .analysis import Analysis, Settings
    from .experiment import Experiment
    from .gamma_events import simulate_gammas

    band = refd["bands"]["closure"]
    exp = Experiment.example(EXAMPLE)
    b_in = exp.excitation.b_up_e2fm
    a = Analysis(exp, simulate_gammas(exp, events, seed=seed))
    rows, worst, ok = [], 0.0, True
    for det in ("Ring2", "Ring3", None):
        res = a.run(Settings(detectors=None if det is None else [det], particle_gate="all"))
        dev = res.b_e2fm4 / b_in - 1
        allowed = band["rel"] + band["mc_sigmas"] * res.monte_carlo
        passed = abs(dev) <= allowed
        rows.append({"rings": det or "both", "b_e2fm4": res.b_e2fm4, "b_in_e2fm4": b_in, "deviation": dev,
                     "monte_carlo": res.monte_carlo, "band": allowed, "passed": passed})
        worst, ok = max(worst, abs(dev) / allowed), ok and passed
    return {"events": events, "seed": seed, "rows": rows, "passed": ok, "worst": worst}


def check_paper(refd: dict, scan: dict) -> dict:
    """The paper's per-ring matrix element through GOSIA's ratio(M) and back through physim's: B(E2) per ring
    against the paper's, within the paper's total relative uncertainty of B(E2)."""
    p = refd["paper"]
    fits = refd["gosia"]["fit"]
    band = refd["bands"]["paper"]["rel"]
    rows, worst, ok = [], 0.0, True
    wsum = bsum = 0.0
    for r in RINGS:
        me_paper = p[REF_RING[r]]["me_eb"]
        me_stat = p[REF_RING[r]]["me_stat_eb"]
        f = fits[REF_RING[r]]
        ratio = ratio_of(f["a"], f["b"], me_paper)
        a, b = ratio_fit({m: v[r] for m, v in scan.items()})
        me_app = me_of(a, b, ratio)
        b_paper, b_app = me_paper**2, me_app**2
        sigma = 2 * me_stat / me_paper * b_paper  # statistical, on M² = B(E2↑)
        dev = b_app / b_paper - 1
        passed = abs(dev) <= band
        rows.append({"ring": r, "me_paper_eb": me_paper, "ratio_from_gosia": ratio, "physim_a": a, "physim_b": b,
                     "me_physim_eb": me_app, "b_paper_e2b2": b_paper, "b_physim_e2b2": b_app, "deviation": dev,
                     "band": band, "pull_stat": (b_app - b_paper) / sigma, "stat_rel": sigma / b_paper,
                     "passed": passed})
        worst, ok = max(worst, abs(dev) / band), ok and passed
        wsum += 1 / sigma**2
        bsum += b_app / sigma**2
    final = p["final"]
    combined = bsum / wsum
    return {"rows": rows, "combined_b_e2b2": combined, "final_b_e2b2": final["b_e2b2"],
            "final_unc_e2b2": final["b_unc_e2b2"], "combined_deviation": combined / final["b_e2b2"] - 1,
            "passed": ok, "worst": worst}


def check_sensitivities(refd: dict, base: dict, q_plus: dict, q_minus: dict, feeding: dict) -> dict:
    g = refd["gosia"]["sensitivity"]
    band = refd["bands"]["sensitivity"]["pp"]
    rows, worst, ok = [], 0.0, True
    for what, ours in (("q_plus", {r: q_plus[r] / base[r] - 1 for r in RINGS}),
                       ("q_minus", {r: q_minus[r] / base[r] - 1 for r in RINGS}),
                       ("feeding", feeding)):
        for r in RINGS:
            want = g[REF_RING[r]][what]
            dev_pp = 100 * (ours[r] - want)
            passed = abs(dev_pp) <= band
            rows.append({"what": what, "ring": r, "physim_percent": 100 * ours[r], "gosia_percent": 100 * want,
                         "difference_pp": dev_pp, "band_pp": band, "passed": passed})
            worst, ok = max(worst, abs(dev_pp) / band), ok and passed
    return {"rows": rows, "passed": ok, "worst": worst}


def feeding_effect(angle_step: float = 2.0, energies: int = 6) -> dict:
    """The change of the 2⁺₁ → 0⁺ γ yield per ring when ENSDF's higher states feed it, against the 0⁺, 2⁺ scheme
    alone (both Q = 0, all orders, :class:`~physim.nuclear.multistep.Multistep` through the target). Relative."""
    from .multistep import Multistep

    out = {}
    got = {}
    for name in ("two", "ensdf"):
        exp = setup(levels=name)
        y = Multistep(exp, exp.levels["beam"], "beam", angle_step=angle_step, energies=energies).yields()
        got[name] = {r: y.detectors[r][(1, 0)] for r in RINGS}
    for r in RINGS:
        out[r] = got["ensdf"][r] / got["two"][r] - 1
    return out


def run(ref=None, events: int = EVENTS, seed: int = SEED) -> Optional[dict]:
    """The five checks (cached per reference folder and run size); None if the reference file is not there."""
    refd = reference(ref)
    if refd is None:
        return None
    key = (str(ref), events, seed)
    if key in _CACHE:
        return _CACHE[key]
    me0 = refd["gosia"]["me_eb"]
    scan = {}
    for m in refd["bands"]["paper"]["physim_scan_me_eb"]:
        scan[m] = ratios(setup(me_eb=m, levels="two"))
    base = scan[me0] if me0 in scan else ratios(setup(me_eb=me0, levels="two"))
    q = refd["gosia"]["sensitivity"]["q_eb"]
    q_plus = ratios(setup(me_eb=me0, q_eb=q, levels="two"))
    q_minus = ratios(setup(me_eb=me0, q_eb=-q, levels="two"))
    checks = {
        "thickness": check_thickness(refd),
        "gosia": check_gosia(refd, base),
        "closure": check_closure(refd, events, seed),
        "paper": check_paper(refd, scan),
        "sensitivities": check_sensitivities(refd, base, q_plus, q_minus, feeding_effect()),
    }
    out = {"checks": checks, "ratios_q0": base, "scan": scan, "passed": all(c["passed"] for c in checks.values()),
           "worst": max(c["worst"] for c in checks.values())}
    out["detail"] = summary(out)
    _CACHE[key] = out
    return out


# -- what to show ------------------------------------------------------------------------------------------------


def summary(r: dict) -> str:
    c = r["checks"]
    th = {x["ion"]: x for x in c["thickness"]["rows"]}
    g = {(x["setting"], x["ring"]): x for x in c["gosia"]["rows"]}
    cl = {x["rings"]: x for x in c["closure"]["rows"]}
    pa = {x["ring"]: x for x in c["paper"]["rows"]}
    se = {(x["what"], x["ring"]): x for x in c["sensitivities"]["rows"]}

    def ok(name):
        return "pass" if c[name]["passed"] else "FAIL"

    return (
        f"1 thickness ({ok('thickness')}): losses 58Ni {th['58Ni']['loss_mev']:.2f}, 60Ni {th['60Ni']['loss_mev']:.2f},"
        f" 62Ni {th['62Ni']['loss_mev']:.2f} MeV (42.7(8), 42.1(8), 42.4(8)); "
        f"2 GOSIA ({ok('gosia')}): ring 2/3 {100 * g['matched', 'Ring2']['deviation']:+.2f}/"
        f"{100 * g['matched', 'Ring3']['deviation']:+.2f} % against nuclear masses, DIP off (0.5 %), "
        f"{100 * g['as_run', 'Ring2']['deviation']:+.2f}/{100 * g['as_run', 'Ring3']['deviation']:+.2f} % against "
        f"GOSIA as run (2 %); "
        f"3 closure ({ok('closure')}): {100 * cl['Ring2']['deviation']:+.2f}/{100 * cl['Ring3']['deviation']:+.2f}/"
        f"{100 * cl['both']['deviation']:+.2f} % (MC {100 * cl['Ring2']['monte_carlo']:.2f}/"
        f"{100 * cl['Ring3']['monte_carlo']:.2f}/{100 * cl['both']['monte_carlo']:.2f} %); "
        f"4 paper ({ok('paper')}, the paper's ratios rebuilt through GOSIA): B(E2) "
        f"{pa['Ring2']['b_physim_e2b2']:.4f}/{pa['Ring3']['b_physim_e2b2']:.4f} e2b2 against "
        f"{pa['Ring2']['b_paper_e2b2']:.4f}/{pa['Ring3']['b_paper_e2b2']:.4f} ({100 * pa['Ring2']['deviation']:+.1f}/"
        f"{100 * pa['Ring3']['deviation']:+.1f} %, band 4.1 %; pulls {pa['Ring2']['pull_stat']:+.2f}/"
        f"{pa['Ring3']['pull_stat']:+.2f} stat), both {c['paper']['combined_b_e2b2']:.4f} against 0.0906(37); "
        f"5 sensitivities ({ok('sensitivities')}): Q +0.1 eb {se['q_plus', 'Ring2']['physim_percent']:+.2f}/"
        f"{se['q_plus', 'Ring3']['physim_percent']:+.2f} % (GOSIA {se['q_plus', 'Ring2']['gosia_percent']:+.2f}/"
        f"{se['q_plus', 'Ring3']['gosia_percent']:+.2f}), Q −0.1 eb {se['q_minus', 'Ring2']['physim_percent']:+.2f}/"
        f"{se['q_minus', 'Ring3']['physim_percent']:+.2f} % (GOSIA {se['q_minus', 'Ring2']['gosia_percent']:+.2f}/"
        f"{se['q_minus', 'Ring3']['gosia_percent']:+.2f}), feeding {se['feeding', 'Ring2']['physim_percent']:+.2f}/"
        f"{se['feeding', 'Ring3']['physim_percent']:+.2f} % (GOSIA {se['feeding', 'Ring2']['gosia_percent']:+.2f}/"
        f"{se['feeding', 'Ring3']['gosia_percent']:+.2f})")


def table(r: dict) -> str:
    """A plain-text table, one line per compared number."""
    c = r["checks"]
    lines = [f"{'check':14s} {'what':40s} {'physim':>12s} {'reference':>12s} {'deviation':>11s} {'band':>10s}  result"]

    def line(check, what, ours, theirs, dev, band, passed):
        lines.append(f"{check:14s} {what:40s} {ours:>12s} {theirs:>12s} {dev:>11s} {band:>10s}  "
                     f"{'pass' if passed else 'FAIL'}")

    for x in c["thickness"]["rows"]:
        line("1 thickness", f"{x['ion']} loss at {x['entry_mev']:.1f} MeV", f"{x['loss_mev']:.2f} MeV",
             f"{x['reference_mev']:.1f} MeV", f"{x['loss_mev'] - x['reference_mev']:+.2f} MeV",
             f"+-{x['band_mev']:.1f} MeV", x["passed"])
    names = {"matched": "nuclear masses, DIP off", "as_run": "as run (A=62/12, DIP 3.8)"}
    for x in c["gosia"]["rows"]:
        line("2 GOSIA", f"{x['ring']} R_el vs {names[x['setting']]}", f"{x['physim']:.4e}", f"{x['gosia']:.4e}",
             f"{100 * x['deviation']:+.2f} %", f"+-{100 * x['band']:.1f} %", x["passed"])
    for x in c["closure"]["rows"]:
        line("3 closure", f"B(E2) {x['rings']} (MC {100 * x['monte_carlo']:.2f} %)", f"{x['b_e2fm4']:.1f}",
             f"{x['b_in_e2fm4']:.1f}", f"{100 * x['deviation']:+.2f} %", f"+-{100 * x['band']:.2f} %", x["passed"])
    for x in c["paper"]["rows"]:
        line("4 paper", f"{x['ring']} B(E2) e2b2 (pull {x['pull_stat']:+.2f})", f"{x['b_physim_e2b2']:.4f}",
             f"{x['b_paper_e2b2']:.4f}", f"{100 * x['deviation']:+.2f} %", f"+-{100 * x['band']:.1f} %", x["passed"])
    p = c["paper"]
    lines.append(f"{'':14s} {'both rings (stat-weighted) e2b2':40s} {p['combined_b_e2b2']:>12.4f} "
                 f"{p['final_b_e2b2']:>12.4f} {100 * p['combined_deviation']:>+9.2f} % {'(reported)':>10s}")
    labels = {"q_plus": "Q_s = +0.1 eb", "q_minus": "Q_s = -0.1 eb", "feeding": "feeding (ENSDF states)"}
    for x in c["sensitivities"]["rows"]:
        line("5 sensitivity", f"{x['ring']} {labels[x['what']]}", f"{x['physim_percent']:+.2f} %",
             f"{x['gosia_percent']:+.2f} %", f"{x['difference_pp']:+.2f} pp", f"+-{x['band_pp']:.1f} pp",
             x["passed"])
    return "\n".join(lines)
