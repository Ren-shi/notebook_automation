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

**To all orders.** Every number that is judged solves the coupled equations of the 0⁺, 2⁺ scheme (Q = 0), as GOSIA
does. What first order gives, the app's default for a setup without a level scheme, is shown beside it as a second
line and not judged: it is 0.7 % (ring 2) and 0.4 % (ring 3) above all orders in the ratio.

**The checks** (the bands, their reasons and which of them were fixed without knowing the answer are in the
reference file; :func:`run`)::

    from physim.nuclear import coulex_reproduction as cr

    r = cr.run()          # the checks; r["passed"], r["checks"][name]["passed"], r["detail"]
    print(cr.table(r))    # one line per number: physim, reference, band, pass/MISS/FAIL

1. *Thickness*: physim's stopping through the example's carbon gives the measured ⁵⁸Ni and ⁶⁰Ni energy losses
   (42.7(8) and 42.1(8) MeV at 1.8 MeV/u) within their errors, and 42.4(8) MeV for ⁶²Ni (SRIM and LISE++ ATIMA,
   each tuned the same way, give 42.39 and 42.42 MeV).
2. *physim against GOSIA*: σ_C/σ_R(elastic) per ring, integrated over the target, at <0⁺‖E2‖2⁺> = 0.29734 e b,
   Q = 0. Against GOSIA as the outside tools ran it (integer mass numbers, DIP 3.8: the E1 polarisation with d =
   0.0038, which physim has too by default) the band is the 1 % that the verification's plan fixed beforehand. Until
   physim had the E1 polarisation (2026-11-25) it missed that band by a known amount, +1.3 % (ring 2) and +1.2 % (ring
   3), and the whole check failed; with it the numbers are inside, +0.8 % and +0.7 %: GOSIA's integer masses (about
   0.6 %) and its non-relativistic kinematics (0.24 %). The reference file keeps the record of the miss. Against GOSIA
   run with the nuclear masses and without E1 polarisation, physim without it too: within 0.5 %. GOSIA run with the
   nuclear masses and the E1 polarisation is shown as a second line, not judged.
3. *Closure*: the example with the 0⁺, 2⁺ scheme, simulated and analysed (gate "all") with the seeds of the
   reference file, chosen before they were run: the mean B(E2) per ring and for both rings within a quarter of the
   paper's statistical error + 3 Monte Carlo errors of the mean. The scatter of the seeds is reported beside it. The
   scheme has the 2⁺ state's half-life, so some of the nuclei decay before they have left the target, with a tail
   below the peak that the setup without a scheme does not have.
4. *The paper*: per ring, the paper's matrix element → GOSIA's ratio → physim's B(E2), within the paper's total
   uncertainty of 4.1 % of its per-ring B(E2); the pulls against the per-ring statistical errors are reported.
5. *Sensitivities*: Q_s(2⁺) = ±0.1 e b and the feeding from ENSDF's higher states (0⁺₂, 2⁺₂, 4⁺₁) move the ratio
   as in GOSIA, within 0.5 percentage points and within 2 % of GOSIA's change (a constant factor on the effect,
   which the percentage points alone let through: ±2.9 % for ±2.7 % is 6 % too much and 0.2 points).
6. *Stored run*: a closure sample as a run stores it (the segment of each particle, not its angle; the stored
   columns), through the app's own loading: the excitations expected per particle against the rates' integral over
   the ring (0.1 %), and B(E2) back.

Results are dicts of plain numbers: ``checks`` (one dict per check with its values, bands and ``passed``),
``passed`` (every number inside its band), ``known_misses`` and ``only_known_misses`` (a miss of check 2 that the
reference file records as known; none is, since 2026-11-25), ``worst`` (the largest deviation in units of its band, the
misses included), ``worst_other`` (the same without the known misses) and ``detail`` (a sentence).
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
#: The rows of the closure and of the stored run, and the detectors each takes.
ROWS = {"Ring2": ["Ring2"], "Ring3": ["Ring3"], "both": None}

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


def setup(me_eb: Optional[float] = None, q_eb: Optional[float] = None, levels: Optional[str] = None,
          e1_polarisation: Optional[float] = None):
    """The example, optionally with <0⁺‖E2‖2⁺> = ``me_eb`` (e b; B(E2↑) = M² in the [reaction]) and a level scheme
    for the beam: ``levels`` "two" (0⁺, 2⁺ with Q_s = ``q_eb``, 0 if not given; all orders) or "ensdf" (up to the
    4⁺₁, Q_s = 0). ``e1_polarisation`` is the setup's (0.0038 when None, as in GOSIA as run)."""
    from .experiment import Experiment

    d = Experiment.example(EXAMPLE).to_dict()
    if me_eb is None:
        me_eb = math.sqrt(Experiment.example(EXAMPLE).excitation.b_up_e2fm) / 100.0
    d["reaction"]["b_up"] = f"{(100 * me_eb) ** 2:.10g} e2fm4"
    if e1_polarisation is not None:
        d["reaction"]["e1_polarisation"] = e1_polarisation
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


def check_gosia(refd: dict, base: dict, off: dict, first: Optional[dict] = None,
                first_off: Optional[dict] = None) -> dict:
    """The ratios against GOSIA's: ``base`` (E1 polarisation on, as in GOSIA as run) and ``off`` (without, as in the
    matched run), both to all orders. ``first`` and ``first_off`` (first order) are the second line, under
    ``second_line``, with ``base`` against GOSIA's run with nuclear masses and the E1 polarisation. A row outside its
    band has ``known_miss`` when the reference file records a miss (``known_miss`` of its bands) and its size is still
    the recorded one."""
    g = refd["gosia"]
    bands = refd["bands"]["gosia"]
    miss = bands.get("known_miss")
    masses = g["app_input"]["nuclear_masses"]
    rows, second, worst, ok = [], [], 0.0, True
    for setting, band, ours, theirs in (("as_run", bands["as_run_rel"], base, first),
                                        ("matched", bands["matched_rel"], off, first_off)):
        for r in RINGS:
            want = g[setting][REF_RING[r]]["ratio_elastic"]
            dev = ours[r] / want - 1
            passed = abs(dev) <= band
            known = (not passed and miss is not None and setting == "as_run"
                     and abs(dev - miss[REF_RING[r]]) <= miss["size_abs"])
            rows.append({"setting": setting, "ring": r, "order": "all", "physim": ours[r], "gosia": want,
                         "deviation": dev, "band": band, "passed": passed, "known_miss": known,
                         # What GOSIA's own run with the nuclear masses changes, and what is then left.
                         "explained": masses[REF_RING[r]] / want - 1 if setting == "as_run" else 0.0})
            worst, ok = max(worst, abs(dev) / band), ok and passed
            if theirs is not None:
                dev = theirs[r] / want - 1
                second.append({"setting": setting, "ring": r, "order": "first", "physim": theirs[r], "gosia": want,
                               "deviation": dev, "band": band, "inside": abs(dev) <= band})
    for r in RINGS:
        want, band = masses[REF_RING[r]], bands["matched_rel"]
        dev = base[r] / want - 1
        second.append({"setting": "nuclear_masses", "ring": r, "order": "all", "physim": base[r], "gosia": want,
                       "deviation": dev, "band": band, "inside": abs(dev) <= band})
    return {"rows": rows, "second_line": second, "passed": ok, "worst": worst,
            "why_missed": miss["why"] if miss is not None else ""}


def _closure_band(refd: dict, row: str) -> float:
    """The part of the closure's band that is not Monte Carlo: a share of the paper's relative statistical error of
    B(E2) = M² in that ring, or of the weighted mean of both."""
    p = refd["paper"]
    stat = {r: 2 * p[REF_RING[r]]["me_stat_eb"] / p[REF_RING[r]]["me_eb"] for r in RINGS}
    both = math.sqrt(1 / sum(1 / (stat[r] * p[REF_RING[r]]["me_eb"] ** 2) ** 2 for r in RINGS))
    stat["both"] = both / p["final"]["b_e2b2"]
    return refd["bands"]["closure"]["share_of_paper_stat"] * stat[row]


def check_closure(refd: dict, events: Optional[int] = None, seeds=None, keep: Optional[dict] = None) -> dict:
    """The mean of B(E2)/input − 1 over the seeds, per ring and for both, events and analysis to all orders. The
    second line (``first_order``) is what the same analyses read in first order. ``keep`` is given the first seed's
    analysis ("analysis") and results ("results"), for the stored run."""
    from .analysis import Analysis, Settings
    from .gamma_events import simulate_gammas

    band = refd["bands"]["closure"]
    events = int(band["events"] if events is None else events)
    seeds = [int(s) for s in (band["seeds"] if seeds is None else seeds)]
    exp = setup(levels="two")
    b_in = exp.excitation.b_up_e2fm
    got: dict = {row: [] for row in ROWS}
    for seed in seeds:
        a = Analysis(exp, simulate_gammas(exp, events, seed=seed))
        for row, dets in ROWS.items():
            got[row].append(a.run(Settings(detectors=dets, particle_gate="all")))
        if keep is not None and "analysis" not in keep:
            keep.update(analysis=a, results={row: got[row][0] for row in ROWS})
    rows, worst, ok = [], 0.0, True
    for row, results in got.items():
        each = np.array([r.b_e2fm4 / b_in - 1 for r in results])
        mc = math.sqrt(sum(r.monte_carlo**2 for r in results)) / len(results)
        scatter = float(np.std(each, ddof=1)) / math.sqrt(len(each)) if len(each) > 1 else math.nan
        dev, allowed = float(each.mean()), _closure_band(refd, row) + band["mc_sigmas"] * mc
        passed = abs(dev) <= allowed
        first = [(r.b_first_order if r.b_first_order is not None else r.b_e2fm4) / b_in - 1 for r in results]
        rows.append({"rings": row, "b_e2fm4": (1 + dev) * b_in, "b_in_e2fm4": b_in, "deviation": dev,
                     "monte_carlo": mc, "scatter": scatter, "each": each.tolist(),
                     "first_order": float(np.mean(first)), "tails": sum(r.fit.tails is not None for r in results),
                     "band": allowed, "passed": passed})
        worst, ok = max(worst, abs(dev) / allowed), ok and passed
    return {"events": events, "seeds": seeds, "rows": rows, "passed": ok, "worst": worst}


def stored_run(a, hours: float, seed: int):
    """A run of ``hours`` made of a weighted sample (an :class:`~physim.nuclear.analysis.Analysis` of one), as the app
    stores and loads a run (:class:`~physim.nuclear.runs.RunData`): the counted particles with the stored columns and
    their precision, at the mean angles of the segments they hit, and the γ rays of each drawn as often as their
    rate says, each counting once, with the Doppler correction worked out again from the stored columns. Only the
    particles differ from a run's: they keep their weights (10⁹ of them would each count once), which no more than
    counts them."""
    import dataclasses

    from . import runs

    ev, g, real = a.events, a.gammas, 3600.0 * hours
    rows = np.flatnonzero(ev.columns["counted"])
    cols = {k: np.asarray(v)[rows] for k, v in ev.columns.items()}
    kept = runs._load_particles(runs._store_particles(cols, 0))
    kept.update(recoil=kept["recoil"].astype(bool), counted=np.ones(len(rows), dtype=bool), weight=cols["weight"])
    kept["theta"], kept["phi"] = runs._segment_directions(a.experiment, kept)
    kept["beam_energy"] = runs._beam_energies(a.experiment, kept["depth"])
    row_of = np.full(len(ev.columns["counted"]), -1)
    row_of[rows] = np.arange(len(rows))
    parent = row_of[g["particle"]]
    times = np.where(parent >= 0, np.random.default_rng(seed).poisson(g["weight"] * real), 0)
    pick = np.repeat(np.arange(len(g)), times)
    arrays = {k: (parent if k == "particle" else np.asarray(g[k]))[pick].astype(kind)
              for k, kind in runs.GAMMA_COLUMNS.items()}
    summary = {"number": 1, "kind": "beam", "status": "complete", "seed": seed, "duration_s": real, "real_s": real,
               "generated": int(ev.n_events), "detectors": list(ev.detectors), "channels": list(ev.channels),
               "particles": len(rows), "gammas": len(pick), "nuclei": 1, "one_particle": True}
    run = runs.RunData(Path("synthetic"), summary, a.experiment.to_dict())
    run._events = dataclasses.replace(ev, columns=kept, beam_time_s=real)
    run._gammas[False] = runs._gamma_events(run, arrays, False)
    return run


def check_stored(refd: dict, kept: dict) -> dict:
    """A closure sample as a stored run (:func:`stored_run`), analysed as the app analyses a run."""
    from .analysis import Analysis, Settings

    band = refd["bands"]["stored"]
    sigmas = refd["bands"]["closure"]["mc_sigmas"]
    a = kept["analysis"]
    run = stored_run(a, band["hours"], band["draw_seed"])
    exp = run.experiment
    b = Analysis(exp, run.gammas(), run=run)
    per_ring = ratios(exp)
    b_in = exp.excitation.b_up_e2fm
    rows, worst, ok = [], 0.0, True
    for row, dets in ROWS.items():
        s = Settings(detectors=dets, particle_gate="all")
        if row in per_ring:
            mine = b._elastic_rows(s)
            w = b.events.columns["weight"][mine]
            share = b._excited_share(mine, mode="forward")
            mean = float(np.average(share, weights=w))
            # The sample's own share of the segments: the error of the mean of the segments' expectations.
            err = math.sqrt(float(np.sum((w * (share - mean)) ** 2))) / float(w.sum()) / mean
            dev, allowed = mean / per_ring[row] - 1, band["probability_rel"] + sigmas * err
            rows.append({"what": "excitations per particle", "rings": row, "physim": mean,
                         "reference": per_ring[row], "deviation": dev, "error": err, "band": allowed,
                         "passed": abs(dev) <= allowed})
        res = b.run(s)
        err = math.hypot(res.statistical, kept["results"][row].monte_carlo)
        dev, allowed = res.b_e2fm4 / b_in - 1, _closure_band(refd, row) + sigmas * err
        first = res.b_first_order if res.b_first_order is not None else res.b_e2fm4
        rows.append({"what": "B(E2)", "rings": row, "physim": res.b_e2fm4, "reference": b_in, "deviation": dev,
                     "error": err, "band": allowed, "passed": abs(dev) <= allowed, "counts": res.area,
                     "weighted": kept["results"][row].b_e2fm4 / b_in - 1, "first_order": first / b_in - 1})
    for x in rows:
        worst, ok = max(worst, abs(x["deviation"]) / x["band"]), ok and x["passed"]
    return {"hours": band["hours"], "rows": rows, "passed": ok, "worst": worst}


def check_paper(refd: dict, scan: dict, first: Optional[dict] = None) -> dict:
    """The paper's per-ring matrix element through GOSIA's ratio(M) and back through physim's to all orders
    (``scan``): B(E2) per ring against the paper's, within the paper's total relative uncertainty of B(E2). The
    same through physim's first-order ratios (``first``) is the second line."""
    out = _paper(refd, scan)
    out["second_line"] = None if first is None else _paper(refd, first)
    return out


def _paper(refd: dict, scan: dict) -> dict:
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
    band, rel = refd["bands"]["sensitivity"]["pp"], refd["bands"]["sensitivity"]["rel"]
    rows, worst, ok = [], 0.0, True
    for what, ours in (("q_plus", {r: q_plus[r] / base[r] - 1 for r in RINGS}),
                       ("q_minus", {r: q_minus[r] / base[r] - 1 for r in RINGS}),
                       ("feeding", feeding)):
        for r in RINGS:
            want = g[REF_RING[r]][what]
            dev_pp = 100 * (ours[r] - want)
            dev_rel = ours[r] / want - 1
            passed = abs(dev_pp) <= band and abs(dev_rel) <= rel
            rows.append({"what": what, "ring": r, "physim_percent": 100 * ours[r], "gosia_percent": 100 * want,
                         "difference_pp": dev_pp, "band_pp": band, "difference_rel": dev_rel, "band_rel": rel,
                         "passed": passed})
            worst, ok = max(worst, abs(dev_pp) / band, abs(dev_rel) / rel), ok and passed
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


def run(ref=None, events: Optional[int] = None, seeds=None) -> Optional[dict]:
    """The six checks (cached per reference folder and closure run); None if the reference file is not there.
    ``events`` and ``seeds`` of the closure are the reference file's by default."""
    refd = reference(ref)
    if refd is None:
        return None
    key = (str(ref), events, None if seeds is None else tuple(seeds))
    if key in _CACHE:
        return _CACHE[key]
    me0 = refd["gosia"]["me_eb"]
    points = [*refd["bands"]["paper"]["physim_scan_me_eb"]]
    scan = {m: ratios(setup(me_eb=m, levels="two")) for m in points}
    first = {m: ratios(setup(me_eb=m)) for m in points}
    base = scan[me0] if me0 in scan else ratios(setup(me_eb=me0, levels="two"))
    base_first = first[me0] if me0 in first else ratios(setup(me_eb=me0))
    off = ratios(setup(me_eb=me0, levels="two", e1_polarisation=0.0))
    first_off = ratios(setup(me_eb=me0, e1_polarisation=0.0))
    q = refd["gosia"]["sensitivity"]["q_eb"]
    q_plus = ratios(setup(me_eb=me0, q_eb=q, levels="two"))
    q_minus = ratios(setup(me_eb=me0, q_eb=-q, levels="two"))
    kept: dict = {}
    checks = {
        "thickness": check_thickness(refd),
        "gosia": check_gosia(refd, base, off, base_first, first_off),
        "closure": check_closure(refd, events, seeds, kept),
        "paper": check_paper(refd, scan, first),
        "sensitivities": check_sensitivities(refd, base, q_plus, q_minus, feeding_effect()),
        "stored": check_stored(refd, kept),
    }
    out = {"checks": checks, "ratios_q0": base, "ratios_first_order": base_first, "ratios_without_e1": off,
           "scan": scan, **outcome(checks)}
    out["detail"] = summary(out)
    _CACHE[key] = out
    return out


def outcome(checks: dict) -> dict:
    """What the checks come to: ``passed`` (every number inside its band), the ``known_misses`` among the rows that
    are not, whether they are all there is (``only_known_misses``), and the worst deviation in units of its band,
    of everything (``worst``) and without the known misses (``worst_other``)."""
    misses = [x for x in checks["gosia"]["rows"] if x["known_miss"]]
    outside = sum(not x["passed"] for c in checks.values() for x in c["rows"])
    others = [abs(x["deviation"]) / x["band"] for x in checks["gosia"]["rows"] if not x["known_miss"]]
    others += [c["worst"] for name, c in checks.items() if name != "gosia"]
    return {"passed": outside == 0, "known_misses": misses, "only_known_misses": outside == len(misses) > 0,
            "worst": max(c["worst"] for c in checks.values()), "worst_other": max(others)}


# -- what to show ------------------------------------------------------------------------------------------------


def verdict(r: dict) -> str:
    """One sentence: passed, or what missed or failed, and the worst deviation in units of its band."""
    worst = f"worst {r['worst']:.2f} of its band"
    if r["passed"]:
        return f"every number inside its band ({worst})"
    other = f"{r['worst_other']:.2f} without the known miss" if r["known_misses"] else ""
    if r["only_known_misses"]:
        sizes = " and ".join(f"{100 * x['deviation']:+.2f} % ({x['ring']})" for x in r["known_misses"])
        return (f"NOT PASSED: a known, explained miss of the 1 % band against GOSIA as run, {sizes}; every other "
                f"number inside its band ({worst}; {other})")
    bad = [name for name, c in r["checks"].items() if any(not x["passed"] and not x.get("known_miss")
                                                          for x in c["rows"])]
    return f"FAILED: {', '.join(bad)} ({worst}{'; ' + other if other else ''})"


def summary(r: dict) -> str:
    c = r["checks"]
    th = {x["ion"]: x for x in c["thickness"]["rows"]}
    g = {(x["setting"], x["ring"]): x for x in c["gosia"]["rows"]}
    g1 = {(x["setting"], x["ring"]): x for x in c["gosia"]["second_line"] if x["order"] == "first"}
    gm = {(x["setting"], x["ring"]): x for x in c["gosia"]["second_line"] if x["order"] == "all"}
    cl = {x["rings"]: x for x in c["closure"]["rows"]}
    pa = {x["ring"]: x for x in c["paper"]["rows"]}
    p1 = c["paper"]["second_line"]
    pa1 = {x["ring"]: x for x in p1["rows"]}
    se = {(x["what"], x["ring"]): x for x in c["sensitivities"]["rows"]}
    sp = {x["rings"]: x for x in c["stored"]["rows"] if x["what"] != "B(E2)"}
    sb = {x["rings"]: x for x in c["stored"]["rows"] if x["what"] == "B(E2)"}

    def ok(name):
        bad = [x for x in c[name]["rows"] if not x["passed"]]
        return "pass" if not bad else ("MISS, known" if all(x.get("known_miss") for x in bad) else "FAIL")

    def pair(rows, key, setting=None):
        return "/".join(f"{100 * (rows[x] if setting is None else rows[setting, x])[key]:+.2f}" for x in RINGS)

    def three(rows, key, sign="+"):
        return "/".join(f"{100 * rows[x][key]:{sign}.2f}" for x in ROWS)

    return (
        f"{verdict(r)}. All orders (first order in brackets, not judged). "
        f"1 thickness ({ok('thickness')}): losses 58Ni {th['58Ni']['loss_mev']:.2f}, 60Ni {th['60Ni']['loss_mev']:.2f},"
        f" 62Ni {th['62Ni']['loss_mev']:.2f} MeV (42.7(8), 42.1(8), 42.4(8)); "
        f"2 GOSIA ({ok('gosia')}): ring 2/3 {pair(g, 'deviation', 'as_run')} % [{pair(g1, 'deviation', 'as_run')} %] "
        f"against GOSIA as run, E1 polarisation in both (1 %), of which {pair(g, 'explained', 'as_run')} % is GOSIA's "
        f"integer masses (its run with nuclear masses: {pair(gm, 'deviation', 'nuclear_masses')} %, not judged); "
        f"{pair(g, 'deviation', 'matched')} % [{pair(g1, 'deviation', 'matched')} %] with nuclear masses and the E1 "
        f"polarisation off in both (0.5 %), mostly GOSIA's non-relativistic kinematics (0.24 %); "
        f"3 closure ({ok('closure')}): mean of {len(c['closure']['seeds'])} seeds {three(cl, 'deviation')} % "
        f"[{three(cl, 'first_order')} %] for ring 2/ring 3/both (MC error {three(cl, 'monte_carlo', '')} %, scatter "
        f"error {three(cl, 'scatter', '')} %, bands {three(cl, 'band', '')} %); "
        f"4 paper ({ok('paper')}, the paper's ratios rebuilt through GOSIA): B(E2) "
        f"{pa['Ring2']['b_physim_e2b2']:.4f}/{pa['Ring3']['b_physim_e2b2']:.4f} "
        f"[{pa1['Ring2']['b_physim_e2b2']:.4f}/{pa1['Ring3']['b_physim_e2b2']:.4f}] e2b2 against "
        f"{pa['Ring2']['b_paper_e2b2']:.4f}/{pa['Ring3']['b_paper_e2b2']:.4f} ({100 * pa['Ring2']['deviation']:+.1f}/"
        f"{100 * pa['Ring3']['deviation']:+.1f} %, band 4.1 %; pulls {pa['Ring2']['pull_stat']:+.2f}/"
        f"{pa['Ring3']['pull_stat']:+.2f} stat), both {c['paper']['combined_b_e2b2']:.4f} "
        f"[{p1['combined_b_e2b2']:.4f}] against 0.0906(37); "
        f"5 sensitivities ({ok('sensitivities')}): Q +0.1 eb {se['q_plus', 'Ring2']['physim_percent']:+.2f}/"
        f"{se['q_plus', 'Ring3']['physim_percent']:+.2f} % (GOSIA {se['q_plus', 'Ring2']['gosia_percent']:+.2f}/"
        f"{se['q_plus', 'Ring3']['gosia_percent']:+.2f}), Q −0.1 eb {se['q_minus', 'Ring2']['physim_percent']:+.2f}/"
        f"{se['q_minus', 'Ring3']['physim_percent']:+.2f} % (GOSIA {se['q_minus', 'Ring2']['gosia_percent']:+.2f}/"
        f"{se['q_minus', 'Ring3']['gosia_percent']:+.2f}), feeding {se['feeding', 'Ring2']['physim_percent']:+.2f}/"
        f"{se['feeding', 'Ring3']['physim_percent']:+.2f} % (GOSIA {se['feeding', 'Ring2']['gosia_percent']:+.2f}/"
        f"{se['feeding', 'Ring3']['gosia_percent']:+.2f}); "
        f"6 stored run ({ok('stored')}): excitations per particle {pair(sp, 'deviation')} % from the rates, B(E2) "
        f"{three(sb, 'deviation')} % [{three(sb, 'first_order')} %], errors {three(sb, 'error', '')} %")


def table(r: dict) -> str:
    """A plain-text table, one line per compared number. "MISS" is outside its band by the known amount; a line
    marked "1st order" is the second line, not judged."""
    c = r["checks"]
    lines = [f"{'check':14s} {'what':46s} {'physim':>12s} {'reference':>12s} {'deviation':>11s} {'band':>10s}  result"]

    def line(check, what, ours, theirs, dev, band, result):
        if not isinstance(result, str):
            result = "pass" if result else "FAIL"
        lines.append(f"{check:14s} {what:46s} {ours:>12s} {theirs:>12s} {dev:>11s} {band:>10s}  {result}")

    second = "(1st order)"
    for x in c["thickness"]["rows"]:
        line("1 thickness", f"{x['ion']} loss at {x['entry_mev']:.1f} MeV", f"{x['loss_mev']:.2f} MeV",
             f"{x['reference_mev']:.1f} MeV", f"{x['loss_mev'] - x['reference_mev']:+.2f} MeV",
             f"+-{x['band_mev']:.1f} MeV", x["passed"])
    names = {"matched": "nuclear masses, DIP off", "as_run": "as run (A=62/12, DIP 3.8)",
             "nuclear_masses": "nuclear masses, DIP 3.8"}
    for x in c["gosia"]["rows"]:
        line("2 GOSIA", f"{x['ring']} R_el vs {names[x['setting']]}", f"{x['physim']:.4e}", f"{x['gosia']:.4e}",
             f"{100 * x['deviation']:+.2f} %", f"+-{100 * x['band']:.1f} %",
             "MISS (known, explained)" if x["known_miss"] else x["passed"])
    for x in c["gosia"]["second_line"]:
        line("2 " + (second if x["order"] == "first" else "(reported)"), f"{x['ring']} R_el vs {names[x['setting']]}",
             f"{x['physim']:.4e}", f"{x['gosia']:.4e}", f"{100 * x['deviation']:+.2f} %", "(reported)",
             "inside the band" if x["inside"] else "outside the band")
    n = len(c["closure"]["seeds"])
    for x in c["closure"]["rows"]:
        line("3 closure", f"B(E2) {x['rings']}, {n} seeds (MC {100 * x['monte_carlo']:.2f}, "
             f"scatter {100 * x['scatter']:.2f} %)", f"{x['b_e2fm4']:.1f}", f"{x['b_in_e2fm4']:.1f}",
             f"{100 * x['deviation']:+.2f} %", f"+-{100 * x['band']:.2f} %", x["passed"])
    for x in c["closure"]["rows"]:
        line("3 " + second, f"B(E2) {x['rings']}, the same events", f"{(1 + x['first_order']) * x['b_in_e2fm4']:.1f}",
             f"{x['b_in_e2fm4']:.1f}", f"{100 * x['first_order']:+.2f} %", "(reported)", "")
    for label, part in (("4 paper", c["paper"]), ("4 " + second, c["paper"]["second_line"])):
        judged = label == "4 paper"
        for x in part["rows"]:
            line(label, f"{x['ring']} B(E2) e2b2 (pull {x['pull_stat']:+.2f})", f"{x['b_physim_e2b2']:.4f}",
                 f"{x['b_paper_e2b2']:.4f}", f"{100 * x['deviation']:+.2f} %",
                 f"+-{100 * x['band']:.1f} %" if judged else "(reported)", x["passed"] if judged else "")
        line(label, "both rings (stat-weighted) e2b2", f"{part['combined_b_e2b2']:.4f}", f"{part['final_b_e2b2']:.4f}",
             f"{100 * part['combined_deviation']:+.2f} %", "(reported)", "")
    labels = {"q_plus": "Q_s = +0.1 eb", "q_minus": "Q_s = -0.1 eb", "feeding": "feeding (ENSDF states)"}
    for x in c["sensitivities"]["rows"]:
        line("5 sensitivity", f"{x['ring']} {labels[x['what']]}", f"{x['physim_percent']:+.2f} %",
             f"{x['gosia_percent']:+.2f} %", f"{x['difference_pp']:+.2f} pp", f"+-{x['band_pp']:.1f} pp",
             x["passed"])
        line("5 sensitivity", "  relative to GOSIA's change", "", "", f"{100 * x['difference_rel']:+.2f} %",
             f"+-{100 * x['band_rel']:.0f} %", x["passed"])
    for x in c["stored"]["rows"]:
        number = "{:.4e}" if x["what"] != "B(E2)" else "{:.1f}"
        line("6 stored run", f"{x['what']} {x['rings']} (+- {100 * x['error']:.2f} %)", number.format(x["physim"]),
             number.format(x["reference"]), f"{100 * x['deviation']:+.2f} %", f"+-{100 * x['band']:.2f} %", x["passed"])
    for x in (x for x in c["stored"]["rows"] if x["what"] == "B(E2)"):
        line("6 " + second, f"B(E2) {x['rings']}, the same run", f"{(1 + x['first_order']) * x['reference']:.1f}",
             f"{x['reference']:.1f}", f"{100 * x['first_order']:+.2f} %", "(reported)", "")
    if r["known_misses"]:  # one the reference file records, while its size stays
        lines.append("")
        lines.append("The known miss: " + c["gosia"]["why_missed"])
    return "\n".join(line.rstrip() for line in lines)
