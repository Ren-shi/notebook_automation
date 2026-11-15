"""A published fusion-evaporation experiment reproduced end to end through the planner (backlog item 75, D10).

**The experiment.** A. R. Usman, M. U. Khandaker, H. Haba, N. Otuka and M. Murakami, Nucl. Instrum. Methods B 469,
42 (2020), Table 5, EXFOR E2650: ¹⁶⁵Ho(α,xn)¹⁶⁸⁻ˣTm, x = 1–4, by stacked-foil activation at the RIKEN AVF cyclotron.
A 50.4 MeV α beam of 5.45 × 10¹¹ α/s (87.3 pnA) ran for 2.0 h through 26 metallic Ho foils (12.29 µm, 10.81 mg/cm²)
interleaved with Cu (8.70 mg/cm²) and Ti (4.69 mg/cm²) degraders and monitors; the Tm residues were counted off
line by their decay γ rays. The measured values are ``tests/reference/nuclear/evaporation_4He_165Ho_literature.csv``
(``ref = Usman2020``); Kormazeva et al. 2022 and Gadkari et al. 1997 in the same file are the cross-checks.

**The path, as a user takes it** (:func:`run`)::

    from physim.nuclear import reproduction

    r = reproduction.run()            # every Ho foil of the stack, and a few energies beyond it
    r.rows                            # per foil and channel: σ filled, residues made, σ inferred, measured
    c = reproduction.compare(r.rows, reproduction.measured())
    c["passed"], c["peaks"], c["strongest"], c["magnitudes"], c["ratio_3n_4n"]

For each foil: the example setup ``fusion_tm165`` (⁴He on 10.81 mg/cm² Ho, 87.3 pnA, 2.0 h) opened in the
:class:`~physim.nuclear.planner.Planner`, the beam energy set so that the mean energy in the foil is the paper's
(:func:`entry_energy`; activation papers tabulate the mean energy in each foil), the channel table filled from
physim's own fusion and evaporation (:meth:`Planner.fill_channels`, the app's "Fill from physim" button, 20 000
evaporation chains), the residues simulated (:func:`physim.nuclear.residues.simulate_residues`), and the production
cross section of each residue inferred from the residues made, as the experimenters did: σ = N / (Φ n t) with the
paper's beam intensity Φ, the foil's atoms n and the irradiation time t (not the setup's own luminosity, so the
setup's beam and target bookkeeping is part of what is checked). The statistical error of σ is physim's: the
binomial error of the channel's share of the evaporation chains and the residue simulation's own, added in
quadrature. (The counting statistics of the residues themselves, ~10⁸ per foil, are negligible.) The off-line decay
counting is outside the planner, so the comparison is at the level of production cross sections.

**The comparisons**, with the tolerances of backlog 75 as the request for this reproduction stated them, fixed
before the comparison was run (:func:`compare`). All are judged only where E_cm ≥ 1.1 V_B (physim's
Bass barrier): nearer the barrier physim's fusion misses the coupled-channels enhancement (``docs/theory/fusion.md``).

1. *Peak energies*: the peak of 2n, 3n and 4n, each the vertex of the parabola through the highest point and its
   two neighbours (physim on its own points, which go beyond the stack to 59 MeV), within the measurement's energy
   spread at its peak (EN-RSL-HW) plus 1 MeV. 1n is reported: its peak (19.6 MeV) is set by the fusion barrier.
2. *Strongest channel*: physim's strongest channel the same as the measured one at each foil, where the
   measurement resolves it (the two strongest measured channels differ by more than their combined error) and
   knows it (below the measured 4n peak, since (α,5n) was not measured).
3. *Magnitudes*: the strongest three channels that exceed 10% of the measured Σσ(1n–4n) at the foil, each within
   30% of the band the three measurements span (each ± its error, over ± its stated energy uncertainty; the file's
   advice: "keep the spread as the measurement's uncertainty"), where the foil is at or below the channel's
   measured peak. Above its peak a channel's tail is fed by pre-equilibrium emission of a fast first neutron, which
   physim (like PACE4) does not model: reported with that limit beside it, not judged.
4. *σ(3n)/σ(4n) at 45–50 MeV*: physim against the three measurements, reported (the most sensitive number to
   pre-equilibrium; both statistical codes are 3–5× low).

**The stack** (:func:`stack_energies`). The paper's foil order is not in EXFOR, nor in the RIKEN progress report
(RIKEN Accel. Prog. Rep. 53, 184), and the paper itself was not at hand. :data:`STACK` is the order inferred from
the tabulated foil energies with physim's stopping: between successive Ho foils the energy drops by one Ho foil
alone, or by one Ho foil and one Cu or one Ti foil, and each gap fits one of these to about 0.1 MeV, which gives
Ti, Ho, Cu, Ho, Ti, Ho, Ho, Cu, Ho, Ho, Ti, … (pairs of Ho foils between alternating Cu and Ti; the Ti in front is
the 50 MeV monitor). :func:`stack_energies` follows the beam from 50.4 MeV through that stack with physim's stopping.
Because the order was read off the same energies, this checks the consistency of physim's stopping with the
authors' (SRIM), not an independent foil order.
"""

from __future__ import annotations

import contextlib
import dataclasses
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

#: The example setup each foil starts from.
EXAMPLE = "fusion_tm165"
#: The paper (EXFOR E2650) and its measured values in the reference folder.
REFERENCE_FILE = "evaporation_4He_165Ho_literature.csv"
REF = "Usman2020"
CROSS_CHECKS = ("Kormazeva2022", "Gadkari1997")
#: The saved result of :func:`run` (written by ``scripts/reproduce_usman2020.py``).
RESULT_FILE = "reproduction_usman2020.csv"
#: PACE4 at the same energies (``scripts/pace4_reference.py --only excitation``), for the figure.
PACE4_FILE = "pace4_excitation_4He_165Ho.csv"

#: The beam: α particles per second (EXFOR INC-SOURCE), irradiation time (s), energy extracted (MeV; RIKEN Accel.
#: Prog. Rep. 53, 184).
BEAM_PER_S = 5.45e11
IRRADIATION_S = 2.0 * 3600.0
BEAM_MEV = 50.4
#: The foils: areal densities, mg/cm² (EXFOR SAMPLE and METHOD).
FOIL_MG_CM2 = {"Ho": 10.81, "Cu": 8.70, "Ti": 4.69}

#: Evaporation chains per foil (the app's default for "Fill from physim") and residue events per foil.
EVENTS = 20_000
RESIDUE_EVENTS = 20_000
SEED = 1
#: Beam energies beyond the stack (MeV, mean in the foil) where physim is also run, to find its peaks.
BEYOND = (51.0, 53.0, 55.0, 57.0, 59.0)

#: The xn channels the paper measured, and the residue of each.
MEASURED = {"1n": "168Tm", "2n": "167Tm", "3n": "166Tm", "4n": "165Tm"}
#: The measured channel each Tm residue stands for (physim labels one neutron "n").
CHANNEL_OF = {residue: channel for channel, residue in MEASURED.items()}
#: Channels whose peak, strongest-channel rank and magnitude are judged (1n peaks at the fusion barrier).
JUDGED = ("2n", "3n", "4n")
#: Judged only where E_cm ≥ ABOVE_BARRIER × V_B (as the other fusion-evaporation checks).
ABOVE_BARRIER = 1.1
#: Peak energies: within the measurement's energy spread at its peak plus this, MeV.
PEAK_MARGIN_MEV = 1.0
#: Magnitudes: channels above this share of the measured Σσ(1n–4n), the strongest this many, within this.
SHARE_CUT = 0.10
STRONGEST = 3
MAGNITUDE_TOLERANCE = 0.30
#: The energy range of the σ(3n)/σ(4n) comparison, MeV.
RATIO_RANGE = (45.0, 50.5)

def _stack() -> list:
    """Ti, then Ho foils with the degraders between them: after the 1st Ho a Cu, after the 2nd a Ti, then pairs of
    Ho foils with a Cu or a Ti (alternating, Cu first) after each pair; 26 Ho foils."""
    gaps = ["Cu", "Ti"] + ["" if k % 2 == 0 else ("Cu" if (k // 2) % 2 == 0 else "Ti") for k in range(23)]
    layers = ["Ti"]
    for g in gaps:
        layers.append("Ho")
        if g:
            layers.append(g)
    layers.append("Ho")
    return layers


#: The stack inferred from the foil energies (see the module docstring): the layers in beam order, "Ho" the 26
#: targets.
STACK = tuple(_stack())


# ---------------------------------------------------------------------------------------------------------------
# Energies


def _stopping(material: str):
    from .stopping import Stopping

    return Stopping("4He", material)


def mean_energy(e_in: float, material: str = "Ho", mg_cm2: Optional[float] = None, points: int = 41) -> float:
    """The mean beam energy through a foil (depth average, as the fusion's thick-target average), MeV."""
    s = _stopping(material)
    t = FOIL_MG_CM2[material] if mg_cm2 is None else mg_cm2
    xs = np.linspace(0.0, t, points)
    return float(np.mean([e_in] + [float(s.energy_after(e_in, float(x))) for x in xs[1:]]))


def entry_energy(e_mean: float, material: str = "Ho") -> float:
    """The beam energy entering a foil whose mean energy is ``e_mean`` (MeV)."""
    lo, hi = e_mean, e_mean + 3.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if mean_energy(mid, material) < e_mean else (lo, mid)
    return 0.5 * (lo + hi)


def stack_energies(beam_mev: float = BEAM_MEV, stack=STACK) -> list:
    """Follow the beam through :data:`STACK` with physim's stopping: per Ho foil (in beam order) the energy in, the
    mean energy and the energy out, MeV."""
    out, e = [], float(beam_mev)
    for layer in stack:
        s = _stopping(layer)
        after = float(s.energy_after(e, FOIL_MG_CM2[layer]))
        if layer == "Ho":
            out.append({"foil": len(out) + 1, "e_in_mev": e, "e_mean_mev": mean_energy(e), "e_out_mev": after})
        e = after
    return out


# ---------------------------------------------------------------------------------------------------------------
# The reproduction


@contextlib.contextmanager
def _cached_evaporation(cache: dict):
    """Within the block, :func:`physim.nuclear.evaporation.evaporate` remembers its results by its inputs (the
    compound nucleus, E*, the shape of σ_ℓ, events, seed). The evaporation is the slow part and does not depend on
    σ_fus's normalisation, so the validation suite reruns the reproduction cheaply, and a broken fusion still shows."""
    from . import evaporation

    original = evaporation.evaporate

    def cached(compound, excitation_mev, *, sigma_l_mb=None, events=10000, seed=None, **kw):
        if kw or sigma_l_mb is None or not isinstance(seed, (int, np.integer)):
            return original(compound, excitation_mev, sigma_l_mb=sigma_l_mb, events=events, seed=seed, **kw)
        w = np.asarray(sigma_l_mb, dtype=float)
        key = (str(compound), round(float(excitation_mev), 9), tuple(np.round(w / w.sum(), 12)), int(events),
               int(seed))
        if key not in cache:
            cache[key] = original(compound, excitation_mev, sigma_l_mb=sigma_l_mb, events=events, seed=seed)
        # The chains depend on σ_ℓ's shape only (the spins are drawn from σ_ℓ / Σσ_ℓ); σ_fus is its sum.
        return dataclasses.replace(cache[key], sigma_fusion_mb=float(w.sum()))

    evaporation.evaporate = cached
    try:
        yield
    finally:
        evaporation.evaporate = original


_CACHE: dict = {}


def atoms_per_cm2() -> float:
    """¹⁶⁵Ho atoms per cm² in one foil (10.81 mg/cm², Avogadro's number, AME2020 atomic mass)."""
    from . import data

    return FOIL_MG_CM2["Ho"] * 1e-3 * 6.02214076e23 / data.nuclide("165Ho").atomic_mass_u


def foil(e_mean: float, events: int = EVENTS, residue_events: int = RESIDUE_EVENTS, seed: int = SEED,
         where: str = "stack") -> list:
    """One foil through the planner (module docstring): one row per channel of the filled table."""
    from .planner import Planner
    from .residues import simulate_residues

    p = Planner.example(EXAMPLE)
    e_in = entry_energy(e_mean)
    if not p.set("beam", "energy", f"{e_in:.6f} MeV"):
        raise ValueError(f"the setup at {e_in:.4g} MeV is invalid: {p.problems}")
    if not p.fill_channels(events, seed):
        raise ValueError(f"filling the channels failed: {p.problems}")
    made = p.physim_channels(events, seed)
    # What the app gives back: the table as the setup now holds it, and the residues simulated from it.
    ev = simulate_residues(p.experiment, events=residue_events, seed=seed)
    norm = BEAM_PER_S * atoms_per_cm2() * IRRADIATION_S * 1e-27   # residues per mb
    rows = []
    for ch in ev.table.channels:
        share = ch.cross_section_mb / made.sigma_fusion_mb
        share_err = ch.cross_section_mb * math.sqrt(max(1.0 - share, 0.0) / (share * events)) if share > 0 else 0.0
        rate, rate_err = ev.rate(ch.label)
        n, n_err = rate * ev.beam_time_s, rate_err * ev.beam_time_s
        sigma = n / norm
        rows.append({"e_lab_mev": float(e_mean), "e_in_mev": e_in, "where": where,
                     "channel": CHANNEL_OF.get(ch.residue, ch.label), "residue": ch.residue,
                     "sigma_table_mb": ch.cross_section_mb, "share_unc_mb": share_err,
                     "residues_made": n, "residues_made_unc": n_err, "sigma_mb": sigma,
                     "sigma_unc_mb": math.hypot(share_err, sigma * n_err / n if n else 0.0),
                     "sigma_fus_mb": made.sigma_fusion_mb, "excitation_mev": made.excitation_mev,
                     "mean_l": made.mean_l})
    return rows


@dataclass
class Reproduction:
    """The result of :func:`run`: one row per foil and channel (see :func:`foil`)."""

    rows: list
    energies: list
    events: int
    residue_events: int
    seed: int
    notes: list = field(default_factory=list)


def measured(reference=None) -> list:
    """The rows of the measured ¹⁶⁵Ho(α,xn) file (all three data sets)."""
    from .validation import read_reference

    return read_reference(_reference(reference) / REFERENCE_FILE)[1]


def foil_energies(rows: list) -> list:
    """The paper's foil energies (mean in each Ho foil), highest first."""
    return sorted({float(r["e_lab_mev"]) for r in rows if r["ref"] == REF}, reverse=True)


def run(events: int = EVENTS, residue_events: int = RESIDUE_EVENTS, seed: int = SEED, energies=None,
        beyond=BEYOND, reference=None, progress=None) -> Reproduction:
    """Every Ho foil of the paper (or ``energies``), then ``beyond``, through the planner."""
    if energies is None:
        energies = foil_energies(measured(reference))
    rows = []
    with _cached_evaporation(_CACHE):
        for where, es in (("stack", energies), ("beyond", beyond)):
            for e in es:
                rows += foil(float(e), events, residue_events, seed, where)
                if progress:
                    progress(e)
    return Reproduction(rows, [float(e) for e in energies], events, residue_events, seed)


# ---------------------------------------------------------------------------------------------------------------
# The comparisons


def _curves(rows: list) -> dict:
    """{(ref, channel): (E, σ, error, energy uncertainty)} sorted by energy."""
    out: dict = {}
    for r in rows:
        out.setdefault((r["ref"], r["channel"]), []).append(
            (float(r["e_lab_mev"]), float(r["sigma_mb"]), float(r["uncertainty_mb"]), float(r["e_lab_unc_mev"])))
    return {k: tuple(np.array(c) for c in zip(*sorted(v))) for k, v in out.items()}


def peak(e, s) -> tuple:
    """(energy, value) of a curve's peak: the vertex of the parabola through its highest point and the two
    neighbours, or the end point where the highest point is at an end (then flagged by energy = that end)."""
    e, s = np.asarray(e, dtype=float), np.asarray(s, dtype=float)
    i = int(np.argmax(s))
    if i == 0 or i == len(s) - 1:
        return float(e[i]), float(s[i]), True
    a, b, c = np.polyfit(e[i - 1:i + 2], s[i - 1:i + 2], 2)
    if a >= 0:
        return float(e[i]), float(s[i]), False
    v = -b / (2 * a)
    return float(v), float(np.polyval([a, b, c], v)), False


def _band(curves: dict, channel: str, e: float):
    """(low, high) of the measurements of one channel at energy e: each data set interpolated (log σ, linear in
    E) over e ± its energy uncertainty, ± its error; None where no set covers e."""
    lo, hi = math.inf, 0.0
    for (_, ch), (en, sig, unc, eu) in curves.items():
        if ch != channel or not en[0] <= e <= en[-1]:
            continue
        d = float(np.interp(e, en, eu))
        grid = np.clip(np.linspace(e - d, e + d, 21), en[0], en[-1])
        s = np.exp(np.interp(grid, en, np.log(sig)))
        rel = np.interp(grid, en, unc / sig)
        lo, hi = min(lo, float(np.min(s * (1 - rel)))), max(hi, float(np.max(s * (1 + rel))))
    return (lo, hi) if hi > 0 else None


def physim_curves(rows: list) -> dict:
    """{channel: (E, σ inferred, error)} of a reproduction, sorted by energy (stack and beyond)."""
    out: dict = {}
    for r in rows:
        out.setdefault(r["channel"], []).append((float(r["e_lab_mev"]), float(r["sigma_mb"]),
                                                 float(r["sigma_unc_mb"])))
    return {k: tuple(np.array(c) for c in zip(*sorted(v))) for k, v in out.items()}


def compare(rows: list, measured_rows: list) -> dict:
    """The four comparisons of the module docstring. ``rows`` from :func:`run` (or the saved result read back);
    returns their entries, whether each passes, the worst deviation as a fraction of its tolerance, and
    ``passed``."""
    from .fusion import Fusion

    f = Fusion("4He", "165Ho")
    vb = f.barrier().height_mev

    def above(e: float) -> bool:
        return float(f.e_cm(e)) / vb >= ABOVE_BARRIER

    curves = _curves(measured_rows)
    ours = physim_curves(rows)
    at: dict = {}  # foil energy -> {channel: σ inferred}
    for r in rows:
        if r["where"] == "stack":
            at.setdefault(float(r["e_lab_mev"]), {})[r["channel"]] = float(r["sigma_mb"])
    usman = {ch: curves[(REF, ch)] for ch in MEASURED if (REF, ch) in curves}
    meas_at: dict = {}
    for ch, (en, sig, unc, eu) in usman.items():
        for e, s, u, d in zip(en, sig, unc, eu):
            meas_at.setdefault(float(e), {})[ch] = (float(s), float(u), float(d))
    peak_e = {ch: float(en[int(np.argmax(sig))]) for ch, (en, sig, _, _) in usman.items()}

    # 1. Peak energies.
    peaks = []
    for ch in MEASURED:
        en, sig, _, eu = usman[ch]
        em, sm, edge_m = peak(en, sig)
        tol = float(eu[int(np.argmax(sig))]) + PEAK_MARGIN_MEV
        pe, ps, pu = ours.get(ch, (np.zeros(0),) * 3)
        if len(pe) >= 3:
            ep, sp, edge_p = peak(pe, ps)
        else:
            ep, sp, edge_p = math.nan, math.nan, True
        judged = ch in JUDGED and above(em) and not edge_m
        peaks.append({"channel": ch, "measured_mev": em, "measured_mb": sm, "physim_mev": ep, "physim_mb": sp,
                      "physim_at_edge": edge_p, "tolerance_mev": tol, "judged": judged,
                      "passed": (abs(ep - em) <= tol and not edge_p) if judged else None})

    # 2. Strongest channel.
    strongest = []
    for e in sorted(meas_at):
        here = meas_at[e]
        if e not in at:
            continue
        order = sorted(here, key=lambda c: -here[c][0])
        top = order[0]
        resolved = len(order) < 2 or (here[top][0] - here[order[1]][0]
                                      > math.hypot(here[top][1], here[order[1]][1]))
        nxt = f"{int(top[:-1]) + 1}n"
        known = nxt in usman or e <= peak_e[top]
        mine = max(at[e], key=at[e].get)
        judged = above(e) and resolved and known
        strongest.append({"e_lab_mev": e, "measured": top, "physim": mine, "resolved": resolved, "known": known,
                          "above_barrier": above(e), "judged": judged, "passed": (mine == top) if judged else None})

    # 3. Magnitudes.
    magnitudes = []
    for e in sorted(meas_at):
        if e not in at:
            continue
        here = meas_at[e]
        total = sum(v[0] for v in here.values())
        big = sorted((c for c in here if here[c][0] > SHARE_CUT * total), key=lambda c: -here[c][0])[:STRONGEST]
        for ch in sorted(here):
            s, u, _ = here[ch]
            mine = at[e].get(ch, 0.0)
            in_tail = e > peak_e[ch]
            judged = ch in big and ch in JUDGED and above(e) and not in_tail
            band = _band(curves, ch, e)
            if band is None:
                dev = math.nan
            else:
                lo, hi = band
                dev = mine / hi - 1 if mine > hi else (mine / lo - 1 if mine < lo else 0.0)
            why = ("judged" if judged else "below 1.1 V_B" if not above(e) else
                   "tail above the measured peak (pre-equilibrium)" if in_tail and ch in big else
                   "below 10% of the measured total" if ch not in big else "1n (not judged)")
            magnitudes.append({"e_lab_mev": e, "channel": ch, "physim_mb": mine, "measured_mb": s,
                               "measured_unc_mb": u, "band_mb": band, "deviation": dev, "ratio": mine / s,
                               "judged": judged, "why": why,
                               "passed": (abs(dev) <= MAGNITUDE_TOLERANCE) if judged else None})

    # 4. σ(3n)/σ(4n) at 45–50 MeV.
    ratios = []
    for name in (REF, *CROSS_CHECKS):
        if (name, "3n") not in curves or (name, "4n") not in curves:
            continue
        e3, s3, _, _ = curves[(name, "3n")]
        e4, s4, _, _ = curves[(name, "4n")]
        for e, s in zip(e4, s4):
            if RATIO_RANGE[0] <= e <= RATIO_RANGE[1] and e in e3:
                p3, p4 = ours["3n"], ours["4n"]
                q = float(np.interp(e, p3[0], p3[1]) / np.interp(e, p4[0], p4[1]))
                ratios.append({"ref": name, "e_lab_mev": float(e), "measured": float(s3[list(e3).index(e)] / s),
                               "physim": q})

    def worst(entries, value, tol):
        vals = [abs(value(x)) / tol for x in entries if x["judged"]]
        return max(vals) if vals else 0.0

    w_peak = worst(peaks, lambda x: (x["physim_mev"] - x["measured_mev"]) / x["tolerance_mev"], 1.0)
    w_mag = worst(magnitudes, lambda x: x["deviation"], MAGNITUDE_TOLERANCE)
    ok = {"peaks": all(x["passed"] for x in peaks if x["judged"]),
          "strongest": all(x["passed"] for x in strongest if x["judged"]),
          "magnitudes": all(x["passed"] for x in magnitudes if x["judged"])}
    # Reported, not judged: physim/Usman by channel where judged, in the tails, and for 1n; and the range of a
    # common scale factor on physim's cross sections that the magnitudes would still pass (the test's power).
    reported: dict = {}
    for x in magnitudes:
        if x["measured_mb"] <= 0 or not above(x["e_lab_mev"]):
            continue
        kind = ("judged" if x["judged"] else "1n" if x["channel"] == "1n" else
                "tail" if x["e_lab_mev"] > peak_e[x["channel"]] else "rising")
        reported.setdefault((x["channel"], kind), []).append((x["e_lab_mev"], x["ratio"]))
    judged_m = [x for x in magnitudes if x["judged"] and x["band_mb"]]
    scale = ((max(x["band_mb"][0] * (1 - MAGNITUDE_TOLERANCE) / x["physim_mb"] for x in judged_m),
              min(x["band_mb"][1] * (1 + MAGNITUDE_TOLERANCE) / x["physim_mb"] for x in judged_m))
             if judged_m else (math.nan, math.nan))
    return {"peaks": peaks, "strongest": strongest, "magnitudes": magnitudes, "ratio_3n_4n": ratios,
            "reported": reported, "scale_passing": scale,
            "ok": ok, "passed": all(ok.values()),
            "worst": {"peaks": w_peak, "magnitudes": w_mag,
                      "strongest": float(sum(1 for x in strongest if x["judged"] and not x["passed"]))},
            "barrier_mev": vb}


# ---------------------------------------------------------------------------------------------------------------
# Saving and reading the result

COLUMNS = ("e_lab_mev", "e_in_mev", "where", "channel", "residue", "sigma_table_mb", "share_unc_mb",
           "residues_made", "residues_made_unc", "sigma_mb", "sigma_unc_mb", "sigma_fus_mb", "excitation_mev",
           "mean_l", "measured_mb", "measured_unc_mb")


def result_lines(rep: Reproduction, measured_rows: list, produced_by: str, comparison: Optional[dict] = None) -> list:
    """The saved result: a provenance header and one CSV row per foil and channel."""
    import physim

    usman = {(float(r["e_lab_mev"]), r["channel"]): (r["sigma_mb"], r["uncertainty_mb"])
             for r in measured_rows if r["ref"] == REF}
    lines = [
        f"# tool: physim {physim.__version__}, the planner's own path for each Ho foil: Planner.example(\"{EXAMPLE}\"), "
        "beam energy set, fill_channels (physim's fusion and evaporation), simulate_residues; σ inferred from the "
        "residues made",
        f"# produced by: {produced_by}",
        "# reference: reproduction of A. R. Usman et al., Nucl. Instrum. Methods B 469, 42 (2020), Table 5 (EXFOR "
        "E2650): 165Ho(alpha,xn)Tm by stacked-foil activation, RIKEN AVF cyclotron",
        f"# settings: 4He on 165Ho, {FOIL_MG_CM2['Ho']} mg/cm2 Ho per foil, {BEAM_PER_S:.3g} alpha/s (87.3 pnA), "
        f"{IRRADIATION_S / 3600:g} h; {rep.events} evaporation chains and {rep.residue_events} residue events per "
        f"foil, seed {rep.seed}",
        "#   e_lab_mev is the paper's mean energy in the foil; e_in_mev the beam energy entering it in the setup (physim's "
        "stopping); where = stack (the paper's foils) or beyond (physim only, to find its peaks)",
        "#   sigma_table_mb is the channel table fill_channels wrote (thick-target fusion x the channel's share of the "
        "chains), share_unc_mb its binomial error; residues_made in the irradiation, from simulate_residues, with its "
        "statistical error",
        "#   sigma_mb = residues_made / (beam/s x Ho atoms/cm2 x time) with the paper's numbers; sigma_unc_mb physim's "
        "statistical error (share and residue simulation in quadrature); measured_mb, measured_unc_mb = Usman2020 "
        "where measured",
    ]
    if comparison is not None:
        lines.append("# comparison: " + summary(comparison))
    lines.append(",".join(COLUMNS))
    for r in rep.rows:
        m = usman.get((float(r["e_lab_mev"]), r["channel"]), ("", ""))
        vals = {**r, "measured_mb": m[0], "measured_unc_mb": m[1]}
        lines.append(",".join(_fmt(vals[c]) for c in COLUMNS))
    return lines


def _fmt(v) -> str:
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.6g}"
    return str(v)


def read_result(path) -> list:
    """The rows of a saved result, numbers as floats."""
    from .validation import read_reference

    rows = read_reference(path)[1]
    num = set(COLUMNS) - {"where", "channel", "residue", "measured_mb", "measured_unc_mb"}
    return [{k: (float(v) if k in num else v) for k, v in r.items()} for r in rows]


def reported(c: dict) -> str:
    """physim / Usman, as ranges: where judged, in the tails above each channel's measured peak (pre-equilibrium), on
    the rising edges below 10% of the total, and for 1n, all above 1.1 V_B; and the common scale factors on physim's
    cross sections that the magnitudes would still pass."""
    parts = []
    for kind, what in (("judged", "where judged"), ("tail", "tails above the measured peaks (pre-equilibrium)"),
                       ("rising", "rising edges under 10% of the total"),
                       ("1n", "1n (its peak is at the fusion barrier; physim has no 1n above 31 MeV)")):
        items = [f"{ch} {min(q for _, q in v):.2f}–{max(q for _, q in v):.2f} ({min(e for e, _ in v):g}–"
                 f"{max(e for e, _ in v):g} MeV)" for (ch, k), v in sorted(c["reported"].items()) if k == kind]
        if items:
            parts.append(f"{what}: " + ", ".join(items))
    lo, hi = c["scale_passing"]
    return ("physim/Usman " + "; ".join(parts)
            + f". A common factor on physim's cross sections between {lo:.2f} and {hi:.2f} would still pass the "
            "magnitudes (the three measurements span a wide band)")


def summary(c: dict) -> str:
    """One line per comparison, for headers and reports."""
    p = "; ".join(f"{x['channel']} {x['physim_mev']:.1f} vs {x['measured_mev']:.1f} MeV" for x in c["peaks"]
                  if x["judged"])
    s = [x for x in c["strongest"] if x["judged"]]
    m = [x for x in c["magnitudes"] if x["judged"]]
    worst = max(m, key=lambda x: abs(x["deviation"])) if m else None
    q = [x for x in c["ratio_3n_4n"]]
    return (f"peaks {'pass' if c['ok']['peaks'] else 'FAIL'} ({p}); strongest channel "
            f"{sum(x['passed'] for x in s)}/{len(s)} foils; magnitudes {'pass' if c['ok']['magnitudes'] else 'FAIL'} "
            f"({len(m)} points" + (f", worst {worst['deviation']:+.0%} outside the band, {worst['channel']} at "
                                   f"{worst['e_lab_mev']:g} MeV" if worst else "") + "); "
            "σ(3n)/σ(4n) at 45–50 MeV physim " + ", ".join(f"{x['physim']:.3f}" for x in q if x["ref"] == REF)
            + " vs measured " + ", ".join(f"{x['measured']:.2f} ({x['ref'][:-4]} {x['e_lab_mev']:g})" for x in q))


# ---------------------------------------------------------------------------------------------------------------
# The figure


def plot(rows: list, measured_rows: list, pace4: Optional[dict] = None, axes=None):
    """The four excitation functions: physim (inferred from its residues, through the planner), the three
    measurements, and PACE4 where given ({channel: (E, σ)}). Returns the figure."""
    from ..plot import _plt

    plt = _plt()
    if axes is None:
        fig, axes = plt.subplots(2, 2, figsize=(9, 6.5), sharex=True)
    else:
        fig = axes.flat[0].figure
    curves = _curves(measured_rows)
    ours = physim_curves(rows)
    marks = {REF: ("o", "#222222", "Usman 2020 (reproduced)"), "Kormazeva2022": ("s", "#1f77b4", "Kormazeva 2022"),
             "Gadkari1997": ("^", "#2ca02c", "Gadkari 1997")}
    for ax, ch in zip(axes.flat, ("1n", "2n", "3n", "4n")):
        for name, (m, col, label) in marks.items():
            if (name, ch) in curves:
                en, sig, unc, eu = curves[(name, ch)]
                ax.errorbar(en, sig, yerr=unc, xerr=eu, fmt=m, ms=4, color=col, mfc="none" if name != REF else col,
                            lw=0.7, label=label)
        if ch in ours:
            pe, ps, pu = ours[ch]
            ax.fill_between(pe, ps - pu, ps + pu, color="#d62728", alpha=0.25, lw=0)
            ax.plot(pe, ps, "-", color="#d62728", lw=1.6, label="physim (planner)")
        if pace4 and ch in pace4:
            ax.plot(*pace4[ch], "--", color="#ff7f0e", lw=1.2, label="PACE4 4.36.0")
        ax.set_yscale("log")
        ax.set_ylim(0.05, 3000)
        ax.set_title(f"¹⁶⁵Ho(α,{ch}){_SUPERSCRIPT[MEASURED[ch]]}", fontsize=9, loc="left")
        ax.grid(alpha=0.3, lw=0.5)
    for ax in axes[-1]:
        ax.set_xlabel("α energy in the foil (MeV, lab)")
    for ax in axes[:, 0]:
        ax.set_ylabel("σ (mb)")
    axes.flat[0].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    return fig


_SUPERSCRIPT = {"168Tm": "¹⁶⁸Tm", "167Tm": "¹⁶⁷Tm", "166Tm": "¹⁶⁶Tm", "165Tm": "¹⁶⁵Tm"}


def _reference(reference=None) -> Path:
    from .validation import reference_dir

    ref = reference_dir() if reference is None else Path(reference)
    if ref is None:
        raise FileNotFoundError("no reference folder (tests/reference/nuclear)")
    return ref


def saved(reference=None) -> list:
    """The saved result (``reproduction_usman2020.csv``), as :func:`read_result` reads it."""
    return read_result(_reference(reference) / RESULT_FILE)


def figure(reference=None):
    """:func:`plot` of the saved result, the measurements and PACE4."""
    ref = _reference(reference)
    pace4 = ref / PACE4_FILE
    return plot(saved(ref), measured(ref), read_pace4(pace4) if pace4.exists() else None)


def read_pace4(path) -> dict:
    """{channel: (E, σ)} for the xn channels from the PACE4 excitation-function file."""
    from .validation import read_reference

    rows = read_reference(path)[1]
    out: dict = {}
    for r in rows:
        z, a = int(r["Z"]), int(r["A"])
        if z == 69 and 165 <= a <= 168:
            out.setdefault(f"{169 - a}n", []).append((float(r["beam_energy_mev"]), float(r["sigma_mb"])))
    return {k: tuple(np.array(c) for c in zip(*sorted(v))) for k, v in out.items()}
