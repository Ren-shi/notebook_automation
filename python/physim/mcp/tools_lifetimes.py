"""The lifetime tools: plan a lifetime measurement (plunger RDDS or DSAM), set up the plunger, predict its spectra,
and analyse the runs taken at its distances.

Each tool is a method of :class:`LifetimeTools`; its docstring is what the model reads. The tools are thin: the
planner works everything out (:meth:`~physim.nuclear.planner.Planner.lifetimes`,
:meth:`~physim.nuclear.planner.Planner.lifetime_spectra`, :meth:`~physim.nuclear.planner.Planner.analyse_lifetimes`
and the draft methods ``set_lifetime``, ``set_stopper``, ``set_degrader2``), and these turn its numbers into text
with units, as the other tools do ("8.63 ps", "4.38 % of c", "25 um"), leaving out the spectra's bins.
"""

from __future__ import annotations

import math
from typing import Any, Optional, Union

from ..nuclear.experiment import Experiment
from ..nuclear.planner import Planner
from .session import Session, ToolProblem
from .shape import plain
from .tools_setup import _changed

#: A quantity as the tools take it: text with a unit, or a number in the unit the argument names.
Amount = Union[str, int, float]


def _num(x: Optional[float], digits: int = 4) -> str:
    if x is None or not math.isfinite(x):
        return "—"
    return f"{x:.{digits}g}"


def _ps(x: Optional[float], digits: int = 4) -> Optional[str]:
    return None if x is None or not math.isfinite(x) else f"{_num(x, digits)} ps"


def _time(seconds: Optional[float]) -> Optional[str]:
    """A beam time as text: s, min or h."""
    if seconds is None or not math.isfinite(seconds):
        return None
    if seconds < 120:
        return f"{seconds:.3g} s"
    if seconds < 7200:
        return f"{seconds / 60:.3g} min"
    return f"{seconds / 3600:.3g} h"


def _um(values) -> list:
    return [f"{v:g} um" for v in values]


def level_plan(row: dict) -> dict:
    """One level of :meth:`Planner.lifetimes` as text with units."""
    out: dict[str, Any] = {
        "level": row["level"], "energy": f"{row['energy_kev']:g} keV",
        "line": None if row["gamma_kev"] is None else f"{row['gamma_kev']:g} keV to {row['final']}",
        "tau": _ps(row["tau_ps"]), "tau_source": row["tau_source"], "method": row["method"],
        "recoil_beta": f"{100 * row['beta']:.3g} % of c",
    }
    if row["recommended_distances_um"]:
        out["recommended_distances"] = _um(row["recommended_distances_um"])
    if row["sensitive_um"]:
        lo, hi = row["sensitive_um"]
        out["sensitive_range"] = f"{lo:.3g} to {hi:.3g} um"
    seps = {}
    for s in row["separation"]:  # one line per angle: the detectors at one angle share it
        key = f"{s['theta_deg']:g} deg"
        if key not in seps:
            seps[key] = (f"{s['separation_kev']:+.3g} keV against FWHM {s['fwhm_kev']:.3g} keV: "
                         + ("resolved" if s["resolved"] else "not resolved"))
    out["peak_separation"] = seps
    if row["attenuation"] is not None:
        out["attenuation_F"] = round(float(row["attenuation"]), 4)
    if row["rate_per_s"] is not None:
        out["line_rate"] = f"{row['rate_per_s']:.4g} /s in the γ-ray detectors (full energy)"
    if row["hours_per_distance"] is not None:
        out["beam_time_per_distance"] = _time(3600 * row["hours_per_distance"])
        out["counts_per_distance"] = round(float(row["counts_per_distance"]))
        out["at_distances"] = _um(row["distances_used_um"])
    if row.get("backing_advice"):
        out["backing_advice"] = backing_view(row["backing_advice"])
    if row["note"]:
        out["note"] = row["note"]
    return out


#: The backings listed per level (the setup's is added when it is not among them).
BACKINGS_SHOWN = 5


def backing_view(advice: list) -> list:
    """The DSAM backing advice of a level, best first, as text with units: the first :data:`BACKINGS_SHOWN` and
    the setup's backing material."""
    out = []
    for rank, c in enumerate(advice, start=1):
        if rank > BACKINGS_SHOWN and not c.get("setup"):
            continue
        um = (lambda x: f" ({x:.2g} um)" if x else "")
        row = {"rank": rank, "material": c["material"],
               "attenuation_F": round(float(c["attenuation"]), 3),
               "sensitivity_dF_dlntau": round(float(c["sensitivity"]), 4),
               "of_best": f"{100 * c['relative_sensitivity']:.0f} %",
               "slowing_down_time": _ps(c["slowing_time_ps"], 3), "time_to_stop": _ps(c["time_to_stop_ps"], 3),
               "thickness_to_stop": f"{c['thickness_to_stop_mg_cm2']:.3g} mg/cm2{um(c.get('thickness_to_stop_um'))}",
               "thickness_advised": f"{c['thickness_mg_cm2']:.2g} mg/cm2{um(c.get('thickness_um'))}"}
        if c.get("setup"):
            row["setup_backing"] = True
        out.append(row)
    return out


def foil_view(f: dict) -> dict:
    """The beam in one of the plunger's foils as text with units."""
    out: dict[str, Any] = {
        "foil": f"{f['thickness_mg_cm2']:.3g} mg/cm2 {f['material']} {f['foil']}",
        "beam": f"{f['beam_in_mev']:.4g} MeV in, " + ("stops in it" if f["beam_out_mev"] <= 0 else
                                                    f"{f['beam_out_mev']:.4g} MeV out"),
        "power": f"{f['power_w']:.3g} W ({f['deposited_mev']:.4g} MeV per particle); limit {f['power_limit_w']:g} W"}
    lines = []
    for x in f.get("lines") or []:
        line = f"{x['nuclide']} {x['gamma_kev']:.1f} keV ({x['transition']})"
        if x.get("rate_per_s"):
            line += f": {x['rate_per_s']:.3g} /s in the γ-ray detectors (singles)"
        lines.append(line)
    if lines:
        out["coulex_lines"] = lines[:8]
    if f.get("missing"):
        out["no_level_data"] = [m["nuclide"] for m in f["missing"]]
    return out


def plan_view(plan: dict) -> dict:
    """:meth:`Planner.lifetimes` as text with units."""
    out: dict[str, Any] = {"available": plan["available"], "reaction": plan["reaction"], "nucleus": plan["nucleus"],
                           "method": plan["method"], "precision_wanted": f"{100 * plan['precision']:g} % on tau"}
    if plan.get("reason"):
        out["reason"] = plan["reason"]
    if plan["distances_um"]:
        out["distances"] = _um(plan["distances_um"])
        if plan.get("offset_um"):
            out["offset"] = f"{plan['offset_um']:g} um"
    if plan.get("recoil"):
        r = plan["recoil"]
        out["recoil"] = (f"{r['nuclide']} made at {r['energy_made_mev']:.3g} MeV ({100 * plan['beta_made']:.3g} % of "
                         f"c), leaving the target at {r['energy_after_target_mev']:.3g} MeV "
                         f"({100 * plan['beta_after_target']:.3g} % of c)")
    if plan["time_to_stop_ps"]:
        out["time_to_stop"] = _ps(plan["time_to_stop_ps"], 3)
    out["levels"] = [level_plan(row) for row in plan["levels"]]
    if plan.get("stopper"):
        out["beam_in_plunger_foils"] = [foil_view(f) for f in plan["stopper"]]
    out["problems"] = list(plan["problems"])
    out["warnings"] = list(plan["warnings"])
    return plain(out)


def spectra_view(sp: dict) -> dict:
    """:meth:`Planner.lifetime_spectra` without the bins: per distance and angle where the peaks are, the line's
    counts and the unshifted share."""
    out: dict[str, Any] = {"available": sp["available"], "method": sp["method"]}
    if sp.get("reason"):
        out["reason"] = sp["reason"]
        return out
    out.update(level=sp["level"], line=f"{sp['gamma_kev']:g} keV", model=sp["model"],
               beam_time_per_distance=_time(sp["seconds_per_distance"]))
    edges = sp["edges_kev"]
    out["energy_range"] = f"{edges[0]:.1f} to {edges[-1]:.1f} keV in {len(edges) - 1} bins"
    rows = []
    for s in sp["spectra"]:
        row = {"angle": f"{s['theta_deg']:g} deg", "shifted_at": f"{s['shifted_kev']:.1f} keV",
               "unshifted_at": f"{s['unshifted_kev']:.1f} keV", "counts_in_line": round(float(sum(s["counts"])))}
        if s.get("distance_um") is not None:
            row = {"distance": f"{s['distance_um']:g} um", **row}
        if s.get("unshifted_share") is not None:
            row["unshifted_share"] = round(float(s["unshifted_share"]), 4)
        if s.get("attenuation") is not None:
            row["attenuation_F"] = round(float(s["attenuation"]), 4)
        rows.append(row)
    out["spectra"] = rows
    out["notes"] = list(sp.get("notes", []))
    return plain(out)


def analysis_view(a: dict) -> dict:
    """:meth:`Planner.analyse_lifetimes` as text with units."""
    out: dict[str, Any] = {"available": a["available"], "level": a["level"], "method": a["method"]}
    if not a["available"]:
        out["reason"] = a["reason"]
        return out
    out["tau"] = (f"{_num(a['tau_ps'])} ± {_num(a['stat_ps'], 2)} (statistics) ps; ± {_num(a['total_ps'], 2)} ps "
                  "with every source")
    out["budget"] = {k: _ps(v, 2) for k, v in a["budget"].items() if v}
    if a.get("beta"):
        out["recoil_beta"] = f"{100 * a['beta']:.3g} % of c (from the shifted peaks)" if a["points"] else \
            f"{100 * a['beta']:.3g} % of c"
    if a.get("rings_deg"):
        out["angles_used"] = [f"{t:g} deg" for t in a["rings_deg"]]
    if a["b_down"]:
        out[f"B({a.get('b_multipole', 'E2')} down)"] = (f"{_num(a['b_down'][0])} ± {_num(a['b_down'][1], 2)} "
                                                       f"{a['b_unit']}")
    out["points"] = [{"distance": f"{p['distance_um']:g} um", "unshifted_share": round(p["unshifted"], 4),
                      "error": round(p["error"], 4),
                      **({"tau": f"{_num(p['tau_ps'])} ± {_num(p['tau_err_ps'], 2)} ps"}
                         if p.get("tau_ps") is not None else {})} for p in a["points"]]
    out["compare"] = [{"source": c["source"], "tau": _ps(c["tau_ps"]),
                       **({"B down": f"{_num(c['b_down'])} {a['b_unit']}"} if c["b_down"] is not None else {})}
                      for c in a["compare"]]
    out["notes"] = list(a["notes"])
    return plain(out)


def _distances(values) -> Optional[list]:
    if values is None:
        return None
    if isinstance(values, (str, int, float)):
        values = [values]
    return [v if isinstance(v, str) else float(v) for v in values]


class LifetimeTools:
    """The lifetime tools, on one :class:`~physim.mcp.session.Session`."""

    def __init__(self, session: Session):
        self.session = session

    def plan_lifetime_measurement(self, level: Optional[Union[str, int]] = None, method: Optional[str] = None,
                                  precision: float = 0.1) -> dict:
        """How to measure the lifetime of a level of the nucleus the open experiment studies (Coulomb
        excitation's excited nucleus, or fusion-evaporation's studied residue). level is its name ("2+1", "2+"),
        its number in the level scheme or its energy ("652 keV"); without one every level with a strong line is
        listed, the studied one first. For each: tau and where it comes from (ENSDF half-life, the setup's B(E2),
        typed), the method that suits it (plunger "rdds" for about 1 ps to 1 ns, "dsam" for 10 fs to 1 ps), the
        recoil's speed, the recommended plunger distances and sensitive range, the shifted-unshifted peak
        separation at each γ-ray detector angle against its resolution, F(tau) for DSAM on the backing, the
        DSAM backing advice (the usual backings ranked by the sensitivity |dF/d ln tau|, with F(tau), the
        slowing-down time and the thickness to stop the recoil), and the beam time per distance for the wanted
        precision on tau (0.1 = 10%, statistics only). For a plunger, the beam in its foils: the power it deposits
        against what a thin foil takes, and the foils' own Coulomb-excitation lines. method "rdds" or
        "dsam" plans as if the setup used that method (nothing is changed; set_plunger changes it). Read the
        problems and warnings."""
        p = self.session.require()
        q = p
        note = None
        if method is not None:
            if method not in ("rdds", "dsam"):
                raise ToolProblem("method is \"rdds\" (a plunger) or \"dsam\" (Doppler-shift attenuation).")
            if p.experiment.lifetime is None or p.experiment.lifetime.method != method:
                q = Planner(Experiment.from_dict(p.experiment.to_dict()))
                if not q.set_lifetime(method=method):
                    out = plan_view(q.lifetimes(level, precision))
                    out["problems"] = list(q.problems) + out["problems"]
                    out["note"] = (f"Planned as if the method were {method}; that setup is not valid (see "
                                   "problems). Nothing was changed.")
                    return out
                note = f"Planned as if the method were {method}; the setup is unchanged (set_plunger changes it)."
        out = plan_view(q.lifetimes(level, precision))
        if note:
            out["note"] = note
        return out

    def set_plunger(self, distances: Optional[list] = None, offset: Optional[Amount] = None,
                    kind: Optional[str] = None, material: Optional[str] = None, thickness: Optional[Amount] = None,
                    density: Optional[Amount] = None, second_foil: Optional[dict] = None,
                    deorientation: Optional[dict] = None, method: Optional[str] = None,
                    remove: bool = False) -> dict:
        """Set up the lifetime measurement of the open experiment (its [lifetime] section; made if there is none).
        A plunger (method "rdds", the default for a new one): distances are the target-to-stopper distances (a list
        of text with units, ["10 um", "50 um"], or numbers in um); offset the distance at electrical contact ("8
        um"); the foil after the gap: kind "stopper" (stops the recoil), "degrader" (it flies through, slower) or
        "differential" (the first of two foils), its material ("Au", "Cu", "Ta"), thickness ("10 mg/cm2", or a
        length with a density); second_foil for a differential plunger: {"material": "Ta", "thickness": "6.5
        mg/cm2", "gap": "100 um"}; deorientation of the recoil's alignment in vacuum: {"alpha2": 0.18, "gamma2":
        "0.142 /ps", "alpha4": 0.24, "gamma4": "0.314 /ps"} or {"g": 0.47, "c2": "3.3 ps", "c4": "1.5 ps",
        "alpha2": ..., "alpha4": ...}, {} to remove it (the static model; the reaction's deorientation in vacuum,
        GOSIA's two-state model, does nothing beside a lifetime measurement). method "dsam" uses the target's backing
        instead (the plunger's fields are kept aside and come back with "rdds"). remove true takes the lifetime
        measurement out.
        Only what is given changes. Returns ok, the lifetime section and the summary with problems and warnings."""
        p = self.session.require()
        if remove:
            out = _changed(p, p.remove_lifetime())
            out["lifetime"] = None
            return out
        if method is not None and method not in ("rdds", "dsam"):
            raise ToolProblem("method is \"rdds\" (a plunger) or \"dsam\" (Doppler-shift attenuation).")
        foil = [kind, material, thickness]
        if any(x is not None for x in foil) and not all(x is not None for x in foil):
            current = (p.draft.get("lifetime") or {}).get("stopper") or {}
            kind = kind or current.get("kind")
            material = material or current.get("material")
            thickness = thickness or current.get("thickness")
            if None in (kind, material, thickness):
                raise ToolProblem("The plunger's foil needs kind, material and thickness (the first time).")
        if second_foil is not None and not {"material", "thickness", "gap"} <= set(second_foil):
            raise ToolProblem("second_foil needs material, thickness and gap, e.g. {\"material\": \"Ta\", "
                              "\"thickness\": \"6.5 mg/cm2\", \"gap\": \"100 um\"}.")
        ok = p.set_lifetime(method=method, distances=_distances(distances), offset=offset,
                            deorientation=deorientation)
        if kind is not None:
            ok = p.set_stopper(kind, material, thickness, density)
        if second_foil is not None:
            ok = p.set_degrader2(second_foil["material"], second_foil["thickness"], second_foil["gap"],
                                 second_foil.get("density"))
        out = _changed(p, ok)
        out["lifetime"] = plain(p.draft.get("lifetime"))
        return out

    def lifetime_spectra(self, distance: Optional[Amount] = None, level: Optional[Union[str, int]] = None) -> dict:
        """The predicted γ-ray spectra of a level's line (the studied one by default) for the open experiment's
        lifetime measurement: at each plunger distance (or only at distance, "25 um" or a number in um) and each
        γ-ray detector angle, where the shifted and unshifted peaks are, the counts in the line in the beam time
        per distance and the unshifted share; for DSAM the line shape's centroid and F at each angle. The bins are
        left out (the app draws them). Says how the spectra were made (model)."""
        p = self.session.require()
        d = None
        if distance is not None:
            if isinstance(distance, str):
                from ..nuclear.quantity import Quantity

                d = float(Quantity.parse(distance).to("um"))
            else:
                d = float(distance)
        return spectra_view(p.lifetime_spectra(d, level))

    def analyse_lifetimes(self, level: Optional[Union[str, int]] = None, method: Optional[str] = None) -> dict:
        """Analyse the open experiment's lifetime runs for a level (the studied one by default) as an
        experimenter would: a plunger from its distance series (one run per distance, taken with
        start_run-like tools or the app), the shifted and unshifted peak areas at the angles away from 90°, then
        method "ddcm" (default; the differential decay-curve method with the feeding line), "ddcm_three_foil"
        (default for a differential plunger) or "decay_curve"; DSAM from the line shapes of its runs ("lineshape").
        Returns tau with its uncertainty budget, the points per distance, B(E2 down) (or the line's multipole) and
        the comparison with ENSDF, the setup's Coulex B(E2), shell-model values and the level scheme's own."""
        p = self.session.require()
        if method is not None and method not in ("ddcm", "decay_curve", "ddcm_three_foil", "lineshape"):
            raise ToolProblem("method is \"ddcm\", \"decay_curve\", \"ddcm_three_foil\" or \"lineshape\".")
        return analysis_view(p.analyse_lifetimes(level, method=method))

    #: The tools in the order they are listed to the model, and whether each writes files.
    TOOLS: tuple = (
        ("plan_lifetime_measurement", False), ("set_plunger", True), ("lifetime_spectra", False),
        ("analyse_lifetimes", False),
    )
