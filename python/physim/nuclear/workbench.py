"""The workbench's structure, without the web page (backlog item 64): the setup as a list of decisions, the status
strip, the checks, the templates, and new experiments.

The app draws what this module returns, so every line it shows can be had from Python too::

    from physim.nuclear.planner import Planner
    from physim.nuclear import workbench

    p = Planner(workbench.new_setup("16O", "64 MeV", "58Ni", measure="coulex-target"))
    for d in workbench.decisions(p):
        print(d.title, "·", d.state, "·", d.consequence)
    workbench.status_strip(p)          # experiment · beam · target · detectors · last run · warnings
    workbench.checks(p)                # warnings, overlaps, "setup changed since run 3"

**The decisions**, in the order they are made: Beam → Target and reaction → Particle detectors → γ-ray detectors
→ Run conditions. Each is a card with its current state in one line and its consequence in a second; its fields,
and the rarely touched ones behind a "more" fold.

**The stages** of the experiment are the app's tabs, left to right: :data:`STAGES`.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Optional

from . import guide

#: The app's tabs, left to right: the stages of an experiment. (key, label)
STAGES = (("setup", "Setup"), ("plan", "Plan"), ("run", "Run"), ("data", "Data"), ("analysis", "Analysis"),
          ("report", "Report"), ("physics", "Physics"))

# -- the fields of each decision: (field, label, placeholder) -------------------------------------------------------
BEAM_FIELDS = [("nuclide", "Nuclide", "4He"), ("energy", "Energy", "5.5 MeV or 4 MeV/u"),
               ("current", "Current", "1 pnA or 10 enA")]
BEAM_MORE = [("charge_state", "Charge state", "6"), ("energy_spread", "Energy spread (FWHM)", "0.1 %"),
             ("spot_size", "Spot size (FWHM)", "2 mm")]
TARGET_FIELDS = [("material", "Material", "Au, 208Pb, CD2"), ("thickness", "Thickness", "0.5 mg/cm2 or 1 um")]
TARGET_MORE = [("position", "Position along the beam", "0 mm"), ("tilt", "Tilt", "0 deg"),
               ("density", "Density", "19.3 g/cm3")]
BACKING_FIELDS = [("material", "Backing material", "C"), ("thickness", "Backing thickness", "20 ug/cm2")]
REACTION_FIELDS = [("energy", "State energy", "1.454 MeV"), ("b_up", "B(Eλ↑)", "0.0695 e2b2 or 695 e2fm4")]
REACTION_TYPES = {"elastic": "Elastic (Rutherford) scattering", "coulex": "Coulomb excitation"}
DETECTOR_FIELDS = [("name", "Name", "D1"), ("theta", "θ", "45 deg"), ("phi", "φ", "0 deg"),
                   ("distance", "Distance", "100 mm"), ("width", "Width", "50 mm"), ("height", "Height", "50 mm"),
                   ("strips_x", "Strips x", "16"), ("strips_y", "Strips y", "16"), ("radius", "Radius", "5 mm"),
                   ("inner_radius", "Inner radius", "9 mm"), ("outer_radius", "Outer radius", "41 mm"),
                   ("rings", "Rings", "16"), ("sectors", "Sectors", "24"), ("thickness", "Thickness", "300 um")]
DETECTOR_MORE = [("dead_layer", "Dead layer", "0.5 um"), ("resolution", "Resolution (FWHM)", "20 keV"),
                 ("threshold", "Threshold", "200 keV"), ("material", "Material", "Si")]
GAMMA_FIELDS = [("name", "Name", "Ge1"), ("theta", "θ", "90 deg"), ("phi", "φ", "90 deg"),
                ("distance", "Distance", "120 mm"), ("radius", "Crystal radius", "35 mm")]
GAMMA_MORE = [("resolution", "Resolution (FWHM)", "2.5 keV"), ("efficiency", "Efficiency (full peak)", "2 %"),
              ("material", "Crystal material", "Ge"), ("absorbers", "Absorbers", "Pb 1 mm, Cu 0.5 mm")]
#: For a γ-ray detector with real crystals (a model of the catalogue), in place of the radius; all behind "more".
CRYSTAL_FIELDS = [("crystals", "Crystals", "4"), ("crystal_diameter", "Crystal diameter", "50 mm"),
                  ("crystal_length", "Crystal length", "70 mm"), ("crystal_pitch", "Crystal pitch", "45 mm"),
                  ("housing_side", "Housing side", "101 mm"), ("window_gap", "Window to crystal", "5 mm")]
CLOVER_FIELDS = [("addback_factor", "Add-back factor at 1332 keV", "1.5"),
                 ("shield_thickness", "Shield thickness", "25 mm"),
                 ("suppression_factor", "Suppression factor at 1332 keV", "3")]
RUN_FIELDS = [("beam_time", "Beam time", "12 h"), ("counts_wanted", "Counts wanted", "5000")]
RUN_MORE = [("coincidence_window", "Coincidence window", "100 ns"), ("dead_time", "Dead time", "5 us"),
            ("room_background", "Room background per crystal", "1 /s"),
            ("shaping_time", "Shaping time (pile-up)", "3 us")]
#: Which size fields each particle-detector shape uses.
SHAPE_FIELDS = {"rectangle": {"width", "height", "strips_x", "strips_y"}, "circle": {"radius"},
                "annular": {"inner_radius", "outer_radius", "rings", "sectors"}}

#: The decision groups, in order.
GROUPS = ("Beam", "Target and reaction", "Particle detectors", "γ-ray detectors", "Run conditions")

#: Starting points for a new experiment's detectors: the examples, named for what they are. (key, label, example)
TEMPLATES = (
    ("coulex-cd-clovers", "Coulomb excitation of a target nucleus with a CD and clovers", "coulex_ni58"),
    ("elastic-dssd", "Elastic scattering with silicon pads", "alpha_on_gold"),
    ("elastic-array", "Elastic scattering with an array of strip detectors", "oxygen_on_lead_array"),
    ("blank", "Blank: one detector, chosen from the kinematics", None),
)

#: The Physics tab's sections, in the order of the physics register: (key, title, register page, the tests that
#: validate it). Each holds panels that need no run.
PHYSICS = (
    ("data", "Nuclear data: masses and level schemes", "physics-register/nuclear-data",
     ("test_nuclear_data.py", "test_nuclear_levels.py")),
    ("kinematics", "Two-body kinematics", "physics-register/kinematics", ("test_nuclear_kinematics.py",)),
    ("stopping", "Energy loss and straggling", "physics-register/stopping", ("test_nuclear_stopping.py",)),
    ("rutherford", "Rutherford scattering: the orbit and the safe distance", "physics-register/rutherford",
     ("test_nuclear_rutherford.py",)),
    ("detectors", "Detectors: solid angles and the γ-ray response", "physics-register/detectors",
     ("test_nuclear_detectors.py", "test_nuclear_response.py")),
    ("rates", "Count rates and events", "physics-register/rates-and-events",
     ("test_nuclear_validation.py", "test_nuclear_events.py")),
    ("coulex", "Coulomb excitation: the excitation, the correlation and the Doppler shift", "theory/coulex",
     ("test_nuclear_coulex.py", "test_nuclear_orientation.py", "test_nuclear_gamma_events.py")),
)

#: Which Physics section shows the physics behind a number: by the key of its explanation (the ? of the Plan, the
#: Analysis and the run record) or of a setup card.
PHYSICS_FOR = {
    "solid_angle": "detectors", "rate": "rates", "counts": "rates", "beam_time": "rates",
    "excitation_probability": "coulex", "gamma_efficiency": "detectors", "coincidences": "coulex",
    "doppler": "coulex", "area": "coulex", "yield": "coulex", "mean_probability": "coulex", "b_e2": "coulex",
    "weisskopf": "data", "beta2": "data", "q0": "data", "lifetime": "data", "uncertainty": "rates",
    "run_rate": "rates", "run_counts": "rates", "run_busy": "rates", "run_live": "rates",
    "run_coincidences": "coulex",
    "beam": "rutherford", "target": "stopping", "run": "rates", "detector": "detectors", "gamma": "detectors",
}


def physics_for(key: str) -> str:
    """The Physics section for an explanation's key, or a setup card's ("detector:2" → "detectors")."""
    base = key.split(":")[0]
    if base.startswith("gate_"):
        return "kinematics"  # a gate sorts the particles by their kinematic lines
    if base.startswith("compare_"):
        return "rates"
    return PHYSICS_FOR.get(base, "coulex")


#: What the runs do not simulate yet, with the backlog item that adds it: the realism is claimed only where it holds.
NOT_SIMULATED = (
    (71, "Cascades in part", "without a level scheme each excitation emits one γ ray; with one, the cascade is "
                             "followed, but the angular correlation between its successive γ rays is not"),
    (72, "Summing and pile-up in part", "γ rays of one cascade sum in a crystal, and pile-up follows the "
                                        "shaping time, but a summed or piled signal is one Gaussian, without the "
                                        "electronics' shapes; dead time is a simple correction"),
    (73, "Lifetimes", "every state decays in flight after leaving the target: no stopped or partly shifted "
                      "components"),
    (74, "Contaminant reactions", "no scattering on carbon, oxygen or other contaminants of the target and backing"),
    (74, "Beam halo", "the beam spot is a Gaussian of the setup's spot size; no halo on a frame or the chamber"),
)

#: What a new experiment can measure. (key, label)
MEASUREMENTS = (("elastic", "Elastic scattering"), ("coulex-target", "Coulomb excitation of the target"),
                ("coulex-beam", "Coulomb excitation of the beam"))


@dataclass
class Decision:
    """One card of the setup panel."""

    #: "beam", "target", "detector:N", "gamma:N" or "run" (N counting from 0, as in the scene).
    key: str
    group: str
    title: str
    #: The current state in one line, and what it leads to in a second.
    state: str
    consequence: str
    #: The planner's section name for edits ("beam", "detector 1", ...).
    section: str
    fields: list = field(default_factory=list)
    more: list = field(default_factory=list)
    values: dict = field(default_factory=dict)


def _num(x: float, digits: int = 3) -> str:
    return guide._fmt(x) if digits == 3 else f"{x:.{digits}g}"


def _deg(a: float) -> str:
    return f"{a:.0f}°"


def _sr(msr: float) -> str:
    return guide._solid_angle(msr)


def _rate(r: float) -> str:
    return f"{guide._fmt(r)}/s" if r >= 0.01 else f"{r:.2g}/s"


def _try(fn, default="—"):
    try:
        return fn()
    except Exception:  # noqa: BLE001 - a consequence line must never break the panel
        return default


def _where(theta_deg: float) -> str:
    if theta_deg >= 170:
        return "behind the target"
    if theta_deg > 100:
        return "backward"
    if theta_deg < 80:
        return "forward"
    return "at the side"


def _beam(p) -> Decision:
    d = p.draft["beam"]
    state = f"{d.get('nuclide', '?')} at {d.get('energy', '?')}" + (f", {d['current']}" if d.get("current") else "")

    def consequence():
        exp = p.experiment
        e = exp.beam.energy_mev
        a = exp.beam.A if hasattr(exp.beam, "A") else None
        t = p.trajectories()
        parts = [f"{guide._fmt(exp.beam.particles_per_second)} particles/s"]
        if a:
            parts.append(f"{e / a:.3g} MeV/u")
        parts.append(f"closest approach {t['d0_fm']:.3g} fm head-on (nuclear range "
                     f"{t['interaction_radius_fm']:.3g} fm)")
        return " · ".join(parts)

    return Decision("beam", "Beam", "Beam", state, _try(consequence), "beam", BEAM_FIELDS, BEAM_MORE, d)


def _target(p) -> Decision:
    d = p.draft
    t = d["target"]
    state = f"{t.get('material', '?')} {t.get('thickness', '')}".strip()
    if t.get("backing"):
        state += f" on {t['backing'].get('thickness', '')} {t['backing'].get('material', '')}".rstrip()
    reaction = d.get("reaction", {"type": "elastic"})
    if reaction.get("type") == "coulex":
        state += f" · Coulomb excitation of the {reaction.get('excite', 'target')} ({reaction.get('energy', '?')})"
    else:
        state += " · elastic scattering"

    def consequence():
        el = p.energy_loss()["layers"]
        lay = el[0]
        parts = [f"the beam loses {guide._energy(lay['loss_mev'])} in the target "
                 f"(straggling {guide._energy(lay['straggling_fwhm_mev'])} FWHM)"]
        if len(el) > 1:
            parts.append(f"{guide._energy(sum(x['loss_mev'] for x in el[1:]))} more in the backing")
        if p.experiment.excitation is not None:
            g = p.gamma()
            if g.get("available"):
                parts.append(f"σ = {g['total_mb']:.3g} mb, safe up to {g['max_safe_angle']:.0f}° (CM)")
        return " · ".join(parts)

    return Decision("target", "Target and reaction", "Target and reaction", state, _try(consequence), "target",
                    TARGET_FIELDS, TARGET_MORE, t)


def _detector(p, i: int, rows: dict) -> Decision:
    det = p.draft["detectors"][i]
    name = det.get("name") or f"D{i + 1}"
    shape = det.get("model") or det.get("shape", "circle")
    try:
        theta = float(str(det.get("theta", "0")).split()[0])
    except ValueError:
        theta = 0.0
    state = f"{name}: {shape} at {det.get('distance', '?')}, {_where(theta)} (θ {det.get('theta', '?')})"

    def consequence():
        row = rows[name]
        lo, hi = row["theta_range"]
        line = f"{_deg(lo)}–{_deg(hi)}, {_sr(row['solid_angle_msr'])}, {_rate(row['rate_per_s'])}"
        if row.get("coincidence_per_s") is not None:
            line += f", {_rate(row['coincidence_per_s'])} in coincidence with a γ ray"
        unsafe = p.safety().get(name)
        if unsafe:
            line += f" · {len(unsafe)} unsafe ring{'s' * (len(unsafe) > 1)}"
        return line

    size = set().union(*SHAPE_FIELDS.values())
    wanted = [f for f in DETECTOR_FIELDS if f[0] not in size or f[0] in SHAPE_FIELDS.get(det.get("shape", "circle"),
                                                                                         ())]
    return Decision(f"detector:{i}", "Particle detectors", name, state, _try(consequence), f"detector {i + 1}",
                    wanted, DETECTOR_MORE, det)


def _gamma_detector(p, i: int) -> Decision:
    gd = p.draft["gamma_detectors"][i]
    name = gd.get("name") or f"γ{i + 1}"
    real = "crystal_diameter" in gd or "model" in gd
    kind = gd.get("model") or ("crystals" if real else "disc")
    extras = [x for x, on in (("add-back", gd.get("addback")), ("BGO shield", gd.get("shield"))) if on]
    state = f"{name}: {kind} at θ {gd.get('theta', '?')}, {gd.get('distance', '?')}" + (
        f" ({', '.join(extras)})" if extras else "")

    def consequence():
        s = p.selection(f"gamma:{i}")
        line = (f"full-energy efficiency {100 * s['peak_efficiency']:.2g}% at {s['gamma_energy_kev']:.0f} keV, "
                f"{_deg(2 * s['half_angle_deg'])} wide")
        if s.get("coincidence_rate_per_s") is not None:
            line += f", {_rate(s['coincidence_rate_per_s'])} particle–γ coincidences"
        return line

    fields = [f for f in GAMMA_FIELDS if not (real and f[0] == "radius")]
    more = GAMMA_MORE + (CRYSTAL_FIELDS if real else []) + (CLOVER_FIELDS if _crystals(gd) == 4 else [])
    return Decision(f"gamma:{i}", "γ-ray detectors", name, state, _try(consequence), f"gamma detector {i + 1}",
                    fields, more, gd)


def _crystals(gd: dict) -> Optional[int]:
    if gd.get("crystals") is not None:
        return gd["crystals"] if isinstance(gd["crystals"], int) else None
    if gd.get("model"):
        from . import catalogue

        try:
            return catalogue.model(gd["model"], "gamma").fields.get("crystals")
        except ValueError:
            return None
    return None


def _run(p) -> Decision:
    r = p.draft["run"]
    state = f"{r.get('beam_time', '?')} of beam" + (f", {r['counts_wanted']} counts wanted"
                                                     if r.get("counts_wanted") else "")

    def consequence():
        t = p.rates()
        what = {"all": "counts", "excitations": "excitations", "coincidences": "coincidences"}[t["measured"]]
        timed = [x for x in t["rows"] if x["beam_time_s"] is not None]
        if t["counts_wanted"] and timed:
            slow = max(timed, key=lambda x: x["beam_time_s"])
            enough = slow["beam_time_s"] <= t["beam_time_s"]
            return (f"{slow['detector']} needs {guide._duration(slow['beam_time_s'])} for {t['counts_wanted']} "
                    f"{what}: the beam time is {'enough' if enough else 'too short'}")
        few = min(t["rows"], key=lambda x: x["counts_in_run"])
        return f"fewest {what}: {guide._fmt(few['counts_in_run'])} in {few['detector']}"

    return Decision("run", "Run conditions", "Run conditions", state, _try(consequence), "run", RUN_FIELDS, RUN_MORE,
                    r)


def decisions(planner) -> list:
    """The setup panel's cards, in the order the decisions are made (see :data:`GROUPS`)."""
    p = planner
    rows = _try(lambda: {r["detector"]: r for r in p.rates()["rows"]}, {})
    out = [_beam(p), _target(p)]
    out += [_detector(p, i, rows) for i in range(len(p.draft["detectors"]))]
    out += [_gamma_detector(p, i) for i in range(len(p.draft.get("gamma_detectors", [])))]
    out.append(_run(p))
    return out


def checks(planner) -> list:
    """Everything the top of the setup panel shows: the setup's problems and warnings, and whether the setup's
    physics changed since the current run. Each is {"level", "text"}."""
    out = [{"level": w.level, "text": w.text} for w in _try(planner.warnings, [])]
    status = planner.run_status()
    if status and status["stale"]:
        out.insert(0, {"level": "warning", "text": f"The setup changed since run {status['number']} "
                       f"({'; '.join(status['changes'][:3])}{'; …' if len(status['changes']) > 3 else ''}). "
                       "Its data stay as taken; take a new run to see the change."})
    return out


def status_strip(planner) -> list:
    """The parts of the status strip, left to right: {"key" (the card it opens), "label", "text"}."""
    p = planner
    d = p.draft
    n_p, n_g = len(d["detectors"]), len(d.get("gamma_detectors", []))
    dets = f"{n_p} particle" + (f" + {n_g} γ" if n_g else "")
    runs = p.runs()
    if p.run is not None:
        last = p.run.describe().split(";")[0] + f" (run {p.run.number}"
        when = p.run.summary.get("finished", "")[:16].replace("T", " ")
        last += f", {when})" if when else ")"
    elif runs:
        last = runs[-1]["describe"].split(";")[0]
    else:
        last = "no run yet"
    ws = _try(p.warnings, [])
    n_warn = sum(w.level in ("error", "warning") for w in ws)
    stale = bool(p.run_status() and p.run_status()["stale"])
    return [
        {"key": None, "label": "Experiment", "text": p.name},
        {"key": "beam", "label": "Beam", "text": f"{d['beam'].get('nuclide', '?')} at {d['beam'].get('energy', '?')}"},
        {"key": "target", "label": "Target", "text": f"{d['target'].get('material', '?')} "
                                                    f"{d['target'].get('thickness', '')}".strip()},
        {"key": "detector:0", "label": "Detectors", "text": dets},
        {"key": "run", "label": "Last run", "text": last + (" · setup changed since" if stale else "")},
        {"key": "checks", "label": "Warnings", "text": f"{n_warn} warning{'s' * (n_warn != 1)}" if n_warn
         else "no warnings"},
    ]


# -- new experiments ------------------------------------------------------------------------------------------------

def template(key: str):
    """The :class:`~physim.nuclear.Experiment` of a template (:data:`TEMPLATES`); "blank" needs a beam and target,
    see :func:`new_setup`."""
    from .experiment import Experiment

    found = next((t for t in TEMPLATES if t[0] == key), None)
    if found is None:
        raise KeyError(f"no template {key!r}; the templates are {[t[0] for t in TEMPLATES]}")
    if found[2] is None:
        return new_setup("4He", "5.5 MeV", "Au")
    return Experiment.example(found[2])


def new_setup(beam: str, energy: str, target: str, thickness: str = "0.5 mg/cm2", measure: str = "elastic",
              current: str = "1 pnA", beam_time: str = "24 h", title: Optional[str] = None):
    """A new experiment's setup from what the "New experiment" dialog asks: the beam and its energy, the target,
    and what to measure (:data:`MEASUREMENTS`). It gets one particle detector where the kinematics send the
    particles (behind the target for a light beam on a heavy target, forward otherwise); more come from the Add
    menus. For Coulomb excitation the first 2⁺ state is taken from the local ENSDF copy when there is one; without
    it the setup stays elastic and :func:`new_setup_notes` says what to enter."""
    from . import data
    from .experiment import SCHEMA, Experiment
    from .levels import LevelScheme
    from .names import parse_nuclide

    z1, a1 = parse_nuclide(beam)
    beam_name = data.nuclide((z1, a1)).name
    try:
        atoms = sorted(data.material(target).atoms, key=lambda x: -x[2])
        heavy = data.nuclide(atoms[0][:2])
    except (KeyError, ValueError, IndexError):
        heavy = None
    light_on_heavy = heavy is None or a1 < heavy.A
    place = guide.placement("ring" if light_on_heavy else "forward")
    d = {"schema": SCHEMA, "title": title or f"{beam_name} on {heavy.name if heavy else target}",
         "beam": {"nuclide": beam, "energy": energy, "current": current},
         "target": {"material": target, "thickness": thickness},
         "run": {"beam_time": beam_time, "counts_wanted": 1000},
         "reaction": {"type": "elastic"},
         "detectors": [place]}
    exp = Experiment.from_dict(copy.deepcopy(d))
    notes = []
    if measure in ("coulex-target", "coulex-beam"):
        role = "target" if measure == "coulex-target" else "beam"
        nuclide = (heavy.name if heavy else target) if role == "target" else beam_name
        try:
            scheme = LevelScheme.from_ensdf(nuclide, max_energy_kev=3000.0)
            twos = [n for n, lev in enumerate(scheme.levels) if n and lev.spin == 2 and scheme.b(0, n, "E2")]
            if not twos:
                raise ValueError(f"no 2+ state with a known B(E2) in {nuclide}")
            n = twos[0]
            b = scheme.b(0, n, "E2")
            d["levels"] = {role: scheme.to_dict()}
            d["reaction"] = {"type": "coulex", "excite": "target" if role == "target" else "projectile",
                             "energy": f"{scheme.levels[n].energy.value:.10g} keV", "multipolarity": "E2",
                             "b_up": f"{b:.6g} e2fm4"}
            exp = Experiment.from_dict(copy.deepcopy(d))
        except Exception as e:  # noqa: BLE001 - no ENSDF copy, or no data for the nuclide
            notes.append(f"Coulomb excitation of {nuclide}: its first 2⁺ state could not be looked up ({e}). Set "
                         "the reaction to Coulomb excitation and enter the state's energy and B(E2) in the Target "
                         "and reaction card.")
    exp._new_setup_notes = notes  # read by new_setup_notes
    return exp


def new_setup_notes(experiment) -> list:
    """What :func:`new_setup` could not fill in."""
    return list(getattr(experiment, "_new_setup_notes", []))


__all__ = ["GROUPS", "MEASUREMENTS", "NOT_SIMULATED", "PHYSICS", "PHYSICS_FOR", "physics_for", "STAGES", "TEMPLATES", "Decision", "checks", "decisions", "new_setup",
           "new_setup_notes", "status_strip", "template"]
