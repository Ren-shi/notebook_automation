"""The planner's guidance: help for every input, the steps of the guided workflow, and a short reading of every
result for the current setup ("how to read this").

Everything here is plain text and data, so the app only lays it out and the tests can check it covers every field and
every tab. The help is practical (what the value is, a typical range, what raising it does); the physics behind each
result is in the "Explain" panels (:data:`physim.nuclear.planner.EXPLAIN`) and the docs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import catalogue


@dataclass(frozen=True)
class Help:
    """Help for one input: what it is, its typical range, and what raising it does to the results."""

    what: str
    typical: str
    effect: str

    def text(self) -> str:
        lines = [self.what, f"Typical: {self.typical}"] + ([f"Raise it: {self.effect}"] if self.effect != "—" else [])
        return "\n".join(lines)


def _h(what, typical, effect):
    return Help(what, typical, effect)


#: Help for every input of the setup panel, by (section, field). Sections: "title", "beam", "target", "backing",
#: "run", "reaction", "detector" (particle detectors) and "gamma" (γ-ray detectors).
HELP = {
    ("title", ""): _h("A name for this plan; it heads the report.", "the reaction and the aim, e.g. "
                      "'Coulomb excitation of 58Ni'", "—"),
    # -- beam
    ("beam", "nuclide"): _h("The ion the accelerator delivers, as mass number and symbol.", "4He, 16O, 12C, 1H",
                            "a heavier, higher-Z beam scatters more strongly (cross section ∝ Z²) and loses energy "
                            "faster in the target."),
    ("beam", "energy"): _h("Kinetic energy of the beam on the target, total (MeV) or per nucleon (MeV/u).",
                           "a few MeV for α particles; 2–5 MeV/u for heavy ions below the Coulomb barrier",
                           "Rutherford rates fall as 1/E², particles reach the detectors with more energy, and above "
                           "the barrier nuclear forces spoil pure Coulomb scattering (a warning says when)."),
    ("beam", "current"): _h("Beam intensity on target: particle nA (pnA) counts ions, electrical nA (enA) counts "
                            "charge.", "0.1–10 pnA", "every rate rises in proportion and the beam time falls; "
                            "detectors near the beam may count too fast (pile-up warning)."),
    ("beam", "charge_state"): _h("Charge of the ions after the accelerator; only needed to convert an electrical "
                                 "current (enA) to particles.", "from the accelerator's tune, e.g. 6 for 16O",
                                 "with a current in enA, fewer particles per nA."),
    ("beam", "energy_spread"): _h("Spread of the beam energy (FWHM), from the accelerator.", "0.05–0.2 %",
                                  "wider peaks in every spectrum."),
    ("beam", "spot_size"): _h("Diameter of the beam spot on the target (FWHM).", "1–3 mm",
                              "a spread of scattering angles in each detector: wider peaks for close detectors."),
    # -- target
    ("target", "material"): _h("What the target is made of: an element, an isotope or a compound.",
                               "Au, 208Pb, 58Ni, CD2, Mylar", "—"),
    ("target", "thickness"): _h("Areal density (mg/cm², µg/cm²) or a length (nm, µm) of the target.",
                                "50 µg/cm² – 1 mg/cm²", "more counts in proportion, but more energy loss and "
                                "straggling: wider peaks and a lower mean energy."),
    ("target", "tilt"): _h("Rotation of the target about the vertical axis, away from facing the beam.", "0–45 deg",
                           "a longer path through the target for the beam (more counts, more loss) and a different "
                           "path out for each detector."),
    ("target", "density"): _h("Mass density, only needed when the thickness is a length.", "looked up for elements",
                              "more mass per area for the same length."),
    ("backing", "material"): _h("A support foil behind the target (beam side last), e.g. carbon.", "C, Au",
                                "—"),
    ("backing", "thickness"): _h("Areal density of the backing.", "10–40 µg/cm² for carbon",
                                 "more scattering off the backing (its own peaks) and more energy loss for particles "
                                 "leaving through it."),
    # -- run
    ("run", "beam_time"): _h("How long the beam is on target.", "8 h – 3 days", "more counts in proportion; the "
                             "statistical error falls as 1/√counts."),
    ("run", "counts_wanted"): _h("Counts you need in each detector for your result, used to work out the beam "
                                 "time.", "1000–10 000 (3% – 1% statistical error)", "a longer beam time "
                                 "needed."),
    # -- reaction
    ("reaction", "type"): _h("What happens in the target: elastic (Rutherford) scattering, or Coulomb excitation of "
                             "one state.", "elastic for target analysis and cross sections; Coulomb excitation to "
                             "measure B(E2)", "—"),
    ("reaction", "excite"): _h("Which nucleus is excited: the target or the beam.", "the target for a stable "
                               "target nucleus; the beam for radioactive beams", "—"),
    ("reaction", "multipolarity"): _h("Character of the transition from the 0⁺ ground state.", "E2 for the first "
                                      "2⁺ state", "—"),
    ("reaction", "energy"): _h("Energy of the excited state, from ENSDF.", "0.1–3 MeV",
                               "a less likely excitation (the adiabatic cutoff) and a larger energy difference "
                               "between elastic and inelastic peaks."),
    ("reaction", "b_up"): _h("Reduced transition probability B(Eλ↑) from the ground state, from ENSDF (e²b² or "
                             "e²fm⁴).", "0.01–1 e²b² for E2", "excitation rates rise in proportion."),
    # -- particle detectors
    ("detector", "shape"): _h("Rectangle: a strip detector (DSSD); circle: a single pad; annular: a CD-type ring "
                              "detector around the beam.", "annular at backward angles, rectangles at the sides",
                              "—"),
    ("detector", "name"): _h("A short name used in tables, plots and the report.", "A30, CD, DSSD-L", "—"),
    ("detector", "theta"): _h("Angle between the beam direction and the detector centre: below 90° is forward, "
                              "above 90° backward (180° is straight back, around the beam).", "20–170 deg",
                              "Rutherford rates fall steeply with angle (1/sin⁴(θ/2)), so backward detectors count "
                              "slowly and the beam time grows; but they see the closest collisions, where Coulomb "
                              "excitation is most likely."),
    ("detector", "phi"): _h("Angle around the beam axis (0 = +x, 90 = +y).", "0, 90, 180, 270 deg",
                            "moves the detector around the beam; rates do not change for an unpolarised beam."),
    ("detector", "distance"): _h("Distance from the target to the detector centre.", "30–200 mm",
                                 "a smaller solid angle (∝ 1/d²): fewer counts, but a narrower angular range and "
                                 "sharper peaks."),
    ("detector", "width"): _h("Width of a rectangular detector.", "50 mm", "more solid angle and counts, a wider "
                              "angular range."),
    ("detector", "height"): _h("Height of a rectangular detector.", "50 mm", "more solid angle and counts."),
    ("detector", "strips_x"): _h("Number of strips across the width (each a separate angle bin).", "16",
                                 "finer angle bins, fewer counts in each."),
    ("detector", "strips_y"): _h("Number of strips across the height.", "16", "finer bins, fewer counts in each."),
    ("detector", "radius"): _h("Radius of a circular detector.", "5–10 mm", "more solid angle (∝ r²), a wider "
                               "angular range."),
    ("detector", "inner_radius"): _h("Inner radius of an annular detector; the beam passes through the hole.",
                                     "9–24 mm", "a smaller angular range, away from the beam."),
    ("detector", "outer_radius"): _h("Outer radius of an annular detector.", "41–48 mm", "more solid angle and a "
                                     "wider angular range."),
    ("detector", "rings"): _h("Number of rings (each a separate polar-angle bin).", "16", "finer angle bins, fewer "
                              "counts in each."),
    ("detector", "sectors"): _h("Number of sectors around the ring.", "16–24", "finer azimuthal bins."),
    ("detector", "thickness"): _h("Thickness of the silicon.", "100–1000 µm", "stops more energetic particles; a "
                                  "particle that is not stopped deposits only part of its energy (warning)."),
    ("detector", "dead_layer"): _h("Inactive layer on the front of the detector (contact and window).",
                                   "0.05–0.5 µm", "a lower measured energy and slightly wider peaks, most for heavy "
                                   "ions."),
    ("detector", "resolution"): _h("Energy resolution (FWHM) of the detector and electronics.", "15–50 keV",
                                   "wider peaks: close peaks may no longer separate."),
    ("detector", "threshold"): _h("Lowest energy the electronics record.", "100–500 keV", "low-energy particles "
                                  "(and noise) are cut; too high a threshold loses slow recoils."),
    ("detector", "material"): _h("Detector material; silicon unless you set another.", "Si", "—"),
    # -- γ-ray detectors
    ("gamma", "name"): _h("A short name used in the Doppler table.", "Ge90, Clover1", "—"),
    ("gamma", "theta"): _h("Angle between the beam and the γ detector.", "25–155 deg", "the Doppler shift is "
                           "largest forward and backward (cos θ) and smallest near 90°."),
    ("gamma", "phi"): _h("Angle around the beam axis.", "keep clear of the particle detectors",
                         "changes the angle to each particle detector, and so the Doppler shift per pair."),
    ("gamma", "distance"): _h("Distance from the target to the crystal face.", "100–250 mm", "a smaller opening "
                              "angle: less Doppler broadening, but less efficiency."),
    ("gamma", "radius"): _h("Radius of the crystal face.", "30–40 mm", "more efficiency, but a wider opening angle "
                            "and more Doppler broadening."),
    ("gamma", "crystals"): _h("Number of crystals: 1, or 4 for a clover.", "4 for a clover",
                              "each crystal is corrected for the Doppler shift on its own."),
    ("gamma", "crystal_diameter"): _h("Diameter of one crystal.", "50 mm for a EUROGAM-type clover",
                                      "more efficiency, but a wider opening angle per crystal."),
    ("gamma", "crystal_length"): _h("Length of one crystal, behind its front face.", "70 mm for a EUROGAM-type "
                                    "clover", "more efficiency at high γ-ray energy."),
    ("gamma", "crystal_pitch"): _h("Distance between the centres of neighbouring crystals of a clover.",
                                   "about 45 mm: check the detector's drawing", "moves each crystal's angle."),
    ("gamma", "housing_side"): _h("Side of the square housing, which stops particles.", "about 100 mm for a clover",
                                  "hides more of whatever lies behind it."),
    ("gamma", "window_gap"): _h("Distance from the housing's front window to the crystals.", "a few mm", "—"),
    ("gamma", "material"): _h("What the crystal is made of: Ge or LaBr3.", "Ge",
                              "LaBr3 is faster and more efficient per volume, with about ten times worse resolution."),
    ("gamma", "absorbers"): _h("Material between the target and the detector, such as a lead or copper sheet "
                               "against X-rays: a material and its thickness, several separated by commas.",
                               "Pb 1 mm, Cu 0.5 mm", "fewer low-energy γ rays and X-rays reach the crystal; "
                               "high-energy γ rays are hardly reduced."),
    ("gamma", "resolution"): _h("Energy resolution (FWHM) of the crystal at 1332 keV; it is scaled to other energies.", "2–3 keV at 1.3 MeV for "
                                "germanium", "wider γ peaks (added to the Doppler broadening)."),
    ("gamma", "efficiency"): _h("Full-energy-peak efficiency of this detector for the γ ray, as a percentage of all "
                                "γ rays emitted (from a source measurement or the array's specification). Left "
                                "empty, the planner uses the typical response of such a crystal at this distance.",
                                "0.5–3 % per germanium crystal at 1.3 MeV",
                                "more particle–γ coincidences: a shorter beam time."),
}


def help_for(section: str, field: str) -> Optional[Help]:
    """Help for a field of a section as the planner names it ("detector 2" and "gamma detector 1" included)."""
    if section.startswith("gamma detector"):
        section = "gamma"
    elif section.startswith("detector"):
        section = "detector"
    return HELP.get((section, field))


# ---------------------------------------------------------------------------------------------------------------
# Where to put a new particle detector


@dataclass(frozen=True)
class Placement:
    """A ready-made particle detector for one region of angles."""

    key: str
    label: str
    #: Why one would put a detector there.
    why: str
    #: Name prefix; a number is added.
    prefix: str
    #: The detector's fields, as in a setup file.
    fields: tuple


PLACEMENTS = (
    Placement("forward", "Forward strip detector (45°)",
              "Many counts: the beam scattered a little, and recoils from the target. Fast, high-energy particles.",
              "F", (("shape", "rectangle"), ("theta", "45 deg"), ("distance", "100 mm"), ("width", "50 mm"),
                    ("height", "50 mm"), ("strips_x", 16), ("strips_y", 16), ("thickness", "500 um"),
                    ("resolution", "30 keV"), ("threshold", "300 keV"))),
    Placement("side", "Side pad (90°)", "A middle angle: moderate rates, a check on the angular distribution.",
              "S", (("shape", "circle"), ("theta", "90 deg"), ("distance", "80 mm"), ("radius", "5 mm"),
                    ("thickness", "300 um"), ("resolution", "20 keV"), ("threshold", "200 keV"))),
    Placement("backward", "Backward pad (150°)",
              "Few counts and a long beam time, but the closest collisions: the strongest test of Rutherford "
              "scattering, and where Coulomb excitation is most likely.",
              "B", (("shape", "circle"), ("theta", "150 deg"), ("distance", "50 mm"), ("radius", "5 mm"),
                    ("thickness", "300 um"), ("resolution", "20 keV"), ("threshold", "200 keV"))),
    Placement("ring", "Backward ring around the beam (CD, 180°)",
              "An annular detector the beam passes through, covering about 125–165°: the most solid angle at "
              "backward angles, the usual choice for Coulomb excitation.",
              "CD", (("shape", "annular"), ("theta", "180 deg"), ("distance", "30 mm"), ("inner_radius", "9 mm"),
                     ("outer_radius", "41 mm"), ("rings", 16), ("sectors", 24), ("thickness", "300 um"),
                     ("resolution", "30 keV"), ("threshold", "300 keV"))),
    Placement("s3", "Micron S3 around the beam (180°)",
              "The real detector: 24 rings and 32 sectors on an active area of 22 to 70 mm in diameter, with its "
              "circuit board. Set the thickness of the one you have.",
              "S3-", (("model", "S3"), ("theta", "180 deg"), ("distance", "30 mm"), ("resolution", "30 keV"),
                      ("threshold", "300 keV"))),
)

#: Ready-made γ-ray detectors for the app's Add menu: (key, label, why, name prefix, fields).
GAMMA_PLACEMENTS = (
    ("disc", "Single germanium crystal (a disc)", "A simple detector: one crystal face of a radius you choose.",
     "Ge", (("theta", "90 deg"), ("phi", "90 deg"), ("distance", "120 mm"), ("radius", "35 mm"),
            ("resolution", "2.5 keV"))),
    ("clover", "Germanium clover, four 50 × 70 mm crystals",
     "The EUROGAM, EUROBALL and AFRODITE type. Each of the four crystals gets its own Doppler correction.",
     "Clover", (("model", "clover"), ("theta", "135 deg"), ("phi", "90 deg"), ("distance", "200 mm"))),
    ("labr3", "LaBr₃(Ce), 2 × 2 inch", "A fast scintillator: good timing, about 2% resolution at 1.3 MeV.",
     "LaBr", (("model", "LaBr3_2x2"), ("theta", "90 deg"), ("phi", "-90 deg"), ("distance", "150 mm"))),
)


def gamma_placement(key: str, taken=()) -> dict:
    """The fields of a new γ-ray detector from :data:`GAMMA_PLACEMENTS`, with its model's values written out and a
    name that does not clash with ``taken``."""
    _, _, _, prefix, fields = next(x for x in GAMMA_PLACEMENTS if x[0] == key)
    k = 1
    while f"{prefix}{k}" in set(taken):
        k += 1
    return dict(catalogue.with_defaults(dict(fields), "gamma"), name=f"{prefix}{k}")


def placement(key: str, taken=()) -> dict:
    """The fields of a new detector from :data:`PLACEMENTS`, named so it does not clash with ``taken``."""
    p = next(x for x in PLACEMENTS if x.key == key)
    k = 1
    while f"{p.prefix}{k}" in set(taken):
        k += 1
    # A model's values are written out, so the setup panel shows them and they can be edited.
    return dict(catalogue.with_defaults(dict(p.fields), "particle"), name=f"{p.prefix}{k}")


# ---------------------------------------------------------------------------------------------------------------
# The guided workflow


@dataclass(frozen=True)
class Step:
    """One step of the guided workflow."""

    key: str
    title: str
    #: What the user decides in this step, in one or two sentences.
    intro: str
    #: Setup sections edited here ("beam", "target", "backing", "run", "reaction", "detectors", "gamma_detectors").
    sections: tuple
    #: Result tabs shown beside the step (see :data:`physim.nuclear.planner.TABS`).
    tabs: tuple
    #: Only for Coulomb excitation.
    coulex_only: bool = False


STEPS = (
    Step("goal", "What do you want to measure?",
         "Choose the kind of experiment. Each starts from a complete example that you then change step by step, so "
         "every result is filled in from the start. For Coulomb excitation, this is also where the excited state "
         "is set.", ("reaction",), ()),
    Step("beam", "Beam",
         "Which ion, at what energy and intensity. The energy sets how the particles scatter and what they arrive "
         "with; the current sets every rate.", ("beam",), ("kinematics",)),
    Step("target", "Target",
         "What the beam hits and how thick it is. A thicker target gives more counts but costs energy and widens "
         "the peaks.", ("target", "backing"), ("energy_loss",)),
    Step("detectors", "Particle detectors",
         "Where the silicon detectors sit and how big they are. Forward angles count fast; backward angles count "
         "slowly but see the closest collisions, which Coulomb excitation needs. Distance and size set the solid "
         "angle and how sharp the peaks are.", ("detectors",), ("geometry", "kinematics")),
    Step("gamma", "γ-ray detectors",
         "Where the germanium detectors sit. Their angle to the beam and to the particle detectors sets the Doppler "
         "shift of the γ ray.", ("gamma_detectors",), ("gamma",), coulex_only=True),
    Step("rates", "Rates and beam time",
         "How many counts each detector collects in your beam time, and how long you need for the counts you want.",
         ("run",), ("rates",)),
    Step("spectra", "Spectra",
         "What the measured energy spectra will look like: which peaks, how wide, and whether they separate.", (),
         ("spectra",)),
    Step("report", "Report",
         "Download the beam-time report: every number above, ready for a proposal or the logbook.", (),
         ("report",)),
)

#: The two kinds of experiment the first step offers, with the example each starts from.
GOALS = {
    "elastic": ("Rutherford (elastic) scattering",
                "Scattered beam particles and recoils: cross sections, target thickness and composition, or "
                "calibrating detectors.", "alpha_on_gold"),
    "coulex": ("Coulomb excitation",
               "Excite a low-lying state electromagnetically, below the Coulomb barrier, and detect its γ ray with "
               "the scattered particle: B(E2) values and collectivity.", "coulex_ni58"),
}


def steps_for(planner) -> list:
    """The steps that apply to the planner's setup, each with its setup problems (empty when the step is fine)."""
    coulex = planner.draft.get("reaction", {}).get("type") == "coulex"
    out = []
    for s in STEPS:
        if s.coulex_only and not coulex:
            continue
        out.append({"step": s, "problems": [p for p in planner.problems if _section_of(p) in s.sections]})
    return out


def _section_of(problem: str) -> str:
    """Which setup section a problem message is about (as named in :attr:`Step.sections`)."""
    head = problem.split(":", 1)[0]
    if head.startswith("gamma detector") or head.startswith("gamma_detectors"):
        return "gamma_detectors"
    if head.startswith("detector") or head.startswith("no detectors"):
        return "detectors"
    if head.startswith("target.backing"):
        return "backing"
    for s in ("beam", "target", "run", "reaction"):
        if head.startswith(s) or f"[{s}]" in head:
            return s
    return "other"


# ---------------------------------------------------------------------------------------------------------------
# How to read each result


def _fmt(x: float) -> str:
    """Three significant figures without exponents: 22,700 rather than 2.27e+04."""
    if x == 0 or not np.isfinite(x):
        return f"{x:g}"
    if abs(x) >= 1000:
        digits = int(np.floor(np.log10(abs(x)))) - 2
        return f"{round(x, -digits):,.0f}"
    return f"{x:.3g}"


def _energy(mev: float) -> str:
    return f"{mev:.3g} MeV" if abs(mev) >= 1 else f"{mev * 1e3:.3g} keV"


def _solid_angle(msr: float) -> str:
    return f"{msr / 1e3:.3g} sr" if msr >= 1000 else f"{msr:.3g} msr"


def _duration(seconds: float) -> str:
    if not np.isfinite(seconds):
        return "never"
    if seconds >= 86400:
        return f"{seconds / 86400:.3g} days"
    if seconds >= 3600:
        return f"{seconds / 3600:.3g} h"
    if seconds >= 60:
        return f"{seconds / 60:.3g} min"
    return f"{seconds:.3g} s"


def reading(planner, tab: str) -> list:
    """One to three sentences on what to look for in a result tab, with this setup's numbers."""
    fn = _READINGS.get(tab)
    return fn(planner) if fn else []


def _geometry(p) -> list:
    g = p.geometry()
    dets = g["detectors"]
    if not dets:
        return []
    big = max(dets, key=lambda d: d["solid_angle_msr"])
    out = [f"The scene shows the experiment to scale: the beam (red) comes in along the axis and meets the target "
           f"at the origin; the particle detectors are blue. {big['name']} covers the most solid angle "
           f"({_solid_angle(big['solid_angle_msr'])})."]
    lo = min(d["theta_range"][0] for d in dets)
    hi = max(d["theta_range"][1] for d in dets)
    out.append(f"Together they cover {lo:.0f}° to {hi:.0f}° from the beam. Turn the view with the mouse to check "
               "nothing blocks another detector or the beam.")
    if g["gamma_detectors"]:
        out.append("The γ-ray detectors are the crystals drawn in bronze.")
    if not any(d["theta_range"][1] > 90 for d in dets):
        out.append("All detectors are forward of 90°. Backward angles count slowly but see the closest collisions: "
                   "add one under Particle detectors (Backward pad, or Backward ring around the beam).")
    if p.experiment.excitation is not None:
        r = p.rates()["rows"]
        share = {x["detector"]: x["excitation_per_s"] / x["rate_per_s"] for x in r if x["rate_per_s"] > 0}
        if share:
            best = max(share, key=share.get)
            out.append(f"{best} has the largest share of excitation events: 1 in {1 / share[best]:,.0f} of its "
                       "particles. Close collisions excite the state: the beam scattered backward, or a target "
                       "recoil sent forward.")
    return out


def _kinematics(p) -> list:
    k = p.kinematics()
    out = ["Each curve is the energy a particle has after scattering into that lab angle; the shaded bands are "
           "the angles your detectors cover, so a band shows the energies that detector will see."]
    ej = [c for c in k["curves"] if c["particle"] == "ejectile" and not c.get("inelastic")]
    if ej and ej[0]["max_angle"] < 179.9:
        out.append(f"The beam is heavier than the target nucleus, so it scatters no further than "
                   f"{ej[0]['max_angle']:.1f}°: detectors beyond that see only recoils.")
    if any(c.get("inelastic") for c in k["curves"]):
        out.append("Dotted curves are the particles after exciting the state: slightly lower in energy at every "
                   "angle.")
    return out


def _rates(p) -> list:
    r = p.rates()
    rows = r["rows"]
    if not rows:
        return []
    top = max(rows, key=lambda x: x["rate_per_s"])
    low = min(rows, key=lambda x: x["rate_per_s"])
    out = [f"{top['detector']} counts fastest ({_fmt(top['rate_per_s'])}/s) and {low['detector']} slowest "
           f"({_fmt(low['rate_per_s'])}/s): scattering falls steeply with angle."]
    need = {x["detector"]: x["beam_time_s"] for x in rows if x["beam_time_s"] is not None}
    if r["measured"] != "all":
        return out + _coulex_rates(r, need)
    if need:
        worst = max(need, key=need.get)
        have = r["beam_time_s"]
        verdict = "" if have is None else (" Your planned beam time is enough." if need[worst] <= have else
                                           f" Your planned {_duration(have)} is not enough.")
        out.append(f"For {r['counts_wanted']:g} counts in every detector you need {_duration(need[worst])} of "
                   f"beam, set by {worst}.{verdict}")
    if top["rate_per_s"] > 5000:
        out.append(f"{top['detector']} is above about 5000/s, where pile-up and dead time start: lower the current "
                   "or move it back.")
    return out


def _coulex_rates(r: dict, need: dict) -> list:
    """The rates reading for Coulomb excitation: what is counted, the beam time it needs, and the efficiency."""
    out = []
    if r["measured"] == "coincidences":
        eff = r["gamma_efficiency"]
        out.append(f"For Coulomb excitation what counts is an excitation event seen in a particle detector together "
                   f"with its γ ray: the γ detectors catch {eff:.1%} of the γ rays"
                   + (" (from the typical response of such crystals: give a measured Efficiency, or an "
                      "efficiency curve, for your own detectors)." if r["gamma_efficiency_typical"] else "."))
        what = "particle–γ coincidences"
    else:
        out.append("For Coulomb excitation what counts is the excitation events; add γ-ray detectors to count the "
                   "particle–γ coincidences the measurement uses.")
        what = "excitation events"
    if need:
        best = min(need, key=need.get)
        have = r["beam_time_s"]
        verdict = "" if have is None else (" Your planned beam time is enough there." if need[best] <= have else
                                           f" Your planned {_duration(have)} is not enough.")
        out.append(f"For {r['counts_wanted']:g} {what}, {best} needs the least beam: {_duration(need[best])}."
                   + verdict)
    return out


def _energy_loss(p) -> list:
    t = p.energy_loss()["layers"][0]
    loss, e0 = t["loss_mev"], t["energy_in_mev"]
    out = [f"The beam loses {_energy(loss)} crossing the target ({loss / e0:.2%} of its energy). "
           "Reactions happen at every depth, so the scattered particles carry this spread into the peaks."]
    if loss / e0 > 0.05:
        out.append("That is a large loss: the cross section changes across the target. A thinner target gives "
                   "sharper peaks.")
    return out


def _spectra(p) -> list:
    rates = p._rates()
    out = ["Each histogram is one detector's measured-energy spectrum from simulated events; lines mark the "
           "expected peaks."]
    for g in rates.array:
        peaks = rates.peaks(g.name)
        if len(peaks) < 2:
            continue
        a, b = sorted(peaks[:2], key=lambda x: -x.mean)
        gap = a.mean - b.mean
        width = max(a.fwhm, b.fwhm)
        verdict = "separate cleanly" if gap > 1.5 * width else "overlap" if gap < 0.7 * width else "just separate"
        out.append(f"In {g.name}, the two strongest peaks ({a.channel}, {a.particle}; {b.channel}, {b.particle}) "
                   f"are {_energy(gap)} apart with widths up to {_energy(width)} FWHM: they {verdict}.")
        # Over a whole detector the change of energy with angle dominates; one strip or ring sees a narrow range.
        sigma = max(a.sigma, b.sigma)
        geo = max(a.components.get("geometry", 0.0), b.components.get("geometry", 0.0))
        if verdict != "separate cleanly" and len(g.segments) > 1 and geo > 0.7 * sigma:
            out.append(f"Most of that width is the spread of angles across {g.name}; each of its "
                       f"{len(g.segments)} segments covers a narrow range, so look at the peaks strip by strip "
                       "(or ring by ring) in the analysis.")
        break
    return out


def _trajectories(p) -> list:
    t = p.trajectories()
    return [f"Coulomb orbits for a range of impact parameters: the closer the beam passes the nucleus, the further "
            f"it is deflected. The closest approach here is {t['d0_fm']:.3g} fm head-on; nuclear forces start at "
            f"about {t['interaction_radius_fm']:.3g} fm."]


def _gamma(p) -> list:
    g = p.gamma()
    if not g["available"]:
        return []
    out = [f"The excitation probability rises toward backward angles (close collisions); in total "
           f"{g['total_mb']:.3g} mb."]
    parts = [r for r in g["particles"] if r["where"] == "centre"]
    if parts:
        r = min(parts, key=lambda x: x["difference_mev"])
        out.append(f"At the centre of {r['detector']}, the {r['nuclide']} from an excitation event has "
                   f"{_energy(r['difference_mev'])} less than an elastic one: compare this with the detector "
                   "resolution to see whether the particle energy alone tells them apart.")
    if g["doppler"]:
        d = max(g["doppler"], key=lambda x: abs(x["shift_kev"]))
        out.append(f"The largest Doppler shift is {d['shift_kev']:+.3g} keV ({d['particle_detector']} with "
                   f"{d['gamma_detector']}). Uncorrected, that pair's γ peak is {d['fwhm_kev']:.3g} keV wide; the "
                   "Doppler correction, using the particle's direction and β, brings it back toward the crystal's "
                   "resolution.")
    return out


def _report(p) -> list:
    ws = [w for w in p.warnings() if w.level in ("error", "warning")]
    if not ws:
        return ["Nothing in the setup needs attention: the report is ready to share."]
    return [f"The report opens with the {len(ws)} warning{'s' if len(ws) != 1 else ''} shown at the top of the "
            "page. Fix them first, or explain them in your proposal."]


_READINGS = {"geometry": _geometry, "kinematics": _kinematics, "rates": _rates, "energy_loss": _energy_loss,
             "spectra": _spectra, "trajectories": _trajectories, "gamma": _gamma, "report": _report}

__all__ = ["GAMMA_PLACEMENTS", "GOALS", "HELP", "Help", "PLACEMENTS", "Placement", "STEPS", "Step", "help_for", "gamma_placement", "placement", "reading",
           "steps_for"]
