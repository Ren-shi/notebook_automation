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

**The steps** of an experiment are the app's step bar, left to right: :data:`STEPS` (Setup → Plan → Run → Data →
Analysis → Report). :func:`steps` gives each step's state and the line under its name, and :func:`next_step` what
the button at the foot of a step leads to. :data:`STAGES` are the steps and the Physics reference, which the app
opens from its header.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Optional

from . import guide

#: Everything the app's main area can show: the steps of an experiment, then the Physics reference. (key, label)
STAGES = (("setup", "Setup"), ("plan", "Plan"), ("run", "Run"), ("data", "Data"), ("analysis", "Analysis"),
          ("report", "Report"), ("physics", "Physics"))
#: The steps of an experiment, in the order they are taken: the app's step bar. (key, label)
STEPS = STAGES[:6]

# -- the fields of each decision: (field, label, placeholder) -------------------------------------------------------
BEAM_FIELDS = [("nuclide", "Nuclide", "4He"), ("energy", "Energy", "5.5 MeV or 4 MeV/u"),
               ("current", "Current", "1 pnA or 10 enA")]
BEAM_MORE = [("charge_state", "Charge state", "6"), ("energy_spread", "Energy spread (FWHM)", "0.1 %"),
             ("spot_size", "Spot size (FWHM)", "2 mm"), ("halo_fraction", "Halo (share of the beam)", "1 %"),
             ("halo_radius", "Halo radius (frame aperture)", "8 mm")]
TARGET_FIELDS = [("material", "Material", "Au, 208Pb, CD2"), ("thickness", "Thickness", "0.5 mg/cm2 or 1 um")]
TARGET_MORE = [("position", "Position along the beam", "0 mm"), ("tilt", "Tilt", "0 deg"),
               ("density", "Density", "19.3 g/cm3"),
               ("contaminants", "Contaminants", "C 10 ug/cm2, O 3 ug/cm2")]
BACKING_FIELDS = [("material", "Backing material", "C"), ("thickness", "Backing thickness", "20 ug/cm2")]
REACTION_FIELDS = [("energy", "State energy", "1.454 MeV"), ("b_up", "B(Eλ↑)", "0.0695 e2b2 or 695 e2fm4"),
                   ("e1_polarisation", "E1 polarisation (GOSIA's d)", "0.0038")]
REACTION_TYPES = {"elastic": "Elastic (Rutherford) scattering", "coulex": "Coulomb excitation",
                  "fusion-evaporation": "Fusion-evaporation (residues from a channel table)"}
#: The fusion-evaporation reaction's own fields, and those of each exit channel (backlog item 75).
FUSION_FIELDS = [("temperature", "Temperature (default from E*)", "1.5 MeV"),
                 ("alignment", "Alignment σ/J", "0.3"), ("side_feeding", "Side feeding", "25 %")]
NEUTRON_WALL_FIELDS = [("cells", "Cells", "50"), ("distance", "Distance", "510 mm"), ("pitch", "Pitch", "130 mm"),
                       ("theta_min", "From θ", "0 deg"), ("radius", "Cell radius", "63.5 mm"),
                       ("depth", "Cell depth", "150 mm"), ("threshold", "Threshold (electron-equivalent)", "0.1 MeV"),
                       ("material", "Scintillator", "EJ-301"), ("layout", "Layout", "hexagonal"),
                       ("gamma_efficiency", "γ rays above threshold", "40 %"),
                       ("gamma_as_neutron", "γ rays called neutrons", "1 %"),
                       ("gamma_multiplicity", "γ rays per event (default from entry spin)", "J/2 + 4"),
                       ("particle_threshold", "Particle-detector threshold (if unset)", "1 MeV")]
CHANNEL_FIELDS = [("particles", "Particles", "4n, p2n, a n"), ("cross_section", "Cross section", "300 mb"),
                  ("entry_spin", "Entry spin (ħ)", "20"), ("source", "Source", "PACE4")]
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
GROUPS = ("Beam", "Target and reaction", "Particle detectors", "γ-ray detectors", "Ancillary detectors",
          "Lifetime measurement", "Run conditions")

#: The lifetime card's choices: the method, and the kind of the plunger's foil after the gap.
LIFETIME_METHODS = {"rdds": "Plunger (RDDS): recoil distance", "dsam": "Line shapes (DSAM): Doppler-shift attenuation"}
FOIL_KINDS = {"stopper": "Stopper: stops the recoil", "degrader": "Degrader: lets it through, slower",
              "differential": "Differential: a degrader, then a second foil that stops it"}
#: The lifetime card's fields, as (section, field, label, placeholder); sections are the setup's own
#: (``lifetime``, ``lifetime.stopper``, ``lifetime.degrader2``, ``lifetime.deorientation``).
PLUNGER_FIELDS = [("lifetime", "distances", "Distances (target to stopper)", "10, 20, 50, 100, 200, 500 um"),
                  ("lifetime", "offset", "Offset at electrical contact", "0 um")]
FOIL_FIELDS = [("lifetime.stopper", "material", "Foil material", "Au"),
               ("lifetime.stopper", "thickness", "Foil thickness", "10 mg/cm2")]
FOIL_MORE = [("lifetime.stopper", "density", "Foil density (for a thickness in µm)", "19.3 g/cm3")]
DEGRADER2_FIELDS = [("lifetime.degrader2", "material", "Second foil material", "Au"),
                    ("lifetime.degrader2", "thickness", "Second foil thickness", "10 mg/cm2"),
                    ("lifetime.degrader2", "gap", "Gap after the first foil", "1000 um")]
#: Deorientation in vacuum, in either of its two forms.
DEORIENTATION_COMMON = [("lifetime.deorientation", "alpha2", "α₂ (hard core)", "0.18"),
                        ("lifetime.deorientation", "alpha4", "α₄ (hard core)", "0.24")]
DEORIENTATION_RATES = [("lifetime.deorientation", "gamma2", "Γ₂", "0.142 /ps"),
                       ("lifetime.deorientation", "gamma4", "Γ₄", "0.314 /ps")]
DEORIENTATION_G = [("lifetime.deorientation", "g", "g factor", "0.47"),
                   ("lifetime.deorientation", "c2", "C₂", "3.3 ps"),
                   ("lifetime.deorientation", "c4", "C₄", "1.5 ps")]
#: The ancillary detectors of a fusion-evaporation setup, as cards: {section: (title, fields, more fields, what
#: it is for)}. (The beam monitor is a particle detector.)
ANCILLARY_CARDS = {
    "neutron_wall": ("Neutron wall", [("cells", "Cells", "50"), ("distance", "Distance", "510 mm"),
                                      ("theta_min", "From θ", "4 deg"), ("material", "Scintillator", "EJ-301")],
                     [f for f in NEUTRON_WALL_FIELDS if f[0] not in ("cells", "distance", "theta_min", "material")],
                     "liquid-scintillator cells at forward angles: tags the xn channels"),
    "charged_particle_ball": ("Charged-particle ball",
                              [("elements", "Elements", "40"), ("radius", "Radius", "60 mm"),
                               ("material", "Material (Si or CsI)", "CsI"), ("threshold", "Threshold", "2 MeV")],
                              [("thickness", "Thickness", "3 mm"), ("theta_min", "From θ", "0 deg"),
                               ("theta_max", "To θ", "180 deg")],
                              "Si or CsI round the target: tags the p and α channels"),
    "recoil_separator": ("Recoil separator",
                         [("kind", "Kind (gas-filled, vacuum, acceptance)", "gas-filled"),
                          ("theta_acceptance", "Angular acceptance (half-angle)", "3 deg"),
                          ("momentum_acceptance", "Momentum acceptance", "10 %"), ("distance", "Entrance at", "600 mm")],
                         [("beam_suppression", "Beam suppression", "1e12"), ("focal_plane", "Focal plane", "DSSD")],
                         "keeps the residues inside its acceptance, with their mass at the focal plane"),
    "gamma_calorimeter": ("γ calorimeter",
                          [("elements", "Elements", "60"), ("radius", "Inner radius", "120 mm"),
                           ("material", "Material (BGO or NaI)", "BGO"), ("fold_min", "Fold at least", "6")],
                          [("thickness", "Thickness", "75 mm"), ("efficiency", "Efficiency per γ ray", "70 %")],
                          "a BGO or NaI shell: keeps the high-fold (high-spin) events"),
}

#: Starting points for a new experiment's detectors: the examples, named for what they are. (key, label, example)
#: The first is what the app opens when there is no experiment yet.
TEMPLATES = (
    ("coulex-beam-backward", "Coulomb excitation of the beam (62Ni on 194Pt) with silicon at backward angles",
     "coulex_ni62_pt194"),
    ("coulex-beam-clovers", "Coulomb excitation of the beam (62Ni on 194Pt) with silicon at backward angles and "
     "four clovers", "coulex_ni62_pt194_clovers"),
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
                             "followed with the angular correlation between its successive γ rays, but a "
                             "conversion electron's own correlation is not"),
    (71, "Excitation in several steps in part", "with a level scheme its levels are populated to all orders (a 4⁺ "
                                                "state through the 2⁺, the reorientation), each with its own "
                                                "orientation; without one, or with order = \"first\", only what one "
                                                "step from the ground state reaches. Every excitation keeps the "
                                                "kinematics of the [reaction] state, and elastic scattering is not "
                                                "reduced by what is excited"),
    (72, "Summing and pile-up in part", "γ rays of one cascade sum in a crystal, and pile-up follows the "
                                        "shaping time, but a summed or piled signal is one Gaussian, without the "
                                        "electronics' shapes; dead time is a simple correction"),
    (73, "Lifetimes in part", "a state with a half-life in the level scheme slows down and decays on its way, "
                              "but without angular straggling or time dilation; a state without one decays after "
                              "the target"),
    (74, "Contaminant reactions in part", "elastic scattering on the backing and on the [target] contaminants is "
                                          "simulated; fusion-evaporation and transfer channels are not"),
    (74, "Beam halo in part", "the halo spreads over the frame's aperture and blurs the lines; particles that hit "
                              "the frame itself are not followed"),
)

#: What a new experiment can measure. (key, label)
MEASUREMENTS = (("elastic", "Elastic scattering"), ("coulex-target", "Coulomb excitation of the target"),
                ("coulex-beam", "Coulomb excitation of the beam"), ("fusion-evaporation", "Fusion-evaporation"))

#: The fusion-evaporation templates: (key, label, what it holds, example the γ-ray array comes from, ancillaries).
FUSION_TEMPLATES = (
    ("fusion-gamma", "γ-ray array only", "5 germanium detectors and a beam monitor, as fusion_dy158",
     "fusion_dy158", ()),
    ("fusion-ball", "γ-ray array + charged-particle ball", "germanium + a CsI ball to tag the p and α channels",
     "fusion_dy158", ("charged_particle_ball",)),
    ("fusion-wall", "γ-ray array + neutron wall", "germanium + a liquid-scintillator wall at forward angles",
     "fusion_dy158", ("neutron_wall",)),
    ("fusion-separator", "γ-ray array + recoil separator", "germanium at the target, residues collected downstream",
     "fusion_dy158", ("recoil_separator",)),
    ("fusion-all", "γ-ray array + every ancillary", "the wall, the ball, the separator and a γ calorimeter",
     "fusion_dy158", ("neutron_wall", "charged_particle_ball", "recoil_separator", "gamma_calorimeter")),
    ("fusion-blank", "Blank", "no detectors yet: the γ-ray array is added on the Setup step", None, ()),
)

#: The lifetime templates: (key, label, what it holds, example the detectors come from, reaction, method). The
#: beam and target are the new experiment's own; the plunger's foil is made thick enough for its recoil.
LIFETIME_TEMPLATES = (
    ("lifetimes-coulex-plunger", "Plunger after Coulomb excitation",
     "ten clovers and a silicon detector at 0°, as coulex_ru98_plunger", "coulex_ru98_plunger", "coulex", "rdds"),
    ("lifetimes-fusion-plunger", "Plunger after fusion-evaporation",
     "eleven germanium detectors at 45° and 142°, as fusion_cr50_plunger", "fusion_cr50_plunger",
     "fusion-evaporation", "rdds"),
    ("lifetimes-fusion-dsam", "Line shapes (DSAM) after fusion-evaporation, backed target",
     "the germanium array and the gold backing of fusion_sm143_dsam", "fusion_sm143_dsam", "fusion-evaporation",
     "dsam"),
)


@dataclass(frozen=True)
class Kind:
    """A kind of experiment the start page offers."""

    key: str
    label: str
    #: What the experiment is, in a sentence or two, and what the planner gives for it.
    what: str
    gets: str
    #: Whether the planner can do it yet (the others are shown greyed, as kinds still to come).
    live: bool = True
    #: Its templates (keys of :data:`TEMPLATES` or :data:`FUSION_TEMPLATES`), its worked examples, its
    #: default beam and target, and the "What to measure" key of :func:`new_setup`.
    templates: tuple = ()
    examples: tuple = ()
    defaults: tuple = ()
    measure: str = ""


#: The kinds of experiment, in the start page's order: the three the planner does, then those still to come.
KINDS = (
    Kind("coulex", "Coulomb excitation",
         "Excite the beam or the target in the Coulomb field, below the barrier. Silicon for the scattered particle, "
         "germanium for the γ rays.", "B(E2)↑, Qs, excitation probability, γ yields per angle, Doppler correction",
         templates=("coulex-beam-backward", "coulex-beam-clovers", "coulex-cd-clovers", "blank"),
         examples=("coulex_ni62_pt194", "coulex_ni62_pt194_clovers", "coulex_ni58"),
         defaults=("62Ni", "237 MeV", "194Pt", "0.5 mg/cm2"), measure="coulex-target"),
    Kind("elastic", "Elastic (Rutherford) scattering",
         "Scatter a beam off a target and count it in silicon: target thickness, beam monitoring, where Rutherford "
         "holds and where nuclear forces begin.", "dσ/dΩ, rates per strip, energy loss, safe distance, beam time",
         templates=("elastic-dssd", "elastic-array", "blank"), examples=("alpha_on_gold", "oxygen_on_lead_array"),
         defaults=("4He", "5.5 MeV", "197Au", "0.5 mg/cm2"), measure="elastic"),
    Kind("fusion-evaporation", "Fusion-evaporation",
         "Fuse beam and target into a hot compound nucleus that boils off neutrons, protons and α particles; study "
         "the residues' γ rays, picked out by what was evaporated.",
         "channel cross sections, residue yields, γ-line rates, channel selection",
         templates=tuple(t[0] for t in FUSION_TEMPLATES), examples=("fusion_dy158", "fusion_tm165"),
         defaults=("12C", "60 MeV", "150Nd", "1 mg/cm2"), measure="fusion-evaporation"),
    Kind("transfer", "Transfer reactions", "(d,p), (t,p) and the like.", "", live=False),
    Kind("lifetimes", "Lifetimes (RDDS, DSAM)",
         "Measure how long a state lives: a plunger (recoil distance) or line shapes in a backed target (Doppler-shift "
         "attenuation), after Coulomb excitation or fusion-evaporation.",
         "τ per level, plunger distances, peak separation, F(τ), beam time per distance, decay curves, B(E2)",
         templates=tuple(t[0] for t in LIFETIME_TEMPLATES),
         examples=("coulex_ru98_plunger", "fusion_cr50_plunger", "fusion_sm143_dsam"),
         defaults=("98Ru", "300 MeV", "24Mg", "0.7 mg/cm2"), measure="lifetimes"),
    Kind("inelastic", "Inelastic scattering", "Nuclear excitation above the barrier.", "", live=False),
)


def kind_of(experiment) -> tuple:
    """The kind of an experiment, from its reaction: (key of :data:`KINDS`, the short label the app shows under
    its name). Coulomb excitation says which nucleus is excited; a setup with a lifetime measurement is of the
    "lifetimes" kind, with its method and reaction."""
    lt = getattr(experiment, "lifetime", None)
    if lt is not None:
        how = "plunger" if lt.method == "rdds" else "DSAM"
        after = "fusion" if experiment.reaction == "fusion-evaporation" else "Coulex"
        return "lifetimes", f"Lifetimes · {how} after {after}"
    if experiment.reaction == "fusion-evaporation":
        return "fusion-evaporation", "Fusion-evaporation"
    if experiment.reaction == "coulex":
        excite = getattr(experiment.excitation, "excite", None) or "target"
        return "coulex", "Coulex · " + ("beam" if excite in ("projectile", "beam") else excite)
    return "elastic", "Elastic"


def nucleus_href(nuclide: str, observable: Optional[str] = None, key: Optional[str] = None) -> str:
    """The address of a nucleus's page in the app ("/?nucleus=62Ni"), at one observable of one level or
    transition if given (``key`` as the results ledger writes it, or as one says it). Every part is
    percent-encoded ("+" as %2B), as the app's :func:`physim.nuclear.app.tab_nucleus.url` writes the address."""
    from urllib.parse import quote

    out = "/?nucleus=" + quote(str(nuclide), safe="")
    if observable is not None:
        out += "&observable=" + quote(str(observable), safe="")
        if key is not None:
            out += "&key=" + quote(str(key), safe="")
    return out


def studied(experiment) -> Optional[dict]:
    """The nucleus an experiment studies, as the rail shows it above the experiment's name: {"nuclide", "label"
    ("⁶²Ni"), "href" (its page)}. None for an experiment that studies none: elastic scattering, and
    fusion-evaporation without a studied channel (:func:`physim.nuclear.results.studied_nuclide`)."""
    from . import results
    from .logbook import nuclide_words

    try:
        nuclide = results.studied_nuclide(experiment)
    except Exception:  # noqa: BLE001 -- an unfinished setup studies no nucleus yet
        nuclide = None
    if nuclide is None:
        return None
    return {"nuclide": nuclide, "label": nuclide_words(nuclide), "href": nucleus_href(nuclide)}


#: The kinds of experiment whose runs physim does not analyse back yet: (what the kind is called in a sentence,
#: what the experiment does have beside its runs, what an analysis of a run would need that physim does not have).
NOT_ANALYSED = {
    "fusion-evaporation": (
        "fusion-evaporation",
        "the plan's rates and the residues",
        "the residues' γ-ray lines fitted in the singles and the γγ-gated spectra, and their areas turned into "
        "the channels' cross sections through the efficiency and the channel selection"),
    "elastic": (
        "elastic scattering",
        "the plan's rates",
        "the counts in each strip fitted to the Rutherford cross section, for the target's thickness, the "
        "normalisation of the beam, or the angle where the cross section leaves Rutherford's"),
}


def not_analysed(experiment) -> Optional[dict]:
    """For an experiment whose runs physim does not analyse back yet (fusion-evaporation and elastic scattering
    without a lifetime measurement): {"kind", "what", "has", "needs", "line"}, where "line" says in a sentence
    what an analysis would need that physim does not have for that kind, and "has" what the experiment does
    have beside its runs (the Plan step shows it, and the report holds it). None for Coulomb excitation and for a
    lifetime measurement, which have an analysis."""
    try:
        kind = kind_of(experiment)[0]
    except Exception:  # noqa: BLE001 -- an unfinished setup
        return None
    if kind not in NOT_ANALYSED:
        return None
    what, has, needs = NOT_ANALYSED[kind]
    return {"kind": kind, "what": what, "has": has, "needs": needs,
            "line": f"physim does not analyse a run of {what} back yet: that would need {needs}."}


def kept_lifetimes(planner) -> Optional[dict]:
    """The last lifetime analysis kept with the experiment (:meth:`physim.nuclear.planner.Planner.kept_lifetimes`:
    there after the app is opened again), read once for each version of its file; None if there is none, the
    setup's physics changed since, or this physim keeps none. Never starts the analysis."""
    import os

    fn, folder = getattr(planner, "kept_lifetimes", None), getattr(planner, "experiment_folder", None)
    if fn is None or folder is None:
        return None
    try:
        stamp = os.stat(os.path.join(str(folder), "lifetime.json")).st_mtime_ns
    except OSError:
        return None
    cache = getattr(planner, "_cache", None)  # the planner's own: it goes when the setup changes
    key = ("workbench", "kept lifetimes")
    if isinstance(cache, dict) and cache.get(key, (None, None))[0] == stamp:
        return cache[key][1]
    try:
        out = fn()
    except Exception:  # noqa: BLE001 -- a file that cannot be read keeps nothing
        out = None
    if not (isinstance(out, dict) and out.get("available") and out.get("tau_ps") is not None):
        out = None
    if isinstance(cache, dict):
        cache[key] = (stamp, out)
    return out


def analysed_level(analysis: dict) -> str:
    """The level a lifetime analysis is of, as the report writes it ("2⁺₁"); the level scheme's own label for an
    analysis that does not say which line it read."""
    from . import results
    from .logbook import level_words

    line = analysis.get("transition")
    if line:
        try:
            x = line["initial"]
            return level_words(results.level_key(line["nuclide"], x["energy_kev"], x.get("spin"), x.get("parity")))
        except (KeyError, TypeError, ValueError):
            pass
    return str(analysis.get("level") or "the studied level")


def runs_not_read(planner, analysis: dict, runs: Optional[list] = None) -> list:
    """The names of this setup's runs that a kept lifetime analysis did not read, because they were taken after
    it (``runs``: :meth:`Planner.runs`, listed if not given). A plunger analysis reads the runs at a distance, a
    line-shape analysis every beam run. Empty for an analysis that does not name the runs it read."""
    read = set(analysis.get("runs") or [])
    if not read:
        return []
    plunger = planner.experiment.lifetime is not None and planner.experiment.lifetime.method == "rdds"
    return [r["label"] for r in (planner.runs() if runs is None else runs)
            if r["label"] not in read and not r.get("stale") and r.get("kind") not in ("source", "alignment")
            and not (plunger and r.get("distance_um") is None)]


def template_label(key: str) -> tuple:
    """A template's (label, what it holds), for :data:`TEMPLATES`, :data:`FUSION_TEMPLATES` and
    :data:`LIFETIME_TEMPLATES`."""
    for t in FUSION_TEMPLATES + LIFETIME_TEMPLATES:
        if t[0] == key:
            return t[1], t[2]
    for k, label, example in TEMPLATES:
        if k == key:
            return label, (example or "one detector, chosen from the kinematics")
    raise KeyError(key)


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
    elif reaction.get("type") == "fusion-evaporation":
        chans = reaction.get("channels", [])
        state += f" · fusion-evaporation, {len(chans)} exit channel{'s' * (len(chans) != 1)}"
    else:
        state += " · elastic scattering"

    def consequence():
        el = p.energy_loss()["layers"]
        lay = el[0]
        parts = [f"the beam loses {guide._energy(lay['loss_mev'])} in the target "
                 f"(straggling {guide._energy(lay['straggling_fwhm_mev'])} FWHM)"]
        if len(el) > 1:
            parts.append(f"{guide._energy(sum(x['loss_mev'] for x in el[1:]))} more in the backing")
        if p.experiment.fusion is not None:
            f = p.fusion()
            lo, hi = f["excitation_range_mev"]
            parts.append(f"σ_fus {guide._fmt(f['sigma_thick_mb'])} mb, E* {lo:.1f}–{hi:.1f} MeV")
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
        if p.experiment.fusion is not None:
            study = p.plan()["study"]
            if study and study["line"] and study["beam_time_s"] is not None:
                enough = study["beam_time_s"] <= p.rates()["beam_time_s"]
                return (f"{study['channel']} → {study['residue']} needs {guide._duration(study['beam_time_s'])} for "
                        f"{r.get('counts_wanted')} counts in its {study['line']['energy_kev']:.1f} keV line: the beam "
                        f"time is {'enough' if enough else 'too short'}")
            return "set the γ-ray detectors and the counts wanted for the beam time"
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


def _ancillary(p, name: str) -> Decision:
    """The card of an ancillary detector: what it is, and what physim makes of it."""
    title, fields, more, _ = ANCILLARY_CARDS[name]
    a = p.draft[name]

    def state() -> str:
        if name == "neutron_wall":
            return (f"{a.get('cells', '?')} {a.get('material', 'EJ-301')} cells at {a.get('distance', '?')}"
                    + (f", from {a['theta_min']}" if a.get("theta_min") else ""))
        if name == "charged_particle_ball":
            return f"{a.get('elements', '?')} {a.get('material', 'CsI')} elements at {a.get('radius', '?')}"
        if name == "recoil_separator":
            return f"{a.get('kind', 'gas-filled')}, ±{a.get('theta_acceptance', '?')} acceptance"
        return f"{a.get('elements', '?')} {a.get('material', 'BGO')} elements, fold ≥ {a.get('fold_min', '?')}"

    def consequence() -> str:
        exp = p.experiment
        if name == "neutron_wall":
            m = exp.neutron_wall.model()
            return (f"{100 * m.geometric_efficiency():.1f} % of 4π, out to θ {m.theta_max():.0f}°; the channel "
                    "selection on the Plan uses it")
        if name == "recoil_separator":
            acc = ancillary_numbers(p).get("recoil_separator", {})
            if acc:
                return (f"{100 * acc['all']:.0f} % of the residues fly inside its acceptance (by angle only; its "
                        "transmission is not computed yet)")
        return "drawn in the scene; the channel selection does not use it yet"

    return Decision(f"anc:{name}", "Ancillary detectors", title, _try(state), _try(consequence), name, fields, more,
                    a)


def parse_distances(text: str) -> list:
    """Plunger distances as typed ("10, 20, 50 um", "10 um; 0.02 mm", "10 20 50") into the setup's list of texts
    with units: a bare number takes the unit of the next number that has one, µm when none has. Text that is not a
    distance is kept as it is, for the setup's checks to report."""
    import re

    tokens = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\s*[a-zA-Zµ]*|[^\s,;]+", str(text or ""))
    out, waiting = [], []
    for tok in (t.strip() for t in tokens):
        if not tok:
            continue
        m = re.fullmatch(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*([a-zA-Zµ]*)", tok)
        if m is None:
            out += [f"{w} um" for w in waiting] + [tok]
            waiting = []
            continue
        number, unit = m.groups()
        if not unit:
            waiting.append(number)
            continue
        unit = "um" if unit in ("µm", "μm") else unit
        out += [f"{w} {unit}" for w in waiting] + [f"{number} {unit}"]
        waiting = []
    return out + [f"{w} um" for w in waiting]


def lifetime_text(text: str) -> bool:
    """Whether a check (a problem or warning) is about the lifetime measurement, for the lifetime card."""
    low = text.lower()
    return any(k in low for k in ("lifetime", "plunger", "stopper", "degrader", "deorientation", "dsam"))


def _tau_text(tau_ps) -> str:
    """A lifetime in its natural unit: fs below 1 ps, ns from 1000 ps."""
    if tau_ps is None:
        return "τ unknown"
    if tau_ps < 1:
        return f"τ {1e3 * tau_ps:.3g} fs"
    return f"τ {tau_ps / 1e3:.3g} ns" if tau_ps >= 1e3 else f"τ {tau_ps:.3g} ps"


def _tau_pm(tau_ps: float, err_ps=None) -> str:
    """A measured lifetime with its uncertainty, in the lifetime's natural unit: "7.6 ± 0.31 ps"."""
    unit, f = ("fs", 1e3) if tau_ps < 1 else (("ns", 1e-3) if tau_ps >= 1e3 else ("ps", 1.0))
    return f"{tau_ps * f:.3g}" + (f" ± {err_ps * f:.2g}" if err_ps else "") + f" {unit}"


def lifetime_plan(planner) -> dict:
    """The lifetime plan of a planner (:meth:`physim.nuclear.planner.Planner.lifetimes`), or ``{"available":
    False, "reason": ...}`` when this physim has none yet or it fails: the app shows the reason, never an error."""
    fn = getattr(planner, "lifetimes", None)
    if fn is None:
        return {"available": False, "missing": True,
                "reason": "The lifetime plan (Planner.lifetimes) is not in this version of physim yet."}
    cache = getattr(planner, "_cache", None)  # the planner's own: it goes when the setup changes
    if isinstance(cache, dict) and ("workbench", "lifetimes") in cache:
        return cache[("workbench", "lifetimes")]
    try:
        out = fn()
    except Exception as err:  # noqa: BLE001 - the card and the panel say why
        out = {"available": False, "reason": f"The lifetime plan could not be worked out: {err}"}
    if isinstance(cache, dict):
        cache[("workbench", "lifetimes")] = out
    return out


def _recoil_line(p) -> str:
    """The recoil as the setup alone gives it (no lifetime plan): what is measured and how fast it flies."""
    import math

    from . import data
    from .names import parse_nuclide

    rec = p.experiment.lifetime_recoil()
    if rec is None:
        return "the emitter does not fly forward out of the target: no recoil to measure"
    nuclide, made, out = rec
    mass = 931.494 * data.nuclide(parse_nuclide(nuclide)).atomic_mass_u
    gamma = 1 + made / mass
    beta = math.sqrt(1 - 1 / gamma ** 2)
    left = "stops in the target and backing" if out <= 0 else f"leaves the target with {out:.3g} MeV"
    return f"{nuclide} made at {made:.3g} MeV (β {beta:.3f}), {left}"


def lifetime_line(planner) -> str:
    """The lifetime card's consequence line: the setup's own lifetime problem if it has one, else the studied
    level's lifetime and what the plan makes of it (peaks resolved, distances, F(τ), beam time per distance)."""
    p = planner
    problems = [t for t in p.problems if t.startswith("lifetime")]
    if problems:
        return problems[0]
    lt = lifetime_plan(p)
    if not lt.get("available"):
        if lt.get("missing"):
            return _recoil_line(p) + " (the plan per level is not in this version yet)"
        return lt.get("reason") or "no lifetime plan"
    levels = lt.get("levels") or []
    if not levels:
        return "no level with a lifetime to measure"
    lev = levels[0]
    parts = [f"{lev.get('level', '?')} → {lev.get('gamma_kev', 0):.0f} keV: {_tau_text(lev.get('tau_ps'))}"
             + (f" ({lev['tau_source']})" if lev.get("tau_source") else "")]
    if lev.get("method") == "dsam":
        if lev.get("attenuation") is not None:
            parts.append(f"F(τ) {lev['attenuation']:.2f}")
    elif lev.get("method") == "rdds":
        sep = lev.get("separation") or []
        if sep:
            parts.append(f"peaks apart at {sum(bool(s.get('resolved')) for s in sep)} of {len(sep)} detectors")
        rec = lev.get("recommended_distances_um") or []
        if rec:
            parts.append(f"distances {min(rec):.3g}–{max(rec):.3g} µm recommended")
    else:
        parts.append("outside what a plunger or DSAM measures")
    if lev.get("hours_per_distance") is not None:
        parts.append(f"{lev['hours_per_distance']:.3g} h per distance")
    n_warn = len(lt.get("warnings") or []) + len(lt.get("problems") or [])
    if n_warn:
        parts.append(f"{n_warn} note{'s' * (n_warn != 1)} in the plan")
    return " · ".join(parts)


def _lifetime(p) -> Decision:
    """The card of the lifetime measurement: a plunger's distances and foils, or DSAM's backing."""
    lt = p.draft["lifetime"]
    if lt.get("method") == "dsam":
        b = p.draft["target"].get("backing")
        title = "Line shapes (DSAM)"
        state = ("the recoil slows down and stops in the " + (f"{b.get('thickness', '')} {b.get('material', '')} "
                                                             "backing".replace("  ", " ") if b else "target itself"))
    else:
        title = "Plunger (RDDS)"
        ds = lt.get("distances") or []
        state = f"{len(ds)} distance{'s' * (len(ds) != 1)}" + (f", {ds[0]} to {ds[-1]}" if ds else "")
        if lt.get("offset"):
            state += f" (+ {lt['offset']} at contact)"
        st = lt.get("stopper") or {}
        foil = f"{st.get('thickness', '?')} {st.get('material', '?')} {st.get('kind', 'stopper')}"
        d2 = lt.get("degrader2")
        if st.get("kind") == "differential" and d2:
            foil += f", then {d2.get('thickness', '?')} {d2.get('material', '?')} {d2.get('gap', '?')} behind"
        state += " · " + foil
        if lt.get("deorientation"):
            state += " · deorientation set"
    return Decision("lifetime", "Lifetime measurement", title, state, _try(lambda: lifetime_line(p)), "lifetime",
                    [], [], lt)


def ancillary_numbers(planner) -> dict:
    """What can be said of the ancillaries from the simulated residues without their own physics: for a recoil
    separator, the share of all residues and of each channel's inside its angular acceptance ({"all": share,
    "channels": {label: share}}). Empty without a recoil separator."""
    import numpy as np

    exp = planner.experiment
    if exp.fusion is None or exp.recoil_separator is None:
        return {}
    from .experiment import _q

    half = _q(exp.recoil_separator.theta_acceptance).to("deg")
    ev = planner.residue_events()
    inside = np.asarray(ev["theta"]) <= half
    weight = np.asarray(ev["weight"])
    channel = np.asarray(ev["channel"])
    total = weight.sum()
    per = {}
    for k, label in enumerate(ev.channels):
        w = weight[channel == k]
        per[label] = float(w[inside[channel == k]].sum() / w.sum()) if w.sum() > 0 else 0.0
    return {"recoil_separator": {"all": float(weight[inside].sum() / total) if total > 0 else 0.0,
                                 "channels": per, "half_angle_deg": half}}


def decisions(planner) -> list:
    """The setup panel's cards, in the order the decisions are made (see :data:`GROUPS`)."""
    p = planner
    rows = _try(lambda: {r["detector"]: r for r in p.rates()["rows"]}, {})
    out = [_beam(p), _target(p)]
    out += [_detector(p, i, rows) for i in range(len(p.draft.get("detectors", [])))]
    out += [_gamma_detector(p, i) for i in range(len(p.draft.get("gamma_detectors", [])))]
    if p.draft.get("reaction", {}).get("type") == "fusion-evaporation":
        out += [_ancillary(p, name) for name in ANCILLARY_CARDS if name in p.draft]
    if "lifetime" in p.draft and p.draft.get("reaction", {}).get("type") in ("coulex", "fusion-evaporation"):
        out.append(_lifetime(p))
    out.append(_run(p))
    return out


class _DraftOnly:
    """A planner's draft and nothing else: its cards say what was chosen, and no consequence is worked out."""

    def __init__(self, planner):
        self.draft = planner.draft

    def __getattr__(self, name):
        raise AttributeError(name)


def draft_decisions(planner) -> list:
    """The cards of the planner's draft with their first line (what was chosen) and their fields, and "—" for the
    consequence: nothing is computed. For a draft whose numbers are not there (yet)."""
    return decisions(_DraftOnly(planner))


def _named(items: list, word: str) -> dict:
    return {(x.get("name") or f"{word} {i + 1}"): x for i, x in enumerate(items)}


def card_setup(setup: dict, key: str, name: Optional[str] = None):
    """The part of a setup (in the structure of a setup file) that the card ``key`` edits: "beam", "target" (with
    the reaction), "run", or a detector, taken by ``name`` if given and by its place in the list ("detector:N",
    "gamma:N") otherwise. None if the setup has no such detector."""
    if key == "target":
        return {"target": setup.get("target"), "reaction": setup.get("reaction", {"type": "elastic"})}
    if key in ("beam", "run"):
        return setup.get(key)
    if key == "lifetime":  # DSAM reads the backing, so a change of the backing changes this card too
        return {"lifetime": setup.get("lifetime"), "backing": (setup.get("target") or {}).get("backing")}
    if key.startswith("anc:"):
        return setup.get(key.split(":", 1)[1])
    kind, i = key.split(":")
    items = setup.get("detectors" if kind == "detector" else "gamma_detectors", [])
    if name is not None:
        return _named(items, "D" if kind == "detector" else "γ").get(name)
    return items[int(i)] if int(i) < len(items) else None


def changes(applied: dict, draft: dict) -> list:
    """What differs between two setups (in the structure of a setup file), one line each: a changed field with
    its old and new value ("Beam: energy 30 MeV → 31 MeV"), a detector added or removed. Detectors are matched by
    name, so one that moved up the list after another was removed is not a change."""
    out = []

    def show(x) -> str:
        return "—" if x is None else str(x)

    def fields(title: str, a: Optional[dict], b: Optional[dict], skip: tuple = ()) -> None:
        a, b = a or {}, b or {}
        for k in dict.fromkeys([*a, *b]):
            if k in skip or a.get(k) == b.get(k):
                continue
            label = k.replace("_", " ")
            if isinstance(a.get(k), (dict, list)) or isinstance(b.get(k), (dict, list)):
                out.append(f"{title}: {label}")
            else:
                out.append(f"{title}: {label} {show(a.get(k))} → {show(b.get(k))}")

    for k in ("title", "description"):
        if applied.get(k) != draft.get(k):
            out.append(f"{k.capitalize()}: {show(draft.get(k))}")
    fields("Beam", applied.get("beam"), draft.get("beam"))
    fields("Target", applied.get("target"), draft.get("target"), skip=("backing",))
    fields("Backing", (applied.get("target") or {}).get("backing"), (draft.get("target") or {}).get("backing"))
    fields("Reaction", applied.get("reaction"), draft.get("reaction"))
    for kind, word in (("detectors", "D"), ("gamma_detectors", "γ")):
        a, b = _named(applied.get(kind, []), word), _named(draft.get(kind, []), word)
        out += [f"{name}: added" for name in b if name not in a]
        out += [f"{name}: removed" for name in a if name not in b]
        for name in b:
            if name in a:
                fields(name, a[name], b[name])
    fields("Run conditions", applied.get("run"), draft.get("run"))
    a_lt, b_lt = applied.get("lifetime"), draft.get("lifetime")
    if a_lt is None and b_lt is not None:
        out.append("Lifetime measurement: added")
    elif a_lt is not None and b_lt is None:
        out.append("Lifetime measurement: removed")
    elif a_lt != b_lt:
        parts = ("stopper", "degrader2", "deorientation")
        fields("Lifetime", {k: v for k, v in a_lt.items() if k not in parts},
               {k: v for k, v in b_lt.items() if k not in parts})
        for part, title in zip(parts, ("Plunger foil", "Second foil", "Deorientation")):
            if a_lt.get(part) is None and b_lt.get(part) is not None:
                out.append(f"{title}: added")
            elif a_lt.get(part) is not None and b_lt.get(part) is None:
                out.append(f"{title}: removed")
            else:
                fields(title, a_lt.get(part), b_lt.get(part))
    known = ("schema", "title", "description", "beam", "target", "reaction", "detectors", "gamma_detectors", "run",
             "lifetime")
    for k in dict.fromkeys([*applied, *draft]):
        if k not in known and applied.get(k) != draft.get(k):
            out.append("Level schemes" if k == "levels" else k.replace("_", " ").capitalize())
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
    n_p, n_g = len(d.get("detectors", [])), len(d.get("gamma_detectors", []))
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


# -- the steps ------------------------------------------------------------------------------------------------------

def _span(seconds: float) -> str:
    if seconds >= 86400:
        return f"{seconds / 86400:.3g} d"
    if seconds >= 3600:
        return f"{seconds / 3600:.3g} h"
    return f"{seconds / 60:.3g} min" if seconds >= 60 else f"{seconds:.3g} s"


def _plan_line(p) -> tuple:
    """(the line under Plan, done, warn) from the rates: the beam time the counts wanted need against the beam
    time planned, or the highest rate when no counts are asked for. For fusion-evaporation, the studied
    channel's beam time."""
    r = p.rates()
    if p.experiment.fusion is not None:
        study = p.plan().get("study")
        needed = study["beam_time_s"] if study else None
        if needed is None:
            return "set the counts wanted", False, False
        if needed != needed or needed == float("inf"):
            return f"{study['channel']}: never reached", False, True
        enough = needed <= r["beam_time_s"]
        return (f"{study['channel']}: {_span(needed)} needed" + ("" if enough else f", {_span(r['beam_time_s'])} "
                                                                  "planned"), enough, not enough)
    rows = r["rows"]
    if not rows:
        return "", False, False
    timed = [x for x in rows if x.get("beam_time_s") is not None]
    if r.get("counts_wanted") and timed:
        needed = max(x["beam_time_s"] for x in timed)
        if needed != needed or needed == float("inf"):
            return "the counts are never reached", False, True
        enough = needed <= r["beam_time_s"]
        return f"{_span(needed)} needed" + ("" if enough else f", {_span(r['beam_time_s'])} planned"), enough, not enough
    fast = max(x["rate_per_s"] for x in rows)
    return f"up to {_rate(fast)}", True, False


def steps(planner) -> list:
    """The step bar, left to right: for each of :data:`STEPS` {"key", "label", "number", "state", "status",
    "warn"}.

    ``state`` is "done" (the step has its result: a setup without errors, a beam time that is enough, a run, a
    run's data, an analysis of the current run), "ready" (it can be worked on) or "waiting" (it needs a run
    first). ``status`` is the line under the step's name, and ``warn`` marks a status that needs attention
    (warnings in the setup, a beam time too short, a run taken with another setup). It reads what the page keeps
    (the rates, the runs' summaries, the analysis kept with the run) and simulates nothing."""
    p = planner
    d = p.draft
    n_p, n_g = len(d.get("detectors", [])), len(d.get("gamma_detectors", []))
    ws = _try(p.warnings, [])
    errors = sum(w.level == "error" for w in ws) + len(p.problems)
    warnings = sum(w.level == "warning" for w in ws)
    setup = f"{n_p} particle" + (f" + {n_g} γ" if n_g else "")
    if d.get("reaction", {}).get("type") == "fusion-evaporation":
        short = {"neutron_wall": "wall", "charged_particle_ball": "ball", "recoil_separator": "RS",
                 "gamma_calorimeter": "cal"}
        setup = " · ".join([f"{n_g} γ"] + [v for k, v in short.items() if k in d] + (["monitor"] if n_p else []))
    if errors:
        setup += f" · {errors} error{'s' * (errors != 1)}"
    elif warnings:
        setup += f" · {warnings} warning{'s' * (warnings != 1)}"
    plan, plan_done, plan_warn = _try(lambda: _plan_line(p), ("", False, False))

    runs = _try(p.runs, [])
    run = p.run
    stale = bool(_try(p.run_status, None) and p.run_status()["stale"])
    if p.running:
        taken, run_state = "a run is being taken", "ready"
    elif runs:
        last = run.describe().split(";")[0] if run is not None else runs[-1]["describe"].split(";")[0]
        taken = f"{len(runs)} run{'s' * (len(runs) != 1)} · " + last[0].lower() + last[1:]
        if stale:
            taken = f"run {run.number}: setup changed since"
        run_state = "done"
    else:
        taken, run_state = "no run yet", "ready"

    if run is None:
        data, data_state = "needs a run", "waiting"
    else:
        gates = len(_try(p.gates, []))
        data = f"run {run.number}" + (f" · {gates} gate{'s' * (gates != 1)}" if gates and run.kind != "source" else "")
        data_state = "done"

    result = _try(p.run_analysis, None) if run is not None and run.kind != "fusion" else None
    coulex = p.experiment.excitation is not None and bool(p.experiment.gamma_detectors)
    lifetime = p.experiment.lifetime
    kept = _try(lambda: kept_lifetimes(p), None) if lifetime is not None else None
    gap = not_analysed(p.experiment)
    if kept is not None:  # a lifetime experiment's result is τ, as a Coulex one's is B(E2)
        later = len(_try(lambda: runs_not_read(p, kept, runs), []))
        analysis = f"τ({analysed_level(kept)}) {_tau_pm(kept['tau_ps'], kept.get('stat_ps'))}" + (
            f" · {later} run{'s' * (later != 1)} since" if later else "")
        analysis_state = "ready" if later else "done"
    elif lifetime is not None and not (result is not None and getattr(result, "b_e2fm4", None) is not None):
        if lifetime.method == "rdds":
            have = {r.get("distance_um") for r in runs if not r.get("stale") and r.get("distance_um") is not None}
            n = len(lifetime.distances_um)
            analysis = (f"lifetimes: {len(have)} of {n} distance{'s' * (n != 1)} taken" if have else
                        "lifetimes: take a run at each distance")
        else:
            analysis = "lifetimes: line shapes of run " + str(run.number) if run is not None else \
                "lifetimes: take a run for the line shapes"
        analysis_state = "ready"
    elif gap is not None and not _coulex_run(p):  # no analysis for this kind yet: the step waits, and says for what
        analysis, analysis_state = "no analysis yet", "waiting"
    elif not coulex and not _coulex_run(p):
        analysis, analysis_state = "needs γ-ray detectors", "waiting"
    elif result is not None and getattr(result, "b_e2fm4", None) is not None:
        analysis, analysis_state = f"run {run.number}: B(E2↑) {_num(result.b_e2fm4, 4)} e²fm⁴", "done"
    elif result is not None:
        analysis, analysis_state = f"run {run.number} analysed", "done"
    elif run is None:
        analysis, analysis_state = "no run: on a sample", "ready"
    elif run.kind == "source":
        analysis, analysis_state = "the current run is a source run", "ready"
    else:
        analysis, analysis_state = f"run {run.number} not analysed yet", "ready"

    report = f"the plan and {len(runs)} run{'s' * (len(runs) != 1)}" if runs else "the plan only"
    out = {"setup": ("ready" if errors else "done", setup, bool(errors or warnings)),
           "plan": ("done" if plan_done and not errors else "ready", plan, plan_warn),
           "run": (run_state, taken, stale),
           "data": (data_state, data, False),
           "analysis": (analysis_state, analysis, False),
           "report": ("ready", report, False)}
    return [{"key": key, "label": label, "number": i, "state": out[key][0], "status": out[key][1],
             "warn": out[key][2]} for i, (key, label) in enumerate(STEPS, 1)]


def _coulex_run(planner) -> bool:
    """Whether the current run is a beam run of a Coulomb excitation with γ-ray detectors, whatever the setup has
    become since: it is analysed."""
    run = planner.run
    if run is None or run.kind == "source":
        return False
    try:
        return run.experiment.excitation is not None and bool(run.experiment.gamma_detectors)
    except Exception:  # noqa: BLE001 -- a setup that cannot be read
        return False


def next_step(planner, key: str, pending: int = 0) -> Optional[dict]:
    """What the foot of the step ``key`` says and where its button leads: {"hint", "label", "to"} (``label`` and
    ``to`` None when there is nothing to go on to yet), or None for a stage that is not a step. ``pending`` is the
    number of changes of the setup that wait for Apply."""
    p = planner
    run = p.run
    # Fusion-evaporation and elastic scattering have no analysis of a run yet: their steps lead from the data to
    # the report, and say what an analysis would need.
    gap = None if _coulex_run(p) else not_analysed(p.experiment)
    if key == "setup":
        if p.problems:
            return {"hint": "The setup has a problem: the checks at the top of the setup panel say what.",
                    "label": None, "to": None}
        hint = (f"{pending} change{'s' * (pending != 1)} wait{'s' * (pending == 1)} for Apply: the plan is for the "
                "setup as it was applied." if pending else "The setup is applied.")
        return {"hint": hint, "label": "Next: check the plan", "to": "plan"}
    if key == "plan":
        return {"hint": "When the rates and the beam time are what you want, take a run.",
                "label": "Next: take a run", "to": "run"}
    if key == "run":
        if p.running:
            return {"hint": "The run is being taken; its data are there when it ends.", "label": None, "to": None}
        plunger = p.experiment.lifetime is not None and p.experiment.lifetime.method == "rdds"
        if run is None:
            return {"hint": "No run yet: choose its length " + ("and the plunger distance " if plunger else "")
                            + "above and press Run.", "label": None, "to": None}
        more = " A plunger measurement takes one run per distance: the Data and Analysis steps read them together." \
            if plunger else ""
        return {"hint": f"Run {run.number} is the current run.{more}",
                "label": f"Next: look at run {run.number}'s data", "to": "data"}
    if key == "data":
        if run is None:
            return {"hint": "The spectra are a run's.", "label": "Take a run", "to": "run"}
        if run.kind == "source":
            return {"hint": f"Run {run.number} is a source run: it measures the γ-ray efficiency. The analysis "
                            "reads a beam run.", "label": "Take a beam run", "to": "run"}
        if gap is not None:
            return {"hint": f"{gap['line'][0].upper()}{gap['line'][1:]} The report holds {gap['has']}, and run "
                            f"{run.number} beside the plan.", "label": "Next: the report", "to": "report"}
        if p.experiment.lifetime is None and not p.experiment.gamma_detectors and not _coulex_run(p):
            return {"hint": "The analysis reads the γ rays: this setup has no γ-ray detector to analyse a run with.",
                    "label": "Next: the report", "to": "report"}
        return {"hint": "A gate saved here can be used by the analysis.", "label": f"Next: analyse run {run.number}",
                "to": "analysis"}
    if key == "analysis":
        if gap is not None:
            return {"hint": f"The report holds {gap['has']}, and the runs beside the plan.",
                    "label": "Next: the report", "to": "report"}
        return {"hint": "The analysis of a run is kept with it, and the report holds it.",
                "label": "Next: the report", "to": "report"}
    if key == "report":
        if gap is not None:
            return {"hint": f"The report is the last step: it holds the setup, {gap['has']}, and the runs beside "
                            "the plan.", "label": None, "to": None}
        return {"hint": "The report is the last step: it holds the setup, the plan, the runs and their analysis.",
                "label": None, "to": None}
    return None


# -- new experiments ------------------------------------------------------------------------------------------------

def template(key: str):
    """The :class:`~physim.nuclear.Experiment` of a template (:data:`TEMPLATES`); "blank" needs a beam and target,
    see :func:`new_setup`."""
    from .experiment import Experiment

    fusion = next((t for t in FUSION_TEMPLATES if t[0] == key), None)
    if fusion is not None:
        return new_setup("12C", "60 MeV", "150Nd", "1 mg/cm2", "fusion-evaporation", detectors_from=key)
    lifetime = next((t for t in LIFETIME_TEMPLATES if t[0] == key), None)
    if lifetime is not None:  # the worked example itself: its beam and target suit its detectors and foils
        return Experiment.example(lifetime[3])
    found = next((t for t in TEMPLATES if t[0] == key), None)
    if found is None:
        known = [t[0] for t in TEMPLATES + FUSION_TEMPLATES + LIFETIME_TEMPLATES]
        raise KeyError(f"no template {key!r}; the templates are {known}")
    if found[2] is None:
        return new_setup("4He", "5.5 MeV", "Au")
    return Experiment.example(found[2])


def new_setup(beam: str, energy: str, target: str, thickness: str = "0.5 mg/cm2", measure: str = "elastic",
              current: str = "1 pnA", beam_time: str = "24 h", title: Optional[str] = None,
              detectors_from: Optional[str] = None):
    """A new experiment's setup from what the form of a kind of experiment on the start page asks (the app's
    ``create_experiment``): the beam and its energy, the target,
    and what to measure (:data:`MEASUREMENTS`). It gets one particle detector where the kinematics send the
    particles (behind the target for a light beam on a heavy target, forward otherwise); more come from the Add
    menus. ``detectors_from`` (a key of :data:`TEMPLATES`) takes a template's particle and γ-ray detectors
    instead; the beam, target and reaction stay the ones given here. For Coulomb excitation the first 2⁺ state is
    taken from the local ENSDF copy when there is one; without it the setup stays elastic and
    :func:`new_setup_notes` says what to enter."""
    from . import data
    from .experiment import SCHEMA, Experiment
    from .names import parse_nuclide
    from .planner import first_state

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
    if detectors_from and detectors_from != "blank" and measure != "fusion-evaporation":
        source = template(detectors_from).to_dict()
        d["detectors"] = source["detectors"]
        if source.get("gamma_detectors"):
            d["gamma_detectors"] = source["gamma_detectors"]
    if measure == "fusion-evaporation":
        return _with_lifetime(_new_fusion_setup(d, detectors_from), detectors_from)
    exp = Experiment.from_dict(copy.deepcopy(d))
    notes = []
    if measure in ("coulex-target", "coulex-beam"):
        role = "target" if measure == "coulex-target" else "beam"
        nuclide = (heavy.name if heavy else target) if role == "target" else beam_name
        try:
            scheme, n, b = first_state(nuclide, "E2")
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
    return _with_lifetime(exp, detectors_from)


def _with_lifetime(exp, key: Optional[str]):
    """A new setup with the lifetime measurement of the lifetime template ``key`` (:data:`LIFETIME_TEMPLATES`;
    any other key leaves it as it is): a plunger with the default distances and a stopper thick enough for this
    setup's recoil, or DSAM with the example's backing when the target has none."""
    from .experiment import Experiment
    from .planner import LIFETIME_REACTIONS, Planner

    found = next((t for t in LIFETIME_TEMPLATES if t[0] == key), None)
    if found is None:
        return exp
    notes = new_setup_notes(exp)
    if exp.reaction not in LIFETIME_REACTIONS:
        notes.append("The lifetime measurement waits for the reaction: set Coulomb excitation in the Target and "
                     "reaction card, then add the plunger or DSAM in the Lifetime measurement card.")
        exp._new_setup_notes = notes
        return exp
    p = Planner(exp)
    if found[5] == "dsam":
        if exp.target.backing is None:
            backing = Experiment.example(found[3]).to_dict()["target"].get("backing") or {"material": "Au",
                                                                                         "thickness": "15 mg/cm2"}
            p.set("backing", "material", backing["material"])
            p.set("backing", "thickness", backing["thickness"])
        p.set_lifetime(method="dsam")
    else:
        p.add_ancillary("plunger")
    if p.problems:
        notes.append("The lifetime measurement of the template does not fit this beam and target as it is: "
                     + "; ".join(p.problems) + ". Change it in the Lifetime measurement card.")
    try:
        out = Experiment.from_dict(copy.deepcopy(p.draft))
    except Exception:  # noqa: BLE001 - keep the setup as made, with the note above
        out = exp
    out._new_setup_notes = notes
    return out


def _new_fusion_setup(d: dict, detectors_from: Optional[str]):
    """A new fusion-evaporation setup: the channel the excitation energy favours as a placeholder (Fill from
    physim in the reaction's card replaces it), and the detectors of a fusion template (:data:`FUSION_TEMPLATES`):
    its γ-ray array and beam monitor, and its ancillaries."""
    from .experiment import Experiment
    from .planner import Planner

    found = next((t for t in FUSION_TEMPLATES if t[0] == detectors_from), None)
    lifetime = next((t for t in LIFETIME_TEMPLATES if t[0] == detectors_from), None)
    if lifetime is not None:  # the γ-ray array (and monitor) of the lifetime template's example, no ancillary
        found = (lifetime[0], lifetime[1], lifetime[2], lifetime[3], ())
    if found is None:
        found = FUSION_TEMPLATES[-1] if detectors_from in (None, "blank", "fusion-blank") else FUSION_TEMPLATES[0]
    _, _, _, example, ancillaries = found
    d.pop("detectors", None)
    if example is not None:
        source = Experiment.example(example).to_dict()
        if source.get("detectors"):
            d["detectors"] = source["detectors"]
        if source.get("gamma_detectors"):
            d["gamma_detectors"] = source["gamma_detectors"]
    d["run"]["counts_wanted"] = 10000
    if not d.get("detectors"):
        # A placeholder monitor while the reaction is still elastic; it goes once the setup is fusion-evaporation.
        d["detectors"] = [guide.placement("forward")]
        drop_monitor = True
    else:
        drop_monitor = False
    p = Planner(Experiment.from_dict(copy.deepcopy(d)))
    p.set("reaction", "type", "fusion-evaporation")
    draft = p.draft
    if drop_monitor:
        draft.pop("detectors", None)
    for name in ancillaries:
        draft[name] = dict(Planner.ANCILLARY_DEFAULTS[name])
    exp = Experiment.from_dict(copy.deepcopy(draft))
    exp._new_setup_notes = []
    return exp


def new_setup_notes(experiment) -> list:
    """What :func:`new_setup` could not fill in."""
    return list(getattr(experiment, "_new_setup_notes", []))


# -- picking nuclides and states ------------------------------------------------------------------------------------

#: The choice of a target made of something other than one element (a compound, or a named material).
COMPOUND = "compound"
#: The mass-number choice of a target of the element in its natural isotopic mix.
NATURAL = "natural"


def elements() -> dict:
    """The elements the beam and target pickers offer, by Z: {symbol: "Ni · nickel (28)"}, those with at least one
    measured mass."""
    from . import data
    from .names import SYMBOLS

    measured = {n.Z for n in data.nuclides() if not n.estimated}
    out = {}
    for z, symbol in enumerate(SYMBOLS):
        if z and z in measured:
            try:
                name = (data.element(symbol).name or "").lower()
            except (KeyError, ValueError):
                name = ""
            out[symbol] = f"{symbol} · {name} ({z})" if name else f"{symbol} ({z})"
    return out


def mass_numbers(symbol: str) -> dict:
    """An element's nuclides with a measured mass, {A: label}: the stable ones with their natural abundance
    ("62 · 3.63 %"), the others as their mass number."""
    from . import data
    from .names import Z_OF

    z = Z_OF[symbol]
    try:
        abundance = dict(data.element(symbol).isotopes)
    except (KeyError, ValueError):
        abundance = {}
    out = {}
    for n in sorted((n for n in data.nuclides() if n.Z == z and not n.estimated), key=lambda n: n.A):
        share = abundance.get(n.A)
        out[n.A] = f"{n.A} · {100 * share:.3g} %" if share else str(n.A)
    return out


def most_abundant(symbol: str) -> Optional[int]:
    """The mass number of an element's most abundant isotope, or None for an element with no stable one."""
    from . import data

    try:
        iso = data.element(symbol).isotopes
    except (KeyError, ValueError):
        return None
    return max(iso, key=lambda x: x[1])[0] if iso else None


def split_nuclide(text: Optional[str]) -> Optional[tuple]:
    """A beam nuclide or target material as (symbol, A), with A None for a natural element ("Pt"); None when it is
    neither (a compound such as "CD2", or nothing)."""
    from .names import Z_OF, parse_nuclide

    text = (text or "").strip()
    if not text:
        return None
    if text in Z_OF:
        return (text, None)
    try:
        z, a = parse_nuclide(text)
    except (KeyError, ValueError):
        return None
    from .names import SYMBOLS

    return (SYMBOLS[z], a)


def states(scheme, multipolarity: str = "E2", limit: int = 30) -> dict:
    """The states of a level scheme that ``multipolarity`` reaches from the ground state with a known B(Eλ↑), for
    the state picker: {level index: "2+ · 1172.9 keV · B(E2↑) 884 e²fm⁴"}, lowest first."""
    if scheme is None:
        return {}
    lam = int(multipolarity[1])
    out = {}
    for n, lev in enumerate(scheme.levels):
        b = scheme.b(0, n, multipolarity) if n else None
        if not b:
            continue
        value = f"{b:,.0f}" if b >= 100 else f"{b:.3g}"
        out[n] = f"{lev.label} · {lev.energy.value:.1f} keV · B({multipolarity}↑) {value} e²fm{_SUPERSCRIPT[2 * lam]}"
        if len(out) >= limit:
            break
    return out


_SUPERSCRIPT = {2: "²", 4: "⁴", 6: "⁶", 8: "⁸"}


__all__ = ["ANCILLARY_CARDS", "ancillary_numbers", "card_setup", "changes", "draft_decisions", "COMPOUND",
           "FUSION_TEMPLATES", "GROUPS", "KINDS", "Kind", "LIFETIME_METHODS", "LIFETIME_TEMPLATES", "FOIL_KINDS",
           "lifetime_line", "lifetime_plan",
           "kind_of", "template_label", "MEASUREMENTS", "NATURAL", "NOT_SIMULATED", "PHYSICS", "PHYSICS_FOR", "physics_for",
           "STAGES", "STEPS", "TEMPLATES", "Decision", "checks", "decisions", "elements", "mass_numbers",
           "most_abundant", "new_setup", "new_setup_notes", "next_step", "split_nuclide", "states", "status_strip",
           "steps", "template"]
