"""The experiment planner app's model: everything the web app shows, without the web app.

The planner app (backlog item 41) is a thin layer over this module. It edits a setup, keeps the last valid
:class:`~physim.nuclear.Experiment`, and computes the content of every result tab as plain data (arrays, tables,
text) that any front end, or a notebook, can draw::

    from physim.nuclear.planner import Planner

    p = Planner.example("alpha_on_gold")
    p.set("beam", "energy", "6 MeV")          # edits are checked; problems go to p.problems
    p.add_detector(name="A75", shape="circle", theta="75 deg", distance="100 mm", radius="2.5 mm",
                   thickness="300 um")
    p.warnings()                              # validity and rate warnings, most important first
    p.rates()                                 # one row per detector
    p.explain("rates")                        # the formula, assumptions, limits, register link
    p.sweep("beam energy", ["4 MeV", "5 MeV", "6 MeV"], "rate", detector="A45")

Tabs: :data:`TABS`. Each ``Planner.<tab>()`` returns a dict; :meth:`Planner.tab` calls one by name.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np

from . import data
from . import data as _data
from . import ensdf as _ensdf
from . import response as _response
from . import scene as _scene
from .detectors import Array, Response
from .experiment import Experiment, SetupError, example_names
from .kinematics import TwoBody
from .levels import LevelScheme
from .names import parse_nuclide
from .quantity import Quantity
from .rates import Rates, beam_energy_at, beam_ion, stack, stopping, tilt_deg
from .rutherford import Rutherford

#: The result tabs, in display order.
TABS = ("geometry", "kinematics", "rates", "energy_loss", "spectra", "trajectories", "gamma", "report")

#: Where each tab's physics is documented and validated (paths under the docs root).
REGISTER = {
    "geometry": "physics-register/detectors",
    "kinematics": "physics-register/kinematics",
    "rates": "physics-register/rates-and-events",
    "energy_loss": "physics-register/stopping",
    "spectra": "physics-register/rates-and-events",
    "trajectories": "physics-register/rutherford",
    "gamma": "theory/coulex",
    "report": "physics-register/README",
}

#: The "Explain" panel of each tab: formula, assumptions, where it stops being valid.
EXPLAIN = {
    "geometry": {
        "title": "Detector geometry",
        "formula": "Ω = ∫ (r̂ · n̂) / r² dA over each face (Gauss–Legendre quadrature).",
        "assumptions": "Flat faces; the beam runs along +z through the target at the origin; point source unless a "
                       "beam spot is given.",
        "limits": "No gaps between strips, no inter-strip effects. A detector in the beam path is flagged.",
    },
    "kinematics": {
        "title": "Two-body kinematics",
        "formula": "Relativistic: E and θ in the lab from the CM angle by a Lorentz boost with β_cm = p₁/(E₁ + m₂).",
        "assumptions": "Nuclear masses from AME2020; the target is at rest.",
        "limits": "Exact for two-body reactions. Where the beam is heavier than the target there are two energies "
                  "at each angle, up to a maximum angle.",
    },
    "rates": {
        "title": "Count rates and beam time",
        "formula": "counts/s = (beam particles/s) × (nuclei/cm²) × ∫ dσ/dΩ dΩ; Rutherford "
                   "dσ/dΩ = (Z₁Z₂e²/4E)² / sin⁴(θ/2). Beam time for N counts = N / rate; relative error 1/√N.",
        "assumptions": "Elastic Rutherford scattering, averaged over the target depth; ejectile and recoil, target "
                       "and backing.",
        "limits": "Rutherford fails above the Coulomb barrier at large angles, at very low energy (screening) and "
                  "for identical nuclei (Mott). The warnings say when.",
    },
    "energy_loss": {
        "title": "Energy loss and straggling",
        "formula": "E_out = R⁻¹(R(E_in) − path), with R the CSDA range from tabulated (ICRU 49/90) and modelled "
                   "stopping powers; straggling from Bohr's formula carried through the layer.",
        "assumptions": "Paths grow as 1/cos of the angle to the layer normal; fully stripped ions at high velocity, "
                       "effective charge below.",
        "limits": "Slow heavy ions (below ~1 MeV/u) are uncertain to 10–50%; straggling below 10 MeV/u to ~2×.",
    },
    "spectra": {
        "title": "Simulated spectra",
        "formula": "Monte Carlo: reaction depth, beam spread and spot, straggling, Rutherford angle, kinematics, "
                   "exit energy loss and detector response, event by event; each event weighted so spectra are "
                   "counts in the planned beam time.",
        "assumptions": "Straight tracks; Gaussian straggling and resolution.",
        "limits": "No multiple scattering, beam divergence, background or pile-up.",
    },
    "trajectories": {
        "title": "Coulomb trajectories",
        "formula": "Classical orbits integrated by physim's engine in the CM frame; deflection 2 arctan(d₀/2b), "
                   "closest approach d₀/2 (1 + 1/sin(θ/2)).",
        "assumptions": "Point charges, classical mechanics (valid when the Sommerfeld parameter η ≫ 1).",
        "limits": "Inside the interaction radius nuclear forces act and the orbit is no longer Coulomb.",
    },
    "gamma": {
        "title": "Coulomb excitation and γ rays",
        "formula": "First-order semiclassical Coulomb excitation: P(θ) from the multipole field integrated along the "
                   "Rutherford orbit (Alder and Winther); dσ/dΩ = P dσ_R/dΩ. γ energy E₀ √(1 − β²)/(1 − β cos α) "
                   "for a nucleus moving at β, emitting at α to its velocity.",
        "assumptions": "One state reached from a 0⁺ ground state; B(Eλ) as given in the setup; decay in flight after "
                       "leaving the target, with the particle–γ angular correlation of first-order excitation "
                       "(or isotropic, by the switch).",
        "limits": "Closer than Cline's safe distance nuclear forces interfere; lifetimes and deorientation are not "
                  "modelled. The rates and the table of all levels are first order; 'Solve with all orders' "
                  "includes multi-step excitation and reorientation for the level scheme.",
    },
    "report": {
        "title": "Beam-time report",
        "formula": "Collects the setup, rates, beam time, peaks and warnings.",
        "assumptions": "As for each tab.",
        "limits": "Only what the tabs compute; check the warnings first.",
    },
}


@dataclass
class Warning:
    """A message for the warnings banner. ``level`` is "error" (the setup cannot be used), "warning" or "note"."""

    level: str
    source: str
    text: str


def _q(x) -> Quantity:
    return Quantity.parse(x) if isinstance(x, str) else x


class Planner:
    """An editable setup and the results computed from it.

    Edits go to a working copy of the setup file's structure. When the copy is valid it becomes :attr:`experiment`
    and the results are recomputed (lazily); otherwise :attr:`problems` lists what is wrong and the results keep
    showing the last valid setup.
    """

    def __init__(self, experiment: Experiment):
        experiment.validate()
        self._draft: dict = experiment.to_dict()
        self.experiment: Experiment = experiment
        #: Problems with the current draft (empty when it is valid).
        self.problems: list = []
        self._cache: dict = {}
        #: The rates before the last edit: detectors the edit left alone keep their rows (see :class:`Rates`).
        self._earlier_rates: Optional[Rates] = None
        #: The Coulomb-excitation settings while the reaction is switched to elastic, so switching back restores them.
        self._stashed_excitation: dict = {}

    # -- loading and saving -----------------------------------------------------------------------------------------

    @classmethod
    def example(cls, name: str) -> Planner:
        return cls(Experiment.example(name))

    @staticmethod
    def examples() -> list:
        return example_names()

    @classmethod
    def load(cls, path: Union[str, Path]) -> Planner:
        return cls(Experiment.load(path))

    def save(self, path: Union[str, Path]) -> None:
        """Write the current (valid) setup to a setup file."""
        if self.problems:
            raise SetupError(self.problems)
        self.experiment.save(path)

    def to_toml(self) -> str:
        return self.experiment.to_toml()

    # -- editing ----------------------------------------------------------------------------------------------------

    @property
    def draft(self) -> dict:
        """A copy of the setup being edited, in the structure of a setup file."""
        return copy.deepcopy(self._draft)

    def _apply(self, draft: dict) -> bool:
        self._draft = draft
        try:
            exp = Experiment.from_dict(copy.deepcopy(draft))
        except SetupError as e:
            self.problems = list(e.problems)
            return False
        self.problems = []
        self.experiment = exp
        self._earlier_rates = self._cache.get("rates", self._earlier_rates)
        self._cache.clear()
        return True

    def place(self, key: str, **fields) -> bool:
        """Set several fields of one detector at once; ``key`` is "detector:N" or "gamma:N", counting from 0 as
        the scene does (:mod:`physim.nuclear.scene`). ``None`` removes a field."""
        d = self.draft
        kind, i = key.split(":")
        entry = d["detectors" if kind == "detector" else "gamma_detectors"][int(i)]
        for k, v in fields.items():
            if v is None:
                entry.pop(k, None)
            else:
                entry[k] = v
        return self._apply(d)

    def move(self, key: str, position, mode: str = "angle") -> bool:
        """Move a detector's front face to ``position`` (x, y, z in mm), as a drag in the scene does: in ``mode``
        "angle" it keeps its distance and takes the new direction, in "distance" it keeps its direction. The
        result is the same as typing the new angles or distance."""
        return self.place(key, **_scene.move_fields(self.experiment, key, position, mode))

    def set(self, section: str, field: str, value: Any) -> bool:
        """Set ``field`` of ``section`` ("beam", "target", "backing", "run", "reaction", "detector N" or
        "gamma detector N", counting from 1); ``None`` removes an optional field. Returns whether the setup is valid
        afterwards.

        Setting the reaction ``type`` to ``"coulex"`` fills in E2 excitation of the target (the state's energy and
        B(E2↑) still have to be given); switching back to ``"elastic"`` keeps the excitation settings for later."""
        d = self.draft
        if section == "reaction" and field == "type":
            return self._set_reaction_type(d, value)
        if section == "backing":
            sec = d["target"].setdefault("backing", {})
        elif section.startswith("gamma detector"):
            sec = d["gamma_detectors"][self._index(section[len("gamma "):], "gamma_detectors")]
        elif section.startswith("detector"):
            sec = d["detectors"][self._index(section)]
        elif section in ("title", "description"):
            d[section] = value
            return self._apply(d)
        else:
            sec = d.setdefault(section, {})
        if value is None:
            sec.pop(field, None)
        else:
            sec[field] = value
        if section == "backing" and not sec:
            d["target"].pop("backing")
        return self._apply(d)

    def _set_reaction_type(self, d: dict, value: str) -> bool:
        reaction = d.setdefault("reaction", {"type": "elastic"})
        if value == "coulex":
            reaction.update(self._stashed_excitation or {"excite": "target", "multipolarity": "E2"})
            reaction["type"] = "coulex"
        else:
            kept = {k: v for k, v in reaction.items() if k != "type"}
            if kept:
                self._stashed_excitation = kept
            d["reaction"] = {"type": value}
        return self._apply(d)

    def _index(self, which, kind: str = "detectors") -> int:
        n = len(self._draft.get(kind, []))
        if isinstance(which, str):
            if which.startswith("detector "):
                i = int(which.split()[1]) - 1
            else:
                names = [det.get("name") for det in self._draft.get(kind, [])]
                if which not in names:
                    raise KeyError(f"no detector named {which!r}")
                i = names.index(which)
        else:
            i = int(which)
        if not 0 <= i < n:
            raise IndexError(f"there are {n} {kind.replace('_', ' ')}")
        return i

    def add_gamma_detector(self, **fields) -> bool:
        """Append a γ-ray detector with the given fields (as in a setup file's ``[[gamma_detectors]]``)."""
        d = self.draft
        d.setdefault("gamma_detectors", []).append(dict(fields))
        return self._apply(d)

    def remove_gamma_detector(self, which) -> bool:
        """Remove a γ-ray detector, by index (from 0) or name."""
        d = self.draft
        d["gamma_detectors"].pop(self._index(which, "gamma_detectors"))
        if not d["gamma_detectors"]:
            d.pop("gamma_detectors")
        return self._apply(d)

    def duplicate_gamma_detector(self, which, **changes) -> bool:
        """Copy a γ-ray detector (by index or name), apply ``changes``, and append it under a new name."""
        return self._duplicate("gamma_detectors", which, changes)

    def add_detector(self, **fields) -> bool:
        """Append a detector with the given fields (as in a setup file)."""
        d = self.draft
        d["detectors"].append(dict(fields))
        return self._apply(d)

    def remove_detector(self, which) -> bool:
        """Remove a detector, by index (from 0) or name."""
        d = self.draft
        d["detectors"].pop(self._index(which))
        return self._apply(d)

    def duplicate_detector(self, which, **changes) -> bool:
        """Copy a detector (by index or name), apply ``changes``, and append it; a copy gets a new name."""
        return self._duplicate("detectors", which, changes)

    def _duplicate(self, kind: str, which, changes: dict) -> bool:
        d = self.draft
        det = copy.deepcopy(d[kind][self._index(which, kind)])
        names = {x.get("name") for x in d[kind]}
        base = det.get("name") or "D"
        k = 2
        while f"{base}-{k}" in names:
            k += 1
        det["name"] = f"{base}-{k}"
        det.update(changes)
        d[kind].append(det)
        return self._apply(d)

    # -- shared results -----------------------------------------------------------------------------------------

    def _rates(self) -> Rates:
        if "rates" not in self._cache:
            self._cache["rates"] = Rates(self.experiment, previous=self._earlier_rates)
            self._earlier_rates = None
        return self._cache["rates"]

    # -- the scene ----------------------------------------------------------------------------------------------

    def live(self, key: str, fields: Optional[dict] = None) -> dict:
        """The quantities that follow a detector while it is dragged, computed quickly: its angles, solid angle,
        how much of it other detectors hide, and (for a particle detector) a coarse estimate of its rate.
        ``fields`` are the setup fields of the trial position (see :func:`physim.nuclear.scene.move_fields`); the
        setup itself is not changed. "problem" says why the position is not allowed, if it is not."""
        try:
            exp = _scene.moved(self.experiment, key, fields) if fields else self.experiment
        except SetupError as e:
            return {"problem": "; ".join(e.problems)}
        found = _scene.problems(exp, only=key)
        i = int(key.split(":")[1])
        if key.startswith("gamma"):
            gd = exp.gamma_detectors[i]
            return {"name": gd.name or f"γ{i + 1}", "theta": _q(gd.theta).to("deg"),
                    "phi": _q(gd.phi).to("deg") if gd.phi is not None else 0.0,
                    "distance_mm": _q(gd.distance).to("mm"), "half_angle_deg": gd.half_angle_deg(),
                    "geometric_efficiency": gd.geometric_efficiency(),
                    "peak_efficiency": gd.peak_efficiency(self.gamma_energy_mev(), exp),
                    "problem": found[0].text if found else None}
        array = Array.from_experiment(exp)
        g = array.geometries[i]
        quick = Rates(exp, depth_points=2, order=3, only=g.name)
        hidden = sum(frac for (_, behind), frac in array.shadowing().items() if behind == g.name)
        x, y, z = g.centre
        return {"name": g.name, "theta": math.degrees(math.atan2(math.hypot(x, y), z)),
                "phi": math.degrees(math.atan2(y, x)), "distance_mm": float(np.linalg.norm(g.centre)),
                "theta_range": g.theta_range(), "solid_angle_msr": g.solid_angle(), "hidden": min(hidden, 1.0),
                "rate_per_s": quick.rate(g.name, counted=False), "problem": found[0].text if found else None}

    def gamma_energy_mev(self) -> float:
        """The γ-ray energy the efficiencies are quoted at: that of the excited state, or 1332 keV without one."""
        exc = self.experiment.excitation
        return exc.energy_mev if exc is not None else _response.REFERENCE_MEV

    def efficiency(self, key: str, points: int = 120) -> dict:
        """The efficiency curve of γ-ray detector ``key`` ("gamma:N", from 0): full-energy-peak and total
        efficiency (fractions of all γ rays emitted at the target) and the resolution against energy, with what
        the response is built from (see :mod:`physim.nuclear.response`)."""
        gd = self.experiment.gamma_detectors[int(key.split(":")[1])]
        r = _response.Response(self.experiment, gd)
        e = np.geomspace(0.03, 3.0, points)
        return {"name": gd.name or f"γ{int(key.split(':')[1]) + 1}", "energy_mev": e,
                "peak": r.peak_efficiency(e), "total": r.total_efficiency(e), "fwhm_kev": 1e3 * r.fwhm(e),
                "transmission": r.transmission(e), **r.describe()}

    def source_run(self, nuclide: str = "152Eu", activity: str = "37 kBq", time: str = "1 h", seed: int = 1):
        """A simulated run with a calibration source at the target position, without beam
        (:func:`physim.nuclear.response.source_run`). Kept until the setup changes."""
        key = ("source", nuclide, activity, time, seed)
        if key not in self._cache:
            self._cache[key] = _response.source_run(self.experiment, nuclide, activity, time, seed)
        return self._cache[key]

    def safety(self) -> dict:
        """For Coulomb excitation, which rings (or strips) of each particle detector see collisions closer than
        Cline's safe distance: {detector: set of first segment indices}. Empty for elastic scattering.

        A ring is marked if the scattered beam or the recoil can reach it from a collision at a CM angle above the
        largest safe one."""
        if "safety" in self._cache:
            return self._cache["safety"]
        from .rates import coulex_for

        r = self._rates()
        out: dict = {}
        ion = beam_ion(self.experiment)
        th_cm = np.linspace(0.5, 180.0, 720)
        for ch in r.channels:
            if ch.excitation is None:
                continue
            safe = coulex_for(self.experiment, ch, r.layers).max_safe_angle()
            if safe >= 180.0:
                continue
            tb = TwoBody(ion, ch.nuclide.name, self.experiment.beam.energy_mev)
            close = th_cm > safe
            lab = np.r_[tb.at_cm(th_cm).theta_lab[close], tb.recoil_for(th_cm).theta_lab[close]]
            for g in r.array:
                for k in sorted({seg[0] for seg in g.segments}):
                    ranges = [g.theta_range(seg) for seg in g.segments if seg[0] == k]
                    lo, hi = min(a for a, _ in ranges), max(b for _, b in ranges)
                    if np.any((lab >= lo) & (lab <= hi)):
                        out.setdefault(g.name, set()).add(k)
        self._cache["safety"] = out
        return out

    def selection(self, key: Optional[str] = None, element: Optional[int] = None) -> dict:
        """What the scene's side panel shows for a selected object: a detector ("detector:N"), one of its rings or
        strips (``element``), a γ-ray detector ("gamma:N"), or with no key (or "target") the whole experiment."""
        r = self._rates()
        what = r.measured
        if key is None or key == "target":
            names = [g.name for g in r.array]
            return {"kind": "experiment", "title": self.experiment.title, "detectors": len(names),
                    "gamma_detectors": len(self.experiment.gamma_detectors),
                    "solid_angle_msr": sum(r.array[n].solid_angle() for n in names),
                    "rate_per_s": sum(r.rate(n) for n in names), "measured": what,
                    "measured_rate_per_s": sum(r.rate(n, what=what) for n in names),
                    "gamma_efficiency": r.gamma_efficiency()[0] if self.experiment.gamma_detectors else None,
                    "target": {"material": self.experiment.target.material,
                               "thickness_um": 1e3 * _scene.target_thickness_mm(self.experiment),
                               "tilt_deg": tilt_deg(self.experiment)},
                    "advice": _scene.advice(self.experiment)}
        i = int(key.split(":")[1])
        if key.startswith("gamma"):
            gd = self.experiment.gamma_detectors[i]
            out = dict(self.live(key), kind="gamma", model=gd.model)
            out["crystals"] = [label or "crystal" for label, _, _ in gd.elements()]
            out["gamma_energy_kev"] = 1e3 * self.gamma_energy_mev()
            out["response"] = _response.Response(self.experiment, gd).describe()
            out["coincidence_rate_per_s"] = None if what == "all" else out["peak_efficiency"] * sum(
                r.rate(g.name, what="excitations") for g in r.array)
            if element is not None and len(out["crystals"]) > 1:
                _, c, rad = gd.elements()[element]
                d = math.hypot(*c)
                out["element"] = {"label": f"crystal {out['crystals'][element]}",
                                  "theta": math.degrees(math.acos(c[2] / d)),
                                  "half_angle_deg": math.degrees(math.atan2(rad, d))}
            return out
        g = r.array.geometries[i]
        hidden = sum(frac for (_, behind), frac in r.array.shadowing().items() if behind == g.name)
        x, y, z = g.centre
        out = {"kind": "detector", "name": g.name, "model": self.experiment.detectors[i].model, "shape": g.shape,
               "theta": math.degrees(math.atan2(math.hypot(x, y), z)), "phi": math.degrees(math.atan2(y, x)),
               "distance_mm": float(np.linalg.norm(g.centre)), "theta_range": g.theta_range(),
               "phi_range": g.phi_range(), "solid_angle_msr": g.solid_angle(), "hidden": min(hidden, 1.0),
               "rate_per_s": r.rate(g.name), "measured": what, "measured_rate_per_s": r.rate(g.name, what=what),
               "counts_in_run": r.counts_in_run(g.name, what=what), "segments": len(g.segments),
               "unsafe": sorted(self.safety().get(g.name, ())), "on_axis": _scene.on_axis(self.experiment, key)}
        if element is not None and len({seg[0] for seg in g.segments}) > 1:
            segs = [seg for seg in g.segments if seg[0] == element]
            ranges = [g.theta_range(seg) for seg in segs]
            omega = g.segment_solid_angles()
            out["element"] = {"label": f"{'strip' if g.shape == 'rectangle' else 'ring'} {element + 1}",
                              "theta_range": (min(a for a, _ in ranges), max(b for _, b in ranges)),
                              "solid_angle_msr": sum(omega[seg] for seg in segs),
                              "rate_per_s": sum(r.rate(g.name, seg) for seg in segs),
                              "measured_rate_per_s": sum(r.rate(g.name, seg, what=what) for seg in segs),
                              "safe": element not in out["unsafe"]}
        return out

    def warnings(self) -> list:
        """Setup problems, then validity and rate warnings: everything the banner should show."""
        out = [Warning("error", "setup", p) for p in self.problems]
        for text in self._rates().warnings():
            level = "warning" if any(k in text for k in ("pile-up", "not Rutherford", "Mott", "stops inside",
                                                         "records nothing", "beam path", "blocks", "hides")) else "note"
            source = "rates" if ("counts" in text or "collects" in text or "records" in text) else "physics"
            out.append(Warning(level, source, text))
        for p in _scene.problems(self.experiment):
            if p.kind == "overlap":  # a detector in the beam is already reported above
                out.append(Warning("warning", "setup", p.text.replace("would overlap", "overlaps")
                                   + " They cannot both stand there."))
        order = {"error": 0, "warning": 1, "note": 2}
        return sorted(out, key=lambda w: order[w.level])

    def explain(self, tab: str) -> dict:
        """The "Explain" panel of a tab, with a link to its physics-register page."""
        if tab not in EXPLAIN:
            raise KeyError(f"unknown tab {tab!r}; tabs are {TABS}")
        return dict(EXPLAIN[tab], register=REGISTER[tab])

    def tab(self, name: str, **kw) -> dict:
        """The content of a result tab by name (see :data:`TABS`)."""
        if name not in TABS:
            raise KeyError(f"unknown tab {name!r}; tabs are {TABS}")
        return getattr(self, name)(**kw)

    # -- tabs -------------------------------------------------------------------------------------------------------

    def geometry(self) -> dict:
        """Beam axis, target and detector outlines in 3D (mm), with each detector's solid angle and angles."""
        array = Array.from_experiment(self.experiment)
        extent = max(float(np.linalg.norm(g.centre)) for g in array) * 1.25 if array.geometries else 100.0
        dets = []
        for g in array:
            lo, hi = g.theta_range()
            dets.append({"name": g.name, "outline": g.outline(), "centre": g.centre, "normal": g.n,
                         "solid_angle_msr": g.solid_angle(), "theta_range": (lo, hi), "phi_range": g.phi_range(),
                         "segments": len(g.segments)})
        gammas = []
        gap = np.full((1, 3), np.nan)
        for i, gd in enumerate(self.experiment.gamma_detectors):
            u = np.array(gd.direction())
            dist = _q(gd.distance).to("mm")
            a, b = (np.array(x) for x in gd.face_axes())
            t = np.linspace(0.0, 2 * np.pi, 49)[:, None]
            # The front face of each crystal, and the crystal's length behind it where the setup gives one.
            length = _q(gd.crystal_length).to("mm") if gd.crystal_length is not None else 0.0
            parts = []
            for _, centre, rad in gd.elements():
                ring = np.array(centre) + rad * (np.cos(t) * a + np.sin(t) * b)
                parts += [ring, gap]
                if length:
                    parts += [ring + length * u, gap]
                    for k in (0, 12, 24, 36):
                        parts += [np.array([ring[k], ring[k] + length * u]), gap]
            housing = gd.housing()
            if housing is not None:
                c, side = np.array(housing[0]), housing[1] / 2
                corners = [c + sa * side * a + sb * side * b for sa, sb in ((1, 1), (-1, 1), (-1, -1), (1, -1), (1, 1))]
                parts += [np.array(corners), gap]
            gammas.append({"name": gd.name or f"γ{i + 1}", "outline": np.concatenate(parts[:-1]),
                           "centre": dist * u, "theta": _q(gd.theta).to("deg"),
                           "half_angle_deg": gd.half_angle_deg(), "distance_mm": dist, "model": gd.model,
                           "crystals": len(gd.elements())})
        if gammas:
            extent = max(extent, 1.25 * max(g["distance_mm"] for g in gammas))
        blocking = [{"name": b.name, "outline": b.outline()} for b in array.blockers if "(housing)" not in b.name]
        return {"beam": np.array([[0.0, 0.0, -extent], [0.0, 0.0, extent]]), "extent": extent,
                "target_size_mm": 0.04 * extent, "detectors": dets, "gamma_detectors": gammas,
                "blocking": blocking}

    def kinematics(self, points: int = 361) -> dict:
        """Lab energy against lab angle for the scattered beam and the recoil of every target nuclide, with the
        angular range each detector covers."""
        beam = self.experiment.beam
        ion = beam_ion(self.experiment)
        th_cm = np.linspace(0.0, 180.0, points)
        curves = []
        for lay in stack(self.experiment):
            for key in lay.nuclides:
                nuc = data.nuclide(key).name
                tb = TwoBody(ion, nuc, beam.energy_mev)
                ej = tb.at_cm(th_cm)
                rec = tb.recoil_for(th_cm)
                label = f"{nuc} ({lay.name})"
                curves.append({"label": f"{ion} on {label}", "particle": "ejectile", "theta": ej.theta_lab,
                               "energy": ej.energy, "max_angle": tb.max_angle()})
                curves.append({"label": f"{nuc} recoil ({lay.name})", "particle": "recoil", "theta": rec.theta_lab,
                               "energy": rec.energy, "max_angle": tb.max_angle("recoil")})
        for ch in self._rates().channels:
            if ch.excitation is None:
                continue
            try:
                tb = TwoBody(ion, ch.nuclide.name, beam.energy_mev, **ch.kinematics_args)
            except ValueError:  # below the excitation threshold
                continue
            ej, rec = tb.at_cm(th_cm), tb.recoil_for(th_cm)
            star = f"{ch.excitation.energy_mev * 1e3:g} keV"
            who = ch.nuclide.name if ch.excitation.excite == "target" else ion
            curves.append({"label": f"{ion} on {ch.nuclide.name}, {who}* {star}", "particle": "ejectile",
                           "theta": ej.theta_lab, "energy": ej.energy, "max_angle": tb.max_angle(),
                           "inelastic": True})
            curves.append({"label": f"{ch.nuclide.name} recoil, {who}* {star}", "particle": "recoil",
                           "theta": rec.theta_lab, "energy": rec.energy, "max_angle": tb.max_angle("recoil"),
                           "inelastic": True})
        cover = [{"name": g.name, "theta_range": g.theta_range()} for g in Array.from_experiment(self.experiment)]
        return {"beam_energy_mev": beam.energy_mev, "curves": curves, "detectors": cover}

    def rates(self) -> dict:
        """One row per detector: solid angle, angles, rate, counts in the run, beam time for the counts wanted,
        relative error; and per-strip rates for each detector."""
        r = self._rates()
        what = r.measured
        rows, strips = [], {}
        for g in r.array:
            rate = r.rate(g.name)
            counts = r.counts_in_run(g.name, what=what)
            row = {"detector": g.name, "solid_angle_msr": g.solid_angle(), "theta_range": g.theta_range(),
                   "rate_per_s": rate, "rate_all_per_s": r.rate(g.name, counted=False),
                   "counts_in_run": counts,
                   "relative_error": r.relative_error(g.name, what=what) if counts > 0 else math.inf,
                   "beam_time_s": (r.beam_time_for(g.name, what=what) if r.counts_wanted is not None else None)}
            if what != "all":
                row["excitation_per_s"] = r.rate(g.name, what="excitations")
                row["coincidence_per_s"] = r.rate(g.name, what="coincidences") if what == "coincidences" else None
                # The γ-ray efficiency in coincidence with this detector, with the angular correlation.
                row["gamma_efficiency"] = r.gamma_efficiency(g.name)[0] if what == "coincidences" else None
            rows.append(row)
            strips[g.name] = r.per_segment(g.name)
        eff, typical = r.gamma_efficiency() if what == "coincidences" else (None, False)
        return {"rows": rows, "strips": strips, "beam_time_s": r.beam_time_s, "counts_wanted": r.counts_wanted,
                "particles_per_second": r.particles_per_second, "measured": what, "gamma_efficiency": eff,
                "gamma_efficiency_typical": typical}

    def energy_loss(self) -> dict:
        """The beam through each layer (energy in and out, loss, straggling), the beam energy through the target,
        and each detector's dead-layer loss and punch-through energy for the scattered beam."""
        exp = self.experiment
        ion = beam_ion(exp)
        layers = stack(exp)
        rows = []
        e = exp.beam.energy_mev
        c = math.cos(math.radians(tilt_deg(exp)))
        profile_z, profile_e, offset = [], [], 0.0
        for i, lay in enumerate(layers):
            st = stopping(ion, lay.material)
            path = lay.thickness / c
            out = st.energy_after(e, path)
            rows.append({"layer": lay.name, "material": lay.material.name, "thickness_mg_cm2": lay.thickness,
                         "energy_in_mev": e, "energy_out_mev": out, "loss_mev": e - out,
                         "straggling_fwhm_mev": 2.3548 * st.straggling(e, path, steps=50) if out > 0 else 0.0})
            z = np.linspace(0.0, lay.thickness, 41)
            profile_z.append(offset + z)
            profile_e.append(beam_energy_at(exp, i, z, layers))
            offset += lay.thickness
            e = out
        main = data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        dets = []
        for det, g in zip(exp.detectors, Array.from_experiment(exp)):
            resp = Response(ion, det)
            k = float(TwoBody(ion, main, exp.beam.energy_mev).at_lab(g.mean_theta())[0].energy)
            dets.append({"detector": g.name, "ejectile_energy_mev": k,
                         "after_dead_layer_mev": float(resp.after_dead_layer(k)) if np.isfinite(k) else math.nan,
                         "punch_through_mev": resp.punch_through_energy()})
        return {"layers": rows, "depth_mg_cm2": np.concatenate(profile_z), "beam_energy_mev": np.concatenate(profile_e),
                "detectors": dets}

    def spectra(self, events: int = 200_000, seed: int = 1, bins: int = 200) -> dict:
        """Simulated measured-energy spectra per detector (counts in the planned beam time per bin)."""
        key = ("events", events, seed)
        if key not in self._cache:
            from .events import simulate

            self._cache[key] = simulate(self.experiment, events, seed)
        ev = self._cache[key]
        out = {}
        for name in ev.detectors:
            m = ev.select(name)
            if not m.any():
                out[name] = {"counts": np.zeros(bins), "edges": np.linspace(0, 1, bins + 1)}
                continue
            hi = float(ev["measured"][m].max()) * 1.03
            h, edges = ev.spectrum(name, bins=bins, range=(0.0, hi))
            out[name] = {"counts": h, "edges": edges}
        return {"events": ev, "spectra": out}

    def gamma_events(self, events: int = 400_000, seed: int = 1):
        """γ rays in coincidence with the detected particles (:func:`physim.nuclear.gamma_events.simulate_gammas`),
        built on the same particle events as :meth:`spectra`. Kept until the setup changes."""
        from .events import simulate
        from .gamma_events import simulate_gammas

        key = ("gammas", events, seed)
        if key not in self._cache:
            ev_key = ("events", events, seed)
            if ev_key not in self._cache:
                self._cache[ev_key] = simulate(self.experiment, events, seed)
            self._cache[key] = simulate_gammas(self.experiment, events, seed, particle_events=self._cache[ev_key])
        return self._cache[key]

    def gamma_spectra(self, events: int = 400_000, seed: int = 1, bins: int = 300) -> dict:
        """The γ-ray side of the Monte Carlo: for each γ-ray detector the raw and Doppler-corrected spectra in
        coincidence with all particle detectors (counts in the run per bin, random coincidences included), the
        particle × γ matrix of true and random coincidences, and the live fraction.
        ``{"available": False}`` without an excited state or γ-ray detectors."""
        exp = self.experiment
        if exp.excitation is None or not exp.gamma_detectors:
            return {"available": False, "reason": "Coulomb excitation with γ-ray detectors is needed."}
        g = self.gamma_events(events, seed)
        emitter = "recoil" if exp.excitation.excite == "target" else "projectile"
        spectra = {}
        for name in g.detector_names():
            total = None
            out = {}
            for key in (None, "projectile", "recoil"):
                parts = [g.spectrum(name, d, corrected=key, bins=bins, range=(0.0, 1.3 * g.energy_mev))
                         for d in g.events.detectors]
                counts = np.sum([h for h, _ in parts], axis=0)
                edges = parts[0][1]
                out[key or "measured"] = counts
                total = edges
            spectra[name] = {"edges": total, **out, "counts": g.counts(name)}
        return {"available": True, "spectra": spectra, "emitter": emitter, "energy_kev": 1e3 * g.energy_mev,
                "coincidences": g.coincidences(), "live_fraction": g.live_fraction, "window_s": g.window_s,
                "gamma_detectors": g.detector_names(), "particle_detectors": list(g.events.detectors),
                "notes": g.notes, "n_gammas": len(g)}

    def analysis(self, events: int = 400_000, seed: int = 1):
        """The :class:`~physim.nuclear.analysis.Analysis` of the simulated γ rays (:meth:`gamma_events`), kept
        with its history until the setup changes."""
        from .analysis import Analysis

        key = ("analysis", events, seed)
        if key not in self._cache:
            self._cache[key] = Analysis(self.experiment, self.gamma_events(events, seed))
        return self._cache[key]

    def analyse(self, settings=None, events: int = 400_000, seed: int = 1) -> dict:
        """Run the automatic analysis with ``settings`` (a :class:`~physim.nuclear.analysis.Settings`, or a dict
        of its fields) and return the result as a dict for display. ``{"available": False}`` without Coulomb
        excitation and γ-ray detectors."""
        from dataclasses import asdict

        from .analysis import Settings

        exp = self.experiment
        if exp.excitation is None or not exp.gamma_detectors:
            return {"available": False, "reason": "Coulomb excitation with γ-ray detectors is needed."}
        if isinstance(settings, dict):
            settings = Settings(**settings)
        try:
            r = self.analysis(events, seed).run(settings)
        except ValueError as e:
            return {"available": False, "reason": str(e)}
        out = asdict(r)
        out.update(available=True, total_unc=r.total_unc, runs=len(self.analysis(events, seed).history))
        return out

    def explanations(self, detector: Optional[str] = None, gamma_detector: Optional[str] = None,
                     result=None) -> list:
        """Every number of the plan (and of the analysis ``result``, if given) with its explanation, as dicts
        (:mod:`physim.nuclear.record`)."""
        from .record import analysis_explanations, explanations

        out = [x.as_dict() for x in explanations(self, detector, gamma_detector)]
        if result is not None:
            out += [x.as_dict() for x in analysis_explanations(self, result)]
        return out

    def record_html(self, result=None, scene_png: Optional[str] = None) -> str:
        """The run record as one HTML page (:func:`physim.nuclear.record.record_html`); ``result`` is the last
        analysis result if None and one has been run."""
        from .record import record_html

        if result is None:
            for key, value in self._cache.items():
                if key[0] == "analysis" and value.history:
                    result = value.history[-1]
        return record_html(self, result, scene_png=scene_png)

    def multistep(self, role: Optional[str] = None, angle_step: float = 4.0, energies: int = 2,
                  shapes: bool = True) -> dict:
        """Multi-step Coulomb excitation of the setup's level scheme (:mod:`physim.nuclear.multistep`): the γ
        yields of every transition in every particle detector with all orders and with first order, the
        populated levels, and (``shapes``) the prolate–zero–oblate comparison of the first 2⁺ state. Uses the
        scheme of ``role`` ("target" or "beam"; the excited nucleus's by default). Kept until the setup changes.
        ``{"available": False}`` without a level scheme."""
        from .multistep import Multistep

        exp = self.experiment
        if role is None:
            wanted = None if exp.excitation is None else ("target" if exp.excitation.excite == "target" else "beam")
            role = wanted if wanted in exp.levels else next(iter(exp.levels), None)
        if role not in exp.levels:
            return {"available": False, "reason": "Look up a level scheme first (Excitation and γ rays)."}
        key = ("multistep", role, angle_step, energies, shapes)
        if key in self._cache:
            return self._cache[key]
        scheme = exp.levels[role]
        ms = Multistep(exp, scheme, role, angle_step=angle_step, energies=energies)
        try:
            full = ms.yields()
            first = Multistep(exp, scheme, role, angle_step=angle_step, energies=energies, first_order=True).yields()
        except ValueError as e:
            return {"available": False, "reason": str(e)}
        labels = {t: f"{scheme.levels[t[0]].label} ({scheme.levels[t[0]].energy.value:g} keV) → "
                     f"{scheme.levels[t[1]].label}" for t in full.transitions}
        out = {"available": True, "role": role, "nuclide": scheme.nuclide, "detectors": list(full.detectors),
               "transitions": [{"transition": t, "label": labels[t], "energy_kev": scheme.levels[t[0]].energy.value
                                - scheme.levels[t[1]].energy.value,
                                "all_orders": {d: full.detectors[d].get(t, 0.0) for d in full.detectors},
                                "first_order": {d: first.detectors[d].get(t, 0.0) for d in first.detectors}}
                               for t in full.transitions],
               "levels": [{"index": n, "label": scheme.levels[n].label, "energy_kev": scheme.levels[n].energy.value,
                           "all_orders": {d: full.levels[d].get(n, 0.0) for d in full.levels},
                           "first_order": {d: first.levels[d].get(n, 0.0) for d in first.levels}}
                          for n in range(1, len(scheme.levels)) if any(full.levels[d].get(n, 0.0) > 0
                                                                       for d in full.levels)],
               "notes": full.notes, "gosia_input": ms.gosia_input(full), "beam_time_s": _q(exp.run.beam_time).to("s")}
        if shapes:
            twos = [n for n, lev in enumerate(scheme.levels) if n and lev.spin == 2 and scheme.b(0, n, "E2")]
            if twos:
                sh = ms.shapes(twos[0])
                out["shapes"] = {"level": twos[0], "q_efm2": sh["q_efm2"], "rings": sh["rings"],
                                 "totals": sh["totals"], "separable": sh["separable"],
                                 "total_difference_sigma": sh["total_difference_sigma"],
                                 "transition": sh["transition"]}
        self._cache[key] = out
        return out

    def fit_matrix_elements(self, measured: dict, free: list, role: Optional[str] = None) -> dict:
        """Fit up to three matrix elements of the level scheme to measured counts (see
        :meth:`physim.nuclear.multistep.Multistep.fit`)."""
        from .multistep import Multistep

        exp = self.experiment
        if role is None:
            role = "target" if exp.excitation is None or exp.excitation.excite == "target" else "beam"
        return Multistep(exp, exp.levels[role], role).fit(measured, free)

    def trajectories(self, impact_parameters: Optional[list] = None, nuclide: Optional[str] = None) -> dict:
        """Coulomb orbits (CM frame, fm) for a range of impact parameters, from physim's engine."""
        exp = self.experiment
        ion = beam_ion(exp)
        target = nuclide or data.nuclide(max(stack(exp)[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        r = Rutherford(ion, target, exp.beam.energy_mev)
        bs = np.asarray(impact_parameters if impact_parameters is not None
                        else r.d0 * np.array([0.1, 0.25, 0.5, 1.0, 2.0, 4.0]), dtype=float)
        orbits = r.trajectories(bs, distance=40 * r.d0)
        return {"target": target, "d0_fm": r.d0, "interaction_radius_fm": r.interaction_radius,
                "orbits": [{"b_fm": o.b, "xy": o.pos, "deflection_deg": o.deflection(),
                            "closest_fm": o.closest_approach()} for o in orbits],
                "grazing_angle_deg": r.grazing_angle()}

    def gamma(self) -> dict:
        """Coulomb excitation: the excited state, its cross section, the excitation probability against angle, and
        the Doppler-shifted γ-ray energy for every particle × γ detector pair. ``{"available": False}`` without an
        excited state."""
        exc = self.experiment.excitation
        if exc is None:
            return {"available": False, "reason": "The setup has no excited state ([reaction] type = \"coulex\")."}
        from .gamma import doppler_table
        from .rates import coulex_for

        r = self._rates()
        ch = next(c for c in r.channels if c.excitation is not None)
        cx = coulex_for(self.experiment, ch, r.layers)
        th = np.linspace(1.0, 180.0, 180)
        rates = {g.name: sum(v for (label, _), v in r.by_channel(g.name).items() if label == ch.label)
                 for g in r.array}
        if "gamma" not in self._cache:
            self._cache["gamma"] = doppler_table(self.experiment) if self.experiment.gamma_detectors else []
        return {"available": True, "state": {"excite": exc.excite, "energy_kev": exc.energy_mev * 1e3,
                                             "multipolarity": exc.multipolarity, "b_up_e2fm": exc.b_up_e2fm},
                "xi": cx.xi, "eta": cx.eta, "safe_distance_fm": cx.safe_distance, "max_safe_angle": cx.max_safe_angle(),
                "total_mb": cx.total(), "theta_cm": th, "probability": cx.probability(th),
                "rates": rates, "particles": self._particle_energies(ch), "doppler": self._cache["gamma"],
                "emission": exc.emission or "correlated", "correlation": self.correlation()}

    def correlation(self) -> list:
        """For every particle detector × γ-ray crystal, the γ rays seen in coincidence relative to isotropic
        emission (:func:`physim.nuclear.gamma.correlation_table`). Empty without γ-ray detectors."""
        if "correlation" not in self._cache:
            from .gamma import correlation_table

            exp = self.experiment
            self._cache["correlation"] = (correlation_table(exp) if exp.excitation is not None
                                          and exp.gamma_detectors else [])
        return self._cache["correlation"]

    def populations(self, role: Optional[str] = None) -> dict:
        """First-order Coulomb excitation of a whole level scheme of the setup (``role`` "target" or "beam"; by
        default the nucleus the reaction excites, if it has a scheme, else whichever has one): for each level the
        cross section for exciting it directly and with feeding from above, and for each transition the cross
        section of its γ ray, integrated over all scattering angles at the mid-target energy.
        ``{"available": False}`` without a level scheme."""
        from .orientation import Excitation

        exp = self.experiment
        if role is None:
            wanted = None if exp.excitation is None else ("target" if exp.excitation.excite == "target" else "beam")
            role = wanted if wanted in exp.levels else next(iter(exp.levels), None)
        if role not in exp.levels:
            return {"available": False, "reason": "Look up a level scheme first."}
        key = ("populations", role)
        if key not in self._cache:
            scheme = exp.levels[role]
            layers = stack(exp)
            e_mid = round(float(beam_energy_at(exp, 0, [layers[0].thickness / 2], layers)[0]), 9)
            target = data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
            if role == "target":
                target = scheme.nuclide
            try:
                ex = Excitation(beam_ion(exp), target, e_mid, scheme,
                                excite="target" if role == "target" else "projectile")
            except ValueError as e:
                return {"available": False, "reason": str(e)}
            cs = ex.cross_sections()
            levels = [{"index": n, "label": lev.label, "energy_kev": lev.energy.value,
                       "direct_mb": cs["direct"].get(n, 0.0), "populated_mb": cs["populated"].get(n, 0.0)}
                      for n, lev in enumerate(scheme.levels) if n and cs["populated"].get(n, 0.0) > 0]
            gammas = [{"initial": i, "final": f, "energy_kev": scheme.levels[i].energy.value
                       - scheme.levels[f].energy.value, "label": f"{scheme.levels[i].label} → "
                       f"{scheme.levels[f].label}", "gamma_mb": v}
                      for (i, f), v in cs["gamma"].items() if v > 0]
            self._cache[key] = {"available": True, "role": role, "nuclide": scheme.nuclide, "levels": levels,
                                "gammas": sorted(gammas, key=lambda g: -g["gamma_mb"]), "notes": list(ex.notes),
                                "beam_energy_mev": e_mid}
        return self._cache[key]

    # -- level schemes ----------------------------------------------------------------------------------------------

    def level_nuclides(self) -> dict:
        """The nuclide whose level scheme each role takes: the beam, and the target's most abundant nuclide (all of
        them are listed under "target_choices")."""
        d = self._draft
        out = {"beam": None, "target": None, "target_choices": []}
        try:
            z, a = parse_nuclide(d["beam"]["nuclide"])
            out["beam"] = _data.nuclide((z, a)).name
        except (KeyError, ValueError):
            pass
        try:
            atoms = sorted(_data.material(d["target"]["material"]).atoms, key=lambda x: -x[2])
            out["target_choices"] = [_data.nuclide((z, a)).name for z, a, _ in atoms]
            out["target"] = out["target_choices"][0]
        except (KeyError, ValueError, IndexError):
            pass
        return out

    def levels(self) -> dict:
        """The level schemes in the setup, as tables for display (:meth:`LevelScheme.table`), with whether a local
        ENSDF copy is present and which nuclide each role would look up."""
        schemes = self.experiment.levels if not self.problems else {}
        return {"ensdf": _ensdf.available(), "ensdf_folder": str(_ensdf.folder()), "nuclides": self.level_nuclides(),
                "beam": schemes["beam"].table() if "beam" in schemes else None,
                "target": schemes["target"].table() if "target" in schemes else None}

    def lookup_levels(self, role: str, max_energy: str = "3 MeV", nuclide: Optional[str] = None) -> bool:
        """Read the level scheme of the beam or the target nuclide from the local ENSDF copy into the setup, up to
        ``max_energy``. Matrix elements the user had set for the same nuclide are kept. Raises
        :class:`~physim.nuclear.ensdf.EnsdfMissing` without a copy or without data for the nuclide."""
        if role not in ("beam", "target"):
            raise ValueError("role must be 'beam' or 'target'")
        nuclide = nuclide or self.level_nuclides()[role]
        if nuclide is None:
            raise ValueError(f"the setup has no {role} nuclide yet")
        scheme = LevelScheme.from_ensdf(nuclide, max_energy_kev=_q(max_energy).to("keV"))
        d = self.draft
        old = d.get("levels", {}).get(role)
        if old and old.get("nuclide") == scheme.nuclide:
            energies = [lv.energy.value for lv in scheme.levels]
            old_energy = [_q(lv["energy"]).to("keV") for lv in old.get("level", [])]
            for m in old.get("matrix_element", []):
                if m.get("source") != "user":
                    continue
                try:  # the same pair of levels, found by energy
                    a, b = (energies.index(old_energy[m[k]]) for k in ("from", "to"))
                except (ValueError, IndexError, KeyError):
                    continue
                scheme.set_matrix_element(a, b, m.get("multipolarity", "E2"), m["value"], m.get("value_unc"),
                                          note=m.get("note", ""))
        d.setdefault("levels", {})[role] = scheme.to_dict()
        return self._apply(d)

    def remove_levels(self, role: str) -> bool:
        """Take a level scheme out of the setup."""
        d = self.draft
        d.get("levels", {}).pop(role, None)
        if not d.get("levels"):
            d.pop("levels", None)
        return self._apply(d)

    def set_matrix_element(self, role: str, a: int, b: int, multipolarity: str, value: str) -> bool:
        """Set a matrix element of a level scheme (text with a unit, e.g. ``"1.28 eb"``); it is marked as the
        user's."""
        scheme = self.experiment.levels[role]
        scheme.set_matrix_element(a, b, multipolarity, value)
        d = self.draft
        d["levels"][role] = scheme.to_dict()
        return self._apply(d)

    def use_state(self, role: str, level: int, multipolarity: str = "E2") -> bool:
        """Plan Coulomb excitation of one state of a level scheme: sets the reaction's excited nucleus, energy,
        multipolarity and B(Eλ↑) from the scheme."""
        scheme = self.experiment.levels[role]
        b_up = scheme.b(0, level, multipolarity)
        if not b_up:
            raise ValueError(f"the scheme has no {multipolarity} matrix element between the ground state and "
                             f"level {level}")
        lam = int(multipolarity[1])
        d = self.draft
        d["reaction"] = {"type": "coulex", "excite": "target" if role == "target" else "projectile",
                         "energy": f"{scheme.levels[level].energy.value:.10g} keV", "multipolarity": multipolarity,
                         "b_up": f"{b_up:.6g} e2fm{2 * lam}"}
        return self._apply(d)

    def _particle_energies(self, channel) -> list:
        """What each particle detector measures in a Coulomb-excitation run: the energy of the beam particle or
        recoil at the detector's smallest, central and largest lab angle, after elastic scattering and after exciting
        the state, and the speed β of the excited nucleus in that event (for the Doppler correction). Energies are at
        the reaction point, before losses in the target."""
        beam = self.experiment.beam
        ion = beam_ion(self.experiment)
        elastic = TwoBody(ion, channel.nuclide.name, beam.energy_mev)
        try:
            inelastic = TwoBody(ion, channel.nuclide.name, beam.energy_mev, **channel.kinematics_args)
        except ValueError:  # below the excitation threshold
            return []
        excited = channel.kinematics_args["excite"]  # "recoil" or "ejectile"
        m_exc = inelastic.m4 if excited == "recoil" else inelastic.m3

        def beta(t, m):
            return float(np.sqrt(t * (t + 2 * m)) / (t + m))

        rows = []
        for g in Array.from_experiment(self.experiment):
            lo, hi = g.theta_range()
            for particle, name in (("ejectile", ion), ("recoil", channel.nuclide.name)):
                reach = min(elastic.max_angle(particle), inelastic.max_angle(particle))
                for where, th in (("min", lo), ("centre", 0.5 * (lo + hi)), ("max", hi)):
                    if not 0 < th <= reach:
                        continue
                    el = elastic.at_lab(th, particle)[0]
                    ex = inelastic.at_lab(th, particle)[0]
                    if not (np.isfinite(el.energy) and np.isfinite(ex.energy)):
                        continue
                    # The excited nucleus is the detected particle itself or its partner in the same event.
                    if particle == excited:
                        t_exc = float(ex.energy)
                    else:
                        partner = "recoil" if particle == "ejectile" else "ejectile"
                        t_exc = float(inelastic.at_cm(180.0 - float(ex.theta_cm), partner).energy)
                    rows.append({"detector": g.name, "particle": particle, "nuclide": name, "where": where,
                                 "theta_lab": float(th), "elastic_mev": float(el.energy),
                                 "excited_mev": float(ex.energy),
                                 "difference_mev": float(el.energy) - float(ex.energy),
                                 "beta_excited": beta(t_exc, m_exc)})
        return rows

    def report(self) -> dict:
        """What the beam-time report will contain (the export itself is item 42)."""
        return {"title": self.experiment.title, "setup": self.experiment.to_dict(), "rates": self.rates()["rows"],
                "peaks": {g.name: self._rates().peaks(g.name) for g in self._rates().array},
                "warnings": self.warnings()}

    # -- sweeps -----------------------------------------------------------------------------------------------------

    def sweep(self, parameter: str, values: list, quantity: str = "rate", detector: Optional[str] = None) -> dict:
        """Vary one parameter and recompute one quantity for a detector.

        ``parameter``: "beam energy", "target thickness", or "detector angle" (of ``detector``). ``quantity``:
        "rate" (counts/s), "beam time" (s for the counts wanted), "peak energy" (MeV, the strongest peak) or
        "peak width" (FWHM, MeV). Values are given with units, as in a setup file.
        """
        detector = detector or self.experiment.detectors[0].name
        out = []
        for v in values:
            d = self.draft
            if parameter == "beam energy":
                d["beam"]["energy"] = v
            elif parameter == "target thickness":
                d["target"]["thickness"] = v
            elif parameter == "detector angle":
                det = d["detectors"][self._index(detector)]
                det["theta"] = v
                det.pop("position", None)
            else:
                raise ValueError("parameter must be 'beam energy', 'target thickness' or 'detector angle'")
            r = Rates(Experiment.from_dict(d))
            if quantity == "rate":
                y = r.rate(detector)
            elif quantity == "beam time":
                y = r.beam_time_for(detector, what="measured")
            elif quantity in ("peak energy", "peak width"):
                pk = r.peaks(detector)
                y = (pk[0].mean if quantity == "peak energy" else pk[0].fwhm) if pk else math.nan
            else:
                raise ValueError("quantity must be 'rate', 'beam time', 'peak energy' or 'peak width'")
            out.append(y)
        return {"parameter": parameter, "values": list(values), "x": [_q(v).value for v in values],
                "quantity": quantity, "detector": detector, "y": np.array(out, dtype=float)}


__all__ = ["EXPLAIN", "REGISTER", "TABS", "Planner", "Warning"]
