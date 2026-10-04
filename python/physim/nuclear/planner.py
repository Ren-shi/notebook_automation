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
from .detectors import Array, Response
from .experiment import Experiment, SetupError, example_names
from .kinematics import TwoBody
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
                       "leaving the target, isotropic.",
        "limits": "Closer than Cline's safe distance nuclear forces interfere; P ≳ 0.1 needs multi-step excitation "
                  "(GOSIA); lifetimes, angular distributions and feeding are not modelled.",
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
        self._cache.clear()
        return True

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
            self._cache["rates"] = Rates(self.experiment)
        return self._cache["rates"]

    def warnings(self) -> list:
        """Setup problems, then validity and rate warnings: everything the banner should show."""
        out = [Warning("error", "setup", p) for p in self.problems]
        for text in self._rates().warnings():
            level = "warning" if any(k in text for k in ("pile-up", "not Rutherford", "Mott", "stops inside",
                                                         "records nothing", "beam path", "blocks")) else "note"
            source = "rates" if ("counts" in text or "collects" in text or "records" in text) else "physics"
            out.append(Warning(level, source, text))
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
        for i, gd in enumerate(self.experiment.gamma_detectors):
            u = np.array(gd.direction())
            dist, rad = _q(gd.distance).to("mm"), _q(gd.radius).to("mm")
            # Two unit vectors across the detector face, for its outline.
            a = np.cross(u, [0.0, 0.0, 1.0] if abs(u[2]) < 0.9 else [1.0, 0.0, 0.0])
            a /= np.linalg.norm(a)
            b = np.cross(u, a)
            t = np.linspace(0.0, 2 * np.pi, 49)[:, None]
            outline = dist * u + rad * (np.cos(t) * a + np.sin(t) * b)
            gammas.append({"name": gd.name or f"γ{i + 1}", "outline": outline, "centre": dist * u,
                           "theta": _q(gd.theta).to("deg"), "half_angle_deg": gd.half_angle_deg(),
                           "distance_mm": dist})
        if gammas:
            extent = max(extent, 1.25 * max(g["distance_mm"] for g in gammas))
        return {"beam": np.array([[0.0, 0.0, -extent], [0.0, 0.0, extent]]), "extent": extent,
                "target_size_mm": 0.04 * extent, "detectors": dets, "gamma_detectors": gammas}

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
            rows.append(row)
            strips[g.name] = r.per_segment(g.name)
        eff, geometric = r.gamma_efficiency() if what == "coincidences" else (None, False)
        return {"rows": rows, "strips": strips, "beam_time_s": r.beam_time_s, "counts_wanted": r.counts_wanted,
                "particles_per_second": r.particles_per_second, "measured": what, "gamma_efficiency": eff,
                "gamma_efficiency_geometric": geometric}

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
                "rates": rates, "particles": self._particle_energies(ch), "doppler": self._cache["gamma"]}

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
