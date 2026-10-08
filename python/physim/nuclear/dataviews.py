"""The Data tab's model (backlog item 67): every spectrum of a run, the operations applied to them, gates and
conditions as named objects, the 2D views, run overlays and exports.

A run's data are views, not pre-baked results: the user picks the Doppler correction, a particle gate, singles or
coincidences with the randoms shown or subtracted, add-back on or off, and the binning::

    from physim.nuclear import dataviews as dv

    run = planner.run                                     # the current run (physim.nuclear.runs.RunData)
    grid = dv.grid(run)                                   # one spectrum per particle and per γ-ray detector
    g = dv.Gate("excited CD", detectors=["CD"], group="excited")
    h = dv.gamma_spectrum(run, "Ge90", correction="recoil", gate=g)
    m = dv.energy_vs_ring(run, "CD")                      # the 2D view, with the kinematic line of each group
    dv.Gate.from_region("my gate", "CD", rings=(3, 12), energy=(25.0, 27.5))

**Gates** are kept with the experiment (``gates.json`` in its folder) and the analysis uses them by name
(:meth:`Gate.settings`). A gate's ``group`` picks the particles by their energy against the kinematic lines of the
elastic and excited groups (what an experimentalist can do), not by the simulation's knowledge of the channel; a
gate drawn on the energy-against-ring view is a ring range and an energy window.
"""

from __future__ import annotations

import csv
import io
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional, Union

import numpy as np

#: The Doppler corrections, with what each assumes.
CORRECTIONS = {"off": "No correction: the energy as measured, shifted and broadened by the motion of the emitter.",
               "projectile": "As if the scattered beam emitted the γ ray: its velocity from the particle's angle "
                             "and two-body kinematics.",
               "recoil": "As if the target recoil emitted the γ ray: the partner of the detected particle, its "
                         "velocity from two-body kinematics."}
#: How random coincidences are treated.
RANDOMS = {"shown": "Random coincidences included, as measured.",
           "subtracted": "Randoms subtracted: the true coincidences plus the randoms' Poisson noise.",
           "none": "True coincidences only (what a simulation can show and a measurement cannot)."}


def correction_reason(experiment, choice: str) -> str:
    """Why a Doppler correction is right or wrong for this setup, in one sentence."""
    exc = experiment.excitation
    if exc is None:
        return "No excited state: there is no γ ray to correct."
    emitter = "recoil" if exc.excite == "target" else "projectile"
    who = "the target nucleus" if emitter == "recoil" else "the beam nucleus"
    if choice == "off":
        return f"The reaction excites {who}, which emits the γ ray in flight: uncorrected, the line is shifted " \
               "by its velocity and spread over the crystals' angles."
    if choice == emitter:
        return f"Right for this setup: the reaction excites {who}, and this corrects for its velocity."
    return f"The reaction excites {who}, not the {choice}: this correction uses the wrong velocity and leaves " \
           "the line shifted and broad."


@dataclass
class Gate:
    """A named condition on the particles: which particle detectors, which rings (or strips, along u) counting from
    0, which group ("any", "excited" or "elastic", by energy against the kinematic lines), and an energy window
    (MeV) drawn on the energy-against-ring view."""

    name: str
    detectors: Optional[list] = None
    rings: Optional[list] = None
    group: str = "any"
    energy: Optional[list] = None
    note: str = ""

    def __post_init__(self):
        if self.group not in ("any", "excited", "elastic"):
            raise ValueError("a gate's group is 'any', 'excited' or 'elastic'")
        if self.rings is not None:
            self.rings = sorted(int(r) for r in self.rings)
        if self.energy is not None:
            lo, hi = (float(x) for x in self.energy)
            self.energy = [min(lo, hi), max(lo, hi)]

    @classmethod
    def from_region(cls, name: str, detector: str, rings: tuple, energy: tuple, note: str = "") -> Gate:
        """A gate drawn on the energy-against-ring view: rings ``rings[0]`` to ``rings[1]`` (from 0, inclusive) of
        ``detector`` and measured energies in ``energy`` (MeV)."""
        lo, hi = sorted(int(r) for r in rings)
        return cls(name, [detector], list(range(lo, hi + 1)), "any", list(energy), note)

    def as_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> Gate:
        return cls(**{k: d.get(k) for k in ("name", "detectors", "rings", "group", "energy", "note") if k in d})

    def describe(self) -> str:
        parts = [", ".join(self.detectors) if self.detectors else "all particle detectors"]
        if self.rings is not None:
            parts.append(f"rings {self.rings[0] + 1}–{self.rings[-1] + 1}" if self.rings == list(
                range(self.rings[0], self.rings[-1] + 1)) else "rings " + ", ".join(str(r + 1) for r in self.rings))
        if self.group != "any":
            parts.append(f"the {self.group} group")
        if self.energy is not None:
            parts.append(f"{self.energy[0]:.4g}–{self.energy[1]:.4g} MeV")
        return "; ".join(parts)

    def mask(self, events, experiment) -> np.ndarray:
        """Which particle rows of ``events`` pass the gate."""
        c = events.columns
        m = np.ones(len(events), dtype=bool)
        if self.detectors is not None:
            m &= np.isin(c["detector"], [events.detectors.index(d) for d in self.detectors if d in events.detectors])
        if self.rings is not None:
            m &= np.isin(c["segment_i"], self.rings)
        if self.energy is not None:
            m &= (c["measured"] >= self.energy[0]) & (c["measured"] <= self.energy[1])
        if self.group != "any":
            m &= group_of(events, experiment) == self.group
        return m

    def settings(self, **overrides):
        """The analysis settings this gate stands for (:class:`physim.nuclear.analysis.Settings`)."""
        from .analysis import Settings

        s = Settings(detectors=self.detectors, rings=self.rings,
                     particle_gate="inelastic" if self.group == "excited" else "all",
                     particle_energy=tuple(self.energy) if self.energy is not None else None)
        for k, v in overrides.items():
            setattr(s, k, v)
        return s


# -- kinematic lines and groups ----------------------------------------------------------------------------------

def kinematic_lines(experiment, detector: str) -> list:
    """The measured energy each group would have in each ring (or strip) of a particle detector: after the way
    out of the middle of the target and the detector's dead layer, at the ring's centre. A list of {"label",
    "group" ("elastic" or "excited"), "particle", "rings": [index], "energy": [MeV]}."""
    from . import data
    from .detectors import Array
    from .kinematics import TwoBody
    from .rates import beam_energy_at, beam_ion, exit_energy, stack

    array = Array.from_experiment(experiment)
    names = [g.name for g in array.geometries]
    g = array.geometries[names.index(detector)]
    layers = stack(experiment)
    e_mid = float(beam_energy_at(experiment, 0, [layers[0].thickness / 2], layers)[0])
    ion = beam_ion(experiment)
    target = data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
    rings = sorted({seg[0] for seg in g.segments})
    # Each ring's direction: the segment of the ring whose angle to the beam is nearest the ring's mean (the
    # sectors of an annular ring all have the same angle; the strips along v of a rectangle do not).
    centres = []
    for k in rings:
        segs = [seg for seg in g.segments if seg[0] == k]
        cs = np.array([g.segment_centre(seg, weighted=True) for seg in segs])
        cs = cs / np.linalg.norm(cs, axis=1, keepdims=True)
        th = np.degrees(np.arccos(np.clip(cs[:, 2], -1, 1)))
        centres.append(cs[int(np.argmin(np.abs(th - th.mean())))])
    dirs = np.array(centres)
    theta = np.degrees(np.arccos(np.clip(dirs[:, 2], -1, 1)))
    kinds = [("elastic", TwoBody(ion, target, e_mid))]
    exc = experiment.excitation
    if exc is not None:
        try:
            kinds.append(("excited", TwoBody(ion, target, e_mid, excitation_mev=exc.energy_mev,
                                              excite="recoil" if exc.excite == "target" else "ejectile")))
        except ValueError:
            pass
    out = []
    # The other layers (the backing and the contaminants): their elastic scattering, a line for each nuclide.
    for j, lay in enumerate(layers):
        for key, share in lay.nuclides.items():
            nuc = data.nuclide(key).name
            if (j == 0 and nuc == target) or share < 0.01 * sum(lay.nuclides.values()):
                continue
            e_layer = float(beam_energy_at(experiment, j, [lay.thickness / 2], layers)[0])
            try:
                tb = TwoBody(ion, nuc, e_layer)
            except ValueError:
                continue
            for particle, species in (("ejectile", ion), ("recoil", nuc)):
                e = np.asarray(tb.at_lab(theta, particle)[0].energy, dtype=float)
                ok = np.isfinite(e) & (e > 0)
                if not ok.any():
                    continue
                e_out = exit_energy(experiment, layers, j, np.full(ok.sum(), lay.thickness / 2), species, e[ok],
                                    dirs[ok])
                resp = array.response(g.name, species)
                cos_i = np.clip(-(dirs[ok] @ g.n), 1e-6, 1.0)
                dep = np.array([float(np.asarray(resp.deposited(np.array([x]),
                                                                float(np.degrees(np.arccos(ci)))))[0])
                                for x, ci in zip(e_out, cos_i)])
                keep = dep > 0
                if keep.any():
                    what = f"scattered {ion} on {nuc}" if particle == "ejectile" else f"{nuc} recoil"
                    out.append({"label": f"{what} ({lay.name})", "group": "elastic", "particle": particle,
                                "layer": lay.name, "rings": [r for r, k in zip(np.array(rings)[ok], keep) if k],
                                "energy": [float(x) for x in dep[keep]]})
    for group, tb in kinds:
        for particle, species in (("ejectile", ion), ("recoil", target)):
            e = np.asarray(tb.at_lab(theta, particle)[0].energy, dtype=float)
            ok = np.isfinite(e) & (e > 0)
            if not ok.any():
                continue
            e_out = exit_energy(experiment, layers, 0, np.full(ok.sum(), layers[0].thickness / 2), species, e[ok],
                                dirs[ok])
            resp = array.response(g.name, species)
            cos_i = np.clip(-(dirs[ok] @ g.n), 1e-6, 1.0)
            dep = np.array([float(np.asarray(resp.deposited(np.array([x]), float(np.degrees(np.arccos(ci)))))[0])
                            for x, ci in zip(e_out, cos_i)])
            keep = dep > 0
            if not keep.any():
                continue
            label = f"{'scattered ' + ion if particle == 'ejectile' else target + ' recoil'}" + (
                ", excited" if group == "excited" else ", elastic")
            out.append({"label": label, "group": group, "particle": particle,
                        "rings": [r for r, k in zip(np.array(rings)[ok], keep) if k],
                        "energy": [float(x) for x in dep[keep]]})
    return out


def group_of(events, experiment) -> np.ndarray:
    """For each particle row: "elastic" or "excited", whichever kinematic line of its ring its measured energy is
    nearer (by energy alone, as an experimentalist would sort them); "elastic" without an excited state."""
    c = events.columns
    out = np.full(len(events), "elastic", dtype=object)
    if experiment.excitation is None or not len(events):
        return out
    for i, name in enumerate(events.detectors):
        m = c["detector"] == i
        if not m.any():
            continue
        lines = kinematic_lines(experiment, name)
        rings = c["segment_i"][m].astype(np.int64)
        meas = c["measured"][m]
        best = {"elastic": np.full(m.sum(), np.inf), "excited": np.full(m.sum(), np.inf)}
        size = int(rings.max()) + 1 if len(rings) else 1
        for line in lines:
            lookup = np.full(max(size, max(line["rings"], default=0) + 1), np.nan)
            lookup[np.asarray(line["rings"], dtype=np.int64)] = line["energy"]
            e = lookup[np.clip(rings, 0, len(lookup) - 1)]
            d = np.abs(meas - e)
            d = np.where(np.isfinite(d), d, np.inf)
            best[line["group"]] = np.minimum(best[line["group"]], d)
        out[np.flatnonzero(m)] = np.where(best["excited"] < best["elastic"], "excited", "elastic")
    return out


# -- spectra ------------------------------------------------------------------------------------------------------

def _edges(lo: float, hi: float, bins: int) -> np.ndarray:
    return np.linspace(lo, hi, int(bins) + 1)


def particle_spectrum(run, detector: str, gate: Optional[Gate] = None, bins: int = 200,
                      range: Optional[tuple] = None, scaled: bool = False) -> dict:
    """Measured-energy spectrum of a particle detector in the run's real part (or, ``scaled``, the whole run):
    {"counts", "edges", "label"}."""
    ev = run.events()
    exp = run.experiment
    m = ev.select(detector)
    if gate is not None:
        m &= gate.mask(ev, exp)
    x = ev["measured"][m]
    if range is None:
        top = float(ev["measured"][ev.select(detector)].max()) * 1.03 if ev.select(detector).any() else 1.0
        range = (0.0, top)
    counts, edges = np.histogram(x, bins=_edges(range[0], range[1], bins))
    counts = counts.astype(float)
    if scaled:
        counts = run.scaled(counts)
    return {"counts": counts, "edges": edges, "label": detector + (f" · {gate.name}" if gate else "")}


def gamma_range(g) -> tuple:
    """The energies (MeV) a γ-ray spectrum is histogrammed over unless a range is given: from 0 to the highest of
    1.3 × the transition energy, the top of the singles spectrum (which holds the room-background and extra lines)
    and the highest energy any crystal recorded in the run, so no line of the run is cut off. A view may zoom in;
    the histogram keeps the whole range."""
    top = 1.3 * g.energy_mev
    if len(g.singles_edges):
        top = max(top, float(g.singles_edges[-1]))
    if len(g):
        top = max(top, 1.05 * float(g["measured"].max()))
    return (0.0, top)


def gamma_spectrum(run, name: str, correction: str = "off", gate: Optional[Gate] = None, mode: str = "coincidence",
                   randoms: str = "shown", addback: bool = True, bins: int = 400, range: Optional[tuple] = None,
                   scaled: bool = False, gamma_gate: Optional[tuple] = None) -> dict:
    """The γ-ray spectrum of a γ-ray detector (or one crystal, "Clover1 A") in the run's real part: in
    coincidence with the particles that pass ``gate`` (all counted particles without one), Doppler-corrected
    (``correction`` "off", "projectile" or "recoil"), with the random coincidences ``randoms`` ("shown",
    "subtracted" or "none"), with or without add-back and suppression (``addback``). ``mode`` "singles" is every
    γ ray the crystal records, from the rates (the run follows only the γ rays in coincidence). ``gamma_gate``
    (MeV, low and high, with the same correction) keeps the γ rays whose cascade has another γ ray in that window:
    a gate on one line selects the lines in coincidence with it (backlog item 71)."""
    exp = run.experiment
    g = run.gammas(plain=not addback)
    if g is None:
        raise ValueError("this run has no γ rays")
    if range is None:
        range = gamma_range(g)
    edges = _edges(range[0], range[1], bins)
    t = g.events.beam_time_s * g.live_fraction
    rng = np.random.default_rng(int(run.summary["seed"]) + 7919)
    if mode == "singles":
        expected = np.zeros(len(edges) - 1)
        centres = (g.singles_edges[:-1] + g.singles_edges[1:]) / 2
        for c in g.crystals:
            if g._in(c, name):
                h, _ = np.histogram(centres, bins=edges, weights=g.singles_spectrum[c.name])
                expected += h
        expected = g.pile_up(expected * t, edges, name, singles=True) / t if abs(edges[0]) < 1e-12 else expected
        counts = rng.poisson(expected * t).astype(float)
        out = {"counts": counts, "edges": edges, "label": f"{name} singles (from the rates)", "expected": expected * t}
    else:
        key = {"off": "measured", "projectile": "corrected_projectile", "recoil": "corrected_recoil"}[correction]
        m = g.select(name)
        dets = list(g.events.detectors)
        if gate is not None:
            passed = gate.mask(g.events, exp)
            m &= passed[g["particle"]]
            dets = gate.detectors or dets
        if gamma_gate is not None:
            lo, hi = sorted(gamma_gate)
            cas = g._cascades()
            hit = g["counted"] & (g[key] >= lo) & (g[key] <= hi)
            gated_cascades, counts_in = np.unique(cas[hit], return_counts=True)
            n_in = dict(zip(gated_cascades.tolist(), counts_in.tolist()))
            # Another γ ray of the cascade in the gate: a γ ray in the window itself needs a second one there.
            others = np.array([n_in.get(int(k), 0) - (1 if h else 0) for k, h in zip(cas, hit)])
            m &= others > 0
        counts, _ = np.histogram(g[key][m], bins=edges, weights=g["weight"][m] * t)
        counts = g.pile_up(counts, edges, name)
        if randoms != "none":
            rnd = np.zeros(len(edges) - 1)
            for d in dets:
                rnd += g.random_spectrum(name, d, edges) * t
            # A ring range or energy window takes a share of each detector's particle singles that the rates do
            # not give: the randoms are those of the whole detectors, an upper limit.
            drawn = rng.poisson(rnd).astype(float)
            counts = counts + drawn - (rnd if randoms == "subtracted" else 0.0)
        out = {"counts": counts, "edges": edges, "label": f"{name} · {correction} correction"
               + (f" · {gate.name}" if gate else "")
               + (f" · γ gate {1e3 * min(gamma_gate):.0f}–{1e3 * max(gamma_gate):.0f} keV" if gamma_gate else "")}
    if scaled:
        out["counts"] = run.scaled(np.maximum(out["counts"], 0.0))
    return out


def grid(run, crystals: bool = False, bins: int = 200, **gamma_options) -> dict:
    """Every spectrum of a run at once: {"particles": {detector: spectrum}, "gammas": {γ detector (or crystal):
    spectrum}}. A source run gives its crystals' spectra under "gammas"."""
    if run.kind == "source":
        sr = run.source()
        return {"particles": {}, "gammas": {n: {"counts": sr.spectra[n], "edges": sr.edges, "label": n}
                                            for n in sr.names()}}
    ev = run.events()
    out = {"particles": {d: particle_spectrum(run, d, bins=bins) for d in ev.detectors}, "gammas": {}}
    g = run.gammas()
    if g is not None:
        names = g.crystal_names if crystals else g.detector_names()
        for n in names:
            out["gammas"][n] = gamma_spectrum(run, n, bins=2 * bins, **gamma_options)
    return out


# -- 2D views -------------------------------------------------------------------------------------------------------

def energy_vs_ring(run, detector: str, bins: int = 150, gate: Optional[Gate] = None) -> dict:
    """Measured energy against ring (or strip along u) of a particle detector, with the kinematic line of each
    group: {"counts" (rings × bins), "rings", "edges", "lines"}. How one sees the elastic and excited groups."""
    ev = run.events()
    exp = run.experiment
    m = ev.select(detector)
    if gate is not None:
        m &= gate.mask(ev, exp)
    rings = np.arange(int(ev["segment_i"][ev.select(detector)].max()) + 1) if ev.select(detector).any() else \
        np.arange(1)
    top = float(ev["measured"][ev.select(detector)].max()) * 1.03 if ev.select(detector).any() else 1.0
    edges = _edges(0.0, top, bins)
    h, _, _ = np.histogram2d(ev["segment_i"][m], ev["measured"][m],
                             bins=[np.arange(len(rings) + 1) - 0.5, edges])
    return {"counts": h, "rings": rings, "edges": edges, "lines": kinematic_lines(exp, detector),
            "detector": detector}


def gamma_vs_crystal(run, correction: str = "off", bins: int = 300, range: Optional[tuple] = None,
                     offset_mm: float = 0.0) -> dict:
    """γ-ray energy against crystal, raw or corrected: {"counts" (crystals × bins), "crystals", "edges"}. How one
    sees a misplaced target: each crystal's corrected line sits elsewhere. ``offset_mm`` corrects with the target
    assumed that far along the beam from where it is (an alignment check)."""
    g = run.gammas()
    if g is None:
        raise ValueError("this run has no γ rays")
    key = {"off": "measured", "projectile": "corrected_projectile", "recoil": "corrected_recoil"}[correction]
    values = g[key]
    if offset_mm and correction != "off":
        from .alignment import with_offset
        from .gamma_events import recorrect

        values = recorrect(g, with_offset(run.experiment, offset_mm))[key]
    if range is None:
        range = (0.9 * g.energy_mev, 1.1 * g.energy_mev)
    edges = _edges(range[0], range[1], bins)
    m = g.select()
    h, _, _ = np.histogram2d(g["crystal"][m], values[m], bins=[np.arange(len(g.crystals) + 1) - 0.5, edges])
    return {"counts": h, "crystals": g.crystal_names, "edges": edges, "correction": correction}


# -- comparing runs -------------------------------------------------------------------------------------------------

def overlay(first, second, spectrum: str, which: str = "particle", **options) -> dict:
    """The same spectrum of two runs, each labelled with its run and setup, and what differs between their
    setups: {"spectra": [spectrum, spectrum], "differences": [...]}. ``which`` is "particle" or "gamma"."""
    fn = particle_spectrum if which == "particle" else gamma_spectrum
    a = fn(first, spectrum, **options)
    opts = dict(options)
    opts.setdefault("range", (a["edges"][0], a["edges"][-1]))
    b = fn(second, spectrum, **opts)
    a["label"] = f"run {first.number}: {a['label']} ({first.describe().split(';')[0]})"
    b["label"] = f"run {second.number}: {b['label']} ({second.describe().split(';')[0]})"
    return {"spectra": [a, b], "differences": first.changes(second.setup)}


# -- exports --------------------------------------------------------------------------------------------------------

def spectrum_csv(spectrum: dict) -> str:
    """A spectrum as CSV: bin low edge, bin high edge, counts."""
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["low", "high", "counts"])
    e = spectrum["edges"]
    for lo, hi, c in zip(e[:-1], e[1:], spectrum["counts"]):
        w.writerow([f"{lo:.6g}", f"{hi:.6g}", f"{c:.6g}"])
    return buf.getvalue()


def cut_expression(gate: Gate, detectors: list) -> str:
    """A gate as a ROOT cut on the ``events`` tree (the elastic/excited group has no simple cut: it is noted)."""
    parts = []
    if gate.detectors is not None:
        parts.append("(" + " || ".join(f"detector == {detectors.index(d)}" for d in gate.detectors
                                       if d in detectors) + ")")
    if gate.rings is not None:
        parts.append("(" + " || ".join(f"segment_i == {r}" for r in gate.rings) + ")")
    if gate.energy is not None:
        parts.append(f"measured >= {gate.energy[0]:.6g} && measured <= {gate.energy[1]:.6g}")
    expr = " && ".join(parts) or "1"
    if gate.group != "any":
        expr += f"  /* and the {gate.group} group by energy, not expressible as a cut */"
    return expr


def write_root(run, path: Union[str, Path], gates: tuple = (), bins: int = 400) -> Path:
    """The run's real part as a ROOT file (:func:`physim.nuclear.rootio.write_root`), with each gate as a cut
    (``cut_<name>``) and the run's summary."""
    from . import rootio

    path = Path(path)
    ev = run.events()
    rootio.write_root(ev, path, run.experiment, bins=bins, gammas=run.gammas())
    import uproot

    with uproot.update(path) as f:
        for gate in gates:
            f[f"cut_{_safe(gate.name)}"] = cut_expression(gate, list(ev.detectors))
        f["run"] = json.dumps({k: v for k, v in run.summary.items() if k not in ("counters_raw", "generator")},
                              default=str)
    return path


def _safe(name: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in name)


# -- gates kept with the experiment -------------------------------------------------------------------------------

@dataclass
class GateBook:
    """The named gates of an experiment, saved in its folder as ``gates.json`` (in memory without a folder)."""

    folder: Optional[Path] = None
    gates: dict = field(default_factory=dict)

    @classmethod
    def load(cls, folder: Optional[Path]) -> GateBook:
        book = cls(folder)
        p = book._path()
        if p is not None and p.is_file():
            for d in json.loads(p.read_text(encoding="utf-8")):
                g = Gate.from_dict(d)
                book.gates[g.name] = g
        return book

    def _path(self) -> Optional[Path]:
        return Path(self.folder) / "gates.json" if self.folder is not None else None

    def save(self) -> None:
        p = self._path()
        if p is not None:
            p.write_text(json.dumps([g.as_dict() for g in self.gates.values()], indent=1), encoding="utf-8")

    def put(self, gate: Gate) -> None:
        self.gates[gate.name] = gate
        self.save()

    def remove(self, name: str) -> None:
        self.gates.pop(name, None)
        self.save()

    def __getitem__(self, name: str) -> Gate:
        if name not in self.gates:
            raise KeyError(f"no gate named {name!r}; the gates are {list(self.gates)}")
        return self.gates[name]


__all__ = ["CORRECTIONS", "RANDOMS", "Gate", "GateBook", "correction_reason", "cut_expression", "energy_vs_ring",
           "gamma_spectrum", "gamma_vs_crystal", "grid", "group_of", "kinematic_lines", "overlay",
           "particle_spectrum", "spectrum_csv", "write_root"]
