"""Simulated events and spectra as a ROOT file, written with ``uproot`` (no ROOT installation needed).

::

    from physim.nuclear import Experiment
    from physim.nuclear.events import simulate
    from physim.nuclear.rootio import write_root

    exp = Experiment.example("alpha_on_gold")
    write_root(simulate(exp, 1_000_000, seed=1), "alpha.root", exp)

The file holds:

- ``events``: a TTree with one entry per particle that reached a detector, one branch per column of
  :class:`~physim.nuclear.events.Events` (``detector``, ``segment_i``, ``measured``, ``theta``, ``weight``, ...);
  units as in :data:`physim.nuclear.events.COLUMNS`. ``weight`` is a rate (1/s); multiply by the beam time for
  counts.
- ``spectrum_<detector>``: a TH1D of the measured energy (MeV) of the counted particles, in counts per bin for the
  planned beam time, with the Monte Carlo statistical error of each bin (Σw²) stored.
- ``setup``: the setup file (TOML) the events came from, as a TObjString; ``info``: seed, number of events, beam
  time, detector and channel names, as a TObjString.

In ROOT: ``events->Draw("measured", "weight*(detector==2 && counted)")``; ``spectrum_A45->Fit("gaus")``.

Requires ``uproot`` (``pip install physim-engine[root]``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Union

import numpy as np


def _uproot():
    try:
        import uproot
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise ImportError("writing ROOT files needs uproot: pip install uproot") from e
    return uproot


def _th1d(name: str, title: str, counts: np.ndarray, sumw2: np.ndarray, edges: np.ndarray, entries: int,
          xs: np.ndarray, ws: np.ndarray):
    from uproot.writing.identify import to_TAxis, to_TH1x

    # ROOT stores underflow and overflow bins at either end.
    data = np.concatenate([[0.0], counts, [0.0]])
    sw2 = np.concatenate([[0.0], sumw2, [0.0]])
    axis = to_TAxis("xaxis", "measured energy (MeV)", len(counts), float(edges[0]), float(edges[-1]))
    return to_TH1x(name, title, data, float(entries), float(ws.sum()), float((ws**2).sum()),
                   float((ws * xs).sum()), float((ws * xs**2).sum()), sw2, axis)


def write_root(events, path: Union[str, Path], experiment=None, bins: int = 400,
               range_mev: Optional[tuple] = None, gammas=None) -> Path:
    """Write ``events`` (from :func:`~physim.nuclear.events.simulate`) to a ROOT file; ``experiment`` is stored as
    its setup file. With ``gammas`` (from :func:`~physim.nuclear.gamma_events.simulate_gammas`) the file also
    holds a ``gammas`` tree, one entry per γ ray in coincidence with a particle (``particle`` is the row of
    ``events`` it belongs to), and ``gamma_<crystal>`` histograms of the Doppler-corrected energy. Returns the
    path."""
    uproot = _uproot()
    path = Path(path)
    cols = events.columns
    tree = {}
    for name, arr in cols.items():
        a = np.asarray(arr)
        tree[name] = a.astype(np.int32) if a.dtype in (np.uint32,) else a
    t = events.beam_time_s
    with uproot.recreate(path) as f:
        f["events"] = tree
        for name in events.detectors:
            m = events.select(name)
            x = np.asarray(cols["measured"])[m]
            w = np.asarray(cols["weight"])[m] * t
            lo, hi = range_mev if range_mev is not None else (0.0, float(x.max()) * 1.03 if len(x) else 1.0)
            counts, edges = np.histogram(x, bins=bins, range=(lo, hi), weights=w)
            sumw2, _ = np.histogram(x, bins=bins, range=(lo, hi), weights=w**2)
            inside = (x >= lo) & (x < hi)
            f[f"spectrum_{name}"] = _th1d(f"spectrum_{name}", f"{name}: measured energy, counts in the run",
                                          counts, sumw2, edges, int(m.sum()), x[inside], w[inside])
        if gammas is not None:
            g = {name: (np.asarray(a).astype(np.int32) if np.asarray(a).dtype in (np.uint32, np.int64)
                        else np.asarray(a)) for name, a in gammas.columns.items()}
            f["gammas"] = g
            key = "corrected_recoil" if gammas.emitter == experiment_target(experiment, gammas) else \
                "corrected_projectile"
            for name in gammas.detector_names():
                counts, edges = gammas.spectrum(name, corrected=key.split("_")[1], bins=bins, randoms=False)
                f[f"gamma_{name}"] = _th1d(f"gamma_{name}", f"{name}: Doppler-corrected γ-ray energy, counts in "
                                           "the run", counts, counts, edges, int(gammas.select(name).sum()),
                                           np.zeros(0), np.zeros(0))
        if experiment is not None:
            f["setup"] = experiment.to_toml()
        info = {"seed": events.seed, "n_events": events.n_events, "beam_time_s": t,
                "detectors": list(events.detectors), "channels": list(events.channels),
                "min_energy_mev": events.min_energy, "weight_unit": "1/s"}
        if gammas is not None:
            info.update(crystals=gammas.crystal_names, gamma_detectors=gammas.detector_names(),
                        window_s=gammas.window_s, live_fraction=gammas.live_fraction)
        f["info"] = json.dumps(info)
    return path


def experiment_target(experiment, gammas) -> str:
    """The nuclide the γ rays come from when the target is excited (for choosing the corrected spectrum)."""
    return gammas.emitter if experiment is not None and experiment.excitation is not None \
        and experiment.excitation.excite == "target" else ""


def read_root(path: Union[str, Path]) -> dict:
    """Read a file written by :func:`write_root`: ``{"events": {column: array}, "spectra": {detector: (counts,
    errors, edges)}, "setup": text or None, "info": dict}``."""
    uproot = _uproot()
    with uproot.open(path) as f:
        events = f["events"].arrays(library="np")
        info = json.loads(str(f["info"]))
        spectra = {}
        for name in info["detectors"]:
            h = f[f"spectrum_{name}"]
            spectra[name] = (h.values(), h.errors(), h.axis().edges())
        setup = str(f["setup"]) if "setup" in f else None
        gammas = f["gammas"].arrays(library="np") if "gammas" in f else None
    return {"events": events, "spectra": spectra, "setup": setup, "info": info, "gammas": gammas}


__all__ = ["read_root", "write_root"]
