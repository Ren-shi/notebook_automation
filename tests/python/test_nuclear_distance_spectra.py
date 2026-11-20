"""dataviews.distance_spectra (lifetime wave 3A): a plunger's distance series as spectra per distance and per ring
of γ-ray detectors, to the frozen shape {edges_kev, by_distance: {distance: {theta_deg, counts}}}. The runs are
synthetic, with the surfaces of a particle run's GammaEvents and a fusion run's FusionGammas."""

from types import SimpleNamespace

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import dataviews as dv
from physim.nuclear.gamma_events import crystals_of
from physim.nuclear.runs import DistanceSeries


class _CoulexGammas:
    live_fraction = 0.5

    def __init__(self, measured_kev, crystal, crystals, weight=1.0):
        n = len(measured_kev)
        m = 1e-3 * np.asarray(measured_kev, dtype=float)
        self.columns = {"measured": m, "crystal": np.asarray(crystal), "weight": np.full(n, weight),
                        "corrected_projectile": m - 0.010, "corrected_recoil": m + 0.010}
        self.crystals = crystals
        self.events = SimpleNamespace(beam_time_s=4.0)

    def select(self):
        return np.ones(len(self.columns["measured"]), dtype=bool)

    def __getitem__(self, key):
        return self.columns[key]


class _FusionGammas:
    def __init__(self, measured_kev, crystal, names):
        m = 1e-3 * np.asarray(measured_kev, dtype=float)
        self.columns = {"measured": m, "crystal": np.asarray(crystal), "corrected": m + 0.005}
        self.crystals = names

    def __len__(self):
        return len(self.columns["measured"])

    def __getitem__(self, key):
        return self.columns[key]


class _Run:
    def __init__(self, experiment, gammas, label, kind="beam", real_s=10.0, duration_s=40.0):
        self.experiment, self._g, self.label, self.kind = experiment, gammas, label, kind
        self.real_s, self.duration_s = real_s, duration_s

    def gammas(self, plain=False):
        return self._g


@pytest.fixture(scope="module")
def coulex():
    return Experiment.example("coulex_ru98_plunger")


def _ring_of(exp):
    crystals = crystals_of(exp, plain=True)[0]
    return crystals, [round(float(exp.gamma_detectors[c.detector].theta.to("deg")), 1) for c in crystals]


def test_spectra_by_distance_and_ring_have_the_frozen_shape(coulex):
    crystals, theta = _ring_of(coulex)
    first = {t: theta.index(t) for t in sorted(set(theta))}
    # Distance 10: two γ rays in the 41.5° ring at 660 keV, one at 90° at 652 keV; distance 50: one at 138.5°.
    a = _Run(coulex, _CoulexGammas([660.2, 660.4, 652.1], [first[41.5], first[41.5], first[90.0]], crystals), "a")
    b = _Run(coulex, _CoulexGammas([640.3], [first[138.5]], crystals), "b")
    series = DistanceSeries(coulex.to_dict(), {10.0: [a], 50.0: [b]})
    out = dv.distance_spectra(series)
    assert set(out) >= {"edges_kev", "by_distance"} and out["doppler"] is True
    edges = out["edges_kev"]
    assert np.all(np.diff(edges) == pytest.approx(dv.DISTANCE_BIN_KEV)) and edges[-1] >= 660.4
    assert sorted(out["by_distance"]) == [10.0, 50.0]
    ten = out["by_distance"][10.0]
    assert ten["theta_deg"] == [41.5, 90.0, 138.5] == dv.gamma_rings(coulex)
    assert ten["counts"].shape == (3, len(edges) - 1)
    # Each γ ray counts weight × beam time × live fraction (here 1 × 4 s × 0.5 = 2) in its ring and bin.
    assert ten["counts"][0].sum() == pytest.approx(4.0) and ten["counts"][1].sum() == pytest.approx(2.0)
    assert ten["counts"][2].sum() == 0
    k = np.searchsorted(edges, 660.2, side="right") - 1
    assert ten["counts"][0][k] == pytest.approx(4.0)
    assert out["by_distance"][50.0]["counts"][2].sum() == pytest.approx(2.0)
    assert ten["scale"] == pytest.approx(4.0) and ten["runs"] == ["a"]


def test_runs_at_one_distance_add_and_the_correction_moves_the_energies(coulex):
    crystals, theta = _ring_of(coulex)
    c0 = theta.index(41.5)
    runs = [_Run(coulex, _CoulexGammas([700.0], [c0], crystals), "a"),
            _Run(coulex, _CoulexGammas([700.0], [c0], crystals), "b")]
    series = DistanceSeries(coulex.to_dict(), {20.0: runs})
    edges = np.arange(600.0, 800.0, 1.0)
    out = dv.distance_spectra(series, edges_kev=edges)
    counts = out["by_distance"][20.0]["counts"][0]
    assert counts.sum() == pytest.approx(4.0) and counts[np.searchsorted(edges, 700.0, side="right") - 1] == 4.0
    assert out["by_distance"][20.0]["runs"] == ["a", "b"]
    # 98Ru is the projectile: Doppler-corrected energies are the projectile's column (10 keV below here).
    corrected = dv.distance_spectra(series, edges_kev=edges, doppler=False)["by_distance"][20.0]["counts"][0]
    assert corrected[np.searchsorted(edges, 690.0, side="right") - 1] == 4.0
    with pytest.raises(ValueError):
        dv.distance_spectra(series, edges_kev=[1.0, 0.5])


def test_a_fusion_series_reads_its_crystals_by_name():
    exp = Experiment.example("fusion_cr50_plunger")
    crystals, theta = _ring_of(exp)
    names = [c.name for c in crystals]
    back, fwd = theta.index(142.3), theta.index(45.0)
    g = _FusionGammas([760.0, 805.0, 805.5], [back, fwd, fwd], names)
    series = DistanceSeries(exp.to_dict(), {17.0: [_Run(exp, g, "f", kind="fusion")]})
    out = dv.distance_spectra(series)
    assert out["kind"] == "fusion"
    d = out["by_distance"][17.0]
    assert d["theta_deg"] == [45.0, 142.3]
    assert d["counts"][0].sum() == 2 and d["counts"][1].sum() == 1  # one count per γ ray


def test_no_gamma_detectors_is_refused(coulex):
    d = coulex.to_dict()
    d["gamma_detectors"] = []
    d["reaction"] = {"type": "elastic"}
    d.pop("lifetime")
    exp = Experiment.from_dict(d)
    with pytest.raises(ValueError):
        dv.distance_spectra(DistanceSeries(exp.to_dict(), {}))
