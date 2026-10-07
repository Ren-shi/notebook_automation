"""Summing and pile-up (backlog item 72): a ⁶⁰Co source close to a crystal shows the 2.5 MeV sum peak with the area
the efficiencies predict, the efficiency points are corrected for summing and say so, and pile-up takes counts
from the peaks into a shoulder above them."""

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.planner import Planner
from physim.nuclear.response import pile_up, source_run


def close_germanium(distance="40 mm", shaping=None) -> Experiment:
    p = Planner.example("coulex_ni58")
    d = p.draft
    d["gamma_detectors"] = [dict(d["gamma_detectors"][0], distance=distance)]
    if shaping:
        d["run"]["shaping_time"] = shaping
    return Experiment.from_dict(d)


def test_the_60co_sum_peak_has_the_area_the_efficiencies_predict():
    exp = close_germanium()
    sr = source_run(exp, "60Co", "50 kBq", "2 h", seed=3, bin_kev=1.0)
    name = sr.names()[0]
    (peak,) = sr.sum_peaks[name]
    assert peak["energy_mev"] == pytest.approx(2.505720, abs=1e-6)
    area, unc = sr.peak_area(name, peak["energy_mev"])
    assert abs(area - peak["expected"]) < 4 * unc + 0.03 * peak["expected"], (area, unc, peak["expected"])
    assert area > 10 * unc  # clearly there
    plain = source_run(exp, "60Co", "50 kBq", "2 h", seed=3, bin_kev=1.0, summing=False)
    assert not plain.sum_peaks[name] and sr.spectra[name].sum() < plain.spectra[name].sum() + 1e9


def test_the_efficiency_points_are_corrected_for_summing_and_say_so():
    exp = close_germanium()
    sr = source_run(exp, "60Co", "50 kBq", "2 h", seed=3, bin_kev=1.0)
    points = sr.efficiency_points(sr.names()[0])
    assert len(points) == 2
    for p in points:
        assert p["summing_correction"] > 1.0
        assert p["uncorrected"] < p["efficiency"]
        assert abs(p["efficiency"] - p["true"]) < 4 * p["uncertainty"] + 0.01 * p["true"], p
        # Without the correction the point would be low by the share summing took.
        assert p["uncorrected"] == pytest.approx(p["efficiency"] / p["summing_correction"])
    far = source_run(close_germanium("250 mm"), "60Co", "50 kBq", "2 h", seed=3, bin_kev=1.0)
    near_corr = max(p["summing_correction"] for p in points)
    far_corr = max(p["summing_correction"] for p in far.efficiency_points(far.names()[0]))
    assert 1.0 < far_corr < near_corr  # summing falls with distance


def test_pile_up_is_a_loss_and_a_shoulder():
    counts = np.zeros(400)
    counts[100] = 1e5  # one line
    piled, p = pile_up(counts, rate=20_000.0, shaping_time_s=3e-6)
    assert p == pytest.approx(1 - np.exp(-2 * 3e-6 * 2e4))
    assert piled[100] == pytest.approx((1 - p) * 1e5)
    assert piled[200] == pytest.approx(0.5 * p * 1e5)       # two signals of the line, piled
    assert piled.sum() == pytest.approx(1e5 * (1 - p / 2))  # a pair counts once
    # In a source run: the shaping time takes counts out of the peaks.
    exp = close_germanium("60 mm")
    calm = source_run(exp, "60Co", "2 MBq", "10 min", seed=1, bin_kev=1.0)
    busy = source_run(close_germanium("60 mm", shaping="6 us"), "60Co", "2 MBq", "10 min", seed=1, bin_kev=1.0)
    name = calm.names()[0]
    assert busy.pile_up[name] > 0.01
    a_calm, _ = calm.peak_area(name, 1.332492)
    a_busy, _ = busy.peak_area(name, 1.332492)
    assert a_busy < a_calm * (1 - 0.5 * busy.pile_up[name])
