"""Contaminant reactions, beam spot and halo (backlog item 74): a carbon contaminant gives its elastic line at the
right energy in every ring, and a 3 mm spot broadens the Doppler-corrected γ peak by what the geometry gives."""

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import dataviews as dv
from physim.nuclear.detectors import Array
from physim.nuclear.events import simulate
from physim.nuclear.gamma import excitation_of
from physim.nuclear.gamma_events import _reconstruct, _unit, simulate_gammas
from physim.nuclear.kinematics import TwoBody
from physim.nuclear.planner import Planner
from physim.nuclear.rates import stack


def with_carbon() -> Experiment:
    p = Planner.example("coulex_ni58")
    d = p.draft
    d["target"]["thickness"] = "0.2 mg/cm2"
    d["target"]["contaminants"] = [["C", "20 ug/cm2"]]
    return Experiment.from_dict(d)


def test_a_carbon_contaminant_has_its_elastic_line_in_every_ring():
    exp = with_carbon()
    layers = stack(exp)
    assert [lay.name for lay in layers] == ["target", "C contaminant"]
    ev = simulate(exp, 400_000, 3)
    carbon = ev.channels.index("12C (C contaminant)")
    for det in ("DSSD-L", "DSSD-R"):
        lines = dv.kinematic_lines(exp, det)
        line = next(x for x in lines if x["label"] == "scattered 16O on 12C (C contaminant)")
        at = dict(zip(line["rings"], line["energy"]))
        m = ev.select(det) & (ev["channel"] == carbon) & ~ev["recoil"]
        seen = sorted(set(ev["segment_i"][m].tolist()))
        assert len(seen) >= 10, det
        for ring in seen:
            here = m & (ev["segment_i"] == ring)
            if here.sum() < 20:
                continue
            median = float(np.median(ev["measured"][here]))
            # The strip spans a range of angles: the line is at its middle, within the spread across the strip.
            spread = float(np.std(ev["measured"][here]))
            assert abs(median - at[ring]) < max(3 * spread / np.sqrt(here.sum()) + 0.03, 0.05), (det, ring)


def _spot_setup(spot):
    p = Planner.example("coulex_ni58")
    d = p.draft
    d["detectors"] = [x for x in d["detectors"] if x["name"] == "CD"]
    d["gamma_detectors"] = [gd for gd in d["gamma_detectors"] if gd["name"] in ("Ge0", "Ge45")]
    for gd in d["gamma_detectors"]:
        gd["distance"] = "60 mm"
    d["target"]["thickness"] = "0.1 mg/cm2"
    d["reaction"]["emission"] = "isotropic"
    if spot:
        d["beam"]["spot_size"] = spot
    else:
        d["beam"].pop("spot_size", None)
    return Experiment.from_dict(d)


def _width(g) -> float:
    x = g["corrected_recoil"]
    m = g["counted"] & (np.abs(x - 1.4542) < 0.03) & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-9)
    return float(np.sqrt(np.cov(x[m], aweights=g["weight"][m]))), int(m.sum())


def test_a_3mm_spot_broadens_the_corrected_peak_by_what_the_geometry_gives():
    sharp, spot = _spot_setup(None), _spot_setup("3 mm")
    g0 = simulate_gammas(sharp, 600_000, 7, gammas_per_event=10)
    g3 = simulate_gammas(spot, 600_000, 7, gammas_per_event=10)
    w0, n0 = _width(g0)
    w3, n3 = _width(g3)
    assert w3 > w0
    measured = np.sqrt(w3**2 - w0**2)

    # What the geometry gives: the same γ rays, corrected as if they had come from points of the spot while the
    # correction assumes the target's centre.
    sigma = 3.0 / 2.3548
    rng = np.random.default_rng(1)
    c = g0.events.columns
    full = g0["counted"] & (np.abs(g0["corrected_recoil"] - 1.4542) < 0.03)
    rows, crystals, measured_e = g0["particle"][full], g0["crystal"][full], g0["measured"][full]
    offset = np.c_[rng.normal(0, sigma, len(rows)), rng.normal(0, sigma, len(rows)), np.zeros(len(rows))]
    array = Array.from_experiment(sharp)
    seg = np.array([array.geometries[int(c["detector"][r])].segment_centre(
        (int(c["segment_i"][r]), int(c["segment_j"][r])), weighted=True) for r in rows])
    centres = np.array([np.array(cn, dtype=float) for gd in sharp.gamma_detectors for _, cn, _ in gd.elements()])
    ex = excitation_of(sharp)
    tb = TwoBody(ex.beam, ex.target, ex.beam_energy, excitation_mev=1.4542, excite="recoil")
    layers = stack(sharp)

    def corrected(shift):
        p_dir = _unit(seg - shift)
        th = np.degrees(np.arccos(np.clip(p_dir[:, 2], -1, 1)))
        d_em, b_em = _reconstruct(sharp, tb, layers, p_dir, th, c["recoil"][rows], "recoil", ex.beam_energy)
        cos_a = np.einsum("ij,ij->i", _unit(centres[crystals] - shift), d_em)
        return measured_e * (1 - b_em * cos_a) / np.sqrt(1 - b_em**2)

    predicted = float(np.std(corrected(offset) - corrected(np.zeros_like(offset))))
    # The statistical error of a width from n counts is about width / √(2n); on the difference of squares it adds.
    err = np.sqrt((w3**2 / np.sqrt(2 * n3)) ** 2 + (w0**2 / np.sqrt(2 * n0)) ** 2) / measured
    assert measured == pytest.approx(predicted, abs=3 * err + 0.1 * predicted), (measured, predicted, err)


def test_the_halo_spreads_the_beam_over_the_frame():
    p = Planner.example("alpha_on_gold")
    d = p.draft
    d["beam"]["halo_fraction"] = "20 %"
    d["beam"]["halo_radius"] = "8 mm"
    d["beam"].pop("spot_size", None)
    ev = simulate(Experiment.from_dict(d), 200_000, 2)
    r = np.hypot(ev["x"], ev["y"])
    assert np.all(r <= 8.0 + 1e-9)
    halo = r > 1e-12
    assert halo.mean() == pytest.approx(0.2, abs=0.03)
    # The halo is even over the disc: half its points lie outside r/√2 of the radius.
    assert np.mean(r[halo] > 8.0 / np.sqrt(2)) == pytest.approx(0.5, abs=0.05)
    plain = simulate(Planner.example("alpha_on_gold").experiment, 50_000, 2)
    assert not np.any(plain["x"])  # no spot, no halo: on the axis
