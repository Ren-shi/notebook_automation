"""The particle–γ correlation table's points over a particle detector (:func:`physim.nuclear.gamma._table_points`): cut and
crowded where a particle's largest lab angle crosses the face, equally spaced round a ring; its factors by segment and
kind of particle (:func:`physim.nuclear.gamma.segment_correlation`); and the analysis's use of them."""

import json
import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.detectors import Array
from physim.nuclear.gamma import _crossings, _line_nodes, _table_points, correlation_table, segment_correlation


def test_the_points_crowded_to_a_cut_integrate_an_inverse_square_root_to_a_part_in_a_million():
    # ∫_0^1 (0.6 − s)^(−1/2) ds over s < 0.6 = 2√0.6, with the share beyond the cut 0, as a particle's cross section
    # per lab solid angle up to its largest angle.
    s, w = _line_nodes(0.0, 1.0, [0.6], 8)
    f = np.where(s < 0.6, 1.0 / np.sqrt(np.abs(0.6 - s)), 0.0)
    assert abs(float(w @ f) / (2 * math.sqrt(0.6)) - 1) < 1e-6
    # Without a cut, the plain Gauss points.
    s, w = _line_nodes(-1.0, 1.0, [], 5)
    x, wx = np.polynomial.legendre.leggauss(5)
    assert np.allclose(s, x) and np.allclose(w, wx)


def test_a_line_crosses_the_cone_of_an_angle_where_its_polar_angle_is_that():
    start, along = np.array([0.0, 0.0, 25.0]), np.array([1.0, 0.0, 0.0])
    (s,) = _crossings(start, along, 1.5, 14.3, [14.0])
    assert math.degrees(math.atan2(s, 25.0)) == pytest.approx(14.0, abs=1e-9)
    assert _crossings(start, along, 1.5, 14.3, [40.0]) == []  # beyond the face
    assert _crossings(start, along, 1.5, 14.3, [180.0 - 14.0]) == []  # the mirror cone, behind the target


def test_the_points_cover_the_solid_angle_of_the_face():
    exp = Experiment.example("coulex_ni58")
    for g in Array.from_experiment(exp):
        whole = g.solid_angle() * 1e-3
        assert _table_points(g, 8)[1].sum() == pytest.approx(whole, rel=1e-6)
        # Cut where no particle turns back, the pieces (5 points each, crowded to their ends) still cover the face.
        cut = (g.theta_range()[0] * 0.6 + g.theta_range()[1] * 0.4,)
        assert _table_points(g, 8, cut)[1].sum() == pytest.approx(whole, rel=1e-4)
        assert sum(_table_points(g, 3, (), seg)[1].sum() for seg in g.segments) == pytest.approx(whole, rel=1e-6)


@pytest.fixture(scope="module")
def ru98():
    return Experiment.example("coulex_ru98_plunger")


def test_the_table_has_converged_at_its_order_where_the_scattered_beam_turns_back_on_the_face(ru98):
    # 98Ru on 24Mg: the scattered 98Ru goes no further than 14° and the detector at 0° reaches 30°. Before the points
    # were cut there, order 8 was 7.6 % off per crystal.
    rows8 = correlation_table(ru98)
    rows24 = correlation_table(ru98, order=24)
    worst = max(abs(a["factor"] / b["factor"] - 1) for a, b in zip(rows8, rows24))
    assert worst < 1e-3


def test_a_segment_holding_the_whole_face_has_the_table_s_factor():
    # coulex_ni58's annular detector made one segment: only the scattered 16O reaches it, so the factor of that kind
    # in that segment is the table's.
    raw = Experiment.example("coulex_ni58").to_dict()
    for d in raw["detectors"]:
        if d["name"] == "CD":  # backwards: only the scattered beam comes there
            d["rings"], d["sectors"] = 1, 1
    exp = Experiment.from_dict(json.loads(json.dumps(raw)))
    names, f = segment_correlation(exp, [("CD", (0, 0), False), ("CD", (0, 0), True)])
    table = {r["gamma_detector"]: r["detector_factor"] for r in correlation_table(exp) if r["particle_detector"] == "CD"}
    assert np.allclose(f[0], [table[n] for n in names], rtol=1e-12)
    assert np.isnan(f[1]).all()  # no recoil goes backwards


def test_the_factor_of_a_segment_and_kind_is_its_own():
    # In the strips of coulex_ni58 both the scattered 16O and the 58Ni recoils arrive; their factors differ, and the
    # detector's lies between them.
    exp = Experiment.example("coulex_ni58")
    names, f = segment_correlation(exp, [("DSSD-L", (8, 8), False), ("DSSD-L", (8, 8), True)])
    assert np.all(np.isfinite(f))
    assert np.max(np.abs(f[0] / f[1] - 1)) > 0.02


def test_a_weighted_sample_keeps_the_gamma_rays_of_counted_particles_only():
    from physim.nuclear.gamma_events import simulate_gammas

    exp = Experiment.example("coulex_ni58")
    g = simulate_gammas(exp, 20_000, seed=3)
    assert g.events.columns["counted"][g["particle"]].all()
    assert (~g.events.columns["counted"]).any()  # there were particles below threshold


def test_a_window_on_the_particle_energy_divides_by_the_share_it_keeps():
    from physim.nuclear.analysis import Analysis, Settings, _narrower
    from physim.nuclear.gamma_events import simulate_gammas

    assert _narrower(Settings()) and _narrower(Settings(particle_gate="all", particle_energy=(8.0, 16.0)))
    assert not _narrower(Settings(particle_gate="all"))
    exp = Experiment.example("coulex_ni58")
    a = Analysis(exp, simulate_gammas(exp, 40_000, seed=5))
    r = a.run(Settings(particle_gate="all", particle_energy=(8.0, 16.0)))
    step = next(x["text"] for x in r.steps if x["step"] == "yield")
    assert 0.0 < a._gate_efficiency < 0.5 and f"gate efficiency {a._gate_efficiency:.3f}" in step
