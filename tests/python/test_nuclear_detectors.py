import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.detectors import Array, Geometry, Response, direction, exit_path


def rect_corner(x, y, d):
    """Solid angle of the rectangle [0, x] × [0, y] in a plane at distance d, seen from the foot of the normal."""
    return math.atan(x * y / (d * math.sqrt(x * x + y * y + d * d)))


# -- solid angles -----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("d, radius", [(100.0, 2.5), (40.0, 30.0), (10.0, 50.0)])
def test_disc_on_axis(d, radius):
    g = Geometry("disc", "circle", [0, 0, d], [0, 0, -1], radius=radius)
    exact = 2 * math.pi * (1 - d / math.hypot(d, radius)) * 1e3
    assert g.solid_angle() == pytest.approx(exact, rel=1e-10)


@pytest.mark.parametrize("d, a, b", [(150.0, 25.0, 25.0), (50.0, 40.0, 10.0), (20.0, 50.0, 50.0)])
def test_rectangle_facing_the_target(d, a, b):
    g = Geometry("r", "rectangle", direction(37.0, 20.0) * d, -direction(37.0, 20.0), width=2 * a, height=2 * b)
    exact = 4 * rect_corner(a, b, d) * 1e3
    assert g.solid_angle() == pytest.approx(exact, rel=1e-10)


def test_off_centre_rectangle():
    """A rectangle in the plane z = d, facing straight back along −z, not centred on the axis: the exact answer adds
    and subtracts rectangles with a corner at the foot of the normal."""
    d, cx, cy, w, h = 80.0, 60.0, -15.0, 40.0, 30.0
    g = Geometry("r", "rectangle", [cx, cy, d], [0, 0, -1], width=w, height=h)
    x1, x2, y1, y2 = cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2
    exact = (math.copysign(1, x2) * math.copysign(1, y2) * rect_corner(abs(x2), abs(y2), d)
             - math.copysign(1, x1) * math.copysign(1, y2) * rect_corner(abs(x1), abs(y2), d)
             - math.copysign(1, x2) * math.copysign(1, y1) * rect_corner(abs(x2), abs(y1), d)
             + math.copysign(1, x1) * math.copysign(1, y1) * rect_corner(abs(x1), abs(y1), d)) * 1e3
    assert g.solid_angle() == pytest.approx(exact, rel=1e-10)


def test_annulus_on_axis_and_its_segments():
    g = Geometry("cd", "annular", [0, 0, -40], [0, 0, 1], inner_radius=9.0, outer_radius=41.0, rings=16, sectors=24)
    exact = 2 * math.pi * (40 / math.hypot(40, 9) - 40 / math.hypot(40, 41)) * 1e3
    assert g.solid_angle() == pytest.approx(exact, rel=1e-10)
    segs = g.segment_solid_angles()
    assert len(segs) == 16 * 24
    assert sum(segs.values()) == pytest.approx(exact, rel=1e-10)
    # Ring k on axis: 2π(cos θ_in − cos θ_out) / sectors for each sector.
    r = 9.0 + 2.0 * np.arange(17)
    for k in (0, 7, 15):
        ring = 2 * math.pi * (40 / math.hypot(40, r[k]) - 40 / math.hypot(40, r[k + 1])) * 1e3 / 24
        assert segs[(k, 5)] == pytest.approx(ring, rel=1e-9)


def test_monte_carlo_agrees_for_tilted_faces():
    """Faces that do not look at the target square-on, checked by an independent Monte Carlo count."""
    cases = [
        Geometry("tilted", "rectangle", direction(60, 30) * 120, -direction(40, 30), width=50, height=30),
        Geometry("cd", "annular", [5, 0, -60], [0.2, 0.1, 1], inner_radius=10, outer_radius=45, rotation_deg=17),
        Geometry("disc", "circle", direction(100, -45) * 70, -direction(110, -40), radius=12),
    ]
    for g in cases:
        mc, err = g.solid_angle_monte_carlo(400_000, seed=3)
        assert abs(mc - g.solid_angle()) < 4 * err, (g.name, mc, err, g.solid_angle())


def test_strips_add_up_and_are_hit_where_expected():
    exp = Experiment.example("oxygen_on_lead_array")
    d = Array.from_experiment(exp)["DSSD1"]
    segs = d.segment_solid_angles()
    assert len(segs) == 256 and sum(segs.values()) == pytest.approx(d.solid_angle(), rel=1e-12)
    for seg in [(0, 0), (3, 12), (15, 15)]:
        du, dv = d.width / 16, d.height / 16
        local = (-d.width / 2 + (seg[0] + 0.5) * du, -d.height / 2 + (seg[1] + 0.5) * dv)
        p = d.centre + local[0] * d.u + local[1] * d.v
        hit = d.hit(p)
        assert hit.segment == seg and hit.distance == pytest.approx(np.linalg.norm(p))
        assert hit.local == pytest.approx(local)
    assert d.hit(direction(150, 0)) is None  # behind the target
    assert d.hit(-d.centre) is None  # opposite direction


def test_angular_coverage():
    exp = Experiment.example("oxygen_on_lead_array")
    array = Array.from_experiment(exp)
    d = array["DSSD1"]  # 50 × 50 mm at 150 mm, θ = 35°, φ = 0, facing the target
    lo, hi = d.theta_range()
    assert lo == pytest.approx(35 - math.degrees(math.atan(25 / 150)), abs=1e-9)
    # The largest θ is at the far corners: check against a dense grid over the face.
    u, v = np.meshgrid(np.linspace(-25, 25, 401), np.linspace(-25, 25, 401))
    pts = d.centre + u.reshape(-1, 1) * d.u + v.reshape(-1, 1) * d.v
    grid_max = np.degrees(np.arccos(pts[:, 2] / np.linalg.norm(pts, axis=1))).max()
    assert hi == pytest.approx(grid_max, abs=1e-6)
    assert d.mean_theta() == pytest.approx(35.0, abs=0.5)
    cd = array["CD"]  # annulus on the axis, upstream: θ from 180 − atan(41/40) to 180 − atan(9/40)
    assert cd.theta_range() == pytest.approx((180 - math.degrees(math.atan(41 / 40)),
                                              180 - math.degrees(math.atan(9 / 40))), abs=1e-6)
    assert cd.phi_range() == (-180.0, 180.0)
    lo, hi = array["DSSD2"].phi_range()  # at φ = 180°: the range is unwrapped around it, symmetric about 180°
    assert lo < 180 < hi and (lo + hi) / 2 == pytest.approx(180.0, abs=1e-9)
    assert array["DSSD1"].phi_range() == pytest.approx((-20.874, 20.874), abs=1e-3)


# -- layout warnings --------------------------------------------------------------------------------------------


def test_shadowing_and_beam_warnings():
    front = Geometry("front", "circle", [0, 30, 50], -_unit([0, 30, 50]), radius=5)
    back = Geometry("back", "circle", [0, 60, 100], -_unit([0, 60, 100]), radius=20)
    array = Array([front, back])
    shade = array.shadowing()
    # The small disc sits entirely inside the big one's cone: it hides exactly its own solid angle.
    assert shade[("front", "back")] == pytest.approx(front.solid_angle() / back.solid_angle(), rel=1e-3)
    assert ("back", "front") not in shade
    assert any("front hides" in w for w in array.warnings())
    exp = Experiment.example("oxygen_on_lead_array")
    assert Array.from_experiment(exp).warnings() == ["CD hides <1% of DSSD4 from the target."]
    in_beam = Array([Geometry("zero", "circle", [0, 0, 100], [0, 0, -1], radius=5),
                     Geometry("back", "circle", [0, 0, -50], [0, 0, 1], radius=5),
                     Geometry("ring", "annular", [0, 0, -50], [0, 0, 1], inner_radius=4, outer_radius=20)])
    w = in_beam.warnings()
    assert any("zero sits in the beam" in x for x in w)
    assert any("back blocks the incoming beam" in x for x in w)
    assert not any(x.startswith("ring") and "beam" in x for x in w)


def test_visible_part_of_a_hidden_detector():
    """A small disc centred in front of a big one hides a cone: what the big one still sees is the difference of
    two cones' solid angles, 2π(cos α_front − cos α_back)."""
    from physim.nuclear.rates import SHADOW_ORDER

    axis = _unit([0, 1, 1])
    front = Geometry("front", "circle", 50 * axis, -axis, radius=2)
    back = Geometry("back", "circle", 100 * axis, -axis, radius=10)
    array = Array([front, back])
    exact = 2 * math.pi * (math.cos(math.atan(2 / 50)) - math.cos(math.atan(10 / 100)))
    dirs, dom = back.directions(None, SHADOW_ORDER)
    assert float(np.sum(dom * array.visible("back", dirs))) == pytest.approx(exact, rel=5e-3)
    # The front disc is not hidden at all, and directions that miss a detector do not reach it.
    dirs, dom = front.directions(None, 12)
    assert array.visible("front", dirs).all()
    assert np.isinf(front.distances(np.array([[1.0, 0.0, 0.0]]))).all()
    # Distances agree with hit().
    d = _unit([0, 1.02, 1])
    assert back.distances(d[None, :])[0] == pytest.approx(back.hit(d).distance)


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def test_exit_path_through_the_target():
    t = 1.0
    assert exit_path(direction(0), t, 0.5) == (0.5, False)
    assert exit_path(direction(60), t, 0.5)[0] == pytest.approx(1.0)
    assert exit_path(direction(180), t, 0.25)[0] == pytest.approx(0.25)  # backwards: out through the front
    assert exit_path(direction(120), t, 0.25)[0] == pytest.approx(0.5)
    # Near 90° the path through a flat target grows as 1/cos; exactly at 90° it is capped.
    assert exit_path(direction(89.9), t, 0.5)[0] == pytest.approx(0.5 / math.cos(math.radians(89.9)))
    assert exit_path(direction(90), t, 0.5) == (1000.0, True)
    # A target tilted by 30° about y: a track along the new normal crosses it square-on.
    assert exit_path(direction(30, 0), t, 0.0, tilt_deg=30)[0] == pytest.approx(1.0)


# -- response ---------------------------------------------------------------------------------------------------


def test_response_dead_layer_punch_through_threshold_resolution():
    exp = Experiment.example("oxygen_on_lead_array")
    det = exp.detectors[0]  # 500 um Si, 0.5 um dead layer, 40 keV FWHM, 300 keV threshold
    r = Response("16O", det)
    path = r.dead_layer + r.thickness
    assert r.stopping.range(r.punch_through_energy()) == pytest.approx(path, rel=1e-6)
    assert r.punch_through_energy(60.0) > r.punch_through_energy()
    e = 50.0  # stops inside: deposits everything that passed the dead layer
    assert r.deposited(e) == pytest.approx(r.after_dead_layer(e), rel=1e-12)
    assert r.after_dead_layer(e) == pytest.approx(r.stopping.energy_after(e, r.dead_layer), rel=1e-12)
    assert r.after_dead_layer(e, 45.0) < r.after_dead_layer(e)
    hi = r.punch_through_energy() * 1.5  # punches through: deposits only part
    assert r.deposited(hi) < r.after_dead_layer(hi)
    assert r.deposited(hi) == pytest.approx(r.after_dead_layer(hi) - r.stopping.energy_after(
        r.after_dead_layer(hi), r.thickness), rel=1e-12)
    assert math.isnan(r.measured_energy(0.2))  # below the 300 keV threshold
    rng = np.random.default_rng(7)
    samples = r.measured_energy(np.full(40000, e), rng=rng)
    assert np.std(samples) == pytest.approx(0.040 / 2.3548, rel=0.02)
    assert np.mean(samples) == pytest.approx(r.deposited(e), abs=3e-4)
    array = Array.from_experiment(exp)
    assert array.response("CD", "4He").thickness == pytest.approx(r.material.areal_density_mg_cm2("300 um"))


# -- pictures ---------------------------------------------------------------------------------------------------


def test_setup_pictures():
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    from physim.nuclear import plot

    exp = Experiment.example("oxygen_on_lead_array")
    ax = plot.coverage(exp)
    lines = {ln.get_label(): ln for ln in ax.get_lines() if not ln.get_label().startswith("_")}
    assert set(lines) == {"DSSD1", "DSSD2", "DSSD3", "DSSD4", "CD"}
    th = lines["DSSD1"].get_ydata()
    lo, hi = Array.from_experiment(exp)["DSSD1"].theta_range()
    assert np.nanmin(th) == pytest.approx(lo, abs=0.05) and np.nanmax(th) == pytest.approx(hi, abs=0.05)
    ax3 = plot.setup_3d(exp)
    assert len(ax3.texts) == 5
