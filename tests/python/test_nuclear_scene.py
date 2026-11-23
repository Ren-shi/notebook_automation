"""The interactive scene (backlog 52): solids to scale, where a detector may be moved, and the quick numbers that
follow a drag."""

import asyncio
import math
import re
from types import SimpleNamespace

import numpy as np
import pytest

from physim.nuclear import Experiment, scene
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates

BASE = {
    "schema": "physim.experiment/1",
    "title": "scene",
    "beam": {"nuclide": "16O", "energy": "60 MeV", "current": "1 pnA"},
    "target": {"material": "194Pt", "thickness": "1 mg/cm2"},
    "run": {"beam_time": "1 h"},
}


def setup(detectors, gammas=(), **extra) -> dict:
    d = dict(BASE, detectors=[dict(x) for x in detectors], **extra)
    if gammas:
        d["gamma_detectors"] = [dict(g) for g in gammas]
    return d


S3 = {"name": "S3", "model": "S3", "theta": "180 deg", "distance": "20 mm"}
PAD = {"name": "Pad", "shape": "rectangle", "width": "20 mm", "height": "20 mm", "thickness": "300 um",
       "theta": "60 deg", "phi": "0 deg", "distance": "80 mm"}
CLOVER = {"name": "Clo", "model": "clover", "theta": "90 deg", "phi": "90 deg", "distance": "200 mm"}


def by_key(exp) -> dict:
    return {s.key: s for s in scene.solids(exp)}


# -- solids, to scale -----------------------------------------------------------------------------------------------


def test_the_s3_is_drawn_with_its_real_dimensions():
    exp = Experiment.from_dict(setup([S3]))
    s = by_key(exp)["detector:0"]
    assert np.allclose(s.centre, [0, 0, -20])
    rings = [p for p in s.parts if p.role == "active" and p.element is not None]
    assert len(rings) == 24
    assert rings[0].r_in == pytest.approx(11.0) and rings[-1].r_out == pytest.approx(35.0)
    assert all(b.r_in == pytest.approx(a.r_out) for a, b in zip(rings, rings[1:]))
    assert sum(p.kind == "line" for p in s.parts) == 32  # the sector boundaries
    board = next(p for p in s.parts if p.role == "board")
    assert board.r_out == pytest.approx(50.0)
    # The face looks at the target: its normal points from the detector to the origin.
    assert np.allclose(s.axes[:, 2], [0, 0, 1])


def test_a_clover_has_four_crystals_in_its_housing():
    exp = Experiment.from_dict(setup([S3], [CLOVER]))
    s = by_key(exp)["gamma:0"]
    crystals = [p for p in s.parts if p.role == "crystal"]
    assert len(crystals) == 4
    assert {(round(abs(p.cx), 3), round(abs(p.cy), 3)) for p in crystals} == {(22.5, 22.5)}
    assert all(p.r_out == pytest.approx(25.0) and p.depth == pytest.approx(70.0) for p in crystals)
    # In the laboratory the crystal faces are where the setup puts them.
    faces = sorted(tuple(np.round(s.world([[p.cx, p.cy, 0.0]])[0], 6)) for p in crystals)
    assert faces == sorted(tuple(np.round(c, 6)) for _, c, _ in exp.gamma_detectors[0].elements())
    assert next(p for p in s.parts if p.role == "housing").w == pytest.approx(101.0)


def test_every_frame_is_a_rotation():
    for name in Planner.examples():
        for s in scene.solids(Experiment.example(name)):
            assert np.allclose(s.axes.T @ s.axes, np.eye(3), atol=1e-12), s.key
            assert np.linalg.det(s.axes) == pytest.approx(1.0), s.key


def test_the_examples_have_no_overlaps():
    for name in Planner.examples():
        assert scene.problems(Experiment.example(name)) == [], name


# -- what a layout may not do ---------------------------------------------------------------------------------------


def test_a_detector_in_the_beam_is_found():
    exp = Experiment.from_dict(setup([dict(PAD, theta="0 deg")]))
    found = scene.problems(exp)
    assert [p.kind for p in found] == ["beam"] and "Pad" in found[0].text
    # An annular detector lets the beam through its hole.
    assert scene.problems(Experiment.from_dict(setup([S3]))) == []


def test_overlaps_are_found_even_between_thin_plates():
    crossing = [PAD, dict(PAD, name="Pad2", theta="62 deg", distance="80 mm")]
    found = scene.problems(Experiment.from_dict(setup(crossing)))
    assert [p.kind for p in found] == ["overlap"]
    apart = [PAD, dict(PAD, name="Pad2", theta="100 deg")]
    assert scene.problems(Experiment.from_dict(setup(apart))) == []
    # A thin detector standing inside a clover's housing.
    inside = scene.problems(Experiment.from_dict(setup([dict(PAD, theta="90 deg", phi="90 deg", distance="220 mm")],
                                                       [CLOVER])))
    assert [p.kind for p in inside] == ["overlap"]


def annulus(name, r_in, r_out, distance=50.0, **extra) -> dict:
    return dict({"name": name, "shape": "annular", "theta": "0 deg", "distance": f"{distance} mm",
                 "inner_radius": f"{r_in} mm", "outer_radius": f"{r_out} mm", "rings": 4, "thickness": "1000 um"},
                **extra)


def kinds(detectors, gammas=()) -> list:
    return [p.kind for p in scene.problems(Experiment.from_dict(setup(detectors, gammas)))]


def stands_there(detectors) -> list:
    return [w.text for w in Planner(Experiment.from_dict(setup(detectors))).warnings() if "both stand" in w.text]


def roof(theta, gap=0.0, tilt=45.0) -> list:
    """Two pads that share a front edge at 80 mm in the direction theta, each tilted by ``tilt`` away from the
    plane square to that direction, and pulled apart by ``gap``."""
    t, k = math.radians(theta), math.radians(tilt)
    u, e = np.array([math.sin(t), 0, math.cos(t)]), np.array([math.cos(t), 0, -math.sin(t)])
    out = []
    for name, s in (("A", 1), ("B", -1)):
        c = 80 * u + s * (gap / 2) * e + 10 * (s * math.cos(k) * e - math.sin(k) * u)
        n = -math.cos(k) * u - s * math.sin(k) * e  # the normal: towards the target, leaning to the shared plane
        d = {"name": name, "shape": "rectangle", "width": "20 mm", "height": "20 mm", "thickness": "300 um"}
        out.append(dict(d, position=[f"{x:.9f} mm" for x in c], facing=[float(x) for x in n]))
    return out


def test_rings_that_share_an_edge_only_touch():
    # Allmond 2014: two rings in one plane that meet at 28 degrees.
    r = 26.5855
    touching = [annulus("In", 12.4664, r), annulus("Out", r, 48.2844)]
    assert kinds(touching) == [] and stands_there(touching) == []
    for gap in (0.01, 0.5):
        assert kinds([annulus("In", 12, r), annulus("Out", r + gap, 48)]) == []
    for depth in (0.1, 0.5, 5):
        over = [annulus("In", 12.4664, r), annulus("Out", r - depth, 48.2844)]
        assert kinds(over) == ["overlap"], depth
        assert stands_there(over) == ["In overlaps Out. They cannot both stand there."]
    # Within the tolerance on one side, beyond it on the other.
    assert kinds([annulus("In", 12, 30), annulus("Out", 30 - 0.9 * scene.TOUCH_MM, 48)]) == []
    assert kinds([annulus("In", 12, 30), annulus("Out", 30 - 1.1 * scene.TOUCH_MM, 48)]) == ["overlap"]


def test_nested_and_identical_detectors_overlap():
    ring, pad = annulus("A", 12, 40), PAD
    assert kinds([ring, dict(ring, name="B")]) == ["overlap"]
    assert kinds([pad, dict(pad, name="B")]) == ["overlap"]
    assert kinds([S3], [CLOVER, dict(CLOVER, name="B")]) == ["overlap"]
    assert kinds([ring, annulus("B", 20, 30)]) == ["overlap"]  # one within the other
    assert kinds([ring, annulus("B", 12, 30)]) == ["overlap"]  # sharing the inner edge
    assert kinds([ring, annulus("B", 20, 40)]) == ["overlap"]  # sharing the outer edge
    # Three rings, each in the hole of the next.
    assert kinds([annulus("A", 10, 20), annulus("B", 20, 30), annulus("C", 30, 40)]) == []


def test_stacked_detectors_touch_face_to_face():
    # The second starts where the first (1 mm thick) ends.
    assert kinds([annulus("dE", 12, 40), annulus("E", 12, 40, 51.0)]) == []
    assert kinds([annulus("dE", 12, 40), annulus("E", 12, 40, 50.5)]) == ["overlap"]
    thin = {"thickness": "20 um"}
    assert kinds([annulus("dE", 12, 40, **thin), annulus("E", 12, 40, 50.02)]) == []
    assert kinds([annulus("dE", 12, 40, **thin), annulus("E", 12, 40, 50.01)]) == ["overlap"]


def test_rectangles_that_share_an_edge_only_touch():
    # Side by side in one plane at 80 mm: their centres are 20 mm apart.
    flat = {k: v for k, v in PAD.items() if k not in ("theta", "phi", "distance")}
    side = lambda x: dict(flat, name=f"P{x}", position=[f"{x} mm", "30 mm", "80 mm"], facing=[0, 0, -1])  # noqa: E731
    assert kinds([side(-10), side(10)]) == []
    assert kinds([side(-10), side(9.5)]) == ["overlap"]
    assert kinds([side(-10), side(9.9)]) == ["overlap"]
    # A corner in the other's corner.
    assert kinds([side(-10), dict(side(10), position=["10 mm", "50 mm", "80 mm"])]) == []
    assert kinds([side(-10), dict(side(10), position=["9.5 mm", "49.5 mm", "80 mm"])]) == ["overlap"]
    for theta, tilt in ((50, 45), (73, 20), (90, 60)):
        assert kinds(roof(theta, tilt=tilt)) == [], (theta, tilt)
        assert kinds(roof(theta, gap=-0.5, tilt=tilt)) == ["overlap"], (theta, tilt)


def test_thin_bars_that_cross_overlap():
    bar = {"shape": "rectangle", "thickness": "300 um", "facing": [0, 0, -1]}
    a = dict(bar, name="A", width="100 mm", height="1 mm", position=["0 mm", "30 mm", "80 mm"])
    b = dict(bar, name="B", width="1 mm", height="100 mm", position=["17 mm", "60 mm", "80 mm"])
    assert kinds([a, b]) == ["overlap"]
    assert kinds([a, dict(b, position=["17 mm", "80.5 mm", "80 mm"])]) == []
    assert kinds([a, dict(b, position=["17 mm", "80 mm", "80 mm"])]) == ["overlap"]


def test_a_detector_can_be_dragged_until_it_touches_another():
    exp = Experiment.from_dict(setup([annulus("dE", 12, 40), annulus("E", 12, 40, 60.0)]))
    lim = scene.limits(exp, "detector:1")
    assert "dE would overlap E" in lim.lo_why and lim.lo == pytest.approx(51.0, abs=0.1)
    p = Planner(exp)
    pos = lim.apply((0.0, 0.0, 20.0))  # dragged right through the first
    fields = scene.move_fields(exp, "detector:1", pos, lim.mode)
    assert scene.blocked(exp, "detector:1", fields) is None and p.move("detector:1", pos)
    assert scene.problems(p.experiment) == []
    assert lim.stopped(pos) == lim.lo_why
    assert "overlap" in scene.blocked(exp, "detector:1", {"distance": "50.5 mm"})
    # Sideways, against a neighbour.
    exp = Experiment.from_dict(setup([PAD, dict(PAD, name="Pad2", theta="100 deg")]))
    touch = 60 + 2 * math.degrees(math.atan(10 / 80))  # where the front edges meet
    assert scene.blocked(exp, "detector:1", {"theta": f"{touch:.9f} deg"}) is None
    assert "Pad would overlap Pad2" in scene.blocked(exp, "detector:1", {"theta": f"{touch - 0.5:.2f} deg"})


def test_a_drop_in_the_view_stops_where_the_detectors_touch():
    """The same through the app's 3D view: its drag and drop handlers in NiceGUI's simulated browser."""
    pytest.importorskip("nicegui")
    from nicegui.testing import user_simulation

    from physim.nuclear.scene_view import SceneView

    p = Planner(Experiment.from_dict(setup([annulus("dE", 12, 40), annulus("E", 12, 40, 60.0)])))
    views, lives = [], []

    def page():
        views.append(SceneView(lambda: p, on_live=lives.append, on_moved=p.move))

    async def main():
        async with user_simulation(page) as user:
            await user.open("/")
            with user:
                view = views[-1]
                view.select("detector:1")
                at = dict(object_name="detector:1", x=0.0, y=0.0, z=20.0)
                view._dragging(SimpleNamespace(args=at))
                view._dropped(SimpleNamespace(**at))

    asyncio.run(main())
    assert lives and lives[-1]["problem"] is None
    assert p.draft["detectors"][1]["distance"] == "51 mm" and scene.problems(p.experiment) == []


def test_the_beam_is_kept_clear_to_its_radius():
    # The overlaps' tolerance is not given to the beam: a hole must hold the whole beam.
    hole = scene.BEAM_RADIUS_MM
    assert kinds([annulus("A", hole + 0.01, 30)]) == []
    assert kinds([annulus("A", hole - 0.01, 30)]) == ["beam"]


def test_a_move_through_the_chamber_wall_is_blocked():
    chamber = {"radius": "100 mm", "wall_thickness": "3 mm"}
    exp = Experiment.from_dict(setup([PAD], [CLOVER], chamber=chamber))
    assert scene.blocked(exp, "detector:0", {"distance": "60 mm"}) is None
    assert "outside the chamber" in scene.blocked(exp, "detector:0", {"distance": "120 mm"})
    assert "inside the chamber wall" in scene.blocked(exp, "gamma:0", {"distance": "90 mm"})


# -- moving ---------------------------------------------------------------------------------------------------------


def test_dragging_the_s3_is_the_same_as_typing_the_distance():
    dragged = Planner(Experiment.from_dict(setup([S3])))
    typed = Planner(Experiment.from_dict(setup([S3])))
    assert dragged.move("detector:0", (0.0, 0.0, -40.0))
    assert typed.set("detector 1", "distance", "40 mm")
    assert dragged.draft["detectors"][0]["distance"] == "40 mm"
    a, b = dragged.rates(), typed.rates()
    assert a["rows"][0]["theta_range"] == b["rows"][0]["theta_range"]
    assert a["rows"][0]["rate_per_s"] == b["rows"][0]["rate_per_s"]
    assert a["strips"] == b["strips"] or all(
        np.allclose(a["strips"][k]["rate_per_s"], b["strips"][k]["rate_per_s"]) for k in a["strips"])
    # The rings have moved to new angles: at 40 mm the outer ring is at 180° − atan(35/40).
    lo, _ = dragged.geometry()["detectors"][0]["theta_range"]
    assert lo == pytest.approx(180 - math.degrees(math.atan(35 / 40)), abs=1e-6)


def test_a_detector_on_the_axis_only_slides_along_it():
    exp = Experiment.from_dict(setup([S3]))
    assert scene.on_axis(exp, "detector:0")
    assert scene.move_fields(exp, "detector:0", (15.0, -8.0, -33.0)) == {"distance": "33 mm"}
    lim = scene.limits(exp, "detector:0")
    assert lim.mode == "distance" and lim.lo == pytest.approx(scene.MIN_DISTANCE_MM)
    assert np.allclose(lim.apply((15.0, -8.0, -33.0)), [0, 0, -33.0])
    assert np.allclose(lim.apply((0.0, 0.0, 5.0)), [0, 0, -lim.lo])  # it does not pass through the target


def test_moving_in_angle_keeps_the_distance_and_in_distance_the_direction():
    exp = Experiment.from_dict(setup([PAD], [CLOVER]))
    f = scene.move_fields(exp, "gamma:0", (100.0, 0.0, 100.0), "angle")
    assert f == {"theta": "45 deg", "phi": "0 deg"}
    f = scene.move_fields(exp, "detector:0", 2 * scene.position_of(exp, "detector:0"), "distance")
    assert f == {"distance": "160 mm"}
    by_position = dict(PAD, position=["40 mm", "0 mm", "40 mm"])
    del by_position["theta"], by_position["phi"], by_position["distance"]
    exp = Experiment.from_dict(setup([by_position]))
    f = scene.move_fields(exp, "detector:0", (0.0, 30.0, 30.0), "angle")
    assert np.allclose([float(c.split()[0]) for c in f["position"]], [0, 40, 40], atol=0.01)


def test_limits_stop_at_the_beam_at_another_solid_and_at_the_wall():
    exp = Experiment.from_dict(setup([S3, PAD], [CLOVER], chamber={"radius": "150 mm"}))
    lim = scene.limits(exp, "detector:1", "angle")
    assert lim.mode == "angle" and 5 < lim.lo < 20 and "beam" in lim.lo_why
    for theta in (lim.lo + 0.2, lim.hi - 0.2):
        t = math.radians(theta)
        f = scene.move_fields(exp, "detector:1", (80 * math.sin(t), 0, 80 * math.cos(t)), "angle")
        assert not [p for p in scene.problems(scene.moved(exp, "detector:1", f), only="detector:1")
                    if p.kind == "beam"]
    assert scene.blocked(exp, "detector:1", {"theta": f"{lim.lo - 1:.2f} deg"})
    out = scene.limits(exp, "detector:1", "distance")
    assert "chamber" in out.hi_why and out.hi < 150
    assert scene.blocked(exp, "detector:1", {"distance": f"{out.hi - 0.2:.2f} mm"}) is None
    assert out.stopped(out.apply(1e3 * np.array(out.direction))) == out.hi_why
    # The clover comes no closer than the chamber wall.
    g = scene.limits(exp, "gamma:0", "distance")
    assert "chamber" in g.lo_why and g.lo >= 150


def _js(expression: str) -> str:
    """A drag constraint (JavaScript) as a Python expression: ``c ? a : b`` becomes a conditional."""
    expression = expression.replace("Math.sqrt", "math.sqrt").replace("Math.abs", "abs")
    depth = 0
    for i, ch in enumerate(expression):
        depth += ch == "("
        depth -= ch == ")"
        if ch == "?" and depth == 0:
            level = 0
            for j in range(i + 1, len(expression)):
                level += expression[j] == "?"
                if expression[j] == ":" and level == 0:
                    return (f"(({_js(expression[i + 1:j])}) if ({expression[:i]}) "
                            f"else ({_js(expression[j + 1:])}))")
                level -= expression[j] == ":"
    out, i = "", 0
    while i < len(expression):  # the same inside brackets
        if expression[i] == "(":
            depth, j = 1, i + 1
            while depth:
                depth += expression[j] == "("
                depth -= expression[j] == ")"
                j += 1
            out += "(" + _js(expression[i + 1:j - 1]) + ")"
            i = j
        else:
            out += expression[i]
            i += 1
    return out


def run_constraints(text: str, position) -> np.ndarray:
    """Apply the scene's drag constraints as the browser does: one assignment after the other."""
    p = dict(zip("xyz", map(float, position)))
    for part in text.split(","):
        name, expression = (s.strip() for s in part.split("=", 1))
        p[name] = eval(_js(expression), {"math": math}, dict(p))  # noqa: S307 -- our own expression
    return np.array([p["x"], p["y"], p["z"]])


@pytest.mark.parametrize("key, mode", [("detector:0", "angle"), ("detector:1", "angle"), ("detector:1", "distance"),
                                       ("gamma:0", "angle"), ("gamma:0", "distance")])
def test_the_browser_constraints_follow_the_same_rule(key, mode):
    exp = Experiment.from_dict(setup([S3, PAD], [CLOVER]))
    lim = scene.limits(exp, key, mode)
    text = lim.constraints()
    # NiceGUI splits the constraints at commas and substitutes every letter x, y or z.
    assert len(text.split(",")) == 4
    bare = text.replace("Math.sqrt", "").replace("Math.abs", "")
    assert not re.search(r"[xyz]", re.sub(r"\b[xyz]\b", "", bare))
    rng = np.random.default_rng(3)
    for p in rng.normal(scale=150.0, size=(40, 3)):
        assert np.allclose(run_constraints(text, p), lim.apply(p), atol=1e-3), (key, mode, p)


def test_no_drag_ends_in_the_beam_or_inside_another_solid():
    """Whatever the mouse does, the position kept is an allowed one (the rule the view applies on release)."""
    exp = Experiment.example("coulex_ni58")
    rng = np.random.default_rng(7)
    for s in scene.solids(exp)[1:]:
        for mode in ("angle", "distance"):
            lim = scene.limits(exp, s.key, mode)
            good = scene.position_of(exp, s.key)
            for p in rng.normal(scale=200.0, size=(12, 3)):
                pos = lim.apply(p)
                fields = scene.move_fields(exp, s.key, pos, lim.mode)
                if scene.blocked(exp, s.key, fields):
                    pos = good  # the view goes back to the last allowed place
                    fields = scene.move_fields(exp, s.key, pos, lim.mode)
                else:
                    good = pos
                assert scene.problems(scene.moved(exp, s.key, fields), only=s.key) == [], (s.key, mode)


# -- the numbers that follow a drag ---------------------------------------------------------------------------------


def test_the_quick_rate_matches_the_full_one():
    p = Planner.example("coulex_ni58")
    full = p._rates()
    for i, det in enumerate(p.experiment.detectors):
        live = p.live(f"detector:{i}")
        assert live["rate_per_s"] == pytest.approx(full.rate(det.name, counted=False), rel=1e-3)
        assert live["solid_angle_msr"] == pytest.approx(p.geometry()["detectors"][i]["solid_angle_msr"])
        assert live["problem"] is None


def test_live_values_follow_a_trial_position_without_changing_the_setup():
    p = Planner(Experiment.from_dict(setup([S3, PAD])))
    before = p.to_toml()
    near, far = p.live("detector:0"), p.live("detector:0", {"distance": "40 mm"})
    assert far["distance_mm"] == pytest.approx(40.0) and far["solid_angle_msr"] < near["solid_angle_msr"]
    assert far["theta_range"][0] > near["theta_range"][0]
    assert p.live("detector:1", {"theta": "0 deg"})["problem"]
    assert p.to_toml() == before


def test_rates_keep_the_rows_of_detectors_that_did_not_change():
    p = Planner.example("oxygen_on_lead_array")
    p.rates()
    assert p.move("detector:0", scene.position_of(p.experiment, "detector:0") * 1.2, "distance")
    r = p._rates()
    # The moved detector is computed again; the others keep their rows, unless something shadows them.
    assert "DSSD1" not in r.reused and {"DSSD2", "DSSD3"} <= r.reused
    fresh = Rates(p.experiment)
    key = lambda row: (row.detector, row.segment, row.channel, row.particle)  # noqa: E731
    assert [key(x) for x in sorted(r.rows, key=key)] == [key(x) for x in sorted(fresh.rows, key=key)]
    assert np.allclose([x.rate for x in sorted(r.rows, key=key)], [x.rate for x in sorted(fresh.rows, key=key)],
                       rtol=1e-12)
    # A change of the beam recomputes everything.
    p.set("beam", "energy", "70 MeV")
    assert p._rates().reused == set()


def test_rings_closer_than_the_safe_distance_are_marked():
    d = Experiment.example("coulex_ni58").to_dict()
    assert Planner(Experiment.from_dict(d)).safety() == {}
    d["beam"]["energy"] = "48 MeV"
    p = Planner(Experiment.from_dict(d))
    unsafe = p.safety()
    cd = p.experiment.detectors[0].name
    assert unsafe.get(cd), "the backward detector sees the closest collisions"
    s = p.selection("detector:0", max(unsafe[cd]))
    assert s["unsafe"] == sorted(unsafe[cd]) and s["element"]["safe"] is False


def test_advice_follows_the_kinematics():
    light = scene.advice(Experiment.from_dict(setup([S3])))
    assert "Backward angles are recommended" in light[0] and len(light) == 1
    d = setup([dict(PAD, theta="120 deg")])
    d["beam"] = {"nuclide": "208Pb", "energy": "900 MeV", "current": "1 pnA"}
    d["target"] = {"material": "12C", "thickness": "1 mg/cm2"}
    heavy = scene.advice(Experiment.from_dict(d))
    assert "Forward angles are recommended" in heavy[0]
    assert "Pad" in heavy[1] and "count nothing" in heavy[1]


def test_what_the_side_panel_shows():
    p = Planner.example("coulex_ni58")
    whole = p.selection()
    assert whole["kind"] == "experiment" and whole["detectors"] == 3 and whole["gamma_detectors"] == 4
    assert whole["rate_per_s"] == pytest.approx(sum(r["rate_per_s"] for r in p.rates()["rows"]))
    det = p.selection("detector:0", 2)
    assert det["kind"] == "detector" and det["on_axis"] and det["element"]["label"] == "ring 3"
    rings = [p.selection("detector:0", k)["element"] for k in range(p.experiment.detectors[0].rings)]
    assert sum(r["rate_per_s"] for r in rings) == pytest.approx(det["rate_per_s"])
    assert sum(r["solid_angle_msr"] for r in rings) == pytest.approx(det["solid_angle_msr"])
    gam = p.selection("gamma:0")
    assert gam["kind"] == "gamma" and gam["coincidence_rate_per_s"] > 0
    clover = Planner(Experiment.from_dict(setup([S3], [CLOVER]))).selection("gamma:0", 1)
    assert clover["crystals"] == ["A", "B", "C", "D"] and clover["element"]["label"] == "crystal B"


def test_scene_names_and_labels():
    from physim.nuclear import scene_view

    assert scene_view.parse_name("detector:2|5") == ("detector:2", 5)
    assert scene_view.parse_name("gamma:0") == ("gamma:0", None)
    assert scene_view.parse_name(None) == (None, None) and scene_view.parse_name("ground") == (None, None)
    solids = by_key(Experiment.from_dict(setup([S3, PAD], [CLOVER])))
    assert scene_view.dimensions(solids["detector:0"]) == "⌀ 22–70 mm, 1000 µm thick"
    assert scene_view.dimensions(solids["detector:1"]) == "20 × 20 mm, 300 µm thick"
    assert scene_view.dimensions(solids["gamma:0"]) == "4 × ⌀ 50 × 70 mm"
    assert set(scene_view.PALETTES["light"]) == set(scene_view.PALETTES["dark"])
