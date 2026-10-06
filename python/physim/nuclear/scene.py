"""The experiment as solids, to scale: what the interactive scene draws, and where a detector may be moved.

Everything here is plain geometry, with no drawing code, so it can be tested and used without the app::

    from physim.nuclear import Experiment, scene

    exp = Experiment.example("coulex_ni58")
    for s in scene.solids(exp):
        print(s.key, s.name, s.centre)
    scene.problems(exp)                              # overlaps and detectors in the beam
    scene.move_fields(exp, "detector:0", (0, 0, -40))  # the setup fields of a detector dragged there

A :class:`Solid` has a centre, three axes and :class:`Part` s given in its own coordinates: the third axis ``n``
points from the front face towards the target, and the material lies behind the face (at negative local z).
Lengths are mm. The beam runs along +z through the target at the origin.

Sizes that the setup file does not give are drawn with the typical values below (the target's diameter, the
thickness of a circuit board, the length of a γ-ray detector without crystal dimensions).
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .detectors import Geometry
from .quantity import Quantity

#: Diameter of the target foil as drawn, mm (a setup file gives only its thickness).
TARGET_DIAMETER_MM = 10.0
#: Thickness of a detector's circuit board as drawn, mm (typical FR4).
BOARD_THICKNESS_MM = 1.6
#: Smallest radius of the beam as drawn and as kept clear, mm.
BEAM_RADIUS_MM = 1.0
#: Closest a dragged detector may come to the target, mm.
MIN_DISTANCE_MM = 2.0


def _q(x) -> Optional[Quantity]:
    if x is None:
        return None
    return x if isinstance(x, Quantity) else Quantity.parse(x)


def _mm(x, default: float = 0.0) -> float:
    try:
        return float(_q(x).to("mm")) if x is not None else default
    except (ValueError, KeyError, AttributeError):
        return default


@dataclass
class Part:
    """One piece of a solid, in the solid's coordinates. ``kind`` is "ring" (an annulus or, with ``r_in`` = 0, a
    cylinder, centred at ``cx``, ``cy``), "box" (``w`` × ``h``) or "line" (``points``, drawn only). The piece lies
    between local z = ``z0`` (its front) and ``z0 − depth``.

    ``role`` says what it is ("active", "board", "crystal", "housing", "foil", "divider"); ``element`` is the ring,
    strip or crystal number for a piece that can be selected on its own; ``hull`` marks the pieces used when
    checking for overlaps (the whole active volume, not each ring of it)."""

    kind: str
    role: str
    z0: float = 0.0
    depth: float = 0.0
    cx: float = 0.0
    cy: float = 0.0
    r_in: float = 0.0
    r_out: float = 0.0
    w: float = 0.0
    h: float = 0.0
    points: Optional[list] = None
    element: Optional[int] = None
    label: str = ""
    hull: bool = False

    def lateral(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Whether local points (x, y) lie within the piece, seen along its axis."""
        if self.kind == "box":
            return (np.abs(x - self.cx) <= self.w / 2) & (np.abs(y - self.cy) <= self.h / 2)
        r = np.hypot(x - self.cx, y - self.cy)
        return (r >= self.r_in) & (r <= self.r_out)

    def reach(self) -> float:
        """Largest distance of the piece from the solid's axis."""
        far = math.hypot(self.w, self.h) / 2 if self.kind == "box" else self.r_out
        return math.hypot(self.cx, self.cy) + far

    def edges(self) -> np.ndarray:
        """Line segments along the piece's edges, as an array (n, 2, 3) in local coordinates."""
        if self.kind == "box":
            a, b = self.w / 2, self.h / 2
            c = np.array([[-a, -b], [a, -b], [a, b], [-a, b]]) + [self.cx, self.cy]
            loops = [np.c_[c, np.full(4, z)] for z in (self.z0, self.z0 - self.depth)]
        else:
            t = np.linspace(0.0, 2 * math.pi, 49)[:-1]
            loops = []
            for r in ([self.r_out, self.r_in] if self.r_in > 0 else [self.r_out]):
                for z in (self.z0, self.z0 - self.depth):
                    loops.append(np.c_[self.cx + r * np.cos(t), self.cy + r * np.sin(t), np.full(len(t), z)])
        segs = [np.stack([lp, np.roll(lp, -1, axis=0)], axis=1) for lp in loops]
        # The edges along the axis join the front loop to the back one (every sixth point of a circle).
        step = 1 if self.kind == "box" else 6
        for front, back in zip(loops[0::2], loops[1::2]):
            segs.append(np.stack([front[::step], back[::step]], axis=1))
        return np.concatenate(segs)


@dataclass
class Solid:
    """A detector, or the target, as it stands in the chamber. ``key`` is "detector:N", "gamma:N" (N from 0, as in
    the setup file) or "target"; the columns of ``axes`` are the solid's u, v and n in the laboratory."""

    key: str
    name: str
    kind: str
    centre: np.ndarray
    axes: np.ndarray
    parts: list = field(default_factory=list)
    model: Optional[str] = None

    @property
    def index(self) -> int:
        return int(self.key.split(":")[1]) if ":" in self.key else 0

    def hull(self) -> list:
        return [p for p in self.parts if p.hull]

    def local(self, points: np.ndarray) -> np.ndarray:
        return (np.asarray(points, dtype=float) - self.centre) @ self.axes

    def world(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(points, dtype=float) @ self.axes.T + self.centre

    def reach(self) -> float:
        """Radius of a sphere about the centre that holds the whole solid."""
        return max((math.hypot(p.reach(), max(abs(p.z0), abs(p.z0 - p.depth))) for p in self.hull()), default=0.0)

    def crossed_by(self, segments: np.ndarray) -> bool:
        """Whether any of the line segments (n, 2, 3), in laboratory coordinates, passes through the solid."""
        if not len(segments):
            return False
        p0, p1 = self.local(segments[:, 0]), self.local(segments[:, 1])
        dz = p1[:, 2] - p0[:, 2]
        safe = np.where(np.abs(dz) < 1e-12, 1e-12, dz)
        for part in self.hull():
            lo, hi = part.z0 - part.depth, part.z0
            # The stretch of each segment between the piece's front and back planes ...
            ta, tb = (lo - p0[:, 2]) / safe, (hi - p0[:, 2]) / safe
            t0, t1 = np.clip(np.minimum(ta, tb), 0, 1), np.clip(np.maximum(ta, tb), 0, 1)
            flat = np.abs(dz) < 1e-12
            inside = (p0[:, 2] >= lo) & (p0[:, 2] <= hi)
            t0, t1 = np.where(flat, 0.0, t0), np.where(flat, np.where(inside, 1.0, -1.0), t1)
            ok = t1 > t0 + 1e-12
            if not ok.any():
                continue
            # ... and whether that stretch lies within the piece's outline, at a few points along it.
            for f in np.linspace(0.0, 1.0, 9):
                t = (t0 + f * (t1 - t0))[:, None]
                p = p0 + t * (p1 - p0)
                if (ok & part.lateral(p[:, 0], p[:, 1])).any():
                    return True
        return False

    def edges(self) -> np.ndarray:
        """The edges of the solid's hull, as line segments (n, 2, 3) in laboratory coordinates."""
        segs = np.concatenate([p.edges() for p in self.hull()])
        return self.world(segs.reshape(-1, 3)).reshape(-1, 2, 3)


def _frame(u, v, n) -> np.ndarray:
    return np.column_stack([np.asarray(u, dtype=float), np.asarray(v, dtype=float), np.asarray(n, dtype=float)])


def _thickness_mm(det) -> float:
    """A detector's thickness as a length, mm (a thickness given in mg/cm² is divided by the density)."""
    q = _q(det.thickness)
    try:
        return float(q.to("mm"))
    except (ValueError, KeyError):
        pass
    try:
        from . import data

        rho = data.material(det.material or "Si", _q(det.density)).density_g_cm3
        return float(q.to("mg/cm2") / float(rho) * 1e-2)
    except Exception:  # noqa: BLE001 -- drawing must not fail on an unusual unit
        return 0.3


def _detector_solid(det, i: int) -> Solid:
    g = Geometry.from_setup(det, i + 1)
    t = max(_thickness_mm(det), 0.02)
    parts = []
    if g.shape == "rectangle":
        parts.append(Part("box", "active", depth=t, w=g.width, h=g.height, hull=True))
        dx = g.width / g.strips_x
        if g.strips_x > 1:
            for k in range(g.strips_x):
                parts.append(Part("box", "active", depth=t, cx=-g.width / 2 + (k + 0.5) * dx, w=dx, h=g.height,
                                  element=k, label=f"strip {k + 1}"))
        dy = g.height / g.strips_y
        for k in range(1, g.strips_y):
            y = -g.height / 2 + k * dy
            parts.append(Part("line", "divider", points=[[-g.width / 2, y, 0.02], [g.width / 2, y, 0.02]]))
    else:
        r_in, r_out = g.inner_radius or 0.0, g.outer_radius
        parts.append(Part("ring", "active", depth=t, r_in=r_in, r_out=r_out, hull=True))
        dr = (r_out - r_in) / g.rings
        if g.rings > 1:
            for k in range(g.rings):
                parts.append(Part("ring", "active", depth=t, r_in=r_in + k * dr, r_out=r_in + (k + 1) * dr,
                                  element=k, label=f"ring {k + 1}"))
        if g.sectors > 1:
            for k in range(g.sectors):
                a = 2 * math.pi * k / g.sectors
                c, s = math.cos(a), math.sin(a)
                parts.append(Part("line", "divider", points=[[r_in * c, r_in * s, 0.02], [r_out * c, r_out * s, 0.02]]))
    for name, b_in, b_out in det.blocking():
        if name == "board":
            parts.append(Part("ring", "board", z0=-0.3, depth=BOARD_THICKNESS_MM, r_in=b_in, r_out=b_out, hull=True,
                              label="circuit board"))
    return Solid(f"detector:{i}", g.name, "detector", g.centre + _origin(det), _frame(g.u, g.v, g.n), parts,
                 det.model)


def _origin(det) -> np.ndarray:
    """The target's place from the chamber's centre: the scene draws in the chamber's coordinates."""
    return np.array(det._origin, dtype=float) if getattr(det, "_origin", None) else np.zeros(3)


def _gamma_solid(gd, i: int) -> Solid:
    u = np.array(gd.direction())
    a, b = (np.array(x) for x in gd.face_axes())
    dist = gd.distance_mm()
    centre = dist * u + _origin(gd)
    r = gd.crystal_radius_mm()
    # Without crystal dimensions, a γ-ray detector is drawn as long as it is wide.
    length = _mm(gd.crystal_length, 2 * r)
    parts = []
    crystals = gd.elements()
    for k, (label, c, rad) in enumerate(crystals):
        off = np.array(c) - centre
        parts.append(Part("ring", "crystal", depth=length, cx=float(off @ b), cy=float(off @ a), r_out=rad,
                          element=k if len(crystals) > 1 else None, label=f"crystal {label}" if label else "crystal",
                          hull=gd.housing_side is None))
    if gd.housing_side is not None:
        gap = _mm(gd.window_gap)
        side = _mm(gd.housing_side)
        parts.append(Part("box", "housing", z0=gap, depth=gap + length + gap, w=side, h=side, hull=True,
                          label="housing"))
    # (b, a, −u) is right-handed, as a rotation needs.
    return Solid(f"gamma:{i}", gd.name or f"γ{i + 1}", "gamma", centre, _frame(b, a, -u), parts, gd.model)


def target_thickness_mm(experiment) -> float:
    """The target's thickness as a length, mm."""
    t = experiment.target
    q = _q(t.thickness)
    try:
        return float(q.to("mm"))
    except (ValueError, KeyError):
        pass
    try:
        return float(q.to("mg/cm2") / float(t.material_data().density_g_cm3) * 1e-2)
    except Exception:  # noqa: BLE001
        return 0.001


def _target_solid(experiment) -> Solid:
    from .rates import tilt_deg

    a = math.radians(tilt_deg(experiment))
    # The foil's normal, turned about the vertical (y) axis by the tilt; its front looks back along the beam.
    n = np.array([-math.sin(a), 0.0, -math.cos(a)])
    v = np.array([0.0, 1.0, 0.0])
    t = max(target_thickness_mm(experiment), 0.02)
    part = Part("ring", "foil", z0=t / 2, depth=t, r_out=TARGET_DIAMETER_MM / 2, hull=True, label="target foil")
    return Solid("target", f"{experiment.target.material} target", "target",
                 np.array([0.0, 0.0, experiment.target.position_mm]), _frame(np.cross(v, n), v, n), [part])


def solids(experiment) -> list:
    """The target and every detector of a setup as :class:`Solid` s."""
    out = [_target_solid(experiment)]
    out += [_detector_solid(d, i) for i, d in enumerate(experiment.detectors)]
    out += [_gamma_solid(g, i) for i, g in enumerate(experiment.gamma_detectors)]
    return out


def beam_radius_mm(experiment) -> float:
    """Radius kept clear around the beam axis, mm: the beam spot's FWHM, and at least ``BEAM_RADIUS_MM``."""
    from .rates import FWHM_PER_SIGMA, spot_sigma_mm

    return max(BEAM_RADIUS_MM, spot_sigma_mm(experiment) * FWHM_PER_SIGMA)


def extent_mm(experiment, all_solids: Optional[list] = None) -> float:
    """Half the size of the scene, mm: a little beyond the farthest solid, or the chamber."""
    far = max((float(np.linalg.norm(s.centre)) + s.reach() for s in all_solids or solids(experiment)), default=100.0)
    ch = experiment.chamber
    if ch is not None:
        far = max(far, _mm(ch.radius) + _mm(ch.wall_thickness))
    return 1.15 * far


def _beam_lines(experiment) -> np.ndarray:
    r = beam_radius_mm(experiment)
    t = np.linspace(0.0, 2 * math.pi, 9)[:-1]
    xy = np.r_[[[0.0, 0.0]], np.c_[r * np.cos(t), r * np.sin(t)]]
    return np.stack([np.c_[xy, np.full(len(xy), -1e4)], np.c_[xy, np.full(len(xy), 1e4)]], axis=1)


@dataclass
class Problem:
    """Something a layout may not do. ``keys`` are the solids involved; ``kind`` is "beam" or "overlap"."""

    kind: str
    keys: tuple
    text: str


def _overlap(a: Solid, b: Solid) -> bool:
    if float(np.linalg.norm(a.centre - b.centre)) > a.reach() + b.reach():
        return False
    return b.crossed_by(a.edges()) or a.crossed_by(b.edges())


def problems(experiment, only: Optional[str] = None, all_solids: Optional[list] = None) -> list:
    """Detectors standing in the beam, and solids that overlap. With ``only`` (a solid's key), just the problems
    that involve that solid."""
    ss = all_solids or solids(experiment)
    beam = _beam_lines(experiment)
    out = []
    for s in ss:
        if s.kind != "target" and (only is None or s.key == only) and s.crossed_by(beam):
            out.append(Problem("beam", (s.key,), f"{s.name} would stand in the beam."))
    for i, a in enumerate(ss):
        for b in ss[i + 1:]:
            if only is not None and only not in (a.key, b.key):
                continue
            if _overlap(a, b):
                out.append(Problem("overlap", (a.key, b.key), f"{a.name} would overlap {b.name}."))
    return out


# -- moving a detector ----------------------------------------------------------------------------------------------


def _entry(experiment_dict: dict, key: str) -> dict:
    kind, i = key.split(":")
    return experiment_dict["detectors" if kind == "detector" else "gamma_detectors"][int(i)]


def _setup_object(experiment, key: str):
    kind, i = key.split(":")
    return (experiment.detectors if kind == "detector" else experiment.gamma_detectors)[int(i)]


def position_of(experiment, key: str) -> np.ndarray:
    """Centre of a detector's front face, mm, from the chamber's centre (the scene's coordinates)."""
    obj = _setup_object(experiment, key)
    if key.startswith("gamma"):
        return np.array(obj.chamber_centre_mm(), dtype=float)
    return np.array(obj.chamber_position_mm(), dtype=float)


def on_axis(experiment, key: str) -> bool:
    """Whether a detector sits on the beam axis (an annular detector around the beam): it then slides along it."""
    p = position_of(experiment, key)
    return math.hypot(p[0], p[1]) < 1e-6 * max(1.0, abs(p[2]))


def _num(x: float, digits: int = 1) -> str:
    return f"{round(float(x), digits):g}"


def move_fields(experiment, key: str, position, mode: str = "angle") -> dict:
    """The setup fields of detector ``key`` when its front face is moved to ``position`` (mm).

    ``mode`` "angle" keeps the distance and takes the direction of ``position``; "distance" keeps the direction and
    takes the distance along it. A detector on the beam axis only slides along the axis. Angles and distances are
    rounded to 0.1° and 0.1 mm. A detector placed by ``position`` in the setup file is moved to the new position
    as it is; a value of ``None`` means the field is to be removed."""
    obj = _setup_object(experiment, key)
    new = np.asarray(position, dtype=float)
    old = position_of(experiment, key)
    r_old = float(np.linalg.norm(old))
    u_old = old / r_old
    if on_axis(experiment, key) or mode == "distance":
        r = max(float(new @ u_old), 0.1)
        if key.startswith("detector") and obj.position is not None:
            return {"position": [f"{_num(c, 2)} mm" for c in r * u_old]}
        return {"distance": f"{_num(r)} mm"}
    u = new / float(np.linalg.norm(new))
    theta = math.degrees(math.acos(max(-1.0, min(1.0, u[2]))))
    phi = math.degrees(math.atan2(u[1], u[0]))
    if key.startswith("detector") and obj.position is not None:
        return {"position": [f"{_num(c, 2)} mm" for c in r_old * u]}
    return {"theta": f"{_num(theta)} deg", "phi": f"{_num(phi)} deg"}


def moved(experiment, key: str, fields: dict):
    """A copy of the experiment with ``fields`` set on detector ``key``. Raises
    :class:`~physim.nuclear.experiment.SetupError` if the setup is then not valid (outside the chamber, say)."""
    from .experiment import Experiment

    d = experiment.to_dict()
    entry = _entry(d, key)
    for k, v in fields.items():
        if v is None:
            entry.pop(k, None)
        else:
            entry[k] = v
    return Experiment.from_dict(copy.deepcopy(d))


def blocked(experiment, key: str, fields: dict) -> Optional[str]:
    """Why detector ``key`` may not be given ``fields`` (it would stand in the beam, overlap another solid or
    cross the chamber wall), or ``None`` if it may."""
    from .experiment import SetupError

    try:
        trial = moved(experiment, key, fields)
    except SetupError as e:
        return "; ".join(e.problems)
    found = problems(trial, only=key)
    return found[0].text if found else None


@dataclass
class Limits:
    """How far a detector can be dragged. ``mode`` is "distance" (along ``direction``, between ``lo`` and ``hi``
    mm from the target) or "angle" (over a sphere of ``radius`` mm, with the polar angle between ``lo`` and ``hi``
    degrees). ``lo_why`` and ``hi_why`` say what stops it at each end."""

    mode: str
    lo: float
    hi: float
    direction: tuple = (0.0, 0.0, 1.0)
    radius: float = 0.0
    lo_why: str = ""
    hi_why: str = ""

    def constraints(self) -> str:
        """The limits as the scene's drag constraints: assignments to x, y and z, applied in order to the dragged
        object's position (NiceGUI evaluates each in the browser).

        The expressions avoid commas and any letter x, y or z other than the coordinates, because of how the scene
        splits and substitutes them."""
        f = lambda v: f"({float(v):.6f})"  # noqa: E731
        if self.mode == "distance":
            u = self.direction
            k = max(range(3), key=lambda i: abs(u[i]))
            names = "xyz"
            dot = f"(x * {f(u[0])} + y * {f(u[1])} + z * {f(u[2])})"
            first = f"{names[k]} = {dot} < {f(self.lo)} ? {f(self.lo)} : {dot} > {f(self.hi)} ? {f(self.hi)} : {dot}"
            rest = [f"{names[i]} = {f(u[i])} * {names[k]}" for i in range(3) if i != k]
            return ", ".join([first, *rest, f"{names[k]} = {f(u[k])} * {names[k]}"])
        c_lo, c_hi = math.cos(math.radians(self.hi)), math.cos(math.radians(self.lo))
        cos = "(z / Math.sqrt(x * x + y * y + z * z + 1e-9))"
        r = f(self.radius)
        return ", ".join([
            f"z = {cos} < {f(c_lo)} ? {f(c_lo)} : {cos} > {f(c_hi)} ? {f(c_hi)} : {cos}",
            f"x = x * {r} * Math.sqrt((1 - z * z) / (x * x + y * y + 1e-9))",
            f"y = (y < 0 ? -1 : 1) * Math.sqrt(Math.abs({r} * {r} * (1 - z * z) - x * x))",
            f"z = {r} * z"])

    def apply(self, position) -> np.ndarray:
        """The position the constraints give for a dragged ``position``: the same rule as :meth:`constraints`."""
        p = np.asarray(position, dtype=float)
        if self.mode == "distance":
            u = np.array(self.direction)
            return min(max(float(p @ u), self.lo), self.hi) * u
        c = p[2] / math.sqrt(float(p @ p) + 1e-9)
        c = min(max(c, math.cos(math.radians(self.hi))), math.cos(math.radians(self.lo)))
        rho = math.hypot(p[0], p[1])
        s = self.radius * math.sqrt(1 - c * c)
        if rho < 1e-9:
            return np.array([s, 0.0, self.radius * c])
        return np.array([s * p[0] / rho, s * p[1] / rho, self.radius * c])

    def stopped(self, position, tolerance: float = 0.05) -> Optional[str]:
        """Why a drag to ``position`` was held back, or ``None`` if it was not."""
        p = np.asarray(position, dtype=float)
        if self.mode == "distance":
            d = float(p @ np.array(self.direction))
        else:
            d = math.degrees(math.acos(max(-1.0, min(1.0, p[2] / (float(np.linalg.norm(p)) or 1.0)))))
        if d <= self.lo + tolerance:
            return self.lo_why or None
        if d >= self.hi - tolerance:
            return self.hi_why or None
        return None


def _scan(test, start: float, stop: float, step: float) -> tuple:
    """Walk from ``start`` towards ``stop`` until ``test(value)`` gives a reason; refine the last allowed value by
    bisection. Returns (last allowed value, reason or "")."""
    sign = 1.0 if stop >= start else -1.0
    good, x = start, start
    while (stop - x) * sign > 1e-9:
        x = x + sign * min(step, abs(stop - x))
        why = test(x)
        if why:
            lo, hi = good, x
            for _ in range(6):
                mid = (lo + hi) / 2
                if test(mid):
                    hi = mid
                else:
                    lo = mid
            return lo, why
        good = x
    return good, ""


def limits(experiment, key: str, mode: str = "angle") -> Limits:
    """How far detector ``key`` can be dragged from where it is, in ``mode`` "angle" or "distance" (a detector on
    the beam axis always moves in distance).

    In distance, the limits are where it would meet the target, another solid, the beam or the chamber wall. In
    angle, they are where it would enter the beam; overlaps with other detectors are checked as it moves (see
    :func:`blocked`), since they depend on both angles."""
    p = position_of(experiment, key)
    r = float(np.linalg.norm(p))
    u = p / r
    if on_axis(experiment, key) or mode == "distance":
        far = r + 200.0
        if experiment.chamber is not None and key.startswith("detector"):
            far = _mm(experiment.chamber.radius)

        def test(d: float) -> Optional[str]:
            return blocked(experiment, key, move_fields(experiment, key, d * u, "distance"))

        lo, lo_why = _scan(test, r, MIN_DISTANCE_MM, 2.0)
        hi, hi_why = _scan(test, r, far, 10.0)
        return Limits("distance", lo, hi, tuple(u), r, lo_why or "It cannot come closer to the target than "
                      f"{MIN_DISTANCE_MM:g} mm.", hi_why or "That is as far as it can be dragged; type a "
                      "larger distance if you need one.")
    theta = math.degrees(math.acos(max(-1.0, min(1.0, u[2]))))
    phi = math.atan2(u[1], u[0])

    def test(th: float) -> Optional[str]:
        t = math.radians(th)
        q = r * np.array([math.sin(t) * math.cos(phi), math.sin(t) * math.sin(phi), math.cos(t)])
        trial = moved(experiment, key, move_fields(experiment, key, q, "angle"))
        found = [x for x in problems(trial, only=key) if x.kind == "beam"]
        return found[0].text if found else None

    lo, lo_why = _scan(test, theta, 0.5, 4.0)
    hi, hi_why = _scan(test, theta, 179.5, 4.0)
    return Limits("angle", lo, hi, tuple(u), r, lo_why or "It would stand in the beam.",
                  hi_why or "It would stand in the beam.")


# -- advice ---------------------------------------------------------------------------------------------------------


def advice(experiment) -> list:
    """Where the kinematics send the particles, in plain words: which side to put particle detectors on, and which
    detectors nothing can reach."""
    from . import data
    from .kinematics import TwoBody
    from .rates import beam_ion, stack

    ion = beam_ion(experiment)
    a_beam = experiment.beam.A
    targets = [data.nuclide(k) for lay in stack(experiment)[:1] for k in lay.nuclides]
    heavy = max(targets, key=lambda n: n.A)
    out = []
    if a_beam < heavy.A:
        out.append(f"{ion} is lighter than {heavy.name}: the beam scatters to every angle. Backward angles are "
                   "recommended for the particle detector: the excitation is strongest there and the elastic rate "
                   "lowest.")
    else:
        tb = TwoBody(ion, heavy.name, experiment.beam.energy_mev)
        out.append(f"{ion} is heavier than {heavy.name}: the scattered beam stays within {tb.max_angle():.0f}° of "
                   "the beam and the recoils within 90°. Forward angles are recommended for the particle detector.")
    reach = 0.0
    for lay in stack(experiment):
        for k in lay.nuclides:
            tb = TwoBody(ion, data.nuclide(k).name, experiment.beam.energy_mev)
            reach = max(reach, tb.max_angle(), tb.max_angle("recoil"))
    for i, det in enumerate(experiment.detectors, start=1):
        lo, _ = Geometry.from_setup(det, i).theta_range()
        if lo > reach + 1e-6:
            out.append(f"{det.name or f'D{i}'} starts at {lo:.0f}°, but no particle of this reaction goes beyond "
                       f"{reach:.0f}°: it will count nothing.")
    return out


__all__ = ["BEAM_RADIUS_MM", "BOARD_THICKNESS_MM", "Limits", "MIN_DISTANCE_MM", "Part", "Problem", "Solid",
           "TARGET_DIAMETER_MM", "advice", "beam_radius_mm", "blocked", "extent_mm", "limits", "move_fields",
           "moved", "on_axis", "position_of", "problems", "solids", "target_thickness_mm"]
