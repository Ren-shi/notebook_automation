"""Detector geometry and response: where each detector sits, what it sees, and how it turns a particle's energy
into a measured signal.

Coordinates as in the setup file: the beam runs along +z through the target at the origin, x horizontal, y up;
lengths in mm, angles in degrees, solid angles in msr::

    from physim.nuclear import Experiment
    from physim.nuclear.detectors import Array

    array = Array.from_experiment(Experiment.example("oxygen_on_lead_array"))
    d = array["DSSD1"]
    d.solid_angle()                 # msr, by numerical integration over the face
    d.segment_solid_angles()        # per strip pair / ring-sector, msr
    d.theta_range()                 # (min, max) polar angle covered, degrees
    d.hit([0.5, 0.0, 0.85])         # which segment a direction from the target hits, or None
    array.warnings()                # detectors shadowing each other or sitting in the beam

    resp = d.response("16O")
    resp.measured_energy(60.0, incidence_deg=10.0)   # after the dead layer, punch-through, threshold
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .quantity import Quantity

#: Gauss–Legendre order used for the face integrals (per dimension, per segment).
QUADRATURE_ORDER = 48


def _q(x) -> Optional[Quantity]:
    if x is None:
        return None
    return Quantity.parse(x) if isinstance(x, str) else x


def _unit(v) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def direction(theta_deg: float, phi_deg: float = 0.0) -> np.ndarray:
    """Unit vector at polar angle θ (from +z) and azimuth φ (from +x towards +y)."""
    t, p = math.radians(theta_deg), math.radians(phi_deg)
    return np.array([math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t)])


def angles(v) -> tuple:
    """(θ, φ) in degrees of a vector; φ in (−180, 180]."""
    v = np.asarray(v, dtype=float)
    r = np.linalg.norm(v, axis=-1)
    return np.degrees(np.arccos(np.clip(v[..., 2] / r, -1, 1))), np.degrees(np.arctan2(v[..., 1], v[..., 0]))


@dataclass
class Hit:
    """Where a straight line from the source crosses a detector face."""

    detector: str
    #: (strip along u, strip along v) for rectangles, (ring, sector) for annular detectors, (0, 0) for circles.
    segment: tuple
    #: Distance from the source, mm.
    distance: float
    #: Angle between the track and the face normal, degrees.
    incidence: float
    #: Position in the face's own coordinates (u, v), mm.
    local: tuple


class Geometry:
    """One detector's face: a flat rectangle, disc or annulus with a centre, a normal and in-plane axes u, v.

    The normal ``n`` points from the face towards the target. ``u`` is horizontal and ``v`` "up" in the face (as
    close to +y as the face allows), turned by ``rotation`` about ``n``. Rectangles are divided into ``strips_x``
    along u and ``strips_y`` along v; annular detectors into ``rings`` (equal radial width) and ``sectors`` (equal
    angle, counted from u towards v).
    """

    def __init__(self, name: str, shape: str, centre, normal, rotation_deg: float = 0.0, *, width=None, height=None,
                 strips_x: int = 1, strips_y: int = 1, radius=None, inner_radius=None, outer_radius=None,
                 rings: int = 1, sectors: int = 1):
        self.name, self.shape = name, shape
        self.centre = np.asarray(centre, dtype=float)
        self.n = _unit(normal)
        ref = np.array([0.0, 1.0, 0.0]) if abs(self.n[1]) < 0.99 else np.array([1.0, 0.0, 0.0])
        v = _unit(ref - ref.dot(self.n) * self.n)
        u = np.cross(v, self.n)
        r = math.radians(rotation_deg)
        self.u = math.cos(r) * u + math.sin(r) * v
        self.v = -math.sin(r) * u + math.cos(r) * v
        self.width, self.height = width, height
        self.strips_x, self.strips_y = int(strips_x), int(strips_y)
        self.radius = radius
        self.inner_radius, self.outer_radius = inner_radius, outer_radius
        self.rings, self.sectors = int(rings), int(sectors)
        if shape == "circle":
            self.inner_radius, self.outer_radius, self.rings, self.sectors = 0.0, radius, 1, 1
        if shape not in ("rectangle", "circle", "annular"):
            raise ValueError(f"unknown shape '{shape}'")

    @classmethod
    def from_setup(cls, det, index: int = 1) -> Geometry:
        """From a :class:`~physim.nuclear.Detector` of a setup file."""
        mm = lambda x: None if x is None else _q(x).to("mm")  # noqa: E731
        centre = det.position_mm()
        rotation = _q(det.rotation).to("deg") if det.rotation is not None else 0.0
        return cls(det.name or f"D{index}", det.shape, centre, det.facing_direction(), rotation,
                   width=mm(det.width), height=mm(det.height), strips_x=det.strips_x or 1,
                   strips_y=det.strips_y or 1, radius=mm(det.radius), inner_radius=mm(det.inner_radius),
                   outer_radius=mm(det.outer_radius), rings=det.rings or 1, sectors=det.sectors or 1)

    # -- shape ------------------------------------------------------------------------------------------------

    @property
    def segments(self) -> list:
        if self.shape == "rectangle":
            return [(i, j) for i in range(self.strips_x) for j in range(self.strips_y)]
        return [(i, j) for i in range(self.rings) for j in range(self.sectors)]

    def _segment_of(self, su: np.ndarray, sv: np.ndarray) -> tuple:
        """Segment indices of local points, and whether they lie on the face."""
        if self.shape == "rectangle":
            a, b = self.width / 2, self.height / 2
            inside = (np.abs(su) <= a) & (np.abs(sv) <= b)
            i = np.clip(((su + a) / self.width * self.strips_x).astype(int), 0, self.strips_x - 1)
            j = np.clip(((sv + b) / self.height * self.strips_y).astype(int), 0, self.strips_y - 1)
            return inside, i, j
        rr = np.hypot(su, sv)
        inside = (rr >= self.inner_radius) & (rr <= self.outer_radius)
        i = np.clip(((rr - self.inner_radius) / (self.outer_radius - self.inner_radius) * self.rings).astype(int),
                    0, self.rings - 1)
        ang = np.mod(np.arctan2(sv, su), 2 * math.pi)
        j = np.clip((ang / (2 * math.pi) * self.sectors).astype(int), 0, self.sectors - 1)
        return inside, i, j

    def _boundary(self, s: np.ndarray, inner: bool = False) -> np.ndarray:
        """Points on the edge of the face for s in [0, 1) (3D, mm): around the rectangle, or around the outer (or
        inner) circle."""
        s = np.asarray(s, dtype=float) % 1.0
        if self.shape == "rectangle":
            a, b = self.width / 2, self.height / 2
            corners = np.array([[-a, -b], [a, -b], [a, b], [-a, b], [-a, -b]])
            k = np.minimum((s * 4).astype(int), 3)
            f = s * 4 - k
            loc = corners[k] + f[:, None] * (corners[k + 1] - corners[k])
        else:
            r = self.inner_radius if inner else self.outer_radius
            loc = np.c_[r * np.cos(2 * math.pi * s), r * np.sin(2 * math.pi * s)]
        return self.centre + loc[:, :1] * self.u + loc[:, 1:] * self.v

    def outline(self, points: int = 120) -> np.ndarray:
        """Points around the edge of the face (3D, mm), closed; annular detectors give the outer edge, a NaN row,
        then the inner edge."""
        s = np.linspace(0.0, 1.0, points + 1)
        s[-1] = 0.0
        edge = self._boundary(s)
        if self.shape == "annular" and self.inner_radius:
            edge = np.r_[edge, np.full((1, 3), np.nan), self._boundary(s, inner=True)]
        return edge

    def _extreme_on_boundary(self, f, sign: float) -> float:
        """sign × max of sign·f over the face edge(s): dense sampling, then golden-section refinement."""
        best = -np.inf
        for inner in ((False, True) if self.shape == "annular" and self.inner_radius else (False,)):
            s = np.linspace(0, 1, 4001)[:-1]
            vals = sign * f(self._boundary(s, inner))
            k = int(np.argmax(vals))
            lo, hi = s[k] - 1 / 4000, s[k] + 1 / 4000
            g = (math.sqrt(5) - 1) / 2
            for _ in range(60):
                m1, m2 = hi - g * (hi - lo), lo + g * (hi - lo)
                if sign * f(self._boundary(np.array([m1]), inner))[0] > sign * f(self._boundary(np.array([m2]), inner))[0]:
                    hi = m2
                else:
                    lo = m1
            best = max(best, float(np.max(vals)), sign * float(f(self._boundary(np.array([(lo + hi) / 2]), inner))[0]))
        return sign * best

    # -- rays -------------------------------------------------------------------------------------------------

    def hit(self, d, source=(0.0, 0.0, 0.0), two_sided: bool = False) -> Optional[Hit]:
        """Where the straight line from ``source`` along direction ``d`` crosses this face, or ``None``.

        Only the front (the side facing the target) counts unless ``two_sided``.
        """
        out = self.hits(np.atleast_2d(d), source, two_sided)
        return out[0]

    def _intersect(self, dirs: np.ndarray, source, two_sided: bool) -> tuple:
        s = np.asarray(source, dtype=float)
        d = np.asarray(dirs, dtype=float)
        d = d / np.linalg.norm(d, axis=1, keepdims=True)
        dn = d @ self.n
        with np.errstate(divide="ignore", invalid="ignore"):
            t = ((self.centre - s) @ self.n) / dn
        with np.errstate(invalid="ignore"):
            p = s + t[:, None] * d
            su, sv = (p - self.centre) @ self.u, (p - self.centre) @ self.v
        finite = np.isfinite(su) & np.isfinite(sv)  # a direction parallel to the face never crosses it
        su, sv = np.where(finite, su, np.inf), np.where(finite, sv, np.inf)
        inside, i, j = self._segment_of(np.where(finite, su, 0.0), np.where(finite, sv, 0.0))
        ok = finite & inside & (t > 0) & ((dn < 0) | two_sided)  # the face looks back towards the source
        return ok, t, dn, su, sv, i, j

    def hits(self, dirs: np.ndarray, source=(0.0, 0.0, 0.0), two_sided: bool = False) -> list:
        """:meth:`hit` for many directions at once (an (N, 3) array)."""
        ok, t, dn, su, sv, i, j = self._intersect(dirs, source, two_sided)
        inc = np.degrees(np.arccos(np.clip(-dn, -1, 1)))
        return [Hit(self.name, (int(i[k]), int(j[k])), float(t[k]), float(inc[k]), (float(su[k]), float(sv[k])))
                if ok[k] else None for k in range(len(ok))]

    def distances(self, dirs: np.ndarray, source=(0.0, 0.0, 0.0)) -> np.ndarray:
        """Distance (mm) from ``source`` along each direction to the front of this face; ``inf`` where it misses."""
        ok, t = self._intersect(dirs, source, False)[:2]
        return np.where(ok, t, np.inf)

    # -- solid angle ------------------------------------------------------------------------------------------

    def _nodes(self, segment: Optional[tuple] = None, order: int = QUADRATURE_ORDER) -> tuple:
        """Quadrature points (local u, v) and weights (area, mm²) over the face or one segment."""
        x, w = np.polynomial.legendre.leggauss(order)
        if self.shape == "rectangle":
            a, b = self.width / 2, self.height / 2
            if segment is None:
                u0, u1, v0, v1 = -a, a, -b, b
            else:
                du, dv = self.width / self.strips_x, self.height / self.strips_y
                u0, v0 = -a + segment[0] * du, -b + segment[1] * dv
                u1, v1 = u0 + du, v0 + dv
            uu = (u1 - u0) / 2 * x + (u1 + u0) / 2
            vv = (v1 - v0) / 2 * x + (v1 + v0) / 2
            U, V = np.meshgrid(uu, vv, indexing="ij")
            W = np.outer(w, w) * (u1 - u0) / 2 * (v1 - v0) / 2
            return U.ravel(), V.ravel(), W.ravel()
        r0, r1, a0, a1 = self.inner_radius, self.outer_radius, 0.0, 2 * math.pi
        if segment is not None:
            dr = (r1 - r0) / self.rings
            r0, r1 = r0 + segment[0] * dr, r0 + (segment[0] + 1) * dr
            da = 2 * math.pi / self.sectors
            a0, a1 = segment[1] * da, (segment[1] + 1) * da
        rr = (r1 - r0) / 2 * x + (r1 + r0) / 2
        aa = (a1 - a0) / 2 * x + (a1 + a0) / 2
        R, A = np.meshgrid(rr, aa, indexing="ij")
        W = np.outer(w, w) * (r1 - r0) / 2 * (a1 - a0) / 2 * R
        return (R * np.cos(A)).ravel(), (R * np.sin(A)).ravel(), W.ravel()

    def segment_centre(self, segment: Optional[tuple] = None, weighted: bool = False) -> np.ndarray:
        """The middle of a segment (or of the face), mm: what a Doppler correction knows of where a particle hit.
        With ``weighted`` it is the mean point of the segment weighted by the solid angle and by Rutherford's
        1/sin⁴(θ/2), nearer to where the particles of a steep distribution hit."""
        if weighted and segment is not None:
            dirs, d_omega = self.directions(segment, 6)
            theta = np.arccos(np.clip(dirs[:, 2], -1.0, 1.0))
            w = d_omega / np.maximum(np.sin(theta / 2), 1e-3) ** 4
            su, sv, _ = self._nodes(segment, 6)
            return self.centre + float(np.average(su, weights=w)) * self.u + float(np.average(sv, weights=w)) * self.v
        if segment is None:
            return self.centre.copy()
        if self.shape == "rectangle":
            du, dv = self.width / self.strips_x, self.height / self.strips_y
            su, sv = -self.width / 2 + (segment[0] + 0.5) * du, -self.height / 2 + (segment[1] + 0.5) * dv
        else:
            dr, da = (self.outer_radius - self.inner_radius) / self.rings, 2 * math.pi / self.sectors
            r, a = self.inner_radius + (segment[0] + 0.5) * dr, (segment[1] + 0.5) * da
            su, sv = r * math.cos(a), r * math.sin(a)
        return self.centre + su * self.u + sv * self.v

    def solid_angle(self, segment: Optional[tuple] = None, source=(0.0, 0.0, 0.0)) -> float:
        """Solid angle seen from ``source``, msr: ∫ (r̂ · n) dA / r² over the face (Gauss–Legendre)."""
        _, d_omega = self.directions(segment, source=source)
        return float(np.sum(d_omega)) * 1e3

    def directions(self, segment: Optional[tuple] = None, order: int = QUADRATURE_ORDER,
                   source=(0.0, 0.0, 0.0)) -> tuple:
        """Quadrature over the face (or one segment) as seen from ``source``: unit directions (N, 3) and the solid
        angle each stands for (sr), so ∫ f dΩ ≈ Σ f(direction) dΩ."""
        su, sv, w = self._nodes(segment, order)
        p = self.centre + su[:, None] * self.u + sv[:, None] * self.v - np.asarray(source, dtype=float)
        r = np.linalg.norm(p, axis=1)
        cos = -(p @ self.n) / r
        return p / r[:, None], w * np.clip(cos, 0, None) / r**2

    def segment_solid_angles(self, source=(0.0, 0.0, 0.0)) -> dict:
        """Solid angle of every segment, msr."""
        return {seg: self.solid_angle(seg, source) for seg in self.segments}

    def solid_angle_monte_carlo(self, n: int = 1_000_000, seed: int = 1) -> tuple:
        """Independent Monte Carlo estimate (msr) and its standard error: isotropic directions inside the cone
        that just contains the face, counted when they hit it."""
        rng = np.random.default_rng(seed)
        axis = _unit(self.centre)
        corners = self.outline(400)
        corners = corners[~np.isnan(corners[:, 0])]
        half = float(np.max(np.arccos(np.clip(_unit_rows(corners) @ axis, -1, 1))))
        half = min(math.pi, half * 1.02 + 1e-6)
        cos_min = math.cos(half)
        c = rng.uniform(cos_min, 1.0, n)
        a = rng.uniform(0, 2 * math.pi, n)
        s = np.sqrt(1 - c**2)
        e1 = _unit(np.cross(axis, [0.0, 0.0, 1.0] if abs(axis[2]) < 0.9 else [1.0, 0.0, 0.0]))
        e2 = np.cross(axis, e1)
        d = c[:, None] * axis + s[:, None] * (np.cos(a)[:, None] * e1 + np.sin(a)[:, None] * e2)
        dn = d @ self.n
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (self.centre @ self.n) / dn
        p = t[:, None] * d
        inside, _, _ = self._segment_of((p - self.centre) @ self.u, (p - self.centre) @ self.v)
        frac = np.mean(inside & (t > 0) & (dn < 0))
        cone = 2 * math.pi * (1 - cos_min)
        return cone * frac * 1e3, cone * math.sqrt(frac * (1 - frac) / n) * 1e3

    # -- angular coverage -------------------------------------------------------------------------------------

    def _face_points(self, segment: Optional[tuple] = None) -> np.ndarray:
        su, sv, _ = self._nodes(segment)
        pts = self.centre + su[:, None] * self.u + sv[:, None] * self.v
        if segment is not None:
            return pts
        edge = self.outline(400)
        return np.r_[pts, edge[~np.isnan(edge[:, 0])]]

    def theta_range(self, segment: Optional[tuple] = None) -> tuple:
        """(smallest, largest) polar angle on the face (or a segment), degrees; 0 or 180 if the beam axis
        crosses it. For the whole face the extremes lie on its edge and are found to ~1e-9°; for a segment they
        are taken over its quadrature points (to a fraction of the segment's size)."""
        theta = lambda p: angles(p)[0]  # noqa: E731
        if segment is None:
            lo = self._extreme_on_boundary(theta, -1.0)
            hi = self._extreme_on_boundary(theta, +1.0)
        else:
            th = theta(self._face_points(segment))
            lo, hi = float(th.min()), float(th.max())
        for axis_dir, val in (((0, 0, 1), 0.0), ((0, 0, -1), 180.0)):
            h = self.hit(np.array(axis_dir, dtype=float))
            if h is not None and (segment is None or h.segment == segment):
                lo, hi = (val, hi) if val == 0.0 else (lo, val)
        return lo, hi

    def phi_range(self, segment: Optional[tuple] = None) -> tuple:
        """(smallest, largest) azimuth on the face, degrees, unwrapped around the centre's azimuth;
        (−180, 180) if the face surrounds the beam axis."""
        if self.hit(np.array([0.0, 0.0, 1.0])) is not None or self.hit(np.array([0.0, 0.0, -1.0])) is not None:
            if self.shape != "annular" or self.inner_radius == 0:
                return -180.0, 180.0
        pts = self._face_points(segment)
        _, ph = angles(pts)
        _, ph0 = angles(self.centre if segment is None else pts.mean(axis=0))
        rel = (ph - ph0 + 180.0) % 360.0 - 180.0
        if np.ptp(rel) > 350:
            return -180.0, 180.0
        return float(ph0 + rel.min()), float(ph0 + rel.max())

    def mean_theta(self, segment: Optional[tuple] = None) -> float:
        """Solid-angle-weighted mean polar angle, degrees."""
        su, sv, w = self._nodes(segment)
        p = self.centre + su[:, None] * self.u + sv[:, None] * self.v
        r = np.linalg.norm(p, axis=1)
        dom = w * np.clip(-(p @ self.n) / r, 0, None) / r**2
        th, _ = angles(p)
        return float(np.sum(th * dom) / np.sum(dom))


def _unit_rows(a: np.ndarray) -> np.ndarray:
    return a / np.linalg.norm(a, axis=1, keepdims=True)


# ---------------------------------------------------------------------------------------------------------------
# Response


class Response:
    """How a detector turns a particle's energy into a measured energy.

    The particle crosses the dead layer, then deposits energy in the active thickness: all of it if it stops,
    the difference if it punches through. A Gaussian resolution (FWHM) and a threshold are applied last. Path
    lengths grow as 1/cos(incidence).
    """

    def __init__(self, ion: str, setup_detector, stopping_cls=None):
        from .stopping import Stopping

        cls = stopping_cls or Stopping
        det = setup_detector
        self.ion = ion
        self.material = det.material_data()
        self.stopping = cls(ion, self.material)
        self.dead_layer = self._areal(det.dead_layer)
        self.thickness = self._areal(det.thickness)
        self.fwhm = _q(det.resolution).to("MeV") if det.resolution is not None else 0.0
        self.threshold = _q(det.threshold).to("MeV") if det.threshold is not None else 0.0

    def _areal(self, x) -> float:
        if x is None:
            return 0.0
        return self.material.areal_density_mg_cm2(_q(x))

    def after_dead_layer(self, energy, incidence_deg: float = 0.0):
        """Energy entering the active volume, MeV."""
        if self.dead_layer == 0:
            return energy
        return self.stopping.energy_after(energy, self.dead_layer, incidence_deg)

    def deposited(self, energy, incidence_deg: float = 0.0):
        """Energy deposited in the active volume, MeV (all of it unless the particle punches through)."""
        e_in = np.asarray(self.after_dead_layer(energy, incidence_deg), dtype=float)
        out = e_in - self.stopping.energy_after(e_in, self.thickness, incidence_deg)
        return float(out) if np.ndim(energy) == 0 else out

    def punch_through_energy(self, incidence_deg: float = 0.0) -> float:
        """Smallest incident energy that leaves the back of the detector, MeV."""
        path = (self.dead_layer + self.thickness) / math.cos(math.radians(incidence_deg))
        rng = self.stopping._range
        return float(np.exp(np.interp(math.log(path), np.log(rng), np.log(self.stopping._e))))

    def measured_energy(self, energy, incidence_deg: float = 0.0, rng: Optional[np.random.Generator] = None):
        """Measured energy, MeV: deposited energy smeared by the resolution (if ``rng`` is given), NaN below the
        threshold."""
        e = np.asarray(self.deposited(energy, incidence_deg), dtype=float)
        if rng is not None and self.fwhm > 0:
            e = e + rng.normal(0.0, self.fwhm / 2.3548, size=e.shape)
        e = np.where(e >= self.threshold, e, np.nan)
        return float(e) if np.ndim(energy) == 0 and rng is None else e


# ---------------------------------------------------------------------------------------------------------------
# A whole array


#: How far behind its active face a detector model's dead material is placed, mm, so that the active face is met
#: first where the two overlap.
BEHIND_MM = 0.01


class Array:
    """All detectors of an experiment."""

    def __init__(self, geometries: Sequence[Geometry], setup_detectors: Optional[Sequence] = None,
                 blockers: Sequence[Geometry] = ()):
        self.geometries = list(geometries)
        self._setup = list(setup_detectors) if setup_detectors is not None else [None] * len(self.geometries)
        #: Faces of dead material: they stop particles and count nothing (circuit boards, γ-detector housings).
        self.blockers = list(blockers)

    @classmethod
    def from_experiment(cls, experiment) -> Array:
        """The particle detectors of a setup, with the dead material that can hide them: the boards of detector
        models, placed just behind their active faces, and the front windows of γ-detector housings."""
        dets = experiment.detectors
        faces = [Geometry.from_setup(d, i) for i, d in enumerate(dets, start=1)]
        blockers = []
        for d, g in zip(dets, faces):
            for name, r_in, r_out in d.blocking():
                blockers.append(Geometry(f"{g.name} ({name})", "annular", g.centre - BEHIND_MM * g.n, g.n,
                                         inner_radius=r_in, outer_radius=r_out))
        for i, gd in enumerate(getattr(experiment, "gamma_detectors", ()), start=1):
            housing = gd.housing()
            if housing is not None:
                centre, side = housing
                blockers.append(Geometry(f"{gd.name or f'G{i}'} (housing)", "rectangle", centre,
                                         [-c for c in centre], width=side, height=side))
        return cls(faces, dets, blockers)

    def __getitem__(self, name: str) -> Geometry:
        for g in self.geometries:
            if g.name == name:
                return g
        raise KeyError(name)

    def __iter__(self):
        return iter(self.geometries)

    def response(self, name: str, ion: str) -> Response:
        i = [g.name for g in self.geometries].index(name)
        if self._setup[i] is None:
            raise ValueError("this array was not built from a setup file, so detector response is unknown")
        return Response(ion, self._setup[i])

    def first_hit(self, d, source=(0.0, 0.0, 0.0)) -> Optional[Hit]:
        """The nearest detector face the line from ``source`` along ``d`` crosses."""
        hits = [h for g in self.geometries if (h := g.hit(d, source)) is not None]
        best = min(hits, key=lambda h: h.distance) if hits else None
        if best is not None and any((h := b.hit(d, source)) is not None and h.distance < best.distance
                                    for b in self.blockers):
            return None  # dead material stops the particle first
        return best

    def visible(self, name: str, dirs: np.ndarray, source=(0.0, 0.0, 0.0)) -> np.ndarray:
        """Whether a particle leaving ``source`` along each direction reaches detector ``name`` before any other
        detector's face (the nearest face a particle meets is the one that stops it, as in the event generator)."""
        g = self[name]
        own = g.distances(dirs, source)
        ok = np.ones(len(own), dtype=bool)
        for other in self.geometries + self.blockers:
            if other is not g:
                ok &= ~(other.distances(dirs, source) < own - 1e-9)
        return ok

    def shadowing(self) -> dict:
        """{(front, behind): fraction of `behind`'s solid angle hidden by `front`} for every pair that overlaps."""
        out = {}
        for g in self.geometries:
            su, sv, w = g._nodes()
            pts = g.centre + su[:, None] * g.u + sv[:, None] * g.v
            r = np.linalg.norm(pts, axis=1)
            dom = w * np.clip(-(pts @ g.n) / r, 0, None) / r**2
            dirs = pts / r[:, None]
            for other in self.geometries + self.blockers:
                if other is g:
                    continue
                hidden = other.distances(dirs) < r - 1e-9
                frac = float(np.sum(dom[hidden]) / np.sum(dom))
                if frac > 1e-4:
                    out[(other.name, g.name)] = frac
        return out

    def warnings(self) -> list:
        """Plain-language warnings about the layout."""
        out = []
        for (front, behind), frac in self.shadowing().items():
            amount = f"{100 * frac:.0f}%" if frac >= 0.01 else "<1%"
            out.append(f"{front} hides {amount} of {behind} from the target.")
        beam = np.array([0.0, 0.0, 1.0])
        for g in self.geometries + self.blockers:
            if g.hit(beam, two_sided=True) is not None:
                out.append(f"{g.name} sits in the beam path downstream of the target (the unscattered beam hits it).")
            h = g.hit(beam, source=(0.0, 0.0, -1e6), two_sided=True)
            if h is not None and h.distance < 1e6:
                out.append(f"{g.name} blocks the incoming beam upstream of the target.")
        return out


def exit_path(direction_vec, thickness: float, depth: float, tilt_deg: float = 0.0,
              max_factor: float = 1e3) -> tuple:
    """Path (same units as ``thickness``) from a reaction at ``depth`` into a flat target to its surface, along a
    track, and whether the track runs so close to the target plane that the path was capped.

    The target normal is +z turned by ``tilt_deg`` about the vertical (y) axis. Tracks going forward through the
    target leave by the back face, backward ones by the front face. Near 90° to the normal the path grows as
    1/cos; it is capped at ``max_factor`` × thickness (a real target is finite).
    """
    d = _unit(direction_vec)
    t = math.radians(tilt_deg)
    normal = np.array([math.sin(t), 0.0, math.cos(t)])
    c = float(d @ normal)
    remaining = (thickness - depth) if c > 0 else depth
    if abs(c) < 1e-300 or remaining / abs(c) > max_factor * thickness:
        return max_factor * thickness, True
    return remaining / abs(c), False
