"""Helpers for setting up boundaries and bodies."""

import math


def box_walls(lo, hi):
    """The six walls of the box ``[lo, hi]`` as ``(normal, offset)`` pairs for
    :meth:`physim.World.set_collisions`, with normals pointing into the box."""
    walls = []
    for axis in range(3):
        n = [0.0, 0.0, 0.0]
        n[axis] = 1.0
        walls.append((list(n), float(lo[axis])))
        n[axis] = -1.0
        walls.append((list(n), -float(hi[axis])))
    return walls


def inertia_box(mass, a, b, c):
    """Principal moments of a solid box with sides ``a``, ``b``, ``c`` along x, y, z."""
    return [mass * (b * b + c * c) / 12, mass * (a * a + c * c) / 12, mass * (a * a + b * b) / 12]


def inertia_ellipsoid(mass, a, b, c):
    """Principal moments of a solid ellipsoid with semi-axes ``a``, ``b``, ``c``."""
    return [mass * (b * b + c * c) / 5, mass * (a * a + c * c) / 5, mass * (a * a + b * b) / 5]


def inertia_sphere(mass, radius):
    """Principal moments of a solid sphere."""
    i = 0.4 * mass * radius * radius
    return [i, i, i]


def inertia_cylinder(mass, radius, length):
    """Principal moments of a solid cylinder with its axis along z."""
    t = mass * (3 * radius * radius + length * length) / 12
    return [t, t, 0.5 * mass * radius * radius]


def quaternion_from_axis_angle(axis, angle):
    """Unit quaternion ``(w, x, y, z)`` rotating by ``angle`` (radians) about ``axis``."""
    n = math.sqrt(sum(float(a) * float(a) for a in axis))
    if n == 0:
        return [1.0, 0.0, 0.0, 0.0]
    s = math.sin(0.5 * angle) / n
    return [math.cos(0.5 * angle)] + [float(a) * s for a in axis]
