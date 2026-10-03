"""Helpers for setting up boundaries."""


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
