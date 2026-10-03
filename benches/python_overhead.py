"""Python-side overhead benchmarks. Run after `maturin develop --release`:

    python benches/python_overhead.py

Reports per-step times for:
- stepping from a Python loop vs letting Rust loop (`w.step(dt, n)`),
- `CustomForce` (numpy, called from Rust every force evaluation) vs the built-in Rust force.
"""

import time

import numpy as np

import physim as ps


def per_step(fn, steps, repeats=5):
    """Best-of-`repeats` wall time per step, in microseconds."""
    best = min(_timed(fn) for _ in range(repeats))
    return best / steps * 1e6


def _timed(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


def cluster(n, custom=False):
    rng = np.random.default_rng(42)
    w = ps.World(integrator="verlet")
    for r in rng.random((n, 3)):
        w.add_particle(r, mass=1.0 / n)
    if custom:
        w.add_force(ps.CustomForce(numpy_gravity, velocity_dependent=False))
    else:
        w.add_force(ps.NewtonianGravity(G=1.0, softening=0.01))
    return w


def numpy_gravity(t, pos, vel, mass, eps2=0.01**2):
    d = pos[None, :, :] - pos[:, None, :]
    inv_r3 = (np.einsum("ijk,ijk->ij", d, d) + eps2) ** -1.5
    np.fill_diagonal(inv_r3, 0.0)
    return np.einsum("ijk,ij,j->ik", d, inv_r3, mass)


def main():
    steps = 2_000
    rows = []

    w = cluster(10)
    rust_loop = per_step(lambda: w.step(1e-4, n=steps), steps)
    py_loop = per_step(lambda: [w.step(1e-4) for _ in range(steps)], steps)
    rows.append(("gravity N=10: Rust loop, w.step(dt, n)", rust_loop))
    rows.append(("gravity N=10: Python loop, w.step(dt)", py_loop))

    for n in (10, 100, 1000):
        s = steps if n < 1000 else 50
        builtin = cluster(n)
        custom = cluster(n, custom=True)
        rows.append((f"gravity N={n}: built-in NewtonianGravity", per_step(lambda: builtin.step(1e-4, n=s), s)))
        rows.append((f"gravity N={n}: CustomForce (numpy)", per_step(lambda: custom.step(1e-4, n=s), s)))

    width = max(len(name) for name, _ in rows)
    print(f"{'benchmark':<{width}}  µs/step")
    for name, us in rows:
        print(f"{name:<{width}}  {us:9.2f}")


if __name__ == "__main__":
    main()
