# Testing a new idea

The workflow the engine is built for: **prototype** new physics in Python, **validate** it against something you
know, **port** it to Rust once it is right, and **benchmark** the result. The running example is a central force
with an extra $1/r^3$ term, $\mathbf a = -k\mathbf r/r^3 + 2\beta\mathbf r/r^4$ (potential $-k/r + \beta/r^2$ per
unit mass), whose orbits precess by exactly $2\pi(1/\gamma - 1)$ per radial period with
$\gamma = \sqrt{1 + 2\beta/L^2}$.

## 1. Prototype in Python

```python
import numpy as np
import physim as ps

k, beta = 1.0, 0.01

def accel(t, pos, vel, mass):          # pos, vel: (N, 3); mass: (N,) -> (N, 3) accelerations
    r = np.linalg.norm(pos, axis=1, keepdims=True)
    return -k * pos / r**3 + 2 * beta * pos / r**4

def potential(t, pos, mass):           # optional: enables energy diagnostics
    r = np.linalg.norm(pos, axis=1)
    return float(np.sum(mass * (-k / r + beta / r**2)))

w = ps.World(integrator="yoshida4")
w.add_particle([1.0, 0, 0], vel=[0, 1.2, 0])
w.add_force(ps.CustomForce(accel, potential=potential, velocity_dependent=False))
```

Declare `velocity_dependent=False` only when the function really ignores `vel`: it lets `verlet` and `yoshida4`
reuse accelerations and keeps the symplectic schemes exact. An exception inside the function propagates unchanged,
with the world rolled back to the last completed step and the frames recorded so far attached as `err.trajectory`.
Each force evaluation is a call into Python (tens of microseconds), which is fine for a prototype.

## 2. Validate

Check the idea against things that must hold, before trusting any new result.

**An exact answer.** Here, the precession rate. Find periapses with an event rather than by scanning frames:

```python
L = np.cross(w.positions[0], w.velocities[0])[2]
tr = w.run(2e-3, 60_000, record_every=10, events=[ps.Event.radial_velocity(0, direction=+1)])
phi = np.unwrap(np.arctan2(tr.event_pos[:, 0, 1], tr.event_pos[:, 0, 0]))
measured = np.mean(np.diff(phi))
predicted = 2 * np.pi * (1 / np.sqrt(1 + 2 * beta / L**2) - 1)
assert abs(measured / predicted - 1) < 1e-6
```

**Conservation laws.** A conservative force with a correct potential keeps the energy error bounded under a
symplectic integrator; a drift means the force and potential disagree (or the force depends on velocity):

```python
assert np.abs(ps.relative_energy_error(tr)).max() < 1e-8
```

**The force is minus the gradient of the potential.** Central differences catch sign and factor errors:

```python
pos, m = w.positions, w.masses
a = accel(0, pos, None, m)
h = 1e-6
for axis in range(3):
    e = np.zeros_like(pos); e[0, axis] = h
    grad = (potential(0, pos + e, m) - potential(0, pos - e, m)) / (2 * h)
    assert abs(m[0] * a[0, axis] + grad) < 1e-6
```

**The convergence order.** Halving the step must divide the error by $2^p$ for an order-$p$ integrator; a wrong
order usually means the force is not smooth or not a pure function of its inputs:

```python
def error(dt, t_end=2.0):
    def run(h):
        w = ps.World("yoshida4"); w.add_particle([1.0, 0, 0], vel=[0, 1.2, 0])
        w.add_force(ps.CustomForce(accel, velocity_dependent=False))
        w.step(h, n=round(t_end / h)); return w.positions[0]
    return np.linalg.norm(run(dt) - run(dt / 64))           # against a much finer reference

print(np.log2(error(0.02) / error(0.01)))                   # ≈ 4
```

## 3. Port it to Rust

Implement the `Force` trait. A new file `src/forces/my_force.rs`:

```rust
use crate::error::Result;
use crate::forces::Force;
use crate::vec3::Vec3;

/// `a = -k r/r³ + 2β r/r⁴` about the origin, per unit mass.
pub struct Precessing { pub k: f64, pub beta: f64 }

impl Force for Precessing {
    fn accumulate(&self, _t: f64, pos: &[Vec3], _vel: &[Vec3], _mass: &[f64], acc: &mut [Vec3]) -> Result<()> {
        for (a, r) in acc.iter_mut().zip(pos) {
            let d = r.norm();
            *a += *r * (-self.k / d.powi(3) + 2.0 * self.beta / d.powi(4));   // add, never overwrite
        }
        Ok(())
    }
    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        Ok(Some(pos.iter().zip(mass).map(|(r, m)| m * (-self.k / r.norm() + self.beta / r.norm_squared())).sum()))
    }
    fn name(&self) -> String { "Precessing".into() }
}
```

Register it with `mod my_force; pub use my_force::Precessing;` in `src/forces/mod.rs` (and re-export it from
`src/lib.rs` if Rust users should see it at the top level). Optional trait methods add more:
`velocity_dependent` (default `false`), `jacobian_vector` (analytic tangent dynamics for Lyapunov exponents,
otherwise finite differences), `accumulate_charged`/`accumulate_full` (forces that need charges or radii),
`virial` (pressure), and `builtin`/`params`/`set_param` (saving in checkpoints, `set_force_params`; this needs a
variant in `BuiltinForce`). For quick experiments in Rust without a struct,
`ClosureForce::per_particle("name", |t, r, v, m| ...)` wraps a closure.

**Expose it to Python** in `src/python.rs`: a `#[pyclass(frozen, name = "Precessing", module = "physim")]` wrapper
with a `#[new]` constructor and an `impl PyForceSpec` that builds the Rust force (copy `PyLinearDrag`), its name in
the `try_build!` list of `build_force`, `m.add_class::<PyPrecessing>()?` in the `#[pymodule]`, the export in
`python/physim/__init__.py`, and a stub in `python/physim/_core.pyi`. Rebuild with `maturin develop --release`.

**Check it against the prototype**: the two must agree to round-off on random configurations.

```python
pos = np.random.default_rng(0).normal(size=(50, 3))
w1, w2 = ps.World(), ps.World()
for w in (w1, w2):
    for p in pos: w.add_particle(p)
w1.add_force(ps.CustomForce(accel, velocity_dependent=False))
w2.add_force(ps.Precessing(k=k, beta=beta))
assert np.allclose(w1.accelerations(), w2.accelerations(), rtol=1e-12)
```

Then turn the validation of step 2 into tests: Rust tests in `tests/*.rs` (run in CI with and without the
`parallel` feature) and Python tests in `tests/python/`.

New integrators work the same way: implement the `Integrator` trait in `src/integrators/` and register a name in
`integrators::by_name` and `NAMES`. Only schemes built from kicks and drifts can simply be composed in Python with
`w.use_composition` or `w.use_splitting`, with no Rust at all.

## 4. Benchmark

```bash
cargo bench --bench engine           # criterion; add a group to benches/engine.rs for the new force
python benches/python_overhead.py    # per-step cost of Python loops and CustomForce vs built-in forces
```

In Python, time whole runs rather than single steps so the call overhead is amortised:

```python
import time
start = time.perf_counter(); w2.step(1e-3, n=10_000)
print(f"{(time.perf_counter() - start) / 1e4 * 1e6:.2f} µs per step")
```

For scale, `benches/python_overhead.py` measures gravity as a built-in Rust force 15–40× faster than the same force
as a vectorised numpy `CustomForce` (N = 10 to 1000); unvectorised Python is far slower still. For large N, check
that the force loop is parallel and deterministic: use `crate::parallel` helpers (blocks that depend only on N), and
confirm the result is identical with `RAYON_NUM_THREADS=1` and without the `parallel` feature.
