# physim

A classical-mechanics engine for point particles: the time stepping and forces run in **Rust**, and
**Python** (via [PyO3](https://pyo3.rs) and [maturin](https://www.maturin.rs)) sets up systems and analyses results.

```python
import physim as ps

w = ps.World(integrator="yoshida4")
w.add_particle([0.5, 0, 0], vel=[0, 0.4, 0], mass=0.5)
w.add_particle([-0.5, 0, 0], vel=[0, -0.4, 0], mass=0.5)
w.add_force(ps.NewtonianGravity(G=1.0))

traj = w.run(dt=1e-3, steps=10_000, record_every=10)
traj.pos.shape                      # (frames, particles, 3)
ps.relative_energy_error(traj)      # (E - E0) / |E0| per frame
```

See [`notebooks/01_getting_started.ipynb`](notebooks/01_getting_started.ipynb) for integrator comparisons,
measured convergence orders, a Python-defined force checked against an analytic precession rate, and a
many-particle wave example.

## Setup

Requires a Rust toolchain and Python ≥ 3.9.

```bash
python -m venv .venv && source .venv/bin/activate
pip install maturin numpy matplotlib pytest jupyter
maturin develop --release        # builds the Rust extension into the venv; rerun after Rust changes
```

## Tests

```bash
cargo test --release             # Rust physics tests (convergence orders, conservation laws)
pytest tests/python              # Python binding tests
```

## What's included

**Integrators** (`World(integrator=...)`, switchable at any time via `w.integrator = "rk4"`):

| name | order | symplectic | notes |
|---|---|---|---|
| `explicit_euler` | 1 | no | baseline; energy drifts |
| `symplectic_euler` | 1 | yes | kick-drift |
| `verlet` | 2 | yes | velocity Verlet / KDK leapfrog, time-reversible |
| `yoshida4` | 4 | yes | Yoshida (1990) triple-jump composition of leapfrog |
| `rk4` | 4 | no | classical Runge-Kutta; consistent for velocity-dependent forces |

The symplectic schemes are only symplectic for velocity-independent forces; prefer `rk4` with drag.

**Forces**: `UniformField(g)`, `NewtonianGravity(G, softening)` (direct O(N²), Plummer softening),
`Spring(i, j, k, rest_length)`, `AnchorSpring(i, anchor, k, rest_length=0)`, `LinearDrag(gamma)` (`a = -γv`),
`QuadraticDrag(c)` (`F = -c|v|v`), and `CustomForce` (below). Conservative forces report a potential, so
`traj.energy` and `w.total_energy()` include them.

**Managing a world**:
- `add_force` returns an id. Use it with `remove_force(id)`, `replace_force(id, force)`, `force_params(id)` and
  `set_force_params(id, G=2.0)` to sweep a parameter without rebuilding the world. `w.forces` is `{id: name}`.
- `remove_particle(i)` shifts higher indices down and renumbers springs; it refuses while a force still uses particle `i`.
- `w.masses = [...]` changes masses. A mass of `0` makes a test particle: it feels gravity but does not source it.
  Springs and quadratic drag reject massless particles.
- `w.pin(i)` fixes a particle in place (it still exerts forces); `w.pin(i, False)` releases it.
- If a step fails (an exception in a `CustomForce`, or the state becoming non-finite), the world rolls back to
  the last completed step. From `run`, the exception carries the frames recorded so far as `err.trajectory`.
- `w.integrator_info` gives the integrator's name, order and whether it is symplectic.

**Diagnostics**: `kinetic_energy()`, `potential_energy()`, `total_energy()`, `momentum()`,
`angular_momentum()`, `center_of_mass()`, `accelerations()`; trajectories record `t`, `pos`, `vel`,
`kinetic`, `potential`, `energy`.

## Testing a new idea

### 1. Prototype in Python

```python
def accel(t, pos, vel, mass):          # pos, vel: (N, 3); mass: (N,) -> return (N, 3) accelerations
    r = np.linalg.norm(pos, axis=1, keepdims=True)
    return -pos / r**3 + 2 * beta * pos / r**4

def potential(t, pos, mass):           # optional; enables energy diagnostics
    r = np.linalg.norm(pos, axis=1)
    return float(np.sum(mass * (-1 / r + beta / r**2)))

w.add_force(ps.CustomForce(accel, potential=potential, velocity_dependent=False))
```

Exceptions raised inside your function propagate unchanged. Each force evaluation is a Python call,
so this is convenient but slow.

### 2. Port it to Rust for speed

Implement the `Force` trait in `src/forces.rs`:

```rust
pub struct MyForce { pub beta: f64 }

impl Force for MyForce {
    fn accumulate(&self, t: f64, pos: &[Vec3], vel: &[Vec3], mass: &[f64], acc: &mut [Vec3]) -> Result<()> {
        for (a, r) in acc.iter_mut().zip(pos) { /* *a += ... */ }
        Ok(())
    }
    fn potential(&self, t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> { Ok(Some(/* ... */)) }
    fn name(&self) -> String { "MyForce".into() }
}
```

Then expose it to Python in `src/python.rs`: add a `#[pyclass]` wrapper with an `impl PyForceSpec`
(copy `PyLinearDrag`), list it in `build_force` and the `#[pymodule]`, and export it from
`python/physim/__init__.py`. Run `maturin develop --release`.

New integrators work the same way: implement `Integrator` in `src/integrators.rs` and register a name in
`integrators::by_name`.

## Roadmap

Planned work lives in [`backlog/`](backlog/README.md), one file per item with priority, scope and acceptance criteria.

## Layout

```
src/
  vec3.rs         3-vector math
  state.rs        positions, velocities, masses, time; conserved quantities
  forces.rs       Force trait, ForceSet, built-in forces
  integrators.rs  Integrator trait and schemes
  world.rs        World (state + forces + integrator) and Trajectory recording
  python.rs       PyO3 bindings (feature "python")
python/physim/    Python package (re-exports the extension, analysis helpers, type stubs)
tests/            Rust tests; tests/python for the bindings
notebooks/        examples
backlog/          planned work
```
