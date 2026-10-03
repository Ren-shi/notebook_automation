# physim

[![CI](https://github.com/Ren-shi/notebook_automation/actions/workflows/ci.yml/badge.svg)](https://github.com/Ren-shi/notebook_automation/actions/workflows/ci.yml)

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
pip install -r requirements-dev.txt jupyter
maturin develop --release        # builds the Rust extension into the venv; rerun after Rust changes
```

## Tests

```bash
cargo test --release             # Rust physics tests (convergence orders, conservation laws)
pytest tests/python              # Python binding tests
python scripts/run_notebooks.py  # execute the example notebooks, fail on any error
```

## Benchmarks

```bash
cargo bench --bench engine           # Rust engine (criterion); reports in target/criterion/
python benches/python_overhead.py    # Python call overhead and CustomForce cost
```

Baseline numbers and what they imply are in [`backlog/02-benchmarks.md`](backlog/02-benchmarks.md).

CI (`.github/workflows/ci.yml`) runs all of the above plus `cargo fmt --check` and
`cargo clippy --all-targets --features python -- -D warnings` on every push and pull request.

## What's included

**Integrators** (`World(integrator=...)`, switchable at any time via `w.integrator = "rk4"`):

| name | order | symplectic | notes |
|---|---|---|---|
| `explicit_euler` | 1 | no | baseline; energy drifts |
| `symplectic_euler` | 1 | yes | kick-drift |
| `verlet` | 2 | yes | velocity Verlet / KDK leapfrog, time-reversible; 1 force evaluation per step*; RATTLE with rods |
| `yoshida4` | 4 | yes | Yoshida (1990) triple-jump composition of leapfrog; 3 force evaluations per step*; supports rods |
| `rk4` | 4 | no | classical Runge-Kutta; consistent for velocity-dependent forces |

The symplectic schemes are only symplectic for velocity-independent forces; prefer `rk4` with drag.

\* `verlet` and `yoshida4` reuse the end-of-step acceleration when no force depends on velocity (otherwise 2 and 6
evaluations). Forces must therefore be pure functions of `t`, positions, velocities and masses: a `CustomForce`
that reads mutable outside state should be declared `velocity_dependent=True` (the default) to disable reuse.

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

**Events**: find the exact moment something happens, instead of scanning recorded frames.
```python
periapsis = ps.Event.radial_velocity(0, 1, direction=+1)       # (r0 - r1)·(v0 - v1) rising through 0
ground = ps.Event.coordinate(0, "y", direction=-1, terminal=True)  # stop when particle 0 reaches y = 0
custom = ps.Event(lambda t, pos, vel, mass: pos[0, 0] - 2.0)     # any g(t, state); zeros are events
traj = w.run(dt, steps, events=[periapsis, ground, custom])
traj.event_t, traj.event_index, traj.event_pos, traj.terminated_by
```
After each step, a sign change of `g` is located by re-stepping the integrator from the start of the step with a
fraction of `dt` (Illinois root finding, to ~1e-12 of a step), so event states are as accurate as the integrator
itself. `direction` is +1 (rising), -1 (falling) or 0 (either); a terminal event stops the run and leaves the world
at the event. Built-in events (`radial_velocity`, `coordinate`, `separation`) run in Rust at no per-step Python cost.
Two crossings of one event within a single step are not detected, so keep `dt` small compared with event spacing.

**Rigid constraints**: rods instead of stiff springs.
```python
bob = w.add_particle([1, 0, 0])
w.add_rod(bob, [0, 0, 0])            # to a fixed point; length = current distance
w.add_rod(2, bob, length=0.5)        # between particles (length must match their distance)
traj = w.run(dt, steps)
traj.tension                          # (frames, rods): force pulling each rod's ends together
w.constraint_tensions(), w.constraints, w.remove_constraint(id), w.constraint_tolerance
```
- Integrated with RATTLE (velocity Verlet plus Lagrange multipliers): `verlet` is second order and `yoshida4`
  composes it to fourth order; both stay symplectic. Other integrators refuse to step a constrained world.
- Lengths hold to `constraint_tolerance` (default 1e-10, relative) at every step; velocity components along rods are
  removed. All rods are solved together (conjugate gradients on `G M⁻¹ Gᵀ`), so chains and closed loops converge:
  a 2-rod double pendulum costs ~1.7 µs/step (Verlet; 5 µs with yoshida4); chains of 10/100/300 links ~6 µs/0.36 ms/3.3 ms per step (cost grows
  as links²).
- The rigid pendulum's period matches the elliptic-integral result to 1e-10 from small angles up to 3 rad;
  tensions agree with m(g + v²/L) to O(dt²). Rods to pinned particles behave like rods to fixed points; massless
  particles cannot be constrained. Constraints are saved in checkpoints and trajectory metadata.

**Saving, loading and streaming**:
```python
ps.save_checkpoint(w, "run.npz")                     # state, forces (with ids), integrator
w2 = ps.load_checkpoint("run.npz", custom_forces={"precession": precession})
# continuing w2 gives exactly the same numbers as continuing w

ps.save_trajectory(traj, "orbit.h5")                 # or .npz; includes events and traj.metadata
traj = ps.load_trajectory("orbit.h5")

with ps.TrajectoryWriter("long.h5") as out:          # stream to disk: memory stays flat
    w.run(dt, 10_000_000, sink=out, energies=False)
```
- Restarting from a checkpoint is bit-for-bit identical to never having stopped. `CustomForce` functions cannot be
  saved: they are stored by name and must be passed again via `custom_forces` (by name or force id).
- `traj.metadata` records the engine version, integrator, `dt`, `record_every`, masses and every force's
  parameters; it is stored with saved trajectories as JSON.
- `run(..., sink=f, chunk_size=1024)` passes recorded frames and events to `f` in chunks (as `Trajectory`
  objects) instead of keeping them; any callable works. A 10⁷-frame run streams to `.h5` or `.npz` at a constant
  ~30–50 MB of memory.
- `run(..., energies=False)` skips the per-frame kinetic and potential energy. The potential is a full force pass
  (O(N²) for gravity), so this makes recording every step nearly free for N-body runs.
- `.npz` needs nothing extra; `.h5` needs `h5py` and can be read in slices (`h5py.File(p)["pos"][a:b]`) when a
  file is larger than memory. In Rust, `World::checkpoint`/`from_checkpoint` and `World::run_into` with a
  `Recorder` are the same features.

**Parallelism**: direct-sum gravity and per-particle forces use all CPU cores (rayon, `parallel` feature, on by
default). Work is split into blocks that depend only on N, so results are bit-identical for any thread count and
with the feature off. Set `RAYON_NUM_THREADS` to limit threads. Small systems (under ~360 bodies for gravity) run
on one thread, avoiding overhead.

**Diagnostics**: `kinetic_energy()`, `potential_energy()`, `total_energy()`, `momentum()`,
`angular_momentum()`, `center_of_mass()`, `accelerations()`; trajectories record `t`, `pos`, `vel`,
`kinetic`, `potential`, `energy` (the last three are `None` with `energies=False`).

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
  world.rs        World (state + forces + integrator), run loop, Recorder, Trajectory
  events.rs       event functions and crossing detection
  constraints.rs  rigid rods and the RATTLE projections
  checkpoint.rs   save/restore a World
  parallel.rs     deterministic block-parallel helpers
  python.rs       PyO3 bindings (feature "python")
scripts/          developer tools (notebook runner)
python/physim/    Python package (re-exports the extension, file I/O, analysis helpers, type stubs)
tests/            Rust tests; tests/python for the bindings
benches/          Rust (criterion) and Python benchmarks
notebooks/        examples
backlog/          planned work
```
