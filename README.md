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

## Install

Prebuilt wheels (Linux x86_64/aarch64, macOS x86_64/arm64, Windows x64; one abi3 wheel per platform covers
Python ≥ 3.9) need no Rust toolchain. The distribution is called `physim-engine` because `physim` is taken on PyPI;
the import name is `physim`:

```bash
pip install physim-engine                 # once a release has been published
pip install "physim-engine[plot]"         # with matplotlib for physim.plot ([plot3d] adds plotly, [io] h5py)
```

Releases: bump `version` in `Cargo.toml` (shared by the crate, the wheel and `physim.__version__`), note the
changes in [`CHANGELOG.md`](CHANGELOG.md), merge, and push a tag `vX.Y.Z`. The
[release workflow](.github/workflows/release.yml) builds and tests the wheels and the sdist and publishes to PyPI by
trusted publishing (configure the PyPI project's trusted publisher for this repository and workflow, environment
`pypi`). Publishing the Rust crate is opt-in (repository variable `PUBLISH_CRATE=true` and secret
`CARGO_REGISTRY_TOKEN`). Licensed under the MIT license ([`LICENSE`](LICENSE)).

## Setup (from source)

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
| `dopri5` | 5 | no | Dormand-Prince 5(4) at a fixed step; 6 force evaluations per step; also the adaptive method below |
| `yoshida6`, `yoshida8` | 6, 8 | yes | Yoshida compositions of 7 and 15 Verlet substeps; support rods |
| `pefrl` | 4 | yes | Omelyan et al. (2002) optimised splitting: 4 evaluations, error 60× below `yoshida4` at equal cost |
| `blanes_moan4` | 4 | yes | Blanes & Moan (2002) 6-stage splitting: 6 evaluations, error 190× below `yoshida4` at equal cost |
| `gauss2`, `gauss4`, `gauss6` | 2, 4, 6 | yes | Gauss-Legendre implicit RK (`gauss2` = implicit midpoint): symplectic and symmetric even with velocity-dependent forces, conserve quadratic invariants exactly; fixed-point iteration, so several evaluations per stage |
| `wisdom_holman` | 2 | yes | Wisdom-Holman for planetary systems (particle 0 is the star; needs one `NewtonianGravity`, no pins): error ∝ planet/star mass ratio |
| `boris` | 2 | no (volume-preserving) | Boris pusher for charged particles: exact rotation about B, so `|v|` is conserved to round-off in a pure magnetic field; non-magnetic forces must not depend on velocity |

The explicit symplectic schemes are only symplectic for velocity-independent forces; with drag prefer `rk4` or
`dopri5`, and for Hamiltonian velocity-dependent forces (e.g. magnetic) the `gauss*` schemes. Aliases:
`forest_ruth` = `yoshida4`, `omelyan` = `pefrl`, `implicit_midpoint` = `gauss2`.

Your own schemes: `w.use_composition(weights, order)` composes Verlet substeps of `weights[k]·dt`, and
`w.use_splitting([("kick", b1), ("drift", a1), ...], order)` applies kicks and drifts in order. Both are saved in
checkpoints. On the outer solar system (Sun to Pluto) at 1/20 of Jupiter's period, `wisdom_holman` holds the energy
error at 4e-6 for a million years, against 2e-3 for `verlet` (`cargo run --release --example outer_solar_system`).

\* `verlet` and `yoshida4` reuse the end-of-step acceleration when no force depends on velocity (otherwise 2 and 6
evaluations). Forces must therefore be pure functions of `t`, positions, velocities and masses: a `CustomForce`
that reads mutable outside state should be declared `velocity_dependent=True` (the default) to disable reuse.

**Adaptive time stepping**: `w.run_adaptive(t_end, rtol=1e-9, atol=1e-12)` integrates to `t_end` with
Dormand-Prince 5(4), choosing each step so the estimated local error of every position and velocity component
stays below `atol + rtol·|value|`. It records every accepted step, or exactly the times you ask for with
`times=[...]` (interpolated by the method's continuous extension, at no extra cost). `events`, `energies` and
`sink` work as in `run`, events are located on the continuous extension, and `t_end` may lie in the past.
`traj.metadata["adaptive"]` reports accepted and rejected steps and force evaluations. Constraints are not supported.

```python
traj = w.run_adaptive(2 * math.pi, rtol=1e-10, atol=1e-13)                    # one orbit, every step
traj = w.run_adaptive(100.0, times=np.linspace(w.t, 100.0, 1001), events=[periapsis])
```

When to use it: on an e = 0.99 Kepler orbit, adaptive stepping reaches a given energy error with about 350 to 16 000 times
fewer force evaluations than any fixed-step integrator here (`cargo run --release --example eccentric_orbit`).
The catch is that it is not symplectic: its energy error grows linearly with time (about 1.6e-10 per orbit at
e = 0.5 and `rtol=1e-10`, so 1.6e-6 after 10 000 orbits), whereas `yoshida4` at the same cost stays bounded
(1e-12). Use adaptive steps for transients, close encounters and moderately long runs with widely varying time
scales; use a fixed-step symplectic integrator for very long integrations of smooth conservative motion.

**Molecular dynamics**: pair potentials with a cutoff — `LennardJones(epsilon, sigma, cutoff=2.5, shift=True,
box=L)`, `Morse(depth, a, r0, cutoff)`, and `TabulatedPair` (any `V(r)`, e.g. from `ps.tabulate_pair(f, r_min, r_max)`,
evaluated in Rust) — in an optional periodic box (minimum image; coordinates are not wrapped, so diffusion is easy to
measure). Neighbours come from cell lists and Verlet lists with a skin; results are identical whether a list is fresh
or reused, so restarts stay bit-exact. Thermostats: `w.use_langevin(T, friction, seed)` (BAOAB, reproducible
counter-based noise) and `w.use_nose_hoover(T, tau)` (`total_energy() + thermostat_energy()` is conserved).
Observables: `w.temperature()`, `w.pressure(volume)` (virial), `ps.radial_distribution(positions, box)`.

**Plotting**: `ps.plot` (matplotlib, imported on first use) — `orbits(traj, labels=...)`, `energy_error({"verlet":
t1, "rk4": t2})`, `phase_space(traj, particle, "x")`, `tidy(ax)` for the house style, and `animate(traj, trail=30)`
(save with `anim.save("orbit.gif", writer="pillow")` or mp4; in a notebook `HTML(anim.to_jshtml())`). `view3d(traj)`
gives an interactive plotly figure (`pip install physim[plot3d]`). They accept `Trajectory` and `RigidTrajectory`.

**Scenarios and units**: `ps.scenarios` builds standard systems as ready-to-run worlds — `two_body(m1, m2, a, e, i,
Omega, omega, nu=... | M=...)` (barycentric), `solar_system(["Jupiter", ...])` (JPL J2000 mean elements, Sun first),
`plummer_sphere(n)` (equilibrium cluster), `figure_eight()` (returns to its start to 3e-8 after one period),
`pendulum(length, theta0)` and `spring_lattice(shape, boundary="fixed")` — plus `elements_to_state` /
`state_to_elements`. The engine has no units; `ps.units` converts between `SI`, `ASTRONOMICAL` (AU, solar mass,
Gaussian year, so `G = 4π²`), `DIMENSIONLESS` and your own `units.scaled(name, length, mass, time)`:
`ps.units.ASTRONOMICAL.to_si(v, "velocity")`, `ps.units.convert(x, "energy", src, dst)`.

**Rigid bodies**: `ps.RigidSystem(order=2 | 4)` holds bodies with a mass, principal moments of inertia, an
orientation (unit quaternion) and angular momentum — free (`position` is the centre of mass) or turning about a fixed
`pivot` with a `center_of_mass` offset. Forces act at body points and give torques: `BodyGravity(g)`,
`BodySpring((body, point), (body, point) | anchor, k, rest_length)`, or any `f(t, positions, orientations) ->
(forces, torques)`. The step splits the free rotation into exact turns about the body axes (symplectic, time
reversible, angular momentum of a free body exact). A box spun about its intermediate axis flips with the
elliptic-integral period to 1e-4, and a fast heavy top precesses within 0.3% of `mgl / (I3 ω3)`
(`cargo run --release --example spinning_tops`). Inertia helpers: `ps.inertia_box`, `inertia_cylinder`,
`inertia_ellipsoid`, `inertia_sphere`; `ps.quaternion_from_axis_angle`.

**Collisions and contact**: particles get a size with `add_particle(..., radius=r)` or `w.radii = [...]`.
- Hard collisions: `w.set_collisions(restitution=1.0, walls=ps.box_walls(lo, hi))` makes spheres bounce off each
  other and off planar walls (`(normal, offset)` pairs). Each contact is found inside the step (swept spheres on a uniform
  grid), the world is moved to that instant, the impulse is applied and the step continues — exact event-driven
  dynamics when no forces act between collisions, and contact times refined on the integrator's trajectory when they do.
  A 500-sphere gas conserves energy to 3e-15 and relaxes to Maxwell-Boltzmann (⟨v⁴⟩/⟨v²⟩² = 1.669 vs 5/3).
- Soft contact: `SoftContact(k, damping=0, law="linear" | "hertz")`, an ordinary force (contact duration matches
  π√(μ/k) and the Hertz law to 6e-5).

**Chaos indicators**: `w.lyapunov(dt, steps, n=1)` integrates the variational equations with the world's own
integrator (tangent vectors ride along as extra pseudo-particles whose acceleration is the Jacobian-vector product
`∂a/∂x·δx + ∂a/∂v·δv`: analytic for most built-in forces, finite differences otherwise) and returns the `n` largest
Lyapunov exponents (Benettin's method, QR re-orthonormalisation), their running estimates, and MEGNO `⟨Y⟩` (→ 2 for
regular motion, grows like `λt/2` for chaos). `ps.poincare_section(w, dt, steps, event)` returns the states at every
crossing of a section, located to ~1e-12 of a step. On Hénon-Heiles (`ps.HenonHeiles()`) at E = 1/6 the chaotic sea
gives λ ≈ 0.1 while regular orbits give λ → 0 and ⟨Y⟩ = 2.00; the full spectrum comes in ± pairs summing to 1e-16.

```python
out = w.lyapunov(0.01, 200_000, n=1)          # out["exponents"], out["running"], out["mean_megno"]
y, py = ps.poincare_section(w, 0.01, 200_000, ps.Event.coordinate(0, "x", direction=1),
                            coordinates=lambda p, v: (p[:, 0, 1], v[:, 0, 1]))
```

**Forces**:
- Basics: `UniformField(g)`, `NewtonianGravity(G, softening)` (direct O(N²), Plummer softening),
  `TreeGravity(G, softening, theta=0.5, quadrupole=False)` (Barnes-Hut, O(N log N): rms force error 9e-4 at θ = 0.5 and
  20× faster than direct summation at N = 10⁵, 69× at θ = 1; `cargo run --release --example tree_gravity`),
  `LinearDrag(gamma)` (`a = -γv`), `QuadraticDrag(c)` (`F = -c|v|v`).
- Springs: `Spring(i, j, k, rest_length)`, `AnchorSpring(i, anchor, k, rest_length=0)`,
  `DampedSpring(i, j, k, rest_length, c)` (dashpot along the bond), `ModulatedSpring(i, to, k, depth, omega,
  phase=0, rest_length=0)` (`k(t) = k(1 + depth·cos(ωt + φ))`, to a particle or a point: parametric driving), and
  `SpringNetwork(i, j, k, rest_length)` (all bonds of a lattice or polymer in one force, from index arrays; ~1.75×
  faster per evaluation than separate springs, identical results).
- External central potentials about `center` (per unit mass, acting on every particle): `PowerLaw(k, n)`
  (`Φ = k rⁿ`, `k ln r` for `n = 0`), `Yukawa(k, length)` (`-k e^{-r/λ}/r`), `PlummerPotential(GM, a)`,
  `HernquistPotential(GM, a)`, `HarmonicTrap(omega)` (scalar or per-axis frequencies).
- `HenonHeiles(lam=1, center)`: the Hénon-Heiles potential in the xy-plane, the standard chaos test bed.
- Electromagnetism (charges set with `add_particle(..., charge=q)` or `w.charges = [...]`): `ElectricField(E)`,
  `MagneticField(B)` (`a = (q/m)(E + v × B)`), `Coulomb(k=1, softening=0)` (pairwise, same parallel pair loop as
  gravity), and `FieldForce(E=f, B=g)` for fields given as Python functions `f(t, pos) -> (N, 3)`.
- Drives: `PeriodicForce(i, amplitude, omega, phase=0)` (`F = A cos(ωt + φ)` on particle `i`).
- Orbital perturbations from a central particle, used alongside `NewtonianGravity`: `PostNewtonian(central, c)`
  (first post-Newtonian GR correction, test-particle limit; gives the `6πGM/(c²a(1−e²))` periapsis advance) and
  `J2Oblateness(central, J2, radius, axis)` (oblate central body; momentum-conserving).
- `CustomForce` (below) for Python prototypes; in Rust, `ClosureForce::per_particle(name, |t, r, v, m| ...)` or
  `ClosureForce::new(...)` turns a closure into a force without writing a struct; `FieldFunctions` does the same for
  E and B fields.

Conservative forces report a potential, so `traj.energy` and `w.total_energy()` include them. Every built-in force
is saved in checkpoints and run metadata, and its scalar and vector parameters can be changed with
`set_force_params`.

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
  adaptive.rs     adaptive Dormand-Prince run loop (error control, dense output)
  chaos.rs        variational equations, Lyapunov exponents, MEGNO
  collisions.rs   event-driven hard collisions and walls (broadphase.rs: uniform grid)
  events.rs       event functions and crossing detection
  constraints.rs  rigid rods and the RATTLE projections
  rigid.rs        rigid bodies: quaternions, torques, splitting step (RigidSystem)
  checkpoint.rs   save/restore a World
  parallel.rs     deterministic block-parallel helpers
  python.rs       PyO3 bindings (feature "python"; python/rigid.rs for RigidSystem)
scripts/          developer tools (notebook runner)
python/physim/    Python package (re-exports the extension, file I/O, analysis helpers, scenarios, units, plotting, type stubs)
tests/            Rust tests; tests/python for the bindings
benches/          Rust (criterion) and Python benchmarks
examples/         standalone Rust programs (cargo run --release --example <name>)
notebooks/        examples
backlog/          planned work
```
