# 19 · More built-in forces

**Priority:** P2 · **Size:** M · **Area:** Physics · **Status: Done**

> **Done** (`src/forces/`: `springs.rs`, `central.rs`, `drives.rs`, `orbital.rs`, `closure.rs`; all exposed to Python
> except the Rust-only `ClosureForce`). New forces: `DampedSpring`, `ModulatedSpring` (to a particle or a point),
> `SpringNetwork`, `PowerLaw` (with `n = 0` as the logarithmic potential), `Yukawa`, `PlummerPotential`,
> `HernquistPotential`, `HarmonicTrap` (anisotropic), `PeriodicForce`, `PostNewtonian` (1PN, test-particle limit),
> `J2Oblateness` (with the reaction on the central body, so momentum is conserved), and
> `ClosureForce::new`/`per_particle`. All built-ins are saved in checkpoints (restarts stay bit-exact), described in
> run metadata, and their scalar/vector parameters work with `set_force_params`.
>
> Results (`tests/forces.rs`, `tests/python/test_more_forces.py`, notebook §8):
> - Driven damped oscillator: steady amplitude and phase match the resonance curve to ~2e-10 at seven drive
>   frequencies from 0.3 ω0 to 3 ω0.
> - 1PN periapsis advance: 6πGM/(c²a(1−e²)) to 3.0e-4 at c = 100 (e = 0.5). The residual is physical: it falls to 7.5e-5
>   and 1.9e-5 at c = 200 and 400, the O(GM/(c²a)) relative correction that 1PN leaves out.
> - J2 nodal regression: −(3/2) n J2 (R/a)² cos i to 1.0e-3 (J2 = 1e-3, so this is the expected O(J2) residual).
> - Mathieu equation (`ModulatedSpring`): Floquet growth rate at the principal resonance matches ε ω0/4 to 4.7e-5.
> - Every conservative force equals −∇U by finite differences; circular orbits in the Plummer, Hernquist,
>   logarithmic and Yukawa potentials stay at constant radius to 1e-9; `PowerLaw(n=−1, k=−GM)` reproduces Kepler.
> - `SpringNetwork` gives bit-identical trajectories to one `Spring` per bond. Per evaluation on a 10 000-bond chain it
>   costs ~98 µs against ~172 µs (≈10 ns against 17 ns per bond; folding the index/mass checks into the main loop
>   saved 20%). Per Verlet step the gain is smaller (235 µs against 351 µs) because the step itself has ~200 µs of
>   per-particle overhead at N = 10 000, mostly the state backup used to roll back failed steps (see the known
>   limits in the backlog README).
>
> Not done: a parallel `SpringNetwork` (bonds scatter into shared particles, so it needs coloured or per-block
> buffers), pairwise Yukawa (belongs with generic pair potentials in item 12), and 1PN back-reaction on the central body.

## Why
Several common forces exist only as Python `CustomForce` code, which is slow. They belong in Rust.

## Scope
- Damped spring (spring + dashpot along the bond).
- Central potentials: power law, Yukawa, Plummer and Hernquist halos, isotropic and anisotropic harmonic traps.
- Time-dependent drives: sinusoidal forcing, parametric modulation of a spring constant.
- Orbital perturbations: 1PN (relativistic) correction, J2 oblateness of the central body.
- A Rust-side "closure force" that takes a function or closure, so new forces in Rust need no boilerplate struct.
- `SpringNetwork`: all bonds of a lattice or polymer in one force, stored as index and parameter arrays.
  Separate `Spring` objects cost ~21 ns per spring per evaluation (item 02 baseline) from dynamic dispatch and
  per-spring checks.
- Expose each one to Python.

## Done when
- Each force has an analytic check: a driven damped oscillator's resonance curve; 1PN perihelion advance of
  6πGM/(c²a(1−e²)) per orbit; J2 nodal regression rate.
