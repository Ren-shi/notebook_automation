# 12 · Molecular dynamics: pair potentials, periodic boxes, thermostats

**Priority:** P2 · **Size:** L · **Area:** Physics · **Status: Done**

> **Done** (`src/forces/pair.rs`, `src/integrators/thermostat.rs`, `src/rng.rs`; Python `LennardJones`, `Morse`,
> `TabulatedPair` + `ps.tabulate_pair`, `w.use_langevin`, `w.use_nose_hoover`, `w.temperature()`, `w.pressure(V)`,
> `w.thermostat_energy()`, `ps.radial_distribution`).
> - `PairPotential { kind: LennardJones | Morse | Table, cutoff, shift, period }`: optional periodic box with the
>   minimum-image convention (coordinates are not wrapped). Cell lists build Verlet neighbour lists with a skin of
>   0.12 × cutoff, rebuilt when any particle moves half the skin. Each particle's list is sorted and pairs are re-filtered
>   by the cutoff every evaluation, so results never depend on when the list was built (restarts stay bit-exact). Small
>   systems (< 8192 particles) use each pair once (serial); larger ones gather per particle in parallel.
> - Tables: cubic Hermite interpolation of V with its derivative (continuous forces), linear continuation below the
>   table; a vectorised Python function is tabulated once, then evaluated in Rust.
> - Pressure from the virial `Σ r·F` (new `Force::virial`), temperature over free particles (k_B = 1).
> - Thermostats as integrators: Langevin by BAOAB with a counter-based normal generator (reproducible, independent of
>   threads; the counter is saved in checkpoints), and Nosé-Hoover by the Martyna et al. Trotter splitting with its
>   extended energy available. Both checkpoint exactly (new `Scheme::Langevin`/`NoseHoover`).
>
> Results (`tests/md.rs`, `tests/python/test_md.py`, notebook §13, `cargo run --release --example lj_liquid`):
> - NVE Lennard-Jones liquid, N = 500 at ρ = 0.8442, T = 0.721, rc = 2.5σ shifted, dt = 0.005: over 10⁵ steps the
>   total energy drifts by −4.8e-5 relative (maximum deviation 1.1e-4; the residual comes from the force jump at the
>   cutoff). Momentum is conserved to 1e-10.
> - g(r) at that state point: first peak 2.99 at r = 1.09σ, minimum 0.57 at 1.57σ, second peak 1.29 at 2.05σ — the
>   classic triple-point structure (Verlet 1968).
> - Langevin, target T = 1 from a cold start: ⟨T⟩ = 1.0015 with relative fluctuation 0.048 (canonical: √(2/3N) =
>   0.051). Nosé-Hoover: ⟨T⟩ = 0.998 with extended energy conserved to 5e-5, but temperature fluctuations of 0.10 —
>   a single Nosé-Hoover thermostat rings slowly after a large temperature jump; Nosé-Hoover chains would fix this.
> - Virial pressure equals the kinetic term minus dU/dV from box scaling to 1e-5; analytic forces equal −∇U for
>   Lennard-Jones and Morse across the periodic boundary; a table built from LJ reproduces it to 1e-6.
> - Cost: ~0.25 ms per Verlet step for 500 atoms between list rebuilds; rebuilds (every ~10 steps in this hot
>   liquid) cost ~2.5 ms because a box this small has only three cells per side.
>
> Not done: Nosé-Hoover chains, mixtures (per-species parameters and mixing rules), long-range corrections and Ewald
> sums, barostats (NPT), and wrapping coordinates for output.

## Why
Statistical-mechanics experiments (phase behaviour, transport, equilibration) need short-range pair
potentials in a periodic box and temperature control.

## Scope
- Generic pair potential with cutoff: Lennard-Jones, Morse, and a user-defined `V(r)` given as a table
  or a vectorised Python function evaluated once per step (avoids per-pair Python calls).
- Periodic boundary conditions with minimum-image convention; cell lists or neighbour lists for O(N).
- Thermostats: Langevin (needs a seeded RNG and a stochastic integrator such as BAOAB), Nosé–Hoover.
- Observables: temperature, pressure (virial), radial distribution function g(r).

## Done when
- NVE Lennard-Jones liquid conserves energy over 10⁵ steps with a shifted potential.
- Langevin and Nosé–Hoover runs reach the target temperature; g(r) matches published LJ data at a
  reference state point.
