# 31 · Seeded, reproducible randomness

**Priority:** P2 · **Size:** S · **Area:** Core · **Status: Done**

> **Done** (`src/rng.rs` now public, `src/world.rs` `World::seed`, `src/checkpoint.rs`, `src/python.rs`,
> `python/physim/random.py`, `tests/rng.rs`, `tests/python/test_random.py`, `docs/theory/integrators.md`).
> - The counter-based generator already used by the Langevin thermostat (SplitMix64 mixing of
>   `(seed, counter, index)`) is now the public `physim::rng` with its key layout documented, rather than a new
>   Philox implementation: it has the same property that matters (every value a pure function of its key).
> - `World::seed` (default 1, Python `World.seed`), saved in and restored from checkpoints (old checkpoints get 1);
>   `use_langevin(seed=None)` takes the world's seed.
> - `physim.random`: `uniform`, `normal` (the same numbers as Rust), `positions_in_box`, `positions_in_sphere`,
>   `maxwell_boltzmann(world, T)` (pinned particles stay at rest, centre-of-mass momentum removed). Setup draws use
>   counters from 2⁶² up (`physim.random.SETUP`), so they never coincide with a run's (counter = step).
>
> Results: a Langevin run of 1500 particles with blocked, multi-threaded gravity is bit-identical with 1 and 4
> threads; runs restart bit for bit from checkpoints (Rust and Python, through `.npz` files); the Python draws are
> exactly the thermostat's (a BAOAB step checked to 1e-12); uniforms pass χ² (100 bins, 10⁶ values) and lag-1 and
> cross-seed correlation tests; Maxwell–Boltzmann variances are k_BT/m within 5% and the ball's radial distribution
> follows r³.
>
> Left out: a stochastic force usable with any integrator (white noise is only well defined inside a stochastic
> integrator such as BAOAB).

## Why
Thermostats, noise forces and random initial conditions all need random numbers. Runs must be reproducible,
including across threads and checkpoint restarts.

## Scope
- One counter-based RNG (e.g. Philox) owned by `World`, seeded explicitly, with independent streams per particle
  so parallel results do not depend on thread count.
- RNG state stored in checkpoints (item 06), so a restart reproduces the original run bit for bit.
- Python helpers for random initial conditions (Maxwell–Boltzmann velocities, random positions in a box or sphere)
  drawing from the same seed.

## Done when
- A noisy run gives bit-identical results with 1 thread and with all threads, and across a checkpoint restart.
