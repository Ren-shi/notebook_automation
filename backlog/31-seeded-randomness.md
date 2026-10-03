# 31 · Seeded, reproducible randomness

**Priority:** P2 · **Size:** S · **Area:** Core

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
