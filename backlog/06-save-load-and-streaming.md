# 06 · Save, load and stream simulation data

**Priority:** P1 · **Size:** M · **Area:** Core · **Status: Done**

> **Done.** Engine: `World::checkpoint`/`World::from_checkpoint` (`src/checkpoint.rs`), a `Recorder` trait and
> `World::run_into` for streaming, `RunOptions { energies }`. Python: `save_checkpoint`/`load_checkpoint` (`.npz`),
> `save_trajectory`/`load_trajectory` and a streaming `TrajectoryWriter` (`.npz` or `.h5`), `run(..., sink=,
> chunk_size=, energies=)`, and `Trajectory.metadata` (engine version, integrator, dt, masses, every force's parameters).
>
> Results:
> - Checkpoint → reload → continue matches an uninterrupted run **bit for bit** for all five integrators, with every
>   built-in force, a pinned particle, a removed force id and a re-supplied custom force, including gravity on the
>   parallel path (N = 500) (`tests/checkpoint.rs`, `tests/python/test_io.py`).
> - A 10⁷-step run recording every step uses flat memory when streamed: peak RSS 27 → 27 MiB with a counting sink,
>   30 MiB to `.npz`, 46 MiB to `.h5` (the same at 10⁶ and 10⁷ frames; a 720 MB file). Kept in memory, 10⁶ frames
>   already take 149 MiB.
> - Per-frame energies: for gravity with N = 1 000, recording all 100 steps takes 191 ms with energies and 118 ms
>   without (recording only the ends: 125 ms). For the spring chain the potential was not the main cost: copying
>   frames was. Pre-sizing the trajectory cut `record_every=1` from 84 ms to 64 ms (55 ms without energies;
>   30 ms recording only the ends).
> - Along the way: force-set versions are now globally unique, so the acceleration cache can never reuse a value
>   from a different `ForceSet` after `world.forces` is replaced.
>
> Not done: `CustomForce` functions are saved by name only (as planned); HDF5 checkpoints (checkpoints are small,
> `.npz` is enough).

## Why
`World::run` keeps every recorded frame in memory and there is no way to checkpoint or reproduce a run
later. Long or large runs need to write to disk as they go.

## Scope
- Save and load full state (time, particles, force descriptions, integrator) so a run can be restarted.
- Save trajectories to disk: `.npz` to start; HDF5 (`h5py` on the Python side) for large runs.
- Streaming: an observer callback or a writer that receives each recorded frame instead of
  accumulating it, so memory stays flat.
- Store metadata with every file: engine version, integrator, dt, force parameters.

- Make per-frame energies optional: recording every step currently more than doubles run time (item 02
  baseline), mostly from evaluating the potential, which is O(N²) for gravity.

## Notes
- `CustomForce` holds a Python function and cannot be serialised; save its name and require it to be
  supplied again on load.

## Done when
- Run, checkpoint, reload, continue gives the same final state as an uninterrupted run, bit for bit.
- A 10⁷-step run with recording every step uses constant memory when streaming.
