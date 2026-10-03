# 06 · Save, load and stream simulation data

**Priority:** P1 · **Size:** M · **Area:** Core

## Why
`World::run` keeps every recorded frame in memory and there is no way to checkpoint or reproduce a run
later. Long or large runs need to write to disk as they go.

## Scope
- Save and load full state (time, particles, force descriptions, integrator) so a run can be restarted.
- Save trajectories to disk: `.npz` to start; HDF5 (`h5py` on the Python side) for large runs.
- Streaming: an observer callback or a writer that receives each recorded frame instead of
  accumulating it, so memory stays flat.
- Store metadata with every file: engine version, integrator, dt, force parameters.

## Notes
- `CustomForce` holds a Python function and cannot be serialised; save its name and require it to be
  supplied again on load.

## Done when
- Run, checkpoint, reload, continue gives the same final state as an uninterrupted run, bit for bit.
- A 10⁷-step run with recording every step uses constant memory when streaming.
