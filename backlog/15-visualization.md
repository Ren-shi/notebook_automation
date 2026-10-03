# 15 · Visualization and animation helpers

**Priority:** P3 · **Size:** S · **Area:** Usability

## Why
Every notebook repeats the same plotting code.

## Scope
- `physim.plot`: `orbits(traj)`, `energy_error(traj)`, `phase_space(traj, particle, axis)`.
- `physim.plot.animate(traj)`: matplotlib animation, saveable as mp4 or gif, with trails.
- Optional interactive 3D view for larger systems (e.g. plotly), behind an optional dependency.

## Done when
- The example notebook uses the helpers and loses most of its plotting boilerplate.
