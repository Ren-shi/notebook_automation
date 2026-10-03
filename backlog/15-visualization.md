# 15 · Visualization and animation helpers

**Priority:** P3 · **Size:** S · **Area:** Usability · **Status: Done**

> **Done** (`python/physim/plot.py`, `tests/python/test_plot.py`, notebook §16).
> - `ps.plot.orbits`, `energy_error` (one run, a list, or a `{label: run}` dict; `time_unit` for e.g. days → years),
>   `phase_space` (velocity or momentum), `tidy` (spines, grid, frameless legend), `color(i)` / `COLORS` palette.
> - `ps.plot.animate`: a `FuncAnimation` with fading trails, frame skipping and a time label; saves to gif (pillow)
>   or mp4 (ffmpeg); tested by saving a gif.
> - `ps.plot.view3d`: interactive plotly figure with paths and a frame slider; `plotly` is an optional extra
>   (`physim[plot3d]`) and is in `requirements-dev.txt` so CI tests it.
> - matplotlib is imported on first use, so `import physim` does not need it.
> - The notebook now uses the helpers: every hand-written spine/legend styling line (10 + 9) and the
>   orbit and energy-error plotting loops are gone (21 helper calls); §16 shows an animation and phase portraits.

## Why
Every notebook repeats the same plotting code.

## Scope
- `physim.plot`: `orbits(traj)`, `energy_error(traj)`, `phase_space(traj, particle, axis)`.
- `physim.plot.animate(traj)`: matplotlib animation, saveable as mp4 or gif, with trails.
- Optional interactive 3D view for larger systems (e.g. plotly), behind an optional dependency.

## Done when
- The example notebook uses the helpers and loses most of its plotting boilerplate.
