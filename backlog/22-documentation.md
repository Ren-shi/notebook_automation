# 22 · Documentation

**Priority:** P3 · **Size:** M · **Area:** Usability · **Status: Done** (Pages publishing needs the owner to enable it, below)

> **Done** (`docs/`, `.github/workflows/docs.yml`, `notebooks/02-05`).
> - Sphinx site (furo theme, MyST Markdown, myst-nb): install, quickstart, user guide (included from the README so the
>   two never drift apart), theory notes, example notebooks executed at build time, Python API reference
>   (autodoc of every class and module; the 53 bindings members that had no docstring got one), Rust API reference
>   (`cargo doc`, linked in), changelog.
> - Theory notes: `theory/integrators.md` (every scheme: equations, order, cost, what it conserves, when to use or
>   avoid it; why symplectic matters; the adaptive controller; thermostats), `theory/forces.md` (the equations of every
>   built-in force, cutoffs, minimum image, neighbour lists, virial), `theory/methods.md` (events, RATTLE/M-SHAKE, hard
>   collisions, Lyapunov/MEGNO, the rigid-body splitting, deterministic parallelism, exact restarts).
> - `new-ideas.md`: prototype → validate (exact answer, conservation, gradient check, convergence order) → port to
>   Rust (trait, registration, Python wrapper, cross-check against the prototype) → benchmark. Its Python snippets
>   were run as written (precession to 9e-10, measured order 4.00).
> - New notebooks: `02_nbody_cluster` (Plummer sphere with tree gravity; Lagrangian radii match the analytic ones),
>   `03_three_body_chaos` (Pythagorean problem; 1e-9 perturbation grows to O(1) by t ≈ 60, while the figure-eight's
>   stays ~1e-7 and its MEGNO → 2), `04_driven_oscillator` (Duffing bifurcation diagram and strange attractor,
>   λ ≈ 0.1), `05_charged_particles` (magnetic bottle with Boris: mirror points match L cot α to 1e-4, μ conserved to
>   5e-5). All run in the notebook CI job (2-7 s each).
> - CI: `docs.yml` builds the site with warnings as errors (`sphinx-build -W`, `RUSTDOCFLAGS=-D warnings`) on every
>   push and pull request, and deploys it to GitHub Pages on version tags or a manual run.
> - The README's stale paths (`src/forces.rs`, `src/integrators.rs`) now point to the module directories.
>
> **Left for the repository owner:** enable GitHub Pages with source "GitHub Actions" (Settings → Pages) so the deploy
> job can publish.

## Why
The README and type stubs cover the basics. As the engine grows, users will need reference docs and the
reasoning behind the numerics.

## Scope
- Python API reference generated from docstrings (Sphinx or MkDocs) and Rust API docs (`cargo doc`).
- Theory notes per integrator and force: equations, order, conserved quantities, when to use or avoid it.
- More example notebooks: N-body cluster, chaotic three-body, driven damped oscillator, charged particles once
  item 09 lands.
- A "testing a new idea" guide expanding on the README: prototype in Python, validate, port to Rust, benchmark.

## Done when
- Docs build in CI and are published (e.g. GitHub Pages) on each release.
