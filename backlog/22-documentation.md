# 22 · Documentation

**Priority:** P3 · **Size:** M · **Area:** Usability

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
