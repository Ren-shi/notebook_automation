# 01 · Continuous integration

**Priority:** P1 · **Size:** S · **Area:** Tooling · **Status: Done**

> **Done.** `.github/workflows/ci.yml` has two jobs on every push and pull request: **Rust** (fmt, clippy with
> `-D warnings`, `cargo test --release`) and **Python** (build with maturin, `pytest tests/python`, execute
> notebooks via `scripts/run_notebooks.py`). Dependencies are pinned in `requirements-dev.txt`; cargo and pip
> are cached. Multi-platform wheels are left to item 21.

## Why
Nothing checks the build automatically. Every later item changes numerics, and a silent regression in
an integrator or force is exactly the kind of bug that is hard to notice by eye.

## Scope
GitHub Actions workflow on push and pull request:
- `cargo fmt --check`, `cargo clippy --all-targets --features python -- -D warnings`
- `cargo test --release`
- Build the extension with `maturin develop --release` in a venv, then `pytest tests/python`
- Execute `notebooks/*.ipynb` (nbclient) so examples cannot rot

## Notes
- Cache `~/.cargo` and `target/`.
- One Linux job is enough to start; add macOS and Windows wheels later (`maturin-action`) if the package
  will be installed elsewhere.

## Done when
- A PR that breaks a physics test shows a red check.
- The workflow finishes in under ~5 minutes on a warm cache.
