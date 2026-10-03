# 21 · Packaging and release

**Priority:** P3 · **Size:** S · **Area:** Tooling · **Status: Done** (publishing needs the owner's PyPI setup, below)

> **Done** (`.github/workflows/release.yml`, `LICENSE`, `CHANGELOG.md`, `pyproject.toml`, `Cargo.toml`,
> `tests/python/test_packaging.py`).
> - Distribution name **`physim-engine`**: `physim` is already taken on PyPI by an unrelated project. The import
>   name stays `physim`. Changing the name is one line in `pyproject.toml`.
> - Release workflow: abi3 wheels for Linux (x86_64, aarch64; manylinux), macOS (x86_64 cross-built and arm64) and
>   Windows (x64), plus an sdist; every wheel is installed on Linux/macOS/Windows with Python 3.9 and 3.13 without a
>   Rust toolchain and the binding tests run against it; the sdist is compiled from scratch and tested. On a tag
>   `vX.Y.Z` (checked against `Cargo.toml`) it publishes to PyPI by trusted publishing; the crate is opt-in.
>   Pull requests touching packaging files run the build and test jobs.
> - One version, in `Cargo.toml`: the wheel takes it through maturin (`dynamic = ["version"]`) and
>   `physim.__version__` from `CARGO_PKG_VERSION`; a test checks all three agree. Semantic versioning and the
>   release steps are described in `CHANGELOG.md`.
> - MIT `LICENSE` (matching `Cargo.toml`), shipped in the wheel and sdist; crate metadata for crates.io; optional
>   extras `plot`, `plot3d`, `io`.
>
> Verified locally: the wheel installs into a clean virtualenv with no `cargo` on `PATH` and passes the binding tests
> (88 passed, 4 skipped for h5py/plotly); the sdist builds a working wheel; `cargo publish --dry-run` passes. Building the
> sdist found that Cargo's `exclude = ["python/"]` also dropped `src/python/`; patterns are now anchored (`/python/`).
>
> **Left for the repository owner:** the first upload needs a PyPI account with a trusted publisher for
> `physim-engine` pointing at this repository, `release.yml` and environment `pypi`; then pushing tag `v0.1.0`
> publishes it and `pip install physim-engine` works on a clean machine.

## Why
Installing physim currently means cloning the repo and having a Rust toolchain. That is fine for one person;
anyone else needs prebuilt wheels.

## Scope
- Build wheels for Linux, macOS and Windows with `maturin-action` (the abi3 build means one wheel per platform
  covers Python 3.9+).
- Publish to PyPI on tagged releases; optionally publish the Rust crate to crates.io.
- Semantic versioning, a CHANGELOG, and a version string exposed as `physim.__version__`.
- Choose and add a LICENSE file (Cargo.toml currently says MIT, but there is no license file).

## Depends on
01.

## Done when
- `pip install physim` works on a clean machine without Rust installed.
