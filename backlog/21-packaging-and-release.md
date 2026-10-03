# 21 · Packaging and release

**Priority:** P3 · **Size:** S · **Area:** Tooling

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
