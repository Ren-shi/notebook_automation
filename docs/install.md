# Installation

## From PyPI

Prebuilt wheels need no Rust toolchain. One abi3 wheel per platform covers Python 3.9 and newer, on Linux
(x86_64, aarch64), macOS (x86_64, arm64) and Windows (x64). The distribution is called `physim-engine` because
`physim` is taken on PyPI; the import name is `physim`.

```bash
pip install physim-engine                 # numpy is the only required dependency
pip install "physim-engine[plot]"         # + matplotlib for physim.plot
pip install "physim-engine[plot3d]"       # + plotly for physim.plot.view3d
pip install "physim-engine[io]"           # + h5py for .h5 trajectory files
```

## From source

Requires a Rust toolchain (<https://rustup.rs>) and Python ≥ 3.9.

```bash
git clone https://github.com/Ren-shi/notebook_automation && cd notebook_automation
python -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt
maturin develop --release        # builds the extension into the venv; rerun after Rust changes
```

Checks a contributor runs:

```bash
cargo fmt --check && cargo clippy --all-targets --features python -- -D warnings
cargo test --release             # Rust tests (convergence orders, conservation laws, restarts)
pytest tests/python              # Python binding tests
python scripts/run_notebooks.py  # execute every example notebook
```

## Building these docs

```bash
pip install -r docs/requirements.txt
sphinx-build -W -b html docs docs/_build/html
cargo doc --no-deps && cp -r target/doc docs/_build/html/rust    # Rust API reference
```

## Using the Rust crate directly

The engine is an ordinary Rust library (`physim`); the Python bindings are behind the `python` feature and
multi-threading behind `parallel` (on by default).

```toml
[dependencies]
physim = { git = "https://github.com/Ren-shi/notebook_automation" }
```
