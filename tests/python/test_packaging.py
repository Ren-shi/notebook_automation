import re
from importlib import metadata
from pathlib import Path

import pytest

import physim as ps


def test_version_is_shared():
    # The distribution is physim-engine (physim is taken on PyPI); the import name is physim.
    assert metadata.version("physim-engine") == ps.__version__
    assert re.fullmatch(r"\d+\.\d+\.\d+([.-]?\w+)?", ps.__version__)
    cargo = Path(__file__).resolve().parents[2] / "Cargo.toml"
    if not cargo.exists():
        pytest.skip("not run from the source tree")
    version = re.search(r'^version = "([^"]+)"', cargo.read_text(), re.M).group(1)
    assert version == ps.__version__


def test_package_ships_type_information():
    root = Path(ps.__file__).parent
    assert (root / "py.typed").exists() and (root / "_core.pyi").exists()
