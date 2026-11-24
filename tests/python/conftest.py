"""``--quick``: leave out the tests of .github/slow-tests.txt, as pull requests do (plain pytest runs everything)."""

from pathlib import Path

import pytest

SLOW_TESTS = Path(__file__).resolve().parents[2] / ".github" / "slow-tests.txt"


def slow_entries() -> set:
    """The lines of the list: ``file::test``, ``file::test[parameters]`` or ``file::fixture:NAME``."""
    if not SLOW_TESTS.is_file():  # an unpacked source distribution
        return set()
    lines = (x.strip() for x in SLOW_TESTS.read_text(encoding="utf-8").splitlines())
    return {x for x in lines if x and not x.startswith("#")}


def is_slow(file: str, name: str, fixtures, entries: set) -> bool:
    """Whether the test ``name`` (with its parameters) of ``file`` is in the list, itself or by a fixture it takes."""
    return (f"{file}::{name}" in entries or f"{file}::{name.partition('[')[0]}" in entries
            or any(f"{file}::fixture:{x}" in entries for x in fixtures))


def pytest_addoption(parser):
    parser.addoption("--quick", action="store_true", help="leave out the tests of .github/slow-tests.txt")


# Last, so that the whole collection is judged by its final order and fixtures.
@pytest.hookimpl(trylast=True)
def pytest_collection_modifyitems(config, items):
    if not config.getoption("--quick"):
        return
    entries = slow_entries()
    out = [x for x in items if is_slow(x.path.name, x.name, getattr(x, "fixturenames", ()), entries)]
    if out:
        config.hook.pytest_deselected(items=out)
        gone = set(map(id, out))
        items[:] = [x for x in items if id(x) not in gone]
