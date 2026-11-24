"""``--quick``: leave out the tests of .github/slow-tests.txt, as pull requests do (plain pytest runs everything)."""

import itertools
from pathlib import Path

import pytest

SLOW_TESTS = Path(__file__).resolve().parents[2] / ".github" / "slow-tests.txt"
_NUCLEI_HOMES = itertools.count()


def slow_entries() -> set:
    """The lines of the list: ``file::test``, ``file::test[parameters]`` or ``file::fixture:NAME``."""
    if not SLOW_TESTS.is_file():  # an unpacked source distribution
        return set()
    lines = (x.strip() for x in SLOW_TESTS.read_text(encoding="utf-8").splitlines())
    return {x for x in lines if x and not x.startswith("#")}


@pytest.fixture(autouse=True)
def nuclei_home_apart(tmp_path_factory, monkeypatch):
    """The results ledger's rows of a nucleus (shell-model results are copied there whenever one is read or
    applied) go to a folder of this test run, never to ``~/.physim/nuclei``. A test that wants its own sets
    ``PHYSIM_NUCLEI`` itself."""
    monkeypatch.setenv("PHYSIM_NUCLEI", str(tmp_path_factory.getbasetemp() / "nuclei" / str(next(_NUCLEI_HOMES))))


def slow_entry(file: str, name: str, fixtures, entries):
    """The first of the list's entries that a test is in, without its file; None if it is in none."""
    mine = (f"{file}::{name}", f"{file}::{name.partition('[')[0]}", *(f"{file}::fixture:{x}" for x in fixtures))
    return next((x.partition("::")[2] for x in sorted(entries) if x in mine), None)


def pytest_addoption(parser):
    parser.addoption("--quick", action="store_true", help="leave out the tests of .github/slow-tests.txt")
    parser.addoption("--slow-entry", metavar="ENTRY",
                     help="run only the tests of one entry of .github/slow-tests.txt (without its file)")


# Last, so that the whole collection is judged by its final order and fixtures.
@pytest.hookimpl(trylast=True)
def pytest_collection_modifyitems(config, items):
    quick, only = config.getoption("--quick"), config.getoption("--slow-entry")
    if not quick and only is None:
        return
    entries = slow_entries()
    out = [x for x in items if slow_entry(x.path.name, x.name, getattr(x, "fixturenames", ()), entries) != only]
    if out:
        config.hook.pytest_deselected(items=out)
        gone = set(map(id, out))
        items[:] = [x for x in items if id(x) not in gone]
