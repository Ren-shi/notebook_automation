"""What CI runs: every test file on exactly one runner, and the slow tests that pull requests leave out
(scripts/ci_shard.py, conftest.py, .github/slow-tests.txt)."""

import importlib.util
import re
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parents[1]

pytestmark = pytest.mark.skipif(not (ROOT / "scripts" / "ci_shard.py").is_file(), reason="scripts are not in this tree")


def _load(path):
    spec = importlib.util.spec_from_file_location(f"_{path.stem}_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def shard():
    return _load(ROOT / "scripts" / "ci_shard.py")


@pytest.mark.parametrize("quick, whole", [(False, ()), (True, ()), (True, ("test_nuclear_validation.py",))])
def test_every_test_file_is_on_one_runner(shard, quick, whole):
    files = sorted(p.name for p in TESTS.glob("test_*.py"))
    for shards in (1, 4, 8):
        plan = shard.plan(shards, quick, whole)
        assert sorted(f for _, part in plan for f in part) == files and len(plan) == shards


def test_a_pull_request_runs_whole_what_it_changes(shard):
    got = shard.slow_in(["tests/python/test_nuclear_validation.py", "python/physim/nuclear/analysis.py", "CHANGELOG.md",
                         "tests/python/test_that_was_deleted.py", "docs/analysis.md"])
    assert "test_nuclear_validation.py" in got and "test_that_was_deleted.py" not in got
    assert {"test_nuclear_analysis.py", "test_nuclear_analysis_closure.py", "test_nuclear_analysis_run.py"} < set(got)
    assert "test_nuclear_report.py" not in got and shard.slow_in(["README.md"]) == []


def test_the_slow_tests_exist():
    """A renamed test would come back to the pull requests without anyone noticing: its line must go with it."""
    entries = _load(TESTS / "conftest.py").slow_entries()
    assert entries
    for x in sorted(entries):
        file, _, name = x.partition("::")
        source = (TESTS / file).read_text(encoding="utf-8")
        name = name.removeprefix("fixture:").partition("[")[0]
        assert re.search(rf"^def {re.escape(name)}\(", source, re.M), x


def test_what_counts_as_slow():
    is_slow = _load(TESTS / "conftest.py").is_slow
    entries = {"test_a.py::test_one", "test_a.py::test_two[b]", "test_a.py::fixture:big"}
    assert is_slow("test_a.py", "test_one", (), entries) and is_slow("test_a.py", "test_one[x]", (), entries)
    assert is_slow("test_a.py", "test_two[b]", (), entries) and not is_slow("test_a.py", "test_two[c]", (), entries)
    assert is_slow("test_a.py", "test_three", ("tmp_path", "big"), entries)
    assert not is_slow("test_b.py", "test_one", ("big",), entries)
