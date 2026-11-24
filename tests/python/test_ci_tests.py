"""What CI runs: the tests that what changed can affect (scripts/ci_tests.py, scripts/ci_record.py), on one runner each, and
the slow tests that pull requests leave out (conftest.py, .github/slow-tests.txt)."""

import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parents[1]

pytestmark = [pytest.mark.skipif(not (ROOT / "scripts" / "ci_tests.py").is_file(), reason="scripts are not in this tree"),
              pytest.mark.skipif(sys.version_info < (3, 12), reason="CI's Python has sys.monitoring")]


def _load(path):
    spec = importlib.util.spec_from_file_location(f"_{path.stem}_under_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def ci():
    return _load(ROOT / "scripts" / "ci_tests.py")


# -- what counts as a change ------------------------------------------------------------------------------------------

BEFORE = '''"""A module."""
import math

SCALE = 2.0


def used(x, y=1):
    """What it does."""
    return helper(x) + y  # a comment


def helper(x):
    return math.sqrt(x)


class Thing:
    """A thing."""
    size = 3

    def area(self):
        return self.size ** 2
'''


@pytest.mark.parametrize("old, new, changes", [
    ("# a comment", "# another, longer one", set()),
    ('"""What it does."""', '"""What it really does, in\n    more lines."""', set()),
    ("helper(x) + y", "helper( x )+y", set()),
    ("math.sqrt(x)", "math.sqrt(x + 1)", {"helper"}),
    ("y=1", "y=2", {"used"}),
    ("SCALE = 2.0", "SCALE = 3.0", {"rest"}),
    ("size = 3", "size = 4", {"rest"}),
    ("self.size ** 2", "self.size * 2", {"Thing.area"}),
    ("class Thing:", "def new(x):\n    return x\n\n\nclass Thing:", set()),  # nobody can call it yet
    ("class Thing:", "@register\ndef new(x):\n    return x\n\n\nclass Thing:", {"rest"}),
    ("class Thing:", "def new(x=used(1)):\n    return x\n\n\nclass Thing:", {"rest"}),
    ("class Thing:", "def len(x):\n    return 1\n\n\nclass Thing:", {"rest"}),
    ("class Thing:", "def math(x):\n    return 1\n\n\nclass Thing:", {"rest"}),
    ("size = 3", "size = 3\n\n    def other(self):\n        return 0", {"rest"}),
])
def test_only_what_runs_is_a_change(ci, old, new, changes):
    assert BEFORE.count(old) == 1
    a, b = ci.Code(BEFORE), ci.Code(BEFORE.replace(old, new))
    got = {u for u in a.units if b.units.get(u) != a.units[u]} | ({"rest"} if a.rest != b.rest else set())
    assert got == changes
    assert (a.whole == b.whole) == (not changes and "def new" not in new)
    assert ci.same("defs:x", " ".join(a.defs), " ".join(b.defs)) and "used" in a.defs


def test_a_function_that_is_gone_is_a_change(ci):
    """Someone imports it, perhaps without calling it."""
    a, b = ci.Code(BEFORE), ci.Code(BEFORE.replace("def helper", "def renamed"))
    assert not ci.same("defs:x", " ".join(a.defs), " ".join(b.defs))
    assert not ci.same("defs:x", "", "-") and not ci.same("unit:x", "a", "a b")


def test_docstrings_that_a_model_reads_are_code(ci):
    assert ci.keeps_docstrings("python/physim/mcp/tools_setup.py") and not ci.keeps_docstrings("python/physim/plot.py")
    new = BEFORE.replace("What it does.", "What else.")
    assert ci.Code(BEFORE, True).units["used"] != ci.Code(new, True).units["used"]


def test_the_place_of_code_is_known(ci):
    code = ci.Code(BEFORE)
    at = {line: code.unit_at(n) for n, line in enumerate(BEFORE.splitlines(), 1)}
    assert at["def used(x, y=1):"] == "used" and at["    def area(self):"] == "Thing.area"
    assert at["class Thing:"] is None and at["SCALE = 2.0"] is None
    twice = ci.Code("class A:\n    @property\n    def x(self):\n        return 1\n\n    @x.setter\n    def x(self, v):\n        pass\n")
    assert list(twice.units) == ["A.x", "A.x#2"] and twice.unit_at(6) == "A.x#2"


# -- records, end to end on a small repository ------------------------------------------------------------------------


@pytest.fixture
def repo(tmp_path, ci):
    """A repository with two test files, and a function that runs, records and chooses in it."""
    def write(path, text):
        (tmp_path / path).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / path).write_text(text, encoding="utf-8")

    def git(*args):
        subprocess.run(["git", "-c", "user.name=x", "-c", "user.email=x@example.org", *args], cwd=tmp_path, check=True,
                       capture_output=True)

    for name in ("ci_tests.py", "ci_record.py"):
        write(f"scripts/{name}", (ROOT / "scripts" / name).read_text(encoding="utf-8"))
    write("tests/python/conftest.py", (TESTS / "conftest.py").read_text(encoding="utf-8"))
    write(".github/slow-tests.txt", "test_a.py::test_slow\n")
    write(".github/test-durations.json", '{"full": {}, "quick": {}}')
    write("python/physim/__init__.py", "")
    write("python/physim/one.py", BEFORE + "\n\ndef read(p):\n    return open(p).read()\n")
    write("python/physim/two.py", "from .one import helper\n\n\ndef other():\n    return 1\n")
    write("data/a.csv", "1\n")
    write("README.md", "text\n")
    start = "import sys\nfrom pathlib import Path\nsys.path.insert(0, str(Path(__file__).parents[2] / 'python'))\n"
    write("tests/python/test_a.py", start + "from physim import one\n\n\ndef test_a():\n    assert one.used(4) == 3\n"
          "    assert one.read(Path(__file__).parents[2] / 'data' / 'a.csv')\n\n\ndef test_slow():\n    one.Thing().area()\n")
    write("tests/python/test_b.py", start + "import subprocess\nfrom physim import two\n\n\ndef test_b():\n"
          "    assert two.other() == 1\n")
    git("init", "-q")
    git("add", "-A")
    git("commit", "-q", "-m", "x")

    class Repo:
        path = tmp_path

        def edit(self, path, old, new):
            text = (tmp_path / path).read_text(encoding="utf-8")
            assert old in text
            write(path, text.replace(old, new))

        def call(self, *args):
            env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PHYSIM_CI_RECORD")}
            return subprocess.run([sys.executable, "scripts/ci_tests.py", *args], cwd=tmp_path, capture_output=True,
                                  text=True, env=env)

        def chosen(self, *stores, event="pull_request"):
            got = self.call("select", "--event", event, "--base", "HEAD", "--store", *stores)
            assert got.returncode == 0, got.stderr
            out = dict(x.split("=", 1) for x in got.stdout.splitlines())
            assert out["python"] == ("true" if json.loads(out["shards"]) else "false")
            return sorted(t for s in json.loads(out["shards"]) for t in s["tests"].split())

    return Repo()


def test_a_run_is_not_repeated_until_what_it_used_changes(repo):
    both = ["test_a.py:quick", "test_b.py:quick"]
    assert repo.chosen() == both
    ran = repo.call("run", "--records", "r.json", *both)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert repo.chosen("r.json") == []
    # Text: a document, a comment, a docstring.
    repo.edit("README.md", "text", "other text")
    repo.edit("python/physim/one.py", "# a comment", "# more")
    repo.edit("python/physim/one.py", "What it does.", "What it did.")
    repo.edit("tests/python/test_b.py", "def test_b():", "# a note\ndef test_b():")
    assert repo.chosen("r.json") == []
    # A function that a test ran, through another: that test alone, and one of its names is enough to run it whole.
    repo.edit("python/physim/two.py", "return 1", "return 2 - 1")
    assert repo.chosen("r.json") == ["test_b.py:quick"]
    repo.edit("python/physim/one.py", "math.sqrt(x)", "math.sqrt(x) + 0")
    assert repo.chosen("r.json") == both
    # The fix that only one needed: after the other passed on it, only that one is left.
    assert repo.call("run", "--records", "r2.json", "test_a.py:quick").returncode == 0
    assert repo.chosen("r.json", "r2.json") == ["test_b.py:quick"]


@pytest.mark.parametrize("path, old, new, runs", [
    ("data/a.csv", "1", "2", ["test_a.py:quick"]),
    ("python/physim/one.py", "SCALE = 2.0", "SCALE = 3.0", ["test_a.py:quick", "test_b.py:quick"]),
    ("python/physim/one.py", "self.size ** 2", "self.size", []),  # only the slow test, which a pull request leaves out
    ("python/physim/one.py", "def read", "def read_it", ["test_a.py:quick", "test_b.py:quick"]),
    ("tests/python/test_a.py", "== 3", "== 4", ["test_a.py:full"]),
    (".github/slow-tests.txt", "test_a.py::test_slow\n", "", ["test_a.py:quick"]),
    ("scripts/ci_record.py", "TOOL = 4", "TOOL = 3", ["test_a.py:quick", "test_b.py:quick"]),
])
def test_what_a_change_runs(repo, path, old, new, runs):
    assert repo.call("run", "--records", "r.json", "test_a.py:quick", "test_b.py:quick").returncode == 0
    repo.edit(path, old, new)
    assert repo.chosen("r.json") == runs


def test_a_failure_is_not_recorded_and_does_not_stop_the_others(repo):
    repo.edit("python/physim/two.py", "return 1", "return 3")
    ran = repo.call("run", "--records", "r.json", "test_b.py:quick", "test_a.py:quick")
    assert ran.returncode == 1 and "failed: test_b.py" in ran.stdout
    assert repo.chosen("r.json") == ["test_b.py:quick"]
    repo.edit("python/physim/two.py", "return 3", "return 1")  # the fix, in what the other test doesn't use
    assert repo.chosen("r.json") == ["test_b.py:quick"]


def test_what_starts_a_process_uses_the_whole_package(repo):
    repo.edit("tests/python/test_b.py", "== 1", "== 1\n    subprocess.run([sys.executable, '-c', 'pass'], check=True)")
    assert repo.call("run", "--records", "r.json", "test_b.py:full").returncode == 0
    assert repo.chosen("r.json") == ["test_a.py:quick"]
    repo.edit("python/physim/one.py", "self.size ** 2", "self.size")
    assert repo.chosen("r.json") == ["test_a.py:quick", "test_b.py:full"]


def test_main_runs_everything_or_nothing(repo):
    full = ["test_a.py:full", "test_b.py:full"]
    assert repo.call("run", "--records", "r.json", *full).returncode == 0
    assert repo.chosen("r.json", event="push") == [] and repo.chosen("r.json") == []  # a whole run counts for a quick one
    assert repo.chosen("r.json", event="workflow_dispatch") == full and repo.chosen("r.json", event="tag") == full
    repo.edit("README.md", "text", "other")
    assert repo.chosen("r.json", event="push") == []
    assert repo.call("code-changed", "HEAD").stdout.strip() == "false"
    repo.edit("data/a.csv", "1", "2")
    assert repo.chosen("r.json", event="push") == full
    assert repo.call("code-changed", "HEAD").stdout.strip() == "true"


def test_stores_are_merged_and_a_broken_one_counts_for_nothing(repo, ci):
    assert repo.call("run", "--records", "r.json", "test_a.py:quick").returncode == 0
    (repo.path / "broken.json").write_text("{")
    assert repo.call("merge", "--out", "m.json", "broken.json", "r.json", "r.json").returncode == 0
    assert len(json.loads((repo.path / "m.json").read_text())["records"]) == 1
    assert repo.chosen("m.json", "broken.json", "missing.json") == ["test_b.py:quick"]
    many = [{"file": "a", "mode": "quick", "used": {"x": str(i)}} for i in range(20)]
    assert ci.merged(many) == many[-ci.KEPT:]


# -- sharing out, and the slow tests ----------------------------------------------------------------------------------


def test_every_test_file_is_on_one_runner(ci):
    tests = [(p.name, mode) for p in sorted(TESTS.glob("test_*.py")) for mode in ("quick",)]
    times = ci.known_times([])
    for shards in (None, 1, 4, 8):
        plan = ci.plan(tests, times, shards)
        assert sorted(f for _, part in plan for f in part) == sorted(tests) and len(plan) <= ci.MAX_SHARDS
    assert ci.plan([], times) == [] and len(ci.plan(tests[:1], times)) == 1


def test_a_pull_request_runs_whole_what_it_changes(ci):
    files = sorted(p.name for p in TESTS.glob("test_*.py"))
    got = ci.slow_in(["tests/python/test_nuclear_validation.py", "python/physim/nuclear/analysis.py", "CHANGELOG.md",
                      "tests/python/test_that_was_deleted.py", "docs/analysis.md"], files)
    assert "test_nuclear_validation.py" in got and "test_that_was_deleted.py" not in got
    assert {"test_nuclear_analysis.py", "test_nuclear_analysis_closure.py", "test_nuclear_analysis_run.py"} < set(got)
    assert "test_nuclear_report.py" not in got and ci.slow_in(["README.md"], files) == []


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
