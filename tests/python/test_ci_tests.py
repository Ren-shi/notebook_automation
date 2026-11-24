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
from __future__ import annotations
import math
from dataclasses import dataclass

SCALE = 2.0
TABLE = {}
TABLE["a"] = 1


def used(x: Thing, y=1):
    """What it does."""
    return helper(x) + y  # a comment


def helper(x):
    return math.sqrt(x)


@dataclass
class Thing:
    """A thing."""
    size: int = 3

    def area(self):
        return self.size ** 2


class Other(Thing):
    pass
'''


@pytest.mark.parametrize("old, new, changes", [
    ("# a comment", "# another, longer one", set()),
    ('"""What it does."""', '"""What it really does, in\n    more lines."""', set()),
    ("helper(x) + y", "helper( x )+y", set()),
    ("math.sqrt(x)", "math.sqrt(x + 1)", {"helper"}),
    ("y=1", "y=2", {"used"}),
    ("SCALE = 2.0", "SCALE = 3.0", {"name SCALE"}),
    ('TABLE["a"] = 1', 'TABLE["a"] = 2', {"rest"}),
    ("import math", "import cmath as math", {"name math"}),
    ("import math", "import math\nimport json", set()),
    ("size: int = 3", "size: int = 4", {"name Thing"}),
    ("self.size ** 2", "self.size * 2", {"Thing.area"}),
    ("SCALE = 2.0", "SCALE = 2.0\nMORE = 1", set()),
    ("def helper(x):", "def new(x):\n    return x\n\n\ndef helper(x):", set()),  # nobody can call it yet
    ("def helper(x):", "@register\ndef new(x):\n    return x\n\n\ndef helper(x):", {"rest"}),
    ("def helper(x):", "@lru_cache\ndef helper(x):", {"helper", "name helper"}),
    ("def helper(x):", "def len(x):\n    return 1\n\n\ndef helper(x):", {"rest"}),
    ("def helper(x):", "def math(x):\n    return 1\n\n\ndef helper(x):", {"rest"}),
    ("self.size ** 2", "self.size ** 2\n\n    def other(self):\n        return 0", set()),  # of a class without bases
    ("self.size ** 2", "self.size ** 2\n\n    def __len__(self):\n        return 0", {"name Thing"}),
    ("class Other(Thing):\n    pass", "class Other(Thing):\n    def area(self):\n        return 0", {"name Other"}),
    ("def helper", "def renamed", {"helper", "name helper"}),
])
def test_only_what_runs_is_a_change(ci, old, new, changes):
    assert BEFORE.count(old) == 1
    a, b = ci.Code(BEFORE), ci.Code(BEFORE.replace(old, new))
    got = {u for u in a.units if b.units.get(u, ("-",))[0] != a.units[u][0]} | ({"rest"} if a.rest != b.rest else set())
    got |= {f"name {n}" for n, (value, _) in a.names.items()
            if not ci.same(f"name:{n}", value, b.names[n][0] if n in b.names else "-")}
    assert got == changes or got - changes <= {"name math"} and "def math" in new


def test_a_class_that_lost_a_method_is_a_change(ci):
    a, b = ci.Code(BEFORE), ci.Code(BEFORE.replace("def area", "def size2"))
    assert not ci.same("name:x", a.names["Thing"][0], b.names["Thing"][0])
    assert not ci.same("unit:x", "a", "a b") and not ci.same("name:x", "a", "-")


def test_what_code_names(ci):
    code = ci.Code(BEFORE)
    assert ("helper",) in code.units["used"][1].names and ("Thing",) not in code.units["used"][1].names  # an annotation
    assert ("math", "sqrt") in code.units["helper"][1].names
    assert ("Thing",) in code.names["Other"][1][0].names
    refs = ci.Code("def f():\n    from .a import X as Y\n    return Y + g('physim.nuclear.data.Z')\n").units["f"][1]
    assert refs.imports == [("Y", ".a", "X")] and ("physim", "nuclear", "data", "Z") in refs.names


def test_docstrings_that_a_model_reads_are_code(ci):
    assert ci.keeps_docstrings("python/physim/mcp/tools_setup.py") and not ci.keeps_docstrings("python/physim/plot.py")
    new = BEFORE.replace("What it does.", "What else.")
    assert ci.Code(BEFORE, True).units["used"][0] != ci.Code(new, True).units["used"][0]


def test_the_place_of_code_is_known(ci):
    code = ci.Code(BEFORE)
    at = {line: code.unit_at(n) for n, line in enumerate(BEFORE.splitlines(), 1)}
    assert at["def used(x: Thing, y=1):"] == "used" and at["    def area(self):"] == "Thing.area"
    assert at["class Thing:"] is None and at["SCALE = 2.0"] is None
    twice = ci.Code("class A:\n    @property\n    def x(self):\n        return 1\n\n    @x.setter\n    def x(self, v):\n        pass\n")
    assert list(twice.units) == ["A.x", "A.x#2"] and twice.unit_at(6) == "A.x#2"


def test_names_are_followed_through_the_package(ci):
    """In physim itself: an import of a package's __init__, a module's attribute, a constant's own constants."""
    tree = ci.Tree()
    assert tree.broken_imports() == ()
    refs = ci.Refs(__import__("ast").parse("from physim.nuclear import Experiment, stopping\nstopping.Stopping\n"))
    keys = tree.follow([("tests/python/test_ci_tests.py", refs)])
    assert "name:python/physim/nuclear/experiment.py::Experiment" in keys
    assert "name:python/physim/nuclear/__init__.py::Experiment" in keys
    assert "name:python/physim/nuclear/stopping.py::Stopping" in keys
    assert "name:python/physim/nuclear/stopping.py::known_shortfall" not in keys
    assert not any(k.startswith("skeleton:") for k in keys)  # no module as a whole


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
    write("python/physim/two.py", "from . import one, three\nfrom .one import helper\n\nOWN = three.DEEP + 1\n\n\n"
          "def other():\n    return OWN - one.SCALE - 4\n")
    write("python/physim/three.py", "DEEP = 6\nUNUSED = 1\n\n\ndef nobody():\n    return 2\n")
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
    # A function that a test ran.
    repo.edit("python/physim/two.py", "return OWN -", "return +OWN -")
    assert repo.chosen("r.json") == ["test_b.py:quick"]
    repo.edit("python/physim/one.py", "math.sqrt(x)", "math.sqrt(x) + 0")
    assert repo.chosen("r.json") == both
    # After one passed on it, only the other is left.
    assert repo.call("run", "--records", "r2.json", "test_a.py:quick").returncode == 0
    assert repo.chosen("r.json", "r2.json") == ["test_b.py:quick"]


A, B = "test_a.py:quick", "test_b.py:quick"


@pytest.mark.parametrize("path, old, new, runs", [
    ("data/a.csv", "1", "2", [A]),
    ("python/physim/one.py", "SCALE = 2.0", "SCALE = 3.0", [B]),  # through the module that has it
    ("python/physim/three.py", "DEEP = 6", "DEEP = 7", [B]),  # through a constant of another module
    ("python/physim/three.py", "UNUSED = 1", "UNUSED = 2", []),
    ("python/physim/three.py", "return 2", "return 3", []),
    ("python/physim/three.py", "DEEP = 6", "DEEP = 6\nif DEEP:\n    DEEP = 7", [B]),
    ("python/physim/one.py", 'TABLE["a"] = 1', 'TABLE["a"] = 2', [A, B]),
    ("python/physim/one.py", "self.size ** 2", "self.size", []),  # only the slow test, which a pull request leaves out
    ("python/physim/one.py", "def read", "def read_it", [A]),
    ("python/physim/one.py", "def helper", "def helper_it", [A, B]),  # an import that fails
    ("tests/python/test_a.py", "== 3", "== 4", [A, "test_a.py:test_slow"]),
    (".github/slow-tests.txt", "test_a.py::test_slow\n", "", [A]),
    ("tests/python/conftest.py", "def slow_entries", "X = 1\n\n\ndef slow_entries", [A, B]),
    ("scripts/ci_record.py", "TOOL = 4", "TOOL = 3", [A, B]),
])
def test_what_a_change_runs(repo, path, old, new, runs):
    ran = repo.call("run", "--records", "r.json", A, B)
    assert ran.returncode == 0, ran.stdout
    repo.edit(path, old, new)
    assert repo.chosen("r.json") == runs


def test_a_failure_is_not_recorded_and_does_not_stop_the_others(repo):
    repo.edit("python/physim/three.py", "DEEP = 6", "DEEP = 8")
    ran = repo.call("run", "--records", "r.json", B, A)
    assert ran.returncode == 1 and "failed: test_b.py::quick" in ran.stdout
    assert repo.chosen("r.json") == [B]
    repo.edit("python/physim/three.py", "DEEP = 8", "DEEP = 6")  # the fix, in what the other test doesn't use
    assert repo.chosen("r.json") == [B]


def test_what_starts_a_program_uses_the_whole_package(repo):
    repo.edit("tests/python/test_b.py", "== 1", "== 1\n    subprocess.run([sys.executable, '-c', 'pass'], check=True)")
    assert repo.call("run", "--records", "r.json", B).returncode == 0
    assert repo.chosen("r.json") == [A]
    repo.edit("python/physim/three.py", "return 2", "return 3")
    assert repo.chosen("r.json") == [A, B]


def test_main_runs_everything_or_nothing(repo):
    every = [A, "test_a.py:test_slow", B]
    ran = repo.call("run", "--records", "r.json", *every)
    assert ran.returncode == 0, ran.stdout
    assert ran.stdout.count("1 passed, 1 deselected") == 2  # the slow test by itself, and the rest without it
    assert repo.chosen("r.json", event="push") == [] and repo.chosen("r.json") == []
    assert repo.chosen("r.json", event="workflow_dispatch") == every and repo.chosen("r.json", event="tag") == every
    repo.edit("README.md", "text", "other")
    assert repo.chosen("r.json", event="push") == []
    assert repo.call("code-changed", "HEAD").stdout.strip() == "false"
    repo.edit("data/a.csv", "1", "2")
    assert repo.chosen("r.json", event="push") == every
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


def test_every_part_is_on_one_runner(ci):
    tree = ci.Tree()
    tests, times = ci.parts(tree, tree.tests()), ci.known_times([], tree)
    assert {f for f, _ in tests} == set(tree.tests()) and len(tests) > len(tree.tests())
    for shards in (None, 1, 4, 8):
        plan = ci.plan(tests, times, shards)
        assert sorted(f for _, part in plan for f in part) == sorted(tests) and len(plan) <= ci.MAX_SHARDS
    assert ci.plan([], times) == [] and len(ci.plan(tests[:1], times)) == 1
    assert ci.token(("test_x.py", "test_a[one thing]")) == "test_x.py:test_a[one%20thing]"


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


def test_a_slow_test_is_in_one_entry():
    entry = _load(TESTS / "conftest.py").slow_entry
    entries = {"test_a.py::test_one", "test_a.py::test_two[b]", "test_a.py::fixture:big"}
    assert entry("test_a.py", "test_one", (), entries) == entry("test_a.py", "test_one[x]", (), entries) == "test_one"
    assert entry("test_a.py", "test_two[b]", (), entries) == "test_two[b]"
    assert entry("test_a.py", "test_two[c]", (), entries) is None
    assert entry("test_a.py", "test_three", ("tmp_path", "big"), entries) == "fixture:big"
    assert entry("test_a.py", "test_one", ("big",), entries) in ("fixture:big", "test_one")  # in one of them
    assert entry("test_b.py", "test_one", ("big",), entries) is None
