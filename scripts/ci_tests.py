"""Which tests CI has to run, and running them.

Usage: python scripts/ci_tests.py select --event EVENT [--base REV] [--store STORE ...]   the outputs of CI's first job
       python scripts/ci_tests.py run --records OUT TEST ...                              one runner's tests
       python scripts/ci_tests.py passed --records OUT NAME                               another job passed ("rust")
       python scripts/ci_tests.py merge --out STORE [--store STORE ...] [OUT ...]         the records after a run
       python scripts/ci_tests.py code-changed REV                                        "true" unless only text did

**Records.** Each test file runs in a process of its own with ``scripts/ci_record.py``, which sees the functions, modules,
files and directory listings of this repository that it used. When it passes, they are kept with a hash of each. A file
whose record still holds (every hash is the same in the tree at hand) would do the same again, so it is not run: after a
test failed and was fixed, what runs is that file and whatever else used the code that the fix changed. Code is compared
without its comments and layout, and without docstrings (except those of ``physim.mcp``, which a model reads), function
by function; the rest of a module (its imports, constants, classes and their methods' signatures) is compared as one. A new
function in a module changes nothing for a test that doesn't reach it, but one that is gone may (it is imported somewhere).
Test files are compared whole. A test that starts another process holds everything under ``python/physim`` as used. The Rust sources, the
requirements and this machinery are used by all.

**Pull requests** are quick besides: they leave out the tests of ``.github/slow-tests.txt``, except in the test files they
change and in those named after a module they change (``analysis.py``: ``test_nuclear_analysis*.py``), which run whole.
**main** runs every test when any file's record no longer holds, and nothing when none is affected (a change of text).
**Tags and runs started by hand** run everything.

**What a record misses:** whether a file exists that was only looked for, new versions of the packages of the requirements,
and what a program outside the repository does. main's full runs and those by hand see those.

The files are shared out over runners (about ``SHARD_S`` seconds of tests each) by the time that they took last.
"""

import argparse
import ast
import builtins
import functools
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests" / "python"
DURATIONS = ROOT / ".github" / "test-durations.json"
SLOW_TESTS = ROOT / ".github" / "slow-tests.txt"
UNKNOWN_S, START_S = 10.0, 2.0  # a file without a time; what a process takes beside its tests
SHARD_S, MAX_SHARDS = 75.0, 10
KEPT = 8  # records of a file and mode

RUST = re.compile(r"(src/|benches/|examples/|tests/[^/]+\.rs$|Cargo\.|\.cargo/|[^/]*\.toml$|rust-toolchain)")
MACHINERY = ("requirements-dev.txt", "scripts/ci_tests.py", "scripts/ci_record.py", "tests/python/conftest.py")


def digest(data) -> str:
    return hashlib.sha256(data.encode("utf-8") if isinstance(data, str) else data).hexdigest()[:16]


# -- code without its text ------------------------------------------------------------------------------------------


def _without_docstrings(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            first = node.body[0] if node.body else None
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                node.body = node.body[1:] or [ast.Pass()]
    return tree


class Code:
    """A Python source: the hash of each of its functions and methods, of what is around them, and of all of it."""

    def __init__(self, text: str, docstrings: bool = False):
        tree = ast.parse(text)
        if not docstrings:
            _without_docstrings(tree)
        self.whole = digest(ast.dump(tree))
        self.units, self.spans = {}, []
        plain = self._plain(tree)
        self._take(tree, "")
        self.defs = sorted(f.name for f in plain)
        tree.body = [x for x in tree.body if x not in plain]
        self.rest = digest(ast.dump(tree))

    @staticmethod
    def _plain(tree) -> list:
        """The module's own functions that do nothing until they are called: without a decorator, with constants as
        defaults, and with a name that nothing else in the module binds."""
        def constant(x):
            return all(isinstance(n, (ast.Constant, ast.Tuple, ast.Load, ast.UnaryOp, ast.USub)) for n in ast.walk(x))

        bound = {}
        for node in ast.walk(tree):
            names = ([node.name] if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) else
                     [node.id] if isinstance(node, ast.Name) and not isinstance(node.ctx, ast.Load) else
                     [(node.asname or node.name).partition(".")[0]] if isinstance(node, ast.alias) else [])
            for n in names:
                bound[n] = bound.get(n, 0) + 1
        return [x for x in tree.body if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef)) and not x.decorator_list
                and bound[x.name] == 1 and not hasattr(builtins, x.name) and not x.name.startswith("__")
                and all(constant(d) for d in [*x.args.defaults, *filter(None, x.args.kw_defaults)])]

    def _take(self, node, prefix):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                name = key = prefix + child.name
                n = 1
                while key in self.units:  # a property's setter, a function defined twice
                    n += 1
                    key = f"{name}#{n}"
                self.units[key] = digest(ast.dump(child))
                first = min([child.lineno] + [d.lineno for d in child.decorator_list])
                self.spans.append((first, child.end_lineno, key))
                child.body = [ast.Pass()]
            elif isinstance(child, ast.ClassDef):
                self._take(child, f"{prefix}{child.name}.")
            elif isinstance(child, ast.stmt) or isinstance(child, (ast.ExceptHandler, ast.match_case)):
                self._take(child, prefix)

    def unit_at(self, line: int):
        return next((key for first, last, key in self.spans if first <= line <= last), None)


def keeps_docstrings(path: str) -> bool:
    return path.startswith("python/physim/mcp/")


# -- a tree ---------------------------------------------------------------------------------------------------------


class Tree:
    """The tracked files as they are in the working directory (``rev``: in a commit)."""

    def __init__(self, rev=None):
        self.rev = rev
        out = _git("ls-files", "-z") if rev is None else _git("ls-tree", "-r", "-z", "--name-only", rev)
        self.files = sorted(x for x in out.split("\0") if x and (rev or (ROOT / x).is_file()))
        self._set = set(self.files)
        self._hashes = {}

    def read(self, path: str) -> bytes:
        if self.rev is None:
            return (ROOT / path).read_bytes().replace(b"\r\n", b"\n")
        return subprocess.run(["git", "show", f"{self.rev}:{path}"], cwd=ROOT, capture_output=True,
                              check=True).stdout.replace(b"\r\n", b"\n")

    @functools.lru_cache(maxsize=None)
    def code(self, path: str):
        """:class:`Code` of a file, or None if it is not here or not Python that can be read."""
        if path not in self._set:
            return None
        try:
            return Code(self.read(path).decode("utf-8"), keeps_docstrings(path))
        except (SyntaxError, UnicodeDecodeError, ValueError):
            return None

    @functools.lru_cache(maxsize=None)
    def bytes_hash(self, path: str) -> str:
        return digest(self.read(path)) if path in self._set else "-"

    def meaning(self, path: str) -> str:
        """The hash of a file, without its text if it is Python."""
        code = self.code(path) if path.endswith(".py") else None
        return code.whole if code else self.bytes_hash(path)

    def _all(self, paths) -> str:
        return digest("\n".join(f"{p} {self.meaning(p)}" for p in paths))

    def hash(self, key: str) -> str:
        if key not in self._hashes:
            self._hashes[key] = self._hash(key)
        return self._hashes[key]

    def _hash(self, key: str) -> str:
        kind, _, arg = key.partition(":")
        if kind == "unit":
            path, _, unit = arg.partition("::")
            code = self.code(path)
            return code.units.get(unit, "-") if code else "-"
        if kind == "rest":
            code = self.code(arg)
            return code.rest if code else self.bytes_hash(arg)
        if kind == "defs":
            code = self.code(arg)
            return " ".join(code.defs) if code else "-"
        if kind == "whole":
            return self.meaning(arg)
        if kind == "file":
            return self.bytes_hash(arg)
        if kind == "list":
            prefix = "" if arg == "." else arg + "/"
            return digest("\n".join(sorted({p[len(prefix):].partition("/")[0] for p in self.files
                                            if p.startswith(prefix)})))
        if kind == "package":
            return self._all(p for p in self.files if p.startswith("python/physim/"))
        if kind == "rust":
            return self._all(p for p in self.files if RUST.match(p))
        if kind == "machinery":
            return digest(f"{sys.version_info[:2]} {self.hash('rust:')} " + self._all(p for p in self.files if p in MACHINERY))
        if kind == "slow":
            text = self.read(".github/slow-tests.txt").decode() if ".github/slow-tests.txt" in self._set else ""
            return digest("\n".join(sorted(x.strip() for x in text.splitlines() if x.startswith(f"{arg}::"))))
        raise ValueError(key)

    def tests(self) -> list:
        return sorted(PurePosixPath(p).name for p in self.files if re.fullmatch(r"tests/python/test_[^/]*\.py", p))


def _git(*args) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True).stdout.decode("utf-8")


def used_by(seen: dict, file: str, mode: str, tree: Tree) -> dict:
    """{key: hash} of what ``ci_record`` saw a test file use."""
    own = f"tests/python/{file}"
    keys = {"machinery:", f"whole:{own}"}
    if mode == "quick":
        keys.add(f"slow:{file}")
    for path, line, name in seen["code"]:
        code = tree.code(path) if path != own else None
        if code is None:
            continue  # the file itself, or not of this repository
        keys.update((f"rest:{path}", f"defs:{path}"))
        unit = code.unit_at(line) if name != "<module>" else None
        if unit:
            keys.add(f"unit:{path}::{unit}")
    known = set(tree.files)
    keys.update(f"file:{p}" for p in seen["files"] if p in known and p != own)
    directories = {str(parent) for p in tree.files for parent in PurePosixPath(p).parents}
    keys.update(f"list:{p}" for p in seen["lists"] if p in directories)
    if seen["spawned"]:
        keys.add("package:")
    if file.startswith("test_ci_"):  # what reads the list of slow tests for pytest is tested there
        keys.update(("file:.github/slow-tests.txt", "list:tests/python"))
    return {k: tree.hash(k) for k in sorted(keys)}


# -- records --------------------------------------------------------------------------------------------------------


def load(paths) -> list:
    """The records of the stores and runs, the oldest first. One that cannot be read counts for nothing."""
    out = []
    for p in paths:
        try:
            out += json.loads(Path(p).read_text(encoding="utf-8"))["records"]
        except (OSError, ValueError, KeyError, TypeError):
            print(f"no records in {p}", file=sys.stderr)
    return out


def save(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps({"records": records}, separators=(",", ":")), encoding="utf-8")


def merged(records) -> list:
    """The newest ``KEPT`` different records of each file and mode."""
    out, seen, count = [], set(), {}
    for r in reversed(records):
        what = digest(json.dumps([r["file"], r["mode"], r["used"]], sort_keys=True))
        of = (r["file"], r["mode"])
        if what not in seen and count.get(of, 0) < KEPT:
            seen.add(what)
            count[of] = count.get(of, 0) + 1
            out.append(r)
    return out[::-1]


def same(key: str, then: str, now: str) -> bool:
    return now == then or key.startswith("defs:") and now != "-" and set(then.split()) <= set(now.split())


def holds(records, file: str, mode: str, tree: Tree) -> bool:
    """Whether a test file passed on what it would use in this tree (a run of all its tests counts for a quick one)."""
    return any(r["file"] == file and r["mode"] in (mode, "full")
               and all(same(k, h, tree.hash(k)) for k, h in r["used"].items()) for r in records)


# -- selection ------------------------------------------------------------------------------------------------------


def changed_in_meaning(tree: Tree, base: Tree) -> list:
    """The paths that differ by more than text: all but comments, layout and docstrings of Python."""
    paths = [p for p in _git("diff", "--name-only", "-z", "--no-renames", base.rev).split("\0") if p]
    return [p for p in paths if tree.meaning(p) != base.meaning(p)]


def is_text(path: str) -> bool:
    return path.endswith((".md", ".rst", ".txt")) and not path.endswith("requirements-dev.txt") or path.startswith(
        ("docs/", "backlog/", "planning/", ".github/"))


def slow_in(changed, files) -> list:
    """The test files that a pull request runs whole, of the paths it changes: those changed, and those named after a
    changed module of physim."""
    out = set()
    for path in map(PurePosixPath, changed):
        if path.suffix != ".py":
            continue
        if path.parts[:2] == ("tests", "python"):
            out.update({path.name} & set(files))
        elif path.parts[:2] == ("python", "physim"):
            about = re.compile(rf"test_(nuclear_)?{re.escape(path.stem)}(_|\.)")
            out.update(f for f in files if about.match(f))
    return sorted(out)


def plan(tests, times, shards=None) -> list:
    """[(seconds, [(file, mode)])] of each runner: the longest first, to the one with the least so far."""
    if not tests:
        return []
    t = {x: times.get(x, UNKNOWN_S) + START_S for x in tests}
    shards = shards or min(max(math.ceil(sum(t.values()) / SHARD_S), 1), MAX_SHARDS, len(tests))
    out = [[0.0, []] for _ in range(shards)]
    for x in sorted(tests, key=lambda x: (-t[x], x)):
        least = min(out, key=lambda b: b[0])
        least[0] += t[x]
        least[1].append(x)
    return [(s, sorted(f)) for s, f in out]


def known_times(records) -> dict:
    """{(file, mode): seconds}: of the newest record, or of ``.github/test-durations.json``."""
    known = json.loads(DURATIONS.read_text(encoding="utf-8"))
    out = {(f, mode): s for mode in ("full", "quick") for f, s in {**known["full"], **known[mode]}.items()}
    out.update({(r["file"], r["mode"]): r["seconds"] for r in records if "seconds" in r})
    return out


def choose(event: str, tree: Tree, base, records) -> dict:
    """What a run of CI does: the test files with their modes, and whether Rust's checks run."""
    files = tree.tests()
    if event not in ("pull_request", "push"):
        return {"tests": [(f, "full") for f in files], "rust": True}
    rust = not holds(records, "(rust)", "full", tree)
    whole = set(slow_in(changed_in_meaning(tree, base), files)) if event == "pull_request" and base else set()
    modes = {f: "full" if event == "push" or f in whole else "quick" for f in files}
    tests = [(f, mode) for f, mode in modes.items() if not holds(records, f, mode, tree)]
    if event == "push" and tests:
        tests = list(modes.items())
    return {"tests": tests, "rust": rust}


def select(args) -> int:
    tree = Tree()
    records = load(args.store)
    got = choose(args.event, tree, Tree(args.base) if args.base else None, records)
    shards = plan(got["tests"], known_times(records))
    print(f"{len(got['tests'])} of {len(tree.tests())} test files on {len(shards)} runners; "
          f"{len(records)} records", file=sys.stderr)
    for i, (seconds, tests) in enumerate(shards, 1):
        print(f"  {i}: {seconds:.0f} s: {' '.join(f if m == 'quick' else f + ' (whole)' for f, m in tests)}",
              file=sys.stderr)
    print("shards=" + json.dumps([{"index": i, "tests": " ".join(f"{f}:{m}" for f, m in tests)}
                                  for i, (_, tests) in enumerate(shards, 1)]))
    print(f"python={'true' if shards else 'false'}")
    print(f"rust={'true' if got['rust'] else 'false'}")
    return 0


def run(args) -> int:
    tree, records, failed = Tree(), [], []
    for x in args.tests:
        file, _, mode = x.partition(":")
        mode = mode or "full"
        command = [sys.executable, "-m", "pytest", f"tests/python/{file}", "-p", "ci_record", "-q", "-rs", "--durations=5",
                   "--durations-min=5", *(["--quick"] if mode == "quick" else [])]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "used.json"
            path = os.pathsep.join(filter(None, [str(ROOT / "scripts"), os.environ.get("PYTHONPATH")]))
            print(f"=== {file} ({mode})", flush=True)
            began = time.perf_counter()
            code = subprocess.run(command, cwd=ROOT, env={**os.environ, "PYTHONPATH": path,
                                                          "PHYSIM_CI_RECORD": str(out)}).returncode
            seconds = time.perf_counter() - began
            if code not in (0, 5):  # 5: no test left, as all of them are slow
                failed.append(file)
                continue
            seen = json.loads(out.read_text(encoding="utf-8"))
        records.append({"file": file, "mode": mode, "seconds": round(seconds, 1), "used": used_by(seen, file, mode, tree)})
    save(args.records, records)
    print(f"{len(records)} of {len(args.tests)} test files passed" + (f"; failed: {' '.join(failed)}" if failed else ""))
    return 1 if failed else 0


def main(argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__.partition("\n")[0])
    sub = parser.add_subparsers(dest="what", required=True)
    p = sub.add_parser("select")
    p.add_argument("--event", required=True)
    p.add_argument("--base")
    p.add_argument("--store", nargs="*", default=[])
    p = sub.add_parser("run")
    p.add_argument("--records", required=True)
    p.add_argument("tests", nargs="+")
    p = sub.add_parser("passed")
    p.add_argument("--records", required=True)
    p.add_argument("name", choices=["rust"])
    p = sub.add_parser("merge")
    p.add_argument("--out", required=True)
    p.add_argument("--store", nargs="*", default=[])
    p.add_argument("runs", nargs="*")
    p = sub.add_parser("code-changed")
    p.add_argument("base")
    args = parser.parse_args(argv)
    if args.what == "select":
        return select(args)
    if args.what == "run":
        return run(args)
    if args.what == "passed":
        save(args.records, [{"file": f"({args.name})", "mode": "full", "used": {f"{args.name}:": Tree().hash(f"{args.name}:")}}])
    elif args.what == "merge":
        save(args.out, merged(load([*args.store, *args.runs])))
    else:
        changed = [p for p in changed_in_meaning(Tree(), Tree(args.base)) if not is_text(p)]
        print("true" if changed else "false")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
