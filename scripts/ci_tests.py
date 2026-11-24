"""Which tests CI has to run, and running them.

Usage: python scripts/ci_tests.py select --event EVENT [--base REV] [--store STORE ...]   the outputs of CI's first job
       python scripts/ci_tests.py run --records OUT TEST ...                              one runner's tests
       python scripts/ci_tests.py passed --records OUT NAME                               another job passed ("rust")
       python scripts/ci_tests.py merge --out STORE [--store STORE ...] [OUT ...]         the records after a run
       python scripts/ci_tests.py code-changed REV                                        "true" unless only text did

**Records.** Each test file runs in a process of its own with ``scripts/ci_record.py``, which sees the functions, modules,
files and directory listings of this repository that it used. When it passes, they are kept with a hash of each. One whose
record still holds (every hash is the same in the tree at hand) would do the same again, so it is not run: after a test
failed and was fixed, what runs is that file and whatever else used the code that the fix changed. Code is compared without
its comments and layout, and without docstrings (except those of ``physim.mcp``, which a model reads), piece by piece: the
functions and methods that ran; the constants, imports and classes that they, the test and each other name (followed from
module to module); and of each module that was loaded the statements that do something else (``TABLE.update(...)``, an
``if``, a decorator that registers). A module that is handed around as a whole counts as a whole. A new function, or a new
method of a class without bases, changes nothing for a test that doesn't reach it. Test files are compared whole. A test that
starts a program (but the ``HARMLESS``) holds everything under ``python/physim`` as used. The Rust sources, the requirements
and this machinery are used by all.

**Pull requests** leave out the tests of ``.github/slow-tests.txt``, except in the test files they change and in those named
after a module they change (``analysis.py``: ``test_nuclear_analysis*.py``). Each entry of that list runs by itself ("part":
``file:quick`` is the rest of a file, ``file:ENTRY`` one entry), so that a long file is not one runner's.
**main** runs every test when any part's record no longer holds, and nothing when none is affected (a change of text).
**Tags and runs started by hand** run everything.

**What a record misses:** names looked up by a string that is put together, whether a file exists that was only looked for, new
versions of the packages of the requirements, and what a program outside the repository does. main's full runs and those by hand see those.

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
import urllib.parse
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests" / "python"
DURATIONS = ROOT / ".github" / "test-durations.json"
SLOW_TESTS = ROOT / ".github" / "slow-tests.txt"
UNKNOWN_S, START_S = 10.0, 2.0  # a file without a time; what a process takes beside its tests
SHARD_S, MAX_SHARDS = 75.0, 14
KEPT = 4  # records of a part
#: Programs that use nothing of this repository but what is given to them: those that libraries ask about the system, git (for
#: the commit in a record) and the browser that prints a report.
HARMLESS = {"fc-list", "ldconfig", "gcc", "cc", "ld", "objdump", "lscpu", "uname", "lsb_release", "file", "git", "google-chrome"}

RUST = re.compile(r"(src/|benches/|examples/|tests/[^/]+\.rs$|Cargo\.|\.cargo/|[^/]*\.toml$|rust-toolchain)")
MACHINERY = ("requirements-dev.txt", "scripts/ci_tests.py", "scripts/ci_record.py")
PACKAGE, CONFTEST = "python/physim/", "tests/python/conftest.py"
#: Decorators that do nothing but to what they decorate.
CALM = {"dataclass", "dataclasses.dataclass", "lru_cache", "functools.lru_cache", "cache", "functools.cache", "contextmanager",
        "contextlib.contextmanager", "runtime_checkable", "staticmethod", "classmethod", "property", "cached_property",
        "functools.cached_property", "functools.total_ordering", "total_ordering"}


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


def _dotted(node):
    """("a", "b", "c") of the expression a.b.c, or None."""
    out = []
    while isinstance(node, ast.Attribute):
        out.append(node.attr)
        node = node.value
    return (node.id, *reversed(out)) if isinstance(node, ast.Name) else None


def _decorator(node) -> str:
    return ".".join(_dotted(node.func if isinstance(node, ast.Call) else node) or ("?",))


FUNCTION = (ast.FunctionDef, ast.AsyncFunctionDef)


class Refs:
    """What a piece of code names: the dotted names it reads, and the imports inside it. Annotations don't count."""

    def __init__(self, *nodes):
        self.names, self.imports = set(), []
        for n in nodes:
            self._visit(n)

    def _visit(self, node):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            self.imports += _aliases(node)
            return
        if isinstance(node, ast.Attribute) and _dotted(node):
            self.names.add(_dotted(node))
            return
        if isinstance(node, ast.Name):
            self.names.add((node.id,))
        elif isinstance(node, ast.Constant) and isinstance(node.value, str) and re.fullmatch(r"physim(\.\w+)+", node.value):
            self.names.add(tuple(node.value.split(".")))  # monkeypatch.setattr("physim.nuclear.data.X", ...)
        for field, value in ast.iter_fields(node):
            if field in ("annotation", "returns"):
                continue
            for child in value if isinstance(value, list) else [value]:
                if isinstance(child, ast.AST):
                    self._visit(child)


def _aliases(node) -> list:
    """[(local name, module with its leading dots, name in it or None)] of an import."""
    if isinstance(node, ast.Import):
        return [(a.asname, a.name, None) if a.asname else (a.name.partition(".")[0], a.name.partition(".")[0], None)
                for a in node.names]
    module = "." * node.level + (node.module or "")
    return [(a.asname or a.name, module, a.name) for a in node.names]


class Code:
    """A Python source in pieces, each with its hash and :class:`Refs`: "units" are the functions and methods; "names" what
    the module binds in plain statements (constants, imports, classes without their methods; a plain function's is "def",
    as its unit has the rest); "rest" the statements that are none of these. "skeleton" is all but the units' bodies."""

    def __init__(self, text: str, docstrings: bool = False):
        tree = ast.parse(text)
        if not docstrings:
            _without_docstrings(tree)
        self.whole, self.refs = digest(ast.dump(tree)), Refs(tree)
        self.units, self.spans, self.names = {}, [], {}
        bound = self._bound(tree)
        top = [x for x in tree.body if isinstance(x, FUNCTION)]
        self._take(tree, "")  # which leaves no function its body
        self.skeleton = digest(ast.dump(tree))
        self.star = [a.module and "." * a.level + a.module or "." * a.level for a in tree.body
                     if isinstance(a, ast.ImportFrom) and any(n.name == "*" for n in a.names)]
        rest, pieces = [], {}
        for node in tree.body:
            got = self._plain(node, bound) if node not in top or self._calm(node) else None
            if got is None:
                rest.append(node)
                continue
            for name, value, refs in got:
                pieces.setdefault(name, []).append((value, refs))
        for name, parts in pieces.items():
            value = parts[0][0] if len(parts) == 1 else digest(" ".join(v for v, _ in parts))
            self.names[name] = (value, [r for _, r in parts])
        self.rest = digest(" ".join(ast.dump(x) for x in rest))
        self.rest_refs = Refs(*rest)
        self.imports = [a for node in tree.body for a in (_aliases(node) if isinstance(node, (ast.Import, ast.ImportFrom))
                                                          else [])] + self.rest_refs.imports
        self.bound = set(bound) | {a[0] for a in self.imports}

    @staticmethod
    def _bound(tree) -> dict:
        out = {}

        def visit(node):
            names = []
            if isinstance(node, (*FUNCTION, ast.ClassDef)):
                names = [node.name]
            elif isinstance(node, ast.Name) and not isinstance(node.ctx, ast.Load):
                names = [node.id]
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                names = [a[0] for a in _aliases(node)]
            for n in names:
                out[n] = out.get(n, 0) + 1
            for child in ast.iter_child_nodes(node):
                if not isinstance(node, (*FUNCTION, ast.ClassDef, ast.Lambda)):
                    visit(child)

        visit(tree)
        return out

    @staticmethod
    def _calm(node) -> bool:
        return all(_decorator(d) in CALM for d in node.decorator_list)

    def _plain(self, node, bound):
        """[(name, value, refs)] of a statement that only binds names, else None."""
        if isinstance(node, FUNCTION):
            outside = [*node.decorator_list, *node.args.defaults, *filter(None, node.args.kw_defaults)]
            if bound[node.name] > 1 or hasattr(builtins, node.name) or node.name.startswith("__"):
                return None
            if all(isinstance(n, (ast.Constant, ast.Tuple, ast.Load, ast.UnaryOp, ast.USub)) for d in outside
                   for n in ast.walk(d)):
                return [(node.name, "def", Refs())]
            return [(node.name, digest(ast.dump(node)), Refs(*outside))]
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            if any(a.name == "*" for a in node.names):
                return None
            return [(a[0], digest(repr(a)), Refs()) for a in _aliases(node)]
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if not all(isinstance(t, ast.Name) for t in targets):
                return None
            refs = Refs(*filter(None, [node.value]))
            return [(t.id, digest(ast.dump(node)), refs) for t in targets]
        if isinstance(node, ast.ClassDef):
            if not self._calm(node) or node.keywords:
                return None
            methods = []
            if not node.bases and not any(isinstance(x, FUNCTION) and x.name == "__getattr__" for x in node.body):
                names = [x.name for x in node.body if isinstance(x, FUNCTION)]
                methods = [x for x in node.body if isinstance(x, FUNCTION) and names.count(x.name) == 1
                           and not (x.name.startswith("_") and x.name.endswith("_"))
                           and all(_decorator(d) in CALM for d in x.decorator_list)]
            outside = Refs(node)
            node.body = [x for x in node.body if x not in methods]
            return [(node.name, " ".join([digest(ast.dump(node)), *sorted(m.name for m in methods)]), outside)]
        return None

    def _take(self, node, prefix):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, FUNCTION):
                name = key = prefix + child.name
                n = 1
                while key in self.units:  # a property's setter, a function defined twice
                    n += 1
                    key = f"{name}#{n}"
                self.units[key] = (digest(ast.dump(child)), Refs(child))
                first = min([child.lineno] + [d.lineno for d in child.decorator_list])
                self.spans.append((first, child.end_lineno, key))
                child.body = [ast.Pass()]
            elif isinstance(child, ast.ClassDef):
                self._take(child, f"{prefix}{child.name}.")
            elif isinstance(child, (ast.stmt, ast.ExceptHandler, ast.match_case)):
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
        if kind in ("unit", "name"):
            path, _, what = arg.partition("::")
            code = self.code(path)
            pieces = {} if code is None else code.units if kind == "unit" else code.names
            return pieces[what][0] if what in pieces else "-"
        if kind in ("rest", "skeleton"):
            code = self.code(arg)
            return getattr(code, kind) if code else self.bytes_hash(arg)
        if kind == "whole":
            return self.meaning(arg)
        if kind == "file":
            return self.bytes_hash(arg)
        if kind == "list":
            prefix = "" if arg == "." else arg + "/"
            return digest("\n".join(sorted({p[len(prefix):].partition("/")[0] for p in self.files
                                            if p.startswith(prefix)})))
        if kind == "package":
            return self._all(p for p in self.files if p.startswith(PACKAGE))
        if kind == "imports":
            return digest("\n".join(self.broken_imports()))
        if kind == "rust":
            return self._all(p for p in self.files if RUST.match(p))
        if kind == "machinery":
            return digest(f"{sys.version_info[:2]} {self.hash('rust:')} " + self._all(p for p in self.files if p in MACHINERY))
        if kind == "slow":
            return digest("\n".join(sorted(self.slow(arg))))
        raise ValueError(key)

    def tests(self) -> list:
        return sorted(PurePosixPath(p).name for p in self.files if re.fullmatch(r"tests/python/test_[^/]*\.py", p))

    def slow(self, file: str) -> list:
        """The entries of the list of slow tests for a test file."""
        text = self.read(".github/slow-tests.txt").decode() if ".github/slow-tests.txt" in self._set else ""
        return [x.strip().partition("::")[2] for x in text.splitlines() if x.startswith(f"{file}::")]

    # What names mean, from module to module of physim.

    def module(self, dotted: str, where: str = ""):
        """The path of a module of physim, by its name in an import in the file ``where``; None if it is not one."""
        name = dotted.lstrip(".")
        dots = len(dotted) - len(name)
        if dots:
            package = PurePosixPath(where).parent
            for _ in range(dots - 1):
                package = package.parent
            base = str(package) + ("/" + name.replace(".", "/") if name else "")
        else:
            base = "python/" + name.replace(".", "/")
        if not (base + "/").startswith(PACKAGE):
            return None
        return next((p for p in (base + "/__init__.py", base + ".py") if p in self._set), None)

    def _submodule(self, path: str, name: str):
        base = path.removesuffix("/__init__.py")
        return next((p for p in (f"{base}/{name}/__init__.py", f"{base}/{name}.py") if p in self._set),
                    None) if base != path else None

    def meant(self, path: str, chain: tuple, imports=(), depth: int = 0) -> list:
        """The keys of what a dotted name in a file stands for, and [(path, name)] of those to follow."""
        code = self.code(path)
        if code is None or depth > 8:
            return []
        head, out = chain[0], []
        local = next((a for a in imports if a[0] == head), None)
        found = local or next((a for a in code.imports if a[0] == head), None)
        if head in code.names and not local:
            out.append((path, head))
        if found:
            target = self.module(found[1], path)
            if target and found[2] is not None:
                inner = self.code(target)
                child = self._submodule(target, found[2])
                if child and not (inner and found[2] in inner.names and inner.names[found[2]][0] == "def"):
                    out += self._in_module(child, chain[1:], depth)
                else:
                    out += self.meant(target, (found[2], *chain[1:]), (), depth + 1)
            elif target:
                out += self._in_module(target, chain[1:], depth)
        elif head not in code.names:
            for star in code.star:
                target = self.module(star, path)
                if target:
                    out += self.meant(target, chain, (), depth + 1)
            if path.startswith("tests/") or head == "physim":
                target = self.module(head)
                if target and not out and len(chain) > 1:
                    out += self._in_module(target, chain[1:], depth)
        return out

    def _in_module(self, path: str, chain: tuple, depth: int) -> list:
        if not chain:
            return [(path, None)]  # the module itself, handed around
        code = self.code(path)
        child = self._submodule(path, chain[0])
        if child and not (code and code.names.get(chain[0], ("",))[0] == "def"):
            return self._in_module(child, chain[1:], depth + 1)
        return self.meant(path, chain, (), depth + 1)

    def follow(self, starts) -> set:
        """The keys of all that pieces of code name, and what that names: ``starts`` are (path, Refs)."""
        keys, seen, queue = set(), set(), []

        def add(path, refs):
            for chain in refs.names:
                queue.extend(self.meant(path, chain, refs.imports))
            for local, module, name in refs.imports:
                if name:  # what it takes from a module; of a module it takes, what its name is used for
                    queue.extend(x for x in self.meant(path, (local,), [(local, module, name)]) if x[1] is not None)

        for path, refs in starts:
            add(path, refs)
        while queue:
            at = queue.pop()
            if at in seen:
                continue
            seen.add(at)
            path, name = at
            if name is None:
                keys.add(f"skeleton:{path}")
                continue
            keys.add(f"name:{path}::{name}")
            for refs in self.code(path).names.get(name, ("", []))[1]:
                add(path, refs)
        return keys

    @functools.lru_cache(maxsize=None)
    def broken_imports(self) -> tuple:
        """The imports between modules of physim that name what isn't there (an import fails for all its users)."""
        out = []
        for path in self.files:
            if not (path.startswith(PACKAGE) and path.endswith(".py")):
                continue
            code = self.code(path)
            if code is None:
                out.append(f"{path} cannot be read")
                continue
            every = code.imports + [a for _, refs in code.units.values() for a in refs.imports]
            for _, module, name in every:
                if not (module.startswith(".") or module.partition(".")[0] == "physim"):
                    continue
                target = self.module(module, path)
                if target is None:
                    missing = not module.rpartition(".")[2].startswith("_core")
                else:
                    inner = self.code(target)
                    missing = (name not in (None, "*") and inner is not None and name not in inner.bound and not inner.star
                               and not self._submodule(target, name) and "__getattr__" not in inner.units
                               and not name.startswith("_core"))
                if missing:
                    out.append(f"{path}: {module} {name}")
        return tuple(sorted(set(out)))


def _git(*args) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True).stdout.decode("utf-8")


def used_by(seen: dict, file: str, mode: str, tree: Tree) -> dict:
    """{key: hash} of what ``ci_record`` saw a part of a test file use."""
    own = f"tests/python/{file}"
    keys = {"machinery:", "imports:", f"whole:{own}", f"slow:{file}"}
    starts = [(own, tree.code(own).refs)] if tree.code(own) else []
    for path, line, name in seen["code"]:
        code = tree.code(path) if path != own else None
        if code is None:
            continue  # the file itself, or not of this repository
        unit = code.unit_at(line) if name != "<module>" else None
        if unit:
            keys.add(f"unit:{path}::{unit}")
        if not path.startswith(PACKAGE):
            keys.add(f"skeleton:{path}")
            continue
        keys.add(f"rest:{path}")
        starts.append((path, code.rest_refs))
        if unit:
            refs = code.units[unit][1]
            starts.append((path, refs))
            if "." in unit and unit.partition(".")[0] in code.names:  # the class that it is a method of
                outer = Refs()
                outer.names.add((unit.partition(".")[0],))
                starts.append((path, outer))
    keys |= tree.follow(starts)
    known = set(tree.files)
    keys.update(f"file:{p}" for p in seen["files"] if p in known and p != own)
    directories = {str(parent) for p in tree.files for parent in PurePosixPath(p).parents}
    keys.update(f"list:{p}" for p in seen["lists"] if p in directories)
    if set(seen["spawned"]) - HARMLESS:
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
    """The newest ``KEPT`` different records of each part."""
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
    """Whether a piece is what it was. A class ("hash method method …") may have got new methods."""
    if now == then:
        return True
    (a, *before), (b, *after) = then.split(), now.split()
    return key.startswith("name:") and a == b and set(before) <= set(after)


def holds(records, file: str, mode: str, tree: Tree) -> bool:
    """Whether a part passed on what it would use in this tree."""
    return any(r["file"] == file and r["mode"] == mode and all(same(k, h, tree.hash(k)) for k, h in r["used"].items())
               for r in records)


# -- selection ------------------------------------------------------------------------------------------------------


def changed_in_meaning(tree: Tree, base: Tree) -> list:
    """The paths that differ by more than text: all but comments, layout and docstrings of Python."""
    paths = [p for p in _git("diff", "--name-only", "-z", "--no-renames", base.rev).split("\0") if p]
    return [p for p in paths if tree.meaning(p) != base.meaning(p)]


def is_text(path: str) -> bool:
    return path.endswith((".md", ".rst", ".txt")) and not path.endswith("requirements-dev.txt") or path.startswith(
        ("docs/", "backlog/", "planning/", ".github/"))


def slow_in(changed, files) -> list:
    """The test files whose slow tests a pull request runs, of the paths it changes: those changed, and those named after a
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
    """[(seconds, [(file, part)])] of each runner: the longest first, to the one with the least so far."""
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


def known_times(records, tree: Tree) -> dict:
    """{(file, part): seconds}: of the newest record, or of ``.github/test-durations.json``."""
    known = json.loads(DURATIONS.read_text(encoding="utf-8"))
    out = {}
    for f, whole in known["full"].items():
        quick, slow = known["quick"].get(f, whole), tree.slow(f)
        out[f, "quick"] = quick
        out.update({(f, entry): (whole - quick) / len(slow) for entry in slow})
    out.update({(r["file"], r["mode"]): r["seconds"] for r in records if "seconds" in r})
    return out


def parts(tree: Tree, files) -> list:
    return [(f, part) for f in files for part in ("quick", *tree.slow(f))]


def token(part) -> str:
    """A part in one word, for a command line."""
    return f"{part[0]}:{urllib.parse.quote(part[1], safe=':[]_-.=+,')}"


def choose(event: str, tree: Tree, base, records) -> dict:
    """What a run of CI does: the parts of the test files, and whether Rust's checks run."""
    files = tree.tests()
    if event not in ("pull_request", "push"):
        return {"tests": parts(tree, files), "rust": True}
    rust = not holds(records, "(rust)", "quick", tree)
    every = parts(tree, files)
    if event == "pull_request":
        whole = set(slow_in(changed_in_meaning(tree, base), files)) if base else set()
        every = [x for x in every if x[1] == "quick" or x[0] in whole]
    tests = [x for x in every if not holds(records, *x, tree)]
    return {"tests": every if event == "push" and tests else tests, "rust": rust}


def select(args) -> int:
    tree = Tree()
    records = load(args.store)
    got = choose(args.event, tree, Tree(args.base) if args.base else None, records)
    shards = plan(got["tests"], known_times(records, tree))
    print(f"{len({f for f, _ in got['tests']})} of {len(tree.tests())} test files on {len(shards)} runners; "
          f"{len(records)} records", file=sys.stderr)
    for i, (seconds, tests) in enumerate(shards, 1):
        print(f"  {i}: {seconds:.0f} s: {' '.join(f if m == 'quick' else f'{f}::{m}' for f, m in tests)}",
              file=sys.stderr)
    print("shards=" + json.dumps([{"index": i, "tests": " ".join(map(token, tests))}
                                  for i, (_, tests) in enumerate(shards, 1)]))
    print(f"python={'true' if shards else 'false'}")
    print(f"rust={'true' if got['rust'] else 'false'}")
    return 0


def run(args) -> int:
    tree, records, failed = Tree(), [], []
    for x in args.tests:
        file, _, mode = x.partition(":")
        mode = urllib.parse.unquote(mode) or "quick"
        command = [sys.executable, "-m", "pytest", f"tests/python/{file}", "-p", "ci_record", "-q", "-rs", "--durations=5",
                   "--durations-min=5", *(["--quick"] if mode == "quick" else ["--slow-entry", mode])]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "used.json"
            path = os.pathsep.join(filter(None, [str(ROOT / "scripts"), os.environ.get("PYTHONPATH")]))
            print(f"=== {file} ({mode})", flush=True)
            began = time.perf_counter()
            code = subprocess.run(command, cwd=ROOT, env={**os.environ, "PYTHONPATH": path,
                                                          "PHYSIM_CI_RECORD": str(out)}).returncode
            seconds = time.perf_counter() - began
            if code not in (0, 5):  # 5: no test left, as all of them are slow
                failed.append(f"{file}::{mode}")
                continue
            seen = json.loads(out.read_text(encoding="utf-8"))
        copies = [p for p, _, _ in seen["code"] if "site-packages/physim/" in p]
        if copies:  # an installed copy, perhaps an old one out of the cache, instead of the code at hand
            print(f"{file} ran {copies[0]}")
            failed.append(f"{file}::{mode}")
            continue
        others = sorted(set(seen["spawned"]) - HARMLESS)
        if others:
            print(f"started: {' '.join(others)}")
        records.append({"file": file, "mode": mode, "seconds": round(seconds, 1), "used": used_by(seen, file, mode, tree)})
    save(args.records, records)
    print(f"{len(records)} of {len(args.tests)} parts of test files passed" + (f"; failed: {' '.join(failed)}" if failed else ""))
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
        save(args.records, [{"file": f"({args.name})", "mode": "quick", "used": {f"{args.name}:": Tree().hash(f"{args.name}:")}}])
    elif args.what == "merge":
        save(args.out, merged(load([*args.store, *args.runs])))
    else:
        changed = [p for p in changed_in_meaning(Tree(), Tree(args.base)) if not is_text(p)]
        print("true" if changed else "false")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
