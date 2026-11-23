"""The Python tests of one of several CI runners, shared out so that each takes about the same time.

Usage: python scripts/ci_shard.py SHARD SHARDS [--quick [--slow-in FILES]]    (SHARD from 1) the arguments for pytest
       python scripts/ci_shard.py --plan SHARDS [--quick [--slow-in FILES]]   what each runner gets
       python scripts/ci_shard.py --changed < CHANGED                         FILES for a pull request

Each file goes whole to one runner (its module fixtures are slow to make), the longest first, to the runner with the
least so far. How long a file takes (seconds, one at a time on GitHub's Linux runner) is in
``.github/test-durations.json``: whole, and "quick", without the tests of ``.github/slow-tests.txt``; a file that is not
there counts as 10 s. The times only balance the runners: every test file is always on exactly one of them.

Pull requests are quick: they leave out the tests of ``.github/slow-tests.txt``, except in the test files they change
and in those named after a module they change (``analysis.py``: ``test_nuclear_analysis*.py``), which run whole. main,
tags and runs started by hand run everything.
"""

import json
import re
import sys
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests" / "python"
DURATIONS = ROOT / ".github" / "test-durations.json"
UNKNOWN_S = 10.0


def slow_in(changed) -> list:
    """The test files that run whole, of the paths a pull request changes: those changed, and those named after a
    changed module of physim."""
    files = sorted(p.name for p in TESTS.glob("test_*.py"))
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


def plan(shards: int, quick: bool = False, whole=()) -> list:
    """[(seconds, files)] of each runner."""
    known = json.loads(DURATIONS.read_text(encoding="utf-8"))
    times = dict(known["full"])
    if quick:
        times.update({f: s for f, s in known["quick"].items() if f not in whole})
    files = sorted(p.name for p in TESTS.glob("test_*.py"))
    out = [[0.0, []] for _ in range(shards)]
    for name in sorted(files, key=lambda f: (-times.get(f, UNKNOWN_S), f)):
        least = min(out, key=lambda b: b[0])
        least[0] += times.get(name, UNKNOWN_S)
        least[1].append(name)
    return [(s, sorted(f)) for s, f in out]


def main(args) -> int:
    if args == ["--changed"]:
        print(",".join(slow_in(sys.stdin.read().split())))
        return 0
    quick = "--quick" in args
    names = args[args.index("--slow-in") + 1:] if "--slow-in" in args else []  # nothing follows if no file applies
    whole = [x for x in ",".join(names[:1]).split(",") if x]
    if args[0] == "--plan":
        for i, (seconds, files) in enumerate(plan(int(args[1]), quick, whole), 1):
            print(f"{i}: {seconds:.0f} s, {len(files)} files")
        return 0
    shard, shards = int(args[0]), int(args[1])
    if not 1 <= shard <= shards:
        raise SystemExit(f"shard {shard} of {shards}")
    out = [f"tests/python/{f}" for f in plan(shards, quick, whole)[shard - 1][1]]
    if quick:
        out.append("--quick")
        if whole:
            out.append("--slow-in=" + ",".join(whole))
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
