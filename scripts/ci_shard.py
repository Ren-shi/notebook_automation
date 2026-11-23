"""The Python test files of one of several CI runners, shared out so that each takes about the same time.

Usage: python scripts/ci_shard.py SHARD SHARDS      (SHARD from 1)   prints the files, one per line
       python scripts/ci_shard.py --plan SHARDS                      prints what each runner gets

Each file goes whole to one runner (its module fixtures are slow to make), the longest first, to the runner with the
least so far. How long a file takes (seconds, one at a time on GitHub's Linux runner) is in
``.github/test-durations.json``; a file that is not there counts as 10 s. The times only balance the runners: every
test file is always on exactly one of them.
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DURATIONS = ROOT / ".github" / "test-durations.json"
UNKNOWN_S = 10.0


def plan(shards: int) -> list:
    """[(seconds, files)] of each runner."""
    known = json.loads(DURATIONS.read_text(encoding="utf-8"))
    files = sorted(p.name for p in (ROOT / "tests" / "python").glob("test_*.py"))
    out = [[0.0, []] for _ in range(shards)]
    for name in sorted(files, key=lambda f: (-known.get(f, UNKNOWN_S), f)):
        least = min(out, key=lambda b: b[0])
        least[0] += known.get(name, UNKNOWN_S)
        least[1].append(name)
    return [(s, sorted(f)) for s, f in out]


def main(args) -> int:
    if args[0] == "--plan":
        for i, (seconds, files) in enumerate(plan(int(args[1])), 1):
            print(f"{i}: {seconds:.0f} s, {len(files)} files")
        return 0
    shard, shards = int(args[0]), int(args[1])
    if not 1 <= shard <= shards:
        raise SystemExit(f"shard {shard} of {shards}")
    print("\n".join(f"tests/python/{f}" for f in plan(shards)[shard - 1][1]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
