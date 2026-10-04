"""Execute notebooks and fail if any cell raises.

Usage: python scripts/run_notebooks.py [notebook ...]   (default: notebooks/*.ipynb)
Notebooks are executed in memory; the files on disk are not modified.
"""

import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parent.parent


def main(paths):
    paths = [Path(p) for p in paths] or sorted((ROOT / "notebooks").glob("*.ipynb"))
    failed = []
    for path in paths:
        nb = nbformat.read(path, as_version=4)
        start = time.perf_counter()
        try:
            NotebookClient(nb, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(path.parent)}}).execute()
        except CellExecutionError as e:
            failed.append(path)
            print(f"FAIL {path}\n{e}", file=sys.stderr)
            continue
        # Errors raised after a cell's code finishes (e.g. while the inline backend draws a
        # figure) are recorded as outputs without failing the cell, so look for them too.
        errors = [(i, o) for i, c in enumerate(nb.cells) if c.cell_type == "code"
                  for o in c.get("outputs", []) if o.output_type == "error"]
        if errors:
            failed.append(path)
            for i, o in errors:
                print(f"FAIL {path} cell {i}: {o.ename}: {o.evalue}", file=sys.stderr)
            continue
        print(f"ok   {path} ({time.perf_counter() - start:.1f}s)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
