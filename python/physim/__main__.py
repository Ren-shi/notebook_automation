"""``physim`` on the command line: ``physim app`` starts the experiment planner, ``physim report SETUP`` writes a
beam-time report. ``python -m physim ...`` does the same."""

from __future__ import annotations

import sys

USAGE = """usage: physim <command> [options]

commands:
  app       start the experiment planner in the browser (physim app --help)
  report    write a beam-time report for a setup file (physim report --help)
"""


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(USAGE)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "app":
        from .nuclear.app import main as app_main

        app_main(rest)
        return 0
    if cmd == "report":
        from .nuclear.report import main as report_main

        return report_main(rest)
    print(f"unknown command {cmd!r}\n\n{USAGE}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
