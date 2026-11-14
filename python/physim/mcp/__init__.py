"""``physim mcp``: the experiment planner as an MCP server, so an LLM client (Claude Code, Claude Desktop, ...)
can open experiments, change the setup and read the results, as a person does in the app.

The server runs over stdin/stdout; the client starts it::

    claude mcp add physim -- physim mcp

It works on the saved experiments under ``$PHYSIM_EXPERIMENTS`` or ``~/.physim/experiments`` (``--experiments``
picks another folder), the same folders the app opens. It needs the ``mcp`` extra
(``pip install "physim-engine[mcp]"``, Python 3.10 or later).

The tools are plain Python (:mod:`physim.mcp.tools_setup`, on a :class:`~physim.mcp.session.Session`) and can be
called without a client.
"""

from __future__ import annotations

import sys
from typing import Optional

USAGE_MISSING = (
    "physim mcp needs the MCP package: pip install \"physim-engine[mcp]\" (Python 3.10 or later)."
)


def main(argv: Optional[list] = None) -> int:
    """``physim mcp [--experiments PATH]``: run the MCP server over stdin/stdout."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="physim mcp",
        description="Run physim's experiment planner as an MCP server over stdin/stdout, for an LLM client. "
                    "Add it to Claude Code with: claude mcp add physim -- physim mcp")
    parser.add_argument("--experiments", metavar="PATH", default=None,
                        help="the experiments folder (default: $PHYSIM_EXPERIMENTS or ~/.physim/experiments)")
    args = parser.parse_args(argv)
    try:
        from .server import build_server
    except ImportError as e:
        if getattr(e, "name", "") and e.name.split(".")[0] in ("mcp", "mcp_types"):
            print(USAGE_MISSING, file=sys.stderr)
            return 2
        raise
    build_server(args.experiments).run()
    return 0
