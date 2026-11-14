"""The MCP server: the session, and the tools registered on an ``MCPServer`` from the official ``mcp`` package.

``build_server`` makes a server that a client talks to; ``physim mcp`` runs it over stdin/stdout. The tools
themselves are plain methods (:mod:`physim.mcp.tools_setup`); this module only registers them, marks which ones
write files, and turns the planner's exceptions into messages the model can act on.
"""

from __future__ import annotations

import functools
from pathlib import Path
from typing import Callable, Optional, Union

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp_types import ToolAnnotations

from .. import __version__
from ..nuclear.experiment import SetupError
from ..nuclear.runs import RunInProgress
from .session import Session, ToolProblem
from .tools_setup import SetupTools

INSTRUCTIONS = """\
physim's experiment planner for nuclear-physics beam time (Rutherford scattering and Coulomb excitation).

One experiment is open at a time. Start with list_experiments and open_experiment, or list_examples and
new_experiment. The other tools act on the open experiment.

Quantities are text with a unit, as in the app and setup files: "30 MeV", "0.5 mg/cm2", "45 deg", "100 mm",
"24 h". Every valid edit is saved to the experiment's folder at once. Tools that change the setup return the
experiment's summary, with problems (the setup is invalid until they are fixed) and warnings (most important
first): read them before going on.

The planner app (physim app) shows the same experiment folders. It does not watch for changes: an experiment
changed here shows the changes in the app once it is opened there again.
"""

READS = ToolAnnotations(readOnlyHint=True, openWorldHint=False)
WRITES = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False, openWorldHint=False)
#: Tools that overwrite what is there (the setup as a whole), for clients that ask before such calls.
REPLACES = ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=True, openWorldHint=False)
_REPLACING = {"replace_setup"}


def message(e: BaseException) -> str:
    """The text of a refusal for the model: the planner's own words, without Python's quoting."""
    if isinstance(e, SetupError):
        return "The setup is not valid: " + "; ".join(str(p) for p in e.problems)
    if isinstance(e, KeyError) and e.args:
        return str(e.args[0])
    if isinstance(e, IndexError) and str(e).startswith("there are"):
        return f"No such detector: {e}. Count detectors from 1, as in the summary."
    return str(e) or type(e).__name__


#: The exceptions that are a refusal (bad input, the wrong moment), not a fault in physim.
_REFUSALS = (ToolProblem, SetupError, RunInProgress, KeyError, IndexError, ValueError, LookupError)


def _refusing(fn: Callable) -> Callable:
    @functools.wraps(fn)
    def tool(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except _REFUSALS as e:
            raise ToolError(message(e)) from e

    return tool


def build_server(root: Union[str, Path, None] = None, session: Optional[Session] = None) -> MCPServer:
    """An MCP server on the experiments folder ``root`` (default ``$PHYSIM_EXPERIMENTS`` or
    ``~/.physim/experiments``), or on an existing ``session``."""
    session = session or Session(root)
    server = MCPServer("physim-planner", title="physim experiment planner", instructions=INSTRUCTIONS,
                       version=__version__)
    tools = SetupTools(session)
    for name, writes in SetupTools.TOOLS:
        annotations = REPLACES if name in _REPLACING else WRITES if writes else READS
        server.add_tool(_refusing(getattr(tools, name)), name=name, annotations=annotations)
    server.session = session  # for tests and scripts: the state behind the tools
    return server
