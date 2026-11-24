"""The MCP server: the session, and the tools registered on an ``MCPServer`` from the official ``mcp`` package.

``build_server`` makes a server that a client talks to; ``physim mcp`` runs it over stdin/stdout. The tools
themselves are plain methods (:mod:`physim.mcp.tools_setup`, :mod:`physim.mcp.tools_shellmodel`,
:mod:`physim.mcp.tools_lifetimes`, :mod:`physim.mcp.tools_results`); this module only
registers them, marks which ones write files, and turns the planner's exceptions into messages the model can act on.
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
from .tools_lifetimes import LifetimeTools
from .tools_results import ResultTools
from .tools_setup import SetupTools
from .tools_shellmodel import ShellModelTools

INSTRUCTIONS = """\
physim's experiment planner for nuclear-physics beam time (Rutherford scattering, Coulomb excitation and
fusion-evaporation from a table of exit channels).

One experiment is open at a time. Start with list_experiments and open_experiment, or list_examples and
new_experiment. The other tools act on the open experiment.

Quantities are text with a unit, as in the app and setup files: "30 MeV", "0.5 mg/cm2", "45 deg", "100 mm",
"24 h". Every valid edit is saved to the experiment's folder at once. Tools that change the setup return the
experiment's summary, with problems (the setup is invalid until they are fixed) and warnings (most important
first): read them before going on.

The planner app (physim app) shows the same experiment folders. It does not watch for changes: an experiment
changed here shows the changes in the app once it is opened there again.

Shell-model calculations: physim drives KSHELL (an outside program) for laptop-sized valence spaces, with no
truncation. shell_model_options lists the interactions for a nucleus, estimate_shell_model says how big a request
is, run_shell_model starts a background job (follow it with shell_model_status), shell_model_results gives the
levels, B(E2), B(M1), Q_s and mu, and import_shell_model reads results run elsewhere. Two rules of the user's:
a run may use at most half of the RAM free when it starts (allow_over_limit only when the user explicitly says
so), and shell-model values never replace ENSDF, derived or user values in a level scheme (apply_shell_model's
replace only when the user explicitly asks).

Lifetime measurements (a plunger, RDDS, for about 1 ps to 1 ns; Doppler-shift attenuation, DSAM, for 10 fs to 1 ps;
after Coulomb excitation or fusion-evaporation): plan_lifetime_measurement says which method suits each level, the
plunger distances, whether the shifted and unshifted peaks separate at each detector angle, F(tau) and the beam time
per distance; set_plunger sets the distances, the stopper or degrader foils and the deorientation (or DSAM);
lifetime_spectra predicts the spectra per distance; analyse_lifetimes gets tau and B(E2) back from the runs taken at
the distances and compares them with ENSDF, the setup's Coulex B(E2) and the shell model.

Results across experiments (these tools act on a nucleus and need no open experiment): every result is a row with
its nucleus, level or transition, method and evidence. list_nuclei gives the nuclei with their experiments and what
is still open; list_results the rows of one nucleus; snapshot_evaluation takes ENSDF's adopted values in as
evaluated rows; add_data_point adds a published value or an own measurement (a value with its uncertainty, never
raw spectra; a correction is a new point that names the one it corrects, and retract_data_point withdraws one:
nothing is deleted); import_data_points adds the points of a JSON or CSV file and says what was added, already
there or refused; compare_results sets every result for one observable of one level or transition beside a
reference, with the routes to it and what the value gives. Evidence is in two groups that are never mixed:
independent (measured, evaluated, theory, outside program) and projected (a simulated run analysed back: it shows
the precision a setup would reach and is compared only with the value its simulation was given, never with a
measurement as if it tested it).
"""

READS = ToolAnnotations(readOnlyHint=True, openWorldHint=False)
WRITES = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False, openWorldHint=False)
#: Tools that overwrite what is there (the setup as a whole), for clients that ask before such calls.
REPLACES = ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=True, openWorldHint=False)
_REPLACING = {"replace_setup"}
#: Tools that can overwrite or stop something (apply_shell_model with replace=true; cancelling a running job).
DESTRUCTIVE = ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=False, openWorldHint=False)
_DESTRUCTIVE = {"apply_shell_model", "cancel_shell_model"}


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
    for tools in (SetupTools(session), ShellModelTools(session), LifetimeTools(session), ResultTools(session)):
        for name, writes in tools.TOOLS:
            annotations = (REPLACES if name in _REPLACING else DESTRUCTIVE if name in _DESTRUCTIVE
                           else WRITES if writes else READS)
            server.add_tool(_refusing(getattr(tools, name)), name=name, annotations=annotations)
    server.session = session  # for tests and scripts: the state behind the tools
    return server
