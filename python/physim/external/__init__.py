"""Outside programs, run out of sight: no windows, no console, no focus taken from the user.

``runner`` holds what every program shares (the hidden desktop, the job object that stops a run and everything it
started, the timeout, one copy at a time where a program needs it). Each program has a thin module of its own that
knows how to feed it, start it and collect what it writes:

- ``gosia``: GOSIA, a console program (input on standard input, results in ``fort.22``);
- ``srim``: SRIM 2013's SR Module (console) and TRIM (a window program run with ``TRIMAUTO`` = 1);
- ``pace4``: PACE4 from LISE++ in its silent mode (``-s``);
- ``gemini``: GEMINI++ from LISE++, a window program with no silent mode, whose Execute button is pressed through
  Windows UI Automation by ``gemini_drive.ps1`` on the hidden desktop.

The modules return the programs' output files as text; reading the numbers out of them is left to the callers.
Nothing is written into a program's install folder: runs happen in scratch folders or copies.
"""

from .runner import HIDDEN_DESKTOP, ProgramFailed, Run, run

__all__ = ["HIDDEN_DESKTOP", "ProgramFailed", "Run", "run"]
