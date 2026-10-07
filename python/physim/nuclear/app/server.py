"""Starting the planner's server: ``physim app``, the desktop launch, ports."""

from __future__ import annotations

import os
import time
from typing import Optional

# -- starting the server --------------------------------------------------------------------------------------------

#: What ``/physim-planner`` answers, so a second launch can recognise a planner already running on the port.
MARKER = "physim experiment planner"


def log_path():
    """Where ``physim app --desktop`` writes its log (there is no console to write to)."""
    from pathlib import Path

    base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), ".local", "state")
    return Path(base) / "physim" / "planner.log"


def planner_at(port: int, host: str = "127.0.0.1") -> bool:
    """Whether a physim planner is already serving at ``host:port``."""
    import json
    import urllib.request

    try:
        with urllib.request.urlopen(f"http://{host}:{port}/physim-planner", timeout=2) as r:
            return json.loads(r.read().decode()).get("app") == MARKER
    except (OSError, ValueError):
        return False


def port_free(port: int, host: str = "127.0.0.1") -> bool:
    import socket

    with socket.socket() as s:
        try:
            s.bind((host, port))
        except OSError:
            return False
    return True


def pick_port(port: int, host: str = "127.0.0.1") -> int:
    """``port`` if it is free, otherwise a free port chosen by the system."""
    import socket

    if port_free(port, host):
        return port
    with socket.socket() as s:
        s.bind((host, 0))
        return s.getsockname()[1]


def main(argv: Optional[list] = None) -> None:
    """``physim app``: start the planner in the browser."""
    import argparse

    ap = argparse.ArgumentParser(prog="physim app", description="Start the experiment planner in the browser.")
    ap.add_argument("--example", default=None,
                    help="start from this example as a template (default: the most recent experiment)")
    ap.add_argument("--port", type=int, default=8080)
    ap.add_argument("--host", default="127.0.0.1",
                    help="address to listen on (default: this computer only; 0.0.0.0 serves the local network)")
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    ap.add_argument("--desktop", action="store_true",
                    help="started from a shortcut: log to a file, reuse a planner that is already running, use "
                         "another port if this one is taken, and stop when the last browser tab has been closed")
    ap.add_argument("--idle-exit", type=float, default=None, metavar="SECONDS",
                    help="stop when no browser tab has been open for this long (default with --desktop: 60)")
    args = ap.parse_args(argv)
    idle_exit = args.idle_exit if args.idle_exit is not None else (60.0 if args.desktop else None)
    local = "127.0.0.1" if args.host in ("0.0.0.0", "::") else args.host

    if args.desktop:
        import sys

        path = log_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        # Started with pythonw.exe there is no console: sys.stdout is None, and printing would fail.
        log = open(path, "a", encoding="utf-8", buffering=1)  # noqa: SIM115 -- open for the life of the server
        sys.stdout = sys.stderr = log
        print(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')} physim app {' '.join(argv or [])}")
        if planner_at(args.port, local):
            print(f"a planner is already running on port {args.port}; opening it")
            if not args.no_browser:
                import webbrowser

                webbrowser.open(f"http://{local}:{args.port}/")
            return
        args.port = pick_port(args.port, args.host)

    try:
        import matplotlib

        matplotlib.use("Agg")  # the report draws its figures in a worker thread, without a screen
    except ImportError:
        pass
    try:
        from nicegui import app, ui
    except ImportError:
        raise SystemExit("The planner app needs NiceGUI: pip install physim-engine[app]") from None

    @ui.page("/")
    def index(example: Optional[str] = None, template: Optional[str] = None, experiment: Optional[str] = None,
              theme: Optional[str] = None):
        from .figures import THEMES
        from .page import build_page

        build_page(example=example or args.example, template=template, experiment=experiment,
                   theme=theme if theme in THEMES else None)

    @app.get("/physim-planner")
    def marker():
        from ... import __version__

        return {"app": MARKER, "version": __version__}

    if idle_exit is not None:
        from nicegui import Client

        idle = {"since": time.monotonic()}

        def check_idle():
            if any(c.has_socket_connection for c in list(Client.instances.values())):
                idle["since"] = time.monotonic()
            elif time.monotonic() - idle["since"] > idle_exit:
                print(f"no browser tab open for {idle_exit:g} s; stopping")
                app.shutdown()

        def start_clock():
            idle["since"] = time.monotonic()
            _start_idle_timer(check_idle, min(5.0, idle_exit / 4))

        app.on_startup(start_clock)

    ui.run(title="physim experiment planner", host=args.host, port=args.port, show=not args.no_browser,
           reload=False, show_welcome_message=True)


def _start_idle_timer(check, every: float) -> None:
    import asyncio

    async def loop():
        while True:
            await asyncio.sleep(every)
            check()

    asyncio.get_running_loop().create_task(loop())
