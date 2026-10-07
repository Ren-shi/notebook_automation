"""The Run tab: take a run of the setup for a duration, and the list of runs taken (backlog items 63 and 66)."""

from .. import runs as _runs


def render(ctx) -> None:
    ctx.run_panel()


def describe_runs(planner) -> list:
    """The run list as table rows."""
    rows = []
    for s in planner.runs():
        rows.append({"n": s["number"], "label": s["label"], "what": s["describe"],
                     "when": str(s.get("finished", ""))[:16].replace("T", " "),
                     "stale": "setup changed since" if s["stale"] else ""})
    return rows


__all__ = ["describe_runs", "render", "_runs"]
