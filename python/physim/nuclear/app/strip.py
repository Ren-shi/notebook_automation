"""The status strip across the top of the main area, on every tab: experiment · beam · target · detectors · last
run · warnings. Clicking a part opens its card in the setup panel."""

from __future__ import annotations

from .. import workbench as wb


def render(ctx) -> None:
    ui = ctx.ui
    try:
        parts = wb.status_strip(ctx.P())
    except Exception:  # noqa: BLE001 - the strip must never break the page
        return
    with ui.row().classes("w-full items-stretch gap-0 ps-strip no-wrap").style("overflow-x: auto"):
        for part in parts:
            cell = ui.element("div").classes("ps-strip-cell")
            if part["key"] is not None:
                cell.classes("cursor-pointer").on("click", lambda k=part["key"]: ctx.open_card(k))
            with cell:
                ui.label(part["label"]).classes("text-xs ps-muted")
                colour = ""
                if part["label"] == "Warnings" and part["text"] != "no warnings":
                    colour = " ps-warn"
                if part["label"] == "Last run" and "changed" in part["text"]:
                    colour = " ps-warn"
                ui.label(part["text"]).classes("text-sm font-medium" + colour)
