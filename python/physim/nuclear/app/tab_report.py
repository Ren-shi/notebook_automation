"""The Report tab (backlog item 69): the report of the experiment, one document with the setup, the plan, every run
and its spectra, the gates, the Doppler correction, the analysis and every number explained
(:mod:`physim.nuclear.logbook`), shown on the page and saved as HTML, as PDF or as a zip with the data."""

from __future__ import annotations

import io
import zipfile
from functools import partial

from .. import logbook, paper

#: Monte Carlo events of the plan's figures in the report.
PLAN_EVENTS = 100_000


def _last_result(p):
    """The analysis of the current run: this session's, or the one kept in the run's folder. The report finds
    the other runs' analyses itself (:func:`physim.nuclear.logbook.run_analyses`)."""
    try:
        return p.run_analysis()
    except Exception:  # noqa: BLE001 - no analysis is not an error here
        return None


def experiment_block(ctx) -> None:
    from nicegui import run

    ui = ctx.ui
    p = ctx.P()
    state = ctx.state
    runs = p.runs()
    name = p.experiment_folder.name if p.folder else "experiment"
    ui.label(f"{p.name}: one report, with everything in it. The setup with its level schemes and warnings; the plan "
             f"(rates, beam time, kinematics, expected peaks, predicted spectra); the {len(runs)} run"
             f"{'s' * (len(runs) != 1)} with their spectra (each particle detector, the energy against ring with "
             "your gates drawn on it, each γ-ray detector as measured and Doppler-corrected to the top of the "
             "level scheme); what each gate does to the γ-ray spectrum; the Doppler correction and its "
             "parameters; the analysis of every run you analysed; every number explained; and the method. The "
             "analysis is kept with its run, so it is still here after the setup changes or the app "
             "closes.").classes("text-sm ps-muted")
    made = {"page": None}  # the report as last built here, for the downloads that follow a look at it

    def options() -> dict:
        return {"plan_events": PLAN_EVENTS, "journal": state["journal"], "width": state["width"],
                "validation": True}

    async def page() -> str:
        if made["page"] is None:
            note = ui.notification("Building the report…", spinner=True, timeout=None)
            try:
                made["page"] = await run.io_bound(partial(logbook.experiment_html, p.snapshot(), _last_result(p),
                                                          **options()))
            finally:
                note.dismiss()
        return made["page"]

    async def show() -> None:
        holder.clear()
        with holder:
            ui.spinner(size="md")
        text = await page()
        holder.clear()
        with holder:
            ui.html('<div class="ps-record">' + text.split("<body>", 1)[1].rsplit("</body>", 1)[0] + "</div>",
                    sanitize=False).classes("w-full ps-plate")

    async def save_html() -> None:
        ui.download.content((await page()).encode("utf-8"), f"{name}-report.html")

    async def save_pdf() -> None:
        text = await page()
        note = ui.notification("Printing the report to PDF…", spinner=True, timeout=None)
        try:
            data = await run.io_bound(logbook.html_to_pdf, text)
        finally:
            note.dismiss()
        if data is None:
            ui.notify("No browser was found to print the PDF (Edge, Chrome or Chromium). Save the report as HTML, "
                      "open it and print it to PDF: it has the same pages.", type="warning", timeout=12000)
            return
        ui.download.content(data, f"{name}-report.pdf")

    async def save_zip() -> None:
        note = ui.notification("Building the report and packing the data…", spinner=True, timeout=None)
        try:
            data = await run.io_bound(partial(logbook.experiment_zip, p.snapshot(), _last_result(p), figures=True,
                                              **options()))
        finally:
            note.dismiss()
        ui.download.content(data, f"{name}-report.zip")

    async def folder() -> None:
        def pack() -> bytes:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
                for f in sorted(p.experiment_folder.rglob("*")):
                    if f.is_file():
                        z.write(f, p.experiment_folder.name + "/" + f.relative_to(p.experiment_folder).as_posix())
            return buf.getvalue()

        data = await run.io_bound(pack)
        ui.download.content(data, f"{p.experiment_folder.name}.zip")

    def pick(key: str, value) -> None:
        state[key] = value
        if state["width"] not in paper.JOURNALS[state["journal"]].widths:
            state["width"] = "single"
        made["page"] = None
        ctx.report_view.refresh()

    with ui.row().classes("gap-2 items-center"):
        ui.button("Show the report", icon="menu_book", on_click=show).props("unelevated no-caps")
        ui.button("Save as PDF", icon="picture_as_pdf", on_click=save_pdf).props("flat no-caps")
        ui.button("Save as HTML", icon="download", on_click=save_html).props("flat no-caps")
        ui.button("Save everything (zip)", icon="folder_zip", on_click=save_zip).props("flat no-caps").tooltip(
            "The report as HTML and PDF, the setup file, the plan's tables (CSV) and figures (PNG and PDF), and "
            "the experiment folder with every run")
        if p.folder is not None:
            ui.button("Save the experiment folder", icon="folder", on_click=folder).props("flat no-caps")
    with ui.expansion("The style of the plan's figures").classes("w-full text-sm"):
        ui.label("The setup, coverage, kinematics and predicted-spectra figures are drawn in a journal's style, at "
                 "its column width.").classes("text-sm ps-muted")
        with ui.row().classes("items-end gap-2"):
            ui.select({k: j.label for k, j in paper.JOURNALS.items()}, value=state["journal"], label="Journal style",
                      on_change=lambda e: pick("journal", e.value)).props("dense outlined").classes("w-80")
            ui.select({k: f"{paper.WIDTHS[k]} · {mm:g} mm"
                       for k, mm in paper.JOURNALS[state["journal"]].widths.items()}, value=state["width"],
                      label="Figure width", on_change=lambda e: pick("width", e.value)).props(
                "dense outlined").classes("w-56")
    holder = ui.column().classes("w-full")  # the report, once shown, below its buttons


def render(ctx) -> None:
    ctx.report_view()


def render_view(ctx) -> None:
    ctx.block("report", "The report")
