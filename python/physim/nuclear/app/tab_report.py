"""The Report tab (backlog item 69): the record of the experiment and its runs, then the beam-time report of the
current setup and the run record (:mod:`physim.nuclear.logbook`, :mod:`physim.nuclear.report`)."""

from __future__ import annotations

import io
import zipfile

from .. import logbook


def _last_result(p):
    try:
        a = p.analysis() if p.run is not None and p.run.kind != "source" and p._data_experiment().excitation \
            is not None else None
    except Exception:  # noqa: BLE001 - no analysis is not an error here
        return None
    return a.history[-1] if a is not None and a.history else None


def experiment_block(ctx) -> None:
    from nicegui import run

    ui = ctx.ui
    p = ctx.P()
    ui.label("The record of the experiment").classes("text-lg font-semibold")
    runs = p.runs()
    ui.label(f"{p.name}: the setup and how it changed between runs, the Plan, the {len(runs)} run"
             f"{'s' * (len(runs) != 1)} with their summaries, the gates, the analysis of the current run "
             "(the last one you ran on the Analysis tab), and every number of the run explained.").classes(
        "text-sm ps-muted")
    holder = ui.column().classes("w-full")

    async def show() -> None:
        holder.clear()
        with holder:
            ui.spinner(size="md")
        page = await run.io_bound(logbook.experiment_html, p, _last_result(p))
        holder.clear()
        with holder:
            ui.html('<div class="ps-record">' + page.split("<body>", 1)[1].rsplit("</body>", 1)[0] + "</div>",
                    sanitize=False).classes("w-full ps-plate")

    async def download() -> None:
        note = ui.notification("Building the record of the experiment…", spinner=True, timeout=None)
        try:
            data = await run.io_bound(logbook.experiment_zip, p, _last_result(p))
        finally:
            note.dismiss()
        ui.download.content(data, f"{p.experiment_folder.name if p.folder else 'experiment'}-record.zip")

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

    with ui.row().classes("gap-2"):
        ui.button("Show the record of the experiment", icon="menu_book", on_click=show).props("flat no-caps")
        ui.button("Download it (zip: record, run record, setup, beam-time report, every run)", icon="download",
                  on_click=download).props("flat no-caps")
        if p.folder is not None:
            ui.button("Download the experiment folder", icon="folder_zip", on_click=folder).props("flat no-caps")


def render(ctx) -> None:
    ctx.report_view()


def render_view(ctx) -> None:
    experiment_block(ctx)
    ctx.ui.separator().classes("my-3")
    ctx.block("report", "The beam-time report of the current setup")
