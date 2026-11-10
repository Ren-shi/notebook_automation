"""The Data tab (backlog item 67): every spectrum of the current run on screen, the operations as the user's
choices, gates as named objects, the 2D views, run overlays and exports (:mod:`physim.nuclear.dataviews`).

Loading a run's events and histogramming them takes seconds for a long run, so none of it is done while the tab
is built: :func:`prepare` computes every view for the current choices on the page's worker thread and keeps each
by what it depends on (``ctx.keep``), and the tab is then built from what is kept. A choice that changes one
kind of view (the Doppler correction, say) leaves the others as they are."""

from __future__ import annotations

import tempfile
from pathlib import Path

from .. import dataviews as dv
from .figures import figure_energy_ring, figure_gamma_crystal, figure_histogram, themed


def options(ctx) -> dict:
    return ctx.data_options


def _gate(ctx, p=None):
    name = options(ctx)["gate"]
    if not name:
        return None
    try:
        return (p or ctx.P()).gate(name)
    except KeyError:
        options(ctx)["gate"] = None
        return None


def _gate_key(gate) -> tuple:
    """What tells one gate from another (a gate can be saved again under its name)."""
    return None if gate is None else (gate.name, gate.describe())


def panels(ctx, p) -> list:
    """The grid for the current choices: (name, spectrum, its small figure) per particle detector and per γ-ray
    detector (or crystal), or per crystal of a source run. Each spectrum and figure is computed once per choice
    it depends on and kept with the run."""
    run = p.run
    o = options(ctx)
    theme = ctx.state["theme"]
    out = []

    def add(name: str, key: tuple, make, view=None) -> None:
        spec = ctx.keep(p, ("data", "spectrum") + key, make, of_run=True)
        if view is not None:
            # "γ from / to" only zooms the view: the histogram keeps the whole range of the run.
            spec = dict(spec, view=view)
        fig = ctx.keep(p, ("data", "small figure", theme, view) + key,
                       lambda: ctx.in_theme(figure_histogram([spec], height=190)), of_run=True)
        out.append((name, spec, fig))

    if run.kind == "source":
        sr = run.source()
        for n in sr.names():
            add(n, ("source", n), lambda n=n: dict(counts=sr.spectra[n], edges=sr.edges, label=n, scale=1e3,
                                                   xlabel="γ-ray energy (keV)"))
        return out
    ev = run.events()
    gate = _gate(ctx, p)
    bins, scaled = int(o["bins"]), o["scaled"]
    for d in ev.detectors:
        add(d, ("particle", d, _gate_key(gate), bins, scaled),
            lambda d=d: dict(dv.particle_spectrum(run, d, gate=gate, bins=bins, scaled=scaled),
                             xlabel="energy (MeV)"))
    g = run.gammas()
    if g is not None:
        correction, mode, randoms, addback = o["correction"], o["mode"], o["randoms"], o["addback"]
        view = (float(o["lo_kev"]), float(o["hi_kev"])) if o["lo_kev"] is not None and o["hi_kev"] else None
        for n in (g.crystal_names if o["crystals"] else g.detector_names()):
            # The histogram always holds the whole range of the run (every line, up to the highest energy
            # recorded), so nothing is lost from the panel or its CSV.
            add(n, ("gamma", n, correction, _gate_key(gate), mode, randoms, addback, bins, scaled),
                lambda n=n: dict(dv.gamma_spectrum(run, n, correction=correction, gate=gate, mode=mode,
                                                   randoms=randoms, addback=addback, bins=bins,
                                                   range=dv.gamma_range(g), scaled=scaled),
                                 scale=1e3, xlabel="γ-ray energy (keV)"), view)
    return out


def ring_view(ctx, p) -> tuple:
    """(detector, figure) of the energy-against-ring view chosen."""
    run = p.run
    ev = run.events()
    det = options(ctx)["view_detector"] or ev.detectors[0]
    view = ctx.keep(p, ("data", "ring", det), lambda: dv.energy_vs_ring(run, det), of_run=True)
    return det, ctx.keep(p, ("data", "ring figure", ctx.state["theme"], det),
                         lambda: ctx.in_theme(figure_energy_ring(view)), of_run=True)


def crystal_view(ctx, p):
    """The γ-ray-energy-against-crystal figure for the correction chosen."""
    run = p.run
    corr = options(ctx)["view_correction"]
    return ctx.keep(p, ("data", "crystal figure", ctx.state["theme"], corr),
                    lambda: ctx.in_theme(figure_gamma_crystal(dv.gamma_vs_crystal(run, corr))), of_run=True)


def compared(ctx, p):
    """The overlay of the current run and the one chosen under "Compare runs": {"figure", "differences"}, or
    {"error"}; None when no run is chosen. The other run's events are loaded for it."""
    run = p.run
    o = options(ctx)
    if o["other"] is None or run.kind == "source":
        return None
    ev = run.events()
    g = run.gammas()
    name = o["compare_spectrum"] or (list(ev.detectors) + (g.detector_names() if g is not None else []))[0]
    which = "particle" if name in ev.detectors else "gamma"
    opts = {"bins": int(o["bins"])} if which == "particle" else {
        "correction": o["correction"], "randoms": o["randoms"], "addback": o["addback"], "bins": int(o["bins"])}

    def make() -> dict:
        try:
            ov = dv.overlay(run, p.folder.load_run(o["other"]), name, which, **opts)
        except (ValueError, KeyError) as err:
            return {"error": str(err)}
        specs = [dict(x, scale=1e3 if which == "gamma" else 1.0,
                      xlabel="γ-ray energy (keV)" if which == "gamma" else "energy (MeV)") for x in ov["spectra"]]
        return {"figure": ctx.in_theme(figure_histogram(specs, height=340)), "differences": ov["differences"]}

    return ctx.keep(p, ("data", "compare", ctx.state["theme"], o["other"], name) + tuple(sorted(opts.items())),
                    make, of_run=True)


def prepare(ctx, p) -> None:
    """Load the current run's events and compute every view of the tab for the current choices, on the page's
    worker thread (``p`` is a snapshot of the planner). What is kept already is not computed again."""
    run = p.run
    if run is None:
        return
    panels(ctx, p)
    if run.kind == "source":
        return
    ring_view(ctx, p)
    if run.gammas() is not None:
        crystal_view(ctx, p)
    compared(ctx, p)


def controls(ctx, run) -> None:
    ui = ctx.ui
    o = options(ctx)
    exp = run.experiment
    redraw = ctx.refresh_data
    with ui.row().classes("items-end gap-2 w-full"):
        if run.gammas() is not None:
            ui.select({k: k.capitalize() for k in dv.CORRECTIONS}, value=o["correction"], label="Doppler correction",
                      on_change=lambda e: (o.update(correction=e.value), redraw())).props(
                "dense outlined").classes("w-40")
            ui.select({"coincidence": "In coincidence", "singles": "Singles"}, value=o["mode"], label="γ rays",
                      on_change=lambda e: (o.update(mode=e.value), redraw())).props("dense outlined").classes("w-36")
            ui.select({k: k.capitalize() for k in dv.RANDOMS}, value=o["randoms"], label="Random coincidences",
                      on_change=lambda e: (o.update(randoms=e.value), redraw())).props(
                "dense outlined").classes("w-40")
            if any(gd.addback or gd.shield for gd in exp.gamma_detectors):
                ui.switch("Add-back and suppression", value=o["addback"],
                          on_change=lambda e: (o.update(addback=e.value), redraw())).props("dense")
            ui.switch("Per crystal", value=o["crystals"],
                      on_change=lambda e: (o.update(crystals=e.value), redraw())).props("dense")
        gates = {"": "No gate"} | {g.name: g.name for g in ctx.P().gates()}
        ui.select(gates, value=o["gate"] or "", label="Particle gate",
                  on_change=lambda e: (o.update(gate=e.value or None), redraw())).props(
            "dense outlined").classes("w-44")
        ui.number("Bins", value=o["bins"], min=20, max=4000, step=50,
                  on_change=lambda e: o.update(bins=int(e.value or 200))).props("dense outlined").classes(
            "w-24").on("blur", redraw)
        if run.gammas() is not None:
            ui.number("γ view from (keV)", value=o["lo_kev"], on_change=lambda e: o.update(lo_kev=e.value)).props(
                "dense outlined").classes("w-32").on("blur", redraw).tooltip(
                "Zooms the γ-ray panels; the spectra keep their whole range (empty: everything the run recorded)")
            ui.number("to (keV)", value=o["hi_kev"], on_change=lambda e: o.update(hi_kev=e.value)).props(
                "dense outlined").classes("w-24").on("blur", redraw)
        if run.scale > 1.000001:
            ui.switch("Whole run (scaled)", value=o["scaled"],
                      on_change=lambda e: (o.update(scaled=e.value), redraw())).props("dense").tooltip(
                "The real part scaled to the full duration, with the full run's Poisson fluctuations")
    if run.gammas() is not None:
        ui.label(dv.correction_reason(exp, o["correction"]) + " " + dv.RANDOMS[o["randoms"]]).classes(
            "text-xs ps-muted")
    g = _gate(ctx)
    if g is not None:
        ui.label(f"Gate {g.name}: {g.describe()}.").classes("text-xs ps-accent")


def grid(ctx, run) -> None:
    ui = ctx.ui
    o = options(ctx)
    with ui.element("div").classes("w-full").style(
            "display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 8px"):
        for name, spec, fig in panels(ctx, ctx.P()):
            with ui.column().classes("ps-plate gap-0").style("min-width: 0"):
                with ui.row().classes("w-full items-center no-wrap"):
                    ui.label(name).classes("text-sm font-medium")
                    ui.label(f"{spec['counts'].sum():,.0f} counts").classes("text-xs ps-muted")
                    ui.space()
                    ui.button(icon="open_in_full", on_click=lambda n=name, s=spec: enlarge(ctx, run, n, s)).props(
                        "dense flat round size=sm").tooltip("Enlarge, and download the spectrum")
                ui.plotly(fig).classes("w-full")
    ui.label(("The whole run: the real part scaled, with the full run's Poisson fluctuations. " if o["scaled"]
              else f"The real part: {run.describe().split(';', 1)[-1].strip()}. ")
             + "Counts per bin.").classes("text-xs ps-muted")


def enlarge(ctx, run, name: str, spec: dict) -> None:
    ui = ctx.ui
    with ui.dialog() as dialog, ui.card().classes("ps-card").style("width: min(1000px, 95vw); max-width: 95vw"):
        with ui.row().classes("w-full items-center"):
            ui.label(spec.get("label", name)).classes("text-base font-medium")
            ui.space()
            log = ui.switch("Log scale", value=False)
            ui.button("CSV", icon="download",
                      on_click=lambda: ui.download.content(dv.spectrum_csv(spec), f"run{run.number}-{name}.csv")).props(
                "flat no-caps")
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")
        box = ui.column().classes("w-full")

        def draw() -> None:
            box.clear()
            with box:
                ui.plotly(themed(figure_histogram([spec], height=460, log=log.value), ctx.state["theme"])).classes(
                    "w-full")

        log.on_value_change(lambda e: draw())
        draw()
    dialog.open()


def gates_block(ctx, run) -> None:
    ui = ctx.ui
    p = ctx.P()
    ui.label("Gates and conditions").classes("text-lg font-semibold mt-4")
    ui.label("Named conditions on the particles: kept with the experiment, used by the spectra above and by the "
             "analysis by name. Drag a box on the energy-against-ring view below to fill in a gate.").classes(
        "text-xs ps-muted")
    for g in p.gates():
        with ui.row().classes("items-center gap-2 no-wrap"):
            ui.label(g.name).classes("text-sm font-medium")
            ui.label(g.describe()).classes("text-sm ps-muted")
            ui.button(icon="delete", on_click=lambda n=g.name: (p.remove_gate(n), ctx.refresh_data())).props(
                "dense flat round size=sm color=negative")
    f = ctx.gate_form
    names = list(run.events().detectors) if run.kind != "source" else []
    with ui.row().classes("items-end gap-2 w-full"):
        ui.input("Name", value=f["name"], on_change=lambda e: f.update(name=e.value)).props("dense outlined").classes(
            "w-36")
        ui.select({"": "All particle detectors"} | {d: d for d in names}, value=f["detector"] or "",
                  label="Detector", on_change=lambda e: f.update(detector=e.value or None)).props(
            "dense outlined").classes("w-44")
        ui.input("Rings (from–to)", value=f["rings"], placeholder="3-12",
                 on_change=lambda e: f.update(rings=e.value)).props("dense outlined").classes("w-32")
        ui.select({"any": "Any group", "excited": "Excited group", "elastic": "Elastic group"}, value=f["group"],
                  label="Group", on_change=lambda e: f.update(group=e.value)).props("dense outlined").classes("w-36")
        ui.input("Energy (MeV, from–to)", value=f["energy"], placeholder="6.2-6.9",
                 on_change=lambda e: f.update(energy=e.value)).props("dense outlined").classes("w-40")
        ui.button("Save gate", icon="save", on_click=lambda: save_gate(ctx)).props("flat no-caps")


def _pair(text: str, cast):
    text = str(text or "").replace("–", "-").replace(",", "-").strip()
    if not text:
        return None
    parts = [x for x in text.split("-") if x.strip()]
    if len(parts) != 2:
        raise ValueError(f"{text!r} is not a range: give two numbers, as 3-12")
    return cast(parts[0]), cast(parts[1])


def save_gate(ctx) -> None:
    f = ctx.gate_form
    try:
        rings = _pair(f["rings"], int)
        energy = _pair(f["energy"], float)
        if not f["name"]:
            raise ValueError("give the gate a name")
        gate = dv.Gate(f["name"], [f["detector"]] if f["detector"] else None,
                       list(range(min(rings) - 1, max(rings))) if rings else None, f["group"],
                       list(energy) if energy else None)
    except ValueError as err:
        ctx.ui.notify(str(err), type="warning")
        return
    ctx.P().put_gate(gate)
    ctx.data_options["gate"] = gate.name
    ctx.refresh_data()


def views_2d(ctx, run) -> None:
    ui = ctx.ui
    o = options(ctx)
    ev = run.events()
    ui.label("2D views").classes("text-lg font-semibold mt-4")
    with ui.row().classes("items-end gap-2"):
        ui.select(list(ev.detectors), value=o["view_detector"] or ev.detectors[0], label="Energy against ring of",
                  on_change=lambda e: (o.update(view_detector=e.value), ctx.refresh_data())).props(
            "dense outlined").classes("w-48")
    det, figure = ring_view(ctx, ctx.P())
    plot = ui.plotly(figure).classes("w-full")

    def selected(e) -> None:
        r = (e.args or {}).get("range") or {}
        if "x" in r and "y" in r:
            f = ctx.gate_form
            f.update(detector=det, rings=f"{max(1, round(min(r['x'])))}-{round(max(r['x']))}",
                     energy=f"{min(r['y']):.4g}-{max(r['y']):.4g}", group="any",
                     name=f["name"] or f"{det} region")
            ctx.refresh_data()

    plot.on("plotly_selected", selected)
    ui.label("Each group of particles follows its kinematic line: the elastic one dotted, the excited one solid. "
             "Drag a box with the mouse to fill in a gate below the spectra.").classes("text-xs ps-muted")
    if run.gammas() is not None:
        corr = o["view_correction"]
        ui.select({k: k.capitalize() for k in dv.CORRECTIONS}, value=corr, label="γ-ray energy against crystal",
                  on_change=lambda e: (o.update(view_correction=e.value), ctx.refresh_data())).props(
            "dense outlined").classes("w-56 mt-2")
        ui.plotly(crystal_view(ctx, ctx.P())).classes("w-full")
        ui.label("With the right correction every crystal's line sits at the transition energy; a misplaced target "
                 "moves each crystal's line by its own amount.").classes("text-xs ps-muted")


def compare(ctx, run) -> None:
    ui = ctx.ui
    p = ctx.P()
    o = options(ctx)
    others = [r for r in p.runs() if r["number"] != run.number and r["kind"] != "source"]
    ui.label("Compare runs").classes("text-lg font-semibold mt-4")
    if not others:
        ui.label("Take another run to overlay it here.").classes("text-sm ps-muted")
        return
    ev = run.events()
    g = run.gammas()
    spectra = list(ev.detectors) + (g.detector_names() if g is not None else [])
    with ui.row().classes("items-end gap-2"):
        ui.select({r["number"]: f"{r['number']}: {r['describe']}" for r in others}, value=o["other"],
                  label="Overlay run", on_change=lambda e: (o.update(other=e.value), ctx.refresh_data())).props(
            "dense outlined").classes("w-96")
        ui.select(spectra, value=o["compare_spectrum"] or spectra[0], label="Spectrum",
                  on_change=lambda e: (o.update(compare_spectrum=e.value), ctx.refresh_data())).props(
            "dense outlined").classes("w-40")
    ov = compared(ctx, p)
    if ov is None:
        return
    if "error" in ov:
        ui.label(ov["error"]).classes("text-sm ps-warn")
        return
    ui.plotly(ov["figure"]).classes("w-full")
    ui.label("Setups differ in: " + ("; ".join(ov["differences"]) if ov["differences"] else "nothing but labels")
             + ".").classes("text-xs ps-muted")


def exports(ctx, run) -> None:
    ui = ctx.ui
    p = ctx.P()
    ui.label("Export").classes("text-lg font-semibold mt-4")

    def root_file() -> None:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                path = dv.write_root(run, Path(tmp) / f"run{run.number}.root", p.gates())
                data = path.read_bytes()
        except ImportError as err:
            ui.notify(f"ROOT export needs uproot: {err}", type="warning")
            return
        ui.download.content(data, f"run{run.number:02d}.root")

    with ui.row().classes("gap-2"):
        if run.kind != "source":
            ui.button("ROOT file (events, spectra, gates as cuts)", icon="download", on_click=root_file).props(
                "flat no-caps")
        ui.button("The run's events (.npz)", icon="download",
                  on_click=lambda: ui.download.content((run.folder / "events.npz").read_bytes(),
                                                       f"run{run.number:02d}-events.npz")).props("flat no-caps")
    ui.label("A spectrum's CSV is in its enlarged view. The .npz holds named NumPy arrays (see the docs).").classes(
        "text-xs ps-muted")


def render_view(ctx) -> None:
    p = ctx.P()
    run = p.run
    if run is None:
        ctx.ui.markdown("No run yet. Take one on the **Run** tab: its spectra appear here, all at once. Below, a "
                        "preview from a simulated sample of the current setup.")
        ctx.block("spectra", "Preview")
        return
    status = p.run_status()
    head = f"Run {run.number}: {run.describe()}"
    ctx.ui.label(head).classes("text-base font-medium")
    if status and status["stale"]:
        ctx.ui.label("Taken with another setup: " + "; ".join(status["changes"][:4]) + ". These are its data as "
                     "taken; take a new run to see the current setup.").classes("text-sm ps-warn")
    if run.kind == "source":
        grid(ctx, run)
        exports(ctx, run)
        return
    controls(ctx, run)
    grid(ctx, run)
    gates_block(ctx, run)
    views_2d(ctx, run)
    compare(ctx, run)
    exports(ctx, run)


def render(ctx) -> None:
    ctx.data_view()
