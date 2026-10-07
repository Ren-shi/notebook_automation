"""The Data tab (backlog item 67): every spectrum of the current run on screen, the operations as the user's
choices, gates as named objects, the 2D views, run overlays and exports (:mod:`physim.nuclear.dataviews`)."""

from __future__ import annotations

import tempfile
from pathlib import Path

from .. import dataviews as dv
from .figures import figure_energy_ring, figure_gamma_crystal, figure_histogram, themed

GAMMA_RANGE_FACTOR = (0.6, 1.3)


def options(ctx) -> dict:
    return ctx.data_options


def _gate(ctx):
    name = options(ctx)["gate"]
    if not name:
        return None
    try:
        return ctx.P().gate(name)
    except KeyError:
        options(ctx)["gate"] = None
        return None


def _gamma_spec(ctx, run, name: str) -> dict:
    o = options(ctx)
    g = run.gammas()
    e0 = g.energy_mev
    rng = (o["lo_kev"] * 1e-3, o["hi_kev"] * 1e-3) if o["lo_kev"] and o["hi_kev"] else (
        GAMMA_RANGE_FACTOR[0] * e0, GAMMA_RANGE_FACTOR[1] * e0)
    s = dv.gamma_spectrum(run, name, correction=o["correction"], gate=_gate(ctx), mode=o["mode"],
                          randoms=o["randoms"], addback=o["addback"], bins=int(o["bins"]), range=rng,
                          scaled=o["scaled"])
    return dict(s, scale=1e3, xlabel="γ-ray energy (keV)")


def _particle_spec(ctx, run, name: str) -> dict:
    o = options(ctx)
    return dict(dv.particle_spectrum(run, name, gate=_gate(ctx), bins=int(o["bins"]), scaled=o["scaled"]),
                xlabel="energy (MeV)")


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
            ui.number("γ from (keV)", value=o["lo_kev"], on_change=lambda e: o.update(lo_kev=e.value)).props(
                "dense outlined").classes("w-28").on("blur", redraw)
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
    theme = ctx.state["theme"]
    if run.kind == "source":
        sr = run.source()
        panels = [(n, dict(counts=sr.spectra[n], edges=sr.edges, label=n, scale=1e3, xlabel="γ-ray energy (keV)"))
                  for n in sr.names()]
    else:
        ev = run.events()
        panels = [(d, _particle_spec(ctx, run, d)) for d in ev.detectors]
        g = run.gammas()
        if g is not None:
            names = g.crystal_names if o["crystals"] else g.detector_names()
            panels += [(n, _gamma_spec(ctx, run, n)) for n in names]
    with ui.element("div").classes("w-full").style(
            "display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 8px"):
        for name, spec in panels:
            with ui.column().classes("ps-plate gap-0").style("min-width: 0"):
                with ui.row().classes("w-full items-center no-wrap"):
                    ui.label(name).classes("text-sm font-medium")
                    ui.label(f"{spec['counts'].sum():,.0f} counts").classes("text-xs ps-muted")
                    ui.space()
                    ui.button(icon="open_in_full", on_click=lambda n=name, s=spec: enlarge(ctx, run, n, s)).props(
                        "dense flat round size=sm").tooltip("Enlarge, and download the spectrum")
                ui.plotly(themed(figure_histogram([spec], height=190), theme)).classes("w-full")
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
    det = o["view_detector"] or ev.detectors[0]
    view = dv.energy_vs_ring(run, det)
    plot = ui.plotly(themed(figure_energy_ring(view), ctx.state["theme"])).classes("w-full")

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
        ui.plotly(themed(figure_gamma_crystal(dv.gamma_vs_crystal(run, corr)), ctx.state["theme"])).classes("w-full")
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
    if o["other"] is None:
        return
    other = p.folder.load_run(o["other"])
    name = o["compare_spectrum"] or spectra[0]
    which = "particle" if name in ev.detectors else "gamma"
    opts = {"bins": int(o["bins"])} if which == "particle" else {
        "correction": o["correction"], "randoms": o["randoms"], "addback": o["addback"], "bins": int(o["bins"])}
    try:
        ov = dv.overlay(run, other, name, which, **opts)
    except (ValueError, KeyError) as err:
        ui.label(str(err)).classes("text-sm ps-warn")
        return
    specs = [dict(s, scale=1e3 if which == "gamma" else 1.0,
                  xlabel="γ-ray energy (keV)" if which == "gamma" else "energy (MeV)") for s in ov["spectra"]]
    ui.plotly(themed(figure_histogram(specs, height=340), ctx.state["theme"])).classes("w-full")
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
