"""The fusion-evaporation screens: what the Setup and Plan steps show for ``[reaction] type = "fusion-evaporation"``,
in place of the particle-detector numbers of Coulomb excitation and elastic scattering.

The Setup step: the strip above the scene (σ_fus, the studied channel's residues, the γ-ray array, an ancillary),
the whole experiment beside it, and one table per group of detectors (the γ-ray array with its Doppler shifts, the
ancillaries, the beam monitor). The Plan step: the channel studied, which sets the beam time; the fusion; the
residues; the studied residue's band with its known B(E2); the channel selection; what the compound nucleus boils
off. The spectra: the residues' γ rays by channel, one channel alone, the evaporated particles, the neutrons' flight
times and the recoils. Every number comes from :class:`physim.nuclear.planner.Planner`; what physim does not
compute yet (the ball's, separator's and calorimeter's selection, γγ coincidences, lifetimes) is said so."""

from __future__ import annotations

import math

from .tab_plan import _fmt, _time

#: How the ancillary sections are named on the page.
ANCILLARY_NAMES = {"neutron_wall": "Neutron wall", "charged_particle_ball": "Charged-particle ball",
                   "recoil_separator": "Recoil separator", "gamma_calorimeter": "γ calorimeter"}


def is_fusion(p) -> bool:
    return p.experiment.fusion is not None


def _study(p):
    try:
        return p.plan().get("study")
    except Exception:  # noqa: BLE001 -- the page shows what it can
        return None


# -- the Setup step ---------------------------------------------------------------------------------------------------

def strip_cells(p, s: dict) -> list:
    """The cells of the strip above the scene for a fusion-evaporation setup: (title, value, unit, sub-line)."""
    cells = []
    try:
        f = p.fusion()
        cells.append(("Fusion cross section", _fmt(f["sigma_thick_mb"]), "mb", "over the target, Bass + Wong"))
    except Exception:  # noqa: BLE001
        pass
    study = _study(p)
    if study:
        cells.append(("Residues studied", _fmt(study["rate_per_s"]), "/s",
                      f"{study['channel']} → {study['residue']}, {100 * study['share']:.0f} % of the table"))
    if s.get("gamma_efficiency") is not None:
        cells.append(("γ-ray array, full energy", f"{100 * s['gamma_efficiency']:.2f}", "%",
                      f"at {1e3 * p.gamma_energy_mev():.0f} keV"))
    exp = p.experiment
    if exp.recoil_separator is not None:
        from .. import workbench as wb

        acc = wb.ancillary_numbers(p).get("recoil_separator")
        if acc:
            cells.append(("Inside the separator's acceptance", f"{100 * acc['all']:.0f}", "%",
                          f"of the residues, within ±{acc['half_angle_deg']:g}°"))
    elif exp.neutron_wall is not None:
        m = exp.neutron_wall.model()
        cells.append(("Neutron wall", f"{100 * m.geometric_efficiency():.1f}", "% of 4π",
                      f"{exp.neutron_wall.cells} cells, out to θ {m.theta_max():.0f}°"))
    elif exp.detectors:
        cells.append(("Beam monitor", _fmt(s["rate_per_s"]), "/s", "elastically scattered beam"))
    return cells


def experiment_card(ctx, s: dict) -> None:
    """Beside the scene, with nothing selected: the fusion in a few lines, the detectors, the advice."""
    ui, p = ctx.ui, ctx.P()
    exp = p.experiment
    ui.label("The whole experiment").classes("ps-section")
    ui.label(s["title"]).classes("text-base font-medium")
    try:
        f = p.fusion()
        res = p.residues()
        lo, hi = f["excitation_range_mev"]
        rows = sorted(res["rows"], key=lambda r: -r["cross_section_mb"])
        main = ", ".join(f"{r['channel']} → {r['residue']} ({100 * r['share']:.0f} %)" for r in rows[:3])
        cone = p.recoil_cone()
        ui.label(f"{f['beam']} + {f['target']} → {f['compound']}* at E* = {lo:.1f}–{hi:.1f} MeV, spins up to "
                 f"{f['max_l']} ħ. Mostly {main}. The residues leave at β = {cone['beta']:.4f}, 90 % of them within "
                 f"{cone['half_angle_deg']:.1f}° of the beam.").classes("text-sm")
    except Exception as err:  # noqa: BLE001 -- an unfinished setup
        ui.label(f"The fusion cannot be worked out yet: {err}").classes("text-sm ps-warn")
    with ui.row().classes("w-full justify-between no-wrap"):
        ui.label("γ-ray detectors").classes("text-sm ps-muted")
        ui.label(str(len(exp.gamma_detectors))).classes("text-sm ps-num")
    present = [name for key, name in ANCILLARY_NAMES.items() if getattr(exp, key) is not None]
    with ui.row().classes("w-full justify-between no-wrap gap-3"):
        ui.label("Ancillary detectors").classes("text-sm ps-muted")
        ui.label(", ".join(present) or "none").classes("text-sm text-right")
    if exp.detectors:
        with ui.row().classes("w-full justify-between no-wrap"):
            ui.label("Beam monitor").classes("text-sm ps-muted")
            ui.label(", ".join(d.name or f"D{i + 1}" for i, d in enumerate(exp.detectors))).classes("text-sm")
    for line in s["advice"]:
        with ui.row().classes("items-start no-wrap gap-2 mt-1"):
            ui.icon("lightbulb").classes("ps-accent mt-0.5")
            ui.label(line).classes("text-sm")
    ui.label("Click a detector in the scene for its numbers.").classes("text-xs ps-muted mt-1")


def geometry_tables(ctx) -> None:
    """Beside the scene: the γ-ray array with the Doppler shift of the studied line, then a table per ancillary,
    then the beam monitor."""
    ui, p = ctx.ui, ctx.P()
    exp = p.experiment
    cols = ctx.columns
    study = _study(p)
    line = study["line"]["energy_kev"] if study and study["line"] else None
    try:
        beta = p.recoil_cone()["beta"]
    except Exception:  # noqa: BLE001
        beta = None
    gammas = ctx.D().geometry()["gamma_detectors"]
    if gammas:
        rows = []
        for g in gammas:
            row = {"detector": g["name"], "theta": f"{g['theta']:.1f}", "distance": f"{g['distance_mm']:g}",
                   "half": f"{g['half_angle_deg']:.1f}"}
            if line and beta:
                row["shift"] = f"{line * beta * math.cos(math.radians(g['theta'])):+.2f}"
            rows.append(row)
        spec = [("detector", "γ-ray detector"), ("theta", "θ (deg)"), ("distance", "Distance (mm)"),
                ("half", "Half-angle (deg)")]
        if line and beta:
            spec.append(("shift", f"Doppler at {line:.0f} keV"))
        ui.table(columns=cols(spec), rows=rows).props("dense flat").classes("w-full")
    if exp.neutron_wall is not None:
        w, m = exp.neutron_wall, exp.neutron_wall.model()
        _facts(ctx, "Neutron wall", [("Cells", str(w.cells)), ("Distance", str(w.distance)),
                                     ("Cone", f"{w.theta_min or '0 deg'} to {m.theta_max():.0f}°"),
                                     ("Of 4π", f"{100 * m.geometric_efficiency():.1f} %")])
    if exp.charged_particle_ball is not None:
        b = exp.charged_particle_ball
        _facts(ctx, "Charged-particle ball", [("Elements", f"{b.elements} {b.material or 'CsI'}"),
                                              ("Radius", str(b.radius)),
                                              ("Threshold", str(b.threshold) if b.threshold else "—"),
                                              ("Selection", "not computed yet")])
    if exp.recoil_separator is not None:
        from .. import workbench as wb

        r = exp.recoil_separator
        acc = wb.ancillary_numbers(p).get("recoil_separator", {})
        ui.table(columns=cols((("residue", "Residue (channel)"), ("inside", "Inside the acceptance"),
                               ("rate", "Into it (1/s)"))),
                 rows=[{"residue": f"{x['residue']} ({x['channel']})",
                        "inside": f"{100 * acc.get('channels', {}).get(x['channel'], 0):.0f} %",
                        "rate": _fmt(x["rate_per_s"] * acc.get("channels", {}).get(x["channel"], 0))}
                       for x in p.residues()["rows"]]).props("dense flat").classes("w-full")
        ui.label(f"Recoil separator: {r.kind or 'gas-filled'}, ±{r.theta_acceptance} by angle; its transmission "
                 "(charge states, momentum) is not computed yet.").classes("text-xs ps-muted")
    if exp.gamma_calorimeter is not None:
        c = exp.gamma_calorimeter
        _facts(ctx, "γ calorimeter", [("Elements", f"{c.elements} {c.material or 'BGO'}"), ("Radius", str(c.radius)),
                                      ("Fold at least", str(c.fold_min or "—")), ("Selection", "not computed yet")])
    if exp.detectors:
        rows = [{"detector": d["name"], "omega": _fmt(d["solid_angle_msr"]),
                 "theta": f"{d['theta_range'][0]:.1f}–{d['theta_range'][1]:.1f}"}
                for d in ctx.D().geometry()["detectors"]]
        ui.table(columns=cols((("detector", "Beam monitor"), ("omega", "Ω (msr)"), ("theta", "θ (deg)"))),
                 rows=rows).props("dense flat").classes("w-full")


def _facts(ctx, title: str, rows: list) -> None:
    ui = ctx.ui
    with ui.column().classes("w-full ps-card p-2 gap-0"):
        ui.label(title).classes("text-sm font-semibold")
        with ui.grid(columns=len(rows)).classes("w-full gap-x-3"):
            for k, v in rows:
                with ui.column().classes("gap-0"):
                    ui.label(k).classes("text-xs ps-muted")
                    ui.label(v).classes("text-sm ps-num")


# -- the Plan step ----------------------------------------------------------------------------------------------------

def choose_channel(ctx, which: str) -> None:
    """Make ``which`` (a channel's particles) the one studied: set at once, as a result's own button does."""
    ctx.direct(lambda p: p.set("run", "channel", which))


def channel_chips(ctx, pl: dict) -> None:
    ui, p = ctx.ui, ctx.P()
    study = pl["study"]
    res = pl["residues"]
    with ui.row().classes("w-full items-center gap-2 ps-card p-2").style("flex-wrap: wrap"):
        ui.label("The channel you study").classes("text-sm ps-muted")
        chans = p.experiment.fusion.channels
        for k, (ch, row) in enumerate(zip(chans, res["rows"])):
            on = study is not None and k == study["index"]
            ui.button(f"{row['channel']} → {row['residue']}", on_click=lambda c=ch.particles: choose_channel(ctx, c)).props(
                "dense no-caps rounded " + ("unelevated" if on else "outline")).classes("text-xs")
        ui.space()
        wanted = pl["counts_wanted"]
        ui.label(f"Counted: {wanted or '—'} counts in the channel's strongest line" if wanted else
                 "Set the counts wanted in Run conditions").classes("text-xs ps-muted")


def headline(ctx, pl: dict) -> None:
    ui = ctx.ui
    study = pl["study"]
    cells = []
    if study:
        line = study["line"]
        cells.append(("Beam time needed", _time(study["beam_time_s"]),
                      (f"for {pl['counts_wanted']} counts in {study['residue']}'s {line['energy_kev']:.1f} keV line"
                       if line and pl["counts_wanted"] else "no level scheme, or no counts wanted"),
                      "" if pl["enough"] is None else ("ps-ok" if pl["enough"] else "ps-bad")))
    cells.append(("Beam time planned", _time(pl["beam_time_s"]),
                  "enough" if pl["enough"] else ("too short" if pl["enough"] is False else ""), ""))
    if study:
        cells.append(("Residues", f"{_fmt(study['rate_per_s'])} /s",
                      f"{study['residue']}, {_fmt(study['cross_section_mb'])} mb of "
                      f"{_fmt(pl['residues']['total_mb'])}", ""))
        if study["line"]:
            cells.append(("Its strongest line", f"{_fmt(study['line']['detected_per_s'])} /s",
                          f"{study['line']['energy_kev']:.1f} keV ({study['line']['label']}), all γ detectors", ""))
    if pl["detectors"]:
        d = max(pl["detectors"], key=lambda x: x["rate_per_s"])
        cells.append(("Beam monitor", f"{_fmt(d['rate_per_s'])} /s", f"{d['detector']}: elastic counts", ""))
    with ui.element("div").classes("ps-readouts"):
        for title, value, sub, colour in cells:
            with ui.element("div").classes("ps-readout"):
                ui.label(title).classes("text-xs ps-muted")
                ui.label(value).classes("ps-readout-value " + colour)
                ui.label(sub).classes("text-xs ps-muted")


def residues(ctx, pl: dict) -> None:
    """The residues table: the studied channel highlighted; a click on a row studies that one."""
    ui = ctx.ui
    res = pl["residues"]
    study = pl["study"]
    lo, hi = res["excitation_range_mev"]
    head = (f"{res['compound']}* at E* = {res['excitation_mev']:.2f} MeV in the middle of the target ({lo:.2f}–"
            f"{hi:.2f} MeV through it); it moves along the beam at β = {res['beta']:.5f}. Σσ = "
            f"{_fmt(res['total_mb'])} mb.")
    chans = ctx.P().experiment.fusion.channels
    with ctx.panel("Residues of the channel table", "click a row to study that channel"):
        ui.label(head).classes("text-sm")
        rows = []
        for k, r in enumerate(res["rows"]):
            rows.append({"k": k, "particles": chans[k].particles if k < len(chans) else r["channel"],
                         "channel": ("● " if study and k == study["index"] else "") + r["channel"],
                         "residue": r["residue"], "sigma": _fmt(r["cross_section_mb"]),
                         "share": f"{100 * r['share']:.1f}%", "rate": _fmt(r["rate_per_s"]),
                         "run": _fmt(r["residues_in_run"]),
                         "lines": "; ".join(f"{x['energy_kev']:.1f} ({_fmt(x['detected_per_s'])}/s)"
                                            for x in r["lines"]) or "no level scheme",
                         "time": _time(r["beam_time_s"]), "open": r["open"]})
        table = ui.table(columns=ctx.columns((("channel", "Channel"), ("residue", "Residue"), ("sigma", "σ (mb)"),
                                              ("share", "Share"), ("rate", "Residues (1/s)"), ("run", "In the run"),
                                              ("lines", "Strongest lines, keV (counts/s)"), ("time", "Beam time"))),
                         rows=rows, row_key="k").props("dense flat").classes("w-full")
        table.on("rowClick", lambda e: choose_channel(ctx, e.args[1]["particles"]))
        ui.label("Line rates: full-energy-peak counts per second in all γ-ray detectors together (no angular "
                 "distribution, Doppler shift or summing yet). Beam time: for the counts wanted in the channel's "
                 "strongest line.").classes("text-xs ps-muted mt-1")
        for note in res["notes"]:
            ui.label(note).classes("text-xs ps-warn")


def band(ctx, pl: dict) -> None:
    """The studied residue's levels, the rates of its lines, and the B(E2) its scheme gives (stage 1), beside the
    lifetime measurement still to come (stage 2)."""
    ui, p = ctx.ui, ctx.P()
    from .figures import figure_band

    with ui.row().classes("w-full items-start gap-3").style("flex-wrap: wrap"):
        with ui.element("div").style("flex: 2 1 520px; min-width: 0"):
            try:
                b = p.band()
            except Exception as err:  # noqa: BLE001
                b = {"available": False, "reason": str(err), "residue": "?"}
            with ctx.panel(f"{b['residue']}: what the run sees of it",
                           "level scheme, γ rates from the plan, known B(E2)"):
                if not b["available"]:
                    ui.label(b["reason"]).classes("text-sm")
                elif not b["rows"]:
                    ui.label("None of its lines is strong enough to show.").classes("text-sm")
                else:
                    with ui.row().classes("w-full items-start gap-3 no-wrap"):
                        with ui.element("div").style("flex: 0 0 200px"):
                            ui.plotly(ctx.in_theme(figure_band(b))).classes("w-full")
                        rows = [{"t": f"{r['initial']} → {r['final']}", "e": f"{r['energy_kev']:.1f}",
                                 "rate": _fmt(r["rate_per_s"]),
                                 "tau": "—" if r["tau_ps"] is None else _fmt(r["tau_ps"]),
                                 "b": ("unknown" if r["b_e2_wu"] is None else
                                       f"{r['b_e2_wu']:.0f}" + (f" ({r['b_e2_unc_wu']:.0f})"
                                                                 if r["b_e2_unc_wu"] else ""))}
                                for r in b["rows"]]
                        ui.table(columns=ctx.columns((("t", "Transition"), ("e", "Eγ (keV)"), ("rate", "Rate (1/s)"),
                                                      ("tau", "τ (ps)"), ("b", "B(E2)↓ (W.u.)"))),
                                 rows=rows).props("dense flat").classes("w-full")
                    ui.label("B(E2)↓ from the level's half-life, branching and conversion in the "
                             + ("setup's level scheme" if b["from_setup"] else "local ENSDF copy")
                             + ": a known value beside the rates, not something this run measures."
                             + (f" {b['reference']}" if b.get("reference") else "")).classes("text-xs ps-muted")
        with ui.element("div").style("flex: 1 1 300px; min-width: 0"):
            with ctx.panel("Measure the lifetimes", "stage 2: plunger (RDDS) or line shapes (DSAM)"):
                ui.label("Later: a plunger or a backed target as an ancillary on the Setup step. The plan would then "
                         "show the shifted and unshifted peaks against the plunger distance, or the line shapes, and "
                         "the beam time for a wanted precision on τ, hence on B(E2). The levels above with no known "
                         "lifetime would be the first to measure.").classes("text-sm")
                ui.button("Add a plunger…").props("dense flat no-caps disable").tooltip("Not in this version")


def selection(ctx) -> None:
    """The channel selection: the neutron wall's matrix (worked out when asked), and what the other ancillaries
    would add."""
    ui, p = ctx.ui, ctx.P()
    exp = p.experiment
    from .tab_plan import selection as wall_selection

    if exp.neutron_wall is not None:
        wall_selection(ctx)
    others = [ANCILLARY_NAMES[k] for k in ("charged_particle_ball", "recoil_separator", "gamma_calorimeter")
              if getattr(exp, k) is not None]
    if exp.neutron_wall is None or others:
        with ctx.panel("Channel selection", "what each gate keeps"):
            if exp.neutron_wall is None:
                ui.label("Add a neutron wall (Setup, Ancillary detectors) to see which channel each event looks "
                         "like.").classes("text-sm")
            if others:
                ui.label(f"{', '.join(others)}: in the setup and the scene; their part in the selection (protons and "
                         "α in the ball, the separator's transmission, the calorimeter's fold) is not computed yet. "
                         "The separator's angular acceptance is on the Setup step.").classes("text-sm ps-muted")


def evaporation(ctx) -> None:
    ui, p = ctx.ui, ctx.P()
    from .figures import figure_evaporation

    try:
        ev = p.evaporation()
    except Exception as err:  # noqa: BLE001
        ui.label(f"The evaporation could not be sampled: {err}").classes("text-sm ps-warn")
        return
    titles = ("Particle energies", "Particles per event", f"Entry states, {ev['channel'] or 'all channels'}",
              f"Residue recoil, {ev['channel'] or 'all channels'}")
    with ctx.panel("What the compound nucleus boils off",
                   "from the channel table's events: particles, entry states, recoil"):
        with ui.element("div").classes("w-full").style(
                "display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 8px"):
            for title, fig in zip(titles, figure_evaporation(ev)):
                with ui.column().classes("gap-0 min-w-0"):
                    ui.label(title).classes("text-xs ps-muted")
                    ui.plotly(ctx.in_theme(fig)).classes("w-full")
        means = ", ".join(f"⟨{k}⟩ {v['mean']:.2f}" for k, v in ev["multiplicity"].items())
        ui.label(f"Mean multiplicities {means}; 90 % of the residues within {ev['cone']['half_angle_deg']:.1f}° "
                 f"at β = {ev['cone']['beta']:.4f}. The particles come from the channel table (their energies from "
                 "an evaporation spectrum at the channel's temperature), not from a full evaporation "
                 "calculation.").classes("text-xs ps-muted")


def render_plan(ctx) -> None:
    """The Plan step's answers for fusion-evaporation (refreshed after every edit)."""
    p = ctx.P()
    try:
        pl = p.plan()
    except Exception as err:  # noqa: BLE001 - say why instead of an empty tab
        ctx.ui.label(f"The plan could not be computed: {err}").classes("ps-bad")
        return
    if not pl.get("residues"):
        ctx.ui.label("The channel table cannot be worked out yet: see the setup's problems.").classes("ps-warn")
        return
    from .tab_plan import fusion as fusion_panel

    channel_chips(ctx, pl)
    headline(ctx, pl)
    if pl.get("fusion"):
        fusion_panel(ctx, pl["fusion"])
    residues(ctx, pl)
    band(ctx, pl)
    selection(ctx)
    evaporation(ctx)


# -- the spectra ----------------------------------------------------------------------------------------------------

def spectra(ctx) -> None:
    """The predicted spectra of a fusion-evaporation setup, in tabs: the residues' γ rays by channel, one channel
    alone, the evaporated particles, the neutrons, the recoils; γγ is not computed yet."""
    ui, p = ctx.ui, ctx.P()
    from nicegui import run

    from .figures import figure_doppler, figure_fusion_gammas, figure_tof

    exp = p.experiment
    ui.label("The residues' spectra").classes("font-semibold mt-3")
    ui.label("Each event picks an exit channel by its cross section; the particles leave with an evaporation "
             "spectrum and the residue recoils from each; its cascade enters the known levels by spin, and the "
             "crystals record the γ rays with their response. The Doppler correction assumes the residue moves "
             "along the beam at the compound nucleus's speed.").classes("text-xs ps-muted")
    holder = ui.column().classes("w-full")

    async def build() -> None:
        holder.clear()
        with holder:
            ui.spinner(size="md")
        try:
            g = await run.io_bound(p.snapshot().residue_gamma_spectra, 100_000, 1)
            p.residue_gamma_spectra(100_000, 1)  # the page's planner keeps them too
        except Exception as err:  # noqa: BLE001
            holder.clear()
            with holder:
                ui.label(str(err)).classes("text-sm ps-warn")
            return
        holder.clear()
        with holder:
            show(g)

    def show(g: dict) -> None:
        with ui.tabs().props("dense no-caps align=left").classes("w-full") as tabs:
            t_single = ui.tab("singles", "γ singles")
            t_one = ui.tab("one", "One channel alone")
            t_gg = ui.tab("gg", "γγ")
            t_part = ui.tab("particles", "Particles")
            t_n = ui.tab("neutrons", "Neutrons")
            t_rec = ui.tab("recoils", "Recoils")
        with ui.tab_panels(tabs, value=t_single).classes("w-full"):
            with ui.tab_panel(t_single):
                if not g["available"]:
                    ui.label(g["reason"]).classes("text-sm")
                else:
                    opts = {"det": "all", "corr": True, "col": True, "log": True}
                    box = ui.column().classes("w-full")

                    def draw() -> None:
                        box.clear()
                        with box:
                            ui.plotly(ctx.in_theme(figure_fusion_gammas(p, opts["det"], opts["corr"], opts["col"],
                                                                        log=opts["log"]))).classes("w-full")

                    with ui.row().classes("items-center gap-2"):
                        ui.select({"all": "all γ-ray detectors", **{n: n for n in g["gamma_detectors"]}}, value="all",
                                  label="Detector", on_change=lambda e: (opts.update(det=e.value), draw())).props(
                            "dense outlined").classes("w-48")
                        ui.toggle({True: "Doppler corrected", False: "as measured"}, value=True,
                                  on_change=lambda e: (opts.update(corr=e.value), draw())).props("dense no-caps")
                        ui.toggle({True: "colour by channel", False: "total only"}, value=True,
                                  on_change=lambda e: (opts.update(col=e.value), draw())).props("dense no-caps")
                        ui.toggle({True: "log", False: "linear"}, value=True,
                                  on_change=lambda e: (opts.update(log=e.value), draw())).props("dense no-caps")
                    draw()
                    ui.label(" ".join(g["notes"])).classes("text-xs ps-muted")
            with ui.tab_panel(t_one):
                if g["available"]:
                    study = _study(p)
                    first = g["channels"][study["index"]] if study and study["index"] < len(g["channels"]) else \
                        g["channels"][0]
                    box1 = ui.column().classes("w-full")

                    def draw1(label: str) -> None:
                        box1.clear()
                        with box1:
                            ui.plotly(ctx.in_theme(figure_fusion_gammas(p, "all", True, only=label))).classes("w-full")

                    ui.select(g["channels"], value=first, label="Channel",
                              on_change=lambda e: draw1(e.value)).props("dense outlined").classes("w-40")
                    draw1(first)
                    ui.label("One channel alone: what a perfect gate on it would keep, beside all channels in grey. "
                             "A real gate keeps less of it and lets some of the others through: the neutron wall's "
                             "share is in the Plan's channel selection; the ball's, separator's and calorimeter's "
                             "are not computed yet.").classes("text-xs ps-muted")
            with ui.tab_panel(t_gg):
                ui.label("γγ coincidences are not computed yet: they need the cascade's coincidences through two "
                         "detectors, which a later version adds.").classes("text-sm ps-muted")
            ev = None
            try:
                ev = p.evaporation()
            except Exception:  # noqa: BLE001
                pass
            with ui.tab_panel(t_part):
                if ev:
                    from .figures import figure_evaporation

                    figs = figure_evaporation(ev)
                    with ui.row().classes("w-full gap-3").style("flex-wrap: wrap"):
                        for title, fig in (("Energies in the compound nucleus's frame", figs[0]),
                                           ("Particles per event", figs[1])):
                            with ui.column().classes("gap-0").style("flex: 1 1 360px; min-width: 0"):
                                ui.label(title).classes("text-xs ps-muted")
                                ui.plotly(ctx.in_theme(fig)).classes("w-full")
                    ui.label("The evaporated particles as they leave the compound nucleus, before any detector: the "
                             "charged-particle ball's response is not computed yet.").classes("text-xs ps-muted")
            with ui.tab_panel(t_n):
                if ev and "tof" in ev:
                    ui.plotly(ctx.in_theme(figure_tof(ev))).classes("w-full")
                    ui.label("Flight times of the evaporated neutrons to the middle of a wall cell, before the "
                             "cells' efficiency and pulse-shape discrimination.").classes("text-xs ps-muted")
                else:
                    ui.label("Add a neutron wall to see the neutrons' flight times.").classes("text-sm ps-muted")
            with ui.tab_panel(t_rec):
                if ev:
                    from .figures import figure_evaporation

                    with ui.row().classes("w-full gap-3").style("flex-wrap: wrap"):
                        with ui.column().classes("gap-0").style("flex: 1 1 360px; min-width: 0"):
                            ui.label(f"Residue recoil, {ev['channel'] or 'all channels'}").classes("text-xs ps-muted")
                            ui.plotly(ctx.in_theme(figure_evaporation(ev)[3])).classes("w-full")
                        study = _study(p)
                        if study and study["line"] and exp.gamma_detectors:
                            with ui.column().classes("gap-0").style("flex: 1 1 360px; min-width: 0"):
                                ui.label("Doppler shift at each γ-ray detector").classes("text-xs ps-muted")
                                ui.plotly(ctx.in_theme(figure_doppler(p, study["line"]["energy_kev"]))).classes(
                                    "w-full")
                    ui.label("The focal-plane mass spectrum needs the separator's transmission, which is not computed "
                             "yet.").classes("text-xs ps-muted")

    if not exp.gamma_detectors:
        ui.label("Add γ-ray detectors to see the residues' γ rays.").classes("text-sm ps-muted")
    with holder:
        ui.button("Simulate the residues' spectra", icon="play_arrow", on_click=build).props("dense flat no-caps")
