"""The Plan tab (backlog item 65): what one checks before asking for beam, on one page, answers first.

The answers: the counts wanted and the beam time they need; per particle detector its rate, dead time, share of the
excitations and safe rings; per γ-ray detector its efficiency at the transition and its coincidence rate; the
particle × γ coincidence rates. Each number has its ? (the run record's explanation). Below them, the checks, in the
order one makes them, as panels that open: the kinematics, the energy loss, the orbit and the safe distance, the
particle–γ correlation and Doppler shifts, the γ-ray efficiency curves, the rates per strip with the sweep, and
the spectra a simulated sample of the setup gives (what a run is predicted to record)."""

from __future__ import annotations

import math

#: The checks below the answers, in the order one makes them: (planner tab, title, only for Coulomb excitation).
CHECKS = (("kinematics", "Kinematics: where the particles go", False),
          ("energy_loss", "Energy loss in the target and the dead layers", False),
          ("trajectories", "The orbit and the safe distance", False),
          ("gamma", "The excitation, the particle–γ correlation and the Doppler shifts", True),
          ("efficiency", "γ-ray efficiency curves", True),
          ("rates", "Rates per strip, and the beam-time sweep", False),
          ("spectra", "Predicted spectra: what each detector will record", False))


def _fmt(x, digits: int = 3) -> str:
    if x is None:
        return "—"
    if isinstance(x, float) and not math.isfinite(x):
        return "never" if x > 0 else "—"
    if abs(x) >= 1000:
        return f"{x:,.0f}"
    return f"{x:.{digits}g}"


def _time(seconds) -> str:
    if seconds is None:
        return "—"
    if not math.isfinite(seconds):
        return "never"
    if seconds >= 86400:
        return f"{seconds / 86400:.3g} d"
    if seconds >= 3600:
        return f"{seconds / 3600:.3g} h"
    return f"{seconds / 60:.3g} min" if seconds >= 60 else f"{seconds:.3g} s"


def _help(ctx, key: str, detector=None, gamma=None) -> None:
    ctx.ui.button(icon="help_outline", on_click=lambda: ctx.open_explanation(key, detector, gamma)).props(
        "dense flat round size=xs").classes("ps-muted").tooltip("How this number is worked out")


def answers(ctx, pl: dict) -> None:
    ui = ctx.ui
    what = {"all": "counts", "excitations": "excitations", "coincidences": "particle–γ coincidences"}[pl["measured"]]
    with ui.element("div").classes("ps-readouts"):
        with ui.element("div").classes("ps-readout"):
            ui.label("Beam time needed").classes("text-xs ps-muted")
            with ui.row().classes("items-center gap-1 no-wrap"):
                colour = "" if pl["enough"] is None else (" ps-ok" if pl["enough"] else " ps-bad")
                ui.label(_time(pl["beam_time_needed_s"])).classes("ps-readout-value" + colour)
                if pl["limiting_detector"]:
                    _help(ctx, "beam_time", pl["limiting_detector"])
            sub = (f"for {pl['counts_wanted']} {what} in {pl['limiting_detector']}" if pl["limiting_detector"]
                   else "set the counts wanted in Run conditions")
            ui.label(sub).classes("text-xs ps-muted")
        with ui.element("div").classes("ps-readout"):
            ui.label("Beam time planned").classes("text-xs ps-muted")
            ui.label(_time(pl["beam_time_s"])).classes("ps-readout-value")
            ui.label("enough" if pl["enough"] else ("too short" if pl["enough"] is False else "")).classes(
                "text-xs " + ("ps-ok" if pl["enough"] else "ps-bad"))
        with ui.element("div").classes("ps-readout"):
            ui.label("Beam").classes("text-xs ps-muted")
            ui.label(f"{_fmt(pl['particles_per_second'])} /s").classes("ps-readout-value")
            ui.label(f"live fraction {pl['live_fraction']:.3f}" + (
                f" (dead time {1e6 * pl['dead_time_s']:.3g} µs per count)" if pl["dead_time_s"] else
                " (no dead time set)")).classes("text-xs ps-muted")

    excite = pl["measured"] != "all"
    # The three tables stand side by side where the window has the width, each as wide as it needs.
    tables = ui.row().classes("w-full items-start gap-3 mt-3").style("flex-wrap: wrap")
    with tables:
        first = ui.element("div").style("flex: 2 1 480px; min-width: 0")
    with first, ctx.panel("Particle detectors"), ui.element("div").classes("w-full").style("overflow-x: auto"):
        with ui.grid(columns=7 if excite else 5).classes("gap-x-6 gap-y-1 items-center ps-num text-sm").style(
                f"grid-template-columns: repeat({7 if excite else 5}, auto); justify-content: space-between"):
            heads = ["Detector", "θ (deg)", "Rate (1/s)", "Busy", "Safe rings"]
            if excite:
                heads[4:4] = ["Excitations", "With γ (1/s)"]
            for h in heads:
                ui.label(h).classes("text-xs ps-muted")
            for d in pl["detectors"]:
                ui.label(d["detector"]).classes("font-medium")
                ui.label(f"{d['theta_range'][0]:.0f}–{d['theta_range'][1]:.0f}")
                with ui.row().classes("items-center gap-0 no-wrap"):
                    ui.label(_fmt(d["rate_per_s"]))
                    _help(ctx, "rate", d["detector"])
                ui.label(f"{100 * d['dead_time_fraction']:.2g}%")
                if excite:
                    ui.label(f"{100 * d['excitation_share']:.0f}%" if d["excitation_share"] is not None else "—")
                    with ui.row().classes("items-center gap-0 no-wrap"):
                        ui.label(_fmt(d["coincidence_per_s"]))
                        if d["coincidence_per_s"] is not None:
                            _help(ctx, "coincidences", d["detector"])
                bad = d["unsafe_rings"]
                ui.label(f"{len(d['safe_rings'])} of {d['rings']}" + (
                    f" (unsafe: {', '.join(str(k + 1) for k in bad[:6])}{'…' if len(bad) > 6 else ''})"
                    if bad else "")).classes("ps-bad" if bad else "")

    other = pl.get("other_peak")
    if other:
        seen = [(d["detector"], d["other_coincidence_per_s"]) for d in pl["detectors"]
                if d.get("other_coincidence_per_s") is not None]
        with first:
            ui.label(f"{other['nuclide']} is excited too (its level scheme is in the setup). Its {other['label']} peak at "
                     f"{other['energy_kev']:.1f} keV, in coincidence: "
                     + "; ".join(f"{name} {_fmt(v)} /s" for name, v in seen) + ". The counts wanted and the beam "
                     "time above are for the nucleus of the reaction.").classes("text-sm ps-muted mt-2")
    if pl["gamma_detectors"]:
        with tables:
            second = ui.element("div").style("flex: 1 1 300px; min-width: 0")
            third = ui.element("div").style("flex: 1 1 340px; min-width: 0")
        with second, ctx.panel("γ-ray detectors"), ui.grid(columns=3).classes(
                "gap-x-6 gap-y-1 items-center ps-num text-sm").style(
                "grid-template-columns: repeat(3, auto); justify-content: space-between"):
            for h in ("Detector", f"Efficiency at {pl['gamma_detectors'][0]['energy_kev']:.0f} keV",
                      "Coincidences (1/s)"):
                ui.label(h).classes("text-xs ps-muted")
            for g in pl["gamma_detectors"]:
                ui.label(g["detector"]).classes("font-medium")
                with ui.row().classes("items-center gap-0 no-wrap"):
                    ui.label(f"{100 * g['efficiency']:.3g}%")
                    _help(ctx, "gamma_efficiency", None, g["detector"])
                ui.label(_fmt(g["coincidence_per_s"]))
        names = [g["detector"] for g in pl["gamma_detectors"]]
        with third, ctx.panel("Particle × γ coincidences per second"), ui.element("div").classes("w-full").style(
                "overflow-x: auto"):
            with ui.grid(columns=len(names) + 1).classes("gap-x-6 gap-y-1 ps-num text-sm").style(
                    f"grid-template-columns: repeat({len(names) + 1}, auto); justify-content: space-between"):
                ui.label("")
                for n in names:
                    ui.label(n).classes("text-xs ps-muted")
                for d, row in pl["coincidences"].items():
                    ui.label(d).classes("font-medium")
                    for n in names:
                        ui.label(_fmt(row.get(n, 0.0)))
        with third:
            ui.label("Each number includes the particle–γ angular correlation; dead time and random coincidences "
                     "are not taken off here (a run measures them).").classes("text-xs ps-muted mt-1")


def answers_monitor(ctx, pl: dict) -> None:
    """For fusion-evaporation: the beam monitor's row (its elastic rate and dead time)."""
    ui = ctx.ui
    with ctx.panel("Beam monitor", "the particle detectors count the elastically scattered beam"), ui.grid(
            columns=4).classes("gap-x-6 gap-y-1 items-center ps-num text-sm").style(
            "grid-template-columns: repeat(4, auto); justify-content: space-between"):
        for h in ("Detector", "θ (deg)", "Rate (1/s)", "Busy"):
            ui.label(h).classes("text-xs ps-muted")
        for d in pl["detectors"]:
            ui.label(d["detector"]).classes("font-medium")
            ui.label(f"{d['theta_range'][0]:.0f}–{d['theta_range'][1]:.0f}")
            ui.label(_fmt(d["rate_per_s"]))
            ui.label(f"{100 * d['dead_time_fraction']:.2g}%")


def residues(ctx, res: dict) -> None:
    """The fusion-evaporation residues (backlog item 75): the compound nucleus, then per exit channel its cross
    section, residues per second and in the run, the beam time for the counts wanted in its strongest line, and
    its strongest lines."""
    ui = ctx.ui
    lo, hi = res["excitation_range_mev"]
    head = (f"{res['compound']}* at E* = {res['excitation_mev']:.2f} MeV in the middle of the target ({lo:.2f}–"
            f"{hi:.2f} MeV through it), E_cm = {res['e_cm_mev']:.2f} MeV, Q = {res['q_value_mev']:.3f} MeV; it "
            f"moves along the beam at β = {res['beta']:.5f}. Σσ = {_fmt(res['total_mb'])} mb.")
    with ctx.panel("Residues of the channel table"), ui.element("div").classes("w-full").style("overflow-x: auto"):
        ui.label(head).classes("text-sm")
        cols = 8
        with ui.grid(columns=cols).classes("gap-x-6 gap-y-1 items-center ps-num text-sm").style(
                f"grid-template-columns: repeat({cols}, auto); justify-content: space-between"):
            heads = ["Channel", "Residue", "σ (mb)", "Share", "Residues (1/s)", "In the run", "Strongest lines",
                     "Beam time"]
            for h in heads:
                ui.label(h).classes("text-xs ps-muted")
            for r in res["rows"]:
                ui.label(r["channel"]).classes("font-medium" + ("" if r["open"] else " ps-bad"))
                ui.label(r["residue"])
                ui.label(_fmt(r["cross_section_mb"]))
                ui.label(f"{100 * r['share']:.1f}%")
                ui.label(_fmt(r["rate_per_s"]))
                ui.label(_fmt(r["residues_in_run"]))
                lines = "; ".join(f"{x['energy_kev']:.1f} keV ({x['label']}, {_fmt(x['detected_per_s'])} /s)"
                                  for x in r["lines"]) or "no level scheme"
                ui.label(lines).classes("text-xs")
                ui.label(_time(r["beam_time_s"]))
        ui.label("Line rates: full-energy-peak counts per second in all γ-ray detectors together, from each "
                 "detector's efficiency at the line's energy (no angular distribution, Doppler shift or summing). "
                 "Beam time: for the counts wanted in the channel's strongest line. A channel in red needs more "
                 "excitation energy than the compound nucleus has.").classes("text-xs ps-muted mt-1")
        for note in res["notes"]:
            ui.label(note).classes("text-xs ps-warn")


def fusion(ctx, f: dict) -> None:
    """The fusion of beam and target (backlog item 75, Stage 2): the barrier, σ_fus at the beam energy and through
    the target, the spins it brings, and the excitation function."""
    ui = ctx.ui
    lo, hi = f["excitation_range_mev"]
    text = (f"{f['beam']} + {f['target']} → {f['compound']}*. Fusion barrier (Bass) {f['barrier_cm_mev']:.1f} MeV "
            f"in the CM, a beam of {f['barrier_lab_mev']:.1f} MeV in the lab, at {f['barrier_radius_fm']:.2f} fm "
            f"(ħω = {f['hbar_omega_mev']:.2f} MeV). σ_fus = {_fmt(f['sigma_mb'])} mb at the beam energy, "
            f"{_fmt(f['sigma_thick_mb'])} mb averaged over the target (the beam leaves it with "
            f"{f['e_out_mev']:.1f} MeV); E* = {lo:.1f}–{hi:.1f} MeV through it. Spins up to {f['max_l']} ħ fuse, "
            f"⟨ℓ⟩ = {f['mean_l']:.1f} ħ (critical ℓ {f['critical_l']} ħ).")
    with ctx.panel("Fusion", "physim's own estimate (Bass barrier, Wong's transmission)"):
        ui.label(text).classes("text-sm")
        try:
            from .figures import figure_fusion

            ui.plotly(ctx.in_theme(figure_fusion(ctx.P()))).classes("w-full")
        except Exception as err:  # noqa: BLE001 - the numbers above still stand
            ui.label(f"The excitation function could not be drawn: {err}").classes("text-xs ps-warn")
        ui.label("The channel table's cross sections are the setup's own; 'Fill from physim' in the reaction's "
                 "card replaces them with physim's fusion and evaporation.").classes("text-xs ps-muted mt-1")
        for note in f.get("warnings", []):
            ui.label(note).classes("text-xs ps-warn")


def selection(ctx) -> None:
    """The channel selection of the neutron wall and the particle detectors (backlog item 75): worked out when asked
    (a few seconds), as two tables: per true channel the share of its events in each detected class, and per class
    what a gate on it contains."""
    ui = ctx.ui
    with ctx.panel("Channel selection", "what the neutron wall and the particle detectors see of each channel"):
        holder = ui.column().classes("w-full")

        async def build() -> None:
            from nicegui import run

            holder.clear()
            with holder:
                ui.spinner(size="md")
            try:
                r = await run.io_bound(ctx.P().snapshot().channel_selection)
            except ValueError as err:
                holder.clear()
                with holder:
                    ui.label(str(err)).classes("text-sm ps-warn")
                return
            holder.clear()
            with holder:
                n = r["neutron"]
                ui.label(f"{n['cells']} cells, {n['solid_angle_msr']:.0f} msr ({100 * n['geometric_efficiency']:.1f}% "
                         f"of 4π, out to θ = {n['theta_max_deg']:.0f}°); a {n['mean_energy_mev']:.2f} MeV neutron "
                         f"(the mean in the lab) is seen with {100 * n['intrinsic_at_mean']:.0f}% of those entering, "
                         f"{100 * n['efficiency_at_mean']:.1f}% of all; cross-talk {100 * n['crosstalk']:.1f}%. "
                         + "; ".join(f"{100 * v:.0f}% of the {k == 'p' and 'protons' or 'α particles'} reach a "
                                     f"particle detector ({100 * r['charged_counted'][k]:.0f}% above its threshold)"
                                     for k, v in r["charged_reached"].items() if v is not None)
                         + f"; {r['gamma']['multiplicity']:.1f} γ rays per event, each a false neutron with "
                           f"{100 * r['gamma']['false_neutron_per_gamma']:.2f}%.").classes("text-sm")
                cols = len(r["assigned"]) + 1
                ui.label("Share of each channel's events in each detected class").classes("text-xs ps-muted mt-2")
                with ui.element("div").classes("w-full").style("overflow-x: auto"):
                    with ui.grid(columns=cols).classes("gap-x-4 gap-y-1 ps-num text-sm").style(
                            f"grid-template-columns: repeat({cols}, auto)"):
                        ui.label("Channel").classes("text-xs ps-muted")
                        for a in r["assigned"]:
                            ui.label(a).classes("text-xs ps-muted")
                        for name, row in zip(r["true"], r["matrix"]):
                            ui.label(name).classes("font-medium")
                            for v in row:
                                ui.label(f"{100 * v:.1f}%" if v > 0 else "·")
                ui.label("What a gate on each class contains (shares of its events, for the table's cross "
                         "sections)").classes("text-xs ps-muted mt-2")
                for cls, c in r["contents"].items():
                    parts = ", ".join(f"{k} {100 * v:.0f}%" for k, v in sorted(c["from"].items(),
                                                                               key=lambda kv: -kv[1]))
                    ui.label(f"{cls}: {_fmt(c['cross_section_mb'])} mb — {parts}").classes("text-xs")
                for note in r["notes"]:
                    ui.label(note).classes("text-xs ps-muted")

        with holder:
            ui.button("Work out the channel selection", icon="play_arrow", on_click=build).props(
                "dense flat no-caps")


def render_answers(ctx) -> None:
    """The answers (refreshed after every edit)."""
    if ctx.P().experiment.fusion is not None:
        from .tab_fusion import render_plan

        render_plan(ctx)
        if ctx.P().experiment.detectors:
            try:
                answers_monitor(ctx, ctx.P().plan())
            except Exception:  # noqa: BLE001 - the monitor's table is not essential
                pass
        return
    try:
        pl = ctx.P().plan()
    except Exception as err:  # noqa: BLE001 - say why instead of an empty tab
        ctx.ui.label(f"The plan could not be computed: {err}").classes("ps-bad")
        return
    answers(ctx, pl)
    if pl.get("fusion"):
        fusion(ctx, pl["fusion"])
    if pl.get("residues"):
        residues(ctx, pl["residues"])
    if ctx.P().experiment.neutron_wall is not None:
        selection(ctx)


def render(ctx) -> None:
    p = ctx.P()
    ctx.readouts()
    ctx.plan_answers()
    ctx.ui.label("The checks").classes("text-lg font-semibold mt-4")
    ctx.ui.label("Open each in turn: what the kinematics, the target and the detectors do to the numbers above. "
                 "The same panels are in Physics, at the foot of the rail on the left.").classes("text-xs ps-muted")
    coulex = p.experiment.excitation is not None and bool(p.experiment.gamma_detectors)
    fusion = p.experiment.fusion is not None
    # For fusion-evaporation the orbit check is Rutherford's (not the point), and the kinematics and per-strip
    # rates are the beam monitor's: shown only with one.
    skip = {"trajectories"} | (set() if p.experiment.detectors else {"kinematics", "rates"}) if fusion else set()
    # A check is built when it is first opened (its figure is computed by then, or a spinner says it is on its
    # way); the ones left open stay open when the tab is rebuilt after an edit.
    opened = ctx.state.setdefault("checks_open", set())
    for tab, title, coulex_only in CHECKS:
        if (coulex_only and not coulex) or tab in skip:
            continue
        if fusion:
            title = {"kinematics": "Kinematics of the beam monitor", "rates": "Rates of the beam monitor",
                     "spectra": "Predicted spectra: the residues' γ rays, the particles, the neutrons, the recoils"
                     }.get(tab, title)
        check = ctx.ui.expansion(title, value=tab in opened).classes("w-full ps-card")
        check.on_value_change(lambda e, tab=tab: (opened.add if e.value else opened.discard)(tab))
        ctx.when_opened(check, (lambda tab=tab: ctx.efficiency_panel() if tab == "efficiency" else ctx.block(tab)),
                        tab)
