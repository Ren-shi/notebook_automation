"""The Plan tab (backlog item 65): what one checks before asking for beam, on one page, answers first.

The answers: the counts wanted and the beam time they need; per particle detector its rate, dead time, share of the
excitations and safe rings; per γ-ray detector its efficiency at the transition and its coincidence rate; the
particle × γ coincidence rates. Each number has its ? (the run record's explanation). Below them, the checks, in the
order one makes them, as panels that open: the kinematics, the energy loss, the orbit and the safe distance, the
particle–γ correlation and Doppler shifts, the γ-ray efficiency curves, and the rates per strip with the sweep."""

from __future__ import annotations

import math

#: The checks below the answers, in the order one makes them: (planner tab, title, only for Coulomb excitation).
CHECKS = (("kinematics", "Kinematics: where the particles go", False),
          ("energy_loss", "Energy loss in the target and the dead layers", False),
          ("trajectories", "The orbit and the safe distance", False),
          ("gamma", "The excitation, the particle–γ correlation and the Doppler shifts", True),
          ("efficiency", "γ-ray efficiency curves", True),
          ("rates", "Rates per strip, and the beam-time sweep", False))


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

    ui.label("Particle detectors").classes("ps-section mt-3")
    excite = pl["measured"] != "all"
    with ui.element("div").classes("w-full").style("overflow-x: auto"):
        with ui.grid(columns=7 if excite else 5).classes("gap-x-4 gap-y-1 items-center ps-num text-sm"):
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

    if pl["gamma_detectors"]:
        ui.label("γ-ray detectors").classes("ps-section mt-3")
        with ui.grid(columns=3).classes("gap-x-4 gap-y-1 items-center ps-num text-sm"):
            for h in ("Detector", f"Efficiency at {pl['gamma_detectors'][0]['energy_kev']:.0f} keV",
                      "Coincidences (1/s)"):
                ui.label(h).classes("text-xs ps-muted")
            for g in pl["gamma_detectors"]:
                ui.label(g["detector"]).classes("font-medium")
                with ui.row().classes("items-center gap-0 no-wrap"):
                    ui.label(f"{100 * g['efficiency']:.3g}%")
                    _help(ctx, "gamma_efficiency", None, g["detector"])
                ui.label(_fmt(g["coincidence_per_s"]))
        ui.label("Particle × γ coincidences per second").classes("ps-section mt-3")
        names = [g["detector"] for g in pl["gamma_detectors"]]
        with ui.element("div").classes("w-full").style("overflow-x: auto"):
            with ui.grid(columns=len(names) + 1).classes("gap-x-4 gap-y-1 ps-num text-sm"):
                ui.label("")
                for n in names:
                    ui.label(n).classes("text-xs ps-muted")
                for d, row in pl["coincidences"].items():
                    ui.label(d).classes("font-medium")
                    for n in names:
                        ui.label(_fmt(row.get(n, 0.0)))
        ui.label("Each number includes the particle–γ angular correlation; dead time and random coincidences are "
                 "not taken off here (the Run tab measures them).").classes("text-xs ps-muted")


def render_answers(ctx) -> None:
    """The answers (refreshed after every edit)."""
    try:
        pl = ctx.P().plan()
    except Exception as err:  # noqa: BLE001 - say why instead of an empty tab
        ctx.ui.label(f"The plan could not be computed: {err}").classes("ps-bad")
        return
    answers(ctx, pl)


def render(ctx) -> None:
    p = ctx.P()
    ctx.readouts()
    ctx.plan_answers()
    ctx.ui.label("The checks").classes("text-lg font-semibold mt-4")
    ctx.ui.label("Open each in turn: what the kinematics, the target and the detectors do to the numbers above. "
                 "The same panels are on the Physics tab.").classes("text-xs ps-muted")
    coulex = p.experiment.excitation is not None and bool(p.experiment.gamma_detectors)
    for tab, title, coulex_only in CHECKS:
        if coulex_only and not coulex:
            continue
        with ctx.ui.expansion(title).classes("w-full ps-card"):
            if tab == "efficiency":
                ctx.efficiency_panel()
            else:
                ctx.block(tab)
