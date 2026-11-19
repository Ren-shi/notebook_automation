"""The Shell model panel beside the level scheme of the target and of the beam.

For one nucleus: choose an interaction, the states and the multipoles; see what the calculation needs before
running it (dimension, memory, disk, the limit of half the RAM free now, a rough time); run KSHELL in the
background with its progress, or import results run elsewhere; compare the shell-model states with the known
levels; and put chosen values into the level scheme through :meth:`physim.nuclear.planner.Planner.apply_shell_model`.

Two rules of the user's hold here. Shell-model values never replace ENSDF, derived or the user's own values unless
the user ticks "Replace existing values" (off by default, with a warning). A run may use at most half of the RAM free
when it starts: over that the Run button stays disabled unless the user ticks an explicit override.

The panel is closed until opened, so the page is not slowed by it. What it holds (the choices, the job, the result)
lives in the page's state, so it comes back as it was when the level scheme is drawn again after an apply."""

from __future__ import annotations

import copy
import datetime
import re
import time
from pathlib import Path
from typing import Optional

from .. import shellmodel as sm
from .. import shellmodel_import as smi
from .. import shellmodel_jobs as jobs
from .. import shellmodel_spaces as spaces
from ..names import parse_nuclide

#: The interaction chosen first when several have the same size (the usual choice for each space first).
PREFERRED = ("ckpot", "usdb", "usda", "gxpf1a", "kb3g", "jun45", "sdpf-mu")

#: What can be applied, as the planner's ``which`` names them, with their labels on the page.
QUANTITIES = (("energy", "Energies"), ("E2", "E2 matrix elements"), ("Q", "Quadrupole moments"))

_STATE = re.compile(r"^\s*(\d+)(/2)?\s*([+-])\s*(?:[x×*]\s*(\d+))?\s*$", re.IGNORECASE)


# -- choices as text ---------------------------------------------------------------------------------------------------

def states_text(states) -> str:
    """``((0, 1, 3), (4, 1, 3))`` → ``"0+ ×3, 2+ ×3"``."""
    return ", ".join(f"{sm.label(j2, p)} ×{n}" for j2, p, n in states)


def parse_states(text: str) -> tuple:
    """``"0+ ×3, 2+ ×3, 3/2- x2"`` → ((2J, parity, how many), ...). ``×n`` may be left out (one state), written
    ``x n`` or ``*n``. ValueError, naming the item, for anything else."""
    out = []
    for item in (text or "").replace(";", ",").split(","):
        if not item.strip():
            continue
        m = _STATE.match(item)
        if not m:
            raise ValueError(f"'{item.strip()}' is not a state: write J and parity with how many, e.g. 2+ ×3 or "
                             f"3/2- ×2")
        j, half, sign, n = m.groups()
        j2 = int(j) if half else 2 * int(j)
        out.append((j2, 1 if sign == "+" else -1, int(n) if n else 1))
    if not out:
        raise ValueError("ask for at least one state, e.g. 0+ ×3, 2+ ×3")
    return tuple(out)


def default_option(options: list) -> Optional[spaces.Option]:
    """The interaction offered first: the smallest calculation, then the usual choice for its space."""
    if not options:
        return None
    rank = {name: i for i, name in enumerate(PREFERRED)}
    return min(options, key=lambda o: (max(o.dimension_by_parity.values(), default=0),
                                       rank.get(o.interaction, len(rank))))


def default_states(option: spaces.Option) -> tuple:
    """The states asked for first: for an even-even nucleus the lowest 3×0⁺, 3×2⁺ and 2×4⁺; for an odd-A one two
    each of 1/2 to 7/2; for an odd-odd one the lowest of J = 0 to 5; in the space's natural parity, and only the
    J that exist there (no more states than exist)."""
    z, n = option.z_val, option.n_val
    parity = option.parities[0] if option.parities else 1
    if (z + n) % 2:
        want = ((1, 2), (3, 2), (5, 2), (7, 2))
    elif z % 2 == 0:
        want = ((0, 3), (4, 3), (8, 2))
    else:
        want = tuple((j2, 1) for j2 in range(0, 12, 2))
    space = option.spec.space
    out = []
    for j2, count in want:
        exist = spaces.j_dimension(space, z, n, parity, j2)
        if exist:
            out.append((j2, parity, min(count, exist)))
    return tuple(out) or ((0 if (z + n) % 2 == 0 else 1, parity, 1),)


def which_of(choice: dict) -> list:
    """The planner's ``which`` from the panel's apply choices: the quantities ticked, the states chosen (none: all
    the matched ones), and ``"new"`` to add the states the scheme has no level for (a new level needs its energy,
    so ``"new"`` alone brings ``"energy"``). Nothing ticked: an empty list, nothing to apply."""
    which = [q for q, _ in QUANTITIES if choice.get(q, True)]
    if choice.get("new"):
        which = (which or ["energy"]) + ["new"]
    if not which:
        return []
    return which + list(choice.get("states") or [])


def preview(result: sm.ShellModelResult, scheme, which, replace: bool) -> list:
    """What applying would change, one line each (and what it would keep), worked out on a copy of the scheme:
    the experiment is not touched. Without a scheme, applying starts one from the result's states."""
    if scheme is None or not scheme.levels:
        return [f"The scheme has no levels yet: applying starts one from the {len(result.states)} shell-model "
                f"states, every value marked \"shell model\"."]
    lines = [plain(line) for line in sm.apply(result, copy.deepcopy(scheme), which, replace=replace)]
    if not any(not kept(line) for line in lines):
        lines.insert(0, "Nothing would change: the scheme holds these values already, or keeps its own.")
    return lines


def kept(line: str) -> bool:
    """Whether a line of :func:`physim.nuclear.shellmodel.apply` says what was kept rather than what changed."""
    return line.startswith("kept ")


def plain(line: str) -> str:
    """A line of :func:`physim.nuclear.shellmodel.apply` in the page's words."""
    return line.replace("apply with replace=True to overwrite them",
                        "tick \"Replace existing values\" to overwrite them")


def elapsed(st: dict) -> Optional[float]:
    """Seconds a job has run (to its end once finished), from its status's ISO times."""
    def when(key):
        try:
            return datetime.datetime.fromisoformat(str(st.get(key)).replace("Z", "+00:00"))
        except (TypeError, ValueError):
            return None
    start = when("started")
    if start is None:
        return None
    end = when("finished") or datetime.datetime.now(datetime.timezone.utc)
    if start.tzinfo is None:
        start = start.replace(tzinfo=datetime.timezone.utc)
    if end.tzinfo is None:
        end = end.replace(tzinfo=datetime.timezone.utc)
    return max(0.0, (end - start).total_seconds())


def _seconds(s: Optional[float]) -> str:
    if s is None:
        return "—"
    if s < 90:
        return f"{s:.0f} s"
    if s < 5400:
        return f"{s / 60:.0f} min"
    if s < 2 * 86400:
        return f"{s / 3600:.1f} h"
    return f"{s / 86400:.0f} days"


def _bytes(n: Optional[int]) -> str:
    if n is None:
        return "—"
    if n < 1024**2:
        return "under 1 MB"
    if n < 1024**3:
        return f"{n / 1024**2:.0f} MB"
    gb = n / 1024**3
    return f"{gb:.2g} GB" if gb < 10 else f"{gb:,.0f} GB".replace(",", " ")


def _int(n: int) -> str:
    return f"{n:,}".replace(",", " ")


def _install():
    from ...external import kshell as _kshell

    try:
        return _kshell.find_install()
    except Exception:  # noqa: BLE001 -- a broken setting reads as no KSHELL; submit says why
        return None


def cached_job(request: sm.ShellModelRequest, install) -> Optional[str]:
    """The job ID of ``request`` if its results are in the cache already (no run needed), else None."""
    if install is None:
        return None
    key = jobs.job_key(request, install.version)
    return key if (jobs.job_folder(key) / "results.json").is_file() else None


# -- the panel --------------------------------------------------------------------------------------------------------

def role_state(ctx, role: str, nuclide: str) -> dict:
    """The panel's state for one role, kept in the page's state; started again when the nuclide changes."""
    store = ctx.state.setdefault("shellmodel", {})
    st = store.get(role)
    if st is None or st.get("nuclide") != nuclide:
        st = store[role] = {"nuclide": nuclide, "open": False, "interaction": None, "states": None,
                            "kinds": ["E2", "M1"], "eff": None, "override": False, "job": None, "result": None,
                            "origin": None, "apply": {"energy": True, "E2": True, "Q": True, "new": False,
                                                      "states": []},
                            "replace": False, "applied": None, "import_path": ""}
    return st


def block(ctx, role: str, nuclide: str) -> None:
    """The Shell model panel of one role (``"target"`` or ``"beam"``) for ``nuclide``: closed until opened."""
    ui = ctx.ui
    st = role_state(ctx, role, nuclide)
    exp = ui.expansion(f"Shell model: calculate or import the levels of {nuclide}", icon="hub",
                       value=st["open"]).classes("w-full ps-plate").props("dense")
    exp.on_value_change(lambda e: st.update(open=bool(e.value)))
    ctx.when_opened(exp, lambda: _Panel(ctx, role, nuclide, st).build())


class _Panel:
    """The panel's content for one role, built when it is opened."""

    def __init__(self, ctx, role: str, nuclide: str, st: dict):
        self.ctx, self.ui, self.role, self.nuclide, self.st = ctx, ctx.ui, role, nuclide, st
        self.install = _install()
        try:
            self.options = spaces.options(nuclide)
        except ValueError:
            self.options = []
        self.by_name = {o.interaction: o for o in self.options}
        self.request: Optional[sm.ShellModelRequest] = None
        self.estimate: Optional[spaces.Estimate] = None
        self.timer = None

    # -- layout
    def build(self) -> None:
        ui = self.ui
        ui.label("KSHELL computes the levels, B(E2), B(M1), quadrupole and magnetic moments of the shell model; "
                 "physim runs it out of sight, or reads results run elsewhere. Applying puts chosen values into the "
                 "level scheme with the source \"shell model\": they fill gaps and never replace ENSDF, derived or "
                 "your own values unless you say so.").classes("text-xs ps-muted")
        with ui.row().classes("w-full items-start gap-6").style("flex-wrap: wrap"):
            with ui.column().classes("gap-1").style("flex: 2 1 420px; min-width: 0"):
                ui.label("Calculate").classes("ps-section")
                self.choose_box = ui.column().classes("w-full gap-1")
                self.estimate_box = ui.column().classes("w-full gap-1")
                self.run_box = ui.column().classes("w-full gap-1")
            with ui.column().classes("gap-1").style("flex: 1 1 300px; min-width: 0"):
                ui.label("Or import results run elsewhere").classes("ps-section")
                self.import_form()
        self.results_box = ui.column().classes("w-full gap-1")
        self.choose()
        self.timer = ui.timer(1.0, self.poll, active=False)
        self.refresh_choice()
        if self.st["job"] and self.st["result"] is None:
            self.poll()  # a job started before the panel was drawn again: follow it
        self.show_results()

    # -- 1. choose
    def choose(self) -> None:
        ui, st = self.ui, self.st
        with self.choose_box:
            if not self.options:
                ui.label(f"No interaction physim knows covers {self.nuclide}: its valence protons and neutrons are "
                         "outside the p, sd, fp, f5pg9 and sd-pf spaces (or it has none). Import results run "
                         "elsewhere instead.").classes("text-sm ps-warn")
                return
            if st["interaction"] not in self.by_name:
                first = default_option(self.options)
                st.update(interaction=first.interaction, states=states_text(default_states(first)),
                          eff=list(first.spec.eff_charge))
            labels = {}
            for o in self.options:
                dims = " / ".join(_int(d) for d in o.dimension_by_parity.values())
                labels[o.interaction] = (f"{o.spec.label}: {o.z_val} valence protons, {o.n_val} neutrons, "
                                         f"M-scheme dimension {dims}")
            with ui.row().classes("w-full items-center gap-2 no-wrap"):
                sel = ui.select(labels, value=st["interaction"], label="Interaction",
                                on_change=lambda e: self.pick(e.value)).props("dense outlined").style(
                    "flex: 1 1 auto; min-width: 0")
                spec = self.by_name[st["interaction"]].spec
                self.cite = ui.icon("info").classes("ps-muted").tooltip(spec.citation + (
                    f". {spec.notes}" if spec.notes else ""))
            self.select = sel
            self.states_in = ui.input("States (J± ×how many)", value=st["states"]).props(
                "dense outlined").classes("w-full").tooltip(
                "The lowest states of each J and parity to compute, e.g. 0+ ×3, 2+ ×3, 4+ ×2; half-integer J as "
                "3/2- ×2.")
            self.states_in.on("blur", lambda: self.set_states(self.states_in.value))
            self.states_in.on("keydown.enter", lambda: self.set_states(self.states_in.value))
            with ui.row().classes("items-center gap-3"):
                ui.label("Transitions").classes("text-sm ps-muted")
                for kind in ("E2", "M1"):
                    ui.checkbox(kind, value=kind in st["kinds"],
                                on_change=lambda e, k=kind: self.set_kind(k, e.value)).props("dense")
                ui.label("Effective charges").classes("text-sm ps-muted ms-3")
                self.ep = ui.number("e_p", value=st["eff"][0], step=0.05, format="%.2f",
                                    on_change=lambda e: self.set_eff(0, e.value)).props("dense outlined").classes("w-20")
                self.en = ui.number("e_n", value=st["eff"][1], step=0.05, format="%.2f",
                                    on_change=lambda e: self.set_eff(1, e.value)).props("dense outlined").classes("w-20")
                self.eff_note = ui.icon("help_outline").classes("ps-muted").tooltip(
                    f"{spec.name.upper()} defaults: {spec.defaults_from}.")

    def pick(self, name: str) -> None:
        o = self.by_name.get(name)
        if o is None:
            return
        self.st.update(interaction=name, states=states_text(default_states(o)), eff=list(o.spec.eff_charge),
                       override=False)
        self.states_in.value = self.st["states"]
        self.ep.value, self.en.value = self.st["eff"]
        self.cite.tooltip(o.spec.citation + (f". {o.spec.notes}" if o.spec.notes else ""))
        self.eff_note.tooltip(f"{o.spec.name.upper()} defaults: {o.spec.defaults_from}.")
        self.refresh_choice()

    def set_states(self, text: str) -> None:
        if (text or "") != self.st["states"]:
            self.st["states"], self.st["override"] = text or "", False
            self.refresh_choice()

    def set_kind(self, kind: str, on: bool) -> None:
        kinds = [k for k in ("E2", "M1") if (k == kind and on) or (k != kind and k in self.st["kinds"])]
        self.st["kinds"] = kinds
        self.refresh_choice()

    def set_eff(self, i: int, value) -> None:
        if value is None:
            return
        self.st["eff"][i] = float(value)
        self.refresh_choice()

    def make_request(self) -> sm.ShellModelRequest:
        st = self.st
        spec = self.by_name[st["interaction"]].spec
        eff = tuple(st["eff"]) if tuple(st["eff"]) != tuple(spec.eff_charge) else None
        return sm.ShellModelRequest(self.nuclide, st["interaction"], parse_states(st["states"]),
                                    tuple(st["kinds"]), eff_charge=eff, allow_over_limit=bool(st["override"]))

    # -- 2. estimate
    def refresh_choice(self) -> None:
        """The estimate for the choices as they stand, then the Run part (and a cached result at once)."""
        ui, st = self.ui, self.st
        self.estimate_box.clear()
        self.request = self.estimate = None
        if not self.options:
            self.run_part()
            return
        with self.estimate_box:
            try:
                self.request = self.make_request()
                self.estimate = spaces.estimate(self.request)
            except ValueError as err:
                ui.label(f"Not a calculation yet: {err}").classes("text-sm ps-warn")
                self.request = self.estimate = None
            if self.estimate is not None:
                self.show_estimate(self.estimate)
        self.run_part()
        if self.request is not None and not self.job_active():
            key = cached_job(self.request, self.install)  # the same choices ran before: their result at once
            if key is not None and (st["result"] is None or st["job"] != key):
                self.take(key, origin="cache")

    def show_estimate(self, e: spaces.Estimate) -> None:
        ui, st = self.ui, self.st
        ui.label("Before running").classes("ps-section mt-1")
        dims = ", ".join(f"{sm.label(j2, p)} {_int(d)}" for (j2, p), d in e.dimension_by_state.items())
        with ui.grid(columns="auto 1fr").classes("gap-x-4 gap-y-0 text-sm"):
            for name, value in (("M-scheme dimension", dims),
                                ("Memory", f"{_bytes(e.memory_bytes)}" + (" (Lanczos vectors on disk)"
                                                                          if e.mode_lv_hdd else "")),
                                ("Disk", _bytes(e.disk_bytes)),
                                ("Limit", f"{_bytes(e.limit_bytes)}: half of the {_bytes(e.free_bytes)} of RAM free "
                                          "now"),
                                ("Rough time", f"{_seconds(e.seconds_guess)} on {e.threads} threads (could be half "
                                               "or twice that)")):
                ui.label(name).classes("ps-muted")
                ui.label(value).classes("ps-num")
        ui.label(e.reason).classes("text-sm " + ("ps-ok" if e.fits else "ps-warn"))
        if not e.fits:
            ui.checkbox("Run anyway, over the limit", value=st["override"],
                        on_change=lambda ev: self.set_override(ev.value)).props("dense").classes("ps-warn")
            ui.label("Over the limit the computer may run out of memory and slow to a crawl, or the run may fail "
                     "part way. Close other programs first, or run it on a cluster and import the results.").classes(
                "text-xs ps-warn")

    def set_override(self, on: bool) -> None:
        self.st["override"] = bool(on)
        try:
            self.request = self.make_request()
        except ValueError:
            pass
        self.run_part()

    # -- 3. run
    def job_active(self) -> bool:
        job = self.st["job"]
        if not job:
            return False
        try:
            return jobs.status(job).get("state") in jobs.ACTIVE
        except (KeyError, ValueError, OSError):
            return False

    def can_run(self) -> tuple:
        """(whether Run is allowed, why not)."""
        if self.install is None:
            return False, ("KSHELL was not found: set PHYSIM_KSHELL to the folder with kshell.exe and transit.exe, "
                           "or install physim with its bundled copy.")
        if self.request is None or self.estimate is None:
            return False, "Choose the interaction and the states first."
        if not self.estimate.fits and not self.st["override"]:
            return False, "Over the memory limit: tick \"Run anyway\" to run it all the same."
        if self.job_active():
            return False, "A calculation is running."
        return True, ""

    def run_part(self) -> None:
        ui = self.ui
        self.run_box.clear()
        with self.run_box:
            ok, why = self.can_run()
            with ui.row().classes("items-center gap-2"):
                button = ui.button("Run KSHELL", icon="play_arrow", on_click=self.start).props("dense no-caps")
                if not ok:
                    button.disable()
                if self.job_active():
                    ui.button("Cancel", icon="stop", on_click=self.cancel).props("dense flat no-caps color=negative")
            if why and not self.job_active():
                ui.label(why).classes("text-xs ps-muted")
            if self.install is not None:
                version = (self.install.version or "unknown").splitlines()[0]
                ui.label(f"KSHELL: {version}").classes("text-xs ps-muted").tooltip(str(self.install.folder))
            self.progress_box = ui.column().classes("w-full gap-0")
        self.show_progress()

    async def start(self) -> None:
        from nicegui import run

        ok, why = self.can_run()
        if not ok:
            self.ui.notify(why, type="warning")
            return
        request = self.request
        try:
            job = await run.io_bound(jobs.submit, request)
        except (FileNotFoundError, ValueError, OSError) as err:
            self.ui.notify(f"Could not start the calculation: {err}", type="negative", multi_line=True)
            return
        self.st.update(job=job, result=None, origin="run", applied=None)
        self.show_results()
        self.poll()

    def cancel(self) -> None:
        job = self.st["job"]
        if not job:
            return
        try:
            jobs.cancel(job)
        except (KeyError, OSError) as err:
            self.ui.notify(f"Could not cancel: {err}", type="negative")
        self.poll()

    def status(self) -> Optional[dict]:
        job = self.st["job"]
        if not job:
            return None
        try:
            return jobs.status(job)
        except (KeyError, ValueError, OSError):
            return None

    def poll(self) -> None:
        """Follow the job: progress while it runs, the result once done, the message if it failed."""
        s = self.status()
        if s is None:
            return
        if s.get("state") in jobs.ACTIVE:
            if self.timer is not None and not self.timer.active:
                self.timer.activate()
                self.run_part()
            else:
                self.show_progress()
            return
        if self.timer is not None:
            self.timer.deactivate()
        if s.get("state") == jobs.DONE and self.st["result"] is None:
            self.take(s.get("job") or self.st["job"], origin=self.st.get("origin") or "run")
            self.ui.notify(f"Shell-model calculation for {self.nuclide} done.", type="positive")
        self.run_part()

    def show_progress(self) -> None:
        ui = self.ui
        box = getattr(self, "progress_box", None)
        if box is None or box.is_deleted:
            return
        box.clear()
        s = self.status()
        if s is None or (s.get("state") == jobs.DONE and self.st["origin"] == "cache"):
            return
        state = s.get("state")
        with box:
            if state in jobs.ACTIVE:
                steps, at = s.get("steps") or 0, s.get("step_index") or 0
                with ui.row().classes("items-center gap-2"):
                    ui.spinner(size="sm")
                    where = f"Step {at} of {steps}: " if steps else ""
                    ui.label(f"{where}{s.get('step', state)}; {_seconds(elapsed(s))} so far").classes("text-sm")
                if steps:
                    ui.linear_progress(value=max(0.0, (at - 1) / steps), show_value=False).props("rounded")
            elif state == jobs.DONE:
                ui.label(f"Done in {_seconds(elapsed(s))}.").classes("text-sm ps-ok")
            elif state == jobs.CANCELLED:
                ui.label("Cancelled.").classes("text-sm ps-muted")
            elif state == jobs.FAILED:
                text = s.get("message") or "The calculation failed."
                if s.get("refused") == "memory":
                    text += " Tick \"Run anyway\" above to run it over the limit."
                ui.label(text).classes("text-sm ps-bad")
                if s.get("log_tail"):
                    with ui.expansion("The end of the log").classes("w-full").props("dense"):
                        ui.code(s["log_tail"], language="text").classes("w-full text-xs")

    # -- 6. import
    def import_form(self) -> None:
        ui, st = self.ui, self.st
        ui.label("A KSHELL output folder (its log files), one log, KSHELL's summary file, or a physim shell-model "
                 "table (CSV) with the numbers of a paper.").classes("text-xs ps-muted")
        path = ui.input("File or folder on this computer", value=st["import_path"]).props(
            "dense outlined clearable").classes("w-full")
        with ui.dialog() as self.upload_dialog, ui.card():
            ui.label(f"Upload shell-model results for {self.nuclide}: the log files of a KSHELL run, its summary, "
                     "or a physim shell-model table. They are kept in the experiment's folder.").classes("text-sm")
            ui.upload(multiple=True, auto_upload=True, on_multi_upload=self.import_uploads).props(
                "accept=.txt,.log,.csv,.dat").classes("w-full")
        with ui.row().classes("items-center gap-2"):
            ui.button("Import", icon="file_open", on_click=lambda: self.import_from(path.value)).props(
                "dense flat no-caps")
            ui.button("Upload files…", icon="upload_file", on_click=self.upload_dialog.open).props(
                "dense flat no-caps")

    def import_from(self, path: str) -> None:
        path = (path or "").strip().strip('"')
        self.st["import_path"] = path
        if not path:
            self.ui.notify("Type the path of a file or a folder first.", type="warning")
            return
        self.read_import(Path(path).expanduser())

    async def import_uploads(self, e) -> None:
        self.upload_dialog.close()
        folder = self.upload_folder()
        for f in e.files:
            await f.save(folder / Path(f.name).name)
        self.read_import(folder)

    def upload_folder(self) -> Path:
        """Where uploaded files are kept: the experiment's folder (so the provenance stays true), else a temporary
        folder."""
        stamp = time.strftime("%Y%m%d-%H%M%S")
        base = self.ctx.P().folder
        if base is not None:
            folder = Path(base) / "shell-model" / f"{self.nuclide}-{stamp}"
        else:
            import tempfile

            folder = Path(tempfile.mkdtemp(prefix="physim-shell-model-")) / f"{self.nuclide}-{stamp}"
        folder.mkdir(parents=True, exist_ok=True)
        return folder

    def read_import(self, path: Path) -> None:
        try:
            result = smi.import_path(path, nuclide=self.nuclide)
        except (FileNotFoundError, ValueError, NotImplementedError, OSError) as err:
            self.ui.notify(f"Could not import: {err}", type="negative", multi_line=True)
            return
        if result.nuclide and parse_nuclide(result.nuclide) != parse_nuclide(self.nuclide):
            self.ui.notify(f"These results are for {result.nuclide}, not {self.nuclide}.", type="negative")
            return
        self.st.update(result=result, origin="import", job=None, applied=None)
        self.ui.notify(f"Imported {len(result.states)} states and {len(result.transitions)} transitions.",
                       type="positive")
        self.show_results()
        self.run_part()

    # -- 4. results
    def take(self, job: str, origin: str) -> None:
        try:
            result = jobs.result(job)
        except (jobs.JobFailed, KeyError, ValueError, OSError):
            return
        if result is None:
            return
        self.st.update(job=job, result=result, origin=origin)
        self.show_results()

    def scheme(self):
        p = self.ctx.P()
        return None if p.problems else p.experiment.levels.get(self.role)

    def show_results(self) -> None:
        ui, st = self.ui, self.st
        self.results_box.clear()
        result: Optional[sm.ShellModelResult] = st["result"]
        if result is None:
            return
        cols = self.ctx.columns
        fmt = self.ctx.fmt
        scheme = self.scheme()
        with self.results_box:
            ui.separator()
            charges = (f"e_p = {result.eff_charge[0]:g}, e_n = {result.eff_charge[1]:g}"
                       if result.eff_charge else "effective charges not given")
            where = {"cache": "from the cache (run before with the same choices)", "run": "run here",
                     "import": "imported"}.get(st["origin"], "")
            ui.label(f"Shell model for {result.nuclide or self.nuclide}: {result.interaction or 'interaction not named'}"
                     f", {charges}, {where}").classes("text-sm font-semibold")
            prov = smi.provenance(result)
            if prov:
                ui.label(f"From {prov['path']} ({prov['format']})" + (f"; {prov['reference']}" if prov.get(
                    "reference") else "")).classes("text-xs ps-muted")
            for note in result.notes[1 if prov else 0:][:4]:
                ui.label(note).classes("text-xs ps-muted")
            pairs = sm.match(result, scheme) if scheme is not None and scheme.levels else [
                sm.Match(s, None) for s in sorted(result.states, key=lambda s: (s.ex_mev, s.j2))]
            rows = []
            for m in pairs:
                s = m.state
                lv = scheme.levels[m.level] if m.level is not None else None
                rows.append({"state": s.label, "known": (f"{lv.energy.value:.1f} ({lv.energy.source})"
                                                         if lv is not None else "—"),
                             "sm": f"{m.shell_kev:.1f}",
                             "delta": f"{int(round(m.delta_kev)):+d}" if m.delta_kev is not None else "",
                             "q": fmt(s.q_s, 3) if s.q_s is not None else "", "mu": fmt(s.mu, 3) if s.mu is not None
                             else ""})
            ui.label("States, beside the known levels (paired by J± and order)").classes("ps-section mt-1")
            with ui.row().classes("w-full items-start gap-4").style("flex-wrap: wrap"):
                ui.table(columns=cols((("state", "State"), ("known", "Known level (keV)"), ("sm", "Shell model (keV)"),
                                       ("delta", "SM − known"), ("q", "Q_s (e fm²)"), ("mu", "μ (μ_N)"))),
                         rows=rows).props("dense flat").classes("max-h-80")
                self.transitions_table(result, scheme, pairs)
            self.apply_part(result, scheme)

    def transitions_table(self, result, scheme, pairs) -> None:
        ui, fmt = self.ui, self.ctx.fmt
        level_of = {m.state.key: m.level for m in pairs if m.level is not None}
        unsigned = smi.unsigned(result)
        rows = []
        for t in sorted((t for t in result.transitions if not t.diagonal),
                        key=lambda t: (t.kind, t.e_initial_mev or 0.0, t.initial, t.final)):
            known = ""
            if scheme is not None and t.kind == "E2" and t.initial in level_of and t.final in level_of:
                a, b = sorted((level_of[t.initial], level_of[t.final]))
                me = scheme.matrix_element(a, b, "E2")
                if me is not None:
                    known = f"{me.value:+.4g} ({me.source})"
            unit = "e²fm⁴" if t.kind == "E2" else ("μ_N²" if t.kind == "M1" else "e²fm²")
            sign = " (sign not known)" if (t.kind, t.initial, t.final) in unsigned else ""
            rows.append({"kind": t.kind, "between": f"{sm.label(*t.initial)} → {sm.label(*t.final)}",
                         "b": f"{fmt(t.b_down, 4)} {unit}" if t.b_down is not None else "",
                         "m": f"{t.mred:+.4g}{sign}", "known": known})
        if rows:
            ui.table(columns=self.ctx.columns((("kind", ""), ("between", "Transition"), ("b", "B↓"),
                                               ("m", "⟨f‖M‖i⟩ (e fm^λ or μ_N)"), ("known", "Known ⟨‖E2‖⟩ (source)"))),
                     rows=rows).props("dense flat").classes("max-h-80")

    # -- 5. apply
    def apply_part(self, result, scheme) -> None:
        ui, st = self.ui, self.st
        choice = st["apply"]
        ui.label(f"Apply to the {self.role}'s level scheme").classes("ps-section mt-2")
        with ui.row().classes("items-center gap-3").style("flex-wrap: wrap"):
            for q, label in QUANTITIES:
                ui.checkbox(label, value=choice.get(q, True),
                            on_change=lambda e, q=q: self.set_apply(q, e.value)).props("dense")
            ui.checkbox("Add the states the scheme lacks as new levels", value=choice.get("new", False),
                        on_change=lambda e: self.set_apply("new", e.value)).props("dense")
        labels = [s.label for s in sorted(result.states, key=lambda s: (s.ex_mev, s.j2))]
        ui.select(labels, value=[x for x in choice.get("states", []) if x in labels], multiple=True,
                  label="Only these states (none chosen: all)",
                  on_change=lambda e: self.set_apply("states", list(e.value or []))).props(
            "dense outlined use-chips").classes("w-full").style("max-width: 640px")
        ui.checkbox("Replace existing values", value=st["replace"],
                    on_change=lambda e: self.set_replace(e.value)).props("dense").classes("ps-warn")
        ui.label("Ticked, the shell model's values overwrite ENSDF, derived and your own values. Not ticked (the "
                 "default), those are kept and only gaps and physim's assumed values are filled.").classes(
            "text-xs " + ("ps-warn" if st["replace"] else "ps-muted"))
        self.preview_box = ui.column().classes("w-full gap-0")
        self.show_preview(result, scheme)
        ui.button("Apply to the level scheme", icon="playlist_add_check", on_click=self.apply).props(
            "dense no-caps")
        if st["applied"]:
            ui.label("Last applied:").classes("text-xs ps-muted mt-1")
            for line in st["applied"][:10]:
                ui.label(line).classes("text-xs")

    def set_apply(self, key: str, value) -> None:
        self.st["apply"][key] = value
        self.show_preview(self.st["result"], self.scheme())

    def set_replace(self, on: bool) -> None:
        self.st["replace"] = bool(on)
        self.show_results()

    def show_preview(self, result, scheme) -> None:
        ui = self.ui
        box = self.preview_box
        if box.is_deleted:
            return
        box.clear()
        with box:
            which = which_of(self.st["apply"])
            if not which:
                ui.label("Tick what to apply.").classes("text-xs ps-muted")
                return
            try:
                lines = preview(result, scheme, which, self.st["replace"])
            except ValueError as err:
                ui.label(str(err)).classes("text-sm ps-warn")
                return
            changes = [x for x in lines if not kept(x)]
            ui.label("Applying would:").classes("text-xs ps-muted")
            for line in changes[:12]:
                ui.label("• " + line).classes("text-xs")
            if len(changes) > 12:
                ui.label(f"… and {len(changes) - 12} more").classes("text-xs ps-muted")
            for line in lines:
                if kept(line):
                    ui.label("• " + line).classes("text-xs ps-warn")

    def apply(self) -> None:
        st, ui = self.st, self.ui
        result = st["result"]
        if result is None:
            return
        which, replace = which_of(st["apply"]), bool(st["replace"])
        if not which:
            ui.notify("Tick what to apply first.", type="warning")
            return
        out: list = []

        def action(planner) -> bool:
            changes = planner.apply_shell_model(self.role, result, which, replace=replace)
            if not out:
                out.append(changes)
            return not planner.problems

        try:
            self.ctx.direct(action)
        except (ValueError, KeyError) as err:
            ui.notify(f"Not applied: {err}", type="negative", multi_line=True)
            return
        if not out:  # the setup is read-only while a run is taken
            ui.notify("A run is being taken: the setup is read-only until it ends.", type="warning")
            return
        lines = [plain(x) for x in out[0]]
        st["applied"] = lines
        if not self.results_box.is_deleted:  # the page draws the level scheme again; this panel follows it
            self.show_results()
        kept_lines = [x for x in lines if kept(x)]
        changed = [x for x in lines if not kept(x)]
        head = (f"{len(changed)} change{'s' if len(changed) != 1 else ''} to the {self.role}'s level scheme, "
                "listed in the Shell model panel." if changed else f"Nothing changed in the {self.role}'s level scheme.")
        ui.notify(head + (" " + kept_lines[0][0].upper() + kept_lines[0][1:] + "." if kept_lines else ""),
                  type="positive" if changed else "info", multi_line=True, close_button=True)
