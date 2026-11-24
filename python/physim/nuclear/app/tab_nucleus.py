"""One nucleus: everything known and projected about it, by method and across its experiments.

The page of one nucleus (``/?nucleus=62Ni``, or a card of the nuclei screen, :mod:`.start`): its observables as
rows; for the chosen row every value side by side in the two groups that are never mixed (independent evidence
against the reference; simulated runs apart, each against the value its simulation was given); the plot of the
values; the routes to the observable with "Plan it" where physim has the kind of experiment; what the value gives;
the form that takes a data point; and the history of the rows.

Nothing is computed here: every number comes from :mod:`physim.nuclear.nucleus` (:func:`~physim.nuclear.nucleus.view`,
:func:`~physim.nuclear.nucleus.observable_view`), read once on a worker thread when the page opens (:func:`load`).
The unit switch of a B(E2) only changes how the canonical value is written, through the registry's units
(:func:`display_factor`).
"""

from __future__ import annotations

import math
import re
from typing import Optional
from urllib.parse import quote

from .. import ensdf, nucleus, observables, result_rows, results, runs
from ..compare import VERDICTS
from ..results import EVIDENCE, METHODS
from .figures import COLORS, THEMES, themed

#: The tables, tags and cards of the nuclei screen and the nucleus page, in the page's colours (both themes).
NUCLEUS_STYLE = """
.ps-nt { width: 100%; border-collapse: collapse; font-size: 13px; }
.ps-nt th { text-align: left; font-weight: 500; color: var(--ps-muted); padding: 6px 10px; white-space: nowrap;
  border-bottom: 1px solid var(--ps-line); }
.ps-nt td { padding: 6px 10px; border-bottom: 1px solid var(--ps-hair); vertical-align: top; }
.ps-nt tr:last-child td { border-bottom: 0; }
.ps-nt .ps-nt-num { font-family: var(--ps-mono); font-variant-numeric: tabular-nums; text-align: right;
  white-space: nowrap; }
.ps-nt tr.ps-nt-group td { background: var(--ps-sunk); font-weight: 600; font-size: 12px; color: var(--ps-muted);
  letter-spacing: .02em; }
.ps-nt tr.ps-nt-click { cursor: pointer; }
.ps-nt tr.ps-nt-click:hover td { background: var(--ps-sunk); }
.ps-nt tr.ps-nt-on td { background: var(--ps-raise); }
.ps-nt-sub { font-size: 11.5px; color: var(--ps-muted); font-family: inherit; white-space: normal; }
.ps-tag { display: inline-block; font-size: 11px; border: 1px solid var(--ps-line); border-radius: 4px;
  padding: 0 6px; color: var(--ps-muted); font-family: var(--ps-mono); white-space: nowrap; }
.ps-tag-done { color: var(--ps-ok); border-color: var(--ps-ok); }
.ps-tag-can { color: var(--ps-accent); border-color: var(--ps-accent); }
.ps-tag-off { color: var(--ps-warn); border-color: var(--ps-warn); }
.ps-open { color: var(--ps-accent); }
.ps-two { display: grid; grid-template-columns: minmax(0, 1.35fr) minmax(0, 1fr); gap: 12px; align-items: start;
  width: 100%; }
.ps-routes { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 10px; width: 100%; }
.ps-route { background: var(--ps-surface); border: 1px solid var(--ps-line); border-radius: 6px;
  padding: 10px 12px; display: flex; flex-direction: column; gap: 5px; font-size: 13px; line-height: 1.45; }
.ps-route-later { border-style: dashed; opacity: .7; }
.ps-nuc-head { display: flex; flex-wrap: wrap; gap: 6px 18px; align-items: baseline; padding: 10px 16px;
  border-bottom: 1px solid var(--ps-line); }
.ps-nuc-body { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1.1fr); }
.ps-nuc-body > div { padding: 10px 16px 12px; min-width: 0; }
.ps-nuc-body > div + div { border-left: 1px solid var(--ps-hair); }
.ps-nuc-h { font-size: 11px; color: var(--ps-muted); letter-spacing: .04em; margin-bottom: 4px; }
.ps-kv { display: grid; grid-template-columns: auto 1fr; gap: 3px 14px; font-size: 13px; }
.ps-kv > :nth-child(odd) { color: var(--ps-muted); }
@media (max-width: 1000px) { .ps-two, .ps-nuc-body { grid-template-columns: 1fr; }
  .ps-nuc-body > div + div { border-left: 0; border-top: 1px solid var(--ps-hair); } }
"""

_SUPERSCRIPT = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
#: The units as the page writes them (the ledger's own spelling otherwise).
UNIT_TEXT = {"e2fm4": "e²fm⁴", "e2b2": "e²b²", "efm2": "e fm²", "eb": "e b", "e2fm2": "e²fm²"}
#: The units the switch of the nucleus page offers for a B(E2), of the registry's: display only.
BE2_UNITS = ("e2fm4", "e2b2")
#: A verdict of the comparison as the page says it, and its colour.
VERDICT_TEXT = {"agree": ("✓ " + VERDICTS["agree"].lower(), "ps-ok"),
                "tension": ("△ " + VERDICTS["tension"].lower(), "ps-warn"),
                "differ": ("✗ " + VERDICTS["differ"].lower(), "ps-bad")}
#: What a simulated run's closure is called: it is no test of the physics, only of the analysis on its own input.
CLOSURE_TEXT = {"agree": ("✓ closes", "ps-ok"), "tension": ("△ closes poorly", "ps-warn"),
                "differ": ("✗ does not close", "ps-bad")}
#: A route that is no experiment of the planner, in words.
MEANS_TEXT = {"calculation": "a calculation: run it in the Shell model panel of an experiment",
              "evaluation": "the evaluation: taken from ENSDF",
              "literature": "known from the literature: its value comes in as a data point"}
#: The order the independent evidence is listed in.
INDEPENDENT_ORDER = ("evaluated", "measured", "outside program", "theory")
#: How many headline rows a card or the readouts show, and how many rows of history the page lists at once.
HEADLINES, HISTORY_ROWS = 4, 25


# ---------------------------------------------------------------------------------------------------------------
# Words: a nuclide, a label, a value in the unit the page shows


def nuclide_text(nuclide: str) -> str:
    """ "62Ni" as it is written: "⁶²Ni"."""
    m = re.fullmatch(r"(\d+)([A-Za-z]+)", str(nuclide))
    return str(nuclide) if m is None else m.group(1).translate(_SUPERSCRIPT) + m.group(2)


def label_text(label: str) -> str:
    """A heading of :func:`physim.nuclear.nucleus.label` with the parities raised: "B(E2; 0⁺₁ → 2⁺₁)↑"."""
    return re.sub(r"(?<=[\d)])([+-])(?=[₀-₉ )]|$)", lambda m: "⁺" if m.group(1) == "+" else "⁻", str(label))


def unit_text(unit: str) -> str:
    return UNIT_TEXT.get(unit, unit or "")


def words(text: str) -> str:
    """A line of the ledger or the registry with its units written as the page writes them."""
    text = str(text or "")
    for unit in ("e2fm4", "e2b2", "e2fm2"):
        text = text.replace(unit, UNIT_TEXT[unit])
    return text


def display_factor(observable: str, display: Optional[str] = None) -> tuple:
    """(factor, unit): what turns a value of ``observable`` in its canonical unit into the unit the page shows it
    in. ``display`` is one of the registry's units of the observable (e²b² for a B(E2)); anything else, and a unit
    that is a convention of its own (Weisskopf units), leaves the canonical unit. Display only: no row changes."""
    o = observables.get(observable)
    size = o.units.get(display) if display else None
    return (1.0 / size, display) if size and display != o.unit else (1.0, o.unit)


def shown(value: float, plus: Optional[float] = None, minus: Optional[float] = None, unit: str = "",
          direction: Optional[str] = None, factor: float = 1.0) -> str:
    """A value with its uncertainty, unit and direction as the page writes it: "884 ± 13 e²fm⁴↑". A value in
    Weisskopf units is always marked ↓."""
    def scaled(x):
        return None if x is None else x * factor

    arrow = "↓" if unit == observables.WEISSKOPF else observables.DIRECTIONS.get(direction, "")
    return f"{nucleus.value_text(value * factor, scaled(plus), scaled(minus))} {unit_text(unit)}".strip() + arrow


def cell_text(cell: dict, display: Optional[str] = None) -> str:
    """A :func:`physim.nuclear.nucleus.cell` as the page writes it, in the unit the page shows."""
    factor, unit = display_factor(cell["observable"], display)
    if cell.get("problem"):  # not in the canonical convention: as given
        factor, unit = 1.0, cell["unit"]
    return shown(cell["value"], cell["plus"], cell["minus"], unit, cell["direction"], factor)


def as_given(cell: dict) -> str:
    """The row's own value, unit and direction, when that is not the convention the page shows: "12.13 W.u.↓"."""
    given = cell.get("as_given") or {}
    if not cell.get("converted") or given.get("value") is None:
        return ""
    return shown(given["value"], unit=given.get("unit") or "", direction=given.get("direction"))


def percent(x: Optional[float]) -> str:
    if x is None or not math.isfinite(x):
        return "—"
    return f"±{x:.1f} %" if x >= 0.095 else f"±{x:.2g} %"


def signed(x: Optional[float], digits: int = 2, unit: str = "") -> str:
    return "—" if x is None or not math.isfinite(x) else f"{x:+.{digits}f}{unit}".replace("-", "−")


def method_words(cell: dict) -> str:
    """The method of a value in words; for a value made from another row, the method of that row and what it was
    made from ("DSAM (line shape), from τ = 2.01 ± 0.12 ps")."""
    via = cell.get("via")
    if via:
        return f"{via['method_label']}, from {label_text(via['text'])}"
    return cell["method_label"]


def stat_sys(value: float, stat: Optional[float], sys: Optional[float], total: Optional[float],
             factor: float = 1.0) -> str:
    """A simulated run's value with its statistical and systematic uncertainty: "919 ± 12 (stat) ± 64 (sys)"."""
    parts = [(x * factor, name) for x, name in ((stat, "stat"), (sys, "sys")) if x is not None]
    if not parts:
        return nucleus.value_text(value * factor, *([None if total is None else total * factor] * 2))
    widths = [x for x, _ in parts if x > 0]
    digits = max(0, 1 - math.floor(math.log10(min(widths)))) if widths else 0
    return f"{value * factor:.{digits}f}" + "".join(f" ± {x:.{digits}f} ({name})" for x, name in parts)


def said(observable: str, key: str) -> str:
    """A key of the ledger as one says it in the form of a data point: "2+1 > 0+1", "2+1"."""
    try:
        if observables.get(observable).belongs_to == "transition":
            k = results.TransitionKey.parse(key)
            return f"{k.initial.level} > {k.final.level}"
        return results.LevelKey.parse(key).level
    except ValueError:
        return ""


def uses(observable: str) -> list:
    """What the registry says a value of the observable is used for, as [{"text", "reference"}]: the registry
    holds one string, or a list of entries each with a text and a reference."""
    held = getattr(observables.get(observable), "used_for", "") or ()
    out = []
    for item in ([held] if isinstance(held, (str, dict)) else list(held)):
        if isinstance(item, dict):
            text, reference = item.get("text"), item.get("reference")
        else:
            text, reference = getattr(item, "text", item), getattr(item, "reference", None)
        if str(text or "").strip():
            out.append({"text": str(text).strip(), "reference": str(reference or "").strip()})
    return out


def url(nuclide: str, observable: Optional[str] = None, key: Optional[str] = None) -> str:
    """The address of a nucleus's page in the app, for a link from an experiment: "/?nucleus=62Ni", with the
    observable and the level or transition to open on. The key is percent-encoded: a bare "+" in an address reads
    as a space."""
    out = "/?nucleus=" + quote(results.nuclide_name(nuclide), safe="")
    for name, value in (("observable", observable), ("key", key)):
        if value:
            out += f"&{name}=" + quote(str(value), safe="")
    return out


# ---------------------------------------------------------------------------------------------------------------
# Reading (on a worker thread) and writing


def backfill() -> list:
    """Rows for the experiments from before the ledger (:func:`physim.nuclear.result_rows.backfill_if_new`), so
    that their kept analyses show. Returns the experiments as :func:`physim.nuclear.runs.list_experiments` lists
    them. A folder that cannot be read is left as it is."""
    found = runs.list_experiments()
    for e in found:
        try:
            result_rows.backfill_if_new(e["path"])
        except Exception:  # noqa: BLE001 -- an experiment that cannot be read shows without rows
            pass
    return found


def load_nuclei() -> dict:
    """What the nuclei screen shows, read once: :func:`physim.nuclear.nucleus.nuclei` after the backfill, and when
    each experiment last changed ({"nuclei", "unassigned", "modified": {path: time}})."""
    listed = backfill()
    return {**nucleus.nuclei(), "modified": {e["path"]: e["modified"] for e in listed}}


def _as_typed(key: Optional[str]) -> Optional[str]:
    """A key from an address: a "+" that was not percent-encoded arrives as a space ("62Ni:2 1>0 1:E2")."""
    return None if not key else re.sub(r"(?<=\d) (?=[\d(@])", "+", str(key))


def load(nuclide: str, observable: Optional[str] = None, key: Optional[str] = None,
         reference: Optional[str] = None, first: bool = True) -> dict:
    """Everything the page of one nucleus shows, read once: {"nuclide", "view", "observable", "notes"}.

    "view" is :func:`physim.nuclear.nucleus.view` and "observable" the
    :func:`~physim.nuclear.nucleus.observable_view` of the row asked for (the first row of the table when none is,
    None for a nucleus without rows). If the evaluation was never taken for the nucleus it is taken first
    (:func:`~physim.nuclear.nucleus.snapshot_evaluation`). "notes" says what could not be done: no local ENSDF
    copy, a row or a reference that is not there. ``first`` is False when the page reads again after a change of
    its own (the experiments from before the ledger have their rows by then). Raises ``ValueError`` for a name that
    is no nuclide."""
    name = results.nuclide_name(nuclide)
    if first:
        backfill()
    notes = []
    view = nucleus.view(name)
    if view["evaluation"]["taken"] is None:
        if not ensdf.available():
            notes.append("There is no local copy of ENSDF, so the adopted values are not here to compare with "
                         "(scripts/fetch_ensdf.py fetches it).")
        else:
            try:
                taken = nucleus.snapshot_evaluation(name)
                if taken["added"] or taken["superseded"]:
                    view = nucleus.view(name)
            except (ensdf.EnsdfMissing, ValueError, OSError) as e:
                notes.append(f"The adopted values of ENSDF could not be taken: {e}")
    return {"nuclide": name, "view": view, "notes": notes,
            "observable": load_observable(name, view, observable, key, reference, notes)}


def load_observable(name: str, view: dict, observable: Optional[str] = None, key: Optional[str] = None,
                    reference: Optional[str] = None, notes: Optional[list] = None) -> Optional[dict]:
    """The :func:`physim.nuclear.nucleus.observable_view` of the row asked for, or of the first row of the
    nucleus's table. What is asked for and is not there is said in ``notes`` and the nearest thing is shown."""
    notes = [] if notes is None else notes
    table = view["table"]
    chosen = None
    if observable and key:
        try:
            chosen = (observable, nucleus.key_of(name, observable, _as_typed(key), check=False).text)
        except ValueError as e:
            notes.append(f"'{key}' of {observable} cannot be shown: {e}")
    if chosen is None:
        row = next((t for t in table if not observable or t["observable"] == observable), table[0] if table else None)
        chosen = None if row is None else (row["observable"], row["key"])
    if chosen is None:
        return None
    try:
        return nucleus.observable_view(name, *chosen, reference=reference)
    except ValueError as e:
        if reference is None:
            notes.append(str(e))
            return None
        notes.append(f"The reference asked for is not there ({e}); the evaluated value is the reference again.")
        return nucleus.observable_view(name, *chosen)


def blank_point(observable: str = "B(E2)", key: str = "") -> dict:
    """The form of a data point, empty, for an observable and a level or transition."""
    o = observables.get(observable)
    return {"observable": observable, "key": said(observable, key), "value": "", "unit": o.unit,
            "direction": o.direction, "form": "total", "uncertainty": "", "stat": "", "sys": "", "plus": "",
            "minus": "", "evidence": "measured", "method": None, "from": "reference", "authors": "", "year": "",
            "doi": "", "where": "", "date": "", "note": "", "detail": "", "corrects": None}


def point_of(row) -> dict:
    """The form filled with a data point of the ledger, to correct it."""
    def text(x) -> str:
        return "" if x is None else f"{x:.12g}"

    origin = row.origin
    form = blank_point(row.observable, row.key)
    form.update(value=text(row.value), unit=row.unit, direction=row.direction, evidence=row.evidence,
                method=row.method, detail=row.detail, corrects=row.id,
                form="asymmetric" if row.plus is not None else "stat-sys" if (
                    row.stat is not None or row.sys is not None) else "total",
                uncertainty=text(row.total), stat=text(row.stat), sys=text(row.sys), plus=text(row.plus),
                minus=text(row.minus))
    form["from"] = "reference" if origin.get("kind") == "reference" else "own"
    form.update({k: str(origin.get(k) or "") for k in ("authors", "year", "doi", "where", "date", "note")})
    return form


def _number(text):
    """A number typed into the form: None for an empty field, the text itself when it is no number (the checks
    of the data point then name it)."""
    text = str("" if text is None else text).strip().replace("−", "-")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return text


def point_arguments(form: dict) -> dict:
    """The form of a data point as the arguments of :func:`physim.nuclear.nucleus.add_data_point`: the numbers
    read, the uncertainty in the one form chosen, and a reference or an own measurement."""
    o = observables.OBSERVABLES.get(form.get("observable"))
    value = _number(form.get("value"))
    args = {"observable": form.get("observable"), "key": str(form.get("key") or "").strip(),
            "value": "" if value is None else value, "unit": form.get("unit") or None,
            "direction": form.get("direction") if o is not None and o.direction else None,
            "evidence": form.get("evidence"), "method": form.get("method"),
            "detail": str(form.get("detail") or "").strip(), "corrects": form.get("corrects") or None}
    kind = form.get("form", "total")
    names = {"total": ("uncertainty",), "stat-sys": ("stat", "sys"), "asymmetric": ("plus", "minus")}[kind]
    args.update({name: _number(form.get(name)) for name in names})
    if form.get("from") == "own":
        args["own"] = {"date": str(form.get("date") or "").strip(), "note": str(form.get("note") or "").strip()}
    else:
        year = _number(form.get("year"))
        args["reference"] = {"authors": str(form.get("authors") or "").strip(),
                             "year": int(year) if isinstance(year, float) and year.is_integer() else year,
                             "doi": str(form.get("doi") or "").strip(), "where": str(form.get("where") or "").strip()}
    return args


def add_point(nuclide: str, form: dict, home=None):
    """Add the data point the form holds to the nucleus's ledger (a correction when the form names the row it
    replaces) and return its row. Raises :class:`physim.nuclear.nucleus.PointProblem` with everything that is
    wrong: the page shows the messages beside the form."""
    args = point_arguments(form)
    return nucleus.add_data_point(nuclide, args.pop("observable"), args.pop("key"), args.pop("value"),
                                  args.pop("unit"), home=home, **args)


# ---------------------------------------------------------------------------------------------------------------
# The plot


def figure_values(view: dict, theme: str = "light", display: Optional[str] = None):
    """The values of one observable of one level or transition (an :func:`~physim.nuclear.nucleus.observable_view`)
    with their uncertainties, one line each: the independent evidence against the reference, whose uncertainty is
    the band; under a rule the simulated runs, open markers, each with the value its simulation was given."""
    import plotly.graph_objects as go

    c, cells = view["comparison"], view["cells"]
    factor, unit = display_factor(view["observable"], display)
    ink, grid = THEMES[theme]["ink"], THEMES[theme]["grid"]
    fig = go.Figure()
    names = []

    def trace(found: list, name: str, colour: str, symbol: str = "circle", truth: bool = False) -> None:
        if not found:
            return
        y = list(range(-len(names), -len(names) - len(found), -1))
        names.extend(cells[e.id]["source"] for e in found)
        widths = [(e.plus, e.minus) if not truth else (e.total, e.total) for e in found]
        fig.add_trace(go.Scatter(
            x=[e.value * factor for e in found], y=y, mode="markers", name=name,
            marker=dict(size=10, color=colour, symbol=symbol, line=dict(width=2, color=colour)),
            error_x=dict(type="data", symmetric=False, array=[(p or 0.0) * factor for p, _ in widths],
                         arrayminus=[(m or 0.0) * factor for _, m in widths], color=colour, thickness=1.6, width=4),
            customdata=[[method_words(cells[e.id]), cell_text(cells[e.id], display)] for e in found],
            hovertemplate="%{customdata[1]}<br>%{customdata[0]}<extra>" + name + "</extra>"))
        if truth:
            fig.add_trace(go.Scatter(
                x=[e.truth * factor for e in found], y=y, mode="markers", name="given to the simulation",
                marker=dict(size=14, color=colour, symbol="line-ns-open", line=dict(width=2, color=colour)),
                hovertemplate="the simulation was given %{x:.5g}<extra></extra>"))

    ref = c.reference
    if ref is not None:
        if ref.plus is not None:
            fig.add_vrect(x0=(ref.value - ref.minus) * factor, x1=(ref.value + ref.plus) * factor,
                          fillcolor=COLORS[4], opacity=0.22, line_width=0, layer="below")
        fig.add_vline(x=ref.value * factor, line=dict(color=COLORS[4], dash="dash", width=1))
    for evidence, colour in zip(INDEPENDENT_ORDER, (ink, COLORS[0], COLORS[3], COLORS[2])):
        trace([e for e in c.entries if e.evidence == evidence], EVIDENCE[evidence].lower(), colour)
    if c.projected and names:
        fig.add_hline(y=0.5 - len(names), line=dict(color=grid, dash="dot", width=1))
    trace(list(c.projected), "projected (simulated runs)", COLORS[1], "circle-open", truth=True)
    fig.update_layout(height=max(230, 120 + 26 * len(names)), margin=dict(l=10, r=14, t=10, b=10),
                      legend=dict(orientation="h", yanchor="bottom", y=1.0, x=0), hovermode="closest")
    fig.update_xaxes(title_text=f"{label_text(view['label'])} ({unit_text(unit)})")
    fig.update_yaxes(tickmode="array", tickvals=list(range(0, -len(names), -1)), ticktext=names, automargin=True,
                     range=[0.6 - len(names), 0.6], showgrid=False)
    return themed(fig, theme)


# ---------------------------------------------------------------------------------------------------------------
# Small parts of both screens


def add_style(ctx) -> None:
    """The style of the nuclei screen and the nucleus page, once per browser page."""
    if not ctx.state.get("nucleus_style"):
        ctx.ui.add_css(NUCLEUS_STYLE)
        ctx.state["nucleus_style"] = True


def later(ui, job) -> None:
    """Start a coroutine of the page on its own (as the page's ``later`` does): the element whose event asked for
    it may be rebuilt before the coroutine ends, and NiceGUI cancels what a deleted element started."""
    from nicegui import background_tasks

    client = ui.context.client

    async def run_() -> None:
        with client:
            await job()

    background_tasks.create(run_(), name="nucleus page")


def waiting(ui, text: str) -> None:
    with ui.row().classes("w-full items-center gap-3 p-4"):
        ui.spinner(size="md")
        ui.label(text).classes("text-sm ps-muted")


def open_nuclei(ctx) -> None:
    """Show the nuclei screen."""
    ctx.state["start"].update(view="nuclei")
    ctx.refresh_start()


def open_nucleus(ctx, nuclide: str, observable: Optional[str] = None, key: Optional[str] = None) -> None:
    """Show the page of one nucleus, on a row of its table if one is named."""
    try:
        name = results.nuclide_name(str(nuclide).strip())
    except (ValueError, KeyError, TypeError, IndexError):
        ctx.ui.notify(f"'{nuclide}' is not a nuclide: write it like 62Ni or Ni-62.", type="warning")
        return
    ctx.state["start"].update(view="nucleus", nucleus=name, observable=observable, key=key, reference=None,
                              point=None, point_problems=[])
    ctx.refresh_start()
    ctx.ui.run_javascript("window.scrollTo({top: 0})")


def plan(ctx, nuclide: str, kind: Optional[str] = None) -> None:
    """Today's new-experiment form for a kind of experiment (the kinds to pick from when none is named), with the
    nucleus filled in where the form has a field for it (:func:`.start.kind_form`). Nothing is computed."""
    st = ctx.state["start"]
    st.update(view="kinds", kind=kind, for_nucleus=nuclide)
    st.pop("values", None)
    ctx.refresh_start()
    ctx.ui.run_javascript("window.scrollTo({top: 0})")


def _table(ui, head: tuple, numeric: tuple = ()):
    """A table with the column heads ``head`` (those at the places ``numeric`` right-aligned, as numbers are); the
    rows are built inside it."""
    table = ui.element("table").classes("ps-nt")
    with table:
        with ui.element("tr"):
            for i, text in enumerate(head):
                with ui.element("th").classes("ps-nt-num" if i in numeric else ""):
                    ui.label(text)
    return table


def _td(ui, text: str = "", sub=(), classes: str = "", colour: str = ""):
    """A cell: its text and, under it, small lines."""
    cell = ui.element("td").classes(classes)
    with cell:
        if text:
            ui.label(text).classes(colour)
        for line in ([sub] if isinstance(sub, str) else sub):
            if line:
                ui.label(line).classes("ps-nt-sub")
    return cell


def _group(ui, text: str, columns: int, colour: str = "") -> None:
    with ui.element("tr").classes("ps-nt-group"):
        with ui.element("td").props(f"colspan={columns}"):
            ui.label(text).classes(colour)


def latest_text(x: dict, display: Optional[str] = None) -> tuple:
    """(text, small line) of an experiment's latest result: its value and, in words, what kind of evidence it is
    and by which method; "not analysed" or "no run" without one."""
    latest = x.get("latest")
    if latest is None:
        return ("no run" if not x["runs"] else "not analysed" if not x["analysed_runs"] else "no result kept"), ""
    return (cell_text(latest, display),
            f"{label_text(latest['label'])} · {EVIDENCE[latest['evidence']].lower()} · {method_words(latest)}")


def experiments_table(ctx, experiments: list, display: Optional[str] = None) -> None:
    """The experiments of a nucleus: name, kind, runs and the latest result or "not analysed". A row opens the
    experiment, as the list of experiments always did."""
    ui = ctx.ui
    with _table(ui, ("Experiment", "Kind", "Runs", "Latest result"), numeric=(2, 3)):
        for x in experiments:
            row = ui.element("tr").classes("ps-nt-click")
            row.on("click", lambda path=x["path"]: ctx.open_experiment(path))
            with row:
                _td(ui, x["name"], "" if x["role"] == "studied" else (
                    "as the partner" if x["role"] == "partner" else "results only"), colour="ps-accent")
                with ui.element("td"):
                    ui.label(x["kind_label"] or "").classes("ps-tag")
                _td(ui, f"{x['runs']} run{'s' * (x['runs'] != 1)}", classes="ps-nt-num")
                text, sub = latest_text(x, display)
                _td(ui, text, sub, classes="ps-nt-num", colour="" if x.get("latest") else "ps-muted")


def summary_text(summary: Optional[dict], display: Optional[str] = None) -> tuple:
    """(text, small line) of the rows of one kind of independent evidence in a row of the observables table: the
    one value, or how many and their range; under it the best one's source and method in words."""
    if summary is None:
        return "—", ""
    best = summary["best"]
    words = f"{best['source']} · {method_words(best)}"
    if summary["count"] == 1:
        return cell_text(best, display), words
    factor, unit = display_factor(best["observable"], display)
    return (f"{summary['count']} values, {summary['min'] * factor:.4g} to {summary['max'] * factor:.4g} "
            f"{unit_text(unit)}",
            f"best {percent(best['uncertainty_percent'])}: {words}" if best["uncertainty_percent"] is not None
            else f"newest: {words}")


def projected_text(summary: Optional[dict]) -> tuple:
    """(text, small line) of the simulated runs in a row of the observables table: the precision the best one
    reaches, how many there are and by which method. Never a value to set beside the independent ones."""
    if summary is None:
        return "—", ""
    best = summary["best"]
    count = summary["count"]
    return (f"reaches {percent(best['uncertainty_percent'])}",
            f"{count} simulated run{'s' * (count != 1)} · {method_words(best)}")


def known_lines(row: dict, display: Optional[str] = None) -> str:
    """One line on a row of the observables table for a card: what is evaluated, measured and calculated, and,
    said apart, what the simulated runs reach."""
    parts = []
    if row["evaluated"] is not None:
        parts.append(f"{cell_text(row['evaluated'], display)} evaluated")
    for name, word in (("measured", "measured"), ("outside_program", "from an outside program"), ("theory", "theory")):
        s = row[name]
        if s is not None:
            parts.append(f"{cell_text(s['best'], display)} {word} ({method_words(s['best'])})" if s["count"] == 1
                         else f"{s['count']} {word}")
    if row["projected"] is not None:
        text, sub = projected_text(row["projected"])
        parts.append(f"projected: {text} ({sub})")
    return " · ".join(parts) or "nothing yet"


# ---------------------------------------------------------------------------------------------------------------
# The page of one nucleus


def render(ctx) -> None:
    """The page of the nucleus ``ctx.state["start"]["nucleus"]``, into the current container. ``ctx`` is the start
    page's (:func:`.start.render`). What it shows is read on a worker thread (:func:`load`), with the evaluation
    taken first if it never was; a spinner stands in meanwhile."""
    ui = ctx.ui
    st = ctx.state["start"]
    name = st["nucleus"]
    add_style(ctx)
    held: dict = {}
    with ui.column().classes("w-full gap-3").style("max-width: 1240px; margin: 0 auto"):
        with ui.row().classes("w-full items-center gap-1"):
            ui.button("Nuclei", on_click=lambda: open_nuclei(ctx)).props("flat no-caps dense").classes("ps-accent")
            ui.label("▸ " + nuclide_text(name)).classes("text-sm ps-muted")
            ui.space()
            if ctx.close_start is not None:
                ui.button("Back to the experiment", icon="arrow_back", on_click=ctx.close_start).props(
                    "flat no-caps").classes("ps-muted")
        body = ui.column().classes("w-full gap-3")
    with body:
        waiting(ui, f"Reading what is known about {nuclide_text(name)}…")

    async def reload() -> None:
        """Read everything again and build the page from it, in place."""
        from nicegui import run

        try:
            data = await run.io_bound(load, name, st.get("observable"), st.get("key"), st.get("reference"),
                                      not held)
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty page
            data = {"error": f"The page of {name} could not be read: {err}"}
        if data is None or body.is_deleted:
            return
        held.clear()
        held.update(data)
        fill()

    async def choose(observable: str, key: str, reference: Optional[str] = None) -> None:
        """Another row of the table, or another reference for the row shown: only that row is read again."""
        from nicegui import run

        notes: list = []
        try:
            found = await run.io_bound(load_observable, name, held["view"], observable, key, reference, notes)
        except Exception as err:  # noqa: BLE001
            ui.notify(str(err), type="negative", multi_line=True)
            return
        if body.is_deleted:
            return
        held.update(observable=found, notes=notes)
        fill()

    def fill() -> None:
        ov = held.get("observable")
        if ov is not None:  # what the page shows now is what a reload shows again
            ref = ov["comparison"].reference if ov["comparison"] is not None else None
            st.update(observable=ov["observable"], key=ov["key"],
                      reference=ref.id if ref is not None and st.get("reference") else None)
        body.clear()
        with body:
            _page(ctx, held, lambda: later(ui, reload), lambda *a: later(ui, lambda: choose(*a)), fill)

    later(ui, reload)


def _page(ctx, held: dict, reload, choose, redraw) -> None:
    """The page from what was read. ``reload()`` reads everything again, ``choose(observable, key, reference)``
    another row or reference, ``redraw()`` builds the page again from what is held (the unit switch)."""
    ui = ctx.ui
    st = ctx.state["start"]
    if held.get("error"):
        ui.label(held["error"]).classes("text-sm ps-warnbox rounded p-3 w-full")
        return
    name, view, ov = held["nuclide"], held["view"], held["observable"]
    display = st.get("unit")
    title = nuclide_text(name)

    # -- the head: the nucleus, what it holds, and the two things one does from here
    studied = [x for x in view["experiments"] if x["role"] == "studied"]
    with ui.row().classes("w-full items-end gap-3"):
        with ui.column().classes("gap-0"):
            ui.label(title).classes("text-2xl font-semibold")
            ui.label(f"Z = {view['z']}, N = {view['n']} · everything known and projected about it, by method · "
                     f"{len(studied)} experiment{'s' * (len(studied) != 1)}, {view['analysed_runs']} analysed "
                     f"run{'s' * (view['analysed_runs'] != 1)}, {view['points']} data "
                     f"point{'s' * (view['points'] != 1)}").classes("text-sm ps-muted")
        ui.space()
        ui.button("Add a data point", icon="add", on_click=lambda: ui.run_javascript(
            "document.querySelector('.ps-point-form')?.scrollIntoView({behavior: 'smooth', block: 'start'})")).props(
            "outline no-caps")
        ui.button(f"New experiment for {title}", icon="add", on_click=lambda: plan(ctx, name)).props(
            "unelevated no-caps")
    for note in held.get("notes") or ():
        ui.label(note).classes("text-sm ps-note p-2 w-full")
    evaluation = view["evaluation"]
    if evaluation["taken"]:
        ui.label(f"The evaluation: {evaluation['reference']}, {evaluation['rows']} values, taken "
                 f"{evaluation['taken'][:10]}.").classes("text-xs ps-muted")

    # -- the readouts: the headline rows, and what is still open
    headline = [t for t in view["table"] if t["headline"]]
    with ui.element("div").classes("ps-readouts"):
        for t in headline[:HEADLINES - 1]:
            best = t["evaluated"] or next((t[k]["best"] for k in ("measured", "outside_program", "theory")
                                           if t[k] is not None), None)
            with ui.element("div").classes("ps-readout"):
                ui.label(label_text(t["label"]) + (", evaluated" if t["evaluated"] else "")).classes(
                    "text-xs ps-muted")
                ui.label(cell_text(best, display) if best else "no independent value").classes(
                    "ps-readout-value" + ("" if best else " ps-muted"))
                ui.label(("" if t["evaluated"] or not best else method_words(best) + " · ")
                         + (f"projected: {projected_text(t['projected'])[0]}" if t["projected"] else
                            "no simulated run")).classes("text-xs ps-muted")
        with ui.element("div").classes("ps-readout"):
            ui.label("Still open").classes("text-xs ps-muted")
            ui.label(str(len(view["open"]))).classes("ps-readout-value" + (" ps-warn" if view["open"] else ""))
            ui.label("routes physim can plan, runs to analyse" if view["open"] else "nothing").classes(
                "text-xs ps-muted")

    _observables(ctx, view, ov, choose, redraw)
    if ov is None:
        ui.label(f"No value of {title} is here yet. Add a data point below, or plan an experiment for it."
                 ).classes("text-sm ps-note p-3 w-full")
    else:
        with ui.element("div").classes("ps-two"):
            _comparison(ctx, held, reload, choose, redraw)
            with ui.column().classes("w-full gap-3 min-w-0"):
                _plot(ctx, ov)
                _projected_note(ctx, ov)
        _routes(ctx, name, ov)
    with ui.element("div").classes("ps-two"):
        with ui.column().classes("w-full gap-3 min-w-0"):
            if ov is not None:
                _gives(ctx, ov)
            _experiments(ctx, view)
        with ui.column().classes("w-full gap-3 min-w-0"):
            _point_form(ctx, name, ov, reload, redraw)
            if ov is not None:
                _history(ctx, ov)


def _panel(ui, title: str, sub: str = ""):
    """A bordered panel with its head line; (the head, to add to; the body, to build in)."""
    with ui.element("div").classes("ps-panel"):
        head = ui.element("div").classes("ps-panel-head")
        with head:
            ui.label(title)
            if sub:
                ui.label(sub).classes("ps-muted text-sm")
        return head, ui.element("div")


def _observables(ctx, view: dict, ov: Optional[dict], choose, redraw) -> None:
    """The observables of the nucleus, one row per level or transition: evaluated, measured, theory and what the
    simulated runs reach in columns of their own, what is still open at the end. A row opens its routes."""
    ui = ctx.ui
    st = ctx.state["start"]
    display = st.get("unit")
    on = (ov["observable"], ov["key"]) if ov is not None else None
    listed = [t for t in view["table"] if t["headline"] or (t["observable"], t["key"]) == on]
    rest = [t for t in view["table"] if not t["headline"]]
    outside = any(t["outside_program"] for t in listed)
    head, body = _panel(ui, "Observables", "one row per level or transition; click a row for its routes")
    with head:
        ui.space()
        if any(t["observable"] == "B(E2)" for t in view["table"]):
            def switch(e) -> None:
                st["unit"] = e.value
                redraw()

            ui.toggle({u: f"B(E2)↑ in {unit_text(u)}" for u in BE2_UNITS}, value=display or BE2_UNITS[0],
                      on_change=switch).props("dense no-caps unelevated").tooltip(
                "How a B(E2) is written here. The rows are not changed; a value in W.u. is always the ↓ value.")
    with body:
        if not view["table"]:
            ui.label("No row yet.").classes("text-sm ps-muted p-3")
            return
        columns = ("Observable", "Evaluated", "Measured") + (("Outside program",) if outside else ()) + (
            "Theory", "Projected (simulated runs)", "Still open")
        if listed:
            with _table(ui, columns):
                for t in listed:
                    row = ui.element("tr").classes(
                        "ps-nt-click" + (" ps-nt-on" if (t["observable"], t["key"]) == on else ""))
                    row.on("click", lambda t=t: choose(t["observable"], t["key"]))
                    with row:
                        _td(ui, label_text(t["label"]), colour="font-medium")
                        _td(ui, cell_text(t["evaluated"], display) if t["evaluated"] else "—",
                            as_given(t["evaluated"]) and "given as " + as_given(t["evaluated"])
                            if t["evaluated"] else "", classes="ps-nt-num")
                        for kind in ("measured",) + (("outside_program",) if outside else ()) + ("theory",):
                            _td(ui, *summary_text(t[kind], display), classes="ps-nt-num")
                        _td(ui, *projected_text(t["projected"]), classes="ps-nt-num")
                        _td(ui, " · ".join(METHODS[m] for m in t["open"]) or "nothing",
                            colour="ps-open" if t["open"] else "ps-muted")
        if rest:
            with ui.expansion(f"{len(rest)} more hold the evaluation only", icon="list").props("dense").classes(
                    "w-full text-sm"):
                rows = [{"id": f"{t['observable']}|{t['key']}", "label": label_text(t["label"]),
                         "value": cell_text(t["evaluated"], display) if t["evaluated"] else "—"} for t in rest]
                more = ui.table(columns=[{"name": "label", "label": "Observable", "field": "label", "align": "left"},
                                         {"name": "value", "label": "Evaluated", "field": "value", "align": "right"}],
                                rows=rows, row_key="id", pagination=12).props("dense flat").classes("w-full")
                more.add_slot("body-cell-label", '<q-td :props="props" style="cursor: pointer; font-family: inherit">'
                                                 '<span class="ps-accent">{{ props.value }}</span></q-td>')
                more.on("rowClick", lambda e: choose(*e.args[1]["id"].split("|", 1)))


def _comparison(ctx, held: dict, reload, choose, redraw) -> None:
    """Every value of the chosen row, side by side, in two tables that are never mixed: the independent evidence
    against the reference (difference, pull, verdict), and the simulated runs, each against the value its
    simulation was given."""
    ui = ctx.ui
    st = ctx.state["start"]
    name, ov = held["nuclide"], held["observable"]
    c, cells = ov["comparison"], ov["cells"]
    display = st.get("unit")
    factor, unit = display_factor(ov["observable"], display)
    heading = f"{label_text(ov['label'])}"
    head, body = _panel(ui, heading, "every value, side by side")
    if c is None:
        with body:
            ui.label("No value yet: add a data point, or analyse a run of an experiment.").classes(
                "text-sm ps-muted p-3")
        return
    ref = c.reference
    with head:
        ui.space()
        ui.label("against").classes("text-xs ps-muted")
        options = {e.id: f"{cells[e.id]['source']} · {cell_text(cells[e.id], display)}" for e in c.entries}
        ui.select(options, value=ref.id if ref is not None else None, label=None if ref is not None else "pick a row",
                  on_change=lambda e: (st.update(reference=e.value), choose(ov["observable"], ov["key"], e.value))
                  ).props("dense outlined options-dense").style("min-width: 220px; max-width: 340px").tooltip(
            "The row every other value is compared with: the evaluated one unless another is picked")
    with body:
        # -- independent: measured, evaluated, theory, outside program, each against the reference
        value_head = f"{c.label} ({unit_text(unit)})"
        with _table(ui, ("Source", "Method", value_head, "Known to", "Difference", "Pull", "Verdict", ""),
                    numeric=(2, 3, 4, 5)):
            _group(ui, "INDEPENDENT · " + ("against the reference: " + cells[ref.id]["source"] if ref is not None
                                             else "no reference yet: pick a row to compare against"), 8)
            for evidence in INDEPENDENT_ORDER:
                found = [e for e in c.entries if e.evidence == evidence]
                if not found:
                    continue
                _group(ui, EVIDENCE[evidence].upper(), 8)
                for e in found:
                    cell = cells[e.id]
                    with ui.element("tr"):
                        _td(ui, cell["source"], ["the reference" if e.is_reference else "",
                                                 "◦ not independent of the reference"
                                                 if e.independent_of_reference is False else ""])
                        given = as_given(cell)
                        _td(ui, method_words(cell), [cell["detail"],
                                                     f"given as {given}: {words(cell['converted'])}"
                                                     if given and not cell.get("via") else "", cell["problem"]])
                        _td(ui, shown(e.value, e.plus, e.minus, factor=factor), classes="ps-nt-num")
                        _td(ui, percent(e.uncertainty_percent) if e.uncertainty_percent is not None
                            else "no uncertainty", classes="ps-nt-num",
                            colour="" if e.uncertainty_percent is not None else "ps-muted")
                        if e.is_reference or ref is None:
                            for _ in range(2):
                                _td(ui, "—", classes="ps-nt-num", colour="ps-muted")
                            _td(ui, "reference" if e.is_reference else "—", colour="ps-muted")
                        else:
                            _td(ui, signed(e.difference_percent, 1, " %"), classes="ps-nt-num")
                            _td(ui, signed(e.pull) if e.pull is not None else (
                                f"× {e.ratio:.2f}" if e.ratio is not None else "—"), classes="ps-nt-num")
                            text, colour = VERDICT_TEXT.get(e.verdict, ("not judged: no uncertainty", "ps-muted"))
                            _td(ui, text, colour=colour)
                        with ui.element("td"):
                            if nucleus.is_data_point(e.row):
                                _row_menu(ctx, name, e.row, cell, reload, redraw)
        ui.label(f"Pull = the difference over the two uncertainties in quadrature (the side of each that faces the "
                 f"other value). Below {c.agree_below:g}: agree; {c.agree_below:g} to {c.differ_above:g}: tension; "
                 f"above {c.differ_above:g}: differ. A value without an uncertainty gets its ratio and no verdict. "
                 f"◦ not independent of the reference: the two share a source, so their agreement is no test."
                 ).classes("text-xs ps-muted").style("padding: 8px 12px; line-height: 1.5")
        # -- projected: simulated runs analysed back, apart, each against what its simulation was given
        with ui.element("div").style("border-top: 1px solid var(--ps-line)"):
            with _table(ui, ("Simulated run", "Method", f"{c.label} ({unit_text(unit)})", "Given to the simulation",
                             "Closure", "Reaches"), numeric=(2, 3, 5)):
                _group(ui, "PROJECTED · simulated runs analysed back. Each is closed against the value its "
                           "simulation was given, never against the reference", 6, "ps-open")
                if not c.projected:
                    with ui.element("tr"):
                        _td(ui, "No simulated run of this is analysed yet.", colour="ps-muted").props("colspan=6")
                for p in c.projected:
                    cell = cells[p.id]
                    with ui.element("tr"):
                        _td(ui, cell["source"])
                        _td(ui, method_words(cell), cell["detail"])
                        _td(ui, stat_sys(p.value, p.stat, p.sys, p.total, factor), classes="ps-nt-num")
                        _td(ui, f"{p.truth * factor:.5g}", classes="ps-nt-num")
                        text, colour = CLOSURE_TEXT.get(p.closure_verdict, ("no uncertainty", "ps-muted"))
                        _td(ui, f"{text} ({signed(p.closure)}σ)" if p.closure is not None else text, colour=colour)
                        _td(ui, percent(p.reach_percent), classes="ps-nt-num")


def _row_menu(ctx, name: str, row, cell: dict, reload, redraw) -> None:
    """Correct and retract, on a data point of the ledger. Neither changes the row: each adds one."""
    ui = ctx.ui
    st = ctx.state["start"]

    def correct() -> None:
        st.update(point=point_of(row), point_problems=[])
        redraw()
        ui.run_javascript("document.querySelector('.ps-point-form')?.scrollIntoView({behavior: 'smooth'})")

    def retract() -> None:
        with ui.dialog() as dialog, ui.card().classes("ps-card").style("min-width: min(480px, 95vw)"):
            ui.label("Retract this data point").classes("text-base font-semibold")
            ui.label(f"{cell['source']} · {method_words(cell)} · {cell_text(cell)}").classes("text-sm ps-muted")
            ui.label("Nothing is deleted: the point leaves what holds now and stays in the history, with the "
                     "reason beside it.").classes("text-xs ps-muted")
            reason = ui.input("Why it is withdrawn").props("dense outlined autofocus").classes("w-full")
            problem = ui.label("").classes("text-sm ps-bad")

            def go() -> None:
                try:
                    nucleus.retract_data_point(name, row.id, reason.value)
                except ValueError as e:
                    problem.set_text(str(e))
                    return
                dialog.close()
                if st.get("reference") == row.id:
                    st["reference"] = None
                ui.notify("The data point is withdrawn; it stays in the history.", type="positive")
                reload()

            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", color=None, on_click=dialog.close).props("flat no-caps").classes("ps-muted")
                ui.button("Retract", on_click=go).props("unelevated no-caps")
        dialog.open()

    with ui.button(icon="more_vert", color=None).props("flat round dense size=sm").classes("ps-muted").tooltip(
            "Correct or retract this data point"):
        with ui.menu():
            ui.menu_item("Correct…", on_click=correct)
            ui.menu_item("Retract…", on_click=retract)


def _plot(ctx, ov: dict) -> None:
    """The values with their uncertainties and the reference band, drawn as the app's other figures are and in the
    page's theme (drawn again when the theme changes)."""
    ui = ctx.ui
    st = ctx.state["start"]
    if ov["comparison"] is None:
        return
    _, body = _panel(ui, "The values", "the band is the reference's uncertainty")
    drawn = {"theme": ctx.state.get("theme", "light")}

    @ui.refreshable
    def figure() -> None:
        ui.plotly(themed_json(ov, drawn["theme"], st.get("unit"))).classes("w-full")

    def follow_theme() -> None:
        if ctx.state.get("theme", "light") != drawn["theme"]:
            drawn["theme"] = ctx.state.get("theme", "light")
            figure.refresh()

    with body:
        with ui.element("div").classes("w-full").style("padding: 6px 8px"):
            figure()
        ui.timer(0.5, follow_theme)


def themed_json(ov: dict, theme: str, display: Optional[str]) -> dict:
    """The figure of the values as the browser takes it."""
    return figure_values(ov, theme, display).to_plotly_json()


def _projected_note(ctx, ov: dict) -> None:
    """What the simulated runs say about the setup, beside what is already known: never a test of the value."""
    c = ov["comparison"]
    reaching = [p for p in (c.projected if c is not None else ()) if p.reach_percent is not None]
    if not reaching:
        return
    ui = ctx.ui
    best = min(reaching, key=lambda p: p.reach_percent)
    cell = ov["cells"][best.id]
    known = (f"the value is known to {percent(best.known_percent)} (the reference)" if best.known_percent is not None
             else "there is no reference with an uncertainty to set it beside")
    measured = (f", and to {percent(best.best_measured_percent)} by the best single measurement"
                if best.best_measured_percent is not None else "")
    share = ""
    if best.sys is not None and best.value:
        share = f", of which {percent(100 * best.sys / abs(best.value))[1:]} is systematic"
    with ui.column().classes("w-full ps-note p-3 gap-1"):
        ui.label("What the projected runs say").classes("text-sm font-semibold")
        ui.label(f"This setup reaches {percent(best.reach_percent)}{share}; {known}{measured}. "
                 f"({cell['source']}, {method_words(cell)}: a simulated run analysed back. It says what the setup "
                 f"would reach, not what the value is.)").classes("text-sm").style("line-height: 1.5")


def _routes(ctx, name: str, ov: dict) -> None:
    """The registry's routes to the observable, each with its state for this level or transition, and "Plan it"
    where physim has the kind of experiment: it opens that kind's new-experiment form for this nucleus."""
    ui = ctx.ui
    title = nuclide_text(name)
    ui.label(f"Routes to {label_text(ov['label'])} in {title}").classes("ps-section").style("margin-top: 6px")
    quiet = []
    with ui.element("div").classes("ps-routes"):
        for r in ov["routes"]:
            status, off = r["status"], r.get("applies") is False
            if status == nucleus.NO_EXPERIMENT and not off:
                quiet.append(r)
                continue
            values = f"{r['independent']} value{'s' * (r['independent'] != 1)}"
            runs_ = f"{r['projected']} simulated run{'s' * (r['projected'] != 1)}"
            tag, colour = (("does not apply", "ps-tag-off") if off else {
                nucleus.DONE: (f"done · {values}" + (f", {runs_}" if r["projected"] else ""), "ps-tag-done"),
                nucleus.ONLY_PROJECTED: (f"projected only · {runs_}", "ps-tag-can"),
                nucleus.IN_PHYSIM: ("in physim", "ps-tag-can"),
                nucleus.LATER: ("later", "")}.get(status, (status, "")))
            with ui.element("div").classes("ps-route" + (" ps-route-later" if off or status == nucleus.LATER else "")):
                with ui.row().classes("w-full items-start no-wrap gap-2"):
                    ui.label(r["label"] if r["method"] != "other" else (r["kind_label"] or r["label"])).classes(
                        "text-sm font-semibold")
                    ui.space()
                    ui.label(tag).classes("ps-tag " + colour)
                how = [r["kind_label"] if r["method"] != "other" else "", r["note"]]
                if any(how):
                    ui.label(": ".join(x for x in how if x)).classes("ps-muted")
                if status == nucleus.ONLY_PROJECTED and not off:
                    ui.label("Simulated runs only: nothing is measured by this route yet.").classes("ps-muted")
                if status == nucleus.LATER:
                    ui.label("A kind of experiment physim does not have yet.").classes("ps-muted")
                if status == nucleus.NO_EXPERIMENT:
                    ui.label(MEANS_TEXT.get(r["means"], "")).classes("ps-muted")
                if r.get("reason"):
                    ui.label(r["reason"]).classes("ps-warn" if off else "")
                if r["live"] and r["kind"] and not off:
                    ui.button(f"Plan it for {title}" if r["open"] else f"Plan another for {title}",
                              icon="arrow_forward", on_click=lambda kind=r["kind"]: plan(ctx, name, kind)).props(
                        "flat dense no-caps").classes("self-start").tooltip(
                        f"Opens the new-experiment form for {r['kind_label']}, with {name} filled in")
    if quiet:
        done = [r for r in quiet if r["rows"]]
        words = [f"{r['label']}" + (f" ({r['rows']} value{'s' * (r['rows'] != 1)} here)" if r["rows"] else "")
                 for r in done + [r for r in quiet if not r["rows"]]]
        ui.label("Not experiments of the planner (the evaluation, a calculation, methods known from the literature: "
                 "their values come in as data points): " + " · ".join(words) + ".").classes("text-xs ps-muted")


def _gives(ctx, ov: dict) -> None:
    """What the reference value gives (the registry's conversions), the checks between observables when there are
    any, and what the value is used for."""
    ui = ctx.ui
    display = ctx.state["start"].get("unit")
    ref = ov["reference"]
    checks = ov.get("checks", []) or []
    used = uses(ov["observable"])
    if ref is None and not checks and not used:
        return
    _, body = _panel(ui, "What this value gives",
                     f"from the reference, {cell_text(ref, display)} ({ref['source']})" if ref else "")
    with body.classes("ps-panel-body"):
        with ui.column().classes("w-full gap-2"):
            if ref is None:
                ui.label("No reference value yet: nothing is converted.").classes("text-sm ps-muted")
            for g in ov["gives"]:
                to = g["to"]
                if g.get("key"):
                    try:
                        to = label_text(nucleus.label(to, g["key"]))
                    except ValueError:
                        pass
                with ui.column().classes("gap-0"):
                    if g["value"] is None:
                        ui.label(f"{to}: not worked out").classes("text-sm ps-muted")
                        ui.label(g["missing"]).classes("text-xs ps-muted")
                        continue
                    factor, unit = (display_factor("B(E2)", display) if g["unit"] == "e2fm4" else (1.0, g["unit"]))
                    ui.label(f"{to} = {shown(g['value'], g['plus'], g['minus'], unit, g['direction'], factor)}"
                             ).classes("text-sm ps-num font-medium")
                    ui.label(words(g["how"])).classes("text-xs ps-muted")
            if checks:
                ui.label("Checks against other observables").classes("ps-section").style("margin-top: 4px")
                for check in checks:
                    with ui.column().classes("gap-0"):
                        ui.label(": ".join(str(x) for x in (check.get("label"), check.get("text")) if x)).classes(
                            "text-sm")
                        if check.get("how"):
                            ui.label(str(check["how"])).classes("text-xs ps-muted")
            if used:
                ui.label("What it is used for").classes("ps-section").style("margin-top: 4px")
                for item in used:
                    with ui.column().classes("gap-0"):
                        ui.label(item["text"]).classes("text-sm").style("line-height: 1.5")
                        if item["reference"]:
                            ui.label(item["reference"]).classes("text-xs ps-muted")


def _experiments(ctx, view: dict) -> None:
    ui = ctx.ui
    found = view["experiments"]
    _, body = _panel(ui, "Experiments", f"{len(found)} with this nucleus; a row opens the experiment")
    with body:
        if not found:
            ui.label("No experiment with this nucleus yet.").classes("text-sm ps-muted p-3")
            return
        experiments_table(ctx, found, ctx.state["start"].get("unit"))


def _history(ctx, ov: dict) -> None:
    """Every row of the chosen observable that was ever written, the newest first: rows are never rewritten."""
    ui = ctx.ui
    rows = ov["history"]
    _, body = _panel(ui, "History", "every row is kept")
    words = {"superseded": "superseded by a newer row", "retracted": "withdrawn", "retraction": "the retraction"}
    with body:
        if not rows:
            ui.label("No row yet.").classes("text-sm ps-muted p-3")
            return
        with ui.element("table").classes("ps-nt"):
            for h in rows[:HISTORY_ROWS]:
                gone = h["status"] != "holds"
                with ui.element("tr"):
                    _td(ui, h["created"][:10], classes="ps-nt-num", colour="ps-muted")
                    _td(ui, h["source"], [f"{EVIDENCE[h['evidence']].lower()} · {h['method_label']}"
                                          + (f" · {h['detail']}" if h["detail"] else ""),
                                          words.get(h["status"], "") + (f": {h['reason']}" if h["reason"] else "")],
                        colour="ps-muted" if gone else "")
                    _td(ui, "" if h["status"] == "retraction" else h["text"].replace(
                        " " + h["unit"], " " + unit_text(h["unit"])), classes="ps-nt-num",
                        colour="ps-muted" if gone else "")
        if len(rows) > HISTORY_ROWS:
            ui.label(f"… and {len(rows) - HISTORY_ROWS} older rows.").classes("text-xs ps-muted p-2")


def _point_form(ctx, name: str, ov: Optional[dict], reload, redraw) -> None:
    """Add a data point: an observable of a level or transition with its value and uncertainty, the method, and a
    reference or "own measurement". What is wrong with it is said beside the form, all at once. The same form
    corrects a point (a new row that supersedes it)."""
    ui = ctx.ui
    st = ctx.state["start"]
    if not st.get("point"):
        st["point"] = blank_point(*((ov["observable"], ov["key"]) if ov is not None else ("B(E2)", "")))
    f = st["point"]

    def set_(key: str, refresh: bool = False):
        def change(e) -> None:
            f[key] = e.value
            if key == "observable":  # the unit and the direction are the observable's
                o = observables.get(e.value)
                f.update(unit=o.unit, direction=o.direction)
            if key == "unit" and e.value == observables.WEISSKOPF:
                f["direction"] = observables.DOWN
            if key == "evidence" and f.get("method") not in nucleus.POINT_METHODS[e.value]:
                f["method"] = None
            if refresh:
                form.refresh()
        return change

    def field(label: str, key: str, placeholder: str = ""):
        return ui.input(label, value=f.get(key) or "", placeholder=placeholder, on_change=set_(key)).props(
            "dense outlined").classes("w-full")

    def submit() -> None:
        try:
            row = add_point(name, f)
        except nucleus.PointProblem as e:
            st["point_problems"] = list(e.problems)
            problems.refresh()
            return
        except ValueError as e:
            st["point_problems"] = [str(e)]
            problems.refresh()
            return
        if st.get("reference") == f.get("corrects"):
            st["reference"] = None
        st.update(point=None, point_problems=[], observable=row.observable, key=row.key)
        ui.notify(f"{'Corrected' if row.supersedes else 'Added'}: {label_text(nucleus.label(row.observable, row.key))} "
                  f"of {nuclide_text(name)}, {nucleus.source(row)}.", type="positive")
        reload()

    def cancel() -> None:
        st.update(point=None, point_problems=[])
        redraw()

    @ui.refreshable
    def form() -> None:
        o = observables.get(f["observable"])
        with ui.grid(columns=2).classes("w-full gap-2"):
            ui.select({k: v.symbol for k, v in observables.OBSERVABLES.items()}, value=f["observable"],
                      label="Observable", on_change=set_("observable", True)).props("dense outlined").classes("w-full")
            field("Transition" if o.belongs_to == "transition" else "Level", "key",
                  "2+1 > 0+1" if o.belongs_to == "transition" else "2+1, or 2+ at 1173 keV")
        with ui.grid(columns=3 if o.direction else 2).classes("w-full gap-2"):
            field("Value", "value")
            ui.select({u: unit_text(u) for u in o.units}, value=f["unit"], label="Unit",
                      on_change=set_("unit", True)).props("dense outlined").classes("w-full")
            if o.direction:
                ui.select({observables.UP: "↑ excitation (lower → upper)", observables.DOWN: "↓ decay"},
                          value=f["direction"], label="Direction", on_change=set_("direction")).props(
                    "dense outlined").classes("w-full")
        with ui.row().classes("w-full items-center gap-2"):
            ui.label("Uncertainty").classes("text-xs ps-muted")
            ui.toggle({"total": "One total", "stat-sys": "Stat and sys", "asymmetric": "Plus and minus"},
                      value=f["form"], on_change=set_("form", True)).props("dense no-caps unelevated")
        with ui.grid(columns=1 if f["form"] == "total" else 2).classes("w-full gap-2"):
            for key, label in {"total": (("uncertainty", "Uncertainty (total)"),),
                               "stat-sys": (("stat", "Statistical"), ("sys", "Systematic")),
                               "asymmetric": (("plus", "Plus"), ("minus", "Minus"))}[f["form"]]:
                field(label, key)
        with ui.grid(columns=2).classes("w-full gap-2"):
            ui.select({k: EVIDENCE[k] for k in nucleus.POINT_EVIDENCE}, value=f["evidence"], label="Evidence",
                      on_change=set_("evidence", True)).props("dense outlined").classes("w-full")
            ui.select({m: METHODS[m] for m in nucleus.POINT_METHODS[f["evidence"]]}, value=f["method"],
                      label="Method", on_change=set_("method")).props("dense outlined").classes("w-full")
        field("Detail of the method (a line, optional)", "detail", "thin target, GOSIA")
        with ui.row().classes("w-full items-center gap-2"):
            ui.label("Where it is from").classes("text-xs ps-muted")
            ui.toggle({"reference": "A reference", "own": "Own measurement"}, value=f["from"],
                      on_change=set_("from", True)).props("dense no-caps unelevated")
        with ui.grid(columns=2).classes("w-full gap-2"):
            if f["from"] == "own":
                field("Date", "date", "2027-06-14")
                field("Note: where and how", "note", "iThemba LABS, runs 12 to 31")
            else:
                field("Authors", "authors", "Allmond et al.")
                field("Year", "year", "2014")
                field("DOI", "doi", "10.1103/PhysRevC.90.034309")
                field("Where in the paper", "where", "Table II")

    @ui.refreshable
    def problems() -> None:
        found = st.get("point_problems") or []
        if not found:
            return
        with ui.column().classes("ps-warnbox rounded p-2 gap-1").style("flex: 1 1 200px; min-width: 0"):
            ui.label("Not a data point yet:").classes("text-sm font-semibold ps-warn")
            for text in found:
                ui.label("· " + text).classes("text-sm")

    correcting = f.get("corrects")
    head, body = _panel(ui, "Correct a data point" if correcting else "Add a data point",
                        f"row {correcting[:8]}… stays in the history, superseded" if correcting
                        else "a value with its uncertainty, its method and where it is from")
    head.classes("ps-point-form")
    with body.classes("ps-panel-body"):
        with ui.row().classes("w-full items-start gap-3"):
            with ui.column().classes("gap-2").style("flex: 2 1 300px; min-width: 0"):
                form()
            problems()
        with ui.row().classes("w-full items-center gap-2").style("margin-top: 8px"):
            ui.label("A point is independent evidence. It is never changed: a correction is a new row that "
                     "supersedes it.").classes("text-xs ps-muted").style("flex: 1 1 200px")
            if correcting:
                ui.button("Cancel", color=None, on_click=cancel).props("flat no-caps").classes("ps-muted")
            ui.button("Save the correction" if correcting else "Add the point", icon="add",
                      on_click=submit).props("unelevated no-caps")
