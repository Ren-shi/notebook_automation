"""What the page of one nucleus (:mod:`.tab_nucleus`), the comparison of one observable (:mod:`.tab_compare`) and the
nuclei screen (:mod:`.start`) share: the words, the reading and writing, going from one view to another, the tables,
the import, the report and the form of a data point.

Nothing is computed here: every number comes from :mod:`physim.nuclear.nucleus` (:func:`~physim.nuclear.nucleus.view`,
:func:`~physim.nuclear.nucleus.observable_view`), read once on a worker thread when a view opens (:func:`load`).
The unit switch of a B(E2) only changes how the canonical value is written, through the registry's units
(:func:`display_factor`).
"""

from __future__ import annotations

import math
import re
import threading
from typing import Optional
from urllib.parse import parse_qs, quote, urlsplit

from .. import ensdf, nucleus, observables, result_rows, results, runs
from ..results import EVIDENCE, METHODS

#: The tables, tags and cards of the nuclei screen and the nucleus page, in the page's colours (both themes).
NUCLEUS_STYLE = """
.ps-nt { width: 100%; border-collapse: collapse; font-size: 13px; }
.ps-nt th { text-align: left; font-weight: 500; color: var(--ps-muted); padding: 6px 10px; white-space: nowrap;
  border-bottom: 1px solid var(--ps-line); }
.ps-nt td { padding: 6px 10px; border-bottom: 1px solid var(--ps-hair); vertical-align: top; }
.ps-nt tr:last-child td { border-bottom: 0; }
.ps-nt th.ps-nt-wrap { white-space: normal; }
.ps-nt .ps-nt-keep .ps-nt-sub { white-space: nowrap; }
.ps-nt-box { overflow-x: auto; }
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
.ps-nuc-body > div { padding: 10px 16px 12px; min-width: 0; overflow-x: auto; }
.ps-nuc-body > div + div { border-left: 1px solid var(--ps-hair); }
.ps-nuc-h { font-size: 11px; color: var(--ps-muted); letter-spacing: .04em; margin-bottom: 4px; }
.ps-kv { display: grid; grid-template-columns: auto 1fr; gap: 3px 14px; font-size: 13px; }
.ps-kv > :nth-child(odd) { color: var(--ps-muted); }
@media (max-width: 1000px) { .ps-two, .ps-nuc-body { grid-template-columns: 1fr; }
  .ps-nuc-body > div + div { border-left: 0; border-top: 1px solid var(--ps-hair); } }
"""

_SUPERSCRIPT = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
#: The units as the page writes them (the ledger's own spelling otherwise).
UNIT_TEXT = {"e2fm4": "e²fm⁴", "e2b2": "e²b²", "efm2": "e fm²", "eb": "e b", "e2fm2": "e²fm²", "mu_N": "μ_N",
             "mu_N2": "μ_N²"}
#: The units the switch of the nucleus page offers for a B(E2), of the registry's: display only.
BE2_UNITS = ("e2fm4", "e2b2")
#: How many headline rows a card or the readouts show.
HEADLINES = 4


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


def value_parts(cell: dict, display: Optional[str] = None) -> tuple:
    """(number, unit) of a :func:`physim.nuclear.nucleus.cell` as the page writes them apart, in the unit the page
    shows: ("884 ± 13", "e²fm⁴↑"). The same words as :func:`cell_text`."""
    factor, unit = display_factor(cell["observable"], display)
    if cell.get("problem"):  # not in the canonical convention: as given
        factor, unit = 1.0, cell["unit"]

    def scaled(x):
        return None if x is None else x * factor

    arrow = "↓" if unit == observables.WEISSKOPF else observables.DIRECTIONS.get(cell["direction"], "")
    return nucleus.value_text(cell["value"] * factor, scaled(cell["plus"]), scaled(cell["minus"])), \
        unit_text(unit) + arrow


def standing(row: dict) -> Optional[dict]:
    """The value a row of the observables table stands for: the evaluated one, else the best independent one
    (measured, from an outside program, theory, in that order), as a cell; None without any."""
    return row["evaluated"] or next((row[k]["best"] for k in ("measured", "outside_program", "theory")
                                     if row[k] is not None), None)


def range_text(summary: Optional[dict], display: Optional[str] = None, word: str = "measured") -> str:
    """How many values of one kind of evidence a row holds and their range, in the unit the page shows: "16
    measured, 1.313 to 2.975 ps"; for one value its method ("1 measured · Moment measurement")."""
    if summary is None:
        return f"No {word} value"
    best = summary["best"]
    if summary["count"] == 1:
        return f"1 {word} · {method_words(best)}"
    factor, unit = display_factor(best["observable"], display)
    return (f"{summary['count']} {word}, {summary['min'] * factor:.4g} to {summary['max'] * factor:.4g} "
            f"{unit_text(unit)}")


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
        return f"{via['method_label']}, from {words(label_text(via['text']))}"
    return cell["method_label"]


def detail_words(cell: dict) -> str:
    """The line beside the method; for a converted value, with what it was converted (the method says what it was
    made from already)."""
    made = re.match(r"from .*? (?=with )", cell["detail"]) if cell.get("via") else None
    return words(cell["detail"][made.end():] if made else cell["detail"])


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
    """The address of a nucleus's page in the app, for a link from an experiment: "/?nucleus=62Ni"; with the
    observable and the level or transition, the address of their comparison. The key is percent-encoded: a bare "+"
    in an address reads as a space."""
    out = "/?nucleus=" + quote(results.nuclide_name(nuclide), safe="")
    for name, value in (("observable", observable), ("key", key)):
        if value:
            out += f"&{name}=" + quote(str(value), safe="")
    return out


def named_in(address: str) -> Optional[tuple]:
    """(nuclide, observable, key) of an address that :func:`url` wrote (the last two None where it names none);
    None for an address that names no nucleus."""
    given = parse_qs(urlsplit(str(address)).query)
    nuclide, observable, key = ((given.get(name) or [None])[0] for name in ("nucleus", "observable", "key"))
    return (nuclide, observable, key) if nuclide else None


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


#: The folders of nuclei whose repeated data points were withdrawn since the app started.
_repeats_withdrawn: set = set()
_repeats_lock = threading.Lock()


def withdraw_repeats() -> None:
    """Once per start of the app, before the nuclei are first read: the data points that repeat an earlier one are
    withdrawn (:func:`physim.nuclear.nucleus.retract_duplicates`: nothing is deleted). A ledger that cannot be read
    or written is left as it is."""
    home = str(results.nuclei_home())
    with _repeats_lock:  # two pages opened at once: the second waits for the first, and finds it done
        if home in _repeats_withdrawn:
            return
        _repeats_withdrawn.add(home)
        try:
            nucleus.retract_duplicates()
        except Exception:  # noqa: BLE001 -- a repeated point shows twice, as it did
            pass


def load_nuclei(fill: bool = True) -> dict:
    """What the nuclei screen shows, read once: :func:`physim.nuclear.nucleus.nuclei` after the backfill and
    :func:`withdraw_repeats`, and when each experiment last changed ({"nuclei", "unassigned", "modified": {path:
    time}}). Without ``fill`` nothing is written: neither of the two."""
    if fill:
        withdraw_repeats()
    listed = backfill() if fill else runs.list_experiments()
    return {**nucleus.nuclei(), "modified": {e["path"]: e["modified"] for e in listed}}


def _as_typed(key: Optional[str]) -> Optional[str]:
    """A key from an address: a "+" that was not percent-encoded arrives as a space ("62Ni:2 1>0 1:E2")."""
    return None if not key else re.sub(r"(?<=\d) (?=[\d(@])", "+", str(key))


def evaluation_taken(name: str) -> bool:
    """Whether the nucleus's ledger holds the evaluation already (a read of its own small file, before the view)."""
    points = results.read(results.nucleus_folder(name) / results.POINTS_FILE)
    return any(nucleus.is_evaluation(r) for r in results.current(points))


def load(nuclide: str, observable: Optional[str] = None, key: Optional[str] = None,
         reference: Optional[str] = None, first: bool = True, compare: bool = True) -> dict:
    """Everything the views of one nucleus show, read once: {"nuclide", "view", "observable", "notes", "rows"}.

    "view" is :func:`physim.nuclear.nucleus.view` and "observable" the
    :func:`~physim.nuclear.nucleus.observable_view` of the row asked for (the first row of the table when none is,
    None for a nucleus without rows, and with ``compare`` off: the nucleus's own page compares nothing). If the
    evaluation was never taken for the nucleus it is taken first
    (:func:`~physim.nuclear.nucleus.snapshot_evaluation`). "notes" says what could not be done: no local ENSDF
    copy, a row or a reference that is not there. ``first`` is False when the page reads again after a change of
    its own (the experiments from before the ledger have their rows by then, and repeated points are withdrawn).
    Raises ``ValueError`` for a name that is no nuclide."""
    name = results.nuclide_name(nuclide)
    if first:
        withdraw_repeats()
        backfill()
    notes = []
    if not evaluation_taken(name):
        if not ensdf.available():
            notes.append("There is no local copy of ENSDF, so the adopted values are not here to compare with "
                         "(scripts/fetch_ensdf.py fetches it).")
        else:
            try:
                taken = nucleus.snapshot_evaluation(name)
                if not taken["added"] and not taken["unchanged"]:
                    notes.append(f"The local copy of ENSDF has no adopted value of {name}.")
            except (ensdf.EnsdfMissing, ValueError, OSError) as e:
                notes.append(f"The adopted values of ENSDF could not be taken: {e}")
    rows = results.of_nuclide(name)  # the ledgers, read once for the view and the observable shown
    view = nucleus.view(name, rows=rows)
    return {"nuclide": name, "view": view, "notes": notes, "rows": rows,
            "observable": load_observable(name, view, observable, key, reference, notes, rows) if compare else None}


def load_landing(nuclide: str, first: bool = True) -> dict:
    """What the landing view of a nucleus shows, read once: :func:`load` without an observable, and with it
    "levels" (:func:`physim.nuclear.nucleus.lowest_levels`), "stability" (the words of
    :func:`~physim.nuclear.nucleus.stability`: "stable, 3.63 % of natural Ni") and "element" ("Nickel-62", or ""
    when the table of elements has no name for it). The routes still open are read apart (:func:`open_routes`):
    their tests take longer."""
    from .. import data

    held = load(nuclide, first=first, compare=False)
    name, view = held["nuclide"], held["view"]
    try:
        element = f"{data.element(view['z']).name}-{view['a']}"
    except (KeyError, ValueError, IndexError, AttributeError):
        element = ""
    try:
        stable = nucleus.stability(name)["words"]
    except (KeyError, ValueError, OSError):
        stable = ""
    return {**held, "levels": nucleus.lowest_levels(name, rows=held["rows"]), "stability": stable,
            "element": element}


def open_routes(nuclide: str, view: dict, rows: list) -> list:
    """The routes physim can plan for the nucleus and that no run or measurement of one's own has taken yet, one
    per method, from the routes of the view's headline rows (:func:`physim.nuclear.nucleus.routes`, with physim's
    tests). Each {"method", "label", "kind", "kind_label", "note", "applies", "reason", "independent", "own"} as
    the route of the first row that has it says ("for": that row's label), and "gives": the rows it gives, each
    {"observable", "key", "label"}. ``rows`` are the nucleus's rows as :func:`load` read them."""
    name = results.nuclide_name(nuclide)
    now, via = nucleus.with_converted(name, results.current([r for r in rows if r.nuclide == name]))
    found: dict = {}
    for t in view["table"]:
        if not t["headline"] or not t["open"]:
            continue
        gives = {"observable": t["observable"], "key": t["key"], "label": t["label"]}
        for r in nucleus.routes(name, t["observable"], t["key"], results.select(now, t["observable"], t["key"]),
                                via, every=now):
            if not r["open"]:
                continue
            if r["method"] in found:
                found[r["method"]]["gives"].append(gives)
                continue
            found[r["method"]] = {**{k: r[k] for k in ("method", "label", "kind", "kind_label", "note", "applies",
                                                       "reason", "independent", "own")},
                                  "for": t["label"], "gives": [gives]}
    return list(found.values())


def load_observable(name: str, view: dict, observable: Optional[str] = None, key: Optional[str] = None,
                    reference: Optional[str] = None, notes: Optional[list] = None,
                    rows: Optional[list] = None) -> Optional[dict]:
    """The :func:`physim.nuclear.nucleus.observable_view` of the row asked for, or of the first row of the
    nucleus's table. What is asked for and is not there is said in ``notes`` and the nearest thing is shown.
    ``rows`` are the nucleus's rows as :func:`load` read them (read again when left out)."""
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
        return nucleus.observable_view(name, *chosen, reference=reference, rows=rows)
    except ValueError as e:
        if reference is None:
            notes.append(str(e))
            return None
        notes.append(f"The reference asked for is not there ({e}); the evaluated value is the reference again.")
        return nucleus.observable_view(name, *chosen, rows=rows)


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


def import_file(nuclide: str, filename: str, data: bytes, home=None) -> dict:
    """The data points of a file sent to the page, added to the nucleus's ledger
    (:func:`physim.nuclear.nucleus.import_points`, which says what was added, already there and refused). The file
    is put in a temporary folder under its own name for the reader, and removed after."""
    import tempfile
    from pathlib import Path

    suffix = Path(str(filename or "points.json")).suffix.lower() or ".json"
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        path = Path(tmp) / f"points{suffix}"
        path.write_bytes(data)
        found = nucleus.import_points(path, nuclide, home)
    found["file"] = str(filename or "")
    return found


def import_lines(found: dict) -> dict:
    """The three lists of an import as the page shows them: {"added", "already", "refused"}, each a list of
    lines ("quantities.BE2_up[0] · B(E2) of 2+1 > gs: 0.0906 ± 0.0037 e2b2↑, coulex-fit, J. M. Allmond et al. 2014
    → Coulex fitted to all orders")."""
    def head(x: dict) -> str:
        return f"{x['where']} · {words(x['what'])}" if x["what"] else x["where"]

    return {"added": [f"{head(x)} → {METHODS[x['row'].method]}" for x in found["added"]],
            "already": [f"{head(x)} (row {x['id'][:8]})" for x in found["already"]],
            "refused": [f"{head(x)}: " + "; ".join(x["reasons"]) for x in found["refused"]]}


def report_html(nuclide: str) -> str:
    """The report of the nucleus (:func:`physim.nuclear.logbook.nucleus_html`): nothing is written."""
    from .. import logbook

    return logbook.nucleus_html(nuclide)


def add_point(nuclide: str, form: dict, home=None):
    """Add the data point the form holds to the nucleus's ledger (a correction when the form names the row it
    replaces) and return its row. Raises :class:`physim.nuclear.nucleus.PointProblem` with everything that is
    wrong: the page shows the messages beside the form."""
    args = point_arguments(form)
    return nucleus.add_data_point(nuclide, args.pop("observable"), args.pop("key"), args.pop("value"),
                                  args.pop("unit"), home=home, **args)


# ---------------------------------------------------------------------------------------------------------------
# Going from one view to another


def add_style(ctx) -> None:
    """The style of the nuclei screen and of the views of a nucleus, once per browser page."""
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


def show_from_top(ctx) -> None:
    """Build the start page again and show it from its top. The client is taken first: the button that asked for
    it is gone with the page it was on."""
    client = ctx.ui.context.client
    ctx.refresh_start()
    client.run_javascript("window.scrollTo({top: 0})")


def open_nuclei(ctx) -> None:
    """Show the nuclei screen."""
    ctx.state["start"].update(view="nuclei")
    show_from_top(ctx)


def open_nucleus(ctx, nuclide: str, observable: Optional[str] = None, key: Optional[str] = None) -> None:
    """Show the page of one nucleus; with an observable, its comparison (:mod:`.tab_compare`)."""
    try:
        name = results.nuclide_name(str(nuclide).strip())
    except (ValueError, KeyError, TypeError, IndexError):
        ctx.ui.notify(f"'{nuclide}' is not a nuclide: write it like 62Ni or Ni-62.", type="warning")
        return
    ctx.state["start"].update(view="observable" if observable else "nucleus", nucleus=name, observable=observable,
                              key=key, reference=None, point=None, point_problems=[], to_point_form=False)
    show_from_top(ctx)


def open_observable(ctx, nuclide: str, observable: str, key: Optional[str] = None, to_point_form: bool = False
                    ) -> None:
    """Show the comparison of one observable of one level or transition (of the first one in the nucleus's table
    without a ``key``), at the form of a data point if asked."""
    ctx.state["start"].update(view="observable", nucleus=nuclide, observable=observable, key=key, reference=None,
                              point=None, point_problems=[], to_point_form=to_point_form)
    show_from_top(ctx)


def open_experiment(ctx, path: str):
    """Open an experiment, as the list of experiments always did. The start page is the nuclei again when it is
    shown next."""
    ctx.state["start"].update(view="nuclei")
    return ctx.open_experiment(path)


def close(ctx) -> None:
    """Back to the open experiment; the start page is the nuclei again when it is shown next."""
    ctx.state["start"].update(view="nuclei")
    ctx.close_start()


def plan(ctx, nuclide: Optional[str], kind: Optional[str] = None) -> None:
    """Today's new-experiment form for a kind of experiment (the kinds to pick from when none is named), with the
    nucleus filled in where the form has a field for it (:func:`.start.kind_form`). Nothing is computed."""
    st = ctx.state["start"]
    st.update(view="kinds", kind=kind, for_nucleus=nuclide)
    st.pop("values", None)
    show_from_top(ctx)


# ---------------------------------------------------------------------------------------------------------------
# Small parts of the screens


def table_head(ui, head: tuple, numeric: tuple = ()):
    """A table with the column heads ``head`` (those at the places ``numeric`` right-aligned, as numbers are); the
    rows are built inside it."""
    table = ui.element("table").classes("ps-nt")
    with table:
        with ui.element("tr"):
            for i, text in enumerate(head):
                with ui.element("th").classes(("ps-nt-num" if i in numeric else "")
                                              + (" ps-nt-wrap" if len(text) > 14 else "")):
                    ui.label(text)
    return table


def td(ui, text: str = "", sub=(), classes: str = "", colour: str = ""):
    """A cell: its text and, under it, small lines."""
    cell = ui.element("td").classes(classes)
    with cell:
        if text:
            ui.label(text).classes(colour)
        for line in ([sub] if isinstance(sub, str) else sub):
            if line:
                ui.label(line).classes("ps-nt-sub")
    return cell


def group_row(ui, text: str, columns: int, colour: str = "") -> None:
    with ui.element("tr").classes("ps-nt-group"):
        with ui.element("td").props(f"colspan={columns}"):
            ui.label(text).classes(colour)


def panel(ui, title: str, sub: str = ""):
    """A bordered panel with its head line; (the head, to add to; the body, to build in)."""
    with ui.element("div").classes("ps-panel"):
        head = ui.element("div").classes("ps-panel-head")
        with head:
            ui.label(title)
            if sub:
                ui.label(sub).classes("ps-muted text-sm")
        return head, ui.element("div").classes("ps-nt-box")


def latest_text(x: dict, display: Optional[str] = None) -> tuple:
    """(text, small line) of an experiment's latest result: its value and, in words, what kind of evidence it is
    and by which method; "not analysed" or "no run" without one."""
    latest = x.get("latest")
    if latest is None:
        return ("no run" if not x["runs"] else "not analysed" if not x["analysed_runs"] else "no result kept"), ""
    return (cell_text(latest, display),
            f"{label_text(latest['label'])} · {EVIDENCE[latest['evidence']].lower()} · {method_words(latest)}")


def experiments_table(ctx, experiments: list, display: Optional[str] = None, opener: bool = False) -> None:
    """The experiments of a nucleus: name, kind, runs and the latest result or "not analysed". A row opens the
    experiment, as the list of experiments always did; with ``opener`` each row ends in "Open →" to say so."""
    ui = ctx.ui
    names = [x["name"] for x in experiments]
    with table_head(ui, ("Experiment", "Kind", "Runs", "Latest result") + (("",) if opener else ()), numeric=(2, 3)):
        for x in experiments:
            row = ui.element("tr").classes("ps-nt-click")
            row.on("click", lambda path=x["path"]: open_experiment(ctx, path))
            with row:
                td(ui, x["name"], [x["folder"] if names.count(x["name"]) > 1 else "",
                                    {"studied": "", "partner": "as the partner"}.get(x["role"], "results only")],
                    colour="ps-accent")
                with ui.element("td"):
                    ui.label(x["kind_label"] or "").classes("ps-tag")
                td(ui, f"{x['runs']} run{'s' * (x['runs'] != 1)}", classes="ps-nt-num")
                text, sub = latest_text(x, display)
                td(ui, text, sub, classes="ps-nt-num", colour="" if x.get("latest") else "ps-muted")
                if opener:
                    td(ui, "Open →", classes="ps-nt-open").props(f'aria-label="Open {x["name"].replace(chr(34), "")}"')


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
# Parts with a form or a file: the import, the report, the data point


def import_button(ctx, name: str, reload, label: str = "Import points from a file") -> None:
    """"Import points from a file" (or ``label``): a window that takes a JSON or CSV file of data points (or a
    literature file of the paper verifications), adds them through :func:`import_file` and shows the three lists:
    added, already there, refused with the reasons. The page reads again when the window closes after an import
    that added something."""
    ui = ctx.ui
    st = ctx.state["start"]

    @ui.refreshable
    def lists() -> None:
        found = st.get("imported")
        if not found:
            return
        said = import_lines(found)
        ui.label(f"{found['file']}: {len(said['added'])} added, {len(said['already'])} already there, "
                 f"{len(said['refused'])} refused").classes("text-sm font-semibold")
        for key, title, colour in (("added", "Added", "ps-ok"), ("already", "Already there: nothing added", ""),
                                   ("refused", "Refused", "ps-warn")):
            if not said[key]:
                continue
            ui.label(f"{title} ({len(said[key])})").classes(f"text-xs font-semibold {colour}").style(
                "margin-top: 6px")
            with ui.column().classes(f"gap-0 ps-import-{key}"):
                for line in said[key]:
                    ui.label(line).classes("text-xs").style("line-height: 1.45")

    async def got(e) -> None:
        from nicegui import run

        try:
            data = await e.file.read()
            found = await run.io_bound(import_file, name, e.file.name, data)
        except (ValueError, OSError) as err:
            ui.notify(f"The file could not be read: {err}", type="negative", multi_line=True)
            return
        st["imported"] = found
        st["import_added"] = st.get("import_added", False) or bool(found["added"])
        lists.refresh()
        upload.reset()

    def closed() -> None:
        """The page reads again when the window closes, if the import added a point (reading again while it is
        open would rebuild the page and the window with it)."""
        st["imported"] = None
        if st.pop("import_added", False):
            reload()

    with ui.dialog().on("hide", closed) as dialog, ui.card().classes("ps-card ps-import").style(
            "width: min(820px, 95vw); max-width: 95vw"):
        with ui.row().classes("w-full items-center"):
            ui.label(f"Import points of {nuclide_text(name)} from a file").classes("text-base font-semibold")
            ui.space()
            ui.button(icon="close", color=None, on_click=dialog.close).props("flat round dense").classes("ps-muted")
        ui.label("JSON or CSV, one data point each with the fields of the form (observable, level or transition, "
                 "value, unit, direction, uncertainty, method, evidence, reference or own measurement), or a "
                 "literature file of a paper verification (nucleus, level, quantities). Every point is tried; one "
                 "that is in the ledger already is not added again, and one that cannot be taken is listed with "
                 "every reason.").classes("text-xs ps-muted").style("line-height: 1.5")
        upload = ui.upload(auto_upload=True, on_upload=got, label="A .json or .csv file").props(
            "accept=.json,.csv max-files=1 flat bordered").classes("w-full")
        with ui.column().classes("w-full gap-0").style("max-height: 55vh; overflow-y: auto"):
            lists()
    ui.button(label, icon="upload_file", on_click=dialog.open).props("outline no-caps")


def report_button(ctx, name: str, label: str = "Report for this nucleus") -> None:
    """"Report for this nucleus" (or ``label``): the report of the nucleus across its experiments
    (:func:`physim.nuclear.logbook.nucleus_html`), saved as HTML or as PDF from its menu, as the Report step of an
    experiment offers its report. Building it writes nothing."""
    ui = ctx.ui
    made: dict = {}

    async def page() -> str:
        from nicegui import run

        note = ui.notification("Building the report of the nucleus…", spinner=True, timeout=None)
        try:
            made["page"] = await run.io_bound(report_html, name)
        finally:
            note.dismiss()
        return made["page"]

    async def save_html() -> None:
        ui.download.content((await page()).encode("utf-8"), f"{name}-nucleus-report.html")

    async def save_pdf() -> None:
        from nicegui import run

        from .. import logbook

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
        ui.download.content(data, f"{name}-nucleus-report.pdf")

    with ui.button(label, icon="description").props("outline no-caps").tooltip(
            "Every result of this nucleus across its experiments, results first: as HTML or PDF"):
        with ui.menu():
            ui.menu_item("Save as HTML", on_click=save_html)
            ui.menu_item("Save as PDF", on_click=save_pdf)


def scroll_to_point_form(ui) -> None:
    """Scroll to the form of a data point where it stands on the page itself: on the landing page of a nucleus
    that has no row yet, beside its experiments. Everywhere else the form is in the window of the comparison
    screen, which :func:`open_observable` opens with ``to_point_form``."""
    ui.run_javascript(
        "document.querySelector('.ps-point-form')?.scrollIntoView({behavior: 'smooth', block: 'start'})")


def point_form(ctx, name: str, ov: Optional[dict], reload, redraw) -> None:
    """Add a data point: an observable of a level or transition with its value and uncertainty, the method, and a
    reference or "own measurement". What is wrong with it is said beside the form, all at once. The same form
    corrects a point (a new row that supersedes it)."""
    ui = ctx.ui
    st = ctx.state["start"]
    blank = blank_point(*((ov["observable"], ov["key"]) if ov is not None else ("B(E2)", "")))
    if not st.get("point") or st["point"] == st.get("point_blank"):  # nothing typed yet: it follows the row shown
        st.update(point=dict(blank), point_blank=blank)
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
    head, body = panel(ui, "Correct a data point" if correcting else "Add a data point",
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
