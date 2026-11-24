"""The results tools: the nuclei and what is known about each, data points, and the comparison of the routes to
one observable.

Each tool is a method of :class:`ResultTools`; its docstring is what the model reads. The tools are thin: the
nucleus layer works everything out (:mod:`physim.nuclear.nucleus`, on the results ledger of
:mod:`physim.nuclear.results`) and these turn its rows into short text with units ("884 ± 13 e2fm4↑"). They act
on the nuclei, not on the open experiment, so none of them needs one.

The experiments are those of the session's folder. A session started on an experiments folder of its own keeps
the nuclei's own rows in the folder ``nuclei`` beside it, whatever ``$PHYSIM_NUCLEI`` says; a session on the
default experiments folder uses the nuclei home (``$PHYSIM_NUCLEI``, or ``~/.physim/nuclei``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

from ..nuclear import compare, observables, results
from ..nuclear import nucleus as _nucleus
from .session import Session
from .shape import plain

#: The observables of the registry, as the tools name them.
OBSERVABLES = tuple(observables.OBSERVABLES)

#: The most rows one answer lists; the rest are counted, and a narrower question shows them.
ROWS_SHOWN = 40

GROUPS = ("Evidence is in two groups that are never mixed. Independent: measured, evaluated, theory, outside "
          "program. Projected: a simulated run analysed back; it shows the precision a setup would reach, is "
          "closed against the value its simulation was given, and is never pulled against a measurement.")


def _percent(x: Optional[float]) -> Optional[str]:
    return None if x is None else f"±{x:.2g} %"


def cell_view(c: dict) -> dict:
    """One row of :func:`physim.nuclear.nucleus.cell` as short text."""
    out: dict[str, Any] = {"id": c["id"], "of": _nucleus.label(c["observable"], c["key"]),
                           "evidence": c["evidence"], "method": c["method_label"], "value": c["text"],
                           "source": c["source"]}
    if c["truth"] is not None:
        out["simulation_was_given"] = f"{c['truth']:.6g} {c['unit']}"
    if c["via"] is not None:
        out["method"] = f"{c['via']['method_label']}, converted"
        out["made_from"] = f"row {c['via']['id']}: {c['via']['text']}"
    given = c["as_given"]
    if c["converted"]:
        out["as_given"] = (f"{given['value']:.6g} {given['unit']}"
                           + {"up": "↑", "down": "↓"}.get(given["direction"], "") + f" ({c['converted']})")
    if c["detail"]:
        out["detail"] = c["detail"]
    if c["problem"]:
        out["problem"] = c["problem"]
    return out


def table_line(t: dict) -> str:
    """One row of :func:`physim.nuclear.nucleus.table` as a line: what each kind of evidence says."""
    parts = []
    if t["evaluated"] is not None:
        parts.append(f"evaluated {t['evaluated']['text']}")
    for name, words in (("measured", "measured"), ("outside_program", "from an outside program"),
                        ("theory", "from theory")):
        s = t[name]
        if s is None:
            continue
        if s["count"] == 1:
            parts.append(f"1 {words}: {s['best']['text']} ({s['best']['source']})")
        else:
            best = _percent(s["best"]["uncertainty_percent"])
            parts.append(f"{s['count']} {words}, {s['min']:.4g} to {s['max']:.4g} {t['unit']}"
                         + (f", the best to {best}" if best else ""))
    if t["projected"] is not None:
        reach = _percent(t["projected"]["best"]["uncertainty_percent"])
        parts.append(f"projected by {t['projected']['count']} simulated run"
                     + ("s" if t["projected"]["count"] != 1 else "") + (f": reaches {reach}" if reach else ""))
    return f"{t['label']} [{t['key']}]: " + ("; ".join(parts) if parts else "no row yet")


def nucleus_view(e: dict) -> dict:
    """One nucleus of :func:`physim.nuclear.nucleus.nuclei` (or :func:`~physim.nuclear.nucleus.view`) as text."""
    experiments = []
    for x in e["experiments"]:
        line = (f"{x['name']} [{x['folder']}], {x['kind_label']}, {x['role']}: {x['runs']} run"
                + ("s" if x["runs"] != 1 else "") + f", {x['analysed_runs']} analysed")
        if x["latest"] is not None:
            line += f"; latest {x['latest']['label']} = {x['latest']['text']} ({x['latest']['evidence']})"
        experiments.append(line)
    headline = [t for t in e["table"] if t["headline"]]
    rest: dict = {}
    for t in e["table"]:
        if not t["headline"]:
            rest[t["observable"]] = rest.get(t["observable"], 0) + 1
    ev = e["evaluation"]
    out: dict[str, Any] = {
        "nucleus": e["nuclide"], "experiments": experiments, "rows_that_hold": e["rows"],
        "data_points": e["points"],
        "evaluation": (f"{ev['rows']} rows from {ev['reference']} (data set of {ev['dataset_date']})" if ev["rows"]
                       else "not taken (snapshot_evaluation takes ENSDF's adopted values)"),
        "observables": [table_line(t) for t in headline]}
    if rest:
        out["evaluated_only"] = (", ".join(f"{n} of {k}" for k, n in rest.items())
                                 + " (list_results with the observable)")
    out["open"] = [o["text"] for o in e["open"]]
    return out


def comparison_view(v: dict) -> dict:
    """:func:`physim.nuclear.nucleus.observable_view` as text: each independent row against the reference, each
    projected row against the value its simulation was given, the routes with whether each applies to this
    nucleus, what the reference gives, the checks against other observables and what a value is used for."""
    c, cells = v["comparison"], v["cells"]
    out: dict[str, Any] = {"nucleus": v["nuclide"], "observable": v["label"], "key": v["key"],
                           "convention": v["convention"]}
    if c is None:
        out["note"] = "No row of this observable holds for this level or transition yet."
    else:
        ref = c.reference
        out["reference"] = (None if ref is None else
                            f"{cells[ref.id]['source']}: {cells[ref.id]['text']} (row {ref.id})")
        if ref is None:
            out["note"] = ("There is no evaluated row to compare against: the rows are listed without "
                           "differences. Name a reference row, or take ENSDF's values with snapshot_evaluation.")
        rows = []
        for e in c.entries:
            row = cell_view(cells[e.id])
            if e.is_reference:
                row["against_reference"] = "the reference"
            elif ref is not None:
                # A reference of 0 has no difference in percent and no ratio: the difference is in the unit then.
                apart = (f"{e.difference:+.3g} {c.unit}" if e.difference_percent is None
                         else f"{e.difference_percent:+.1f} %")
                if e.pull is not None:
                    row["against_reference"] = (f"{apart}, pull {e.pull:+.2f}: "
                                                f"{compare.VERDICTS[e.verdict].lower()}")
                else:
                    row["against_reference"] = ((apart if e.ratio is None else f"× {e.ratio:.3g}")
                                                + " (no uncertainty: not judged)")
                if e.independent_of_reference is False:
                    row["not_independent_of_reference"] = True
            rows.append(row)
        out["independent"] = rows
        projected = []
        for p in c.projected:
            row = cell_view(cells[p.id])
            row["closure"] = (None if p.closure is None else
                              f"{p.closure:+.2f} standard deviations from the value the simulation was given: "
                              f"{compare.VERDICTS[p.closure_verdict].lower()}")
            row["setup_reaches"] = _percent(p.reach_percent)
            row["known_to"] = _percent(p.known_percent)
            row["best_single_measurement"] = _percent(p.best_measured_percent)
            projected.append(row)
        out["projected"] = projected
        out["rule"] = (f"Pull = difference over the two uncertainties in quadrature. Below {c.agree_below:g}: "
                       f"agree; up to {c.differ_above:g}: tension; above: differ. " + GROUPS)
    out["routes"] = [{"route": r["label"], "experiment": r["kind_label"], "status": r["status"],
                      **({"independent": r["independent"], "projected": r["projected"], "own": r["own"]}
                         if r["rows"] else {}),
                      **({"open": True} if r["open"] else {}), **({"note": r["note"]} if r["note"] else {}),
                      "applies": {True: "yes", False: "no", None: "not tested"}[r["applies"]],
                      **({"reason": r["reason"]} if r["reason"] else {})}
                     for r in v["routes"]]
    out["reference_gives"] = [{"gives": g["to"], **({"value": g["text"], "how": g["how"]} if g["value"] is not None
                                                    else {"missing": g["missing"]})} for g in v["gives"]]
    out["checks"] = [{"check": c["label"], "says": c["text"], "how_to_read_it": c["how"], "rows": c["rows"],
                      **({"reference": c["reference"]} if c["reference"] else {})} for c in v["checks"]]
    out["used_for"] = v["used_for"]
    return out


def row_line(row) -> str:
    return (f"{_nucleus.label(row.observable, row.key)} = "
            f"{_nucleus.value_text(row.value, row.unc('plus'), row.unc('minus'))} {row.unit}"
            + {"up": "↑", "down": "↓"}.get(row.direction, ""))


class ResultTools:
    """The results tools, on one :class:`~physim.mcp.session.Session` (for its experiments folder)."""

    def __init__(self, session: Session):
        self.session = session

    @property
    def home(self) -> Path:
        """Where the nuclei's own rows live. A session started on an experiments folder of its own keeps them in
        ``nuclei`` beside that folder, whatever ``$PHYSIM_NUCLEI`` says, so that a server on a scratch folder never
        writes to the real nuclei home. A session on the default experiments folder uses
        :func:`physim.nuclear.results.nuclei_home`."""
        own = getattr(self.session, "_root", None)
        return own.parent / "nuclei" if own is not None else results.nuclei_home()

    def list_nuclei(self) -> dict:
        """The nuclei physim holds something about, lightest first, each with its experiments (the role is
        "studied", "partner" for the other nucleus of a Coulomb-excitation experiment, or "results"), their runs
        and analysed runs, how many result rows hold now, what is known per observable (evaluated, measured,
        theory, and projected by simulated runs) and what is still open (routes physim can plan that no run of
        one's own has taken, runs not analysed, ENSDF's values not taken). not_in_a_nucleus_yet lists the
        experiments that study no nucleus: elastic scattering, and fusion-evaporation without a studied channel.
        No experiment needs to be open."""
        found = _nucleus.nuclei(self.session.root, self.home)
        apart = [{"experiment": f"{u['name']} [{u['folder']}]", "kind": u["kind_label"], "reason": u["reason"]}
                 for u in found["unassigned"]]
        return plain({"nuclei": [nucleus_view(e) for e in found["nuclei"]], "not_in_a_nucleus_yet": apart,
                      "note": GROUPS})

    def list_results(self, nucleus: str, observable: Optional[str] = None,
                     key: Optional[Union[str, dict]] = None, history: bool = False) -> dict:
        """The results of a nucleus ("62Ni"). Without an observable: its experiments, one line per observable and
        level or transition with what each kind of evidence says, and what is still open. With an observable
        ("B(E2)", "tau", "Qs", "B(M1)", "mu", "E_level") the rows that hold now, each with its id, evidence,
        method, value in the observable's convention (B(E2) up in e2fm4, tau in ps, Qs in efm2, B(M1) up in mu_N2,
        mu in mu_N, E in keV) and source; key narrows to one level or transition, as one says it: "2+1" (ENSDF's
        first 2+), "2+ at 1173 keV", "2+1 > 0+1". B(E2) and tau values made from each other through the level's
        one gamma ray are listed too and say what they were made from. history true lists every row ever written,
        the newest first, with its status (holds, superseded, retracted, retraction): rows are never deleted."""
        name = results.nuclide_name(nucleus)
        if key is not None and observable is None:
            raise ValueError("Name the observable the key belongs to: " + ", ".join(OBSERVABLES) + ".")
        if history:
            key = None if key is None else _nucleus.key_of(name, observable, key, check=False).text
            past = _nucleus.history(results.of_nuclide(name, self.session.root, self.home), observable, key)
            rows = [{"id": h["id"], "written": h["created"], "of": _nucleus.label(h["observable"], h["key"]),
                     "evidence": h["evidence"], "method": h["method_label"], "value": h["text"],
                     "source": h["source"], "status": h["status"],
                     **({"superseded_by": h["superseded_by"]} if h["superseded_by"] else {}),
                     **({"reason": h["reason"]} if h["reason"] else {})} for h in past]
            return plain({"nucleus": name, "history": rows[:ROWS_SHOWN], **_more(len(rows))})
        if observable is None:
            return plain(nucleus_view(_nucleus.view(name, self.session.root, self.home)))
        cells = _nucleus.holding(name, observable, key, self.session.root, self.home)
        return plain({"nucleus": name, "rows": [cell_view(c) for c in cells[:ROWS_SHOWN]], **_more(len(cells)),
                      "note": GROUPS})

    def snapshot_evaluation(self, nucleus: str) -> dict:
        """Take ENSDF's adopted values of a nucleus into its results as evaluated rows: level energies, mean
        lives, quadrupole and magnetic moments, B(E2) and B(M1) values, each with its uncertainty and a reference
        that names ENSDF and the data set's date. Taking it again adds nothing; a value a newer ENSDF copy changed
        gets a new row that supersedes the old one. Needs the local ENSDF copy (scripts/fetch_ensdf.py)."""
        s = _nucleus.snapshot_evaluation(nucleus, self.home)
        return {"nucleus": s["nuclide"], "reference": s["reference"], "dataset_date": s["dataset_date"],
                "added": len(s["added"]), "superseded": len(s["superseded"]), "unchanged": s["unchanged"]}

    def add_data_point(self, nucleus: str, observable: str, key: Union[str, dict], value: float, method: str,
                       unit: Optional[str] = None, direction: Optional[str] = None,
                       uncertainty: Optional[float] = None, stat: Optional[float] = None,
                       sys: Optional[float] = None, plus: Optional[float] = None, minus: Optional[float] = None,
                       evidence: str = "measured", reference: Optional[dict] = None,
                       own_measurement: Optional[dict] = None, detail: str = "",
                       corrects: Optional[str] = None) -> dict:
        """Add a data point to a nucleus's results: a value of an observable with its uncertainty, as published or
        as one's own measurement (never raw spectra). observable is "B(E2)", "tau", "Qs", "B(M1)", "mu" or
        "E_level"; key the level or transition as one says it ("2+1", "2+ at 1173 keV", "2+1 > 0+1", {"energy":
        "1173 keV", "jpi": "2+"}); value a number in unit (B(E2): e2fm4, e2b2 or W.u.; tau: fs, ps, ns, s; Qs: efm2
        or eb; B(M1): mu_N2 or W.u.; mu: mu_N; E_level: keV or MeV; the first of each when left out, ps for tau). A
        B(E2) or B(M1) needs direction "up" (from the lower level, as Coulomb excitation gives it) or "down" (the
        decay; W.u. is always down). The uncertainty is absolute, in the same unit, in one form: uncertainty (one
        total), stat and sys, or plus and minus; a measured point needs one. method is "coulex-yield", "coulex-fit", "rdds", "dsam", "nrf", "inelastic-neutron",
        "electron-scattering", "inelastic-scattering", "transfer", "moment" (a static moment measured: laser
        spectroscopy, transient field, NMR) or "other" ("shell-model" for theory). evidence is "measured"
        (default), "outside program" (GOSIA and the like) or "theory". Say where it is from: reference {"authors",
        "year", "doi", "where" in the paper} or own_measurement {"date": "2027-06-14", "note"}, one of the two.
        detail is a free line beside the method. corrects is the id of the data point this one replaces (give the whole point
        again): the old row is kept and superseded. Everything wrong is said at once."""
        row = _nucleus.add_data_point(
            nucleus, observable, key, value, unit, direction=direction, uncertainty=uncertainty, stat=stat,
            sys=sys, plus=plus, minus=minus, method=method, evidence=evidence, reference=reference,
            own=own_measurement, detail=detail, corrects=corrects, home=self.home)
        out = {"added": row_line(row), "id": row.id, "key": row.key, "evidence": row.evidence,
               "source": _nucleus.source(row)}
        if row.supersedes:
            out["supersedes"] = row.supersedes
        return plain(out)

    def retract_data_point(self, nucleus: str, row: str, reason: str) -> dict:
        """Withdraw a data point of a nucleus (row is its id, or the first characters of it; list_results shows
        the ids). Nothing is deleted: a row that says why supersedes the point, which leaves what holds now and
        stays in the history. Only data points can be withdrawn: the evaluation is replaced by a new snapshot,
        an experiment's result by analysing again."""
        made = _nucleus.retract_data_point(nucleus, row, reason, self.home)
        return {"retracted": made.supersedes, "what": row_line(made), "reason": made.retracted, "id": made.id}

    def compare_results(self, nucleus: str, observable: str, key: Union[str, dict],
                        reference: Optional[str] = None) -> dict:
        """Compare every result for one observable of one level or transition of a nucleus ("62Ni", "B(E2)",
        "2+1 > 0+1"). The reference is the evaluated row, or the row whose id is given. Each independent row
        (measured, evaluated, theory, outside program) gets its difference from the reference, its pull and the
        verdict by a rule fixed beforehand (pull below 2 agree, up to 3 tension, above differ; a row without an
        uncertainty gets its ratio only), and is flagged when it is not independent of the reference. Projected
        rows (simulated runs analysed back) are apart: each is closed against the value its simulation was given,
        beside the precision the setup reaches and how well the value is already known. Then the routes to the
        observable (done: an independent result by that method holds; projected: only simulated runs by it,
        nothing measured; in physim: an experiment of that kind can be planned; later; no experiment: a
        calculation, the evaluation or a literature method), each with whether it applies to this nucleus (yes,
        no, or not tested) and the reason: the nucleus as a target or beam, the level's lifetime against the range
        of a plunger or DSAM, whether a shell-model calculation fits this computer now; "not tested" says what is
        not checked and is no "yes". Then what the reference value gives (B down, Weisskopf units, the matrix
        element, tau or B(E2) through the level's gamma ray); the checks against other observables and the
        calculations (plain ratios of independent rows, such as E(4+)/E(2+), and the rotor reading of B(E2) beside
        Q_s), each with how to read it and what it does not show: pass that on with the number; and what a value of
        the observable is used for, with a reference for each."""
        return plain(comparison_view(_nucleus.observable_view(nucleus, observable, key, reference,
                                                              self.session.root, self.home)))

    #: The tools in the order they are listed to the model, and whether each writes files.
    TOOLS: tuple = (
        ("list_nuclei", False), ("list_results", False), ("snapshot_evaluation", True), ("add_data_point", True),
        ("retract_data_point", True), ("compare_results", False),
    )


def _more(n: int) -> dict:
    if n <= ROWS_SHOWN:
        return {}
    return {"more": f"{n - ROWS_SHOWN} more rows: name an observable and a level or transition"}
