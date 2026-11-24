"""Ledger rows from what the analyses already return: the planner writes a row whenever it produces a result.

Nothing is computed here that the analyses do not compute already. Each **builder** takes what an analysis returns
and makes the :class:`~physim.nuclear.results.Row` objects that say it; each **writer** is the one line the planner
calls after the analysis, which appends the rows to the experiment's ledger (``results.jsonl`` in its folder)::

    from physim.nuclear import result_rows

    row = result_rows.coulex_row(result, experiment, "62ni-on-194pt", "03-beam-24h")   # a Row, not written
    result_rows.record(folder, row)                                                    # appended; see below
    result_rows.copy_shell_model(shell_result)                                         # into the nucleus's ledger

=========================================  ======================================================================
What the planner does                      The rows
=========================================  ======================================================================
:meth:`Planner.analyse` (Coulex)           B(E2)↑ of the ``[reaction]`` transition, "coulex-yield"
                                           (:func:`coulex_row`, as :func:`results.backfill` builds it)
:meth:`Planner.analyse_lifetimes`          τ of the level, "rdds" or "dsam", and the B(E2)↓ derived from it, which
                                           depends on the τ row (:func:`lifetime_rows`)
:meth:`Planner.fit_matrix_elements`        B(E2)↑ of each fitted E2 matrix element and Q_s of each fitted diagonal
                                           one, "coulex-fit" (:func:`fit_rows`)
a shell-model result read or applied       level energies, B(E2)↑ and Q_s as theory, in the nucleus's own ledger
                                           (:func:`shell_model_rows`, :func:`copy_shell_model`)
an old experiment opened                   :func:`results.backfill`, once (:func:`backfill_if_new`)
=========================================  ======================================================================

**Everything a simulated run gives is projected**, with ``truth`` the value the setup held: never measured.

**Rows are appended, never rewritten** (:func:`record`). A new analysis of the same thing (the same observable,
level or transition and method; for a Coulex analysis also the same run) supersedes the row that stood. An
analysis that says exactly what the standing row says adds nothing: pressing "Analyse" twice leaves one row.

**A ledger that cannot be written never fails an analysis.** The writers catch whatever goes wrong (a folder that
cannot be written, a level that cannot be named), say it in the log and return no rows; the analysis's own result
is untouched.

What has no observable in the registry yet is left out: a Coulex transition that is not E2, B(λ) from a lifetime
for a multipole other than E2, fitted matrix elements other than E2, and the shell model's B(M1), B(E1) and μ.
"""

from __future__ import annotations

import copy
import functools
import json
import logging
import math
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Union

from . import observables, results
from .results import PROJECTED, Row

log = logging.getLogger(__name__)

#: How each lifetime analysis (:data:`physim.nuclear.lifetime_plan.ANALYSES`) is called in a row: the ledger's
#: method, and the analysis in words for the line of detail.
LIFETIME_METHODS = {"ddcm": ("rdds", "DDCM"), "ddcm_three_foil": ("rdds", "DDCM, three foils"),
                    "decay_curve": ("rdds", "decay curve"), "lineshape": ("dsam", "line shape")}

#: One analysis writes at a time: a row is compared with the ledger and then appended, as one step.
_LOCK = threading.Lock()
#: The shell-model jobs already copied in this session, by nuclei home and job, so that reading a result again
#: (every time its page is drawn) does not read the nucleus's ledger again.
_COPIED: set = set()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _quietly(what: str):
    """A writer that never fails the analysis it follows: whatever goes wrong is said in the log, and no rows
    come back."""
    def wrap(f):
        @functools.wraps(f)
        def inner(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except Exception as e:  # noqa: BLE001 - the ledger is beside the result, never in its way
                log.warning("no ledger row for %s: %s", what, e)
                return []
        return inner
    return wrap


def _number(x) -> Optional[float]:
    """``x`` as a number, or None if it is not one (None, text, NaN, infinite)."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


# ---------------------------------------------------------------------------------------------------------------
# Appending: supersede the row that stood, or leave it when nothing new is said


def _says(row: Row) -> dict:
    """What a row says, apart from when and by which code it was made and which row it replaced."""
    d = row.to_dict()
    for k in ("id", "created", "code", "supersedes"):
        d.pop(k, None)
    d["origin"] = {k: v for k, v in d.get("origin", {}).items() if k not in ("analysed", "backfilled")}
    return d


def _same_thing(old: Row, new: Row, same_run: bool) -> bool:
    return (old.observable == new.observable and old.key == new.key and old.method == new.method
            and old.origin.get("kind") == "run"
            and (not same_run or old.origin.get("run") == new.origin.get("run")))


def record(folder: Union[str, Path], row: Row, same_run: bool = True) -> Row:
    """Append ``row`` to the ledger of the experiment in ``folder`` as the newest analysis of what it is about,
    and return the row that now stands for it.

    The rows about the same thing are those of the same observable, level or transition and method, from a run
    (with ``same_run``, from the same run: a Coulex analysis is of one run, a lifetime analysis of a series). If
    the one that stands says exactly what ``row`` says, nothing is written and that one is returned; otherwise
    ``row`` is appended, superseding it."""
    with _LOCK:
        earlier = results.current([r for r in results.of_experiment(folder) if _same_thing(r, row, same_run)])
        if earlier:
            if _says(earlier[-1]) == _says(row):
                return earlier[-1]
            row = Row.from_dict({**row.to_dict(), "id": "", "supersedes": earlier[-1].id})
        results.append(folder, row)
        return row


# ---------------------------------------------------------------------------------------------------------------
# Coulomb excitation: B(E2) from the yield


def coulex_row(result, experiment, experiment_name: str, run_name: str) -> Optional[Row]:
    """The row of a Coulomb-excitation analysis (a :class:`physim.nuclear.analysis.Result`, or its plain form as
    ``analysis.json`` keeps it): B(E2)↑ of the transition the ``[reaction]`` of ``experiment`` names, with the
    statistical and systematic uncertainties and the budget made absolute, and the value the setup held as
    ``truth``. Built as :func:`physim.nuclear.results.backfill` builds it from the kept file, so that a row
    written at once and a backfilled one look alike. None for a reaction that is not an E2 Coulomb excitation:
    the registry has no observable for it yet."""
    if experiment.excitation is None or experiment.excitation.multipolarity != "E2":
        return None
    d = result if isinstance(result, dict) else result.to_dict()
    key = results.key_of_reaction(experiment)
    value = float(d["b_e2fm4"])
    budget = {k: abs(value) * float(v) for k, v in (d.get("budget") or {}).items() if _number(v) is not None}
    detail = ", ".join(x for x in (f"{d.get('normalisation')} normalisation" if d.get("normalisation") else "",
                                   f"{d.get('gate')} gate" if d.get("gate") else "") if x)
    return Row(key.nuclide, "B(E2)", key, value, "e2fm4", "coulex-yield", PROJECTED, direction=observables.UP,
               stat=abs(value) * float(d["statistical"]), sys=abs(value) * float(d["systematic"]), budget=budget,
               detail=detail, truth=float(d["truth_e2fm4"]),
               origin=results.from_run(experiment_name, run_name, d.get("settings"),
                                       analysed=d.get("saved") or _now()))


def _kept_analysis(run_folder, result) -> dict:
    """The analysis as the run's folder keeps it (with the time it was saved), if what is kept there is this
    result; else the result's own plain form."""
    from .analysis import RESULT_FILE

    try:
        d = json.loads((Path(run_folder) / RESULT_FILE).read_text(encoding="utf-8"))
        if d["b_e2fm4"] == float(result.b_e2fm4) and d["statistical"] == float(result.statistical):
            return d
    except (OSError, ValueError, TypeError, KeyError):
        pass
    return result.to_dict()


@_quietly("the Coulex analysis")
def write_coulex(planner, result) -> list:
    """After :meth:`Planner.analyse` of a beam run: the analysis's row, appended to the experiment's ledger.
    Returns the row that stands (in a list), or none without a run or for a reaction that has no row."""
    run, folder = planner._beam_run(), planner.experiment_folder
    if run is None or folder is None:
        return []
    row = coulex_row(_kept_analysis(run.folder, result), planner._data_experiment(), Path(folder).name, run.label)
    return [] if row is None else [record(folder, row)]


@_quietly("the analyses kept with the runs")
def backfill_if_new(folder) -> list:
    """When an experiment is opened: rows from the analyses its runs already keep
    (:func:`physim.nuclear.results.backfill`), if the experiment has kept analyses and no ledger yet."""
    from .analysis import RESULT_FILE
    from .runs import ExperimentFolder

    folder = Path(folder)
    if (folder / results.RESULTS_FILE).exists():
        return []
    if not any((run / RESULT_FILE).is_file() for run in ExperimentFolder(folder).run_folders()):
        return []
    return results.backfill(folder)


# ---------------------------------------------------------------------------------------------------------------
# Lifetimes: τ, and the B(E2) derived from it


def _held(analysis: dict) -> dict:
    """What the setup's level scheme holds for the analysed level (the first line of the analysis's comparison):
    {"tau_ps", "b_down"}; empty when the scheme has no lifetime for it."""
    for c in analysis.get("compare") or []:
        if str(c.get("source", "")).startswith("level scheme"):
            return c
    return {}


def lifetime_rows(analysis: dict, experiment_name: str, depends_on: Optional[Row] = None) -> list:
    """The rows of a lifetime analysis (the dict of :meth:`Planner.analyse_lifetimes`): the level's τ, by "rdds"
    or "dsam", and (given the τ row that stands as ``depends_on``) B(E2)↓ of the analysed line, which the
    analysis derives from τ through the γ-ray energy, the branching and the conversion and which therefore
    depends on that row. Both are projected; ``truth`` is what the setup's level scheme holds.

    Called without ``depends_on`` it gives the τ row alone; called again with the τ row that stands in the
    ledger, the B(E2) row alone (:func:`write_lifetimes` does both). No rows for an analysis that is not
    available, without the line it is of (``analysis["transition"]``), or when the level scheme holds no lifetime
    to be the truth. No B row for a multipole other than E2: the registry has no observable for it yet."""
    line, held = analysis.get("transition"), _held(analysis)
    if not analysis.get("available") or not line or analysis.get("method") not in LIFETIME_METHODS:
        return []
    method, name = LIFETIME_METHODS[analysis["method"]]
    nuclide = line["nuclide"]
    upper, lower = (results.level_key(nuclide, x["energy_kev"], x.get("spin"), x.get("parity"))
                    for x in (line["initial"], line["final"]))
    runs = list(analysis.get("runs") or [])
    origin = results.from_run(experiment_name, runs[-1] if runs else "",
                              {"method": analysis["method"], "level": analysis.get("level")}, runs=runs,
                              analysis="lifetime")
    detail = f"{name} of the {line['gamma_kev']:.1f} keV line"
    if depends_on is None:
        if _number(held.get("tau_ps")) is None:
            return []
        budget = {k: abs(float(v)) for k, v in (analysis.get("budget") or {}).items()
                  if k != "statistics" and _number(v) is not None}
        return [Row(nuclide, "tau", upper, float(analysis["tau_ps"]), "ps", method, PROJECTED,
                    stat=abs(float(analysis["stat_ps"])), sys=math.sqrt(sum(v * v for v in budget.values())),
                    budget=budget, detail=detail, truth=float(held["tau_ps"]), origin=origin)]
    b = analysis.get("b_down")
    if analysis.get("b_multipole") != "E2" or not b or _number(held.get("b_down")) is None:
        return []
    return [Row(nuclide, "B(E2)", results.transition_key(upper, lower, "E2"), float(b[0]),
                analysis.get("b_unit") or "e2fm4", method, PROJECTED, direction=observables.DOWN,
                total=abs(float(b[1])), truth=float(held["b_down"]), depends_on=(depends_on.id,), origin=origin,
                detail=f"from τ ({detail}) through the γ-ray energy, the branching and the conversion")]


@_quietly("the lifetime analysis")
def write_lifetimes(planner, analysis: dict) -> list:
    """After :meth:`Planner.analyse_lifetimes`: the τ row and the B(E2) row that depends on it, appended to the
    experiment's ledger. Returns the rows that stand."""
    folder = planner.experiment_folder
    if folder is None or not analysis.get("available"):
        return []
    name = Path(folder).name
    out = [record(folder, row, same_run=False) for row in lifetime_rows(analysis, name)]
    if out:
        out += [record(folder, row, same_run=False) for row in lifetime_rows(analysis, name, depends_on=out[0])]
    return out


# ---------------------------------------------------------------------------------------------------------------
# The fit to all orders: matrix elements


def fit_rows(fit: dict, scheme, experiment_name: str, run_name: str, settings: Optional[dict] = None) -> list:
    """The rows of a fit of matrix elements to a run's yields, to all orders (the dict of
    :meth:`physim.nuclear.multistep.Multistep.fit`; ``scheme`` is the level scheme the run was taken with).

    An E2 element between two levels gives B(E2)↑ = ⟨f‖M‖i⟩² / (2J_i + 1) from the lower level
    (:meth:`LevelScheme.b <physim.nuclear.levels.LevelScheme.b>` with the fitted element in place), a diagonal
    one the level's Q_s (:func:`physim.nuclear.observables.quadrupole_moment`). The fit's uncertainty, from the
    curvature at the minimum, is the statistical one: twice the element's relative uncertainty for a B, the same
    factor as the value for Q_s. ``truth`` is the same quantity from the element the scheme holds. The element as
    fitted is kept in the line of detail and in the origin. Elements of other multipoles are left out."""
    rows = []
    for free, value in fit["values"].items():
        a, b, multipole = free
        held = scheme.matrix_element(a, b, multipole) if multipole == "E2" else None
        if held is None:
            continue
        unc = _number((fit.get("uncertainties") or {}).get(free))
        names = [scheme.levels[n].label or scheme.levels[n].jpi for n in (b, a)]
        text = (f"⟨{names[0]}‖E2‖{names[1]}⟩ = {value:.4g}" + (f" ± {unc:.2g}" if unc is not None else "")
                + f" e fm² fitted to the run's yields, χ² {fit.get('chi2', math.nan):.3g}")
        origin = results.from_run(experiment_name, run_name, settings, chi2=fit.get("chi2"),
                                  degrees_of_freedom=fit.get("degrees_of_freedom"),
                                  matrix_element={"from": a, "to": b, "value": float(value), "uncertainty": unc})
        try:
            if a == b:
                j = scheme.levels[a].spin
                stat = None if unc is None else abs(observables.quadrupole_moment(unc, j))
                rows.append(Row(scheme.nuclide, "Qs", results.key_of_level(scheme, a),
                                observables.quadrupole_moment(float(value), j), "efm2", "coulex-fit", PROJECTED,
                                stat=stat, truth=observables.quadrupole_moment(held.value, j), detail=text,
                                origin=origin))
                continue
            low, up = sorted((a, b), key=lambda n: scheme.levels[n].energy.value)
            fitted = copy.deepcopy(scheme)
            fitted.set_matrix_element(a, b, "E2", float(value), source="user")
            b_up = fitted.b(low, up, "E2")
            stat = None if unc is None or not value else b_up * 2 * unc / abs(value)
            rows.append(Row(scheme.nuclide, "B(E2)", results.key_of_transition(scheme, a, b, "E2"), b_up, "e2fm4",
                            "coulex-fit", PROJECTED, direction=observables.UP, stat=stat,
                            truth=scheme.b(low, up, "E2"), detail=text, origin=origin))
        except (ValueError, TypeError) as e:  # an element whose levels cannot be named: the others still count
            log.warning("no ledger row for the fitted element %s: %s", free, e)
    return rows


@_quietly("the fit to all orders")
def write_fit(planner, fit: dict, role: str, gate: Optional[str] = None) -> list:
    """After :meth:`Planner.fit_matrix_elements` on the current run's yields: its rows, appended to the
    experiment's ledger. Returns the rows that stand."""
    run, folder = planner._beam_run(), planner.experiment_folder
    if run is None or folder is None:
        return []
    settings = {"role": role, "gate": gate, "free": [list(k) for k in fit["values"]]}
    rows = fit_rows(fit, planner._data_experiment().levels[role], Path(folder).name, run.label, settings)
    return [record(folder, row) for row in rows]


# ---------------------------------------------------------------------------------------------------------------
# The shell model: theory rows in the nucleus's own ledger


def shell_model_rows(result, job: Optional[str] = None) -> list:
    """A shell-model calculation (a :class:`physim.nuclear.shellmodel.ShellModelResult`) as theory rows: each
    state's excitation energy (the lowest state, at 0 by definition, aside), each state's Q_s and B(E2)↑ of each
    E2 transition, named through :func:`physim.nuclear.results.key_of_shell_state`. The origin is the calculation
    (interaction, valence space, effective charges; the ``job`` in the cache if known). A state or transition
    that cannot be named is passed over. B(M1), B(E1) and μ have no observable in the registry yet."""
    nuclide = result.nuclide or (result.request.nuclide if result.request is not None else "")
    if not nuclide:
        return []
    charges = dict(zip(("proton", "neutron"), result.eff_charge)) if result.eff_charge else {}
    more = {"program": f"KSHELL {result.kshell_version}"} if result.kshell_version else {}
    if job:
        more["job"] = job
    origin = results.shell_model_job(result.interaction, result.space, charges, **more)
    energies = {s.key: 1e3 * s.ex_mev for s in result.states}
    rows = []

    def key(state) -> results.LevelKey:
        return results.key_of_shell_state(nuclide, state, energies.get(tuple(state)))

    def add(observable, where, value, unit, **more):
        try:
            rows.append(Row(nuclide, observable, where(), value, unit, "shell-model", "theory", origin=origin,
                            **more))
        except (ValueError, TypeError, KeyError) as e:
            log.warning("no ledger row for a shell-model %s of %s: %s", observable, nuclide, e)

    lowest = min(result.states, key=lambda s: s.ex_mev, default=None)
    moments = {t.initial: t.moment for t in result.transitions if t.kind == "E2" and t.diagonal}
    for s in result.states:
        if s is not lowest:
            add("E_level", lambda s=s: key(s.key), 1e3 * s.ex_mev, "keV")
        q = s.q_s if s.q_s is not None else moments.get(s.key)
        if q is not None:
            add("Qs", lambda s=s: key(s.key), q, "efm2")
    for t in result.transitions:
        if t.kind != "E2" or t.diagonal:
            continue
        try:
            named = results.transition_key(key(t.initial), key(t.final), "E2")
        except (ValueError, TypeError, KeyError) as e:
            log.warning("no ledger row for a shell-model B(E2) of %s: %s", nuclide, e)
            continue
        # ↑ is from the key's lower level to its upper one. The calculation's upper state is the key's upper level
        # unless ENSDF has the two levels the other way round, and then its B↓ is the key's B↑.
        b = t.b_up if named.initial.text == key(t.initial).text else t.b_down
        if b is not None:
            add("B(E2)", lambda named=named: named, b, "e2fm4", direction=observables.UP)
    return rows


def _calculation(row: Row) -> tuple:
    """What a shell-model row is of and which calculation it is from, the value aside."""
    o = row.origin
    return (row.observable, row.key, row.unit, row.direction, o.get("interaction"), o.get("space"),
            json.dumps(o.get("effective_charges"), sort_keys=True))


@_quietly("the shell-model results")
def copy_shell_model(result, job: Optional[str] = None, home: Union[str, Path, None] = None) -> list:
    """Copy a shell-model calculation's level energies, B(E2) and Q_s into the nucleus's own ledger
    (:func:`physim.nuclear.results.add_point`), so that cleaning the cache of jobs cannot lose them. Returns the
    rows added: none when they are all there already, however often and by whichever way (a job read, a result
    applied or imported) the same calculation comes by. A value that differs from the row of the same
    calculation (run again with another build, say) supersedes it."""
    rows = shell_model_rows(result, job)
    if not rows:
        return []
    added = []
    with _LOCK:
        have = results.read(results.nucleus_folder(rows[0].nuclide, home) / results.POINTS_FILE)
        standing = {_calculation(r): r for r in results.current([r for r in have if r.method == "shell-model"])}
        for row in rows:
            old = standing.get(_calculation(row))
            if old is not None:
                if old.value == row.value:
                    continue
                row = Row.from_dict({**row.to_dict(), "id": "", "supersedes": old.id})
            if results.add_point(row, home):
                standing[_calculation(row)] = row
                added.append(row)
    return added


@_quietly("the shell-model results")
def read_shell_model(result, job: str) -> list:
    """After a finished job's result is read from the cache (:func:`physim.nuclear.shellmodel_jobs.result`): copy
    it into the nucleus's ledger, once per session and nuclei home."""
    mark = (str(results.nuclei_home()), job)
    if mark in _COPIED:
        return []
    _COPIED.add(mark)
    return copy_shell_model(result, job)


__all__ = ["LIFETIME_METHODS", "backfill_if_new", "copy_shell_model", "coulex_row", "fit_rows", "lifetime_rows",
           "read_shell_model", "record", "shell_model_rows", "write_coulex", "write_fit", "write_lifetimes"]
