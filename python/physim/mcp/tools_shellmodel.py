"""The shell-model tools: which interactions apply, how big a calculation is, running KSHELL in the background,
reading the results, importing results run elsewhere, and putting chosen values into the open experiment.

physim does not solve the shell model itself: KSHELL does, as an outside program
(:mod:`physim.nuclear.shellmodel_jobs` runs it in a background process and keeps a cache of results). These tools
only turn the model's arguments into a :class:`~physim.nuclear.shellmodel.ShellModelRequest`, and the results into
numbers with units. Each tool is a method of :class:`ShellModelTools`; its docstring is what the model reads.

Two rules of the user's are carried by the tools and said in their docstrings:

- shell-model values never replace ENSDF, derived or user values in a level scheme unless the user explicitly asks
  (``replace=True`` in :meth:`ShellModelTools.apply_shell_model`);
- a calculation may use at most half of the RAM free when it starts; a bigger one is refused unless the user
  explicitly allows it (``allow_over_limit=True`` in :meth:`ShellModelTools.run_shell_model`).
"""

from __future__ import annotations

import copy
import math
import re
from pathlib import Path
from typing import Any, Optional, Union

from ..external import kshell as _kshell
from ..nuclear import shellmodel as _sm
from ..nuclear import shellmodel_import as _import
from ..nuclear import shellmodel_jobs as _jobs
from ..nuclear import shellmodel_spaces as _spaces
from ..nuclear.names import SYMBOLS, parse_nuclide
from .session import Session, ToolProblem
from .shape import plain
from .tools_setup import summary

#: States as the tools take them: J (or J/2 for half-integer spin), parity, and how many of the lowest.
_STATE = re.compile(r"^(\d+)(?:/(\d+))?\s*([+-])\s*(?:[x×:*]\s*)?(\d+)?$")
_JOB = re.compile(r"^[0-9a-f]{16}$")

#: Units of the numbers in the results, by kind.
_UNITS = {"E1": ("e fm", "e2fm2"), "E2": ("e fm2", "e2fm4"), "M1": ("mu_N", "mu_N2")}


# ---------------------------------------------------------------------------------------------------------------
# Reading the arguments


def parse_states(states: Union[str, list, tuple]) -> tuple:
    """((2J, parity, how many), ...) from ``"0+3, 2+3, 4+2"`` (the lowest three 0⁺, three 2⁺, two 4⁺),
    ``"3/2-2"``, ``"0+"`` (one state) or a list of such items. ToolProblem with what is wrong."""
    if isinstance(states, (list, tuple)):
        pieces = [str(s) for s in states]
    elif isinstance(states, str):
        pieces = re.split(r"[,;]", states)
    else:
        raise ToolProblem(f"states is text such as \"0+3, 2+3, 4+2\" (the lowest three 0+, three 2+ and two 4+ "
                          f"states), not {states!r}")
    items = []
    for piece in pieces:  # "0+3 2+3" without commas is read too, and "0+ 3" as one item
        piece = piece.strip()
        items.extend([piece] if _STATE.match(piece) or not piece else piece.split())
    out, seen = [], set()
    for raw in items:
        text = raw.strip()
        if not text:
            continue
        m = _STATE.match(text)
        if not m:
            raise ToolProblem(f"cannot read {text!r} as states: write J, the parity and how many of the lowest, "
                              f"e.g. \"0+3\" (the three lowest 0+ states), \"2+\" (the lowest 2+) or \"3/2-2\"; "
                              f"separate them with commas")
        num, den, sign, count = m.groups()
        if den is not None:
            if den != "2" or int(num) % 2 == 0:
                raise ToolProblem(f"{text!r}: a half-integer J is written with an odd numerator over 2, e.g. 3/2")
            j2 = int(num)
        else:
            j2 = 2 * int(num)
        n = 1 if count is None else int(count)
        if n < 1:
            raise ToolProblem(f"{text!r}: ask for at least one state")
        parity = 1 if sign == "+" else -1
        if (j2, parity) in seen:
            raise ToolProblem(f"{_sm.label(j2, parity)} is asked for twice in {states!r}; give each J and parity "
                              f"once, with how many states")
        seen.add((j2, parity))
        out.append((j2, parity, n))
    if not out:
        raise ToolProblem("no states asked for: give e.g. \"0+3, 2+3, 4+2\"")
    return tuple(out)


def states_text(states: Any) -> str:
    """``"0+3, 2+3, 4+2"`` from ((2J, parity, how many), ...)."""
    return ", ".join(f"{_sm.label(int(j2), int(p))}{int(n)}" for j2, p, n in states)


def _transitions(transitions: Union[str, list, tuple, None]) -> tuple:
    if transitions is None:
        return ()
    items = re.split(r"[,;\s]+", transitions) if isinstance(transitions, str) else [str(t) for t in transitions]
    items = [t.strip().upper() for t in items if t.strip()]
    if items in (["NONE"], []):
        return ()
    bad = [t for t in items if t not in _sm.KINDS]
    if bad:
        raise ToolProblem(f"transitions are any of {', '.join(_sm.KINDS)} (e.g. \"E2, M1\"), or \"none\"; not "
                          f"{', '.join(bad)}")
    return tuple(items)


def _pair(x: Any, name: str) -> Optional[tuple]:
    if x is None:
        return None
    if isinstance(x, str):
        x = [v for v in re.split(r"[,;\s]+", x.strip()) if v]
    try:
        a, b = (float(v) for v in x)
    except (TypeError, ValueError):
        raise ToolProblem(f"{name} is two numbers, proton then neutron, e.g. [1.5, 0.5] or \"1.5, 0.5\"; "
                          f"not {x!r}") from None
    return a, b


def _which(which: Union[str, list, tuple, None]) -> Union[str, list]:
    if which is None:
        return "all"
    if isinstance(which, str):
        items = [w.strip() for w in re.split(r"[,;]", which) if w.strip()]
    else:
        items = [str(w).strip() for w in which if str(w).strip()]
    return "all" if items in ([], ["all"]) else items


def _size(n: float) -> str:
    if n >= 1024**3:
        gb = n / 1024**3
        return f"{gb:.2g} GB" if gb < 10 else f"{gb:.0f} GB"
    return f"{max(n, 0) / 1024**2:.0f} MB"


def _time(seconds: float) -> str:
    if seconds < 90:
        return f"about {max(1, round(seconds))} s"
    if seconds < 5400:
        return f"about {round(seconds / 60)} min"
    return f"about {seconds / 3600:.1f} h"


def _canonical(nuclide: str) -> str:
    try:
        z, a = parse_nuclide(nuclide.strip())
    except (ValueError, KeyError, AttributeError):
        raise ToolProblem(f"cannot read {nuclide!r} as a nuclide; write it as e.g. \"28Si\"") from None
    return f"{a}{SYMBOLS[z]}"


# ---------------------------------------------------------------------------------------------------------------
# What goes back


def estimate_view(est: _spaces.Estimate) -> dict:
    """An :class:`~physim.nuclear.shellmodel_spaces.Estimate` for the model: dimensions, sizes as text, fits."""
    return {
        "fits": est.fits,
        "reason": est.reason,
        "dimensions": {_sm.label(j2, p): d for (j2, p), d in est.dimension_by_state.items()},
        "memory": _size(est.memory_bytes),
        "disk": _size(est.disk_bytes),
        "lanczos_vectors": "in memory" if est.mode_lv_hdd == 0 else "on disk",
        "limit": _size(est.limit_bytes),
        "free_ram_now": _size(est.free_bytes),
        "rough_time": _time(est.seconds_guess) + f" on {est.threads} threads (a rough guess: off by 2x either way)",
        **({"warnings": list(est.warnings)} if est.warnings else {}),
    }


def _state_row(s: _sm.ShellState) -> dict:
    row = {"state": s.label, "ex_kev": round(1000.0 * s.ex_mev, 1), "energy_mev": round(s.energy_mev, 4)}
    if s.t2 is not None:
        row["T"] = s.t2 / 2 if s.t2 % 2 else s.t2 // 2
    if s.q_s is not None:
        row["q_s_efm2"] = s.q_s
    if s.mu is not None:
        row["mu_nm"] = s.mu
    return row


def result_view(result: _sm.ShellModelResult) -> dict:
    """A :class:`~physim.nuclear.shellmodel.ShellModelResult` as numbers with units: what produced it, the states
    (excitation energy, Q_s, μ) and the transitions (B down and up, signed ⟨final‖M‖initial⟩)."""
    labels = {s.key: s.label for s in result.states}
    unsigned = _import.unsigned(result)
    transitions = []
    for t in result.transitions:
        if t.diagonal:
            continue
        me_unit, b_unit = _UNITS.get(t.kind, ("", ""))
        row = {"kind": t.kind, "from": labels.get(t.initial, _sm.label(*t.initial)),
               "to": labels.get(t.final, _sm.label(*t.final)),
               f"b_down_{b_unit}": t.b_down, f"b_up_{b_unit}": t.b_up, f"matrix_element_{me_unit.replace(' ', '')}":
               t.mred, "sign_known": (t.kind, t.initial, t.final) not in unsigned}
        transitions.append({k: v for k, v in row.items() if v is not None})
    out = {
        "nuclide": result.nuclide,
        "interaction": result.interaction,
        "space": result.space,
        "program": result.kshell_version,
        "effective_charges": None if result.eff_charge is None else
        {"e_p": result.eff_charge[0], "e_n": result.eff_charge[1], "unit": "e"},
        "g_factors": None if result.gl is None and result.gs is None else
        {"g_l": result.gl, "g_s": result.gs, "order": "(proton, neutron)"},
        "states": [_state_row(s) for s in sorted(result.states, key=lambda s: (s.ex_mev, s.j2))],
        "transitions": transitions,
        "units": {"ex_kev": "keV above the lowest state", "energy_mev": "MeV (KSHELL's total energy)",
                  "q_s_efm2": "e fm2", "mu_nm": "mu_N", "b_down": "B(from -> to)", "b_up": "B(to -> from)",
                  "matrix_element": "signed <to||M||from>, Bohr-Mottelson convention"},
        "notes": list(result.notes),
    }
    if result.eff_charge is not None and result.nuclide and result.interaction:
        try:
            warn = _spaces.charge_warnings(result.nuclide, result.interaction, result.eff_charge, how="used")
        except ValueError:  # an interaction physim does not know (imported results): nothing to check against
            warn = []
        if warn:
            out["warnings"] = warn
    where = _import.provenance(result)
    if where is not None:
        out["imported_from"] = where
        out["notes"] = out["notes"][1:]  # the first note is the provenance, given above
    if result.dimensions:
        out["largest_m_scheme_dimension"] = max(result.dimensions.values())
    if result.elapsed_s:
        out["elapsed_s"] = round(result.elapsed_s, 1)
    return out


def _status_view(st: dict) -> dict:
    keep = ("job", "state", "step", "step_index", "steps", "message", "submitted", "started", "finished", "refused",
            "log_tail", "build")
    out = {k: st[k] for k in keep if st.get(k) not in (None, "")}
    req = st.get("request")
    if req:
        out["request"] = {"nuclide": req.get("nuclide"), "interaction": req.get("interaction"),
                          "states": states_text(req.get("states", ())),
                          "transitions": ", ".join(req.get("transitions", ())),
                          "effective_charges": req.get("eff_charge")}
    est = st.get("estimate")
    if est:
        out["estimate"] = {"fits": est.get("fits"), "reason": est.get("reason"),
                           "memory": _size(est.get("memory_bytes", 0)), "limit": _size(est.get("limit_bytes", 0))}
    return out


# ---------------------------------------------------------------------------------------------------------------
# The tools


class ShellModelTools:
    """The shell-model tools, on one :class:`~physim.mcp.session.Session` (for the open experiment's nuclides and
    level schemes). Results imported with :meth:`import_shell_model` are kept for the session, by path."""

    def __init__(self, session: Session):
        self.session = session
        self._imported: dict = {}

    # -- helpers ---------------------------------------------------------------------------------------------------

    def _nuclide(self, nuclide: Optional[str], role: Optional[str] = "target") -> str:
        """``nuclide``, or the open experiment's beam or target nuclide (``nuclide`` may also be "beam" or
        "target")."""
        if nuclide is not None and nuclide.strip().lower() in ("beam", "target"):
            role, nuclide = nuclide.strip().lower(), None
        if nuclide:
            return _canonical(nuclide)
        if role not in ("beam", "target"):
            raise ToolProblem("role is \"beam\" or \"target\"")
        if self.session.planner is None:
            raise ToolProblem("Name the nuclide (e.g. \"28Si\"), or open an experiment so its target or beam can be "
                              "used.")
        found = self.session.planner.level_nuclides()[role]
        if found is None:
            raise ToolProblem(f"The open experiment has no {role} nuclide yet; name the nuclide.")
        return _canonical(found)

    def _request(self, nuclide: str, interaction: str, states, transitions="E2, M1", eff_charge=None,
                 threads: Optional[int] = None, allow_over_limit: bool = False) -> _sm.ShellModelRequest:
        nuclide = self._nuclide(nuclide)
        try:
            spec = _spaces.interaction(interaction)
        except ValueError as e:
            raise ToolProblem(f"{e}. shell_model_options lists those that apply to {nuclide}.") from None
        try:
            _spaces.nuclide_valence(nuclide, spec)
        except ValueError as e:
            raise ToolProblem(f"{e}. shell_model_options lists the interactions that apply.") from None
        if threads is not None and (isinstance(threads, bool) or int(threads) < 1):
            raise ToolProblem("threads is 1 or more (leave it out to use every core)")
        try:
            return _sm.ShellModelRequest(nuclide, spec.name, parse_states(states), _transitions(transitions),
                                         _pair(eff_charge, "eff_charge"), threads=threads,
                                         allow_over_limit=bool(allow_over_limit))
        except ValueError as e:
            raise ToolProblem(str(e)) from None

    @staticmethod
    def _estimate(request) -> _spaces.Estimate:
        try:
            return _spaces.estimate(request)
        except OSError as e:
            raise ToolProblem(f"Cannot read the free memory on this computer ({e}), so the size limit cannot be "
                              f"checked.") from None
        except ValueError as e:
            raise ToolProblem(str(e)) from None

    @staticmethod
    def _install() -> _kshell.KshellInstall:
        install = _kshell.find_install()
        if install is None:
            raise ToolProblem("KSHELL was not found: set the environment variable PHYSIM_KSHELL to the folder "
                              "holding kshell.exe and transit.exe (or install physim with its bundled copy). "
                              "Results run elsewhere can still be read with import_shell_model.")
        return install

    def _import(self, path: str, nuclide: Optional[str] = None,
                interaction: Optional[str] = None) -> tuple[str, _sm.ShellModelResult]:
        p = Path(str(path)).expanduser()
        try:
            result = _import.import_path(p, nuclide=nuclide, interaction=interaction)
        except FileNotFoundError as e:
            raise ToolProblem(f"No such file or folder: {e}") from None
        except NotImplementedError as e:
            raise ToolProblem(str(e)) from None
        except ValueError as e:
            raise ToolProblem(f"Cannot import {p}: {e}") from None
        handle = str(p.resolve())
        self._imported[handle] = result
        return handle, result

    def _load(self, job_or_path: str) -> tuple[str, _sm.ShellModelResult]:
        """A result by job ID (from run_shell_model) or by the path of results imported or to import."""
        key = str(job_or_path).strip()
        if _JOB.match(key):
            try:
                result = _jobs.result(key)
            except KeyError:
                raise ToolProblem(f"No shell-model job {key}; list_shell_model_jobs lists them.") from None
            except _jobs.JobFailed as e:
                st = e.status or {}
                tail = st.get("log_tail")
                raise ToolProblem(f"Job {key} has no results: it {st.get('state', 'failed')} ({st.get('message', '')})"
                                  + (f". End of the log:\n{tail}" if tail else "")) from None
            if result is None:
                st = _jobs.status(key)
                raise ToolProblem(f"Job {key} is still {st.get('state')} ({st.get('step') or st.get('message')}); "
                                  f"ask shell_model_status (with wait_seconds) until it is done.")
            return key, result
        p = Path(key).expanduser()
        handle = str(p.resolve()) if p.exists() else key
        if handle in self._imported:
            return handle, self._imported[handle]
        if not p.exists():
            raise ToolProblem(f"{key!r} is neither a shell-model job ID (16 hexadecimal digits, from "
                              f"run_shell_model or list_shell_model_jobs) nor a file or folder of results.")
        return self._import(key)

    def _comparison(self, result: _sm.ShellModelResult, role: str) -> dict:
        """The result against the open experiment's level scheme for ``role``: pairs, and what apply would do."""
        if role not in ("beam", "target"):
            raise ToolProblem("role is \"beam\" or \"target\"")
        p = self.session.require()
        nuclide = p.level_nuclides()[role]
        if nuclide is None:
            raise ToolProblem(f"The open experiment has no {role} nuclide yet.")
        if result.nuclide and parse_nuclide(result.nuclide) != parse_nuclide(nuclide):
            raise ToolProblem(f"The shell-model result is for {result.nuclide}, but the {role} nuclide is "
                              f"{nuclide}.")
        scheme = p.experiment.levels.get(role)
        if scheme is None or not scheme.levels:
            return {"role": role, "nuclide": nuclide, "level_scheme": None,
                    "note": f"The setup has no level scheme for the {role}. apply_shell_model would start one from "
                            f"the shell-model states (lowest as the ground state); lookup_levels reads ENSDF's "
                            f"first, so that shell-model values only fill its gaps."}
        rows, deltas = [], []
        for m in _sm.match(result, scheme):
            row = {"state": m.state.label, "shell_kev": round(m.shell_kev, 1)}
            if m.level is None:
                row["level"] = None
                row["note"] = "no known level of this J and parity to pair with: apply would add it with \"new\""
            else:
                lv = scheme.levels[m.level]
                row.update(level=m.level, known_kev=m.known_kev, known_source=lv.energy.source,
                           delta_kev=round(m.delta_kev, 1))
                if m.level != 0:
                    deltas.append(m.delta_kev)
            rows.append(row)
        preview = copy.deepcopy(scheme)
        try:
            would = _sm.apply(result, preview, "all", replace=False)
        except ValueError as e:
            would = [f"apply would refuse: {e}"]
        out = {"role": role, "nuclide": nuclide, "pairs": rows,
               "pairing": "within each J and parity, the n-th shell-model state with the n-th known level by energy",
               "apply_would": would or ["nothing to change"]}
        if deltas:
            out["rms_deviation_kev"] = round(math.sqrt(sum(d * d for d in deltas) / len(deltas)), 1)
            out["rms_over_states"] = len(deltas)
        return out

    # -- the tools -------------------------------------------------------------------------------------------------

    def shell_model_options(self, nuclide: Optional[str] = None, role: str = "target") -> dict:
        """List the shell-model interactions that apply to a nucleus, for KSHELL calculations run by physim.
        nuclide is e.g. "28Si"; by default the open experiment's target nucleus (role "target"), or its beam
        nucleus with role "beam". For each interaction: the valence space and core, valence protons and neutrons,
        the parities they can make, the M-scheme dimension, the default effective charges and g-factors with where
        they come from, the paper to cite, and a size estimate for a typical request (the lowest three 0+, three
        2+ and two 4+ states for an even nucleus; two each of 1/2, 3/2, 5/2 for an odd one) against the memory
        limit now (half of the RAM free). Laptop-sized spaces only (p, sd, fp, f5pg9, sd-pf), with no truncation."""
        nuclide = self._nuclide(nuclide, role)
        z, a = parse_nuclide(nuclide)
        try:
            free = _spaces.available_memory_bytes()
        except OSError:
            free = None
        found = []
        for opt in _spaces.options(nuclide):
            spec = opt.spec
            parity = opt.parities[0] if opt.parities else 1
            typical = ((0, parity, 3), (4, parity, 3), (8, parity, 2)) if a % 2 == 0 else \
                ((1, parity, 2), (3, parity, 2), (5, parity, 2))
            row = {
                "interaction": spec.name, "label": spec.label, "space": spec.space.name,
                "core": spec.space.core_label, "valence_protons": opt.z_val, "valence_neutrons": opt.n_val,
                "parities": ["+" if p > 0 else "-" for p in opt.parities],
                "m_scheme_dimension": {("+" if p > 0 else "-"): d for p, d in opt.dimension_by_parity.items()},
                "effective_charges": {"e_p": spec.eff_charge[0], "e_n": spec.eff_charge[1], "unit": "e"},
                "g_factors": {"g_l": spec.gl, "g_s": spec.gs, "order": "(proton, neutron)"},
                "defaults_from": spec.defaults_from, "citation": spec.citation,
            }
            if spec.notes:
                row["notes"] = spec.notes
            warn = _spaces.charge_warnings(nuclide, spec)
            if warn:
                row["warnings"] = warn
            if free is not None:
                try:
                    req = _sm.ShellModelRequest(nuclide, spec.name, typical)
                    est = _spaces.estimate(req, free_bytes=free)
                    row["typical_request"] = {"states": states_text(typical), **estimate_view(est)}
                except ValueError as e:
                    row["typical_request"] = {"states": states_text(typical), "fits": None, "reason": str(e)}
            found.append(row)
        out = {"nuclide": nuclide, "interactions": found,
               "memory_limit": None if free is None else f"{_size(free // 2)} (half of the {_size(free)} of RAM free "
                                                        f"now)",
               "kshell_found": _kshell.find_install() is not None}
        if not found:
            out["note"] = (f"No interaction physim knows applies to {nuclide}: it has no valence nucleons in, or does "
                           f"not fit, the p, sd, fp, f5pg9 or sd-pf spaces. Results run elsewhere can be read with "
                           f"import_shell_model.")
        return plain(out)

    def estimate_shell_model(self, nuclide: str, interaction: str, states: str,
                             threads: Optional[int] = None) -> dict:
        """How big a shell-model calculation is, before running it: the M-scheme dimension of each J and parity,
        the memory and disk it needs, whether the Lanczos vectors stay in memory, the limit (half of the RAM free
        now), whether it fits, why, and a rough time. nuclide e.g. "28Si" (or "target"/"beam" for the open
        experiment's); interaction a name from shell_model_options (e.g. "usdb"); states is J, parity and how many
        of the lowest states, separated by commas: "0+3, 2+3, 4+2" means the lowest three 0+, three 2+ and two 4+
        states; "3/2-2" two 3/2- states. threads: how many cores (default all)."""
        req = self._request(nuclide, interaction, states, threads=threads)
        est = self._estimate(req)
        return plain({"nuclide": req.nuclide, "interaction": req.interaction, "states": states_text(req.states),
                      **estimate_view(est)})

    def run_shell_model(self, nuclide: str, interaction: str, states: str, transitions: str = "E2, M1",
                        eff_charge: Optional[Union[str, list]] = None, allow_over_limit: bool = False) -> dict:
        """Start a KSHELL shell-model calculation in the background and return its job ID at once (or, when the
        same calculation is in the cache, its results straight away). nuclide e.g. "28Si" (or "target"/"beam");
        interaction from shell_model_options (e.g. "usdb"); states is J, parity and how many of the lowest, e.g.
        "0+3, 2+3, 4+2" (the lowest three 0+, three 2+ and two 4+); transitions "E2, M1" (any of E1, E2, M1, or
        "none"), computed between all the states asked for; eff_charge (e_p, e_n) e.g. "1.5, 0.5", by default the
        interaction's (see shell_model_options).
        A calculation may use at most half of the RAM free when it starts: one over that limit is refused with the
        estimate. Set allow_over_limit=true ONLY when the user has explicitly said to run it anyway; never on your
        own. Follow the job with shell_model_status, then read it with shell_model_results. Only one KSHELL job
        runs at a time; others wait queued."""
        req = self._request(nuclide, interaction, states, transitions, eff_charge,
                            allow_over_limit=allow_over_limit)
        install = self._install()
        key = _jobs.job_key(req, install.version)
        if (_jobs.job_folder(key) / "results.json").is_file():
            job = _jobs.submit(req, install=install)  # a cache hit: marks the result as used, starts nothing
            out = {"job": job, "state": _jobs.DONE, "cached": True,
                   "message": "The same calculation is in the cache; its results are below."}
            out.update(result_view(_jobs.result(job)))
            return plain(out)
        est = self._estimate(req)
        if not est.fits and not allow_over_limit:
            raise ToolProblem(f"Refused: {est.reason} Nothing was started. Ask for fewer or lower states, or run it "
                              f"only if the user explicitly says to go over the limit (allow_over_limit=true).")
        try:
            job = _jobs.submit(req, install=install)
        except FileNotFoundError as e:
            raise ToolProblem(str(e)) from None
        out = _status_view(_jobs.status(job))
        out["estimate"] = estimate_view(est)
        if not est.fits:
            out["warning"] = "Over the memory limit, run because allow_over_limit was set (the user's explicit choice)."
        out["next"] = (f"Follow it with shell_model_status(\"{job}\", wait_seconds=...), then read it with "
                       f"shell_model_results(\"{job}\").")
        return plain(out)

    def shell_model_status(self, job: str, wait_seconds: float = 0) -> dict:
        """The state of a shell-model job: queued, running (with the step, e.g. "kshell 2+ (3 states)", and step
        k of n), done, failed (with the end of KSHELL's log) or cancelled (refused "memory" when it was too big
        when it started). wait_seconds (up to 60) waits that long for the job to finish before answering, so a
        short job can be followed without many calls."""
        key = str(job).strip()
        try:
            if wait_seconds and wait_seconds > 0:
                try:
                    st = _jobs.wait(key, timeout_s=min(float(wait_seconds), 60.0))
                except TimeoutError:
                    st = _jobs.status(key)
            else:
                st = _jobs.status(key)
        except KeyError as e:
            raise ToolProblem(f"{e.args[0] if e.args else e}; list_shell_model_jobs lists the jobs.") from None
        out = _status_view(st)
        if st.get("state") == _jobs.DONE:
            out["next"] = f"Read the numbers with shell_model_results(\"{key}\")."
        return plain(out)

    def list_shell_model_jobs(self) -> dict:
        """List the shell-model jobs in the cache, most recent first: job ID, state, nuclide, interaction, the
        states asked for, and the step or message. A done job's results can be read again at any time."""
        rows = []
        for st in _jobs.list_jobs():
            v = _status_view(st)
            row = {"job": v.get("job"), "state": v.get("state"), "message": v.get("message")}
            if "request" in v:
                row.update({k: v["request"][k] for k in ("nuclide", "interaction", "states")})
            if v.get("state") in _jobs.ACTIVE and v.get("step"):
                row["step"] = v["step"]
            row["updated"] = st.get("updated") or st.get("finished")
            rows.append(row)
        return plain({"cache": _jobs.cache_root(), "jobs": rows})

    def cancel_shell_model(self, job: str) -> dict:
        """Stop a queued or running shell-model job (and its KSHELL program). A finished job is left as it is.
        Returns its status."""
        try:
            return plain(_status_view(_jobs.cancel(str(job).strip())))
        except KeyError as e:
            raise ToolProblem(f"{e.args[0] if e.args else e}; list_shell_model_jobs lists the jobs.") from None

    def shell_model_results(self, job: str, role: Optional[str] = None) -> dict:
        """The numbers of a finished shell-model job (job ID from run_shell_model), or of results imported with
        import_shell_model (its handle, the path): what produced them (interaction, effective charges, program),
        the states (state labels such as "2+1", the first 2+; excitation energy in keV, Q_s in e fm2, mu in mu_N)
        and the transitions (B(E2) in e2fm4 and B(M1) in mu_N2, down and up, and the signed reduced matrix element
        <to||M||from>). With role "target" or "beam": compared with the open experiment's level scheme for that
        nucleus, state by state (known energy and its source beside the shell-model one, the rms deviation) and
        what apply_shell_model would change and keep (ENSDF, derived and user values are kept)."""
        handle, result = self._load(job)
        out = {"handle": handle}
        out.update(result_view(result))
        if role is not None:
            out["comparison"] = self._comparison(result, role)
        return plain(out)

    def apply_shell_model(self, job_or_path: str, role: str, which: Union[str, list] = "all",
                          replace: bool = False) -> dict:
        """Put shell-model values into the open experiment's level scheme for role "target" or "beam", marked
        with the source "shell model" and a note naming the interaction and effective charges, and save it.
        job_or_path is a job ID from run_shell_model or the handle (path) from import_shell_model. which is "all"
        or a comma-separated choice of quantities ("energy", "E2", "E1", "Q" for quadrupole moments), states
        ("2+1", the first 2+) and "new" (add the states the scheme has no level for); see shell_model_results for
        a preview. Shell-model values only fill gaps: ENSDF, derived and user values are kept.
        Set replace=true ONLY when the user explicitly asks to overwrite those values with shell-model ones;
        never on your own. B(M1) and mu stay in the results (a level scheme keeps electric matrix elements).
        Returns the changes, ok, the summary and the level schemes."""
        if role not in ("beam", "target"):
            raise ToolProblem("role is \"beam\" or \"target\"")
        p = self.session.require()
        handle, result = self._load(job_or_path)
        try:
            changes = p.apply_shell_model(role, result, _which(which), replace=bool(replace))
        except ValueError as e:
            raise ToolProblem(str(e)) from None
        out = {"changes": changes or ["nothing changed"], "source": handle}
        unsigned = _import.unsigned(result)
        if unsigned and (_which(which) == "all" or any(w.upper() in ("E2", "E1") for w in _which(which))):
            out["warning"] = ("Some transitions came as B values only, so the signs of their matrix elements are not "
                              "known; they were stored as positive roots: "
                              + ", ".join(f"{k} {_sm.label(*i)}->{_sm.label(*f)}" for k, i, f in sorted(unsigned)))
        out["ok"] = not p.problems
        out.update(summary(p))
        out["levels"] = plain(p.levels())
        return plain(out)

    def import_shell_model(self, path: str, role: Optional[str] = None, nuclide: Optional[str] = None,
                           interaction: Optional[str] = None) -> dict:
        """Read shell-model results run elsewhere: a folder of KSHELL logs (or one log, or KSHELL's summary_*.txt),
        or a CSV table in physim's format (first line "# format: physim shell-model table 1", then # nuclide:,
        # interaction:, # effective charges:, # reference: lines and kind,state,final,value,unit rows such as
        "B(E2),2+1,0+1,78.2,e2fm4"). nuclide and interaction fill in what the files do not say. Returns the
        results as shell_model_results does, with a handle (the path) that apply_shell_model and
        shell_model_results take; with role "target" or "beam", compared with the open experiment's level
        scheme. Nothing goes into the experiment until apply_shell_model."""
        handle, result = self._import(path, nuclide=None if nuclide is None else self._nuclide(nuclide),
                                      interaction=interaction)
        return self.shell_model_results(handle, role)

    #: The tools in the order they are listed to the model, and whether each writes files.
    TOOLS: tuple = (
        ("shell_model_options", False), ("estimate_shell_model", False), ("run_shell_model", True),
        ("shell_model_status", False), ("list_shell_model_jobs", False), ("cancel_shell_model", True),
        ("shell_model_results", False), ("apply_shell_model", True), ("import_shell_model", False),
    )
