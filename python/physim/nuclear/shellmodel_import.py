"""Shell-model results run elsewhere, read into a :class:`~physim.nuclear.shellmodel.ShellModelResult`.

A collaborator's KSHELL run on a cluster, KSHELL's own summary of one, or numbers from a paper or an email are as
usable as a calculation run here: :func:`import_path` reads them into the same result, which goes into an
experiment's level scheme only through :func:`physim.nuclear.shellmodel.apply` (so it never replaces ENSDF,
derived or user values unless the user says so)::

    from physim.nuclear import shellmodel, shellmodel_import

    result = shellmodel_import.import_path("from_alice/Si28_usdb")         # a folder of KSHELL logs
    result = shellmodel_import.import_path("summary_Si28_usdb.txt")        # KSHELL's collect_logs.py summary
    result = shellmodel_import.import_path("si28_paper.csv")               # physim's table format (below)
    shellmodel_import.provenance(result)       # {'path': ..., 'format': 'physim shell-model table', 'reference': ...}
    shellmodel_import.to_table(result, "si28.csv")                         # and back to a table

What :func:`import_path` recognises:

- **KSHELL logs**: a folder, one log file or a list of them (``kshell.exe`` and ``transit.exe`` logs, told apart by
  their contents, read by :func:`physim.nuclear.shellmodel.read_logs`). The nuclide comes from the logs' core and
  valence nucleons, the interaction from the interaction file they name; ``nuclide=`` and ``interaction=`` check
  or override them. Logs of runs that are not J-projected (an M = 0 run such as KSHELL's ⁴⁸Cr test) are read as
  they are. Logs from MPI builds have not been seen yet; they are read if they are written the same way.
- **KSHELL's summary** (``summary_*.txt`` from ``collect_logs.py``) when the logs themselves are not there: level
  energies to 1 keV and B(E2)/B(M1) to KSHELL's printed rounding. The summary has no signs, no Q_s or μ, and no
  effective charges; the result says so in its notes, and its matrix elements are the positive roots of the
  B values (see :func:`unsigned`).
- **physim's table format** (``# format: physim shell-model table 1``, below) for numbers typed from a paper.
- **NuShellX output** is not read yet (:class:`NotImplementedError`): no real NuShellX output files were available
  to test a reader against, and physim does not guess a file format. Type the numbers into a table meanwhile.

Every imported result's first note says where it came from (path, format and the reference given in the file);
:func:`provenance` reads it back.

The table format
----------------

A CSV file (as the other reference tables in physim) in "long" form: one value per row, every value with its unit
on the same row, so nothing depends on column positions or on units implied by a heading. The header lines start
with ``#`` and give what produced the numbers::

    # format: physim shell-model table 1
    # nuclide: 28Si
    # interaction: USDB
    # effective charges: 1.31, 0.46
    # reference: A. Author et al., Journal 1, 2 (2024), Table 3
    kind,state,final,value,unit
    Ex,0+1,,0.0,MeV
    Ex,2+1,,1.779,MeV
    Q,2+1,,+18.5,e fm2
    B(E2),2+1,0+1,78.3,e2fm4
    M(E2),4+1,2+1,-31.47,e fm2

Header entries: ``format`` (required, the first line), ``nuclide`` and ``interaction`` (required here or as
arguments; the interaction may be any model name, as the paper gives it), and optionally ``space``,
``effective charges`` (e_p, e_n: needed to make sense of E2 values), ``g_l`` and ``g_s`` (proton, neutron),
``program`` (e.g. ANTOINE, NuShellX, KSHELL), ``reference``, ``notes`` (may be repeated) and ``dimensions``
(``2M parity dimension`` triplets separated by ``;``, e.g. ``0 + 93710; 4 + 81122``). An entry continues on the
next line that starts with ``#`` and two or more spaces; ``##`` starts a comment; an unknown entry is an error
(so a misspelt one is not silently dropped).
``state`` (and ``final``) are written J, parity, index: ``2+1``, ``5/2-2``; ``2+`` alone means the only 2⁺ state
in the file. Kinds and the units they take (spaces, ``^`` and ``_`` in units do not matter):

=============== ======================================================= ==============================
kind            value                                                   units
=============== ======================================================= ==============================
``Ex``          excitation energy of ``state``                          ``MeV``, ``keV``
``E``           absolute (binding) energy, as KSHELL prints it          ``MeV``, ``keV``
``T``           isospin (``0``, ``1``, ``1/2``)                         none
``Q``           spectroscopic quadrupole moment                         ``e fm2``, ``e b``
``mu``          magnetic dipole moment                                  ``mu_N`` (``nm``)
``B(E2)``       B(E2; state → final)                                    ``e2fm4``, ``e2b2``, ``W.u.``
``B(M1)``       B(M1; state → final)                                    ``mu_N2``, ``W.u.``
``B(E1)``       B(E1; state → final)                                    ``e2fm2``, ``e2b``, ``W.u.``
``M(E2)``       signed reduced matrix element ⟨final‖E2‖state⟩          ``e fm2``, ``e b``
``M(M1)``       signed ⟨final‖M1‖state⟩                                 ``mu_N``
``M(E1)``       signed ⟨final‖E1‖state⟩                                 ``e fm``
=============== ======================================================= ==============================

Conventions (KSHELL's and physim's): Bohr–Mottelson reduced matrix elements, B(i → f) = ⟨f‖M‖i⟩² / (2J_i + 1),
so B(E2; 0⁺ → 2⁺) is five times B(E2; 2⁺ → 0⁺); a transition may be written either way round and is stored from
the higher state down. Weisskopf units are (1/4π)(3/(λ+3))²(1.2 A^⅓)^{2λ} e² fm^{2λ} and 1.790 μ_N² for M1. A
diagonal ``M(E2)`` (state = final) is read as Q_s. Every state a transition names needs an ``Ex`` or ``E`` row.
Where a B value is given without the matching ``M`` row the sign is not known: the result keeps the positive root
and lists such transitions in a note (:func:`unsigned`). A pair of states may have one ``M`` row and one ``B`` row
in each direction (B↓ and B↑); they must agree to their printed rounding, and B values are kept as given.

:func:`to_table` writes a result in this format (E, Ex, T, Q, μ, M and B rows); reading it back gives the same
states, transitions, charges and notes (not the proton and neutron parts of Q_s and μ, nor the time taken: the
result's ``to_dict()`` keeps everything).
"""

from __future__ import annotations

import csv
import io
import math
import re
from dataclasses import replace as _replace
from pathlib import Path
from typing import Iterable, Optional, Union

from .levels import quadrupole_factor, weisskopf_unit
from .names import SYMBOLS, parse_nuclide
from .shellmodel import KINDS, ShellModelResult, ShellState, ShellTransition, label, parse_label, read_logs
from .shellmodel_spaces import INTERACTIONS

PathLike = Union[str, Path]

#: The format line a table starts with (any line of the header will do).
TABLE_FORMAT = "physim shell-model table"
TABLE_VERSION = 1
#: Names of the formats, as the provenance note gives them.
FORMATS = ("KSHELL logs", "KSHELL summary", "physim shell-model table")
#: The start of the note listing transitions whose sign is not known.
UNSIGNED_PREFIX = "signs not known (given as B only; ⟨f‖M‖i⟩ is kept as the positive root) for: "
#: B_W(M1) in μ_N²: (10/π)(3/4)².
WEISSKOPF_M1 = 10.0 / math.pi * (3.0 / 4.0) ** 2


# ---------------------------------------------------------------------------------------------------------------
# What a file is, and where a result came from


def import_path(path: Union[PathLike, Iterable[PathLike]], nuclide: Optional[str] = None,
                interaction: Optional[str] = None) -> ShellModelResult:
    """A shell-model result from ``path``: a folder of KSHELL logs, one KSHELL log or a list of them, KSHELL's
    ``collect_logs.py`` summary, or a table in physim's format (see the module's documentation). ``nuclide`` and
    ``interaction`` fill in what the files do not say, and are checked against what they do say (a different
    nuclide is an error; a different interaction name replaces the file's, with a note). The result's first note
    names the path, the format and the reference (:func:`provenance`).

    Raises ``FileNotFoundError`` for a path that is not there, ``ValueError`` for files that cannot be read (with
    what is wrong), and ``NotImplementedError`` for NuShellX output."""
    if isinstance(path, (str, Path)):
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"{p} does not exist")
        if p.is_dir():
            return _import_files(sorted(q for q in p.iterdir() if q.is_file()), p, nuclide, interaction)
        return _import_files([p], p, nuclide, interaction)
    paths = [Path(q) for q in path]
    if not paths:
        raise ValueError("no files given")
    for q in paths:
        if not q.is_file():
            raise FileNotFoundError(f"{q} is not a file")
    where = paths[0].parent if len(paths) > 1 else paths[0]
    return _import_files(paths, where, nuclide, interaction)


def _kind_of(text: str) -> Optional[str]:
    """``"kshell"``, ``"transit"``, ``"summary"``, ``"table"`` or ``None`` for a file's text."""
    head = text[:4000]
    if re.search(r"^\s*#\s*format\s*:\s*" + re.escape(TABLE_FORMAT), head, re.M | re.I):
        return "table"
    if "FN_LOAD_WAVE_L" in text:
        return "transit"
    if "<H>:" in text or re.search(r"^\s*FN_PTN\s*=", text, re.M):
        return "kshell"
    if re.search(r"^\s*Energy levels\s*$", head, re.M) and re.search(r"\bN_Jp\b.*\bEx\(MeV\)", head):
        return "summary"
    return None


def _read_text(path: Path) -> str:
    return path.read_bytes().decode("utf-8", errors="replace")


def _import_files(paths: list, where: Path, nuclide: Optional[str], interaction: Optional[str]) -> ShellModelResult:
    single = len(paths) == 1 and where == paths[0]
    nushellx = [p for p in paths if p.suffix.lower() == ".lpt"]
    kinds: dict = {"kshell": [], "transit": [], "summary": [], "table": []}
    unknown = []
    for p in paths:
        if p.suffix.lower() in (".lpt",):
            continue
        if not single and p.suffix.lower() not in (".txt", ".log", ".csv"):
            continue  # a run folder also holds inputs, partitions, wavefunctions, scripts
        text = _read_text(p)
        kind = _kind_of(text)
        if kind is None:
            unknown.append(p)
        else:
            kinds[kind].append((p, text))
    if kinds["kshell"] or kinds["transit"]:
        if kinds["table"]:
            raise ValueError(f"{where} holds both KSHELL logs and a physim table; import them one at a time")
        version = "unknown"
        if where.is_dir() and (where / "VERSION").is_file():
            version = (where / "VERSION").read_text(encoding="utf-8").strip() or "unknown"
        result = read_logs([t for _, t in kinds["kshell"]], [t for _, t in kinds["transit"]], None, version)
        counts = f"{len(kinds['kshell'])} kshell.exe and {len(kinds['transit'])} transit.exe logs"
        _settle(result, nuclide, interaction, f"the KSHELL logs in {where}")
        extra = [f"read {counts}"]
        for p, text in kinds["summary"]:
            extra.append(_compare_summary(result, p, text))
        _mark(result, where, "KSHELL logs", None, extra)
        return result
    if kinds["summary"]:
        if len(kinds["summary"]) > 1:
            names = ", ".join(p.name for p, _ in kinds["summary"])
            raise ValueError(f"{where} holds several KSHELL summaries ({names}); import one at a time")
        p, text = kinds["summary"][0]
        return parse_summary(text, nuclide, interaction, source=p)
    if kinds["table"]:
        if len(kinds["table"]) > 1:
            names = ", ".join(p.name for p, _ in kinds["table"])
            raise ValueError(f"{where} holds several shell-model tables ({names}); import one at a time")
        p, text = kinds["table"][0]
        return parse_table(text, nuclide, interaction, source=p)
    if nushellx:
        read_nushellx(nushellx[0])
    if single:
        raise ValueError(f"cannot tell what {where} is: it is not a KSHELL log (kshell.exe or transit.exe), "
                         f"KSHELL's collect_logs.py summary, or a physim shell-model table (a CSV file with the line "
                         f"'# format: {TABLE_FORMAT} {TABLE_VERSION}')")
    raise ValueError(f"{where} holds no KSHELL logs, KSHELL summary or physim shell-model table"
                     + (f" (files not recognised: {', '.join(p.name for p in unknown[:5])})" if unknown else ""))


def read_nushellx(path: PathLike) -> ShellModelResult:
    """Not done yet: NuShellX output (``.lpt`` level files and its transition files) needs real NuShellX output to
    write and test a reader against, and none was available. Raises :class:`NotImplementedError`."""
    raise NotImplementedError(
        f"{path}: NuShellX output is not read yet. A reader needs real NuShellX output files to be tested against "
        f"(the .lpt level list and the transition output for one calculation, with the inputs that made them); "
        f"physim does not guess the format. Meanwhile, type the numbers into a physim shell-model table (see "
        f"tests/reference/nuclear/shellmodel/import/table_template.csv) and import that.")


def _interaction_name(text: str) -> str:
    """Known interactions by physim's name (``"USDB"`` → ``"usdb"``); anything else as written."""
    t = text.strip()
    key = t.lower()[:-4] if t.lower().endswith(".snt") else t.lower()
    return key if key in INTERACTIONS else t


def _canonical(nuclide: str) -> str:
    z, a = parse_nuclide(nuclide.strip())
    return f"{a}{SYMBOLS[z]}"


def _settle(result: ShellModelResult, nuclide: Optional[str], interaction: Optional[str], what: str) -> None:
    """Check or fill in the nuclide and interaction; check 2J against A."""
    if nuclide:
        given = _canonical(nuclide)
        if result.nuclide and parse_nuclide(result.nuclide) != parse_nuclide(given):
            raise ValueError(f"{what} are for {result.nuclide}, not {given}")
        result.nuclide = given
    if not result.nuclide:
        raise ValueError(f"cannot tell from {what} which nuclide they are for; pass nuclide=, e.g. '28Si'")
    if interaction and interaction.strip():
        name = _interaction_name(interaction)
        if result.interaction and result.interaction.lower() != name.lower():
            result.notes.append(f"the files name the interaction '{result.interaction}'; called '{name}' as asked")
        result.interaction = name
    if not result.interaction:
        raise ValueError(f"cannot tell from {what} which interaction was used; pass interaction=, e.g. 'usdb'")
    a = parse_nuclide(result.nuclide)[1]
    odd = sorted({s.label for s in result.states if (s.j2 + a) % 2})
    if odd:
        raise ValueError(f"{what}: states {', '.join(odd)} cannot belong to {result.nuclide} (A = {a} needs "
                         f"{'half-integer' if a % 2 else 'integer'} J); is the nuclide right?")


def _mark(result: ShellModelResult, where: Optional[Path], fmt: str, reference: Optional[str],
          extra: Iterable[str] = ()) -> None:
    text = f"imported from {'text given directly' if where is None else Path(where).resolve()} ({fmt})"
    if reference:
        text += f"; reference: {reference}"
    result.notes[:0] = [text] + [e for e in extra if e]


_PROVENANCE = re.compile(r"^imported from (?P<path>.+) \((?P<format>" + "|".join(re.escape(f) for f in FORMATS)
                         + r")\)(?:; reference: (?P<reference>.*))?$", re.S)


def provenance(result: ShellModelResult) -> Optional[dict]:
    """Where an imported result came from: ``{"path", "format", "reference"}`` (reference ``None`` when the file
    gave none), from its first note; ``None`` for a result that was not imported."""
    for note in result.notes[:1]:
        m = _PROVENANCE.match(note)
        if m:
            return {"path": m.group("path"), "format": m.group("format"), "reference": m.group("reference")}
    return None


def unsigned(result: ShellModelResult) -> set:
    """The transitions whose sign is not known (read from a B value only), as (kind, initial, final)."""
    out = set()
    for note in result.notes:
        if note.startswith(UNSIGNED_PREFIX):
            for item in note[len(UNSIGNED_PREFIX):].split(", "):
                kind, _, pair = item.partition(" ")
                a, _, b = pair.partition("→")
                ja, jb = parse_label(a), parse_label(b)
                out.add((kind, (ja[0], ja[1], ja[2] or 1), (jb[0], jb[1], jb[2] or 1)))
    return out


def _unsigned_note(transitions: list, items: list) -> Optional[str]:
    """The note listing ``items`` ((kind, initial, final)), in the order of ``transitions``."""
    chosen = set(items)
    ordered = [(t.kind, t.initial, t.final) for t in transitions if (t.kind, t.initial, t.final) in chosen]
    if not ordered:
        return None
    return UNSIGNED_PREFIX + ", ".join(f"{k} {label(*i)}→{label(*f)}" for k, i, f in ordered)


# ---------------------------------------------------------------------------------------------------------------
# KSHELL's summary (collect_logs.py)


_SUMMARY_LEVEL = re.compile(r"^\s*(\d+)\s+(\d+(?:/2)?)\s+([+-])\s+(\d+)\s+(-?\d+(?:/2)?)\s+(-?\d+\.\d+)\s+"
                            r"(-?\d+\.\d+)(?:\s+(\S+))?\s*$")
_SUMMARY_SECTION = re.compile(r"^\s*B\(([EM]\d)\)\s+\(.*mass\s*=\s*(\d+)")
_SUMMARY_STATE = r"(\d+(?:/2)?)([+-])\(\s*(\d+)\)"
_SUMMARY_B = re.compile(r"^\s*" + _SUMMARY_STATE + r"\s+(-?\d+\.\d+)\s+" + _SUMMARY_STATE + r"\s+(-?\d+\.\d+)\s+"
                        r"(-?\d+\.\d+)\s+(\S+)\(\s*\S+\)\s+(\S+)\(\s*\S+\)\s*$")
_LOG_NAME = re.compile(r"(?:log|summary)_([A-Z][a-z]?)(\d+)_(.+?)(?:_tr)?(?:_[jm]\d+[pn])*\.txt$")


def _j2(text: str) -> int:
    return int(text[:-2]) if text.endswith("/2") else 2 * int(text)


def _decimals(text: str) -> int:
    t = text.strip().lower().split("e")[0]
    return len(t.split(".")[1]) if "." in t else 0


def _half_unit(text: str) -> float:
    """Half the last printed digit of a number written as ``text`` (its rounding)."""
    t = text.strip().lower()
    exp = int(t.split("e")[1]) if "e" in t else 0
    return 0.5 * 10.0 ** (exp - _decimals(t))


def _space_text(name: str) -> str:
    inter = INTERACTIONS.get(name)
    if inter is None:
        return ""
    sp = inter.space
    return " ".join(o.label[2:] for o in sp.proton_orbits) + f" outside {sp.core_label}"


def parse_summary(text: str, nuclide: Optional[str] = None, interaction: Optional[str] = None,
                  source: Optional[PathLike] = None) -> ShellModelResult:
    """A result from KSHELL's ``collect_logs.py`` summary: the energy-level table and the B(E2)/B(M1) tables.
    The nuclide and interaction come from the log-file names it lists (``log_Si28_usdb_j0p.txt``) or the
    summary's own name unless given. Matrix elements are the positive roots of the B values (the summary prints
    no signs), from whichever B column is printed more precisely; Q_s, μ and the effective charges are not in the
    summary."""
    lines = text.splitlines()
    name = Path(source).name if source is not None else "the summary"
    states: dict = {}
    order: list = []
    found_nuclide = found_inter = None
    mass = None
    raw: list = []  # (kind, initial key, final key, ex_i, ex_f, b_lr text, b_rl text, line number)
    section = None
    for n, line in enumerate(lines, start=1):
        if m := _SUMMARY_SECTION.match(line):
            section, mass = m.group(1), int(m.group(2))
            continue
        if m := _SUMMARY_LEVEL.match(line):
            key = (_j2(m.group(2)), 1 if m.group(3) == "+" else -1, int(m.group(4)))
            if key in states:
                raise ValueError(f"{name}, line {n}: {label(*key)} is listed twice")
            states[key] = (float(m.group(6)), float(m.group(7)), _j2(m.group(5)) if m.group(5) != "-1" else None)
            order.append(key)
            if m.group(8) and (lm := _LOG_NAME.search(m.group(8))):
                found_nuclide = found_nuclide or f"{lm.group(2)}{lm.group(1)}"
                found_inter = found_inter or lm.group(3)
            section = None
            continue
        if section and (m := _SUMMARY_B.match(line)):
            ki = (_j2(m.group(1)), 1 if m.group(2) == "+" else -1, int(m.group(3)))
            kf = (_j2(m.group(5)), 1 if m.group(6) == "+" else -1, int(m.group(7)))
            raw.append((section, ki, kf, float(m.group(4)), float(m.group(8)), m.group(10), m.group(11), n))
    if not states:
        raise ValueError(f"{name}: no energy-level table found (the 'N J prty N_Jp T E(MeV) Ex(MeV) log-file' "
                         f"lines of KSHELL's collect_logs.py summary)")
    if source is not None and not found_nuclide and (lm := _LOG_NAME.search(Path(source).name)):
        found_nuclide, found_inter = f"{lm.group(2)}{lm.group(1)}", lm.group(3)
    nuc = None
    if found_nuclide:
        try:
            nuc = _canonical(found_nuclide)
        except ValueError:
            nuc = None
    if nuc and mass is not None and parse_nuclide(nuc)[1] != mass:
        raise ValueError(f"{name}: the log files name {nuc} but the B tables say mass = {mass}")

    shell_states = [ShellState(k[0], k[1], k[2], e, ex, t2) for k, (e, ex, t2) in
                    sorted(states.items(), key=lambda kv: (kv[1][0], order.index(kv[0])))]
    by_key = {s.key: s for s in shell_states}
    transitions, signless = [], []
    for kind, ki, kf, ex_i, ex_f, b_lr, b_rl, n in raw:
        for k in (ki, kf):
            if k not in by_key:
                raise ValueError(f"{name}, line {n}: {label(*k)} is not in the energy-level table")
        if kind not in KINDS:
            raise ValueError(f"{name}, line {n}: B({kind}) is not a multipole physim reads ({', '.join(KINDS)})")
        si, sf = by_key[ki], by_key[kf]
        b_i, b_f = float(b_lr), float(b_rl)  # B(i → f) and B(f → i) as printed
        # ⟨f‖M‖i⟩² from the more precise of the two columns
        err_i, err_f = _half_unit(b_lr) * (si.j2 + 1), _half_unit(b_rl) * (sf.j2 + 1)
        m2 = b_i * (si.j2 + 1) if err_i <= err_f else b_f * (sf.j2 + 1)
        if si.energy_mev < sf.energy_mev:  # listed upwards: turn it round
            si, sf, b_i, b_f = sf, si, b_f, b_i
        t = ShellTransition(kind, si.key, sf.key, math.sqrt(max(m2, 0.0)), b_i, b_f, None, si.energy_mev,
                            sf.energy_mev)
        transitions.append(t)
        signless.append((kind, t.initial, t.final))
    transitions.sort(key=lambda t: (KINDS.index(t.kind), -t.e_initial_mev, -t.e_final_mev))
    inter = _interaction_name(found_inter) if found_inter else ""
    result = ShellModelResult(None, _space_text(inter.lower()), inter, "unknown", None, None, None, shell_states,
                              transitions, {}, 0.0, [], nuc or "")
    _settle(result, nuclide, interaction, name)
    if mass is not None and parse_nuclide(result.nuclide)[1] != mass:
        raise ValueError(f"{name}: the B tables say mass = {mass}, not {result.nuclide}")
    if not result.space:
        result.space = _space_text(result.interaction.lower())
    result.notes.append("KSHELL's summary gives energies to 1 keV and B values to its printed rounding; it has no "
                        "signs, no Q_s or μ, and does not say which effective charges or g-factors were used")
    if note := _unsigned_note(result.transitions, signless):
        result.notes.append(note)
    _mark(result, None if source is None else Path(source), "KSHELL summary", None)
    return result


def _compare_summary(result: ShellModelResult, path: Path, text: str) -> str:
    """A note on a summary found beside the logs: not used, and whether it agrees with them."""
    try:
        summary = parse_summary(text, result.nuclide or None, None, source=path)
    except ValueError as e:
        return f"{path.name} (a KSHELL summary) was not used, and could not be read: {e}"
    off = [s.label for s in summary.states
           if (r := result.state(*s.key)) is None or abs(r.energy_mev - s.energy_mev) > 0.0015]
    if off:
        return (f"{path.name} (a KSHELL summary) was not used; it disagrees with the logs for "
                f"{', '.join(off[:8])} (from another run?)")
    return f"{path.name} (a KSHELL summary) was not used: the logs hold more; it agrees with them"


# ---------------------------------------------------------------------------------------------------------------
# physim's table format


_HEADER_KEYS = ("format", "nuclide", "interaction", "space", "effective charges", "g_l", "g_s", "program",
                "reference", "notes", "dimensions")
_KEY_ALIASES = {"model": "interaction", "eff_charge": "effective charges", "effective charge": "effective charges",
                "gl": "g_l", "gs": "g_s", "note": "notes", "source": "reference"}
_COLUMNS = ("kind", "state", "final", "value", "unit")

_UNITS = {
    "energy": {"mev": 1.0, "kev": 1e-3},
    "Q": {"efm2": 1.0, "eb": 100.0},
    "mu": {"mun": 1.0, "nm": 1.0},
    ("B", "E2"): {"e2fm4": 1.0, "e2b2": 1e4, "w.u.": "wu"},
    ("B", "E1"): {"e2fm2": 1.0, "e2b": 100.0, "w.u.": "wu"},
    ("B", "M1"): {"mun2": 1.0, "nm2": 1.0, "w.u.": "wu"},
    ("M", "E2"): {"efm2": 1.0, "eb": 100.0},
    ("M", "E1"): {"efm": 1.0},
    ("M", "M1"): {"mun": 1.0, "nm": 1.0},
}
_UNIT_NAMES = {"energy": "MeV or keV", "Q": "e fm2 or e b", "mu": "mu_N",
               ("B", "E2"): "e2fm4, e2b2 or W.u.", ("B", "E1"): "e2fm2, e2b or W.u.", ("B", "M1"): "mu_N2 or W.u.",
               ("M", "E2"): "e fm2 or e b", ("M", "E1"): "e fm", ("M", "M1"): "mu_N"}


def _unit_key(text: str) -> str:
    t = text.strip().lower().replace("μ", "mu").replace("µ", "mu").replace("²", "2").replace("⁴", "4")
    t = re.sub(r"[\s^_*]", "", t)
    return "w.u." if t in ("wu", "w.u.", "w.u") else t


def _kind(text: str) -> Optional[tuple]:
    t = re.sub(r"\s", "", text).lower()
    if t in ("ex", "e", "t", "q", "mu"):
        return {"ex": ("Ex",), "e": ("E",), "t": ("T",), "q": ("Q",), "mu": ("mu",)}[t]
    m = re.fullmatch(r"([bm])\(([em]\d)\)", t)
    if m and m.group(2).upper() in KINDS:
        return m.group(1).upper(), m.group(2).upper()
    return None


def _state_ref(text: str) -> tuple:
    t = re.sub(r"\s", "", text).replace("−", "-")
    t = re.sub(r"^(.+?[+-])[(_]?(\d+)\)?$", r"\1\2", t)
    return parse_label(t)


def _pair_numbers(text: str, what: str) -> tuple:
    parts = [p for p in re.split(r"[,;\s]+", re.sub(r"[a-z_]+\s*=", " ", text.strip(), flags=re.I)) if p]
    try:
        values = [float(p) for p in parts]
    except ValueError:
        values = []
    if len(values) != 2:
        raise ValueError(f"{what} must be two numbers (proton, neutron), got '{text.strip()}'")
    return values[0], values[1]


def _float(text: str, what: str) -> float:
    t = text.strip().replace("−", "-")
    try:
        if re.fullmatch(r"[+-]?\d+/2", t):
            return int(t.split("/")[0]) / 2
        return float(t)
    except ValueError:
        raise ValueError(f"{what}: '{text.strip()}' is not a number") from None


def _header(lines: list, name: str) -> tuple:
    """The header entries (key → list of texts) and the remaining (line number, text) lines."""
    entries: dict = {}
    body: list = []
    last = None
    for n, line in enumerate(lines, start=1):
        s = line.strip()
        if not s or s.startswith("##"):
            continue
        if not s.startswith("#"):
            body.append((n, line))
            last = None
            continue
        if re.match(r"^#\s{2,}\S", line.lstrip()) and last is not None:
            entries[last][-1] += " " + s[1:].strip()
            continue
        m = re.match(r"^#\s*([^:]+?)\s*:\s*(.*)$", s)
        if not m:
            last = None
            continue  # a plain comment
        key = m.group(1).strip().lower()
        key = _KEY_ALIASES.get(key, key)
        if key not in _HEADER_KEYS:
            raise ValueError(f"{name}, line {n}: '{m.group(1).strip()}' is not a header entry; use "
                             f"{', '.join(_HEADER_KEYS)} (start a comment with ##)")
        if key in entries and key != "notes":
            raise ValueError(f"{name}, line {n}: '{key}' is given twice")
        entries.setdefault(key, []).append(m.group(2).strip())
        last = key
    return entries, body


def parse_table(text: str, nuclide: Optional[str] = None, interaction: Optional[str] = None,
                source: Optional[PathLike] = None) -> ShellModelResult:
    """A result from a table in physim's format (see the module's documentation). ``nuclide`` and
    ``interaction`` are used when the header does not give them, and checked against it when it does. Raises
    ``ValueError`` naming the line for anything malformed."""
    name = Path(source).name if source is not None else "the table"
    entries, body = _header(text.splitlines(), name)
    fmt = entries.get("format", [""])[0]
    m = re.fullmatch(re.escape(TABLE_FORMAT) + r"\s+(\d+)", fmt, re.I)
    if not m:
        raise ValueError(f"{name}: the first header line must be '# format: {TABLE_FORMAT} {TABLE_VERSION}'")
    if int(m.group(1)) > TABLE_VERSION:
        raise ValueError(f"{name}: written in version {m.group(1)} of the table format; this physim reads up to "
                         f"{TABLE_VERSION}")
    header_nuc = entries.get("nuclide", [""])[0]
    if header_nuc and nuclide and parse_nuclide(header_nuc.strip()) != parse_nuclide(nuclide.strip()):
        raise ValueError(f"{name} is for {_canonical(header_nuc)}, not {_canonical(nuclide)}")
    nuc_text = header_nuc or nuclide
    if not nuc_text:
        raise ValueError(f"{name}: no nuclide; add '# nuclide: 28Si' or pass nuclide=")
    try:
        nuc = _canonical(nuc_text)
    except ValueError as e:
        raise ValueError(f"{name}: {e}") from None
    a_mass = parse_nuclide(nuc)[1]
    inter = entries.get("interaction", [""])[0] or (interaction or "")
    if not inter.strip():
        raise ValueError(f"{name}: no interaction or model; add '# interaction: USDB' or pass interaction=")
    eff = _pair_numbers(entries["effective charges"][0], f"{name}: effective charges") \
        if "effective charges" in entries else None
    gl = _pair_numbers(entries["g_l"][0], f"{name}: g_l") if "g_l" in entries else None
    gs = _pair_numbers(entries["g_s"][0], f"{name}: g_s") if "g_s" in entries else None
    dims = {}
    for part in (entries.get("dimensions", [""])[0]).split(";"):
        if not part.strip():
            continue
        dm = re.fullmatch(r"\s*(-?\d+)\s+([+-])1?\s+(\d+)\s*", part)
        if not dm:
            raise ValueError(f"{name}: dimensions are '2M parity dimension' separated by ';', got '{part.strip()}'")
        dims[(int(dm.group(1)), 1 if dm.group(2) == "+" else -1)] = int(dm.group(3))

    # -- the column row, then the values
    if not body:
        raise ValueError(f"{name}: no values; after the header comes the line '{','.join(_COLUMNS)}'")
    n0, head_line = body[0]
    columns = [c.strip().lower() for c in next(csv.reader([head_line]))]
    missing = [c for c in _COLUMNS if c not in columns]
    extra = [c for c in columns if c not in _COLUMNS + ("comment",)]
    if missing or extra:
        raise ValueError(f"{name}, line {n0}: the columns must be {','.join(_COLUMNS)} (and optionally comment); "
                         + (f"missing {', '.join(missing)}" if missing else f"'{', '.join(extra)}' not known"))
    col = {c: columns.index(c) for c in columns}

    level_values: dict = {}   # (j2, parity, index or None) -> {"Ex": (value, line), ...}
    trans_rows: list = []     # the B and M rows
    for n, line in body[1:]:
        cells = next(csv.reader([line]))
        if len(cells) < len(columns):
            cells += [""] * (len(columns) - len(cells))
        if len(cells) > len(columns) and any(c.strip() for c in cells[len(columns):]):
            raise ValueError(f"{name}, line {n}: more cells than columns")
        get = {c: cells[i].strip() for c, i in col.items()}
        where = f"{name}, line {n}"
        kind = _kind(get["kind"])
        if kind is None:
            raise ValueError(f"{where}: '{get['kind']}' is not a kind; use Ex, E, T, Q, mu, B(E2), B(M1), B(E1), "
                             f"M(E2), M(M1) or M(E1)")
        if not get["state"]:
            raise ValueError(f"{where}: no state")
        try:
            ref = _state_ref(get["state"])
            fin = _state_ref(get["final"]) if get["final"] else None
        except ValueError as e:
            raise ValueError(f"{where}: {e}") from None
        for r in (ref, fin):
            if r is not None and (r[0] + a_mass) % 2:
                raise ValueError(f"{where}: {label(r[0], r[1])} cannot be a state of {nuc} (A = {a_mass} needs "
                                 f"{'half-integer' if a_mass % 2 else 'integer'} J)")
            if r is not None and r[2] is not None and r[2] < 1:
                raise ValueError(f"{where}: state indices start at 1")
        if not get["value"]:
            raise ValueError(f"{where}: no value")
        value = _float(get["value"], where)
        unit = _unit_key(get["unit"])
        if len(kind) == 1:
            if fin is not None:
                raise ValueError(f"{where}: {kind[0]} belongs to one state; leave 'final' empty")
            if kind[0] == "T":
                if unit:
                    raise ValueError(f"{where}: isospin has no unit")
                if value < 0 or abs(2 * value - round(2 * value)) > 1e-9:
                    raise ValueError(f"{where}: isospin is 0, 1/2, 1, …, not {get['value']}")
                stored = round(2 * value)
            else:
                group = "energy" if kind[0] in ("Ex", "E") else kind[0]
                factor = _UNITS[group].get(unit)
                if factor is None:
                    raise ValueError(f"{where}: '{get['unit']}' is not a unit for {kind[0]}; use {_UNIT_NAMES[group]}")
                stored = value * factor
            slot = level_values.setdefault(ref, {})
            if kind[0] in slot:
                raise ValueError(f"{where}: {kind[0]} of {get['state']} is given twice (line {slot[kind[0]][1]})")
            slot[kind[0]] = (stored, n, get["value"], factor if kind[0] not in ("T",) else 1.0)
            continue
        typ, mult = kind
        if fin is None:
            raise ValueError(f"{where}: {get['kind']} needs a final state")
        factor = _UNITS[(typ, mult)].get(unit)
        if factor is None:
            raise ValueError(f"{where}: '{get['unit']}' is not a unit for {typ}({mult}); use "
                             f"{_UNIT_NAMES[(typ, mult)]}")
        if factor == "wu":
            lam = int(mult[1])
            factor = WEISSKOPF_M1 if mult == "M1" else weisskopf_unit(lam, a_mass)
        level_values.setdefault(ref, {})
        level_values.setdefault(fin, {})
        trans_rows.append({"typ": typ, "mult": mult, "value": value * factor, "from": ref, "to": fin, "line": n,
                           "half": _half_unit(get["value"]) * abs(factor)})

    # -- states: indices, energies
    refs = list(level_values)
    for r in refs:
        if r[2] is None and any(o[:2] == r[:2] and o[2] not in (None, 1) for o in refs):
            raise ValueError(f"{name}: {label(r[0], r[1])} without an index means the only {label(r[0], r[1])} state, "
                             f"but there are more; number them ({label(r[0], r[1], 1)}, {label(r[0], r[1], 2)}, …)")

    def key_of(r: tuple) -> tuple:
        return (r[0], r[1], 1 if r[2] is None else r[2])

    values: dict = {}
    for r, v in level_values.items():
        slot = values.setdefault(key_of(r), {})
        for k, item in v.items():
            if k in slot:
                raise ValueError(f"{name}, line {item[1]}: {k} of {label(*key_of(r))} is given twice (line "
                                 f"{slot[k][1]})")
            slot[k] = item
    if not values:
        raise ValueError(f"{name}: no states")
    no_energy = [label(*k) for k, v in values.items() if "Ex" not in v and "E" not in v]
    if no_energy:
        raise ValueError(f"{name}: no energy for {', '.join(sorted(no_energy))}; give an Ex (or E) row for every "
                         f"state a value names")
    both = [(v["E"][0] - v["Ex"][0], v["E"][3] * _half_unit(v["E"][2]) + v["Ex"][3] * _half_unit(v["Ex"][2]))
            for v in values.values() if "E" in v and "Ex" in v]
    notes: list = []
    offset = None
    if both:
        offset = both[0][0]
        for d, h in both[1:]:
            if abs(d - offset) > h + both[0][1] + 1e-9:
                raise ValueError(f"{name}: E - Ex is not the same for every state ({offset:g} and {d:g} MeV); the "
                                 f"absolute and excitation energies disagree")
    elif any("E" in v for v in values.values()) and any("E" not in v for v in values.values()):
        raise ValueError(f"{name}: some states have only E (absolute energy) and others only Ex; give the same kind "
                         f"for all, or both for at least one state")
    elif all("E" in v for v in values.values()):
        offset = min(v["E"][0] for v in values.values())
    else:
        notes.append("the table gives excitation energies only; energy_mev holds the excitation energy")
    states = []
    for k, v in values.items():
        ex = v["Ex"][0] if "Ex" in v else v["E"][0] - offset
        e = v["E"][0] if "E" in v else (ex + offset if offset is not None else ex)
        states.append(ShellState(k[0], k[1], k[2], e, ex, v.get("T", (None,))[0],
                                 v.get("Q", (None,))[0], v.get("mu", (None,))[0]))
    states.sort(key=lambda s: (s.ex_mev, s.j2, s.index))
    by_key = {s.key: s for s in states}
    q_half = {k: v["Q"][3] * _half_unit(v["Q"][2]) for k, v in values.items() if "Q" in v}

    # -- transitions
    groups: dict = {}  # (multipole, the two states) -> {"M": row, "B": [rows, one per direction]}
    for row in trans_rows:
        row["from"], row["to"] = key_of(row["from"]), key_of(row["to"])
        slot = groups.setdefault((row["mult"], frozenset((row["from"], row["to"]))), {"M": None, "B": []})
        twin = slot["M"] if row["typ"] == "M" else next((r for r in slot["B"] if r["from"] == row["from"]), None)
        if twin is not None:
            raise ValueError(f"{name}, line {row['line']}: {row['typ']}({row['mult']}) between {label(*row['from'])} "
                             f"and {label(*row['to'])} is given twice (line {twin['line']})")
        if row["typ"] == "M":
            slot["M"] = row
        else:
            slot["B"].append(row)
    transitions, signless = [], []
    for (mult, _), slot in groups.items():
        first = slot["M"] or slot["B"][0]
        a, b = first["from"], first["to"]
        sa, sb = by_key[a], by_key[b]
        where = f"{name}, line {first['line']}"
        lam = int(mult[1])
        if a == b:
            if mult != "E2" or slot["B"]:
                raise ValueError(f"{where}: a state with itself is read only as a diagonal M(E2) (its Q_s); give "
                                 f"Q or mu rows for moments")
            fq = quadrupole_factor(sa.spin)
            if fq == 0:
                raise ValueError(f"{where}: {sa.label} has no quadrupole moment (J < 1)")
            q = first["value"] * fq
            if sa.q_s is not None and abs(sa.q_s - q) > (first["half"] * fq + q_half.get(a, 0.0)) * 1.01 + 1e-9:
                raise ValueError(f"{where}: the diagonal M(E2) gives Q_s = {q:+.5g} e fm², the Q row "
                                 f"{sa.q_s:+.5g}")
            if sa.q_s is None:
                by_key[a] = _replace(sa, q_s=q)
            continue
        # selection rules
        parity_ok = (sa.parity == sb.parity) == (mult[0] == "E" and lam % 2 == 0 or mult[0] == "M" and lam % 2 == 1)
        if not parity_ok or abs(sa.j2 - sb.j2) > 2 * lam or sa.j2 + sb.j2 < 2 * lam:
            raise ValueError(f"{where}: {mult} between {sa.label} and {sb.label} breaks the selection rules")
        # ⟨f‖M‖i⟩² from each row, with its rounding: M exactly (to 2|M|h), B(a → b) × (2J_a + 1)
        estimates = []  # (M², its rounding, row)
        m_row = slot["M"]
        if m_row is not None:
            estimates.append((m_row["value"] ** 2, 2 * abs(m_row["value"]) * m_row["half"], m_row))
        for b_row in slot["B"]:
            if b_row["value"] < 0:
                raise ValueError(f"{name}, line {b_row['line']}: a B value cannot be negative")
            ja = by_key[b_row["from"]].j2 + 1
            estimates.append((b_row["value"] * ja, b_row["half"] * ja, b_row))
        for m2_x, h_x, row_x in estimates[1:]:
            m2_0, h_0, row_0 = estimates[0]
            if abs(m2_x - m2_0) > (h_x + h_0) * 1.01 + 1e-9 * max(m2_x, m2_0):
                def show(row: dict) -> str:
                    value = f"{row['value']:+.5g}" if row["typ"] == "M" else f"{row['value']:.5g}"
                    return f"{row['typ']}({mult}; {label(*row['from'])} → {label(*row['to'])}) = {value}"
                raise ValueError(f"{name}, lines {row_0['line']} and {row_x['line']}: {show(row_0)} and "
                                 f"{show(row_x)} disagree (they give M² = {m2_0:.5g} and {m2_x:.5g})")
        m2 = estimates[0][0] if m_row is not None else min(estimates, key=lambda e: e[1])[0]
        mred_ba = m_row["value"] if m_row is not None else math.sqrt(m2)  # ⟨b‖M‖a⟩
        # orient from the higher state down: initial i, final f, mred = ⟨f‖M‖i⟩
        if sa.ex_mev >= sb.ex_mev:
            si, sf, mred = sa, sb, mred_ba
        else:
            si, sf = sb, sa
            mred = mred_ba * (-1.0 if ((sa.j2 - sb.j2) // 2) % 2 else 1.0)
        b_down, b_up = m2 / (si.j2 + 1), m2 / (sf.j2 + 1)
        for b_row in slot["B"]:  # keep B values as given
            if b_row["from"] == si.key:
                b_down = b_row["value"]
            else:
                b_up = b_row["value"]
        transitions.append(ShellTransition(mult, si.key, sf.key, mred, b_down, b_up, None, si.energy_mev,
                                           sf.energy_mev))
        if m_row is None:
            signless.append((mult, si.key, sf.key))
    states = sorted(by_key.values(), key=lambda s: (s.ex_mev, s.j2, s.index))
    transitions.sort(key=lambda t: (KINDS.index(t.kind), -t.e_initial_mev, -t.e_final_mev, t.initial, t.final))

    program = entries.get("program", ["unknown"])[0] or "unknown"
    all_notes = list(entries.get("notes", []))
    all_notes += [n for n in notes if n not in all_notes]
    result = ShellModelResult(None, entries.get("space", [""])[0], _interaction_name(inter), program, eff, gl, gs,
                              states, transitions, dims, 0.0, all_notes, nuc)
    if interaction and entries.get("interaction") and \
            _interaction_name(interaction).lower() != result.interaction.lower():
        result.notes.append(f"the table names the interaction '{result.interaction}'; called "
                            f"'{_interaction_name(interaction)}' as asked")
        result.interaction = _interaction_name(interaction)
    if note := _unsigned_note(result.transitions, signless):
        result.notes.append(note)
    reference = entries.get("reference", [None])[0] or None
    _mark(result, None if source is None else Path(source), "physim shell-model table", reference)
    return result


def _num(x: float, digits: int = 10) -> str:
    text = f"{x:.{digits}g}"
    return "0" if text in ("-0", "0") else text


def to_table(result: ShellModelResult, path: Optional[PathLike] = None, reference: Optional[str] = None) -> str:
    """``result`` in physim's table format (see the module's documentation), written to ``path`` when given. The
    reference is ``reference`` or, for an imported result, the one it was imported with. Transitions with a known
    sign get an ``M`` row and a ``B`` row (B↓), others B↓ and B↑ rows (no ``M``, so the sign stays unknown, and
    both B values come back as they were); moments of states with themselves are written as ``Q`` and ``mu``
    rows."""
    if not result.nuclide:
        raise ValueError("the result does not say which nuclide it is for; set result.nuclide first")
    prov = provenance(result)
    reference = reference or (prov or {}).get("reference")
    out = [f"# format: {TABLE_FORMAT} {TABLE_VERSION}", f"# nuclide: {result.nuclide}",
           f"# interaction: {result.interaction or 'not named'}"]
    if result.space:
        out.append(f"# space: {result.space}")
    if result.eff_charge is not None:
        out.append(f"# effective charges: {_num(result.eff_charge[0])}, {_num(result.eff_charge[1])}")
    if result.gl is not None:
        out.append(f"# g_l: {_num(result.gl[0])}, {_num(result.gl[1])}")
    if result.gs is not None:
        out.append(f"# g_s: {_num(result.gs[0])}, {_num(result.gs[1])}")
    if result.kshell_version and result.kshell_version != "unknown":
        out.append(f"# program: {result.kshell_version}")
    if reference:
        out.append(f"# reference: {reference}")
    if result.dimensions:
        out.append("# dimensions: " + "; ".join(f"{m2} {'+' if p > 0 else '-'} {d}"
                                               for (m2, p), d in sorted(result.dimensions.items())))
    for note in result.notes:
        if not note.startswith(UNSIGNED_PREFIX):
            out.append(f"# notes: {' '.join(note.split())}")
    rows = [list(_COLUMNS)]
    absolute = any(abs(s.energy_mev - s.ex_mev) > 1e-12 for s in result.states)
    for s in sorted(result.states, key=lambda s: (s.ex_mev, s.j2, s.index)):
        if absolute:
            rows.append(["E", s.label, "", _num(s.energy_mev), "MeV"])
        rows.append(["Ex", s.label, "", _num(s.ex_mev), "MeV"])
        if s.t2 is not None:
            rows.append(["T", s.label, "", str(s.t2 // 2) if s.t2 % 2 == 0 else f"{s.t2}/2", ""])
        if s.q_s is not None:
            rows.append(["Q", s.label, "", _num(s.q_s), "e fm2"])
        if s.mu is not None:
            rows.append(["mu", s.label, "", _num(s.mu), "mu_N"])
    signless = unsigned(result)
    units_m = {"E2": "e fm2", "M1": "mu_N", "E1": "e fm"}
    units_b = {"E2": "e2fm4", "M1": "mu_N2", "E1": "e2fm2"}
    for t in result.transitions:
        if t.diagonal:
            continue
        i, f = label(*t.initial), label(*t.final)
        known_sign = (t.kind, t.initial, t.final) not in signless
        if known_sign:
            rows.append([f"M({t.kind})", i, f, _num(t.mred), units_m[t.kind]])
        b = t.b_down if t.b_down is not None else t.mred ** 2 / (t.initial[0] + 1)
        rows.append([f"B({t.kind})", i, f, _num(b), units_b[t.kind]])
        if not known_sign:  # B↑ too: with no M row, the more precise of the two gives |M|
            b = t.b_up if t.b_up is not None else t.mred ** 2 / (t.final[0] + 1)
            rows.append([f"B({t.kind})", f, i, _num(b), units_b[t.kind]])
    buf = io.StringIO()
    csv.writer(buf, lineterminator="\n").writerows(rows)
    text = "\n".join(out) + "\n" + buf.getvalue()
    if path is not None:
        Path(path).write_text(text, encoding="utf-8", newline="\n")
    return text


def read_table(path: PathLike, nuclide: Optional[str] = None, interaction: Optional[str] = None) -> ShellModelResult:
    """:func:`parse_table` on a file."""
    p = Path(path)
    return parse_table(_read_text(p), nuclide, interaction, source=p)


__all__ = ["FORMATS", "TABLE_FORMAT", "TABLE_VERSION", "UNSIGNED_PREFIX", "WEISSKOPF_M1", "import_path",
           "parse_summary", "parse_table", "provenance", "read_nushellx", "read_table", "to_table", "unsigned"]
