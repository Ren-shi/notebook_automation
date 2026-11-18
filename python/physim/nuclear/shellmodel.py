"""Shell-model results: KSHELL's logs read into numbers, and chosen values put into a level scheme.

physim does not solve the shell model itself. KSHELL (N. Shimizu et al., Comput. Phys. Commun. 244, 372 (2019))
does, as an outside program: ``kshell.exe`` finds the states of each J± and ``transit.exe`` the E2 and M1
matrix elements between them. This module reads their logs into a :class:`ShellModelResult` and puts the values a
user picks into an experiment's :class:`~physim.nuclear.levels.LevelScheme` with the source ``"shell model"``::

    from physim.nuclear import shellmodel
    from physim.nuclear.levels import LevelScheme

    result = shellmodel.read_folder("runs/Si28_usdb")       # the kshell and transit logs of one calculation
    result.state(4, +1, 1)                                  # the 2⁺₁: energy, Ex, Q_s, μ
    result.transition("E2", (4, 1, 1), (0, 1, 1)).b_down    # B(E2; 2⁺₁ → 0⁺₁) in e² fm⁴
    scheme = LevelScheme.from_ensdf("28Si")
    shellmodel.match(result, scheme)                        # which shell-model state is which known level
    shellmodel.apply(result, scheme, ["Q"])                 # fill in the quadrupole moments ENSDF lacks

Conventions (KSHELL's): a state is labelled by **2J**, **parity ±1** and its **index** (1 for the lowest state of
its 2J±); energies are in **MeV** here and become keV only in a level scheme. Electric matrix elements are in
e fm^λ and B(Eλ) in e² fm^{2λ}; M1 matrix elements in μ_N and B(M1) in μ_N².

How the numbers are read (worked out from the logs and KSHELL's source, which is only read, never copied):

- ``kshell.exe`` prints, for every state, ``<H>`` (the energy), ``<JJ>`` and J, ``<TT>`` and T, the proton and
  neutron quadrupole moments ``<Qp>``, ``<Qn>`` and ``<eQ>`` = e_p·Q_p + e_n·Q_n (the spectroscopic moment Q_s,
  already turned to M = J), and the orbital and spin parts ``<Lp> <Ln> <sp> <sn>`` with
  ``<gm>`` = g_l,p·L_p + g_l,n·L_n + g_s,p·s_p + g_s,n·s_n, which is the magnetic moment μ. The M1 parts are not
  printed for an M = 0 run (they vanish there), so μ is then unknown.
- ``transit.exe`` prints one table per multipole. ``Mred.`` is the signed reduced matrix element ⟨L‖M‖R⟩ between
  the left (first, "f") and right (second, "i") state, in the Bohr–Mottelson convention, so that
  B(R → L) = Mred² / (2J_R + 1) is the ``B(EM)<-`` column and B(L → R) = Mred² / (2J_L + 1) the ``B(EM)->``
  column. ``Mom.`` on a diagonal line is the moment of that state (Q_s or μ).

Here a transition always runs from the higher state (``initial``) to the lower one (``final``), and ``mred`` is
⟨final‖M‖initial⟩, turned round with ⟨a‖M‖b⟩ = (−1)^(J_a − J_b) ⟨b‖M‖a⟩ when KSHELL printed it the other way. The
signs are real relative signs (unlike ENSDF's, which are unknown), so they are kept.

Precedence in a level scheme: shell-model values fill gaps. They never replace an ENSDF, derived or user value
unless the caller asks with ``replace=True`` (the user's explicit choice).
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field, replace as _replace
from pathlib import Path
from typing import Any, Iterable, Optional, Union

from .levels import ELECTRIC, Level, LevelScheme, MatrixElement, Value, quadrupole_factor
from .names import SYMBOLS, parse_nuclide

#: The source of every value this module writes into a level scheme.
SOURCE = "shell model"
#: Sources a shell-model value never replaces unless asked to.
PROTECTED = ("ensdf", "derived", "user")
#: Multipoles transit.exe computes.
KINDS = ("E1", "E2", "M1")
#: How close (MeV) a transit log's energy (3 decimals) must be to a kshell log's (5 decimals) to be the same state.
ENERGY_TOLERANCE_MEV = 0.0015

State = tuple  # (2J, parity, index)


# ---------------------------------------------------------------------------------------------------------------
# The request, the states, the transitions, the result


def _pair(x: Any, name: str) -> Optional[tuple]:
    if x is None:
        return None
    try:
        a, b = x
        return float(a) + 0.0, float(b) + 0.0
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be two numbers (proton, neutron), got {x!r}") from None


def label(j2: int, parity: int, index: Optional[int] = None) -> str:
    """``"2+1"`` for (4, +1, 1), ``"3/2-2"`` for (3, −1, 2); without an index ``"2+"``."""
    j = str(j2 // 2) if j2 % 2 == 0 else f"{j2}/2"
    return j + ("+" if parity > 0 else "-") + ("" if index is None else str(index))


_LABEL = re.compile(r"^(\d+)(/2)?([+-])(\d+)?$")


def parse_label(text: str) -> tuple:
    """(2J, parity, index or None) from ``"2+1"``, ``"3/2-"``, ``"0+3"``."""
    m = _LABEL.match(text.strip())
    if not m:
        raise ValueError(f"cannot read '{text}' as a state; write J, parity and index, e.g. '2+1' or '3/2-2'")
    j2 = int(m.group(1)) if m.group(2) else 2 * int(m.group(1))
    if m.group(2) and j2 % 2 == 0:
        raise ValueError(f"'{text}': a half-integer J has an odd numerator")
    return j2, 1 if m.group(3) == "+" else -1, None if m.group(4) is None else int(m.group(4))


@dataclass(frozen=True)
class ShellModelRequest:
    """What to compute: the nuclide, the interaction, how many states of each 2J± (``states`` holds
    (2J, parity, how many)), the multipoles, and the effective charges and g-factors (``None``: the interaction's
    defaults). ``threads`` and ``allow_over_limit`` change how it runs, not the result."""

    nuclide: str
    interaction: str
    states: tuple = ()
    transitions: tuple = ("E2", "M1")
    eff_charge: Optional[tuple] = None
    gl: Optional[tuple] = None
    gs: Optional[tuple] = None
    threads: Optional[int] = None
    allow_over_limit: bool = False

    def __post_init__(self):
        z, a = parse_nuclide(self.nuclide)
        object.__setattr__(self, "nuclide", f"{a}{SYMBOLS[z]}")
        if not isinstance(self.interaction, str) or not self.interaction.strip():
            raise ValueError("name the interaction, e.g. 'usdb'")
        object.__setattr__(self, "interaction", self.interaction.strip().lower())
        states, seen = [], set()
        for s in self.states:
            try:
                j2, parity, n = (int(x) for x in s)
            except (TypeError, ValueError):
                raise ValueError(f"each state is (2J, parity, how many), got {s!r}") from None
            if j2 < 0 or parity not in (1, -1) or n < 1:
                raise ValueError(f"state {s!r}: 2J must be 0 or more, parity +1 or -1, and at least one state")
            if (j2 + a) % 2:
                raise ValueError(f"state {label(j2, parity)}: 2J must be {'odd' if a % 2 else 'even'} for "
                                 f"{self.nuclide}")
            if (j2, parity) in seen:
                raise ValueError(f"{label(j2, parity)} is asked for twice")
            seen.add((j2, parity))
            states.append((j2, parity, n))
        object.__setattr__(self, "states", tuple(sorted(states, key=lambda s: (-s[1], s[0]))))
        kinds = []
        for k in ([self.transitions] if isinstance(self.transitions, str) else self.transitions):
            k = str(k).upper()
            if k not in KINDS:
                raise ValueError(f"transitions are {', '.join(KINDS)}, not '{k}'")
            if k not in kinds:
                kinds.append(k)
        object.__setattr__(self, "transitions", tuple(kinds))
        for name in ("eff_charge", "gl", "gs"):
            object.__setattr__(self, name, _pair(getattr(self, name), name))
        if self.threads is not None and (isinstance(self.threads, bool) or int(self.threads) < 1):
            raise ValueError("threads must be 1 or more")

    def key(self, build: str = "") -> str:
        """A stable hash of everything that changes the result: nuclide, interaction, states, multipoles,
        effective charges and g-factors (not threads or the size override). ``build`` names the KSHELL build, so
        a job cache can tell results of different builds apart. ``None`` charges stay ``None`` here; resolve
        them to the interaction's defaults first if those may change."""
        payload = {"nuclide": self.nuclide, "interaction": self.interaction,
                   "states": [list(s) for s in self.states], "transitions": sorted(self.transitions),
                   "eff_charge": self.eff_charge, "gl": self.gl, "gs": self.gs, "build": build}
        text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]

    def to_dict(self) -> dict:
        return {"nuclide": self.nuclide, "interaction": self.interaction, "states": [list(s) for s in self.states],
                "transitions": list(self.transitions),
                "eff_charge": None if self.eff_charge is None else list(self.eff_charge),
                "gl": None if self.gl is None else list(self.gl), "gs": None if self.gs is None else list(self.gs),
                "threads": self.threads, "allow_over_limit": self.allow_over_limit}

    @classmethod
    def from_dict(cls, d: dict) -> ShellModelRequest:
        return cls(d["nuclide"], d["interaction"], tuple(tuple(s) for s in d.get("states", ())),
                   tuple(d.get("transitions", ("E2", "M1"))), d.get("eff_charge"), d.get("gl"), d.get("gs"),
                   d.get("threads"), bool(d.get("allow_over_limit", False)))


@dataclass
class ShellState:
    """One shell-model state. ``index`` is 1 for the lowest of its 2J±; ``energy_mev`` is KSHELL's total energy
    and ``ex_mev`` the excitation energy above the lowest state of the result. ``q_s`` (e fm²) and ``mu`` (μ_N) are
    ``None`` where not computed (J = 0, J = 1/2 for Q_s, or an M = 0 run for μ). The proton and neutron parts
    (``q_p``, ``q_n`` without effective charges; ``l_p``, ``l_n``, ``s_p``, ``s_n``) let Q_s and μ be worked out
    again for other charges and g-factors."""

    j2: int
    parity: int
    index: int
    energy_mev: float
    ex_mev: float = 0.0
    t2: Optional[int] = None
    q_s: Optional[float] = None
    mu: Optional[float] = None
    q_p: Optional[float] = None
    q_n: Optional[float] = None
    l_p: Optional[float] = None
    l_n: Optional[float] = None
    s_p: Optional[float] = None
    s_n: Optional[float] = None

    @property
    def key(self) -> State:
        return self.j2, self.parity, self.index

    @property
    def spin(self) -> float:
        return self.j2 / 2

    @property
    def label(self) -> str:
        return label(self.j2, self.parity, self.index)


@dataclass
class ShellTransition:
    """A transition from ``initial`` (the higher state) to ``final``, each (2J, parity, index). ``mred`` is the
    signed ⟨final‖M‖initial⟩ (e fm^λ or μ_N); ``b_down`` = B(initial → final) = mred² / (2J_i + 1) and ``b_up``
    = B(final → initial) = mred² / (2J_f + 1), in e² fm^{2λ} or μ_N². For a state with itself (initial = final),
    ``moment`` is its Q_s (E2, e fm²) or μ (M1, μ_N).

    Straight out of :func:`parse_transit_log`, the indices count states within their wavefunction file
    (``waves`` names the files of initial and final) and ``e_initial_mev``/``e_final_mev`` hold the energies
    printed; :func:`read_logs` replaces them with indices among each 2J± and clears ``waves``."""

    kind: str
    initial: State
    final: State
    mred: float
    b_down: Optional[float] = None
    b_up: Optional[float] = None
    moment: Optional[float] = None
    e_initial_mev: Optional[float] = None
    e_final_mev: Optional[float] = None
    waves: tuple = ()

    @property
    def diagonal(self) -> bool:
        return self.initial == self.final


@dataclass
class ShellModelResult:
    """The states and transitions of one calculation, with what produced them: the nuclide, interaction and
    valence space, the effective charges (e_p, e_n) and g-factors (g_l and g_s, each (proton, neutron)) used,
    the M-scheme dimensions (``dimensions[(2M, parity)]``, with 2M = 2J for KSHELL's J-projected runs), and the
    KSHELL version and time taken. ``notes`` say in plain words what the reader had to work around."""

    request: Optional[ShellModelRequest]
    space: str
    interaction: str
    kshell_version: str
    eff_charge: Optional[tuple]
    gl: Optional[tuple]
    gs: Optional[tuple]
    states: list = field(default_factory=list)
    transitions: list = field(default_factory=list)
    dimensions: dict = field(default_factory=dict)
    elapsed_s: float = 0.0
    notes: list = field(default_factory=list)
    nuclide: str = ""

    def state(self, j2: int, parity: int, index: int = 1) -> Optional[ShellState]:
        for s in self.states:
            if s.key == (j2, parity, index):
                return s
        return None

    def transition(self, kind: str, a: State, b: State) -> Optional[ShellTransition]:
        """The ``kind`` transition between states ``a`` and ``b``, in either order (it is stored from the higher
        state to the lower one)."""
        a, b = tuple(a), tuple(b)
        for t in self.transitions:
            if t.kind == kind and {t.initial, t.final} == {a, b} and (a == b) == t.diagonal:
                return t
        return None

    def to_dict(self) -> dict:
        """The result as JSON-safe values (lists, numbers, text, ``None``)."""

        def st(s: ShellState) -> dict:
            return dict(vars(s))

        def tr(t: ShellTransition) -> dict:
            return {"kind": t.kind, "initial": list(t.initial), "final": list(t.final), "mred": t.mred,
                    "b_down": t.b_down, "b_up": t.b_up, "moment": t.moment, "e_initial_mev": t.e_initial_mev,
                    "e_final_mev": t.e_final_mev}

        return {"nuclide": self.nuclide, "request": None if self.request is None else self.request.to_dict(),
                "space": self.space, "interaction": self.interaction, "kshell_version": self.kshell_version,
                "eff_charge": None if self.eff_charge is None else list(self.eff_charge),
                "gl": None if self.gl is None else list(self.gl), "gs": None if self.gs is None else list(self.gs),
                "states": [st(s) for s in self.states], "transitions": [tr(t) for t in self.transitions],
                "dimensions": [[m2, p, d] for (m2, p), d in sorted(self.dimensions.items())],
                "elapsed_s": self.elapsed_s, "notes": list(self.notes)}

    @classmethod
    def from_dict(cls, d: dict) -> ShellModelResult:
        def tup(x):
            return None if x is None else tuple(x)

        return cls(
            None if d.get("request") is None else ShellModelRequest.from_dict(d["request"]),
            d.get("space", ""), d.get("interaction", ""), d.get("kshell_version", "unknown"),
            tup(d.get("eff_charge")), tup(d.get("gl")), tup(d.get("gs")),
            [ShellState(**s) for s in d.get("states", [])],
            [ShellTransition(t["kind"], tuple(t["initial"]), tuple(t["final"]), t["mred"], t.get("b_down"),
                             t.get("b_up"), t.get("moment"), t.get("e_initial_mev"), t.get("e_final_mev"))
             for t in d.get("transitions", [])],
            {(int(m2), int(p)): int(n) for m2, p, n in d.get("dimensions", [])},
            float(d.get("elapsed_s", 0.0)), list(d.get("notes", [])), d.get("nuclide", ""))


# ---------------------------------------------------------------------------------------------------------------
# Reading the logs


def _float(text: str) -> Optional[float]:
    """A Fortran number, or ``None`` for a field it overflowed (``*****``) or anything else unreadable."""
    try:
        return float(text.replace("D", "E").replace("d", "e"))
    except (ValueError, AttributeError):
        return None


def _namelist(lines: list) -> dict:
    """The ``&INPUT … /`` block KSHELL echoes at the top of a log, keys upper case, values as text (strings
    without their quotes and padding)."""
    out, inside = {}, False
    for line in lines:
        s = line.strip()
        if s.upper().startswith("&INPUT"):
            inside = True
            continue
        if inside and (s == "/" or s.upper().startswith("&END")):
            break
        if inside:
            m = re.match(r"([A-Z_0-9]+)\s*=\s*(.*?),?\s*$", s, re.I)
            if m:
                value = m.group(2).strip()
                if value.startswith('"'):
                    value = value.strip('"').strip()
                out[m.group(1).upper()] = value
    return out


def _numbers(text: str) -> list:
    """The numbers in a namelist value such as ``1.5000000000000000     , 0.50000000000000000``; Fortran's
    repeat form ``2*0.0`` is expanded."""
    out = []
    for part in text.replace(",", " ").split():
        if "*" in part:
            n, _, x = part.partition("*")
            out += [_float(x)] * int(n)
        else:
            out.append(_float(part))
    return out


class KshellLog(list):
    """What :func:`parse_kshell_log` read: the states (it is a list of :class:`ShellState`) and, as attributes,
    the M-scheme ``dimension``, ``mtot`` (2M), ``parity``, ``eff_charge``, ``gl``, ``gs``, the interaction file
    ``fn_int``, the wavefunction file ``wave``, the valence nucleons ``z_val``/``n_val`` and core
    ``core_z``/``core_n``, the model-space ``orbits`` (names like ``"0d5/2"``), ``elapsed_s``, and ``notes``."""

    dimension: Optional[int] = None
    mtot: Optional[int] = None
    parity: Optional[int] = None
    eff_charge: Optional[tuple] = None
    gl: Optional[tuple] = None
    gs: Optional[tuple] = None
    fn_int: str = ""
    wave: str = ""
    z_val: Optional[int] = None
    n_val: Optional[int] = None
    core_z: Optional[int] = None
    core_n: Optional[int] = None
    elapsed_s: float = 0.0

    def __init__(self, states=()):
        super().__init__(states)
        self.orbits: list = []
        self.notes: list = []


_STATE_HEAD = re.compile(r"^\s*(\d+)\s+<H>:\s*(\S+)\s+<JJ>:\s*(\S+)\s+J:\s*(-?\d+)/2\s+prty\s+(-?\d+)")
_EXPECT = re.compile(r"<([A-Za-z]+)>:?\s+(\S+)")
_PARTS = {"Qp": "q_p", "Qn": "q_n", "eQ": "q_s", "Lp": "l_p", "Ln": "l_n", "sp": "s_p", "sn": "s_n", "gm": "mu"}


def parse_kshell_log(text: str) -> KshellLog:
    """The states in one ``kshell.exe`` log: energy, 2J, parity, 2T, Q_s, μ and their proton and neutron parts.
    Each state's index counts the states of its 2J± in this log (logs that are not J-projected, such as an M = 0
    run, mix J values), and ``ex_mev`` is measured from the lowest state of this log; :func:`read_logs` puts the
    states of several logs together. The M-scheme dimension and the header values are attributes of the list
    returned (see :class:`KshellLog`)."""
    lines = text.splitlines()
    out = KshellLog()
    nl = _namelist(lines)
    out.fn_int = nl.get("FN_INT", "")
    out.wave = nl.get("FN_SAVE_WAVE", "")
    if "MTOT" in nl:
        out.mtot = int(_numbers(nl["MTOT"])[0])
    for key, attr in (("EFF_CHARGE", "eff_charge"), ("GL", "gl"), ("GS", "gs")):
        if key in nl:
            values = _numbers(nl[key])
            if len(values) >= 2 and None not in values[:2]:
                setattr(out, attr, (values[0], values[1]))
    rows, current = [], None
    for line in lines:
        if m := re.search(r"total m-scheme dimension\s+(\d+)", line):
            out.dimension = int(m.group(1))
        elif m := re.search(r"parity\s*=\s*([+-])", line):
            out.parity = 1 if m.group(1) == "+" else -1
        elif m := re.search(r"valence protons and neutrons\s*=\s*(\d+)\s+(\d+).*n,z-core\s+(\d+)\s+(\d+)", line):
            out.z_val, out.n_val = int(m.group(1)), int(m.group(2))
            out.core_n, out.core_z = int(m.group(3)), int(m.group(4))
        elif m := re.match(r"\s*\d+\s+\d+\s+\d+\s+\d+\s+-?1\s+-?1\s+\d+\s+([pn])\s+(\d+)([a-z])_(\d+/2)", line):
            name = f"{m.group(2)}{m.group(3)}{m.group(4)}"
            if m.group(1) == "p":
                out.orbits.append(name)
        elif m := re.match(r"\s*effective charges\s+(\S+)\s+(\S+)", line):
            out.eff_charge = (_float(m.group(1)), _float(m.group(2)))
        elif m := re.match(r"\s*gl,gs\s*=\s*(.*)$", line):
            values = [_float(x) for x in m.group(1).split()]
            if len(values) == 4 and None not in values:
                out.gl, out.gs = (values[0], values[1]), (values[2], values[3])
        elif m := re.search(r"total time it took was:\s*(\S+)\s*sec", line):
            out.elapsed_s = _float(m.group(1)) or 0.0
        if m := _STATE_HEAD.match(line):
            energy = _float(m.group(2))
            jj, j2, parity = _float(m.group(3)), int(m.group(4)), int(m.group(5))
            if energy is None:
                out.notes.append(f"state {m.group(1)}: the energy was not printed (overflow); left out")
                current = None
                continue
            current = {"n": int(m.group(1)), "energy": energy, "j2": j2, "parity": parity, "jj": jj}
            rows.append(current)
            continue
        if current is None:
            continue
        if line.strip().startswith("---"):
            current = None
            continue
        if m := re.search(r"<TT>:\s*\S+\s+T:\s*(-?\d+)/2", line):
            current["t2"] = int(m.group(1))
            continue
        for key, value in _EXPECT.findall(line):
            if key in _PARTS:
                current[_PARTS[key]] = _float(value)
    if out.parity is None and rows:
        out.parity = rows[0]["parity"]
    lowest = min((r["energy"] for r in rows), default=0.0)
    count: dict = {}
    for r in sorted(rows, key=lambda r: (r["energy"], r["n"])):
        k = (r["j2"], r["parity"])
        count[k] = count.get(k, 0) + 1
        r["index"] = count[k]
        j = r["j2"] / 2
        if r["jj"] is not None and abs(r["jj"] - j * (j + 1)) > 0.01:
            out.notes.append(f"state {r['n']} ({label(r['j2'], r['parity'], r['index'])}): <J²> = {r['jj']:g} is "
                             f"not J(J+1) = {j * (j + 1):g}; it may not be converged or of good J")
    for r in sorted(rows, key=lambda r: r["n"]):
        out.append(ShellState(r["j2"], r["parity"], r["index"], r["energy"], r["energy"] - lowest, r.get("t2"),
                              *(r.get(a) for a in ("q_s", "mu", "q_p", "q_n", "l_p", "l_n", "s_p", "s_n"))))
    return out


class TransitLog(list):
    """What :func:`parse_transit_log` read: the transitions (it is a list of :class:`ShellTransition`) and, as
    attributes, the wavefunction files ``wave_l`` and ``wave_r``, their parities ``parity_l``/``parity_r``,
    ``eff_charge``, ``gl``, ``gs``, ``fn_int``, ``z``, ``n``, ``elapsed_s`` and ``notes``."""

    wave_l: str = ""
    wave_r: str = ""
    parity_l: Optional[int] = None
    parity_r: Optional[int] = None
    eff_charge: Optional[tuple] = None
    gl: Optional[tuple] = None
    gs: Optional[tuple] = None
    fn_int: str = ""
    z: Optional[int] = None
    n: Optional[int] = None
    elapsed_s: float = 0.0

    def __init__(self, transitions=()):
        super().__init__(transitions)
        self.notes: list = []
        self.kinds: list = []


_TABLE_HEAD = re.compile(r"^\s*([EM]\d)\s+transition\b(.*)$")
#: Column edges of one table line, Fortran format (2(i2,"(",i4,")",f9.3), f8.3, 5f10.4).
_COLUMNS = {"j2_l": (0, 2), "n_l": (3, 7), "e_l": (8, 17), "j2_r": (17, 19), "n_r": (20, 24), "e_r": (25, 34),
            "ex": (34, 42), "mred": (42, 52), "b_lr": (52, 62), "b_rl": (62, 72), "mom": (72, 82)}


def parse_transit_log(text: str) -> TransitLog:
    """The E1, E2 and M1 tables of one ``transit.exe`` log. Each line becomes a :class:`ShellTransition` from the
    higher state to the lower one, with ⟨final‖M‖initial⟩ (turned round from KSHELL's ⟨left‖M‖right⟩ where
    needed), B↓ and B↑ (KSHELL's printed values; worked out from Mred where a column overflowed) and, on a
    diagonal line, the moment. A pair printed twice (both orders, when left and right are the same file) is
    kept once. Indices count within each wavefunction file until :func:`read_logs` renumbers them."""
    lines = text.splitlines()
    out = TransitLog()
    nl = _namelist(lines)
    out.fn_int = nl.get("FN_INT", "")
    out.wave_l, out.wave_r = nl.get("FN_LOAD_WAVE_L", ""), nl.get("FN_LOAD_WAVE_R", "")
    for key, attr in (("EFF_CHARGE", "eff_charge"), ("GL", "gl"), ("GS", "gs")):
        if key in nl:
            values = _numbers(nl[key])
            if len(values) >= 2 and None not in values[:2]:
                setattr(out, attr, (values[0], values[1]))
    seen: dict = {}
    i = 0
    while i < len(lines):
        line = lines[i]
        if m := re.match(r"\s*(left|right)\s+Z,N,A,M,prty:\s*(\d+)\s+(\d+)\s+\d+\s+-?\d+\s+(-?\d+)", line):
            setattr(out, f"parity_{m.group(1)[0]}", int(m.group(4)))
            out.z, out.n = int(m.group(2)), int(m.group(3))
        elif m := re.search(r"total elapsed time:\s*(\S+)", line):
            out.elapsed_s = _float(m.group(1)) or 0.0
        m = _TABLE_HEAD.match(line)
        if not m:
            i += 1
            continue
        kind, rest = m.group(1), m.group(2)
        parities = re.search(r"parity\s+(-?\d+)\s+(-?\d+)", rest)
        p_l, p_r = (int(parities.group(1)), int(parities.group(2))) if parities else (out.parity_l, out.parity_r)
        if kind.startswith("E") and (q := re.search(r"eff_charge=\s*(\S+)\s+(\S+)", rest)):
            charges = (_float(q.group(1)), _float(q.group(2)))
            if kind == "E2":
                out.eff_charge = charges
        if kind == "M1" and (g := re.search(r"gl,gs=\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", rest)):
            out.gl, out.gs = (_float(g.group(1)), _float(g.group(2))), (_float(g.group(3)), _float(g.group(4)))
        header = lines[i + 1] if i + 1 < len(lines) else ""
        i += 2
        if not header.strip().startswith("2xJf"):
            out.notes.append(f"a {kind} table with columns '{header.strip()}' was not read")
            continue
        out.kinds.append(kind)
        while i < len(lines) and lines[i].strip():
            row = lines[i].ljust(82)
            i += 1
            c = {k: row[a:b] for k, (a, b) in _COLUMNS.items()}
            try:
                j2_l, n_l, j2_r, n_r = (int(c[k]) for k in ("j2_l", "n_l", "j2_r", "n_r"))
            except ValueError:
                out.notes.append(f"{kind}: a line that could not be read was left out: '{row.strip()}'")
                continue
            e_l, e_r, mred = _float(c["e_l"]), _float(c["e_r"]), _float(c["mred"])
            if mred is None or e_l is None or e_r is None:
                out.notes.append(f"{kind}: a number overflowed KSHELL's column; left out: '{row.strip()}'")
                continue
            b_lr, b_rl, mom = _float(c["b_lr"]), _float(c["b_rl"]), _float(c["mom"])
            if b_lr is None:
                b_lr = mred**2 / (j2_l + 1)
            if b_rl is None:
                b_rl = mred**2 / (j2_r + 1)
            same = out.wave_l == out.wave_r and n_l == n_r and j2_l == j2_r and p_l == p_r and e_l == e_r
            left, right = (j2_l, p_l, n_l), (j2_r, p_r, n_r)
            if same or e_r >= e_l:  # right is the higher (initial) state: ⟨L‖M‖R⟩ is already ⟨f‖M‖i⟩
                t = ShellTransition(kind, right, left, mred, b_rl, b_lr, mom if same else None, e_r, e_l,
                                    (out.wave_r, out.wave_l))
            else:
                phase = -1.0 if ((j2_l - j2_r) // 2) % 2 else 1.0
                t = ShellTransition(kind, left, right, phase * mred, b_lr, b_rl, None, e_l, e_r,
                                    (out.wave_l, out.wave_r))
            key = (kind, t.initial, t.final, t.waves)
            if key in seen:
                if abs(seen[key].mred - t.mred) > 2e-4:
                    out.notes.append(f"{kind} {label(*t.initial)} → {label(*t.final)}: printed twice with "
                                     f"different values ({seen[key].mred:g}, {t.mred:g}); the first kept")
                continue
            seen[key] = t
            out.append(t)
    return out


# ---------------------------------------------------------------------------------------------------------------
# Putting the logs of one calculation together


def _same(a: Optional[tuple], b: Optional[tuple]) -> bool:
    return a is None or b is None or all(abs(x - y) < 5e-4 for x, y in zip(a, b))


def read_logs(kshell_logs: Iterable[str], transit_logs: Iterable[str] = (), request: Optional[ShellModelRequest] = None,
              kshell_version: str = "unknown") -> ShellModelResult:
    """One result from the texts of kshell.exe and transit.exe logs. States are numbered within each 2J± by
    energy across all kshell logs; each transit line is tied to its states through the wavefunction file it
    names (``FN_LOAD_WAVE_*`` against ``FN_SAVE_WAVE``) or, failing that, by 2J, parity and energy. A state found
    only in a transit log is kept with the energy printed there (3 decimals).

    Q_s and μ come from transit's diagonal lines when there are some (4 decimals, with the transit's charges);
    otherwise from the kshell log, worked out again from the proton and neutron parts if the transit logs used
    other effective charges or g-factors. Raises ``ValueError`` if the transit logs disagree on the charges."""
    klogs = [parse_kshell_log(t) for t in kshell_logs]
    tlogs = [parse_transit_log(t) for t in transit_logs]
    notes: list = []
    for log in klogs + tlogs:
        notes += log.notes

    # -- the charges and g-factors the transitions were computed with
    def agreed(name: str, logs: list) -> Optional[tuple]:
        values = [getattr(log, name) for log in logs if getattr(log, name) is not None]
        for v in values[1:]:
            if not _same(values[0], v):
                raise ValueError(f"the transit logs were run with different {name.replace('_', ' ')}s "
                                 f"({values[0]} and {v}); read each set of logs on its own")
        return values[0] if values else None

    e_logs = [log for log in tlogs if "E2" in log.kinds]
    m_logs = [log for log in tlogs if "M1" in log.kinds]
    eff = agreed("eff_charge", e_logs) or next((k.eff_charge for k in klogs if k.eff_charge), None)
    gl = agreed("gl", m_logs) or next((k.gl for k in klogs if k.gl), None)
    gs = agreed("gs", m_logs) or next((k.gs for k in klogs if k.gs), None)
    if eff is None and request is not None:
        eff = request.eff_charge
    if request is not None and request.eff_charge is not None and not _same(request.eff_charge, eff):
        notes.append(f"the request asked for effective charges {request.eff_charge}, the logs used {eff}")

    # -- the states, numbered within each 2J± by energy
    unique: list = []  # (state, wavefunction file, position in that file: KSHELL counts every state, any J)
    for k, log in enumerate(klogs):
        for n, s in enumerate(log, start=1):
            if any(u[0].j2 == s.j2 and u[0].parity == s.parity and abs(u[0].energy_mev - s.energy_mev) < 1e-5
                   and u[3] != k for u in unique):
                notes.append(f"{label(s.j2, s.parity)} at {s.energy_mev:.5f} MeV is in two kshell logs; kept once")
                continue
            unique.append((s, log.wave, n, k))

    by_file = {(u[1], u[2]): u[0] for u in unique if u[1]}

    def find(j2: int, parity: int, energy: float, wave: str, n: int) -> ShellState:
        s = by_file.get((wave, n))
        if s is not None and s.j2 == j2 and s.parity == parity and abs(s.energy_mev - energy) <= ENERGY_TOLERANCE_MEV:
            return s
        near = [u[0] for u in unique if u[0].j2 == j2 and u[0].parity == parity
                and abs(u[0].energy_mev - energy) <= ENERGY_TOLERANCE_MEV]
        if near:
            return min(near, key=lambda u: abs(u.energy_mev - energy))
        made = ShellState(j2, parity, 0, energy)
        unique.append((made, wave, n, -1))
        if wave:
            by_file[(wave, n)] = made
        transit_only.append(f"{label(j2, parity)} at {energy:.3f} MeV")
        return made

    raw, transit_only = [], []
    for log in tlogs:
        for t in log:
            si = find(t.initial[0], t.initial[1], t.e_initial_mev, t.waves[0], t.initial[2])
            sf = si if t.diagonal else find(t.final[0], t.final[1], t.e_final_mev, t.waves[1], t.final[2])
            raw.append((t, si, sf))

    if transit_only:
        notes.append("known only from the transit logs (no kshell log for them was given), with the energy printed "
                     "there: " + ", ".join(transit_only))
    states = [u[0] for u in unique]
    lowest = min((s.energy_mev for s in states), default=0.0)
    count: dict = {}
    ordered = sorted(states, key=lambda s: s.energy_mev)
    renumber: dict = {}
    for s in ordered:
        k = (s.j2, s.parity)
        count[k] = count.get(k, 0) + 1
        renumber[id(s)] = count[k]
    final_states = []
    for s in ordered:
        final_states.append(_replace(s, index=renumber[id(s)], ex_mev=s.energy_mev - lowest))
    new_state = {id(s): f for s, f in zip(ordered, final_states)}

    transitions, seen = [], {}
    for t, si, sf in raw:
        a, b = new_state[id(si)], new_state[id(sf)]
        nt = _replace(t, initial=a.key, final=b.key, e_initial_mev=a.energy_mev, e_final_mev=b.energy_mev, waves=())
        if not nt.diagonal and a.energy_mev < b.energy_mev:  # the kshell energies put them the other way round
            phase = -1.0 if ((a.j2 - b.j2) // 2) % 2 else 1.0
            nt = _replace(nt, initial=b.key, final=a.key, mred=phase * nt.mred, b_down=t.b_up, b_up=t.b_down,
                          e_initial_mev=b.energy_mev, e_final_mev=a.energy_mev)
        key = (nt.kind, nt.initial, nt.final)
        if key in seen:
            continue
        seen[key] = nt
        transitions.append(nt)
    transitions.sort(key=lambda t: (KINDS.index(t.kind) if t.kind in KINDS else 9, -t.e_initial_mev,
                                    -t.e_final_mev))

    # -- moments: transit's diagonal lines first, else the kshell log (with the transit's charges)
    moments = {(t.kind, t.initial): t.moment for t in transitions if t.diagonal and t.moment is not None}
    k_eff = next((k.eff_charge for k in klogs if k.eff_charge), None)
    k_g = next(((k.gl, k.gs) for k in klogs if k.gl), None)
    redo_q = eff is not None and k_eff is not None and not _same(eff, k_eff)
    redo_mu = k_g is not None and gl is not None and gs is not None and not (_same(gl, k_g[0]) and _same(gs, k_g[1]))
    recomputed = set()
    for i, s in enumerate(final_states):
        q, mu = s.q_s, s.mu
        if ("E2", s.key) in moments:
            q = moments[("E2", s.key)]
        elif redo_q and s.q_p is not None and s.q_n is not None:
            q = eff[0] * s.q_p + eff[1] * s.q_n
            recomputed.add("Q_s")
        elif redo_q:
            q = None
        if ("M1", s.key) in moments:
            mu = moments[("M1", s.key)]
        elif redo_mu and None not in (s.l_p, s.l_n, s.s_p, s.s_n):
            mu = gl[0] * s.l_p + gl[1] * s.l_n + gs[0] * s.s_p + gs[1] * s.s_n
            recomputed.add("μ")
        elif redo_mu:
            mu = None
        final_states[i] = _replace(s, q_s=q, mu=mu)
    if recomputed:
        what = " and ".join(sorted(recomputed))
        notes.append(f"{what} of states without a diagonal transit line: worked out again from the proton and "
                     f"neutron parts in the kshell logs, with the charges and g-factors of the transit logs")

    # -- what produced it
    first = klogs[0] if klogs else None
    fn_int = next((log.fn_int for log in klogs + tlogs if log.fn_int), "")
    interaction = request.interaction if request is not None else Path(fn_int).stem.lower()
    nuclide = request.nuclide if request is not None else ""
    if not nuclide and first is not None and None not in (first.z_val, first.n_val, first.core_z, first.core_n):
        z, n = first.core_z + first.z_val, first.core_n + first.n_val
        nuclide = f"{z + n}{SYMBOLS[z]}" if 0 < z < len(SYMBOLS) else ""
    if not nuclide and tlogs and tlogs[0].z:
        z = tlogs[0].z
        nuclide = f"{z + tlogs[0].n}{SYMBOLS[z]}"
    space = ""
    if first is not None and first.orbits:
        space = " ".join(first.orbits)
        if first.core_z is not None and first.core_n is not None:
            cz, cn = first.core_z, first.core_n
            space += f" outside {cz + cn}{SYMBOLS[cz]}" if 0 < cz < len(SYMBOLS) else ""
    dims = {(log.mtot, log.parity): log.dimension for log in klogs
            if log.dimension is not None and log.mtot is not None and log.parity is not None}
    elapsed = sum(log.elapsed_s for log in klogs + tlogs)
    return ShellModelResult(request, space, interaction, kshell_version, eff, gl, gs, final_states, transitions,
                            dims, elapsed, notes, nuclide)


def read_folder(folder: Union[str, Path], request: Optional[ShellModelRequest] = None,
                kshell_version: Optional[str] = None) -> ShellModelResult:
    """The result of the KSHELL logs in ``folder`` (``*.txt`` and ``*.log``; kshell.exe and transit.exe logs are
    told apart by their contents, anything else such as ``collect_logs.py``'s summary is skipped). The KSHELL
    version is read from a ``VERSION`` file in the folder when not given."""
    folder = Path(folder)
    kshell, transit = [], []
    for path in sorted(folder.iterdir()):
        if not path.is_file() or path.suffix.lower() not in (".txt", ".log"):
            continue
        text = path.read_bytes().decode("utf-8", errors="replace")
        if "FN_LOAD_WAVE_L" in text:
            transit.append(text)
        elif "<H>:" in text or re.search(r"^\s*FN_PTN\s*=", text, re.M):
            kshell.append(text)
    if not kshell and not transit:
        raise ValueError(f"{folder} has no kshell.exe or transit.exe logs")
    if kshell_version is None:
        version_file = folder / "VERSION"
        kshell_version = version_file.read_text(encoding="utf-8").strip() if version_file.is_file() else "unknown"
    return read_logs(kshell, transit, request, kshell_version)


# ---------------------------------------------------------------------------------------------------------------
# Into a level scheme


@dataclass
class Match:
    """A shell-model state and the known level it was paired with (``level``, an index into the scheme's
    levels), or ``None`` when the scheme has no level for it: such a state can be added as a new level."""

    state: ShellState
    level: Optional[int]
    known_kev: Optional[float] = None

    @property
    def shell_kev(self) -> float:
        return 1000.0 * self.state.ex_mev

    @property
    def delta_kev(self) -> Optional[float]:
        """Shell model minus known, keV."""
        return None if self.known_kev is None else self.shell_kev - self.known_kev


def match(result: ShellModelResult, scheme: LevelScheme) -> list:
    """Pair shell-model states with the scheme's levels: within each 2J±, the n-th shell-model state with the
    n-th known level of that J and parity in energy order (levels whose J or parity is not known are not paired).
    One :class:`Match` per state, in order of excitation energy."""
    known: dict = {}
    for i, lv in enumerate(scheme.levels):
        if lv.spin is None or lv.parity is None:
            continue
        known.setdefault((round(2 * lv.spin), lv.parity), []).append(i)
    for k in known:
        known[k].sort(key=lambda i: scheme.levels[i].energy.value)
    out = []
    for s in sorted(result.states, key=lambda s: (s.ex_mev, s.j2)):
        levels = known.get((s.j2, s.parity), [])
        i = levels[s.index - 1] if s.index <= len(levels) else None
        out.append(Match(s, i, None if i is None else scheme.levels[i].energy.value))
    return out


#: Quantities ``apply`` understands, with the spellings accepted.
_QUANTITIES = {"energy": "energy", "energies": "energy", "e1": "E1", "e2": "E2", "q": "Q", "moments": "Q",
               "qs": "Q", "m1": "M1", "mu": "mu"}


def _choice(which: Union[str, Iterable[str], None]) -> tuple:
    """(quantities, states or None for all, add every unmatched state) from ``which``."""
    items = [which] if isinstance(which, str) or which is None else list(which)
    quantities, states, new = set(), set(), False
    for item in items:
        if item is None or str(item).strip().lower() == "all":
            continue
        text = str(item).strip()
        if text.lower() == "new":
            new = True
        elif text.lower() in _QUANTITIES:
            quantities.add(_QUANTITIES[text.lower()])
        else:
            try:
                j2, parity, index = parse_label(text)
            except ValueError:
                raise ValueError(f"cannot apply '{text}': name a state such as '2+1', a quantity "
                                 f"({', '.join(sorted(set(_QUANTITIES.values())))}), 'new' or 'all'") from None
            states.add((j2, parity, 1 if index is None else index))
    if not quantities:
        quantities = {"energy", "E1", "E2", "Q"}
    return quantities, (states or None), new


def _renumber(scheme: LevelScheme) -> dict:
    """Sort the scheme's levels by energy, renumbering transitions and matrix elements; returns old → new."""
    order = sorted(range(len(scheme.levels)), key=lambda i: (scheme.levels[i].energy.value, i))
    new = {old: k for k, old in enumerate(order)}
    if order == list(range(len(order))):
        return new
    scheme.levels = [scheme.levels[i] for i in order]
    scheme.transitions = [_replace(t, initial=new[t.initial], final=new[t.final]) for t in scheme.transitions]
    scheme.matrix_elements = [MatrixElement(min(new[m.a], new[m.b]), max(new[m.a], new[m.b]), m.multipolarity,
                                            m.value) for m in scheme.matrix_elements]
    return new


def _note(result: ShellModelResult, charges: bool) -> str:
    text = f"KSHELL shell model, {result.interaction or 'interaction not named'}"
    if charges and result.eff_charge is not None:
        text += f", e_p = {result.eff_charge[0]:g} e, e_n = {result.eff_charge[1]:g} e"
    return text


def apply(result: ShellModelResult, scheme: LevelScheme, which: Union[str, Iterable[str], None] = "all",
          replace: bool = False) -> list:
    """Put chosen shell-model values into ``scheme`` (changed in place) with the source ``"shell model"`` and a
    note naming the interaction and effective charges. Returns what changed, one line each.

    ``which`` is ``"all"`` or a list of any of: quantities — ``"energy"`` (level energies, from the shell-model
    excitation energy), ``"E2"`` and ``"E1"`` (matrix elements between levels), ``"Q"`` (diagonal E2 matrix
    elements, i.e. quadrupole moments); states such as ``"2+1"`` (only values that involve these states); and
    ``"new"`` (add every state the scheme has no level for). Naming an unmatched state adds it as a new level
    too. With no quantity named, all are applied; with no state named, all matched states.

    Precedence: an ENSDF, derived or user value is never replaced unless ``replace=True``; values marked
    ``"assumed"`` or ``"shell model"`` are. The ground state keeps 0 keV. B(M1) and μ are reported but not stored,
    as a level scheme keeps E1, E2 and E3 matrix elements only.
    """
    if result.nuclide and scheme.nuclide:
        if parse_nuclide(result.nuclide) != parse_nuclide(scheme.nuclide):
            raise ValueError(f"the shell-model result is for {result.nuclide}, the level scheme for {scheme.nuclide}")
    quantities, chosen, add_all = _choice(which)
    changes: list = []
    kept: list = []
    pairs = match(result, scheme)

    def wanted(s: ShellState) -> bool:
        return chosen is None or s.key in chosen

    # -- energies of matched levels, then new levels
    level_of: dict = {}  # state key -> level index; levels are identified by object until renumbered
    for p in pairs:
        if p.level is not None:
            level_of[p.state.key] = scheme.levels[p.level]
    if "energy" in quantities:
        for p in pairs:
            if p.level is None or not wanted(p.state):
                continue
            lv = scheme.levels[p.level]
            e_kev = p.shell_kev
            if p.level == 0 or e_kev <= 0:
                continue  # the ground state stays at 0 keV, and no other level goes there
            if lv.energy.source in PROTECTED and not replace:
                kept.append((lv.energy.source, "energy", f"level {p.level}"))
                continue
            if abs(lv.energy.value - e_kev) < 1e-9 and lv.energy.source == SOURCE:
                continue
            old = lv.energy
            lv.energy = Value(e_kev, None, SOURCE, f"{_note(result, False)}: {p.state.label}")
            changes.append(f"level {p.level} ({lv.label}): energy {old.value:g} keV ({old.source}) → {e_kev:.1f} keV "
                           f"(shell model, {p.state.label})")
    for p in pairs:
        if p.level is not None or not (add_all or (chosen is not None and p.state.key in chosen)):
            continue
        e_kev = p.shell_kev
        if e_kev <= 0:
            changes.append(f"{p.state.label} was not added: it is the shell model's lowest state, but the scheme's "
                           f"ground state has another J±")
            continue
        if any(abs(lv.energy.value - e_kev) < 1e-6 for lv in scheme.levels):
            changes.append(f"{p.state.label} was not added: a level at {e_kev:.1f} keV is already there")
            continue
        lv = Level(Value(e_kev, None, SOURCE, f"{_note(result, False)}: {p.state.label}"), p.state.spin,
                   p.state.parity, label(p.state.j2, p.state.parity))
        scheme.levels.append(lv)
        level_of[p.state.key] = lv
        changes.append(f"new level {lv.label} at {e_kev:.1f} keV (shell model, {p.state.label})")
    _renumber(scheme)
    index = {id(lv): i for i, lv in enumerate(scheme.levels)}
    level_index = {k: index[id(lv)] for k, lv in level_of.items()}

    # -- matrix elements
    note = _note(result, True)
    states = {s.key: s for s in result.states}
    done = set()
    for t in result.transitions:
        if t.kind not in ELECTRIC:
            continue
        if t.diagonal and not (t.kind == "E2" and "Q" in quantities):
            continue
        if not t.diagonal and t.kind not in quantities:
            continue
        si, sf = states.get(t.initial), states.get(t.final)
        if si is None or sf is None or not (wanted(si) or wanted(sf)):
            continue
        if t.initial not in level_index or t.final not in level_index:
            continue
        a, b = level_index[t.final], level_index[t.initial]
        # levels.py keeps ⟨higher level‖M‖lower level⟩; mred is ⟨final‖M‖initial⟩
        value = t.mred if a > b or t.diagonal else t.mred * (-1.0 if ((si.j2 - sf.j2) // 2) % 2 else 1.0)
        done.add((t.kind, min(a, b), max(a, b)))
        what = f"{t.kind} {sf.label} ← {si.label}" if not t.diagonal else f"Q_s of {si.label}"
        changes += _set(scheme, min(a, b), max(a, b), t.kind, value, f"{note}: {what}", replace, kept)
    if "Q" in quantities:  # states with Q_s from the kshell log but no diagonal transit line
        for s in result.states:
            k = s.key
            if s.q_s is None or k not in level_index or not wanted(s) or quadrupole_factor(s.spin) == 0:
                continue
            i = level_index[k]
            if ("E2", i, i) in done:
                continue
            changes += _set(scheme, i, i, "E2", s.q_s / quadrupole_factor(s.spin),
                            f"{note}: Q_s of {s.label} = {s.q_s:+.3f} e fm²", replace, kept)
    if quantities & {"M1", "mu"}:
        changes.append("B(M1) and magnetic moments stay in the shell-model result: a level scheme keeps E1, E2 and "
                       "E3 matrix elements only")
    if kept:
        counts: dict = {}
        for source, what, _ in kept:
            counts[(source, what)] = counts.get((source, what), 0) + 1
        parts = [f"{n} {source} {what if n == 1 else ('energies' if what == 'energy' else what + 's')}"
                 for (source, what), n in counts.items()]
        detail = "" if len(kept) > 6 else " (" + ", ".join(f"{w} of {where}" for _, w, where in kept) + ")"
        changes.append(f"kept {', '.join(parts)}{detail}; apply with replace=True to overwrite them")
    return changes


def _set(scheme: LevelScheme, a: int, b: int, kind: str, value: float, note: str, replace: bool,
         kept: list) -> list:
    existing = scheme.matrix_element(a, b, kind)
    lam = int(kind[1])
    where = f"⟨{b}‖{kind}‖{a}⟩" if a != b else f"⟨{a}‖{kind}‖{a}⟩"
    if existing is not None and existing.source in PROTECTED and not replace:
        kept.append((existing.source, f"{kind} matrix element", where))
        return []
    if existing is not None and existing.source == SOURCE and existing.note == note \
            and abs(existing.value - value) < 1e-12:
        return []
    scheme.set_matrix_element(a, b, kind, value, None, SOURCE, note)
    before = "" if existing is None else f" (was {existing.value:+.4g}, {existing.source})"
    return [f"{where} = {value:+.4g} e fm^{lam}{before}: {note.split(': ', 1)[-1]}"]


def to_level_scheme(result: ShellModelResult, max_energy_kev: Optional[float] = None) -> LevelScheme:
    """A level scheme made of the shell-model values alone (for a nucleus ENSDF knows little about): every state
    up to ``max_energy_kev`` as a level, with its E1/E2 matrix elements and quadrupole moments."""
    if not result.states:
        raise ValueError("the shell-model result has no states")
    lowest = min(result.states, key=lambda s: s.energy_mev)
    scheme = LevelScheme(result.nuclide or "", reference=_note(result, True))
    scheme.levels.append(Level(Value(0.0, None, SOURCE, f"{_note(result, False)}: {lowest.label}"), lowest.spin,
                               lowest.parity, label(lowest.j2, lowest.parity)))
    keep = [s.label for s in result.states if s is not lowest
            and (max_energy_kev is None or 1000.0 * s.ex_mev <= max_energy_kev)]
    apply(result, scheme, ["energy", "E1", "E2", "Q", lowest.label] + keep)
    return scheme


__all__ = ["ENERGY_TOLERANCE_MEV", "KINDS", "KshellLog", "Match", "PROTECTED", "SOURCE", "ShellModelRequest",
           "ShellModelResult", "ShellState", "ShellTransition", "TransitLog", "apply", "label", "match",
           "parse_kshell_log", "parse_label", "parse_transit_log", "read_folder", "read_logs", "to_level_scheme"]
