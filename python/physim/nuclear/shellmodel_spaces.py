"""Shell-model valence spaces and interactions, partition files, M-scheme dimensions and memory estimates.

physim runs KSHELL (N. Shimizu, T. Mizusaki, Y. Utsuno and Y. Tsunoda, Comput. Phys. Commun. 244, 372 (2019)) as
an outside program. This module holds what physim needs to know before a run::

    from physim.nuclear import shellmodel_spaces as sms

    sms.options("28Si")                       # [Option("usda", 6, 6, (1,), ...), Option("usdb", ...), ...]
    usdb = sms.INTERACTIONS["usdb"]
    sms.dimensions(usdb.space, 6, 6, +1)      # {0: 93710, 2: 90338, ..., 28: 1}  (2M -> M-scheme dimension)
    sms.j_dimension(usdb.space, 6, 6, +1, 0)  # 3372 J = 0 states
    text = sms.partitions(usdb.space, 6, 6, +1, snt_name="usdb.snt")   # the .ptn file KSHELL reads
    sms.estimate(request)                     # memory, disk, vector storage mode, whether it fits, rough time

Conventions (KSHELL's): orbits carry ``tz`` −1 for protons and +1 for neutrons, as in ``.snt`` files; angular
momenta are doubled (``j2`` = 2j, 2M, 2J); parity is +1 or −1.

The partition files and dimensions are physim's own (KSHELL's ``gen_partition.py`` is Python 2 only in the
release physim bundles), written in the format KSHELL reads and checked against ``gen_partition.py`` and
``count_dim.exe`` output in ``tests/reference/nuclear/shellmodel/partitions/``. Dimensions are counted by
combinatorics over orbit occupations, not by listing states, so even ⁵⁶Ni in the full fp shell (1 087 455 228
states at M = 0) takes milliseconds.

Memory rule (the user's): a run may use at most **half of the RAM free when it starts**, not of the RAM installed.
16 GB installed with 10 GB in use leaves 6 GB free, so the limit is 3 GB.
"""

from __future__ import annotations

import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Optional, Union

from .names import SYMBOLS, parse_nuclide

L_LETTERS = "spdfghiklmno"

#: KSHELL's default number of Lanczos vectors kept (``max_lanc_vec``) and restart vectors (``n_restart_vec``).
MAX_LANC_VEC = 200
N_RESTART_VEC = 10

#: Bytes per vector component: KSHELL keeps wave functions in double precision.
BYTES_PER_COMPONENT = 8

#: Ad hoc factor on top of the vectors, as KSHELL's ``count_dim`` uses for its own estimate.
MEMORY_FACTOR = 1.2

#: Fixed memory for the program, the interaction and the partition tables (KSHELL's logs show a few MB for sd and
#: fp spaces; 100 MB is a generous round number).
FIXED_OVERHEAD_BYTES = 100 * 1024**2

#: Wave 0 timing used to scale the rough time guess: ⁴⁸Cr with KB3, ground state, M-scheme dimension 1 963 461,
#: 72 s on a Ryzen 5 3600 with 12 threads (400 s on 1 thread: speed-up ≈ threads^0.69).
_TIME_REF_SECONDS = 72.0
_TIME_REF_DIMENSION = 1_963_461
_TIME_REF_THREADS = 12
_THREAD_EXPONENT = 0.69
_SECONDS_PER_RUN = 2.0        # start-up, reading the interaction, writing the wave functions
_SECONDS_PER_TRANSIT = 1.0    # transit.exe for one pair of J± sets (0.06-0.09 s for ²⁸Si; generous)


# ---------------------------------------------------------------------------------------------------------------
# Orbits and valence spaces


@dataclass(frozen=True)
class Orbit:
    """One single-particle orbit: radial quantum number ``n`` (from 0), ``l``, ``j2`` = 2j and ``tz`` (−1 p, +1 n)."""

    n: int
    l: int
    j2: int
    tz: int

    @property
    def capacity(self) -> int:
        """How many nucleons fit: 2j + 1."""
        return self.j2 + 1

    @property
    def parity(self) -> int:
        return -1 if self.l % 2 else 1

    @property
    def label(self) -> str:
        """``"p 0d5/2"``, as KSHELL prints it (without its underscore)."""
        return f"{'p' if self.tz < 0 else 'n'} {self.n}{L_LETTERS[self.l]}{self.j2}/2"


@dataclass(frozen=True)
class ValenceSpace:
    """A valence space: its orbits in ``.snt`` order (protons first) above an inert core of ``core_z``, ``core_n``.

    The order matters: partition files list occupations orbit by orbit in this order, so a space belongs to the
    ``.snt`` files that list their orbits the same way.
    """

    name: str
    core_z: int
    core_n: int
    orbits: tuple[Orbit, ...]

    def species(self, tz: int) -> tuple[Orbit, ...]:
        """The proton (``tz`` −1) or neutron (+1) orbits."""
        return tuple(o for o in self.orbits if o.tz == tz)

    @property
    def proton_orbits(self) -> tuple[Orbit, ...]:
        return self.species(-1)

    @property
    def neutron_orbits(self) -> tuple[Orbit, ...]:
        return self.species(1)

    def capacity(self, tz: int) -> int:
        """How many protons (``tz`` −1) or neutrons (+1) the space holds."""
        return sum(o.capacity for o in self.species(tz))

    def valence(self, z: int, a: int) -> tuple[int, int]:
        """Valence protons and neutrons of the nucleus (Z, A); either may be negative or too many (not checked)."""
        return z - self.core_z, (a - z) - self.core_n

    def holds(self, z_val: int, n_val: int) -> bool:
        return 0 <= z_val <= self.capacity(-1) and 0 <= n_val <= self.capacity(1)

    @property
    def core_label(self) -> str:
        """The core nucleus, as ``"16O"``."""
        return f"{self.core_z + self.core_n}{SYMBOLS[self.core_z]}"


def _orbits(core_z: int, core_n: int, name: str, spec: str) -> ValenceSpace:
    """A space from ``"0d3 0d5 1s1"`` (n, l letter, 2j), the same orbits for protons and neutrons."""
    single = []
    for word in spec.split():
        single.append((int(word[0]), L_LETTERS.index(word[1]), int(word[2:])))
    orbits = tuple(Orbit(n, l, j2, tz) for tz in (-1, 1) for n, l, j2 in single)
    return ValenceSpace(name, core_z, core_n, orbits)


def _data_lines(text: str) -> Iterator[list[str]]:
    """The words of each line of a ``.snt`` file, comments (``!`` or ``#`` onwards) removed, empty lines skipped."""
    for line in text.splitlines():
        words = []
        for word in line.split():
            if word[0] in "!#":
                break
            words.append(word)
        if words:
            yield words


def read_space(snt_path: Union[str, Path], name: Optional[str] = None) -> ValenceSpace:
    """The model-space block of a KSHELL ``.snt`` file; ``name`` defaults to the file's stem (``"usdb"``)."""
    path = Path(snt_path)
    lines = _data_lines(path.read_text(encoding="utf-8", errors="replace"))
    try:
        n_p, n_n, core_z, core_n = (int(w) for w in next(lines)[:4])
        orbits = []
        for k in range(n_p + n_n):
            words = next(lines)
            index, n, l, j2, tz = (int(w) for w in words[:5])
            if index != k + 1:
                raise ValueError(f"orbit {k + 1} is numbered {index}")
            if tz not in (-1, 1) or (tz == -1) != (k < n_p):
                raise ValueError("proton orbits (tz = -1) must come first, then neutron orbits (tz = +1)")
            orbits.append(Orbit(n, l, j2, tz))
    except (StopIteration, ValueError) as exc:
        raise ValueError(f"{path.name}: cannot read the model space: {exc}") from None
    return ValenceSpace(name or path.stem, core_z, core_n, tuple(orbits))


def read_snt_defaults(snt_path: Union[str, Path]) -> dict[str, Any]:
    """The ``!namelist key = value`` lines of a ``.snt`` file: the interaction's own recommended KSHELL inputs.

    Lists come back as tuples of floats (``eff_charge = 1.35, 0.35`` → ``(1.35, 0.35)``), single numbers as int or
    float.
    """
    out: dict[str, Any] = {}
    for line in Path(snt_path).read_text(encoding="utf-8", errors="replace").splitlines():
        m = re.match(r"^\s*!namelist\s+(\w+)\s*=\s*(.*?)\s*$", line)
        if not m:
            continue
        key, value = m.group(1), m.group(2)
        parts = [p for p in re.split(r"[,\s]+", value) if p]
        try:
            numbers = [int(p) if re.fullmatch(r"[+-]?\d+", p) else float(p) for p in parts]
        except ValueError:
            out[key] = value
            continue
        out[key] = tuple(float(x) for x in numbers) if len(numbers) > 1 else numbers[0]
    return out


def needs_cm_term(space: ValenceSpace) -> bool:
    """Whether the space can excite the centre of mass (KSHELL's helper then sets ``beta_cm = 10``).

    That happens when a positive- and a negative-parity orbit differ by at most one unit of j (an ħω excitation of
    the centre of mass connects them), as for sd + pf. Spaces with a core given as holes are left alone.
    """
    if space.core_z < 0 or space.core_n < 0:
        return False
    for tz in (-1, 1):
        orbits = space.species(tz)
        for a in orbits:
            for b in orbits:
                if a.parity == 1 and b.parity == -1 and abs(a.j2 - b.j2) <= 2:
                    return True
    return False


def default_hw_type(space: ValenceSpace) -> int:
    """KSHELL's ``hw_type`` as its helper chooses it: 2 for an ¹⁶O core, else 1.

    ``hw_type`` picks KSHELL's formula for the oscillator ħω (used for E2 operators): 1 is ħω = 41 A^(-1/3) MeV,
    2 is 45 A^(-1/3) − 25 A^(-2/3) MeV.
    """
    return 2 if (space.core_z, space.core_n) == (8, 8) else 1


# ---------------------------------------------------------------------------------------------------------------
# Interactions

#: The valence spaces of the interactions below, in the orbit order of their ``.snt`` files.
P_SPACE = _orbits(2, 2, "p", "0p1 0p3")
SD_SPACE = _orbits(8, 8, "sd", "0d3 0d5 1s1")
FP_SPACE = _orbits(20, 20, "fp", "0f7 1p3 0f5 1p1")
F5PG9_SPACE = _orbits(28, 28, "f5pg9", "1p3 0f5 1p1 0g9")
SDPF_SPACE = _orbits(8, 8, "sdpf", "0d5 0d3 1s1 0f7 0f5 1p3 1p1")

FREE_GL = (1.0, 0.0)
FREE_GS = (5.585, -3.826)
DEFAULT_EFF_CHARGE = (1.5, 0.5)


@dataclass(frozen=True)
class Interaction:
    """An effective interaction KSHELL can use, with physim's default inputs for it.

    ``eff_charge`` (e_p, e_n), ``gl`` and ``gs`` (proton, neutron) are the values used unless a request gives its
    own; ``defaults_from`` says where they come from (the ``.snt`` file's ``!namelist`` lines, the interaction's
    recommendation, or physim's generic defaults 1.5/0.5 and free g-factors). ``hw_type`` and ``beta_cm`` are
    the values KSHELL's own helper would choose for this space.
    """

    name: str
    snt: str
    space: ValenceSpace
    citation: str
    eff_charge: tuple[float, float] = DEFAULT_EFF_CHARGE
    gl: tuple[float, float] = FREE_GL
    gs: tuple[float, float] = FREE_GS
    hw_type: int = 1
    beta_cm: float = 0.0
    defaults_from: str = "physim's generic defaults: e_p 1.5, e_n 0.5, free-nucleon g-factors"
    notes: str = ""

    @property
    def label(self) -> str:
        """``"USDB (sd shell, 16O core)"``."""
        return f"{self.name.upper()} ({self.space.name} shell, {self.space.core_label} core)"


def _interaction(name: str, snt: str, space: ValenceSpace, citation: str, **kw: Any) -> Interaction:
    return Interaction(name, snt, space, citation, hw_type=default_hw_type(space),
                       beta_cm=10.0 if needs_cm_term(space) else 0.0, **kw)


_GENERIC = "physim's generic defaults (e_p 1.5, e_n 0.5, free-nucleon g-factors); the .snt file recommends none"

#: The interactions physim knows, by name. Citations are the ones in the ``.snt`` file headers (KSHELL 2019).
INTERACTIONS: dict[str, Interaction] = {i.name: i for i in (
    _interaction(
        "ckpot", "ckpot.snt", P_SPACE,
        "S. Cohen and D. Kurath, Nucl. Phys. 73, 1 (1965)",
        defaults_from=_GENERIC,
    ),
    _interaction(
        "usda", "usda.snt", SD_SPACE,
        "W. A. Richter, S. Mkhize and B. A. Brown, Phys. Rev. C 78, 064302 (2008)",
        defaults_from=_GENERIC,
        notes="Citation as in the usda.snt header (a paper on sd-shell observables with USDA and USDB). The fit "
              "itself is usually cited as B. A. Brown and W. A. Richter, Phys. Rev. C 74, 034315 (2006) "
              "(not checked here).",
    ),
    _interaction(
        "usdb", "usdb.snt", SD_SPACE,
        "W. A. Richter, S. Mkhize and B. A. Brown, Phys. Rev. C 78, 064302 (2008)",
        defaults_from=_GENERIC,
        notes="The usdb.snt header says 'USDA' (copied from usda.snt), but its numbers are USDB's "
              "(single-particle energies 2.1117, -3.9257, -3.2079 MeV). The citation is the header's (a paper on "
              "sd-shell observables with USDA and USDB). The fit itself is usually cited as B. A. Brown and "
              "W. A. Richter, Phys. Rev. C 74, 034315 (2006) (not checked here).",
    ),
    _interaction(
        "kb3g", "kb3g.snt", FP_SPACE,
        "A. Poves, J. Sanchez-Solano, E. Caurier and F. Nowacki, Nucl. Phys. A 694, 157 (2001)",
        defaults_from=_GENERIC,
    ),
    _interaction(
        "gxpf1a", "gxpf1a.snt", FP_SPACE,
        "M. Honma, T. Otsuka, B. A. Brown and T. Mizusaki, Eur. Phys. J. A 25, Suppl. 1, 499 (2005)",
        defaults_from="effective charges 1.5, 0.5 from the gxpf1a.snt file; free-nucleon g-factors (physim's "
                      "default)",
    ),
    _interaction(
        "jun45", "jun45.snt", F5PG9_SPACE,
        "M. Honma, T. Otsuka, T. Mizusaki and M. Hjorth-Jensen, Phys. Rev. C 80, 064323 (2009)",
        eff_charge=(1.5, 1.1), gs=(3.910, -2.678),
        defaults_from="effective charges 1.5, 1.1 and g-factors gl 1.0, 0.0, gs 3.910, -2.678 (0.7 times free) "
                      "from the jun45.snt file",
    ),
    _interaction(
        "sdpf-mu", "sdpf-mu.snt", SDPF_SPACE,
        "Y. Utsuno, T. Otsuka, B. A. Brown, M. Honma, T. Mizusaki and N. Shimizu, Phys. Rev. C 86, 051301(R) (2012)",
        eff_charge=(1.35, 0.35),
        defaults_from="effective charges 1.35, 0.35 and hw_type 2 from the sdpf-mu.snt file; free-nucleon "
                      "g-factors (physim's default)",
        notes="sd + pf across N or Z = 20: KSHELL's Lawson term (beta_cm = 10) removes centre-of-mass excitations.",
    ),
)}


def interaction(name_or_interaction: Union[str, Interaction]) -> Interaction:
    """The interaction by name (case does not matter), or the interaction itself."""
    if isinstance(name_or_interaction, Interaction):
        return name_or_interaction
    key = str(name_or_interaction).strip().lower()
    if key.endswith(".snt"):
        key = key[:-4]
    if key not in INTERACTIONS:
        raise ValueError(f"unknown interaction {name_or_interaction!r}; physim knows {', '.join(INTERACTIONS)}")
    return INTERACTIONS[key]


# ---------------------------------------------------------------------------------------------------------------
# Options for a nucleus


def parities(space: ValenceSpace, z_val: int, n_val: int) -> tuple[int, ...]:
    """The parities the valence nucleons can make, (+1,), (−1,) or (+1, −1)."""
    if not space.holds(z_val, n_val):
        return ()
    odd_counts = []
    for tz, nf in ((-1, z_val), (1, n_val)):
        cap_odd = sum(o.capacity for o in space.species(tz) if o.parity < 0)
        cap_even = space.capacity(tz) - cap_odd
        odd_counts.append(range(max(0, nf - cap_even), min(nf, cap_odd) + 1))
    found = {(-1) ** (kp + kn) for kp in odd_counts[0][:2] for kn in odd_counts[1][:2]}
    return tuple(p for p in (1, -1) if p in found)


@dataclass(frozen=True)
class Option:
    """An interaction that applies to a nucleus: valence protons and neutrons, reachable parities, and the M-scheme
    dimension at the lowest M (0 or 1/2) for each parity (``dimension_by_parity``)."""

    interaction: str
    z_val: int
    n_val: int
    parities: tuple[int, ...]
    dimension_by_parity: dict[int, int] = field(default_factory=dict, compare=False, hash=False)

    @property
    def spec(self) -> Interaction:
        return INTERACTIONS[self.interaction]


def nuclide_valence(nuclide: str, inter: Union[str, Interaction]) -> tuple[int, int]:
    """Valence protons and neutrons of ``nuclide`` in the interaction's space; ValueError if it does not fit."""
    z, a = parse_nuclide(nuclide)
    spec = interaction(inter)
    z_val, n_val = spec.space.valence(z, a)
    if not spec.space.holds(z_val, n_val):
        raise ValueError(
            f"{nuclide} does not fit in the {spec.space.name} space of {spec.name} ({spec.space.core_label} core, "
            f"room for {spec.space.capacity(-1)} protons and {spec.space.capacity(1)} neutrons): it would have "
            f"{z_val} valence protons and {n_val} valence neutrons")
    return z_val, n_val


def options(nuclide: str) -> list[Option]:
    """The interactions that apply to ``nuclide``: its valence protons and neutrons fit in the space, with at least
    one valence nucleon. In the order of :data:`INTERACTIONS`."""
    z, a = parse_nuclide(nuclide)
    out = []
    for spec in INTERACTIONS.values():
        z_val, n_val = spec.space.valence(z, a)
        if not spec.space.holds(z_val, n_val) or z_val + n_val == 0:
            continue
        pars = parities(spec.space, z_val, n_val)
        dims = {}
        for p in pars:
            d = dimensions(spec.space, z_val, n_val, p)
            dims[p] = d[min(d)] if d else 0
        out.append(Option(spec.name, z_val, n_val, pars, dims))
    return out


# ---------------------------------------------------------------------------------------------------------------
# Dimensions


def _orbit_counts(j2: int) -> dict[tuple[int, int], int]:
    """(particles, 2M) → number of Slater determinants in one orbit of 2j = ``j2``."""
    table: dict[tuple[int, int], int] = {(0, 0): 1}
    for m2 in range(-j2, j2 + 1, 2):
        new = dict(table)
        for (k, m), c in table.items():
            key = (k + 1, m + m2)
            new[key] = new.get(key, 0) + c
        table = new
    return table


def _species_counts(orbits: tuple[Orbit, ...], nf: int) -> dict[tuple[int, int], int]:
    """(2M, parity) → number of Slater determinants of ``nf`` like nucleons in ``orbits``."""
    table: dict[tuple[int, int, int], int] = {(0, 0, 1): 1}
    for orbit in orbits:
        counts = _orbit_counts(orbit.j2)
        new: dict[tuple[int, int, int], int] = defaultdict(int)
        for (k, m, p), c in table.items():
            for (dk, dm), dc in counts.items():
                if k + dk > nf:
                    continue
                new[(k + dk, m + dm, p * orbit.parity ** dk)] += c * dc
        table = new
    return {(m, p): c for (k, m, p), c in table.items() if k == nf}


def dimensions(space: ValenceSpace, z_val: int, n_val: int, parity: int) -> dict[int, int]:
    """2M → M-scheme dimension for 2M ≥ 0 (the dimension is the same for −M), no truncation.

    Empty when the space has no states of that parity (or cannot hold the nucleons).
    """
    if parity not in (1, -1):
        raise ValueError(f"parity is +1 or -1, got {parity!r}")
    if not space.holds(z_val, n_val):
        return {}
    protons = _species_counts(space.proton_orbits, z_val)
    neutrons = _species_counts(space.neutron_orbits, n_val)
    out: dict[int, int] = defaultdict(int)
    for (mp, pp), cp in protons.items():
        for (mn, pn), cn in neutrons.items():
            if pp * pn == parity and mp + mn >= 0:
                out[mp + mn] += cp * cn
    return {m: out[m] for m in sorted(out) if out[m]}


def j_dimension(space: ValenceSpace, z_val: int, n_val: int, parity: int, j2: int) -> int:
    """How many states of 2J = ``j2`` and the given parity: dim(M = J) − dim(M = J + 1)."""
    dims = dimensions(space, z_val, n_val, parity)
    return dims.get(j2, 0) - dims.get(j2 + 2, 0)


# ---------------------------------------------------------------------------------------------------------------
# Partition files


def _occupations(capacities: list[int], nf: int) -> Iterator[tuple[int, ...]]:
    """Every way to put ``nf`` nucleons in orbits of the given capacities."""
    if not capacities:
        if nf == 0:
            yield ()
        return
    rest = sum(capacities[1:])
    for k in range(max(0, nf - rest), min(capacities[0], nf) + 1):
        for tail in _occupations(capacities[1:], nf - k):
            yield (k,) + tail


def partition_lists(space: ValenceSpace, z_val: int, n_val: int, parity: int
                    ) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], list[tuple[int, int]]]:
    """The proton partitions, neutron partitions (orbit occupations) and allowed (proton, neutron) pairs (0-based)
    of a partition file with no truncation, in KSHELL's order. Partitions with no partner are left out."""
    if parity not in (1, -1):
        raise ValueError(f"parity is +1 or -1, got {parity!r}")
    if not space.holds(z_val, n_val):
        raise ValueError(f"the {space.name} space cannot hold {z_val} valence protons and {n_val} neutrons")
    lists, pars = [], []
    for tz, nf in ((-1, z_val), (1, n_val)):
        orbits = space.species(tz)
        occ = sorted(_occupations([o.capacity for o in orbits], nf))
        lists.append(occ)
        pars.append([math.prod(o.parity ** k for o, k in zip(orbits, arr)) for arr in occ])
    pairs = [(i, j) for i, pi in enumerate(pars[0]) for j, pj in enumerate(pars[1]) if pi * pj == parity]
    if not pairs:
        raise ValueError(f"no {'+' if parity > 0 else '-'} parity states for {z_val} protons and {n_val} neutrons "
                         f"in the {space.name} space")
    used_p = sorted({i for i, _ in pairs})
    used_n = sorted({j for _, j in pairs})
    new_p = {old: new for new, old in enumerate(used_p)}
    new_n = {old: new for new, old in enumerate(used_n)}
    return ([lists[0][i] for i in used_p], [lists[1][j] for j in used_n],
            [(new_p[i], new_n[j]) for i, j in pairs])


def partitions(space: ValenceSpace, z_val: int, n_val: int, parity: int, snt_name: Optional[str] = None) -> str:
    """The partition (``.ptn``) file KSHELL reads, with no truncation, laid out as KSHELL's ``gen_partition.py``
    writes it. ``snt_name`` goes in the comment on the first line (default ``<space name>.snt``).

    Lines end with ``\\n``; write the text in text mode (KSHELL reads either line ending).
    """
    prot, neut, pairs = partition_lists(space, z_val, n_val, parity)
    snt_name = snt_name or f"{space.name}.snt"
    out = [f"# partition file of {snt_name}  Z={z_val}  N={n_val}  parity={parity:+d}",
           f" {z_val} {n_val} {parity}",
           "# num. of  proton partition, neutron partition",
           f" {len(prot)} {len(neut)}"]
    for title, occ in (("# proton partition", prot), ("# neutron partition", neut)):
        out.append(title)
        for i, arr in enumerate(occ):
            out.append(f" {i + 1:5d}   " + "".join(f" {k:2d}" for k in arr))
    out.append("# partition of proton and neutron")
    out.append(f"{len(pairs)}")
    out.extend(f"{i + 1:5d} {j + 1:5d}" for i, j in pairs)
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------------------------------------------
# Memory and the estimate


def available_memory_bytes() -> int:
    """RAM free for a new program now, in bytes, from the operating system.

    Windows: ``GlobalMemoryStatusEx`` (``ullAvailPhys``); Linux: ``MemAvailable`` in ``/proc/meminfo``; elsewhere
    the available physical pages from ``os.sysconf``. Raises OSError when none of these can be read.
    """
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [("dwLength", wintypes.DWORD), ("dwMemoryLoad", wintypes.DWORD),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

        status = MEMORYSTATUSEX()
        status.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            raise OSError("GlobalMemoryStatusEx failed")
        return int(status.ullAvailPhys)
    try:
        with open("/proc/meminfo", encoding="ascii") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    try:
        return int(os.sysconf("SC_AVPHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, ValueError, OSError):
        raise OSError("cannot read the free memory on this system") from None


def lanczos_parameters(n_eigen: int) -> tuple[int, int]:
    """(``n_restart_vec``, ``max_lanc_vec``) for ``n_eigen`` states, as KSHELL's helper sets them: at least 10
    restart vectors and 1.5 per state; at least 200 Lanczos vectors and 50 more than the restart vectors."""
    n_restart = max(int(n_eigen * 1.5), N_RESTART_VEC)
    return n_restart, max(MAX_LANC_VEC, n_restart + 50)


def vectors_in_memory(n_eigen: int, mode_lv_hdd: int) -> int:
    """How many vectors of the full dimension KSHELL holds in memory for one J± run.

    With the Lanczos vectors on disk (``mode_lv_hdd = 1``) KSHELL's logs show n_eigen + 3 (6 for three states, 5
    for two); in memory (``0``) the ``max_lanc_vec`` Lanczos vectors come on top of those.
    """
    if mode_lv_hdd == 0:
        return lanczos_parameters(n_eigen)[1] + n_eigen + 3
    return n_eigen + 3


def _memory(dimension: int, n_vectors: int) -> int:
    return int(dimension * BYTES_PER_COMPONENT * n_vectors * MEMORY_FACTOR) + FIXED_OVERHEAD_BYTES


def _gb(n: float) -> str:
    gb = n / 1024**3
    return f"{gb:.2g} GB" if gb < 10 else f"{gb:.0f} GB"


def _state_label(j2: int, parity: int) -> str:
    j = f"{j2 // 2}" if j2 % 2 == 0 else f"{j2}/2"
    return j + ("+" if parity > 0 else "-")


@dataclass
class Estimate:
    """What a shell-model request needs here, before it runs.

    - ``dimension_by_state``: (2J, parity) → M-scheme dimension at M = J (KSHELL's J-projected runs work there);
      0 when no state of that 2J exists.
    - ``memory_bytes``: peak memory in the chosen mode (the largest J± run, or transit.exe if larger);
      ``disk_bytes``: wave-function files kept for the transitions, plus the Lanczos vectors when on disk.
    - ``mode_lv_hdd``: 0 (Lanczos vectors in memory, fastest) when that fits the limit, else 1 (on disk).
    - ``limit_bytes``: half of ``free_bytes``, the RAM free when the estimate was made.
    - ``fits`` and ``reason``: whether it can run here, in plain language.
    - ``seconds_guess``: a **rough** time, scaled from two wave 0 runs on a 6-core desktop (Ryzen 5 3600): expect
      it to be off by a factor of two either way, more for spaces much larger than ⁴⁸Cr in the fp shell.
    """

    dimension_by_state: dict[tuple[int, int], int]
    memory_bytes: int
    disk_bytes: int
    mode_lv_hdd: int
    limit_bytes: int
    fits: bool
    reason: str
    seconds_guess: float
    free_bytes: int = 0
    memory_bytes_in_memory_mode: int = 0
    memory_bytes_on_disk_mode: int = 0
    threads: int = 1
    interaction: str = ""
    z_val: int = 0
    n_val: int = 0

    def to_dict(self) -> dict[str, Any]:
        """JSON-safe: states as ``[{"j2", "parity", "label", "dimension"}]``."""
        return {
            "states": [{"j2": j2, "parity": p, "label": _state_label(j2, p), "dimension": d}
                       for (j2, p), d in self.dimension_by_state.items()],
            "memory_bytes": self.memory_bytes, "disk_bytes": self.disk_bytes, "mode_lv_hdd": self.mode_lv_hdd,
            "limit_bytes": self.limit_bytes, "free_bytes": self.free_bytes, "fits": self.fits,
            "reason": self.reason, "seconds_guess": self.seconds_guess, "seconds_guess_is_rough": True,
            "memory_bytes_in_memory_mode": self.memory_bytes_in_memory_mode,
            "memory_bytes_on_disk_mode": self.memory_bytes_on_disk_mode, "threads": self.threads,
            "interaction": self.interaction, "z_val": self.z_val, "n_val": self.n_val,
        }


def seconds_guess(dimension: int, n_eigen: int, threads: int) -> float:
    """A rough time for one kshell.exe run: linear in the dimension, √(states), threads^0.69 speed-up.

    Scaled from wave 0 (Ryzen 5 3600): ⁴⁸Cr KB3 ground state, dimension 1 963 461, 72 s on 12 threads. The same
    rule gives 7.9, 7.1 and 4.7 s for the ²⁸Si USDB 0⁺ (3 states), 2⁺ (3) and 4⁺ (2) runs, measured 7.6, 13.2 and
    4.8 s, and 405 s for 32 ⁴⁸Cr states, measured 402 s.
    """
    if dimension <= 0:
        return 0.0
    threads = max(1, int(threads))
    speed = (_TIME_REF_THREADS / threads) ** _THREAD_EXPONENT
    return _SECONDS_PER_RUN + _TIME_REF_SECONDS * dimension / _TIME_REF_DIMENSION * math.sqrt(max(1, n_eigen)) * speed


def estimate(request: Any, free_bytes: Optional[int] = None) -> Estimate:
    """Memory, disk, vector storage and a rough time for a shell-model request, and whether it fits here.

    ``request`` is anything with ``nuclide``, ``interaction`` (a name or :class:`Interaction`) and ``states``
    ((2J, parity, how many) tuples), such as :class:`physim.nuclear.shellmodel.ShellModelRequest`; an optional
    ``threads`` (None: all logical cores). The limit is half of ``free_bytes`` (default: the RAM free now). Each
    J± is a separate kshell.exe run, so the peak is the largest run; transitions between every pair of requested
    J± sets are counted for transit.exe.
    """
    spec = interaction(request.interaction)
    z_val, n_val = nuclide_valence(request.nuclide, spec)
    if free_bytes is None:
        free_bytes = available_memory_bytes()
    free_bytes = max(0, int(free_bytes))
    limit = free_bytes // 2
    threads = getattr(request, "threads", None) or os.cpu_count() or 1

    states: dict[tuple[int, int], int] = {}
    for item in request.states:
        j2, parity, n = (int(x) for x in item)
        if parity not in (1, -1):
            raise ValueError(f"parity is +1 or -1, got {parity!r} in {item!r}")
        if j2 < 0 or n < 1:
            raise ValueError(f"a state is (2J >= 0, parity, how many >= 1), got {item!r}")
        if (j2 + z_val + n_val) % 2:
            kind = "odd" if (z_val + n_val) % 2 else "even"
            raise ValueError(f"2J = {j2} is impossible for {request.nuclide}: A is {kind}, so 2J must be {kind}")
        states[(j2, parity)] = states.get((j2, parity), 0) + n
    if not states:
        raise ValueError("the request asks for no states")

    dims_by_parity = {p: dimensions(spec.space, z_val, n_val, p) for p in {p for _, p in states}}
    dimension_by_state = {(j2, p): dims_by_parity[p].get(j2, 0) for (j2, p) in states}

    # J± sets with fewer states than asked: KSHELL finds at most j_dimension of them.
    missing = []
    for (j2, p), n in states.items():
        jdim = j_dimension(spec.space, z_val, n_val, p, j2)
        if jdim < n:
            missing.append(f"{_state_label(j2, p)} ({jdim} exist{'s' if jdim == 1 else ''}, {n} asked)")

    def peak(mode: int) -> int:
        runs = [_memory(dimension_by_state[k], vectors_in_memory(n, mode)) for k, n in states.items()
                if dimension_by_state[k]]
        keys = [k for k in states if dimension_by_state[k]]
        # transit.exe holds both sets of wave functions and one work vector for each.
        transit = [FIXED_OVERHEAD_BYTES + int(((states[a] + 1) * dimension_by_state[a]
                                               + (states[b] + 1) * dimension_by_state[b])
                                              * BYTES_PER_COMPONENT * MEMORY_FACTOR)
                   for ia, a in enumerate(keys) for b in keys[ia:]]
        return max(runs + transit + [FIXED_OVERHEAD_BYTES])

    mem0, mem1 = peak(0), peak(1)
    wave_files = sum(n * dimension_by_state[k] * BYTES_PER_COMPONENT for k, n in states.items())
    largest = max(states, key=lambda k: dimension_by_state[k])
    lanczos_disk = lanczos_parameters(states[largest])[1] * dimension_by_state[largest] * BYTES_PER_COMPONENT

    if mem0 <= limit:
        mode, memory, disk, fits = 0, mem0, wave_files, True
    else:
        mode, memory, disk, fits = 1, mem1, wave_files + lanczos_disk, mem1 <= limit

    limit_text = f"the limit of {_gb(limit)} (half of the {_gb(free_bytes)} of RAM free now)"
    if fits and mode == 0:
        reason = f"Fits: about {_gb(memory)} with the Lanczos vectors in memory, under {limit_text}."
    elif fits:
        reason = (f"Fits with the Lanczos vectors on disk: about {_gb(memory)} of memory and {_gb(disk)} of disk, "
                  f"under {limit_text}. Keeping them in memory would need about {_gb(mem0)}.")
    else:
        reason = (f"Too big for this computer now: it needs about {_gb(mem1)} even with the Lanczos vectors on "
                  f"disk, over {limit_text}. Close other programs, ask for fewer or smaller J, or run it on a "
                  f"cluster.")
    big = f"{max(dimension_by_state.values()):,}".replace(",", " ")
    reason += f" Largest M-scheme dimension {big}."
    if missing:
        reason += " Fewer states exist than asked for " + ", ".join(missing) + "."

    seconds = sum(seconds_guess(dimension_by_state[k], n, threads) for k, n in states.items())
    n_sets = sum(1 for k in states if dimension_by_state[k])
    seconds += _SECONDS_PER_TRANSIT * n_sets * (n_sets + 1) / 2

    return Estimate(dimension_by_state, memory, disk, mode, limit, fits, reason, round(seconds, 1),
                    free_bytes=free_bytes, memory_bytes_in_memory_mode=mem0, memory_bytes_on_disk_mode=mem1,
                    threads=int(threads), interaction=spec.name, z_val=z_val, n_val=n_val)
