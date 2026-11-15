"""Neutron cross sections of carbon and hydrogen from ENDF/B-VIII.0, written to
tests/reference/nuclear/neutron_carbon_endf.csv (the check of physim's carbon table and of Gammel's n-p formula).

Usage: python -I scripts/endf_neutron_reference.py <folder> [--toml]

<folder> holds the ENDF-6 files n-006_C_012.endf, n-006_C_013.endf and n-001_H_001.endf of the ENDF/B-VIII.0
neutron sublibrary (NNDC, https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_neutrons.zip; the IAEA mirrors the
same files at https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.0/n/). The files are only read: the script
parses MF=3 (pointwise cross sections, linearly interpolable in ENDF/B-VIII.0 for these light nuclei) itself, with
no ENDF library. Natural carbon is 98.93% 12C and 1.07% 13C (IUPAC atom fractions).

The reference file lists, at each comparison energy, the pointwise value and the averages over windows of +-2.5 %
and +-10 % in energy (WINDOW_PERCENT, WIDE_PERCENT) for: the total cross section of natural carbon (the one physim uses, to attenuate the
neutrons in a cell and in cross-talk), its elastic and non-elastic parts and 12C(n,alpha) (context: the limits of
docs/theory/neutrons.md), and the n-p total cross section (H-1, MT=1).

With --toml, also prints the window-averaged natural-carbon total on physim's table grid, the [carbon_total] arrays
of python/physim/nuclear/data/neutron_scintillators.toml.
"""
import datetime
import hashlib
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "tests" / "reference" / "nuclear" / "neutron_carbon_endf.csv"
FILES = {"C12": "n-006_C_012.endf", "C13": "n-006_C_013.endf", "H1": "n-001_H_001.endf"}
#: IUPAC atom fractions of 12C and 13C in natural carbon.
ABUNDANCE = {"C12": 0.9893, "C13": 0.0107}
#: Half-widths of the averaging windows, percent of the energy: about one bin of physim's table, and a wide window
#: (about the energy resolution of time of flight over half a metre) in which the resonances average out.
WINDOW_PERCENT = 2.5
WIDE_PERCENT = 10.0
#: The comparison energies, MeV: 0.5-20 MeV, 12 points per decade, plus the neutron-standard energies.
COMPARE = sorted(set([round(0.5 * 10 ** (k / 12), 4) for k in range(20)] + [1.0, 2.0, 5.0, 10.0, 14.1, 20.0]))
#: physim's table grid, MeV: 40 points per decade from 10 keV to 30 MeV. Each entry is the average over its own
#: bin, from half a step below to half a step above (a factor 10^(1/80) = 1.029 each way), so the bins tile the
#: energy axis and a resonance narrower than a bin keeps its area.
PER_DECADE = 40
TABLE = [round(0.01 * 10 ** (k / PER_DECADE), 5) for k in range(int(PER_DECADE * math.log10(3000)) + 1)]


def _float(text: str) -> float:
    """An ENDF-6 number: Fortran style without the E ("1.234567+6", "-2.5-3") or plain."""
    s = text.strip()
    if not s:
        return 0.0
    for i in range(len(s) - 1, 0, -1):
        if s[i] in "+-" and s[i - 1] not in "eE":
            return float(s[:i] + "e" + s[i:])
    return float(s)


def read_mf3(path: Path) -> dict:
    """{MT: (energies in MeV, cross sections in b, interpolation laws)} for every MF=3 section of an ENDF-6 file."""
    sections: dict = {}
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if len(line) < 75:
            continue
        try:
            mf, mt = int(line[70:72]), int(line[72:75])
        except ValueError:
            continue
        if mf == 3 and mt > 0:
            sections.setdefault(mt, []).append(line[:66])
    out = {}
    for mt, lines in sections.items():
        # HEAD record, then the TAB1 record: C1 C2 L1 L2 NR NP, NR (NBT, INT) pairs, NP (E, sigma) pairs
        nr, np_ = int(_float(lines[1][44:55])), int(_float(lines[1][55:66]))
        fields = [_float(ln[i:i + 11]) for ln in lines[2:] for i in range(0, 66, 11)]
        interp = [(int(fields[2 * k]), int(fields[2 * k + 1])) for k in range(nr)]
        start = 6 * math.ceil(2 * nr / 6)
        pairs = fields[start:start + 2 * np_]
        out[mt] = ([pairs[2 * k] * 1e-6 for k in range(np_)], [pairs[2 * k + 1] for k in range(np_)], interp)
    return out


def _interp(xs: list, ys: list, x: float) -> float:
    """Linear interpolation (lin-lin, the law of these sections), 0 outside the table."""
    if x < xs[0] or x > xs[-1]:
        return 0.0
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if xs[mid] <= x else (lo, mid)
    if xs[hi] == xs[lo]:
        return ys[hi]
    return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo])


def _average(xs: list, ys: list, a: float, b: float) -> float:
    """Exact average of the piecewise-linear table over [a, b]."""
    pts = [a] + [x for x in xs if a < x < b] + [b]
    total = sum(0.5 * (_interp(xs, ys, p) + _interp(xs, ys, q)) * (q - p) for p, q in zip(pts[:-1], pts[1:]))
    return total / (b - a)


def check_linear(name: str, mf3: dict, mts, above_mev: float = 1e-3) -> None:
    """Stop unless every interpolation region reaching above ``above_mev`` is lin-lin (law 2), the only law the
    script interpolates (H-1 is log-log at thermal energies only)."""
    for mt in mts:
        xs, _, regions = mf3[mt]
        for nbt, law in regions:
            if law != 2 and xs[nbt - 1] > above_mev:
                raise SystemExit(f"{name} MT={mt}: interpolation law {law} up to {xs[nbt - 1]:g} MeV; this script "
                                 "handles lin-lin (2) only")


class Curve:
    """A cross section as a sum of weighted MF=3 sections: point values and window averages."""

    def __init__(self, parts):
        self.parts = parts  # [(weight, energies, sigmas)]

    def at(self, e: float) -> float:
        return sum(w * _interp(xs, ys, e) for w, xs, ys in self.parts)

    def mean(self, e: float, half_percent: float = WINDOW_PERCENT) -> float:
        """Average over E (1 +- half_percent / 100)."""
        return self.between(e * (1 - half_percent / 100), e * (1 + half_percent / 100))

    def between(self, a: float, b: float) -> float:
        """Average over [a, b], cut to the energies every part covers (H-1 stops at 20 MeV)."""
        top = min(xs[-1] for _, xs, _ in self.parts)
        a, b = min(a, top), min(b, top)
        if b <= a:
            return self.at(top)
        return sum(w * _average(xs, ys, a, b) for w, xs, ys in self.parts)


def curves(folder: Path) -> dict:
    data = {k: read_mf3(folder / f) for k, f in FILES.items()}
    check_linear("C12", data["C12"], (1, 2, 107))
    check_linear("C13", data["C13"], (1, 2))
    check_linear("H1", data["H1"], (1,))

    def part(iso, mt, w=1.0):
        xs, ys, _ = data[iso][mt]
        return (w, xs, ys)

    # The 13C evaluation stops at 20 MeV, the 12C one at 150 MeV: above 20 MeV, 13C takes 12C's values (1% of
    # natural carbon, and the two are within about 10% of each other there).
    c12, c13 = data["C12"], data["C13"]
    for mt in (1, 2):
        xs, ys, regions = c13[mt]
        more = [(x, y) for x, y in zip(c12[mt][0], c12[mt][1]) if x > xs[-1]]
        c13[mt] = (xs + [x for x, _ in more], ys + [y for _, y in more], regions)
    nat = [(iso, ABUNDANCE[iso]) for iso in ("C12", "C13")]
    return {
        "natC_total": Curve([part(i, 1, w) for i, w in nat]),
        "natC_elastic": Curve([part(i, 2, w) for i, w in nat]),
        "natC_nonelastic": Curve([part(i, 1, w) for i, w in nat] + [part(i, 2, -w) for i, w in nat]),
        "C12_total": Curve([part("C12", 1)]),
        "C12_n_alpha": Curve([part("C12", 107)]),
        "H1_total": Curve([part("H1", 1)]),
    }


def main(argv: list) -> None:
    folder = Path(argv[1])
    cs = curves(folder)
    sums = {f: hashlib.sha256((folder / f).read_bytes()).hexdigest()[:16] for f in FILES.values()}
    lines = [
        "# tool: ENDF/B-VIII.0 neutron sublibrary (NNDC, ENDF-B-VIII.0_neutrons.zip), MF=3 pointwise cross sections; "
        "files " + ", ".join(f"{f} (sha256 {s}...)" for f, s in sums.items()),
        f"# produced by: scripts/endf_neutron_reference.py, {datetime.date.today().isoformat()}",
        "# settings: 0 K pointwise data, lin-lin interpolation (the files' own law); natural carbon = 0.9893 12C + "
        "0.0107 13C (IUPAC); nonelastic = total - elastic; "
        f"point_b is the value at E, mean_b the exact average over E*(1 +- {WINDOW_PERCENT / 100:g}), wide_b over "
        f"E*(1 +- {WIDE_PERCENT / 100:g}) (windows cut at 20 MeV for H1); "
        "barns; H1_total is the n-p total (MT=1)",
        "energy_mev,quantity,point_b,mean_b,wide_b",
    ]
    for e in COMPARE:
        for q, c in cs.items():
            lines.append(f"{e:g},{q},{c.at(e):.5g},{c.mean(e):.5g},{c.mean(e, WIDE_PERCENT):.5g}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({len(COMPARE)} energies x {len(cs)} quantities)")
    if "--toml" in argv:
        tot = cs["natC_total"]
        print("energy_mev = [" + ", ".join(f"{e:g}" for e in TABLE) + "]")
        half = 10 ** (0.5 / PER_DECADE)
        print("sigma_b = [" + ", ".join(f"{tot.between(e / half, e * half):.4g}" for e in TABLE) + "]")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv)
