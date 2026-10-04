"""Convert downloaded NIST tables into the CSV files shipped in python/physim/nuclear/data/.

The AME2020 mass table (mass_1.mas20.txt) is shipped unchanged and needs no conversion. The NIST pages are
HTML or a text listing, so they are reduced here to plain CSV; the values are copied as published, only the layout
changes. Usage:

    python scripts/convert_nuclear_data.py <download dir>

where the download directory holds:
    nist_isotopic_compositions.txt   https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all
    nist_xraymasscoef_tab1.html      https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
    nist_xraymasscoef_tab2.html      https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html
"""

from __future__ import annotations

import csv
import html
import re
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "python" / "physim" / "nuclear" / "data"


def _cells(row_html: str) -> list[str]:
    cells = re.findall(r"<TD[^>]*>(.*?)</TD>", row_html, re.S | re.I)
    return [html.unescape(re.sub(r"<[^>]+>", " ", c)).strip() for c in cells]


def _rows(path: Path) -> list[list[str]]:
    text = path.read_text(encoding="latin-1")
    return [_cells(r) for r in re.findall(r"<TR[^>]*>(.*?)</TR>", text, re.S | re.I)]


def isotopes(src: Path) -> tuple[list[dict], dict[int, str]]:
    """Natural isotopes (those with an isotopic composition) and the standard atomic weight of each element."""
    text = src.read_text(encoding="latin-1")
    blocks = re.findall(r"Atomic Number = (\d+)\s*\nAtomic Symbol = (\w+)\s*\nMass Number = (\d+)\s*\n"
                        r"Relative Atomic Mass = ([^\n]*)\nIsotopic Composition = ([^\n]*)\n"
                        r"Standard Atomic Weight = ([^\n]*)\n", text)
    natural, weights = [], {}
    for z, symbol, a, mass, composition, weight in blocks:
        weights.setdefault(int(z), weight.strip())
        if composition.strip():
            natural.append({"Z": int(z), "A": int(a), "relative_atomic_mass": mass.strip(),
                            "isotopic_composition": composition.strip()})
    if len(weights) < 100:
        raise SystemExit(f"{src}: only {len(weights)} elements found; has the format changed?")
    return natural, weights


def elements(src: Path, weights: dict[int, str]) -> list[dict]:
    out = []
    for cells in _rows(src):
        cells = [c for c in cells if c]  # the first row also holds empty spacer cells spanning the table
        if len(cells) == 6 and cells[0].isdigit():
            z, symbol, name, z_over_a, i_ev, density = cells
            out.append({"Z": int(z), "symbol": symbol, "name": name, "standard_atomic_weight": weights[int(z)],
                        "Z_over_A": z_over_a, "I_eV": i_ev, "density_g_cm3": density})
    if len(out) != 92:
        raise SystemExit(f"{src}: expected 92 elements, found {len(out)}")
    return out


def compounds(src: Path) -> list[dict]:
    out = []
    for cells in _rows(src):
        cells = [c for c in cells if c]  # spacer cells, as in table 1
        if len(cells) == 5 and re.match(r"^[\d.]+$", cells[1]):
            name, z_over_a, i_ev, density, composition = cells
            parts = re.findall(r"(\d+)\s*:\s*([\d.]+)", composition)
            out.append({"name": re.sub(r"\s+", " ", name), "Z_over_A": z_over_a, "I_eV": i_ev,
                        "density_g_cm3": density, "mass_fractions": " ".join(f"{z}:{f}" for z, f in parts)})
    if len(out) != 48:
        raise SystemExit(f"{src}: expected 48 materials, found {len(out)}; has the format changed?")
    return out


def write(name: str, rows: list[dict]) -> None:
    with open(OUT / name, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT / name} ({len(rows)} rows)")


def main(download_dir: str) -> None:
    d = Path(download_dir)
    natural, weights = isotopes(d / "nist_isotopic_compositions.txt")
    write("nist_isotopes.csv", natural)
    write("nist_elements.csv", elements(d / "nist_xraymasscoef_tab1.html", weights))
    write("nist_compounds.csv", compounds(d / "nist_xraymasscoef_tab2.html"))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
