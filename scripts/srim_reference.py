"""Write SRIM/TRIM reference values for the physim tests (the cases of the former pending/stopping_srim.md).

    python scripts/srim_reference.py "D:/Projects/srim/SRIM2013" [--sr-only | --trim-only]

Runs only on Windows with SRIM 2013 installed (and its VB5 runtime registered). The installed folder is never
written to: ``SR Module`` and the files TRIM needs are copied to a temporary folder and run there.

- **SR Module** (``SRModule.exe``, no window): one run per ion, target and energy, with E-min = E-max, which gives a
  single table row with stopping, projected range and straggling. Writes ``stopping_srim.csv``.
- **TRIM** (``TRIM.exe`` with ``TRIMAUTO`` = 1, which needs no keyboard or mouse input and exits when all ions are
  done): each ``tests/reference/nuclear/trim/TRIM_n.IN`` is copied byte for byte to ``TRIM.IN``, so a file that runs here opened
  without changes. Writes ``straggling_trim.csv`` from ``SRIM Outputs/TRANSMIT.txt``.

CI reads the CSV files and never needs SRIM.
"""

from __future__ import annotations

import csv
import re
import statistics
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from physim.external.srim import Srim  # noqa: E402

REF = ROOT / "tests" / "reference" / "nuclear"

# --- SR Module ----------------------------------------------------------------------------------------------------

#: Ion, Z, mass entered in SRIM (u, the isotope's atomic mass), mass number (energies are E/u × A).
IONS = [("4He", 2, 4.0026, 4), ("12C", 6, 12.000, 12), ("16O", 8, 15.995, 16), ("40Ar", 18, 39.962, 40),
        ("132Xe", 54, 131.904, 132)]
#: Target, density (g/cm³), elements (Z, name, stoichiometry, mass in u: SRIM's natural masses; D for CD₂).
TARGETS = [
    ("Au", 19.32, [(79, "Gold", 1, 196.967)]),
    ("Si", 2.33, [(14, "Silicon", 1, 28.086)]),
    ("Mylar", 1.38, [(1, "Hydrogen", 8, 1.008), (6, "Carbon", 10, 12.011), (8, "Oxygen", 4, 15.999)]),
    ("CD2", 1.06, [(6, "Carbon", 1, 12.011), (1, "Deuterium", 2, 2.014)]),
    ("CsI", 4.51, [(55, "Cesium", 1, 132.905), (53, "Iodine", 1, 126.904)]),
]
ENERGIES = [0.1, 0.3, 1.0, 3.0, 10.0]  # MeV/u
_LEN = {"A": 1e-4, "um": 1.0, "mm": 1e3, "m": 1e6, "km": 1e9}
_EN = {"eV": 1e-3, "keV": 1.0, "MeV": 1e3, "GeV": 1e6}
_ROW = re.compile(r"^\s*([\d.]+)\s+(eV|keV|MeV|GeV)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.]+)\s+(A|um|mm|m|km)\s+"
                  r"([\d.]+)\s+(A|um|mm|m|km)\s+([\d.]+)\s+(A|um|mm|m|km)")


def sr_in(z: int, mass: float, density: float, elements: list, e_min: float, e_max: float) -> str:
    """An SR.IN control file; compound correction 1 (plain Bragg additivity), stopping in MeV/(mg/cm²)."""
    lines = ["---Stopping/Range Input Data (Number-format: Period = Decimal Point)", "---Output File Name",
             '"out.txt"', "---Ion(Z), Ion Mass(u)", f"{z}   {mass}",
             "---Target Data: (Solid=0,Gas=1), Density(g/cm3), Compound Corr.", f"0    {density}    1",
             "---Number of Target Elements", f" {len(elements)} ",
             "---Target Elements: (Z), Target name, Stoich, Target Mass(u)"]
    lines += [f'{zt}   "{name}"   {st}   {mt}' for zt, name, st, mt in elements]
    lines += ["---Output Stopping Units (1-8)", " 5", "---Ion Energy : E-Min(keV), E-Max(keV)", f" {e_min:g}    {e_max:g}"]
    return "\r\n".join(lines) + "\r\n"


def run_sr(srim: Srim, text: str) -> tuple:
    """Run SR Module on one SR.IN; (output text, [(E keV, Se, Sn, projected range µm, long. straggling µm)])."""
    body = srim.stopping(text, timeout_s=60)
    rows = []
    for line in body.splitlines():
        m = _ROW.match(line)
        if m:
            g = m.groups()
            rows.append((float(g[0]) * _EN[g[1]], float(g[2]), float(g[3]), float(g[4]) * _LEN[g[5]],
                         float(g[6]) * _LEN[g[7]]))
    return body, rows


def stopping_tables(srim: Path) -> None:
    with Srim(srim) as work:
        version = None
        results = []
        for ion, z, mass, a in IONS:
            for target, density, elements in TARGETS:
                for e_u in ENERGIES:
                    e_kev = e_u * a * 1000
                    body, rows = run_sr(work, sr_in(z, mass, density, elements, e_kev, e_kev))
                    if len(rows) != 1 or abs(rows[0][0] / e_kev - 1) > 1e-3:
                        raise RuntimeError(f"unexpected SR Module output for {ion} in {target}:\n{body}")
                    version = version or " / ".join(ln.strip() for ln in body.splitlines()[:4] if "SRIM" in ln)
                    results.append((ion, target, density, e_u, *rows[0][1:]))
                # A full table over the same span must give the same values at its grid points.
                _, table = run_sr(work, sr_in(z, mass, density, elements, 10, ENERGIES[-1] * a * 1000))
                for r in results[-len(ENERGIES):]:
                    for t in table:
                        if abs(t[0] / (r[3] * a * 1000) - 1) < 1e-6 and abs(t[3] / r[6] - 1) > 2e-3:
                            raise RuntimeError(f"table and single-energy runs disagree: {t} {r}")
    path = REF / "stopping_srim.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(f"# tool: SRIM 2013 (installer SRIM-2013-Std.exe), SR Module (SRModule.exe); its output header reads "
                f"\"{version}\"\n")
        f.write(f"# produced by: scripts/srim_reference.py, {date.today().isoformat()}\n")
        f.write("# settings: one SR.IN per row with E-min = E-max = E/u × A keV; stopping units 5 (MeV/(mg/cm2));\n"
                "#   solid targets; compound correction 1 for every target (plain Bragg additivity: SR Module takes the\n"
                "#   factor as input, and SRIM's Core-and-Bond factor for Mylar is computed only by the GUI);\n"
                "#   ion masses 4.0026, 12.000, 15.995, 39.962, 131.904 u; Mylar H8 C10 O4; CD2 C1 + H(Z=1) 2 at\n"
                "#   2.014 u; CsI Cs1 I1; ranges are SRIM's projected ranges, straggling longitudinal, in um.\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["ion", "target", "density_g_cm3", "e_mev_u", "electronic_mev_mg_cm2", "nuclear_mev_mg_cm2",
                    "projected_range_um", "longitudinal_straggling_um"])
        for ion, target, density, e_u, se, sn, rng, strag in results:
            w.writerow([ion, target, density, e_u, f"{se:.4g}", f"{sn:.4g}", f"{rng:.5g}", f"{strag:.4g}"])
    print(f"wrote {path} ({len(results)} rows)")


# --- TRIM ---------------------------------------------------------------------------------------------------------

#: The request's runs: file, then the columns of straggling_trim.csv as physim reads them.
TRIM_RUNS = [
    ("TRIM_1.IN", "4He", "Mylar", "", 5.5, "10 um"),
    ("TRIM_2.IN", "16O", "Au", "", 64.0, "0.5 um"),
    ("TRIM_3.IN", "1H", "Si", "", 10.0, "300 um"),
    ("TRIM_4.IN", "132Xe", "C", "2.0", 660.0, "2 um"),
]
#: Multi-layer exports, run only to see that they open.
TRIM_OPEN_ONLY = ["TRIM_alpha_on_gold.IN", "TRIM_oxygen_on_lead_array.IN"]


def run_trim(srim: Srim, source: Path) -> list:
    """Run TRIM unattended, out of sight, on an unchanged copy of ``source``; transmitted energies (MeV), or []."""
    out = srim.trim(source.read_bytes().decode("latin-1"), timeout_s=3 * 3600).get("TRANSMIT.txt", "")
    return [float(ln.split()[3]) / 1e6 for ln in out.splitlines() if ln.startswith("T ")]


def trim_runs(srim: Path) -> None:
    pending = REF / "trim"
    with Srim(srim) as work:
        rows = []
        for name, ion, material, density, energy, thickness in TRIM_RUNS:
            e = run_trim(work, pending / name)
            ions = int((pending / name).read_text(encoding="latin-1").splitlines()[2].split()[4])
            ok = len(e) > 0.9 * ions  # it ran to the end on the file as written (a few ions may stop or scatter)
            print(f"{name}: {len(e)} of {ions} ions transmitted")
            rows.append([ion, material, density, energy, thickness, f"{statistics.fmean(e):.5g}" if e else "",
                         f"{statistics.stdev(e):.4g}" if len(e) > 1 else "", "yes" if ok else "no"])
        for name in TRIM_OPEN_ONLY:
            e = run_trim(work, pending / name)
            print(f"{name}: {'ran' if e else 'NO TRANSMITTED IONS'}, {len(e)} transmitted")
    path = REF / "straggling_trim.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write("# tool: SRIM 2013 (installer SRIM-2013-Std.exe), TRIM.exe; output headers read \"SRIM-2013.00\"\n")
        f.write(f"# produced by: scripts/srim_reference.py, {date.today().isoformat()}\n")
        f.write("# settings: the request's TRIM_n.IN files exactly as written by physim.nuclear.export.trim_in (2000\n"
                "#   ions, quick Kinchin-Pease, compound correction 1, transmitted ions saved), run with TRIMAUTO = 1;\n"
                "#   mean and sample standard deviation of the energies of the transmitted ions (TRANSMIT.txt); density in g/cm3;\n"
                "#   opened_unchanged = yes if TRIM ran the unedited file to the end (TRIM_3: 1998 of 2000 protons\n"
                "#   transmitted in the first run, none backscattered). Random seed 0 as exported, so a rerun differs\n"
                "#   by about the standard error (sigma / sqrt(2000)).\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["ion", "material", "density", "energy_mev", "thickness", "mean_mev", "sigma_mev", "opened_unchanged"])
        w.writerows(rows)
    print(f"wrote {path}")


if __name__ == "__main__":
    srim = Path(sys.argv[1])
    if "--trim-only" not in sys.argv:
        stopping_tables(srim)
    if "--sr-only" not in sys.argv:
        trim_runs(srim)
