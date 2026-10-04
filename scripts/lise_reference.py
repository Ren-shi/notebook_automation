"""Write LISE++ energy-loss reference values for the physim tests, using LISE++'s Excel library (LISE_excel64.dll).

    python scripts/lise_reference.py "D:/Program Files/LISEcute/LISEcute" [--masses-only]

Runs only on Windows with LISE++ installed. Writes tests/reference/nuclear/stopping_lise.csv and masses_lise.csv; CI
reads those files and never needs LISE++. The library takes energies in MeV/u and thicknesses in mg/cm², and only elemental targets.

The library's energy-loss model must be selected before use: until set_loss is called, residual energies come out
about eight times too low (seen with LISEcute, October 2026). Every value here is computed after set_loss.

straggling_energy also depends on a straggling setting the library keeps between sessions (it follows the GUI): a
second run on 2026-10-04 gave different straggling values, with everything else identical. The committed file is
the first run; regenerate only the masses (--masses-only) unless the straggling setting is known.
"""

from __future__ import annotations

import csv
import ctypes
import os
import sys
from datetime import date
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "tests" / "reference" / "nuclear" / "stopping_lise.csv"
MASSES_OUT = OUT.with_name("masses_lise.csv")

#: LISE++ energy-loss options.
OPTIONS = {0: "Hubert", 1: "Ziegler", 4: "ATIMA 1.4"}
IONS = [(1, 1), (2, 4), (6, 12), (8, 16), (18, 40), (36, 84), (54, 132)]
#: Elemental targets: tabulated by NIST (C, Al, Si, Cu, Ag, Au) and not (Ni, Ta), to test the interpolation.
TARGETS = [6, 13, 14, 28, 29, 47, 73, 79]
ENERGIES = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]  # MeV/u
#: Nuclides for the mass comparison: light, stable, beams and targets in common use, and neutron-rich ones (not the
#: neutron: isotope_mass(0, 1) returns 0).
NUCLIDES = [(1, 1), (1, 2), (1, 3), (2, 3), (2, 4), (2, 6), (3, 6), (3, 7), (4, 9), (5, 10), (5, 11), (6, 12),
            (6, 13), (6, 14), (7, 14), (8, 16), (8, 18), (9, 19), (10, 20), (11, 23), (12, 24), (13, 27), (14, 28),
            (16, 32), (18, 40), (20, 40), (20, 48), (26, 56), (28, 58), (28, 78), (29, 63), (36, 84), (40, 90),
            (47, 107), (50, 100), (50, 120), (50, 132), (54, 132), (56, 138), (62, 152), (73, 181), (79, 197),
            (82, 208), (83, 209), (90, 232), (92, 235), (92, 238)]


def load(folder: str):
    os.add_dll_directory(folder)
    os.chdir(folder)
    lib = ctypes.CDLL(os.path.join(folder, "LISE_excel64.dll"))
    c = ctypes
    sig = {
        "stopping_power": ([c.c_long, c.c_long, c.c_double, c.c_long, c.c_long], c.c_double),
        "trange_option": ([c.c_long, c.c_long, c.c_double, c.c_long, c.c_long], c.c_double),
        "treste_option": ([c.c_long, c.c_long, c.c_double, c.c_long, c.c_double, c.c_long], c.c_double),
        "straggling_energy": ([c.c_long, c.c_long, c.c_double, c.c_long, c.c_double], c.c_double),
        "set_loss": ([c.c_long], None),
        "isotope_mass": ([c.c_long, c.c_long], c.c_double),
        "set_straggling": ([c.c_long], None),
    }
    for name, (args, res) in sig.items():
        f = getattr(lib, name)
        f.argtypes, f.restype = args, res
    return lib


def main(folder: str, masses_only: bool = False) -> None:
    lib = load(folder)
    if not masses_only:
        write_stopping(lib)
    write_masses(lib)


def write_stopping(lib) -> None:
    rows = []
    for opt, label in OPTIONS.items():
        lib.set_loss(opt)
        for z, a in IONS:
            for zt in TARGETS:
                for e in ENERGIES:
                    s = lib.stopping_power(z, a, e, zt, opt)
                    r = lib.trange_option(z, a, e, zt, opt)
                    thick = 0.2 * r  # a layer a fifth of the range: a real slowing-down, well short of stopping
                    after = lib.treste_option(z, a, e, zt, thick, opt)
                    strag = lib.straggling_energy(z, a, e, zt, thick)
                    rows.append({"option": opt, "model": label, "Z": z, "A": a, "Zt": zt, "e_mev_u": e,
                                 "stopping_mev_mg_cm2": f"{s:.6g}", "range_mg_cm2": f"{r:.6g}",
                                 "thickness_mg_cm2": f"{thick:.6g}", "e_after_mev_u": f"{after:.6g}",
                                 "straggling_mev_u": f"{strag:.6g}"})
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        f.write("# tool: LISE++ (LISEcute, Qt build), LISE_excel64.dll\n")
        f.write(f"# produced by: scripts/lise_reference.py, {date.today().isoformat()}\n")
        f.write("# settings: set_loss(option) before each block; stopping_power(Z,A,E,Zt,option) [MeV/(mg/cm2)],\n")
        f.write("#   trange_option [mg/cm2], treste_option over thickness_mg_cm2 [MeV/u], straggling_energy [MeV/u].\n")
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} rows)")


def write_masses(lib) -> None:
    with open(MASSES_OUT, "w", newline="", encoding="utf-8") as f:
        f.write("# tool: LISE++ (LISEcute, Qt build), LISE_excel64.dll\n")
        f.write(f"# produced by: scripts/lise_reference.py, {date.today().isoformat()}\n")
        f.write("# settings: isotope_mass(Z, A), the nuclear (bare) mass in u, from LISE++'s default mass table.\n")
        f.write("Z,A,nuclear_mass_u\n")
        for z, a in NUCLIDES:
            f.write(f"{z},{a},{lib.isotope_mass(z, a):.12f}\n")
    print(f"wrote {MASSES_OUT} ({len(NUCLIDES)} rows)")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--masses-only"]
    if len(args) != 1:
        raise SystemExit(__doc__)
    main(args[0], masses_only="--masses-only" in sys.argv)
