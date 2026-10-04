"""Write LISE++ energy-loss reference values for the physim tests, using LISE++'s Excel library (LISE_excel64.dll).

    python scripts/lise_reference.py "D:/Program Files/LISEcute/LISEcute"

Runs only on Windows with LISE++ installed. Writes tests/reference/nuclear/stopping_lise.csv; CI reads that file and
never needs LISE++. The library takes energies in MeV/u and thicknesses in mg/cm², and only elemental targets.

The library's energy-loss model must be selected before use: until set_loss is called, residual energies come out
about eight times too low (seen with LISEcute, October 2026). Every value here is computed after set_loss.
"""

from __future__ import annotations

import csv
import ctypes
import os
import sys
from datetime import date
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "tests" / "reference" / "nuclear" / "stopping_lise.csv"

#: LISE++ energy-loss options.
OPTIONS = {0: "Hubert", 1: "Ziegler", 4: "ATIMA 1.4"}
IONS = [(1, 1), (2, 4), (6, 12), (8, 16), (18, 40), (36, 84), (54, 132)]
#: Elemental targets: tabulated by NIST (C, Al, Si, Cu, Ag, Au) and not (Ni, Ta), to test the interpolation.
TARGETS = [6, 13, 14, 28, 29, 47, 73, 79]
ENERGIES = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]  # MeV/u


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
        "set_straggling": ([c.c_long], None),
    }
    for name, (args, res) in sig.items():
        f = getattr(lib, name)
        f.argtypes, f.restype = args, res
    return lib


def main(folder: str) -> None:
    lib = load(folder)
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


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
