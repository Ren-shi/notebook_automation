"""Write TRIM references for the recoil transport (``physim.nuclear.transport``) in
``tests/reference/nuclear/lifetimes/trim/``.

    python scripts/transport_reference.py "D:/Projects/srim/SRIM2013" [--ions 5000] [--write-only]

Runs only on Windows with SRIM 2013 installed (the folder may also come from ``physim.external.programs.find("srim")``
or ``PHYSIM_SRIM``). The installed folder is never written to (``physim.external.srim.Srim`` copies what TRIM needs to
a temporary folder). ``--write-only`` writes the TRIM.IN files and runs nothing.

Two recoils of the lifetime experiments (``lifetimes/published/setups.toml``), at the speeds the papers give:

- ⁹⁸Ru at β = 4.5% (its speed after the 0.7 mg/cm² Mg target, Radeck et al. 2012) into Cu, the plunger's stopper;
- ¹⁴³Sm at β = 1.6% (Rajbanshi et al. 2018) into Au, the DSAM backing.

For each, TRIM is run three times (TRIM.IN files written by ``physim.nuclear.export.trim_in``, then their disk-file
line changed): a layer thicker than the range with the stopped ions' positions saved (``RANGE_3D.txt``: projected
range, its straggling, the lateral spread), and two layers of about a quarter and a half of the range with the
transmitted ions saved (``TRANSMIT.txt``: their energies and directions, the angular spread on the way). The
summaries go to ``transport_trim.csv``.

SRIM's own stopping powers differ from physim's by several per cent for these slow heavy ions, which moves the range
by as much. So SR Module also writes SRIM's electronic and nuclear stopping for both recoils (``srim_stopping.csv``,
``--stopping-only`` runs just that): with them in place of physim's, physim's Monte Carlo and TRIM share the stopping
and differ only in how they scatter, which is what the transport test compares. CI reads only the two CSV files.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from datetime import date
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from physim.nuclear import data  # noqa: E402
from physim.nuclear.export import trim_in  # noqa: E402

OUT = ROOT / "tests" / "reference" / "nuclear" / "lifetimes" / "trim"

#: (case, recoil, β, material, [(run, thickness mg/cm², what TRIM saves)]).
CASES = [
    ("ru98_cu", "98Ru", 0.045, "Cu", [("range", 15.0, "ranges"), ("quarter", 1.5, "transmit"),
                                      ("half", 3.0, "transmit")]),
    ("sm143_au", "143Sm", 0.016, "Au", [("range", 8.0, "ranges"), ("quarter", 0.75, "transmit"),
                                        ("half", 1.5, "transmit")]),
]
SEED = 1
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")


def energy_mev(ion: str, beta: float) -> float:
    """Kinetic energy at speed β (the nuclear mass, as physim's transport uses)."""
    m = data.nuclide(ion).nuclear_mass_mev
    return m * (1 / math.sqrt(1 - beta * beta) - 1)


def trim_file(ion: str, energy: float, material: str, thickness: float, save: str, ions: int) -> str:
    text = trim_in(ion, energy, [(material, f"{thickness} mg/cm2")], ions=ions, seed=SEED,
                   title=f"{ion} beta into {material} {thickness} mg/cm2 ({save})")
    lines = text.split("\r\n")
    k = next(i for i, ln in enumerate(lines) if ln.startswith("Diskfiles")) + 1
    ranges, transmit = (1, 0) if save == "ranges" else (0, 1)
    lines[k] = f"{ranges:27d} {0:7d} {transmit:11d} {0:7d} {0:15d} {0:31d}"
    return "\r\n".join(lines)


def files(ions: int) -> list:
    out = []
    for case, ion, beta, material, runs in CASES:
        e = energy_mev(ion, beta)
        for run, thickness, save in runs:
            name = f"TRIM_{case}_{run}.IN"
            out.append((case, ion, beta, e, material, run, thickness, save, name,
                        trim_file(ion, e, material, thickness, save, ions)))
    return out


def _rows(text: str, first: str) -> list:
    """Numeric rows of a TRIM output file: lines starting with ``first`` (``"T"`` in TRANSMIT.txt) or, for an empty
    ``first``, lines of numbers only (RANGE_3D.txt)."""
    rows = []
    for ln in text.splitlines():
        s = ln.strip()
        if first and not s.startswith(first + " "):
            continue
        body = s[len(first):] if first else s
        if not first and (not body or not body[0].isdigit()):
            continue
        nums = [float(x) for x in _NUMBER.findall(body)]
        if nums:
            rows.append(nums)
    return rows


def range_stats(text: str) -> dict:
    """From RANGE_3D.txt (ion number, depth X, lateral Y, lateral Z in Å): µm."""
    a = np.array([r[:4] for r in _rows(text, "") if len(r) >= 4])
    x, y, z = a[:, 1] * 1e-4, a[:, 2] * 1e-4, a[:, 3] * 1e-4
    lateral = np.sqrt(0.5 * (y**2 + z**2).mean())
    return {"ions": len(a), "projected_range_um": x.mean(), "range_straggling_um": x.std(ddof=1),
            "lateral_sigma_um": lateral, "radial_mean_um": np.sqrt(y**2 + z**2).mean()}


def transmit_stats(text: str, ions: int) -> dict:
    """From TRANSMIT.txt (ion, atom, energy eV, depth, lateral Y, Z in Å, cos X, cos Y, cos Z)."""
    a = np.array([r[:9] for r in _rows(text, "T") if len(r) >= 9])
    e = a[:, 2] * 1e-6
    cos_x = np.clip(a[:, 6], -1, 1)
    angle = np.degrees(np.arccos(cos_x))
    return {"ions": ions, "transmitted": len(a) / ions, "mean_energy_mev": e.mean(), "sigma_energy_mev": e.std(ddof=1),
            "rms_angle_deg": math.sqrt((angle**2).mean()), "median_angle_deg": float(np.median(angle)),
            "mean_cos": cos_x.mean(), "angle_90_deg": float(np.percentile(angle, 90)),
            "lateral_sigma_um": math.sqrt(0.5 * ((a[:, 4] ** 2 + a[:, 5] ** 2).mean())) * 1e-4}


#: SR Module: (recoil, Z, mass in u, material, target elements (Z, name, stoichiometry, mass in u), density).
SR_RUNS = [("98Ru", 44, 97.9053, "Cu", [(29, "Copper", 1, 63.546)], 8.96),
           ("143Sm", 62, 142.9146, "Au", [(79, "Gold", 1, 196.967)], 19.32)]


def srim_stopping(folder) -> None:
    """SRIM's stopping tables for the two recoils, 10 keV to above their starting energies (SR Module)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from srim_reference import run_sr, sr_in  # noqa: E402

    from physim.external.srim import Srim

    rows = []
    version = ""
    with Srim(folder) as srim:
        for ion, z, mass, material, elements, density in SR_RUNS:
            e_top = 1.5 * energy_mev(ion, dict((c[1], c[2]) for c in CASES)[ion]) * 1e3
            body, table = run_sr(srim, sr_in(z, mass, density, elements, 10, e_top))
            version = version or " / ".join(ln.strip() for ln in body.splitlines()[:4] if "SRIM" in ln)
            if len(table) < 20:
                raise RuntimeError(f"unexpected SR Module output for {ion} in {material}:\n{body}")
            for e_kev, se, sn, rng, strag in table:
                rows.append([ion, material, density, f"{e_kev * 1e-3:.6g}", f"{se:.5g}", f"{sn:.5g}", f"{rng:.5g}"])
            print(f"{ion} in {material}: {len(table)} rows")
    path = OUT / "srim_stopping.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(f"# tool: SRIM 2013 (installer SRIM-2013-Std.exe), SR Module (SRModule.exe); its output header reads "
                f"\"{version}\"\n")
        f.write(f"# produced by: scripts/transport_reference.py --stopping-only, {date.today().isoformat()}\n")
        f.write("# settings: one SR.IN per recoil, E-min 10 keV to E-max 1.5 x the starting energy; stopping units 5\n"
                "#   (MeV/(mg/cm2)); solid targets, compound correction 1; ion masses 97.9053 and 142.9146 u; the\n"
                "#   projected range in um. The stopping TRIM used for transport_trim.csv.\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["ion", "material", "density_g_cm3", "energy_mev", "electronic_mev_mg_cm2", "nuclear_mev_mg_cm2",
                    "projected_range_um"])
        w.writerows(rows)
    print(f"wrote {path}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("srim", nargs="?")
    ap.add_argument("--ions", type=int, default=5000)
    ap.add_argument("--write-only", action="store_true")
    ap.add_argument("--stopping-only", action="store_true")
    args = ap.parse_args()
    if args.stopping_only:
        from physim.external import programs

        srim_stopping(args.srim or os.environ.get("PHYSIM_SRIM") or programs.find("srim"))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    todo = files(args.ions)
    for *_, name, text in todo:
        (OUT / name).write_bytes(text.encode("latin-1"))
    print(f"wrote {len(todo)} TRIM.IN files in {OUT}")
    if args.write_only:
        return
    from physim.external import programs
    from physim.external.srim import Srim

    folder = args.srim or os.environ.get("PHYSIM_SRIM") or programs.find("srim")
    if not folder:
        raise SystemExit("give the SRIM 2013 folder")
    results = []
    with Srim(folder) as srim:
        for case, ion, beta, e, material, run, thickness, save, name, text in todo:
            print(f"running {name} ...", flush=True)
            out = srim.trim((OUT / name).read_bytes().decode("latin-1"))
            if save == "ranges":
                body = out.get("RANGE_3D.txt", "")
                if not body:
                    raise RuntimeError(f"{name}: TRIM wrote no RANGE_3D.txt ({sorted(out)})")
                stats = range_stats(body)
            else:
                body = out.get("TRANSMIT.txt", "")
                if not body:
                    raise RuntimeError(f"{name}: TRIM wrote no TRANSMIT.txt ({sorted(out)})")
                stats = transmit_stats(body, args.ions)
            print(f"  {stats}", flush=True)
            for key, value in stats.items():
                results.append([case, ion, beta, f"{e:.6g}", material, data.material(material).density_g_cm3, run,
                                thickness, key, f"{value:.6g}"])
    path = OUT / "transport_trim.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write("# tool: SRIM 2013 (installer SRIM-2013-Std.exe), TRIM.exe; output headers read \"SRIM-2013.00\"\n")
        f.write(f"# produced by: scripts/transport_reference.py, {date.today().isoformat()}\n")
        f.write(f"# settings: the TRIM_<case>_<run>.IN files beside this one ({args.ions} ions each, quick Kinchin-Pease,\n"
                f"#   compound correction 1, random seed {SEED}), written by physim.nuclear.export.trim_in with the disk-file\n"
                "#   line changed: 'range' runs save RANGE_3D.txt (stopped ions' depth and lateral position), 'quarter'\n"
                "#   and 'half' runs save TRANSMIT.txt (energies and direction cosines of transmitted ions). Energies\n"
                "#   from beta with the nuclear mass; ion masses in TRIM are atomic masses (u). Lengths in um, angles\n"
                "#   in degrees from the layer normal; lateral_sigma_um is the per-axis sigma, sqrt(<y^2 + z^2>/2).\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["case", "ion", "beta", "energy_mev", "material", "density_g_cm3", "run", "thickness_mg_cm2",
                    "quantity", "value"])
        w.writerows(results)
    print(f"wrote {path}")
    srim_stopping(folder)


if __name__ == "__main__":
    main()
