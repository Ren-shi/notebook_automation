"""Measured electronic stopping of heavy ions in hydrogen-rich targets, from the IAEA Stopping Power Database,
written to tests/reference/nuclear/stopping_measured_hydrogen.csv with SRIM 2013's value for each point.

Usage: python scripts/iaea_stopping_reference.py "D:/Projects/srim/SRIM2013"

Keeps measurements only (the database also serves SRIM 2013 curve values under the reference id "SRIM2013"; they
are dropped), ions C, O, Ar, Kr, Xe, 0.05-15 MeV/u, MeV/(mg/cm2). SRIM values come from SR Module run in a copy of
the install (see srim_reference.py), compound correction 1, gas flag set for H2 and D2.
"""
import csv
import datetime
import json
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import srim_reference as sr  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "tests" / "reference" / "nuclear" / "stopping_measured_hydrogen.csv"
API = "https://www-nds.iaea.org/stopping/api/data/ion:{ion}/target:{target}"
A = {"C": 12, "O": 16, "Ar": 40, "Kr": 84, "Xe": 132}
Z = {"C": 6, "O": 8, "Ar": 18, "Kr": 36, "Xe": 54}
MASS = {12: 12.000, 13: 13.003, 16: 15.995, 40: 39.962, 84: 83.912, 86: 85.911, 129: 128.905, 132: 131.904,
        136: 135.907}
# database target -> (SRIM gas flag, density g/cm3, elements as (Z, name, stoichiometry, mass))
TARGETS = {
    "H2": (1, 8.99e-5, [(1, "Hydrogen", 1, 1.008)]),
    "D2": (1, 1.8e-4, [(1, "Deuterium", 1, 2.014)]),
    "Mylar": (0, 1.397, [(1, "Hydrogen", 8, 1.008), (6, "Carbon", 10, 12.011), (8, "Oxygen", 4, 15.999)]),
    "Polypropylene": (0, 0.9, [(6, "Carbon", 1, 12.011), (1, "Hydrogen", 2, 1.008)]),
    "Polyethylene": (0, 0.93, [(6, "Carbon", 1, 12.011), (1, "Hydrogen", 2, 1.008)]),
    "CH2": (0, 0.93, [(6, "Carbon", 1, 12.011), (1, "Hydrogen", 2, 1.008)]),
}
#: Data sets not used by the check, with the reason (also written to the file header).
EXCLUDED = {"Diw2003": "polyethylene, 7 points 20-40% below both codes and below every other (CH2)n data set"}


def fetch() -> list:
    rows = []
    for ion in A:
        for target in TARGETS:
            req = urllib.request.Request(API.format(ion=ion, target=target), headers={"User-Agent": "curl/8.0"})
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = json.load(r)
            except urllib.error.HTTPError as err:
                if err.code == 404:  # no data for this pair
                    continue
                raise
            for x in data:
                if x["ref_id"] == "SRIM2013":
                    continue
                iso = int(x["ion_isotope"]) if x["ion_isotope"] else A[ion]
                e, unit = float(x["energy"]), x["energy_unit"]
                assert unit in ("MeV", "keV/u", "MeV/u") and x["stopping_unit_converted"] == "MeV/(mg/cm2)", x
                e_u = e / iso if unit == "MeV" else (e / 1000 if unit == "keV/u" else e)
                if 0.05 <= e_u <= 15:
                    rows.append({"ion": ion, "isotope": iso, "target": target, "phase": x["phase"],
                                 "e_mev_u": f"{e_u:.5g}", "electronic_mev_mg_cm2": x["stopping_power_converted"],
                                 "error_percent": x["percentage_error"] or "", "ref_id": x["ref_id"],
                                 "year": x["publication_year"], "doi": x["doi"] or ""})
    return rows


def srim_value(work, r) -> float:
    gas, dens, els = TARGETS[r["target"]]
    iso = int(r["isotope"])
    e_kev = float(r["e_mev_u"]) * iso * 1000
    text = sr.sr_in(Z[r["ion"]], MASS.get(iso, float(iso)), dens, els, e_kev, e_kev)
    text = text.replace(f"\r\n0    {dens}    1\r\n", f"\r\n{gas}    {dens}    1\r\n")
    return sr.run_sr(work, text)[1][0][1]


def main():
    rows = fetch()
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / "SR Module"
        shutil.copytree(Path(sys.argv[1]) / "SR Module", work)  # never run in (or write to) the install itself
        cache = {}
        for r in rows:
            key = (r["ion"], r["isotope"], r["target"], r["e_mev_u"])
            if key not in cache:
                cache[key] = srim_value(work, r)
            r["srim_electronic_mev_mg_cm2"] = cache[key]
            r["used"] = "no" if r["ref_id"] in EXCLUDED else "yes"
    cols = ["ion", "isotope", "target", "phase", "e_mev_u", "electronic_mev_mg_cm2", "error_percent", "ref_id", "year",
            "doi", "used", "srim_electronic_mev_mg_cm2"]
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        f.write("# tool: IAEA Stopping Power Database (www-nds.iaea.org/stopping, database version 2026-09), "
                "read through its API; measurements only (its \"SRIM2013\" curve rows dropped); cite the original "
                "papers by ref_id/doi\n")
        f.write(f"# produced by: scripts/iaea_stopping_reference.py, {datetime.date.today()}\n")
        f.write("# settings: ions C, O, Ar, Kr, Xe; targets H2, D2 (gas), Mylar, polypropylene, polyethylene, CH2; "
                "0.05-15 MeV/u; electronic stopping in MeV/(mg/cm2) as converted by the database;\n")
        f.write("#   srim_electronic_mev_mg_cm2 is SRIM 2013 SR Module at the same ion, isotope and energy, compound "
                "correction 1, gas flag for H2 and D2 (context only, not judged)\n")
        f.write("# excluded: " + "; ".join(f"{k}: {v}" for k, v in EXCLUDED.items()) + " (used = no)\n")
        w = csv.DictWriter(f, cols, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(len(rows), "rows ->", OUT)


if __name__ == "__main__":
    main()
