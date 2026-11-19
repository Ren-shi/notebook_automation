"""Measured electronic stopping of heavy ions (Ar and heavier) in elemental solids, from the IAEA Stopping Power
Database, written to tests/reference/nuclear/stopping_measured_heavy.csv with SRIM 2013's value for each point.

Usage: python scripts/iaea_heavy_ion_reference.py "D:/Projects/srim/SRIM2013" [--no-srim]

Keeps measurements only (the database also serves SRIM 2013 curve values under the reference id "SRIM2013"; they
are dropped), solid targets, 0.005-10 MeV/u, MeV/(mg/cm2). SRIM's value comes from one SR Module table per ion and
target (10 keV/u to 12 MeV/u, see srim_reference.py), interpolated log-log to each point. This is the reference for
physim's stopping of slow heavy ions (below about 0.5 MeV/u), where recoils of Coulomb excitation and
fusion-evaporation slow down in targets, backings and plunger stoppers.
"""
import argparse
import csv
import datetime
import json
import math
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "tests" / "reference" / "nuclear" / "stopping_measured_heavy.csv"
API = "https://www-nds.iaea.org/stopping/api/data/ion:{ion}/target:{target}"
#: Ions: symbol -> (Z, the most common mass number, used when the database gives no isotope).
IONS = {"Ar": (18, 40), "Ti": (22, 48), "Fe": (26, 56), "Ni": (28, 58), "Cu": (29, 63), "Ge": (32, 74),
        "Br": (35, 79), "Kr": (36, 84), "Nb": (41, 93), "Mo": (42, 98), "Ru": (44, 102), "Ag": (47, 107),
        "I": (53, 127), "Xe": (54, 132), "Sm": (62, 152), "Er": (68, 166), "Ta": (73, 181), "Au": (79, 197),
        "Pb": (82, 208), "Bi": (83, 209), "U": (92, 238)}
#: Elemental solid targets: symbol -> (Z, SRIM name, mass in u, density g/cm3).
TARGETS = {"C": (6, "Carbon", 12.011, 2.253), "Mg": (12, "Magnesium", 24.305, 1.74), "Al": (13, "Aluminum", 26.982, 2.70),
           "Si": (14, "Silicon", 28.086, 2.33), "Ti": (22, "Titanium", 47.867, 4.51), "Ni": (28, "Nickel", 58.693, 8.90),
           "Cu": (29, "Copper", 63.546, 8.96), "Ge": (32, "Germanium", 72.63, 5.32), "Ag": (47, "Silver", 107.87, 10.5),
           "Sn": (50, "Tin", 118.71, 7.31), "Ta": (73, "Tantalum", 180.95, 16.65), "Au": (79, "Gold", 196.97, 19.32),
           "Pb": (82, "Lead", 207.2, 11.35)}
E_MIN, E_MAX = 0.005, 10.0  # MeV/u


def fetch(cache: Path) -> list:
    rows = []
    for ion in IONS:
        for target in TARGETS:
            path = cache / f"{ion}_{target}.json"
            if path.exists():
                data = json.loads(path.read_text(encoding="utf-8"))
            else:
                req = urllib.request.Request(API.format(ion=ion, target=target), headers={"User-Agent": "curl/8.0"})
                for attempt in range(12):
                    try:
                        with urllib.request.urlopen(req, timeout=60) as r:
                            data = json.load(r)
                        break
                    except urllib.error.HTTPError as err:
                        if err.code == 404:  # no data for this pair
                            data = []
                            break
                        if err.code != 429 or attempt == 11:
                            raise
                        time.sleep(15 * (attempt + 1))  # the server limits the request rate
                time.sleep(2)
                path.write_text(json.dumps(data), encoding="utf-8")
            for x in data:
                if x["ref_id"] == "SRIM2013" or x.get("phase") not in ("s", None):
                    continue
                unit = x["energy_unit"]
                if x["stopping_unit_converted"] != "MeV/(mg/cm2)" or unit not in ("MeV/u", "keV/u", "MeV"):
                    continue
                a = int(x["ion_isotope"]) if x["ion_isotope"] else IONS[ion][1]
                e = float(x["energy"])
                e_u = e if unit == "MeV/u" else e / 1000 if unit == "keV/u" else e / a
                if not E_MIN <= e_u <= E_MAX:
                    continue
                rows.append({"ion": ion, "z1": IONS[ion][0], "a1": a, "target": target, "z2": TARGETS[target][0],
                             "e_mev_u": e_u, "stopping": float(x["stopping_power_converted"]),
                             "error_percent": x["percentage_error"] or "", "ref_id": x["ref_id"],
                             "doi": x["doi"] or ""})
    return rows


def srim_values(folder: str, rows: list) -> None:
    import srim_reference as sr
    from physim.external.srim import Srim

    pairs = sorted({(r["ion"], r["a1"], r["target"]) for r in rows})
    with Srim(folder) as srim:
        for ion, a, target in pairs:
            z1 = IONS[ion][0]
            z2, name, mass, density = TARGETS[target]
            text = sr.sr_in(z1, float(a), density, [(z2, name, 1, mass)], 10 * a, 12000 * a)
            _, table = sr.run_sr(srim, text)
            es = [math.log(t[0] / a / 1000) for t in table]
            ss = [math.log(t[1]) for t in table]
            for r in rows:
                if (r["ion"], r["a1"], r["target"]) != (ion, a, target):
                    continue
                x = math.log(r["e_mev_u"])
                if x < es[0]:  # below the table: SRIM's S_e ∝ v there
                    r["srim"] = math.exp(ss[0] + 0.5 * (x - es[0]))
                    continue
                k = max(i for i in range(len(es) - 1) if es[i] <= x) if x < es[-1] else len(es) - 2
                f = (x - es[k]) / (es[k + 1] - es[k])
                r["srim"] = math.exp(ss[k] + f * (ss[k + 1] - ss[k]))
            print(f"{ion}{a} in {target}: {sum(1 for r in rows if (r['ion'], r['a1'], r['target']) == (ion, a, target))} points")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("srim", nargs="?")
    ap.add_argument("--no-srim", action="store_true")
    ap.add_argument("--cache", default=str(REPO / "tools" / "iaea-cache"))
    args = ap.parse_args()
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    rows = fetch(cache)
    rows.sort(key=lambda r: (r["z1"], r["z2"], r["e_mev_u"], r["ref_id"]))
    if not args.no_srim:
        sys.path.insert(0, str(REPO / "python"))
        srim_values(args.srim, rows)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        f.write("# tool: IAEA Stopping Power Database (https://www-nds.iaea.org/stopping/), measured values only;\n"
                "#   SRIM 2013 SR Module for the srim column (compound correction 1, one table per ion and target)\n")
        f.write(f"# produced by: scripts/iaea_heavy_ion_reference.py, {datetime.date.today().isoformat()}\n")
        f.write(f"# settings: ions {', '.join(IONS)} in solid {', '.join(TARGETS)}; {E_MIN}-{E_MAX} MeV/u; electronic\n"
                "#   stopping as the database converts it, MeV/(mg/cm2); a1 is the database's isotope (the most common\n"
                "#   one where it gives none). Cite the original papers by ref_id and doi.\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["ion", "z1", "a1", "target", "z2", "e_mev_u", "stopping_mev_mg_cm2", "error_percent", "srim_mev_mg_cm2",
                    "ref_id", "doi"])
        for r in rows:
            w.writerow([r["ion"], r["z1"], r["a1"], r["target"], r["z2"], f"{r['e_mev_u']:.6g}", f"{r['stopping']:.5g}",
                        r["error_percent"], f"{r['srim']:.5g}" if "srim" in r else "", r["ref_id"], r["doi"]])
    print(f"wrote {OUT} ({len(rows)} points)")


if __name__ == "__main__":
    main()
