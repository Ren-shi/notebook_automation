"""Download the NIST PSTAR and ASTAR stopping-power tables and write them as CSV for physim.

    python scripts/fetch_star_tables.py <download dir>

Each material is requested from https://physics.nist.gov/cgi-bin/Star/ap_table.pl on NIST's default energy grid
(protons: 1 keV to 10 GeV; alphas: 1 keV to 1 GeV). The raw HTML answers are kept in the
download directory; python/physim/nuclear/data/nist_pstar.csv and nist_astar.csv get the values, copied as
published. Only the layout changes.
"""

from __future__ import annotations

import csv
import html
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

URL = "https://physics.nist.gov/cgi-bin/Star/ap_table.pl"
HEADERS = {"User-Agent": "physim-fetch-star-tables/1.0 (+https://github.com/Ren-shi/notebook_automation)"}
FORM = "https://physics.nist.gov/PhysRefData/Star/Text/{prog}.html"
OUT = Path(__file__).resolve().parents[1] / "python" / "physim" / "nuclear" / "data"
COLUMNS = ["energy_mev", "electronic_mev_cm2_g", "nuclear_mev_cm2_g", "total_mev_cm2_g", "csda_range_g_cm2",
           "projected_range_g_cm2", "detour_factor"]


def materials(prog: str, download: Path) -> list[tuple[str, str]]:
    page = download / f"{prog.lower()}_form.html"
    if not page.exists():
        req = urllib.request.Request(FORM.format(prog=prog), headers=HEADERS)
        page.write_bytes(urllib.request.urlopen(req, timeout=60).read())
    text = page.read_text(encoding="latin-1")
    found = re.findall(r'<option value="(\d+)">\s*([^\n<]+)', text)
    return [(no, re.sub(r"\s+", " ", name).strip()) for no, name in found]


def fetch(prog: str, matno: str, download: Path) -> str:
    path = download / f"{prog.lower()}_{matno}.html"
    if not path.exists():
        body = urllib.parse.urlencode({"prog": prog, "matno": matno, "ShowDefault": "on", "GraphType": "None"})
        req = urllib.request.Request(URL, data=body.encode(), method="POST", headers=HEADERS)
        path.write_bytes(urllib.request.urlopen(req, timeout=60).read())
        time.sleep(1.0)  # be gentle with NIST's server
    return path.read_text(encoding="latin-1")


def rows(page: str) -> list[list[str]]:
    out = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", page, re.S | re.I):
        cells = [html.unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)]
        if len(cells) == 7 and re.match(r"^\d\.\d+E[+-]\d+$", cells[0]):
            out.append(cells)
    return out


def main(download_dir: str) -> None:
    download = Path(download_dir)
    download.mkdir(parents=True, exist_ok=True)
    for prog in ("PSTAR", "ASTAR"):
        table, grid = [], None
        mats = materials(prog, download)
        if len(mats) != 74:
            raise SystemExit(f"{prog}: expected 74 materials, found {len(mats)}")
        for matno, name in mats:
            data = rows(fetch(prog, matno, download))
            energies = [cells[0] for cells in data]
            grid = grid or energies
            if energies != grid:  # every material is tabulated on the same default grid
                raise SystemExit(f"{prog} {matno} ({name}): energy grid differs ({len(energies)} energies)")
            for cells in data:
                table.append({"material_no": matno, "material": name, **dict(zip(COLUMNS, cells))})
        path = OUT / f"nist_{prog.lower()}.csv"
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["material_no", "material", *COLUMNS], lineterminator="\n")
            w.writeheader()
            w.writerows(table)
        print(f"wrote {path} ({len(mats)} materials x {len(grid)} energies)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
