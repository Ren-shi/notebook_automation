"""Download the current ENSDF distribution from the National Nuclear Data Center for physim's level schemes.

    python scripts/fetch_ensdf.py            # into ~/.physim/ensdf, or $PHYSIM_ENSDF
    python scripts/fetch_ensdf.py <folder>

The file list is read from https://www.nndc.bnl.gov/ensdfarchivals/distributions/files.json and the latest
distribution (one zip file, about 40 MB) is saved unchanged, with a README.md that records where and when it came
from. physim reads the zip as it is (``physim.nuclear.ensdf``).

The copy stays on this machine. It is not part of the repository or of the installer: permission to redistribute
ENSDF has not been confirmed (backlog item 50).
"""

from __future__ import annotations

import datetime
import json
import os
import sys
import urllib.request
from pathlib import Path

BASE = "https://www.nndc.bnl.gov/ensdfarchivals/distributions"
HEADERS = {"User-Agent": "physim-fetch-ensdf/1.0 (+https://github.com/Ren-shi/notebook_automation)"}

README = """# Local copy of ENSDF

- **File:** `{name}`, saved unchanged.
- **Source:** National Nuclear Data Center, Brookhaven National Laboratory, <{url}>.
- **Retrieved:** {today} by `scripts/fetch_ensdf.py`.
- **Cite:** for up to ten nuclides, the Nuclear Data Sheets evaluation named at the top of each dataset; otherwise
  "Evaluated Nuclear Structure Data File (ENSDF), National Nuclear Data Center, doi:10.18139/nndc.ensdf/1845010",
  with the date of this copy.
- **Do not redistribute:** this copy is for use on this machine. Permission to redistribute ENSDF has not been
  confirmed.
"""


def get(url: str) -> bytes:
    return urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=300).read()


def main(argv: list) -> None:
    if len(argv) > 1:
        out = Path(argv[1])
    else:
        out = Path(os.environ["PHYSIM_ENSDF"]) if os.environ.get("PHYSIM_ENSDF") else Path.home() / ".physim" / "ensdf"
    latest = json.loads(get(f"{BASE}/files.json"))["latest"]
    out.mkdir(parents=True, exist_ok=True)
    for name in latest["files"]:
        url = f"{BASE}/dist{latest['year'][2:]}/{name}"
        target = out / name
        if target.exists():
            print(f"{target} is already there")
        else:
            print(f"downloading {url}")
            data = get(url)
            target.write_bytes(data)
            print(f"saved {target} ({len(data) / 1e6:.1f} MB)")
        (out / "README.md").write_text(README.format(name=name, url=url, today=datetime.date.today().isoformat()),
                                       encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv)
