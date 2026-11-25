"""Reproduce Allmond et al., Phys. Rev. C 90, 034309 (2014), ⁶²Ni Coulomb excitation on carbon, with physim.

    python scripts/reproduce_allmond2014.py [--events N] [--seeds S1 S2 ...]

The example setup ``coulex_ni62_c_allmond`` (62Ni at 111.6 MeV on 1.037 mg/cm² of carbon, carbon recoils in two
rings at 14–28° and 28–44°, nine clovers) is checked six ways (:mod:`physim.nuclear.coulex_reproduction` says how;
the paper's numbers, the GOSIA references, the bands and which of them were fixed without knowing the answer are in
``tests/reference/nuclear/coulex/allmond2014.toml``). What is judged is to all orders; first order is printed as a
second line.

1. the thickness: physim's stopping gives the measured ⁵⁸Ni and ⁶⁰Ni energy losses, and 42.4(8) MeV for ⁶²Ni;
2. σ_C/σ_R per ring against GOSIA as run (E1 polarisation in both; 1 %, missed by 1.3 % before physim had it) and with
   the nuclear masses and the E1 polarisation off in both (0.5 %);
3. closure: simulated and analysed with several seeds, the mean B(E2) back per ring within a quarter of the paper's
   statistical error + 3 Monte Carlo errors of the mean;
4. the paper: its per-ring matrix element through GOSIA's ratio(M) and back through physim's, B(E2) within the
   paper's 4.1 %. The paper gives no yields, so this tests whether physim extracts the same B(E2) from the same
   data as GOSIA did; it is not an independent measurement;
5. Q_s = ±0.1 e b and the feeding from ENSDF's higher states change the ratio as in GOSIA (0.5 percentage points);
6. a closure sample as a stored run (segments, not angles), B(E2) back.

Prints the table and returns (from :func:`main` with ``return_result``) the dict of
:func:`physim.nuclear.coulex_reproduction.run` that the validation check ``reproduction_allmond2014`` shows. The exit
code is 0 only if every number is inside its band.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from physim.nuclear import coulex_reproduction as cr  # noqa: E402

REF = ROOT / "tests" / "reference" / "nuclear"


def reproduce(events=None, seeds=None) -> dict:
    r = cr.run(REF, events=events, seeds=seeds)
    if r is None:
        raise SystemExit(f"{REF / cr.REFERENCE_FILE} not found")
    return r


def main(argv=None, return_result: bool = False):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--events", type=int, default=None, help="events per seed of the closure (the reference file's)")
    ap.add_argument("--seeds", type=int, nargs="+", default=None, help="seeds of the closure (the reference file's)")
    args = ap.parse_args(argv)
    r = reproduce(args.events, args.seeds)
    print("Allmond et al., PRC 90, 034309 (2014): 62Ni on natC, through physim\n")
    print(cr.table(r))
    print("\nCheck 4 rebuilds the paper's ratios from its matrix elements through GOSIA's ratio(M): the paper gives no "
          "yields, so it tests whether physim extracts the same B(E2) from the same data as GOSIA did.")
    print("\n" + cr.verdict(r))
    if return_result:
        return r
    return 0 if r["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
