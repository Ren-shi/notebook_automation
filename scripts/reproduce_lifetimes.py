"""Rerun the three published lifetime reproductions of the lifetime validation, optionally with larger runs.

    python scripts/reproduce_lifetimes.py                 # the validation suite's sizes (about 35 s)
    python scripts/reproduce_lifetimes.py --scale 3       # three times the events (statistics √3 better)
    python scripts/reproduce_lifetimes.py ru98 --seed 2   # one of them, another seed

- ``ru98``: Radeck et al., PRC 85, 014301 (2012), RDDS after Coulex, τ(2⁺₁) with the deorientation correction.
- ``cr50``: Beckers et al., NIM A 1042, 167418 (2022), three-foil RDDS after fusion, τ(2⁺₁) and τ(4⁺₁) by γγ DDCM.
- ``sm143``: Rajbanshi et al., PLB 782, 143 (2018), DSAM after fusion on a Au backing, line-shape fits.

The logic is in ``physim.nuclear.lifetime_reproduction``; the published values and setups are in
``tests/reference/nuclear/lifetimes/published/``. Each prints the comparison and whether it is within the
validation check's tolerance. ``--side-steps 5`` runs ¹⁴³Sm with LINESHAPE's five-step side feeding (slow).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))


def main(argv=None) -> int:
    from physim.nuclear import lifetime_reproduction as lr

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("which", nargs="*", help="ru98, cr50, sm143 (default all three)")
    ap.add_argument("--scale", type=float, default=1.0, help="multiply the number of events")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--side-steps", type=int, default=None, help="¹⁴³Sm side-feeding steps (default 1)")
    args = ap.parse_args(argv)
    which = args.which or ["ru98", "cr50", "sm143"]
    unknown = set(which) - {"ru98", "cr50", "sm143"}
    if unknown:
        ap.error(f"unknown reproduction(s): {', '.join(sorted(unknown))}")
    ok = True
    for name in which:
        t0 = time.time()
        if name == "ru98":
            r = lr.ru98(particles=int(2_000_000 * args.scale), seed=args.seed)
        elif name == "cr50":
            r = lr.cr50(residues=int(400_000 * args.scale), seed=args.seed)
            r["detail"] += (f". Singles with the 4⁺ feeder's correction (closure, simulated 12.20 ps): "
                            f"{r['singles'][0]:.2f} ± {r['singles'][1]:.2f} ps")
        else:
            r = lr.sm143(residues=int(20_000 * args.scale), n=int(3000 * args.scale), seed=args.seed,
                         side_steps=args.side_steps)
        print(f"{name} ({time.time() - t0:.0f} s): {'within' if r['passed'] else 'OUTSIDE'} the tolerance "
              f"(worst {r['worst']:.0%})")
        print("  " + r["detail"])
        ok &= bool(r["passed"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
