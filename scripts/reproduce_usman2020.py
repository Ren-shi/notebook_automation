"""Reproduce Usman et al., Nucl. Instrum. Methods B 469, 42 (2020) end to end through the planner (backlog 75, D10).

    python scripts/reproduce_usman2020.py [--no-save]

For every Ho foil of the paper's stack (and a few energies beyond it), the example setup ``fusion_tm165`` is opened
in the planner, its beam energy set to the foil's, its channel table filled from physim's own fusion and
evaporation, the residues simulated, and the production cross sections inferred from the residues made
(:mod:`physim.nuclear.reproduction` says how, and fixes the comparisons). Writes:

- ``tests/reference/nuclear/reproduction_usman2020.csv``: one row per foil and channel, with the provenance header;
- ``docs/images/reproduction-usman2020.png``: the four excitation functions, physim, the three measurements and
  PACE4 (from ``pace4_excitation_4He_165Ho.csv`` when it is there);

and prints the comparisons and the beam's energy at each Ho foil through the inferred stack. The validation check
``reproduction_usman2020`` reruns the same reproduction (same seed, same numbers) and judges it.
"""

from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from physim.nuclear import reproduction as rp  # noqa: E402

REF = ROOT / "tests" / "reference" / "nuclear"
FIGURE = ROOT / "docs" / "images" / "reproduction-usman2020.png"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--no-save", action="store_true", help="print the comparisons only")
    args = ap.parse_args(argv)

    measured = rp.measured(REF)
    rep = rp.run(reference=REF, progress=lambda e: print(f"  foil at {e:g} MeV", flush=True))
    c = rp.compare(rep.rows, measured)

    print("\nPeak energies (MeV): measured, physim, tolerance")
    for x in c["peaks"]:
        print(f"  {x['channel']}: {x['measured_mev']:.2f}  {x['physim_mev']:.2f}{' (edge)' if x['physim_at_edge'] else ''}"
              f"  ±{x['tolerance_mev']:.1f}  {'pass' if x['passed'] else 'FAIL' if x['judged'] else 'reported'}")
    print("\nStrongest channel per foil: measured, physim")
    for x in c["strongest"]:
        tag = ("pass" if x["passed"] else "FAIL") if x["judged"] else (
            "below 1.1 V_B" if not x["above_barrier"] else "unresolved" if not x["resolved"] else "5n unmeasured")
        print(f"  {x['e_lab_mev']:5.1f}: {x['measured']:3s} {x['physim']:4s} {tag}")
    print("\nMagnitudes: physim/Usman, deviation from the band of three measurements")
    for x in c["magnitudes"]:
        band = x["band_mb"]
        print(f"  {x['e_lab_mev']:5.1f} {x['channel']}: physim {x['physim_mb']:8.1f}  Usman {x['measured_mb']:8.1f}"
              f"  ratio {x['ratio']:6.2f}  band {band[0] if band else float('nan'):7.1f}–"
              f"{band[1] if band else float('nan'):7.1f}  dev {x['deviation']:+.0%}  "
              f"{('pass' if x['passed'] else 'FAIL') if x['judged'] else x['why']}")
    print("\nσ(3n)/σ(4n) at 45–50 MeV")
    for x in c["ratio_3n_4n"]:
        print(f"  {x['ref']:14s} {x['e_lab_mev']:5.1f}: measured {x['measured']:.3f}  physim {x['physim']:.3f}")
    print("\nBeam through the inferred stack from 50.4 MeV (mean energy in each Ho foil): physim, paper")
    for s, e in zip(rp.stack_energies(), rp.foil_energies(measured)):
        print(f"  foil {s['foil']:2d}: {s['e_mean_mev']:6.2f}  {e:5.1f}  {s['e_mean_mev'] - e:+.2f}")
    print("\n" + rp.summary(c))
    print("passed" if c["passed"] else "FAILED", c["ok"])

    if args.no_save:
        return 0
    lines = rp.result_lines(rep, measured, f"scripts/reproduce_usman2020.py (Claude Code), {date.today().isoformat()}",
                            c)
    out = REF / rp.RESULT_FILE
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", out)
    pace4 = REF / rp.PACE4_FILE
    fig = rp.plot(rep.rows, measured, rp.read_pace4(pace4) if pace4.exists() else None)
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=130)
    print("wrote", FIGURE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
