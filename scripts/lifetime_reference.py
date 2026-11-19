"""Write the reference files for the lifetime checks (plunger RDDS and DSAM) in ``tests/reference/nuclear/lifetimes/``.

    python scripts/lifetime_reference.py analytic
    python scripts/lifetime_reference.py ahkin <AhKin build folder>
    python scripts/lifetime_reference.py apcad <APCADStopSim folder>
    python scripts/lifetime_reference.py napytau <NaPyTau checkout> --python <python with numpy and scipy>

``analytic`` needs nothing but numpy and imports nothing from physim, so its numbers are independent of the code
they check. The other three take an outside program's files (never shipped with physim; see the README in the
reference folder for where each comes from, its licence and how it was built) and copy or convert what the checks
need. CI only reads the files written here.

- ``analytic/rdds_one_level.csv``: one level, no feeding, a recoil at one speed along the beam. The share of decays
  after a flight of d is exp(-d / (γ β c τ)) (time dilation included; the column ``unshifted_no_dilation`` leaves it
  out).
- ``analytic/rdds_two_levels.csv``: a level fed from one level above (Bateman's equations): N₁ = P₁ e^(−t/τ₁),
  N₂ = P₂ e^(−t/τ₂) + P₁ τ₂/(τ₁ − τ₂) (e^(−t/τ₁) − e^(−t/τ₂)). The lower transition's unshifted share is
  (N₁ + N₂)/(P₁ + P₂), the upper one's N₁/P₁ (all of level 1 decays through level 2).
- ``analytic/dsam_exponential.csv``: stopping proportional to the speed, so v(t) = v₀ e^(−t/α); the mean v/v₀ at
  emission (the attenuation F) is α/(α + τ), and v/v₀ is distributed as (α/τ) x^(α/τ − 1) on (0, 1].
- ``ahkin/y83_in_ta/``: AhKin's own example (⁸³Y residues in Ta, the 1193 keV line), its inputs and outputs, and the
  lifetime it finds (0.1313 +0.0339 −0.0327 ps).
- ``apcad/``: APCAD's two example spectra (continuous angles 20°–160°, τ = 1 ps, and the same with two contaminant
  lines) converted to CSV, with the setup that made them.
- ``napytau/``: NaPyTau's DDCM answers on ``rdds_one_level.csv`` (central values; its uncertainty propagation gives
  meaningless numbers and is not used).
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "reference" / "nuclear" / "lifetimes"
C_UM_PS = 299.792458  # speed of light in µm/ps


def _write_csv(path: Path, header: list[str], comments: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for c in comments:
            f.write(f"# {c}\n")
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(f"{v:.10g}" if isinstance(v, float) else str(v) for v in r) + "\n")


# -- analytic -------------------------------------------------------------------------------------------------------

DISTANCES_UM = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200, 300, 500, 700, 1000, 2000, 5000]


def analytic() -> None:
    today = date.today().isoformat()
    rows = []
    for tau in (1.0, 10.0, 100.0):
        for beta in (0.02, 0.05, 0.10):
            gamma = 1 / np.sqrt(1 - beta**2)
            v = beta * C_UM_PS
            for d in DISTANCES_UM:
                t = d / v
                rows.append((tau, beta, float(d), t, float(np.exp(-t / (gamma * tau))), float(np.exp(-t / tau))))
    _write_csv(OUT / "analytic" / "rdds_one_level.csv",
               ["tau_ps", "beta", "distance_um", "flight_time_ps", "unshifted", "unshifted_no_dilation"],
               [f"Written by scripts/lifetime_reference.py analytic on {today}.",
                "One level, no feeding, recoil along the beam at speed beta; flight time = d / (beta c) in the lab.",
                "unshifted = exp(-t / (gamma tau)) (decays after the stopper: time dilation included); shifted = 1 - unshifted."],
               rows)

    rows = []
    for tau1, tau2 in ((2.0, 10.0), (10.0, 10.5), (30.0, 5.0)):
        for p1, p2 in ((1.0, 0.0), (0.4, 0.6)):
            beta = 0.05
            v = beta * C_UM_PS
            for d in DISTANCES_UM:
                t = d / v  # no time dilation here: the Bateman part alone
                n1 = p1 * np.exp(-t / tau1)
                n2 = p2 * np.exp(-t / tau2) + p1 * tau2 / (tau1 - tau2) * (np.exp(-t / tau1) - np.exp(-t / tau2))
                rows.append((tau1, tau2, p1, p2, beta, float(d), t, float(n1 / p1), float((n1 + n2) / (p1 + p2))))
    _write_csv(OUT / "analytic" / "rdds_two_levels.csv",
               ["tau1_ps", "tau2_ps", "p1", "p2", "beta", "distance_um", "flight_time_ps", "unshifted_upper",
                "unshifted_lower"],
               [f"Written by scripts/lifetime_reference.py analytic on {today}.",
                "Level 1 (populated with p1) decays only to level 2 (populated directly with p2), which decays to the ground state.",
                "Bateman: N1 = p1 exp(-t/tau1); N2 = p2 exp(-t/tau2) + p1 tau2/(tau1 - tau2) (exp(-t/tau1) - exp(-t/tau2)).",
                "unshifted_upper = N1/p1; unshifted_lower = (N1 + N2)/(p1 + p2). No time dilation (see rdds_one_level.csv for it)."],
               rows)

    rows = []
    for alpha in (0.2, 0.5, 1.0):
        for tau in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0):
            k = alpha / tau
            edges = np.linspace(0.0, 1.0, 21)
            share = edges[1:] ** k - edges[:-1] ** k  # CDF of x = v/v0 is x^k
            rows.append((alpha, tau, alpha / (alpha + tau), ";".join(f"{s:.8g}" for s in share)))
    _write_csv(OUT / "analytic" / "dsam_exponential.csv",
               ["alpha_ps", "tau_ps", "attenuation_F", "share_v_over_v0_in_20_bins"],
               [f"Written by scripts/lifetime_reference.py analytic on {today}.",
                "Stopping proportional to speed: v(t) = v0 exp(-t/alpha). Decay at an exponential time with mean tau.",
                "F = <v/v0> at emission = alpha/(alpha + tau). x = v/v0 has density (alpha/tau) x^(alpha/tau - 1) on (0, 1];",
                "the last column is its share in 20 equal bins of x from 0 to 1 (semicolon-separated)."],
               rows)
    print("analytic: written to", OUT / "analytic")


# -- AhKin ----------------------------------------------------------------------------------------------------------

def ahkin(build: Path) -> None:
    """Run AhKin's example (prob_cos, stopping, AhKin_B) in a fresh copy of its in_files and keep the results."""
    src = build / "extracto_AhKin"
    exe = {p: src / f"{p}.exe" for p in ("prob_cos", "stopping", "AhKin_B")}
    work = OUT / "ahkin" / "y83_in_ta"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(src / "in_files", work)
    run = lambda name, text: subprocess.run([str(exe[name])] + (["-i"] if name == "AhKin_B" else []), input=text,
                                            cwd=work, capture_output=True, text=True, check=True).stdout
    run("prob_cos", "in_cos.dat\n100\n1e6\n")
    run("stopping", "in_stp.dat\n")
    out = run("AhKin_B", "in_AhKin.dat\nin_casc.dat\n0.1\n0.2\n")
    lines = [ln for ln in out.splitlines() if ln.strip()]
    i = next(k for k, ln in enumerate(lines) if ln.lstrip().startswith("t a u:"))
    tau = float(lines[i].split()[3])
    up, down = float(lines[i + 1].split()[1]), float(lines[i + 2].split()[1])
    best = f"ech{lines[i].split()[3]}.dat"
    for f in work.glob("ech*.dat"):  # keep only the best-fitting line shape
        if f.name != best:
            f.unlink()
    (work / "result.json").write_text(json.dumps(
        {"tau_ps": tau, "plus_ps": up, "minus_ps": down, "best_lineshape": best,
         "published": "Garzón, Rodríguez, Cristancho, Tao, Comput. Phys. Commun. (2020), AhKin; readme example: "
                      "0.131308 +0.0339269 -0.0326754 ps", "run_on": date.today().isoformat()}, indent=2) + "\n",
        encoding="utf-8")
    print("ahkin: tau", tau, "+", up, "-", down, "->", work)


# -- APCAD ----------------------------------------------------------------------------------------------------------

def apcad(folder: Path) -> None:
    ex = folder / "Example"
    out = OUT / "apcad"
    out.mkdir(parents=True, exist_ok=True)
    for name in ("OneShape_example", "ThreeShapes_example"):
        data = np.loadtxt(ex / "Input" / f"{name}.spc")
        _write_csv(out / f"{name}.csv", ["angle_deg", "energy_kev", "counts"],
                   [f"APCAD example spectrum Example/Input/{name}.spc (Stahl et al., Mendeley Data 4hzzcvx2p7, CC BY 4.0),",
                    "converted by scripts/lifetime_reference.py apcad; the fourth column (sqrt(counts)) is dropped."],
                   [(float(a), float(e), int(round(c))) for a, e, c, _ in data])
    for f in ("Input/Detector_example.det", "Stop/Macro_StopSim_example.mac", "Stop/XSec_example.dat", "README.txt"):
        shutil.copy(ex / f, out / Path(f).name)
    print("apcad: written to", out)


# -- NaPyTau --------------------------------------------------------------------------------------------------------

NAPYTAU_DRIVER = r'''
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
from napytau.core import core
from napytau.core.tau import calculate_tau_i_values
from napytau.core.polynomials import calculate_polynomial_coefficients_for_fit
from napytau.import_export.model.dataset import DataSet
from napytau.import_export.model.datapoint import Datapoint
from napytau.import_export.model.datapoint_collection import DatapointCollection
from napytau.import_export.model.relative_velocity import RelativeVelocity
from napytau.util.model.value_error_pair import ValueErrorPair as V
cases = json.loads(sys.stdin.read())
out = []
for c in cases:
    beta = c["beta"]
    # NaPyTau computes t = distance / (beta * c) with c in m/s; distances are passed scaled so that t is in ps.
    scale = 299792458.0 / 299.792458
    pts = [Datapoint(distance=V(d * scale, 0.0), calibration=V(1.0, 0.0), shifted_intensity=V(s, 0.01 * s + 1e-4),
                     unshifted_intensity=V(u, 0.01 * u + 1e-4)) for d, s, u in zip(c["d"], c["s"], c["u"])]
    ds = DataSet(relative_velocity=V(RelativeVelocity(beta), 0.0), datapoints=DatapointCollection(pts))
    for deg in (2, 3, 4):
        coef = calculate_polynomial_coefficients_for_fit(ds, deg)
        tau_i = calculate_tau_i_values(ds, coef)
        tau, _ = core.calculate_lifetime_for_fit(ds, deg)
        out.append({"tau_ps": c["tau"], "beta": beta, "degree": deg, "tau_final_ps": float(tau),
                    "tau_i_ps": [float(x) for x in tau_i]})
print(json.dumps(out))
'''


def napytau(repo: Path, python: str) -> None:
    import csv
    rows = list(csv.DictReader(ln for ln in (OUT / "analytic" / "rdds_one_level.csv").open(encoding="utf-8")
                               if not ln.startswith("#")))
    cases = []
    for tau, beta in ((10.0, 0.05), (100.0, 0.05), (10.0, 0.10)):
        sel = [r for r in rows if float(r["tau_ps"]) == tau and float(r["beta"]) == beta]
        # the region of sensitivity: flight times from about 0.05 to 3 tau
        sel = [r for r in sel if 0.05 * tau <= float(r["flight_time_ps"]) <= 3 * tau]
        u = [float(r["unshifted_no_dilation"]) for r in sel]
        cases.append({"tau": tau, "beta": beta, "d": [float(r["distance_um"]) for r in sel], "u": u,
                      "s": [1 - x for x in u]})
    res = subprocess.run([python, "-I", "-c", NAPYTAU_DRIVER, str(repo)], input=json.dumps(cases),
                         capture_output=True, text=True, check=True)
    results = json.loads(res.stdout)
    commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    _write_csv(OUT / "napytau" / "ddcm_one_level.csv",
               ["tau_ps", "beta", "polynomial_degree", "tau_final_ps", "tau_i_ps"],
               [f"NaPyTau (github.com/h-mayr/NaPyTau_Public, commit {commit}, EUPL-1.2) on analytic/rdds_one_level.csv",
                "(unshifted_no_dilation, distances with flight times 0.05-3 tau, 1% intensity errors), singles DDCM with one",
                "polynomial over all distances; tau_i per distance (semicolon-separated). Its uncertainties are not usable.",
                f"Written by scripts/lifetime_reference.py napytau on {date.today().isoformat()}."],
               [(r["tau_ps"], r["beta"], r["degree"], r["tau_final_ps"], ";".join(f"{x:.6g}" for x in r["tau_i_ps"]))
                for r in results])
    print("napytau: written to", OUT / "napytau")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="what", required=True)
    sub.add_parser("analytic")
    a = sub.add_parser("ahkin"); a.add_argument("build", type=Path)
    a = sub.add_parser("apcad"); a.add_argument("folder", type=Path)
    a = sub.add_parser("napytau"); a.add_argument("repo", type=Path); a.add_argument("--python", required=True)
    args = ap.parse_args(argv)
    if args.what == "analytic":
        analytic()
    elif args.what == "ahkin":
        ahkin(args.build)
    elif args.what == "apcad":
        apcad(args.folder)
    else:
        napytau(args.repo, args.python)
    return 0


if __name__ == "__main__":
    sys.exit(main())
