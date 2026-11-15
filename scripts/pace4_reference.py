"""Write PACE4 reference values for the fusion-evaporation checks in ``physim.nuclear.validation``.

    python scripts/pace4_reference.py "D:/Program Files/LISEcute/LISEcute/PACE4.exe" [--only evaporation|fusion]

Runs only on Windows with PACE4 (installed with LISE++). The install is never written to: each deck is written to a
fresh temporary folder and PACE4 runs there in silent mode (``PACE4.exe -s deck.pace``), which leaves
``deck.html``, ``deck.rtf`` and ``deck.cs4`` next to the deck. **Runs are strictly one at a time**: two PACE4
processes share a temporary file and corrupt each other's output.

For each of the five systems of ``docs/theory/evaporation.md`` it writes two files to ``tests/reference/nuclear/``:

- ``pace4_<beam>_<target>_<E>MeV.csv``: the evaporation at one beam energy, 20 000 cascades. The header holds E*,
  σ_fus, σ_J, the mean CM energies and the multiplicities of n, p and α; the rows are the residues (Z, A, events,
  σ in mb). ``evaporation_channels_vs_pace4`` and ``evaporation_spectra_vs_pace4`` read these.
- ``pace4_fusion_<beam>_<target>.csv``: the excitation function, one row per beam energy from near the barrier
  upward: E_cm, E*, PACE4's Bass barrier and σ_fus, ℓ_grazing, and σ_J (semicolon-separated, J = 0, 1, ...).
  1000 cascades each (σ_fus and σ_J do not depend on the cascades). ``fusion_vs_pace4`` reads these.

The beam energies of the excitation functions were chosen once from physim's Bass barrier V_B (E_cm/V_B = 1.0, 1.05,
1.1, 1.2, 1.3, 1.5, 1.75, 2.0, rounded to whole MeV in the lab) and are fixed here, so a rerun gives the same rows.

PACE4 settings: its defaults as in the example deck ``files/examples/48Ca_124Sn_e300.pace`` (fusion cross section
and σ_J from the Bass model, L-diffuseness 2, level-density parameter a = A/10), with the ground-state spins of beam
and target set to 0, as physim does not couple them. CI reads the CSV files and never needs PACE4.
"""

from __future__ import annotations

import argparse
import html
import re
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path

REF = Path(__file__).resolve().parents[1] / "tests" / "reference" / "nuclear"

#: (beam, Z, A), (target, Z, A), the evaporation's beam energy (lab MeV), the excitation function's beam energies.
SYSTEMS = [
    (("4He", 2, 4), ("165Ho", 67, 165), 50.0, [18, 19, 20, 22, 24, 27, 31, 35, 40, 50]),
    (("12C", 6, 12), ("150Nd", 60, 150), 70.0, [49, 51, 53, 58, 63, 70, 85, 97]),
    (("16O", 8, 16), ("144Sm", 62, 144), 90.0, [69, 72, 75, 82, 90, 102, 120, 136]),
    (("32S", 16, 32), ("92Mo", 42, 92), 150.0, [115, 120, 126, 137, 150, 172, 200, 229]),
    (("48Ca", 20, 48), ("124Sn", 50, 124), 200.0, [160, 167, 175, 191, 200, 207, 239, 279, 318]),
]
EVAPORATION_CASCADES = 20000
FUSION_CASCADES = 1000


def deck(cascades: int, zp: int, ap: int, zt: int, at: int, e_lab: float) -> str:
    """A PACE4 input deck: the example deck's switches with this system, energy and number of cascades."""
    return (f"{cascades:5d}    1    1    0    0    1    0   10    0    0    0\n"
            f"{zp:5d}{ap:5d}{zt:5d}{at:5d}  0.0  0.0   0.00000\n"
            f"{e_lab:10.5f}   0.00000    0    2    0    0 0 3.00\n")


def run(pace4: Path, text: str) -> tuple:
    """Run one deck in a fresh temporary folder; (deck text, HTML output)."""
    with tempfile.TemporaryDirectory(prefix="pace4_") as tmp:
        work = Path(tmp)
        (work / "run.pace").write_text(text, encoding="ascii", newline="\r\n")
        subprocess.run([str(pace4), "-s", "run.pace"], cwd=work, check=True, timeout=1800)
        out = work / "run.html"
        if not out.exists():
            raise RuntimeError(f"PACE4 wrote no output for\n{text}")
        return text, out.read_text(encoding="latin-1")


# --- reading PACE4's HTML output ------------------------------------------------------------------------------------


def _cells(row: str) -> list:
    return [html.unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row,
                                                                                  re.S)]


def _rows(text: str) -> list:
    return [_cells(r) for r in re.findall(r"<tr[^>]*>(.*?)</tr>", text, re.S)]


def _after(text: str, label: str) -> str:
    """The value cell of a two-cell table row whose first cell contains ``label``."""
    for cells in _rows(text):
        if len(cells) >= 2 and label in cells[0]:
            return cells[1]
    raise KeyError(label)


def parse(text: str) -> dict:
    """The numbers the checks use, from one PACE4 HTML output."""
    plain = html.unescape(re.sub(r"<[^>]+>", " ", text))
    out = {"version": re.search(r"Version\s+([\d.]+)", plain).group(1)}
    m = re.search(r"Bass fusion xsection for E =\s*([\d.]+)\s*MeV is\s*([\d.eE+-]+)\s*mb", plain)
    out["sigma_fus_mb"] = float(m.group(2))
    m = re.search(r"Fusion radius =\s*([\d.]+)\s*fm\. Barrier height is\s*([\d.]+)", plain)
    out["fusion_radius_fm"], out["barrier_mev"] = float(m.group(1)), float(m.group(2))
    out["e_cm_mev"] = float(_after(text, "Center of mass energy"))
    out["excitation_mev"] = float(_after(text, "Compound nucleus excitation energy"))
    out["beta_cn"] = float(_after(text, "Compound nucleus velocity/c"))
    out["l_grazing"] = float(_after(text, "Fusion L-grazing"))
    # Partial cross sections: a table of (J, σ_J) pairs, five pairs per row, J running down the columns.
    part = text.split("Partial cross sections (mb)", 1)[1].split("<p>", 1)[0]
    sig = {}
    for cells in _rows(part)[1:]:
        for j, s in zip(cells[0::2], cells[1::2]):
            if j and s:
                sig[int(j)] = float(s)
    out["sigma_j_mb"] = [sig.get(j, 0.0) for j in range(max(sig) + 1)] if sig else []
    # Residues.
    res = text.split("Yields of residual nuclei", 1)[1].split("</table>", 1)[0]
    out["residues"] = []
    for cells in _rows(res)[1:]:
        if cells and cells[0].isdigit():
            z, n = int(cells[0]), int(cells[1])
            out["residues"].append((z, z + n, int(cells[4]), float(cells[6])))
    # Mean CM energies (the "C.M. spectra of emitted particles" table) and multiplicities.
    for cells in _rows(text):
        if cells and cells[0] == "Average Energy":
            out["mean_cm_energy_mev"] = dict(zip(("n", "p", "a"), (float(c) for c in cells[1:4])))
    mult = dict(re.findall(r"Mode = (NEUT|PROT|ALPH)\b.*?Multiplicity =\s*([\d.eE+-]+)", plain))
    out["multiplicity"] = {"n": float(mult.get("NEUT", 0)), "p": float(mult.get("PROT", 0)),
                           "a": float(mult.get("ALPH", 0))}
    return out


# --- writing the reference files ----------------------------------------------------------------------------------


def _header(version: str, system: str, settings: list, deck_text: str) -> list:
    lines = [f"# tool: PACE4 {version} in LISE++ (LISEcute), run in silent mode: PACE4.exe -s deck.pace",
             f"# produced by: scripts/pace4_reference.py (Claude Code), {date.today().isoformat()}",
             f"# system: {system}",
             f"# settings: {settings[0]}"]
    lines += [f"#   {s}" for s in settings[1:]]
    lines.append("# input deck (first three lines): " + " | ".join(ln.strip() for ln in deck_text.splitlines()))
    return lines


SETTINGS = ["PACE4 defaults as in its example deck; fusion cross section and partial waves from the Bass model",
            "with L-diffuseness 2; level-density parameter a = A/10; ground-state spins of beam and target 0"]


def _cn(beam, target) -> str:
    z, a = beam[1] + target[1], beam[2] + target[2]
    return f"Z = {z}, A = {a}"


def evaporation_file(pace4: Path, beam, target, e_lab: float) -> Path:
    text, out = run(pace4, deck(EVAPORATION_CASCADES, beam[1], beam[2], target[1], target[2], e_lab))
    p = parse(out)
    total = sum(r[2] for r in p["residues"])
    lines = _header(p["version"], f"{beam[0]} + {target[0]} at {e_lab:g} MeV (lab) -> compound nucleus "
                    f"{_cn(beam, target)}", SETTINGS[:1] + [SETTINGS[1] + f"; {EVAPORATION_CASCADES} cascades"], text)
    lines += [f"# excitation_mev: {p['excitation_mev']:g}",
              f"# sigma_fus_mb: {p['sigma_fus_mb']:g}",
              "# sigma_j_mb: " + ",".join(f"{s:g}" for s in p["sigma_j_mb"]),
              "# mean_cm_energy_mev: " + ",".join(f"{k}={v:g}" for k, v in p["mean_cm_energy_mev"].items()),
              "# multiplicity: " + ",".join(f"{k}={v:g}" for k, v in p["multiplicity"].items()),
              f"# cascades: {total}",
              "Z,A,events,sigma_mb"]
    lines += [f"{z},{a},{n},{s:g}" for z, a, n, s in sorted(p["residues"], key=lambda r: -r[2])]
    path = REF / f"pace4_{beam[0]}_{target[0]}_{e_lab:g}MeV.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def fusion_file(pace4: Path, beam, target, energies: list) -> Path:
    rows, version, first = [], None, None
    for e in energies:
        text, out = run(pace4, deck(FUSION_CASCADES, beam[1], beam[2], target[1], target[2], float(e)))
        first = first or text
        p = parse(out)
        version = p["version"]
        rows.append(f"{e:g},{p['e_cm_mev']:g},{p['excitation_mev']:g},{p['barrier_mev']:g},{p['fusion_radius_fm']:g},"
                    f"{p['sigma_fus_mb']:g},{p['l_grazing']:g},{';'.join(f'{s:g}' for s in p['sigma_j_mb'])}")
        print(f"  {beam[0]} + {target[0]} at {e} MeV: sigma_fus = {p['sigma_fus_mb']:g} mb", flush=True)
    lines = _header(version, f"{beam[0]} + {target[0]} -> compound nucleus {_cn(beam, target)}, excitation function",
                    SETTINGS[:1] + [SETTINGS[1] + f"; {FUSION_CASCADES} cascades per energy (σ_fus and σ_J do not "
                                    "depend on them)"], first)
    lines.append("# columns: beam energy (lab MeV), CM energy and E* (MeV) as PACE4 prints them, PACE4's Bass barrier "
                 "(CM MeV)")
    lines.append("#   and fusion radius (fm), σ_fus (mb), ℓ_grazing, σ_J in mb for J = 0, 1, 2, ... (semicolons)")
    lines.append("beam_energy_mev,e_cm_mev,excitation_mev,barrier_mev,fusion_radius_fm,sigma_fus_mb,l_grazing,"
                 "sigma_j_mb")
    lines += rows
    path = REF / f"pace4_fusion_{beam[0]}_{target[0]}.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pace4", type=Path, help="path to PACE4.exe")
    ap.add_argument("--only", choices=("evaporation", "fusion"))
    args = ap.parse_args(argv)
    for beam, target, e_lab, energies in SYSTEMS:  # strictly one PACE4 process at a time
        if args.only in (None, "evaporation"):
            print(f"{beam[0]} + {target[0]} at {e_lab:g} MeV, {EVAPORATION_CASCADES} cascades", flush=True)
            print("  wrote", evaporation_file(args.pace4, beam, target, e_lab).name, flush=True)
        if args.only in (None, "fusion"):
            print("  wrote", fusion_file(args.pace4, beam, target, energies).name, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
