"""What PACE4 does inside its cascades: the numbers behind "PACE4's internals" in ``docs/theory/evaporation.md``.

    python scripts/pace4_internals.py "D:/Program Files/LISEcute/LISEcute/PACE4.exe" [--physim]

Runs only on Windows with PACE4 (installed with LISE++), strictly one PACE4 process at a time, each deck in a fresh
temporary folder (as ``scripts/pace4_reference.py``, whose decks and HTML reader it reuses; nothing is written to the
repository). For the five systems of the theory note, 20 000 cascades each, it reads from PACE4's own decay
summaries:

- the spin J of the compound nucleus at its first decay (the top excitation-energy row of the neutron, proton and α
  tables), against the mean of the σ_J that PACE4 prints;
- the angular momentum its cascades remove: Σ multiplicity × "average spin removed" over n, p, α, E1 and E2, plus the
  spin left in the cascades it reports as trapped;
- its CM neutron spectrum in 1 MeV bins: the share below 1 MeV, the mean, and the slope temperature of the tail
  (4–10 MeV, N/ε ∝ exp(−ε/T));
- its electromagnetic strengths (Weisskopf units).

Then the ⁴He + ¹⁶⁵Ho excitation function, 26–50 MeV in 2 MeV steps: the beam energies where 3n overtakes 2n and 4n
overtakes 3n, against the three measured data sets of ``evaporation_4He_165Ho_literature.csv``.

``--physim`` adds physim's numbers for the same E* and σ_J (20 000 events, seed 1; the excitation function seed 5).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pace4_reference as pr  # noqa: E402

LITERATURE = pr.REF / "evaporation_4He_165Ho_literature.csv"
HO_ENERGIES = list(range(26, 51, 2))
CASCADES = 20000


def _plain(html_text: str) -> str:
    """PACE4's HTML as text: one table row per line, cells separated by tabs."""
    t = re.sub(r"</tr>", "\n", html_text)
    t = re.sub(r"<br\s*/?>|<p>", "\n", t)
    t = re.sub(r"</t[dh]>", "\t", t)
    import html
    return re.sub(r"\n\s*\n+", "\n", html.unescape(re.sub(r"<[^>]+>", "", t)))


_ROW = re.compile(r"(\d+)\s*-\s*(\d+)\t((?:\s*(?:\d+|\.\.)\s*\t){16})\s*(\d+)\t")


def decay_summaries(text: str) -> dict:
    """{mode: {"multiplicity", "spin_removed", "rows": [(Ex_lo, counts per 5-ħ J bin)]}} for NEUT, PROT, ALPH, G-E1,
    G-E2."""
    out = {}
    parts = re.split(r"Decay summary\.\s+Mode = (\S+)", text)
    for name, block in zip(parts[1::2], parts[2::2]):
        mult = float(re.search(r"Multiplicity =\s*([\d.eE+-]+)", block).group(1))
        removed = float(re.search(r"Average spin removed =\s*([\d.eE+-]+)", block).group(1))
        rows = []
        for m in _ROW.finditer(block.split("Ex / J", 1)[0]):
            rows.append((int(m.group(1)), [0 if c.strip() == ".." else int(c)
                                           for c in m.group(3).split("\t") if c.strip()]))
        out[name] = {"multiplicity": mult, "spin_removed": removed, "rows": rows}
    return out


def neutron_cm_spectrum(text: str) -> np.ndarray:
    """Counts of the CM neutron spectrum in 1 MeV bins from 0 MeV."""
    block = text.split("C.M. spectra  of  emitted  particles", 1)[1].split("Average Energy", 1)[0]
    counts = np.zeros(40)
    for m in re.finditer(r"^\s*(\d+)\s*-\s*(\d+)\t([^\n]*)", block, re.M):
        v = m.group(3).split("\t")[0].strip()
        counts[int(m.group(1))] = int(v) if v else 0
    return counts


def tail_temperature(counts: np.ndarray, lo: int = 4, hi: int = 10) -> float:
    """Slope temperature of N(ε)/ε over the 1 MeV bins lo..hi."""
    c = np.arange(lo, hi) + 0.5
    return float(-1.0 / np.polyfit(c, np.log(counts[lo:hi] / c), 1)[0])


def internals(html_text: str) -> dict:
    text = _plain(html_text)
    p = pr.parse(html_text)
    ds = decay_summaries(text)
    first = np.zeros(16)
    for mode in ("NEUT", "PROT", "ALPH"):
        if ds.get(mode, {}).get("rows"):
            first += np.array(max(ds[mode]["rows"], key=lambda r: r[0])[1], float)
    centres = 5 * np.arange(16) + 2.0
    sj = np.array(p["sigma_j_mb"])
    trapped = re.search(r"([\d.]+)\s+percent of cascades trapped.*?average spin =\s*([\d.]+)", text, re.S)
    removed = {m: d["multiplicity"] * d["spin_removed"] for m, d in ds.items()}
    left = float(trapped.group(1)) / 100 * float(trapped.group(2)) if trapped else 0.0
    spec = neutron_cm_spectrum(text)
    strengths = dict(re.findall(r"\b(E1|M1|E2|M2) = ([\d.]+)", text))
    return {"parsed": p, "first_j_bins": first, "first_mean_j": float((first * centres).sum() / first.sum()),
            "sigma_j_mean": float((np.arange(sj.size) * sj).sum() / sj.sum()), "removed": removed,
            "trapped_spin": left, "budget": sum(removed.values()) + left, "spectrum": spec,
            "below_1mev": float(spec[0] / spec.sum()), "tail_t": tail_temperature(spec), "strengths": strengths}


def crossing(es, low, high) -> float:
    """First energy where ``high`` overtakes ``low`` (linear interpolation), nan if it never does."""
    d = np.asarray(low, float) - np.asarray(high, float)
    for i in range(len(d) - 1):
        if np.isfinite(d[i]) and np.isfinite(d[i + 1]) and d[i] > 0 >= d[i + 1]:
            return float(es[i] + (es[i + 1] - es[i]) * d[i] / (d[i] - d[i + 1]))
    return float("nan")


def measured_crossings() -> dict:
    """{("2n", "3n"): [per data set], ("3n", "4n"): [...]}, each data set interpolated on its own energies."""
    rows = [dict(zip(("beam", "target", "residue", "channel", "e", "de", "s", "ds", "ref", "sub"), ln.split(",")))
            for ln in LITERATURE.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith(("#", "beam"))]
    out = {}
    for lo, hi in (("2n", "3n"), ("3n", "4n")):
        vals = []
        for ref in sorted({r["ref"] for r in rows}):
            pts: dict = {}
            for r in rows:
                if r["ref"] == ref and r["channel"] in (lo, hi):
                    pts.setdefault(float(r["e"]), {})[r["channel"]] = float(r["s"])
            es = sorted(k for k, v in pts.items() if lo in v and hi in v)
            vals.append((ref, crossing(es, [pts[k][lo] for k in es], [pts[k][hi] for k in es])))
        out[(lo, hi)] = vals
    return out


def _physim(z, a, e_star, sigma_j, seed):
    from physim.nuclear import evaporation as ev
    return ev.evaporate((z, a), e_star, sigma_l_mb=sigma_j, events=CASCADES, seed=seed)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pace4", type=Path, help="path to PACE4.exe")
    ap.add_argument("--physim", action="store_true", help="also run physim on the same E* and σ_J")
    args = ap.parse_args(argv)
    print("system | <J> sigma_J / first decay | removed: n p a E1 E2 trapped = total | n below 1 MeV, "
          "<e_n>, tail T | strengths")
    for beam, target, e_lab, _ in pr.SYSTEMS:  # strictly one PACE4 process at a time
        _, out = pr.run(args.pace4, pr.deck(CASCADES, beam[1], beam[2], target[1], target[2], e_lab))
        d = internals(out)
        r = d["removed"]
        print(f"{beam[0]}+{target[0]} | {d['sigma_j_mean']:.1f} / {d['first_mean_j']:.1f} | "
              + " ".join(f"{r.get(m, 0):.1f}" for m in ("NEUT", "PROT", "ALPH", "G-E1", "G-E2"))
              + f" {d['trapped_spin']:.1f} = {d['budget']:.1f} | {d['below_1mev']:.1%}, "
              f"{d['parsed']['mean_cm_energy_mev']['n']:.2f}, {d['tail_t']:.2f} MeV | {d['strengths']}", flush=True)
        if args.physim:
            x = _physim(beam[1] + target[1], beam[2] + target[2], d["parsed"]["excitation_mev"],
                        d["parsed"]["sigma_j_mb"], 1)
            en = x.energy_mev[x.kind == 0]
            h = np.histogram(en, np.arange(0, 41))[0].astype(float)
            ok = ~x.fission
            print(f"   physim: <J> {x.initial_spin.mean():.1f}, entry <J> {x.entry_spin[ok].mean():.1f}, n below 1 MeV "
                  f"{h[0] / h.sum():.1%}, <e_n> {en.mean():.2f}, tail T {tail_temperature(h):.2f} MeV", flush=True)
    print("\n4He + 165Ho: shares of 2n 3n 4n, PACE4" + (" | physim" if args.physim else ""))
    shares = {"pace4": {c: [] for c in ("2n", "3n", "4n")}, "physim": {c: [] for c in ("2n", "3n", "4n")}}
    za = {"2n": (69, 167), "3n": (69, 166), "4n": (69, 165)}
    for e in HO_ENERGIES:
        _, out = pr.run(args.pace4, pr.deck(CASCADES, 2, 4, 67, 165, float(e)))
        p = pr.parse(out)
        tot = sum(row[2] for row in p["residues"])
        sh = {(z, a): n / tot for z, a, n, _ in p["residues"]}
        line = f"  {e} MeV: " + " ".join(f"{sh.get(za[c], 0):.3f}" for c in za)
        for c in za:
            shares["pace4"][c].append(sh.get(za[c], 0.0))
        if args.physim:
            x = _physim(69, 169, p["excitation_mev"], p["sigma_j_mb"], 5)
            mine: dict = {}
            for row in x.channel_table():
                if row["Z"] is not None:
                    mine[(row["Z"], row["A"])] = mine.get((row["Z"], row["A"]), 0) + row["share"]
            for c in za:
                shares["physim"][c].append(mine.get(za[c], 0.0))
            line += " | " + " ".join(f"{mine.get(za[c], 0):.3f}" for c in za)
        print(line, flush=True)
    meas = measured_crossings()
    print("\ncrossing (lab MeV): measured per data set | PACE4" + (" | physim" if args.physim else ""))
    for lo, hi in (("2n", "3n"), ("3n", "4n")):
        line = f"  {lo}/{hi}: " + ", ".join(f"{ref} {v:.1f}" for ref, v in meas[(lo, hi)])
        line += f" | {crossing(HO_ENERGIES, shares['pace4'][lo], shares['pace4'][hi]):.1f}"
        if args.physim:
            line += f" | {crossing(HO_ENERGIES, shares['physim'][lo], shares['physim'][hi]):.1f}"
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
