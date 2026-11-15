"""The evaporation Monte Carlo (backlog item 75, Stage 3): each chain conserves energy and momentum exactly; with
only neutrons allowed the first neutron follows Weisskopf's spectrum ε σ_inv(ε) exp(−ε/T); seeded runs repeat bit
for bit; the channel shares sum to one; fission takes over in heavy nuclei; the xn thresholds of 4He + 165Ho fall
where they are measured; and the channels and particle energies of 12C + 150Nd at 70 MeV come close to PACE4's."""

import math
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import data
from physim.nuclear import evaporation as ev

REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "nuclear"


def nuclear_mass(z, a):
    return data.nuclide((z, a)).nuclear_mass_mev


@pytest.fixture(scope="module")
def run():
    """A medium system with all decays open: 162Dy at 59 MeV with a fusion-like spin distribution."""
    sigma_l = np.array([(2 * l + 1) / (1 + math.exp((l - 30) / 2.0)) for l in range(50)])
    return ev.evaporate("162Dy", 59.3, sigma_l_mb=sigma_l, events=4000, seed=7)


# ---------------------------------------------------------------------------------------------------------------
# Conservation and bookkeeping


def test_energy_and_momentum_are_conserved_in_every_event(run):
    z, a = run.compound
    total = nuclear_mass(z, a) + run.excitation_mev
    m_res = np.array([nuclear_mass(zz, aa) for zz, aa in zip(run.residue_z, run.residue_a)])
    e_res = np.sqrt((run.recoil_momentum_mev_c**2).sum(axis=1) + (m_res + run.entry_energy_mev) ** 2)
    e_emit = np.bincount(run.event, weights=run.energy_mev + run.mass_mev, minlength=run.events)
    assert np.abs(e_emit + e_res - total).max() < 1e-4  # MeV, of 1.5e5 MeV
    p_sum = np.stack([np.bincount(run.event, weights=run.momentum_mev_c[:, i], minlength=run.events)
                      for i in range(3)], axis=1)
    assert np.abs(p_sum + run.recoil_momentum_mev_c).max() < 1e-6  # MeV/c


def test_residues_match_the_emitted_particles(run):
    z, a = run.compound
    for i in range(200):
        kinds = run.particles(i)["kind"]
        n, p, al = kinds.count("n"), kinds.count("p"), kinds.count("a")
        if run.fission[i]:
            continue
        assert run.residue_z[i] == z - p - 2 * al
        assert run.residue_a[i] == a - n - p - 4 * al
        assert run.label_of_event[i] == ev.channel_label(n, p, al)


def test_entry_states_are_physical(run):
    ok = ~run.fission
    assert (run.entry_spin[ok] >= 0).all()
    # Above the yrast line of the residue and below where another neutron could leave easily.
    yrast = np.array([ev.yrast_energy_mev(aa, j) for aa, j in zip(run.residue_a[ok], run.entry_spin[ok])])
    assert (run.entry_energy_mev[ok] >= yrast - 1e-9).all()
    assert run.entry_energy_mev[ok].mean() < 20.0
    # The rotational part of the entry energy is that yrast energy.
    assert np.allclose(run.entry_rotational_energy_mev[ok], yrast, rtol=1e-12)
    assert np.isnan(run.entry_rotational_energy_mev[run.fission]).all()


def test_shares_sum_to_one_and_cross_sections_to_sigma_fus(run):
    assert sum(run.shares().values()) == pytest.approx(1.0, abs=1e-12)
    table = run.channel_table()
    assert sum(r["share"] for r in table) == pytest.approx(1.0, abs=1e-12)
    assert sum(r["sigma_mb"] for r in table) == pytest.approx(run.sigma_fusion_mb, rel=1e-12)
    assert table[0]["share"] >= table[-1]["share"]


def test_seeded_runs_are_bit_identical():
    kw = dict(spin_weights=np.arange(40.0), events=1500)
    a = ev.evaporate("150Gd", 45.0, seed=11, **kw)
    b = ev.evaporate("150Gd", 45.0, seed=11, **kw)
    c = ev.evaporate("150Gd", 45.0, seed=12, **kw)
    for name in ("initial_spin", "residue_z", "residue_a", "entry_energy_mev", "entry_spin",
                 "recoil_momentum_mev_c", "offsets", "kind", "energy_mev", "momentum_mev_c"):
        assert np.array_equal(getattr(a, name), getattr(b, name)), name
    assert a.labels == b.labels
    assert not np.array_equal(a.energy_mev[:100], c.energy_mev[:100])


def test_inputs():
    with pytest.raises(ValueError):
        ev.evaporate("162Dy", 40.0, events=10, seed=1)  # no spin distribution
    with pytest.raises(ValueError):
        ev.evaporate("162Dy", 40.0, spin=3, spin_weights=[1, 1], events=10, seed=1)
    with pytest.raises(ValueError):
        ev.evaporate("162Dy", 40.0, spin=3, events=10, seed=1, decays=("n", "x"))
    r = ev.evaporate((66, 162), 40.0, spin=5, events=50, seed=1)
    assert r.sigma_fusion_mb is None and (r.initial_spin == 5).all()
    # A spin above the yrast line at this energy is lowered to the line.
    r = ev.evaporate("162Dy", 10.0, spin=200, events=20, seed=1)
    assert r.spin_clamped == 20 and r.initial_spin.max() < 50


def test_lab_boost_and_events_list(run):
    beta = 0.03
    e_lab, p_lab = run.lab(beta)
    m = run.mass_mev
    # Invariant mass is unchanged by the boost.
    inv_cm = (run.energy_mev + m) ** 2 - (run.momentum_mev_c**2).sum(axis=1)
    inv_lab = (e_lab + m) ** 2 - (p_lab**2).sum(axis=1)
    assert np.allclose(inv_cm, inv_lab, rtol=1e-9, atol=1e-6)
    # Forward-going particles gain energy on average.
    assert e_lab.mean() > run.energy_mev.mean()
    events = run.events_list(beta)
    assert len(events) == run.events
    first = events[0]["particles"]
    assert len(first) == run.offsets[1] - run.offsets[0]
    kind, energy, direction = first[0]
    assert kind in ev.KINDS and energy > 0 and np.linalg.norm(direction) == pytest.approx(1.0)


# ---------------------------------------------------------------------------------------------------------------
# Physics


def _first_neutrons(z, a, e, lo, hi):
    """Histogram of the first neutron's energy (only neutrons allowed, spin 0), divided by ε σ_inv(ε)."""
    r = ev.evaporate((z, a), e, spin=0, events=40000, seed=3, decays=("n",))
    assert set(np.unique(r.kind)) == {0}
    first = r.energy_mev[r.offsets[:-1]]
    edges = np.arange(lo, hi + 1e-9, 0.25)
    counts, _ = np.histogram(first, edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    # ε σ_inv(ε) ∝ Σ (2ℓ+1) T_ℓ(ε), up to a constant.
    t = ev.neutron_transmission(mid, a - 1, 20)
    g = ((2 * np.arange(21) + 1)[:, None] * t).sum(axis=0)
    return mid, np.log(counts / g), counts  # counts are the weights 1/σ² of ln N


def test_neutron_spectrum_is_weisskopf():
    """Only neutrons, spin 0: the first neutron's spectrum is Weisskopf's ε σ_inv(ε) ρ(U − ε), and over a few MeV
    the single-temperature form ε σ_inv(ε) exp(−ε/T)."""
    z, a, e = 70, 168, 30.0
    sep = nuclear_mass(z, a - 1) + nuclear_mass(0, 1) - nuclear_mass(z, a)
    # The full form, 0.5–8 MeV: ln N − ln(ε σ_inv) − ln ρ is flat.
    mid, y, w = _first_neutrons(z, a, e, 0.5, 8.0)
    lnrho = np.log(ev.level_density(z, a - 1, e - sep - mid, 0))
    c = np.average(y - lnrho, weights=w)
    assert (w * (y - lnrho - c) ** 2).sum() / (mid.size - 1) < 2.0
    # One temperature, 0.5–4 MeV.
    mid, y, w = _first_neutrons(z, a, e, 0.5, 4.0)
    slope, intercept = np.polyfit(mid, y, 1, w=np.sqrt(w))
    assert (w * (y - intercept - slope * mid) ** 2).sum() / (mid.size - 2) < 2.0
    t_fit = -1 / slope
    lnrho = np.log(ev.level_density(z, a - 1, e - sep - mid, 0))
    assert t_fit == pytest.approx(-1 / np.polyfit(mid, lnrho, 1)[0], rel=0.05)
    # The Fermi-gas temperature √(U/a), U at the top of the spectrum less the pairing shift.
    u_top = e - sep - ev.pairing_shift_mev(z, a - 1)
    assert t_fit == pytest.approx(math.sqrt(u_top / ev.level_density_parameter(a - 1)), rel=0.15)


def test_fission_competes_in_heavy_nuclei():
    heavy = ev.evaporate("224Th", 50.0, spin=0, events=1500, seed=5)
    spun = ev.evaporate("224Th", 50.0, spin=20, events=1500, seed=5)
    light = ev.evaporate("172Yb", 50.0, spin=20, events=500, seed=5)
    assert 0.2 < heavy.fission.mean() < 0.9
    assert spun.fission.mean() > heavy.fission.mean() + 0.1  # the barrier falls with spin
    assert light.fission.mean() == 0.0
    assert "fission" in heavy.labels and heavy.channel_table()[0]["residue"] in ("fission", "219Th", "220Th")


def test_fission_barrier_is_liquid_drop():
    # Sierk's barrier for 172Yb at L = 0 is 30.04 MeV (as printed by PACE4); the liquid drop gives 30.3 MeV.
    assert ev.fission_barrier_mev(70, 172) == pytest.approx(30.04, rel=0.03)
    assert ev.fission_barrier_mev(90, 224) < 8.0
    assert ev.fission_barrier_mev(98, 252) < ev.fission_barrier_mev(90, 224)


def test_transmission_coefficients():
    eps = np.array([0.1, 1.0, 5.0, 20.0])
    tn = ev.neutron_transmission(eps, 160, 6)
    assert ((tn >= 0) & (tn <= 1)).all()
    assert tn[0, 0] < tn[0, 1] and tn[0, 0] == pytest.approx(4 * 0.069 / 1.0, rel=0.25)  # 4k/K at low energy
    # Charged particles: ½ at the barrier top, tiny far below it, near 1 well above.
    vb, hw, rb = ev.barriers("a", 64, 156, 0)
    t = ev.charged_transmission("a", 64, 156, [vb[0] - 8, vb[0], vb[0] + 8], 0)[0]
    assert t[1] == pytest.approx(0.5, abs=0.01)
    assert t[0] < 1e-6 and t[2] > 0.99
    assert 15.0 < vb[0] < 21.0 and 9.0 < rb[0] < 12.0  # α on Gd: about 17 MeV at 10.5 fm


def test_gdr_strength_gives_radiative_widths_of_the_right_size():
    """Γ_γ at the neutron separation energy of 168Yb from the E1 strength: about 0.1 eV, as measured for
    rare-earth nuclei (RIPL: 50–100 meV)."""
    e, j = 9.06, 1
    eps = np.linspace(1e-3, e, 3000)
    t = ev.e1_transmission(eps, 70, 168)
    w = sum((t * ev.level_density(70, 168, e - eps, jp)).sum() * (eps[1] - eps[0]) for jp in (0, 1, 2))
    gamma_ev = w / (2 * math.pi * ev.level_density(70, 168, e, j)) * 1e6
    assert 0.03 < gamma_ev < 0.4


# ---------------------------------------------------------------------------------------------------------------
# Measured data


def _crossing(es, low, high):
    """First energy where ``high`` overtakes ``low`` (linear interpolation)."""
    d = np.asarray(low, float) - np.asarray(high, float)
    i = int(np.nonzero((d[:-1] > 0) & (d[1:] <= 0))[0][0])
    return es[i] + (es[i + 1] - es[i]) * d[i] / (d[i] - d[i + 1])


def test_xn_thresholds_match_measured_4he_165ho():
    """⁴He + ¹⁶⁵Ho: the beam energies where 3n overtakes 2n and 4n overtakes 3n, from physim's own fusion σ_ℓ,
    against three stacked-foil measurements (Usman 2020, Kormazeva 2022, Gadkari 1997; EXFOR).

    Each crossing sits where the energy left after the neutrons can just pay for one more, so it tests the energy the
    neutrons carry, which the absolute cross sections (normalised by σ_fus, with pre-equilibrium losses above 45 MeV)
    do not. Measured: 2n/3n 30.5–31.8 MeV, 3n/4n 41.8–42.2 MeV. Tolerance 1 MeV from the mean of the three, the
    foils' own energy uncertainty (0.8–2.4 MeV near the crossings). PACE4 4.36 puts 3n/4n at 43.9 MeV, its harder
    neutrons (docs/theory/evaporation.md).
    """
    from physim.nuclear.fusion import Fusion

    rows = [ln.split(",") for ln in (REFERENCE / "evaporation_4He_165Ho_literature.csv").read_text(
        encoding="utf-8").splitlines() if ln and not ln.startswith(("#", "beam"))]
    fusion = Fusion("4He", "165Ho")
    for low, high, energies in (("2n", "3n", (29.0, 30.0, 31.0, 32.0, 33.0)),
                                ("3n", "4n", (40.0, 41.0, 42.0, 43.0, 44.0))):
        measured = []
        for ref in sorted({r[8] for r in rows}):
            pts: dict = {}
            for r in rows:
                if r[8] == ref and r[3] in (low, high):
                    pts.setdefault(float(r[4]), {})[r[3]] = float(r[6])
            es = sorted(e for e, v in pts.items() if low in v and high in v)
            measured.append(_crossing(es, [pts[e][low] for e in es], [pts[e][high] for e in es]))
        shares = {low: [], high: []}
        for e in energies:
            cn = fusion.compound(e)
            r = ev.evaporate("169Tm", cn.excitation_mev, sigma_l_mb=cn.sigma_l_mb, events=10000, seed=5)
            s = r.shares()
            shares[low].append(s.get(low, 0.0))
            shares[high].append(s.get(high, 0.0))
        ours = _crossing(energies, shares[low], shares[high])
        assert ours == pytest.approx(np.mean(measured), abs=1.0), (low, high, ours, measured)


# ---------------------------------------------------------------------------------------------------------------
# PACE4


def _pace(name):
    path = REFERENCE / name
    header, rows = {}, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, _, value = line[1:].partition(":")
            header[key.strip()] = value.strip()
        elif line and not line.startswith("Z,"):
            z, a, n, s = line.split(",")
            rows.append((int(z), int(a), int(n), float(s)))
    return header, rows


def test_channels_and_energies_close_to_pace4():
    """12C + 150Nd at 70 MeV (E* = 59.3 MeV in 162Dy) with PACE4's own spin distribution.

    Tolerances: the strongest residue the same and within 20%, the second within 30% and in order; mean proton and α
    energies within 10%. Mean neutron energies are 10–25% below PACE4's, within 25% here: the spectra share their
    slope above 4 MeV, but PACE4 has only 6–11% of its neutrons below 1 MeV (physim 17–25%), and the measured xn
    thresholds of 4He + 165Ho side with physim (docs/theory/evaporation.md, "Inside PACE4"). The softer neutrons
    leave room for a sixth neutron in more chains, so the third channel (6n) is not compared.
    """
    header, rows = _pace("pace4_12C_150Nd_70MeV.csv")
    sigma_j = [float(x) for x in header["sigma_j_mb"].split(",")]
    e_star = float(header["excitation_mev"])
    mean = dict(kv.split("=") for kv in header["mean_cm_energy_mev"].split(","))
    total = sum(r[2] for r in rows)
    pace = sorted(((n / total, (z, a)) for z, a, n, _ in rows), reverse=True)

    r = ev.evaporate("162Dy", e_star, sigma_l_mb=sigma_j, events=20000, seed=1)
    assert r.sigma_fusion_mb == pytest.approx(sum(sigma_j))
    mine = {(row["Z"], row["A"]): row["share"] for row in r.channel_table() if row["Z"] is not None}
    ranked = sorted(mine, key=lambda k: -mine[k])
    (s1, k1), (s2, k2) = pace[0], pace[1]
    assert ranked[0] == k1
    assert mine[k1] == pytest.approx(s1, rel=0.20)
    assert mine[k2] == pytest.approx(s2, rel=0.30)
    assert mine[k1] > mine[k2]
    assert r.mean_energy_mev("p") == pytest.approx(float(mean["p"]), rel=0.10)
    assert r.mean_energy_mev("a") == pytest.approx(float(mean["a"]), rel=0.10)
    assert r.mean_energy_mev("n") == pytest.approx(float(mean["n"]), rel=0.25)
