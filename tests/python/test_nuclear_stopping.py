import csv
import math
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import Experiment, data
from physim.nuclear import stopping as st
from physim.nuclear.stopping import Stopping, angular_straggling_mrad
from physim.nuclear.validation import LISE_SYMBOLS, STOPPING_TOLERANCE

REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "nuclear"
LISE = REFERENCE / "stopping_lise.csv"


# -- the NIST tables --------------------------------------------------------------------------------------------


def test_tables_are_used_as_published():
    e, se, sn, _ = st._table("pstar", st._elements_in_tables()[79])
    s = Stopping("1H", "Au")
    pick = e[(e > 2e-4 * 1.007) & (e < 9e3)][::7]
    want = np.interp(pick, e, se)
    # The internal grid (801 points) adds < 0.05%, below the tables' 4 significant figures.
    np.testing.assert_allclose(s.stopping_power(pick, "electronic") * 1e3, want, rtol=5e-4)
    np.testing.assert_allclose(s.stopping_power(pick, "nuclear") * 1e3, np.interp(pick, e, sn), rtol=5e-3)
    assert sorted(st._elements_in_tables()) == [1, 2, 4, 6, 7, 8, 10, 13, 14, 18, 22, 26, 29, 32, 36, 42, 47, 50, 54,
                                                64, 74, 78, 79, 82, 92]
    assert "Graphite" in st._elements_in_tables()[6][1]
    assert st._star_compound("Mylar") is not None and st._star_compound("Water") is not None


@pytest.mark.parametrize("ion, prog", [("1H", "pstar"), ("4He", "astar")])
@pytest.mark.parametrize("material", ["C", "Si", "Cu", "Au", "Mylar", "Water", "Polyethylene"])
def test_integrated_range_matches_nist_csda_range(ion, prog, material):
    """Our range is ∫dE/S over our interpolated stopping; NIST tabulates its own CSDA range independently."""
    key = st._star_compound(material) or st._elements_in_tables()[data.element(material).Z]
    e, _, _, r = st._table(prog, key)
    s = Stopping(ion, material)
    energies = np.array([0.002, 0.01, 0.1, 1.0, 10.0, 100.0])
    nist = np.exp(np.interp(np.log(energies), np.log(e), np.log(r))) * 1e3
    np.testing.assert_allclose(s.range(energies), nist, rtol=0.01)


# -- range-energy consistency -----------------------------------------------------------------------------------


def test_energy_after_uses_the_range():
    s = Stopping("16O", "Au")
    e0, e1 = 64.0, 40.0
    path = s.range(e0) - s.range(e1)
    assert s.energy_after(e0, path) == pytest.approx(e1, rel=1e-6)
    assert s.energy_loss(e0, path) == pytest.approx(e0 - e1, rel=1e-6)
    assert s.energy_after(e0, s.range(e0) * 1.01) == 0.0
    # A tilted layer is crossed along a longer path.
    assert s.energy_after(e0, 1.0, tilt_deg=60.0) == pytest.approx(s.energy_after(e0, 2.0), rel=1e-9)
    # Thicknesses as lengths use the density: 1 um of gold is 1.932 mg/cm².
    assert s.energy_after(e0, "1 um") == pytest.approx(s.energy_after(e0, 1.932), rel=1e-9)
    energies = np.array([10.0, 30.0, 64.0])
    np.testing.assert_allclose(s.energy_after(energies, 1.0), [s.energy_after(x, 1.0) for x in energies])
    assert s.thickness_to_stop(e0) == s.range(e0)


def test_isotopes_scale_with_velocity():
    # Electronic stopping depends on velocity: a deuteron at 2E stops like a proton at E (per unit thickness).
    p, d = Stopping("1H", "Si"), Stopping("2H", "Si")
    for e in (0.5, 5.0, 50.0):
        ratio = d.mass_u / p.mass_u
        assert d.stopping_power(e * ratio, "electronic") == pytest.approx(p.stopping_power(e, "electronic"), rel=1e-3)
    he3, he4 = Stopping("3He", "Au"), Stopping("4He", "Au")
    assert he3.stopping_power(3.0 * 3.0, "electronic") == pytest.approx(he4.stopping_power(4.0 * 3.0, "electronic"),
                                                                       rel=5e-3)


def test_heavy_ions_are_fully_stripped_at_high_velocity():
    for ion, z in (("12C", 6), ("40Ar", 18)):
        assert st.effective_charge(z, 500.0) == pytest.approx(1.0, abs=2e-3)
        heavy, p = Stopping(ion, "Si"), Stopping("1H", "Si")
        e_u = 200.0
        expect = z**2 * p.stopping_power(e_u * p.mass_u, "electronic")
        assert heavy.stopping_power(e_u * heavy.mass_u, "electronic") == pytest.approx(expect, rel=1e-2)
    # Slow ions are partly neutralised.
    assert st.effective_charge(54, 0.1) < 0.5


def test_bragg_additivity_for_a_formula():
    cd2 = Stopping("1H", "CD2")
    # A 2 MeV proton is the tabulated particle at 2 MeV: per-electron stopping of C and H, summed by atoms.
    c, d_per_e = (st._per_electron("pstar", z, np.array([2.0]))[0] for z in (6, 1))
    m = data.material("CD2")
    expect = (6 * c + 2 * 1 * d_per_e) * data.AVOGADRO / m.molar_mass * 1e-3
    assert cd2.stopping_power(2.0, "electronic") == pytest.approx(expect, rel=5e-4)  # internal grid


@pytest.mark.parametrize("z", [13, 29, 47, 79, 82])
def test_interpolation_for_untabulated_elements(z):
    """Leave a tabulated element out and rebuild it from its neighbours: above 0.5 MeV within 7% (measured worst
    case 6.8%); near and below the Bragg peak the per-element shell structure is not captured (see the register)."""
    e = np.geomspace(0.5, 1000, 30)
    rebuilt = st._per_electron("pstar", z, e, exclude=frozenset({z}))
    np.testing.assert_allclose(rebuilt, st._per_electron("pstar", z, e), rtol=0.07)


def test_lindhard_scharff_low_velocity_stopping():
    # 10 keV protons in gold: 11.2 eV/(1e15 atoms/cm²) from LS; NIST gives 15.2 (LS is a ~30% estimate there).
    per_atom = st.lindhard_scharff(1, 79, 0.010 / 1.007276467)
    assert per_atom * 1e21 == pytest.approx(11.2, rel=0.02)


# -- straggling -------------------------------------------------------------------------------------------------


def test_bohr_straggling_thin_layer():
    """Thin layer of a fast proton: σ² = 4π r_e² (m c²)² N_A z² (Z/A) t × (1 − β²/2)/(1 − β²)."""
    s = Stopping("1H", "Si")
    t = 1.0  # mg/cm²
    e = 100.0
    beta2 = 1 - (1 / (1 + e / 938.272)) ** 2
    z_over_a = 14 / data.element("Si").molar_mass
    bohr = math.sqrt(0.1535e0 * 1.0 * z_over_a * t * 1e-3 * 4 * math.pi / (4 * math.pi) / 0.1535 * st.BOHR_K
                     * (1 - beta2 / 2) / (1 - beta2))
    assert s.straggling(e, t) == pytest.approx(bohr, rel=2e-3)
    assert st.BOHR_K == pytest.approx(0.1569, rel=1e-3)  # MeV² cm²/g, the textbook constant


def test_straggling_grows_through_a_thick_layer_above_the_bragg_peak():
    s = Stopping("4He", "Au")
    e = 20.0
    t = 0.5 * s.range(e)
    assert s.straggling(e, t) > s.straggling(e, t, propagate=False)
    assert s.straggling(e, s.range(e) * 1.1) == 0.0


def test_highland_multiple_scattering():
    # 10 MeV protons through 100 um of silicon, by hand from the Highland formula.
    mat = data.material("Si")
    x0 = st.radiation_length_g_cm2(mat)
    assert x0 == pytest.approx(21.82, rel=0.025)  # PDG: 21.82 g/cm²; Dahl's fit is good to 2.5%
    m, e = 938.272, 10.0
    p = math.sqrt(e * (e + 2 * m))
    beta = p / (e + m)
    x = 23.3e-3 / x0
    want = 1e3 * 13.6 / (beta * p) * math.sqrt(x) * (1 + 0.038 * math.log(x / beta**2))
    assert angular_straggling_mrad("1H", e, "Si", "100 um") == pytest.approx(want, rel=1e-3)


def test_inputs():
    s = Stopping("4He", "Au")
    with pytest.raises(ValueError, match="part must be"):
        s.stopping_power(5.0, "both")
    with pytest.raises(ValueError, match="density of CH2 is not known"):
        Stopping("1H", "CH2").energy_after(5.0, "10 um")
    assert Stopping("1H", "CH2", density="0.93 g/cm3").energy_after(5.0, "10 um") > 0
    exp = Experiment.example("alpha_on_gold")
    assert Stopping("4He", exp.target.material_data()).range(5.5) == pytest.approx(s.range(5.5))
    assert repr(s) == "Stopping(4He in Au)"


# -- comparison with LISE++ (scripts/lise_reference.py) ---------------------------------------------------------

SYM = LISE_SYMBOLS
#: Worst deviation from LISE++ (ATIMA 1.4) allowed, by quantity, ion class and energy band (MeV/u).
TOLERANCE = STOPPING_TOLERANCE


def _lise_rows():
    with open(LISE, encoding="utf-8") as f:
        return [r for r in csv.DictReader(line for line in f if not line.startswith("#")) if r["model"] == "ATIMA 1.4"]


@pytest.mark.skipif(not LISE.exists(), reason="no LISE++ reference file")
def test_against_lise_atima():
    worst = {}
    cache = {}
    for r in _lise_rows():
        z, a, zt, e_u = int(r["Z"]), int(r["A"]), int(r["Zt"]), float(r["e_mev_u"])
        if e_u < 0.1:
            continue  # below 0.1 MeV/u the codes themselves differ by tens of percent; not tested
        ion = data.nuclide((z, a)).name if z > 1 else "1H"
        s = cache.setdefault((ion, zt), Stopping(ion, SYM[zt]))
        cls = "p,a" if z <= 2 else ("C-Ar" if z <= 18 else "Kr,Xe")
        band = "0.1-1" if e_u < 1 else ("1-10" if e_u < 10 else ">=10")
        e = e_u * a
        got = {
            "stopping": s.stopping_power(e) / float(r["stopping_mev_mg_cm2"]),
            "range": s.range(e) / float(r["range_mg_cm2"]),
            "energy_after": s.energy_after(e, float(r["thickness_mg_cm2"])) / a / float(r["e_after_mev_u"]),
        }
        for q, ratio in got.items():
            key = (q, cls, band)
            worst[key] = max(worst.get(key, 0.0), abs(ratio - 1))
    failures = {k: round(v, 4) for k, v in worst.items() if v > TOLERANCE[k]}
    assert not failures, failures


@pytest.mark.skipif(not LISE.exists(), reason="no LISE++ reference file")
def test_straggling_against_lise_at_high_energy():
    """LISE++ reports the FWHM. Above 30 MeV/u physim agrees within 25%; below, the codes' slow-ion straggling
    theories differ (up to 2-4x at 1 MeV/u), so it is not tested there."""
    for r in _lise_rows():
        z, a, zt, e_u = int(r["Z"]), int(r["A"]), int(r["Zt"]), float(r["e_mev_u"])
        if e_u < 30 or zt not in (6, 79):
            continue
        ion = data.nuclide((z, a)).name if z > 1 else "1H"
        fwhm = 2.3548 * Stopping(ion, SYM[zt]).straggling(e_u * a, float(r["thickness_mg_cm2"]), steps=40) / a
        assert fwhm == pytest.approx(float(r["straggling_mev_u"]), rel=0.25), r
