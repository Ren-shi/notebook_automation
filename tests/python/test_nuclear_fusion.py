import math

import numpy as np
import pytest

from physim.nuclear import data
from physim.nuclear.fusion import Fusion
from physim.nuclear.quantity import Quantity
from physim.nuclear.stopping import Stopping

# Bass 1980 barrier heights V_B (MeV, CM) and positions R_B (fm): I. Dutt and R. K. Puri, Phys. Rev. C 81, 064609
# (2010), Table II, column "Bass 80" (arXiv:1005.5213).
BASS80 = [
    ("12C", "17O", 7.79, 8.13),
    ("4He", "209Bi", 20.30, 11.00),
    ("48Ca", "48Ca", 51.39, 10.40),
    ("16O", "144Sm", 61.34, 10.82),
    ("28Si", "92Zr", 71.21, 10.51),
    ("16O", "208Pb", 75.92, 11.60),
    ("40Ca", "90Zr", 99.28, 10.71),
    ("86Kr", "208Pb", 298.65, 13.15),
]


@pytest.fixture(scope="module")
def opb():
    return Fusion("16O", "208Pb")


@pytest.mark.parametrize("beam, target, vb, rb", BASS80)
def test_bass_barrier_matches_published_values(beam, target, vb, rb):
    b = Fusion(beam, target).barrier()
    assert b.height_mev == pytest.approx(vb, abs=0.03)
    assert b.radius_fm == pytest.approx(rb, abs=0.03)
    assert 2.0 < b.hbar_omega_mev < 6.0


def test_barrier_is_the_maximum_of_the_potential(opb):
    b = opb.barrier()
    r = np.linspace(b.radius_fm - 2, b.radius_fm + 2, 4001)
    v = opb.potential_mev(r)
    assert v.max() == pytest.approx(b.height_mev, abs=1e-6)
    assert r[np.argmax(v)] == pytest.approx(b.radius_fm, abs=2e-3)
    # The curvature against a finite difference of the potential.
    h = 1e-3
    d2 = (opb.potential_mev(b.radius_fm + h) - 2 * b.height_mev + opb.potential_mev(b.radius_fm - h)) / h**2
    assert b.hbar_omega_mev == pytest.approx(197.3269804 * math.sqrt(-d2 / opb.mu), rel=1e-5)


def test_l_dependent_barriers_rise_and_move_in(opb):
    bs = opb.barriers()
    assert len(bs) == opb.critical_l + 1
    v = np.array([b.height_mev for b in bs])
    r = np.array([b.radius_fm for b in bs])
    assert np.all(np.diff(v) > 0) and np.all(np.diff(r) < 0)
    # Low l: V_l ≈ V_0 + ħ²l(l+1)/(2μR_0²).
    l = 5
    expected = bs[0].height_mev + 197.3269804**2 * l * (l + 1) / (2 * opb.mu * bs[0].radius_fm**2)
    assert v[l] == pytest.approx(expected, rel=2e-3)


@pytest.mark.parametrize("factor", [1.5, 2.0, 3.0])
def test_classical_limit_far_above_the_barrier(opb, factor):
    """σ → πR_B²(1 − V_B/E) once E − V_B ≫ ħω, for Wong's formula and the fixed-barrier partial-wave sum."""
    b = opb.barrier()
    e = factor * b.height_mev
    geometric = math.pi * b.radius_fm**2 * (1 - b.height_mev / e) * 10.0
    assert opb.classical_cm(e) == pytest.approx(geometric, rel=1e-12)
    assert opb.wong_cm(e) == pytest.approx(geometric, rel=1e-6)
    assert opb.cross_section_cm(e, l_dependent=False) == pytest.approx(geometric, rel=2e-3)


def test_l_dependent_sum_saturates_at_the_critical_angular_momentum(opb):
    """Far above the barrier every wave up to l_crit fuses: σ → πƛ²(l_crit + 1)²."""
    e = 4 * opb.barrier().height_mev
    lam2 = 197.3269804**2 / (2 * opb.mu * e)
    assert opb.cross_section_cm(e) == pytest.approx(math.pi * lam2 * (opb.critical_l + 1) ** 2 * 10, rel=0.01)
    # Just above the barrier the l-dependence matters little.
    e = opb.barrier().height_mev + 3
    assert opb.cross_section_cm(e) == pytest.approx(opb.wong_cm(e), rel=0.02)


@pytest.mark.parametrize("system", [("16O", "208Pb"), ("48Ca", "48Ca"), ("28Si", "92Zr")])
def test_wong_closed_form_equals_the_partial_wave_sum(system):
    f = Fusion(*system)
    b = f.barrier()
    for de in (-6.0, -3.0, 0.0, 2.0, 5.0, 15.0, 40.0):
        e = b.height_mev + de
        total = f.partial_waves_cm(e, l_dependent=False).sum()
        assert total == pytest.approx(f.wong_cm(e), rel=0.01), de


def test_cross_section_falls_steeply_below_the_barrier(opb):
    b = opb.barrier()
    s = opb.cross_section_cm(b.height_mev - np.arange(0, 11, 1.0))
    assert np.all(np.diff(s) < 0)
    assert s[5] / s[0] < 1e-2 and s[10] / s[0] < 1e-5
    # Deep below, each MeV costs a factor exp(2π/ħω) (the parabola's tunnelling slope).
    assert s[9] / s[8] == pytest.approx(math.exp(-2 * math.pi / b.hbar_omega_mev), rel=0.1)
    assert opb.cross_section_cm(0.0) == 0.0


def test_thick_target_average_equals_the_integral_over_the_energy_loss(opb):
    e_in, t = 90.0, 2.0  # MeV, mg/cm²
    tt = opb.thick_target(e_in, t)
    stop = Stopping("16O", "208Pb")
    assert tt.e_out_mev == pytest.approx(stop.energy_after(e_in, t), rel=1e-9)
    # (1/t) ∫ σ(E) / S(E) dE from E_out to E_in, on a fine energy grid.
    e = np.linspace(tt.e_out_mev, e_in, 2001)
    integrand = opb.cross_section(e) / stop.stopping_power(e)
    expected = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(e)) / t
    assert tt.sigma_mb == pytest.approx(expected, rel=0.01)
    assert tt.sigma_l_mb.sum() == pytest.approx(tt.sigma_mb, rel=1e-12)
    assert opb.cross_section(tt.e_out_mev) < tt.sigma_mb < opb.cross_section(e_in)
    # A thin target gives the thin-target cross section; text thicknesses work.
    thin = opb.thick_target("90 MeV", "1 ug/cm2")
    assert thin.sigma_mb == pytest.approx(opb.cross_section(90.0), rel=1e-3)
    assert opb.excitation_function([90.0], thickness=t)[0] == pytest.approx(tt.sigma_mb)


def test_thick_target_where_the_beam_stops(opb):
    tt = opb.thick_target(90.0, 100.0)  # 100 mg/cm² of lead stops 90 MeV ¹⁶O
    assert tt.e_out_mev == 0.0
    assert tt.sigma_mb < opb.cross_section(90.0)
    assert any("stops in the target" in w for w in opb.warnings(90.0, thickness=100.0))


def test_compound_nucleus(opb):
    cn = opb.compound("90 MeV")
    assert cn.nuclide.name == "224Th"
    # E* = E_cm + Q with nuclear AME2020 masses.
    m1, m2 = data.nuclide("16O").nuclear_mass_mev, data.nuclide("208Pb").nuclear_mass_mev
    q = m1 + m2 - data.nuclide("224Th").nuclear_mass_mev
    e_cm = math.sqrt((m1 + m2) ** 2 + 2 * m2 * 90.0) - m1 - m2
    assert cn.e_cm_mev == pytest.approx(e_cm, abs=1e-9)
    assert cn.excitation_mev == pytest.approx(e_cm + q, abs=1e-6)
    assert cn.excitation_mev == pytest.approx(90.0 * m2 / (m1 + m2) + q, abs=0.05)  # non-relativistic: 17 keV off
    # The CN moves with the CM velocity.
    p = math.sqrt(90.0 * (90.0 + 2 * m1))
    assert cn.beta == pytest.approx(p / (90.0 + m1 + m2), rel=1e-12)
    # σ_l: indexed by l, sums to σ_fus, spin distribution normalised.
    assert cn.sigma_l_mb.ndim == 1 and cn.sigma_l_mb[-1] > 0
    assert cn.sigma_mb == pytest.approx(opb.cross_section(90.0), rel=1e-12)
    assert cn.spin_distribution().sum() == pytest.approx(1.0)
    assert 0 < cn.mean_l < cn.max_l <= opb.critical_l


def test_energies_in_any_form(opb):
    a = opb.cross_section(90.0)
    assert opb.cross_section("90 MeV") == a
    assert opb.cross_section(Quantity(90000, "keV")) == pytest.approx(a)
    arr = opb.cross_section(np.array([85.0, 90.0]))
    assert arr.shape == (2,) and arr[1] == pytest.approx(a)
    assert opb.beam_energy_for_cm(opb.e_cm(90.0)) == pytest.approx(90.0)
    assert opb.barrier_lab_mev == pytest.approx(opb.barrier().height_mev * 224.0 / 208.0, rel=1e-3)
    with pytest.raises(ValueError):
        opb.cross_section("90 mm")
    with pytest.raises(ValueError):
        Fusion("16O", "208Pb", potential="woods-saxon")


def test_bass77_is_available_and_close():
    f77, f80 = Fusion("16O", "208Pb", "bass77"), Fusion("16O", "208Pb")
    assert f77.barrier().height_mev != f80.barrier().height_mev
    assert f77.barrier().height_mev == pytest.approx(f80.barrier().height_mev, rel=0.03)


def test_warnings(opb):
    assert opb.warnings(100.0) == []
    below = opb.warnings(75.0)
    assert len(below) == 1 and "below the fusion barrier" in below[0]
    assert any("near the fusion barrier" in w for w in opb.warnings(84.0))
    assert any("incomplete fusion" in w for w in opb.warnings(200.0))
    assert any("critical angular momentum" in w for w in opb.warnings(250.0))
    heavy = Fusion("86Kr", "208Pb").warnings("450 MeV")
    assert any("quasi-fission" in w for w in heavy)
    assert not any("quasi-fission" in w for w in Fusion("40Ca", "90Zr").warnings("200 MeV"))
