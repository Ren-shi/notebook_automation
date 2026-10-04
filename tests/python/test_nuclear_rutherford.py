import math

import numpy as np
import pytest

from physim.nuclear.rutherford import E2_MEV_FM, Rutherford

# Geiger and Marsden, Phil. Mag. 25, 604 (1913), Table II, gold, first series: lab angle (deg) and number of
# scintillations N in equal times. Transcribed from https://www.chemteam.info/Chem-History/GeigerMarsden-1913/
# GeigerMarsden-1913.html; every row's N / (1/sin⁴(φ/2)) reproduces the paper's last column.
GEIGER_MARSDEN_GOLD = [(150, 33.1), (135, 43.0), (120, 51.9), (105, 69.5), (75, 211), (60, 477), (45, 1435),
                       (37.5, 3300), (30, 7800), (22.5, 27300), (15, 132000)]


def test_cross_section_formula():
    r = Rutherford("4He", "197Au", "20 MeV")
    k = 2 * 79 * E2_MEV_FM
    for th in (10.0, 50.0, 90.0, 179.0):
        expect = (k / (4 * r.e_cm)) ** 2 / math.sin(math.radians(th) / 2) ** 4 * 10
        assert r.cross_section_cm(th) == pytest.approx(expect, rel=1e-12)
    np.testing.assert_allclose(r.cross_section_cm(np.array([30.0, 60.0])), [r.cross_section_cm(30.0),
                                                                            r.cross_section_cm(60.0)])
    # 180°: σ = (k / 4E)² exactly, i.e. (d₀/4)².
    assert r.cross_section_cm(180.0) == pytest.approx((r.d0 / 4) ** 2 * 10, rel=1e-12)


def test_lab_cross_section_is_the_cm_one_times_the_jacobian():
    r = Rutherford("16O", "208Pb", "64 MeV")
    # Independent: count solid angle numerically. A ring dθ_lab maps to a CM ring dθ*; σ_lab = σ_cm dΩ*/dΩ_lab.
    for th_lab in (20.0, 60.0, 120.0, 170.0):
        h = 1e-4
        lo, hi = r.kinematics.at_lab(th_lab - h)[0], r.kinematics.at_lab(th_lab + h)[0]
        d_cos_cm = math.cos(math.radians(lo.theta_cm)) - math.cos(math.radians(hi.theta_cm))
        d_cos_lab = math.cos(math.radians(th_lab - h)) - math.cos(math.radians(th_lab + h))
        mid = r.kinematics.at_lab(th_lab)[0]
        expect = r.cross_section_cm(mid.theta_cm) * d_cos_cm / d_cos_lab
        first, second = r.cross_section_lab(th_lab)
        assert first == pytest.approx(expect, rel=1e-6) and second is None
    # The recoil uses its partner's CM angle: at 0° lab it is the ejectile at 180°, Jacobian (γ(1+g))² = 4 (NR).
    au = Rutherford("4He", "197Au", "20 MeV")
    rec, _ = au.cross_section_lab(1e-6, "recoil")
    assert rec == pytest.approx(au.cross_section_cm(180.0) * 4, rel=1e-3)


def test_inverse_kinematics_gives_two_lab_cross_sections():
    r = Rutherford("208Pb", "1H", "1000 MeV")
    first, second = r.cross_section_lab(0.2)
    assert first > 0 and second > 0


def test_integrated_cross_section():
    r = Rutherford("4He", "197Au", "20 MeV")
    th = np.linspace(20.0, 160.0, 200001)
    trapezoid = getattr(np, "trapezoid", None) or np.trapz  # NumPy 1.x (Python 3.9) has only trapz
    numeric = trapezoid(r.cross_section_cm(th) * np.sin(np.radians(th)), np.radians(th)) * 2 * math.pi
    assert r.integrated(20.0, 160.0) == pytest.approx(numeric, rel=1e-8)
    assert r.integrated(20.0, 160.0, dphi_deg=90.0) == pytest.approx(numeric / 4, rel=1e-8)
    with pytest.raises(ValueError, match="infinite"):
        r.integrated(0.0, 10.0)


def test_orbit_geometry():
    r = Rutherford("4He", "197Au", "20 MeV")
    assert r.d0 == pytest.approx(2 * 79 * E2_MEV_FM / r.e_cm)
    assert r.closest_approach(180.0) == pytest.approx(r.d0)
    for th in (10.0, 50.0, 120.0):
        b = r.impact_parameter(th)
        assert r.angle_for_impact(b) == pytest.approx(th, rel=1e-12)
        assert r.closest_approach_for_impact(b) == pytest.approx(r.closest_approach(th), rel=1e-12)
    # The values read off LISE++'s plots for this case (see the physics register): d₀ ≈ 11.6 fm, b(60°) ≈ 10.1 fm.
    assert r.d0 == pytest.approx(11.6, abs=0.05) and r.impact_parameter(60.0) == pytest.approx(10.05, abs=0.05)


@pytest.mark.parametrize("beam, target, energy", [("4He", "197Au", 20.0), ("16O", "208Pb", 64.0), ("1H", "12C", 2.0)])
def test_integrated_orbits_match_the_analytic_deflection(beam, target, energy):
    """physim's engine integrates the orbits (adaptive Dormand–Prince); the deflection is read from the
    asymptotes of the hyperbola through the first and last points and must equal 2 arctan(d₀/2b)."""
    r = Rutherford(beam, target, energy)
    bs = r.d0 * np.array([0.05, 0.3, 1.0, 3.0, 20.0])
    for orbit in r.trajectories(bs, distance=500 * r.d0):
        assert orbit.deflection() == pytest.approx(r.angle_for_impact(orbit.b), abs=1e-6)
        # The recorded orbit never comes closer than the analytic distance of closest approach.
        assert orbit.closest_approach() >= r.closest_approach_for_impact(orbit.b) * (1 - 1e-9)


def test_lise_value_for_alpha_on_gold():
    """LISE++'s kinematic calculator (user's screenshot, 2026-10-04): 4He + 197Au at 5 MeV/u, 50° CM,
    2640.611 mb/sr CM and 2711.374 mb/sr lab. LISE++ evaluated it with the beam slowed to mid-target
    (4.997 MeV/u); physim at the full 20 MeV is 0.03% below, at mid-target 0.07% above. Tolerance 0.1%."""
    full = Rutherford("4He", "197Au", 20.0)
    assert full.cross_section_cm(50.0) == pytest.approx(2640.611, rel=1e-3)
    lab, _ = full.cross_section_lab(full.kinematics.at_cm(50.0).theta_lab)
    assert lab == pytest.approx(2711.374, rel=1e-3)
    assert full.grazing_angle() == 180.0  # LISE++: grazing angle 180° (below the barrier)


def test_geiger_marsden_angular_distribution():
    """The 1913 gold counts follow 1/sin⁴(θ/2): data/theory stays within ±25% of its mean from 15° to 150° (the
    paper's own N sin⁴(φ/2) column runs from 27.5 to 39.6, ±18%), while a 1/sin²(θ/2) law would be off by a factor
    of 16 over the same range."""
    r = Rutherford("4He", "197Au", 7.7)  # RaC' α particles
    angles = np.array([a for a, _ in GEIGER_MARSDEN_GOLD], dtype=float)
    counts = np.array([n for _, n in GEIGER_MARSDEN_GOLD])
    theory = np.array([r.cross_section_lab(a)[0] for a in angles])
    ratio = counts / theory
    assert np.all(np.abs(ratio / ratio.mean() - 1) < 0.25), ratio / ratio.mean()
    assert ratio.max() / ratio.min() < 1.5
    wrong = counts * np.sin(np.radians(angles) / 2) ** 2
    assert wrong.max() / wrong.min() > 16


def test_validity_checks():
    below = Rutherford("4He", "197Au", 20.0)
    assert below.grazing_angle() == 180.0 and below.warnings() == []
    assert below.interaction_radius == pytest.approx(1.2 * (4 ** (1 / 3) + 197 ** (1 / 3)) + 2.0)
    gamma = 1 + 20.0 / below.kinematics.m1
    assert below.sommerfeld == pytest.approx(2 * 79 / 137.035999 / math.sqrt(1 - 1 / gamma**2), rel=1e-6)
    assert below.coulomb_barrier_cm == pytest.approx(2 * 79 * E2_MEV_FM / below.interaction_radius)
    above = Rutherford("4He", "197Au", 40.0)
    g = above.grazing_angle()
    assert 0 < g < 180
    assert above.closest_approach(g) == pytest.approx(above.interaction_radius, rel=1e-12)
    assert any("not Rutherford" in w for w in above.warnings())
    assert above.warnings(theta_cm_max=g * 0.9) == []
    slow = Rutherford("4He", "197Au", 0.1)
    assert any("screening" in w for w in slow.warnings())
    assert slow.screening_correction() == pytest.approx(1 - 0.049 * 2 * 79 ** (4 / 3) / (slow.e_cm * 1e3))
    assert any("Mott" in w for w in Rutherford("12C", "12C", 20.0).warnings())
    light = Rutherford("1H", "1H", 30.0)
    assert any("Sommerfeld" in w for w in light.warnings())
