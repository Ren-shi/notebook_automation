import csv
import math
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import Experiment, data
from physim.nuclear.kinematics import TwoBody, elastic, kinematic_factor, lab_points

REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "nuclear"

CASES = [  # (beam, target, energy MeV, ejectile, excitation)
    ("4He", "197Au", 5.5, None, 0.0),
    ("16O", "208Pb", 64.0, None, 0.0),
    ("16O", "208Pb", 64.0, None, 2.614),  # inelastic: 208Pb 3- state
    ("208Pb", "1H", 1000.0, None, 0.0),  # inverse kinematics
    ("2H", "12C", 10.0, "1H", 0.0),  # 12C(d,p)13C
    ("1H", "3H", 1.05, "n", 0.0),  # just above threshold: double-valued
    ("12C", "12C", 2000.0, None, 0.0),  # equal masses, relativistic (167 MeV/u)
]


def reactions():
    for beam, target, e, ejectile, ex in CASES:
        yield TwoBody(beam, target, e, ejectile=ejectile, excitation_mev=ex)


def lab_from_four_vectors(r, theta_cm_deg, particle="ejectile"):
    """Independent route: build the CM four-momentum and apply an explicit Lorentz boost along z."""
    m = r.m3 if particle == "ejectile" else r.m4
    th = math.radians(theta_cm_deg)
    p = r.p_cm
    e = math.hypot(p, m)
    px, pz = p * math.sin(th), p * math.cos(th)
    b, g = r.beta_cm, r.gamma_cm
    pz_lab = g * (pz + b * e)
    e_lab = g * (e + b * pz)
    return math.degrees(math.atan2(px, pz_lab)), e_lab - m


# -- independent checks -----------------------------------------------------------------------------------------


@pytest.mark.parametrize("r", list(reactions()), ids=str)
def test_matches_explicit_lorentz_boost(r):
    # Not 180 deg: in elastic scattering the recoil is then at rest in the lab, with no direction.
    for th in (0.0, 10.0, 45.0, 90.0, 137.0, 179.0):
        for particle in ("ejectile", "recoil"):
            point = r.at_cm(th, particle)
            angle, energy = lab_from_four_vectors(r, th, particle)
            assert point.theta_lab == pytest.approx(angle, abs=1e-9)
            assert point.energy == pytest.approx(energy, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("r", list(reactions()), ids=str)
def test_energy_and_momentum_are_conserved(r):
    t1 = r.beam_energy_mev
    p_beam = math.sqrt(t1 * (t1 + 2 * r.m1))
    for th in (5.0, 60.0, 120.0, 175.0):
        ej, rec = r.at_cm(th), r.recoil_for(th)
        assert ej.energy + rec.energy == pytest.approx(t1 + r.q_value_mev, rel=1e-12, abs=1e-9)
        p3 = math.sqrt(ej.energy * (ej.energy + 2 * r.m3))
        p4 = math.sqrt(rec.energy * (rec.energy + 2 * r.m4))
        a3, a4 = math.radians(ej.theta_lab), math.radians(rec.theta_lab)
        # p from the kinetic energy loses digits for slow heavy particles: T = E − m carries a rounding error of
        # about ε m, so δp ≈ ε m² / p, which platforms round differently (e.g. 4e-7 MeV/c for a 0.8 keV ¹⁹⁷Au).
        eps = np.finfo(float).eps
        tol = 20 * eps * (r.m3**2 / max(p3, 1e-300) + r.m4**2 / max(p4, 1e-300))
        assert p3 * math.cos(a3) + p4 * math.cos(a4) == pytest.approx(p_beam, rel=1e-9, abs=tol)
        assert p3 * math.sin(a3) == pytest.approx(p4 * math.sin(a4), rel=1e-9, abs=tol)


@pytest.mark.parametrize("beam, target", [("4He", "197Au"), ("1H", "12C"), ("12C", "1H"), ("16O", "4He")])
def test_non_relativistic_limit_is_the_kinematic_factor(beam, target):
    r = TwoBody(beam, target, 0.001)  # 1 keV: relativistic corrections ~1e-8
    th = np.array([1.0, 10.0, 30.0, 60.0, 90.0, 120.0, 150.0, 179.0])
    th = th[th < r.max_angle() - 1e-6]
    first, _ = r.at_lab(th)
    k = kinematic_factor(r.m1, r.m2, th)
    np.testing.assert_allclose(first.energy / 0.001, k, rtol=0, atol=1e-6)


def test_non_relativistic_recoil_angle():
    # Elastic, non-relativistic: the recoil leaves at (180° - θ*)/2.
    r = TwoBody("4He", "12C", 0.001)
    for th in (20.0, 70.0, 150.0):
        assert r.recoil_for(th).theta_lab == pytest.approx((180.0 - th) / 2.0, abs=1e-5)


@pytest.mark.parametrize("r", list(reactions()), ids=str)
def test_jacobian_and_broadening_match_numerical_derivatives(r):
    h = 1e-5
    for particle in ("ejectile", "recoil"):
        for th in (15.0, 50.0, 100.0, 160.0):
            p = r.at_cm(th, particle)
            lo, hi = r.at_cm(th - h, particle), r.at_cm(th + h, particle)
            dlab = math.radians(hi.theta_lab - lo.theta_lab)
            dcm = math.radians(2 * h)
            numerical = abs(math.sin(math.radians(th)) * dcm / (math.sin(math.radians(p.theta_lab)) * dlab))
            assert p.jacobian == pytest.approx(numerical, rel=1e-6)
            # A wider step for the energy: E = γ(E* + βp* cos θ*) − m loses digits to the subtraction.
            lo2, hi2 = r.at_cm(th - 1e-3, particle), r.at_cm(th + 1e-3, particle)
            de = (hi2.energy - lo2.energy) / (hi2.theta_lab - lo2.theta_lab)
            assert p.de_dtheta == pytest.approx(de, rel=1e-5, abs=1e-12)
    # At 0 deg the limit (dθ*/dθ_lab)^2 is used; it joins smoothly onto the general formula.
    assert r.at_cm(0.0).jacobian == pytest.approx(r.at_cm(1e-4).jacobian, rel=1e-6)


def test_jacobian_non_relativistic_formula():
    r = TwoBody("4He", "12C", 0.001)
    g = r.m1 / r.m3 * r.m3 / r.m4  # elastic: g = m1/m2
    for th in (10.0, 90.0, 170.0):
        c = math.cos(math.radians(th))
        expected = (1 + g * g + 2 * g * c) ** 1.5 / abs(1 + g * c)
        assert r.at_cm(th).jacobian == pytest.approx(expected, rel=1e-6)


# -- lab angles, both branches, maximum angle -------------------------------------------------------------------


@pytest.mark.parametrize("r", list(reactions()), ids=str)
def test_lab_to_cm_inverts_cm_to_lab(r):
    for particle in ("ejectile", "recoil"):
        cm = np.linspace(0.5, 179.5, 60)
        lab = r.at_cm(cm, particle).theta_lab
        first, second = r.theta_cm(lab, particle)
        found = np.where(np.abs(first - cm) < np.abs(np.nan_to_num(second, nan=1e9) - cm), first, second)
        np.testing.assert_allclose(found, cm, atol=1e-7)
        if not r.double_valued(particle):
            assert np.all(np.isnan(second))


def test_double_valued_cases():
    inv = TwoBody("208Pb", "1H", 1000.0)
    assert inv.double_valued() and not inv.double_valued("recoil")
    # Non-relativistically sin θ_max = m_p/m_Pb; relativity changes it slightly at 4.8 MeV/u.
    assert inv.max_angle() == pytest.approx(math.degrees(math.asin(inv.m2 / inv.m1)), rel=2e-3)
    cm = np.linspace(0, 180, 20001)
    assert inv.at_cm(cm).theta_lab.max() == pytest.approx(inv.max_angle(), abs=1e-6)
    fwd, bwd = inv.at_lab(0.2)
    assert fwd.energy > bwd.energy and fwd.theta_cm < bwd.theta_cm
    assert fwd.theta_lab == pytest.approx(0.2) and bwd.theta_lab == pytest.approx(0.2)
    beyond = inv.at_lab(0.3)
    assert math.isnan(beyond[0].energy) and math.isnan(beyond[1].energy)
    assert TwoBody("16O", "208Pb", 64.0).max_angle() == 180.0
    near = TwoBody("1H", "3H", 1.05, ejectile="n")
    assert near.double_valued()
    # Equal masses, elastic: g = 1 exactly (also relativistically), and the maximum is 90 deg.
    for e in (0.001, 2000.0):
        cc = TwoBody("12C", "12C", e)
        assert cc.max_angle() == 90.0 and not cc.double_valued()
    # Every elastic recoil has g = 1 too.
    assert TwoBody("4He", "197Au", 5.5).max_angle("recoil") == 90.0


def test_q_value_and_threshold():
    tpn = TwoBody("1H", "3H", 2.0, ejectile="n")
    assert tpn.recoil.name == "3He"
    assert tpn.q_value_mev == pytest.approx(-0.7638, abs=1e-3)
    assert tpn.threshold_mev == pytest.approx(1.019, abs=5e-4)  # tabulated 3H(p,n)3He threshold, 1.019 MeV
    assert tpn.threshold_mev == pytest.approx(-tpn.q_value_mev * (1 + tpn.m1 / tpn.m2), rel=1e-3)
    with pytest.raises(ValueError, match="below the threshold of 1.019"):
        TwoBody("1H", "3H", 1.0, ejectile="n")
    dp = TwoBody("2H", "12C", 10.0, ejectile="p")
    assert dp.q_value_mev == pytest.approx(data.q_value(["2H", "12C"], ["1H", "13C"]), abs=1e-4)
    assert dp.threshold_mev == 0.0
    inel = TwoBody("16O", "208Pb", 64.0, excitation_mev=2.614)
    assert inel.q_value_mev == pytest.approx(-2.614, abs=1e-9) and inel.label == "208Pb(16O, 16O)208Pb*"
    ej = TwoBody("16O", "208Pb", 64.0, excitation_mev=6.13, excite="ejectile")
    assert ej.label == "208Pb(16O, 16O*)208Pb" and ej.m3 == pytest.approx(ej.ejectile.nuclear_mass_mev + 6.13)
    assert inel.at_cm(90.0).energy < TwoBody("16O", "208Pb", 64.0).at_cm(90.0).energy


def test_input_errors():
    with pytest.raises(ValueError, match="does not conserve"):
        TwoBody("2H", "12C", 10.0, ejectile="p", recoil="12C")
    with pytest.raises(ValueError, match="no recoil conserves"):
        TwoBody("1H", "1H", 10.0, ejectile="4He")
    with pytest.raises(ValueError, match="excite must be"):
        TwoBody("1H", "12C", 10.0, excite="target")
    with pytest.raises(ValueError, match="total energy"):
        TwoBody("16O", "208Pb", "4 MeV/u")
    with pytest.raises(ValueError, match="particle must be"):
        TwoBody("16O", "208Pb", 64.0).at_cm(10.0, "beam")
    assert TwoBody("16O", "208Pb", "64 MeV").beam_energy_mev == 64.0


def test_elastic_for_an_experiment():
    exp = Experiment.example("oxygen_on_lead_array")
    (key, (r, fraction)), = elastic(exp).items()
    assert key == (82, 208) and fraction == 1.0 and r.beam_energy_mev == pytest.approx(64.0)
    carbon = Experiment.from_dict({**exp.to_dict(), "target": {"material": "C", "thickness": "20 ug/cm2"}})
    parts = elastic(carbon)
    assert sorted(parts) == [(6, 12), (6, 13)]
    assert parts[(6, 13)][1] == pytest.approx(0.0107)


# -- comparison with LISE++ (reference files produced by hand, see tests/reference/nuclear/README.md) -----------


def _reference_rows():
    for path in sorted(REFERENCE.glob("kinematics_*.csv")):
        with open(path, encoding="utf-8") as f:
            for row in csv.DictReader(line for line in f if not line.startswith("#")):
                yield path.name, row


@pytest.mark.skipif(not list(REFERENCE.glob("kinematics_*.csv")), reason="no LISE++ reference files yet")
def test_against_lise_reference_files():
    for name, row in _reference_rows():
        r = TwoBody(row["beam"], row["target"], float(row["beam_energy_mev"]), ejectile=row["ejectile"] or None,
                    excitation_mev=float(row["excitation_mev"] or 0))
        first, second = r.at_lab(float(row["theta_lab_deg"]), row["particle"])
        point = first if row["branch"] in ("", "1") else second
        # 1 keV, or the precision the tool displayed if that is coarser.
        tol = max(1e-3, float(row.get("energy_tol_mev") or 0))
        assert point.energy == pytest.approx(float(row["energy_mev"]), abs=tol), (name, row)
        if row.get("theta_cm_deg"):
            assert point.theta_cm == pytest.approx(float(row["theta_cm_deg"]), abs=0.01), (name, row)


def test_lab_points_for_several_energies_equal_the_one_body_methods():
    """The batched transformation (one row per beam energy, as the rates use it through the target) gives what
    TwoBody.at_lab gives body by body: both branches, the NaNs beyond the maximum angle, single-valued recoils."""
    bodies = [TwoBody("16O", "58Ni", e, excitation_mev=1.454) for e in (20.0, 35.0, 50.0)]
    bodies += [TwoBody("58Ni", "12C", 200.0), TwoBody("4He", "197Au", 16.0)]
    th = np.linspace(0.0, 180.0, 361)
    for particle in ("ejectile", "recoil"):
        first, second = lab_points(bodies, th, particle)
        assert first.energy.shape == (len(bodies), len(th))
        for i, b in enumerate(bodies):
            p1, p2 = b.at_lab(th, particle)
            for got, want in ((first, p1), (second, p2)):
                if want is None:
                    assert np.all(np.isnan(got.theta_cm[i]))
                    continue
                # Where the particle has no energy left (θ* at 180° for a recoil with g = 1) the lab angle is
                # numerically undefined; everything else agrees to rounding.
                live = np.nan_to_num(np.asarray(want.energy)) > 1e-9
                for name in ("theta_cm", "theta_lab", "energy", "jacobian", "de_dtheta"):
                    a, c = getattr(got, name)[i], np.asarray(getattr(want, name))
                    assert np.array_equal(np.isnan(a), np.isnan(c)), (particle, i, name)
                    m = ~np.isnan(c) & (live | (name != "theta_lab"))
                    assert np.allclose(a[m], c[m], rtol=1e-12, atol=1e-12), (particle, i, name)
