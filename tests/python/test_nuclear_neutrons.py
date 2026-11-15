"""Neutron detectors and the channel-selection matrix (physim.nuclear.neutrons, backlog item 75)."""

import csv
import math
from dataclasses import replace
from pathlib import Path
from statistics import NormalDist
from types import SimpleNamespace

import numpy as np
import pytest

from physim.nuclear import neutrons as N


@pytest.fixture(scope="module")
def cell():
    return N.NeutronDetector()


# -- cross sections and light output ------------------------------------------------------------------------------


@pytest.mark.parametrize("energy, sigma_b", [(1.0, 4.26), (14.1, 0.689)])
def test_np_cross_section_matches_the_standard(energy, sigma_b):
    # The n-p total cross section, the neutron standard: about 4.26 b at 1 MeV and 0.689 b at 14.1 MeV
    # (e.g. Knoll, Radiation Detection and Measurement, 4th ed., fig. 15.1; ENDF/B-VIII.0 H-1, MT=1).
    assert float(N.hydrogen_cross_section(energy)) == pytest.approx(sigma_b, rel=0.01)


def test_carbon_cross_section_is_interpolated_from_the_table():
    t = N._table()["carbon_total"]
    assert t["typical"] is False and "ENDF/B-VIII.0" in t["source"]
    assert float(N.carbon_cross_section(t["energy_mev"][5])) == pytest.approx(t["sigma_b"][5])
    assert float(N.carbon_cross_section(14.0)) == pytest.approx(1.33, rel=0.02)


# -- the cross sections against ENDF/B-VIII.0 (scripts/endf_neutron_reference.py) ---------------------------------

ENDF = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "neutron_carbon_endf.csv"


def _endf() -> dict:
    """{quantity: [(energy, point, mean ±2.5 %, wide ±10 %)]} from the reference file."""
    if not ENDF.exists():
        pytest.skip("tests/reference/nuclear/neutron_carbon_endf.csv not made yet (scripts/endf_neutron_reference.py)")
    text = ENDF.read_text(encoding="utf-8").splitlines()
    head = [ln for ln in text if ln.startswith("#")]
    assert any(ln.startswith("# tool: ENDF/B-VIII.0") for ln in head) and any("produced by" in ln for ln in head)
    out: dict = {}
    for row in csv.DictReader(ln for ln in text if not ln.startswith("#")):
        out.setdefault(row["quantity"], []).append(
            tuple(float(row[k]) for k in ("energy_mev", "point_b", "mean_b", "wide_b")))
    return out


def test_carbon_table_is_the_endf_evaluation():
    """physim's natural-carbon total cross section (the attenuation in a cell and in cross-talk) against
    ENDF/B-VIII.0 at 0.5–20 MeV. The table holds averages over bins of ±2.9 %, so: within 0.2 % where the cross
    section is smooth (below 1.8 MeV), within 7 % of the ±2.5 % average everywhere (the resonances between 2 and 11
    MeV are narrower than a bin), and within 5 % of the ±10 % average (the narrow 2.08 MeV peak is spread over the
    neighbouring bins)."""
    rows = _endf()["natC_total"]
    assert min(r[0] for r in rows) <= 0.5 and max(r[0] for r in rows) >= 20.0
    for e, point, mean, wide in rows:
        ours = float(N.carbon_cross_section(e))
        x = np.linspace(0.9 * e, 1.1 * e, 4001)
        ours_wide = float(np.trapezoid(N.carbon_cross_section(x), x)) / (0.2 * e)
        assert ours == pytest.approx(mean, rel=0.002 if e < 1.8 else 0.07), e
        assert ours_wide == pytest.approx(wide, rel=0.002 if e < 1.7 else 0.05), e


def test_np_formula_agrees_with_endf():
    """Gammel's closed form against the n–p total cross section of ENDF/B-VIII.0 (H-1, MT=1), 0.5–20 MeV."""
    for e, point, mean, wide in _endf()["H1_total"]:
        assert float(N.hydrogen_cross_section(e)) == pytest.approx(point, rel=0.01), e


def test_light_output_and_its_inverse(cell):
    lo = cell.light_output
    # Cecil et al. (1979) for NE213: a 1 MeV proton gives about 0.21 MeVee
    assert float(lo(1.0)) == pytest.approx(0.83 - 2.82 * (1 - math.exp(-0.25)), rel=1e-12)
    e_th = cell.proton_threshold()
    assert float(lo(e_th)) == pytest.approx(cell.threshold_mevee, abs=1e-9)
    assert 0.5 < e_th < 0.8  # 100 keVee is about 0.6 to 0.7 MeV of proton energy in NE213


def test_materials_and_aliases():
    assert N.material("NE213") == N.material("EJ-301")
    assert N.material("EJ-301")["hydrogen"] == pytest.approx(4.82e22)
    with pytest.raises(ValueError, match="unknown scintillator"):
        N.material("water")


# -- efficiency: the textbook closed form ---------------------------------------------------------------------------


@pytest.mark.parametrize("energy", [1.0, 2.5, 8.0])
def test_efficiency_equals_the_closed_form_of_knoll(cell, energy):
    """Knoll, Radiation Detection and Measurement, 4th ed. (2010), chapter 15 (recoil-proton scintillators):
    ε = n_H σ_H / (n_H σ_H + n_C σ_C) · (1 − exp(−(n_H σ_H + n_C σ_C) d)) · (1 − E_B / E_n)."""
    s = N.material("EJ-301")
    d_cm = 15.0
    sh = s["hydrogen"] * float(N.hydrogen_cross_section(energy)) * 1e-24
    sc = s["carbon"] * float(N.carbon_cross_section(energy)) * 1e-24
    e_b = cell.proton_threshold()
    expected = sh / (sh + sc) * (1 - math.exp(-(sh + sc) * d_cm)) * (1 - e_b / energy)
    assert cell.intrinsic_efficiency(energy) == pytest.approx(expected, rel=1e-12)


def test_thin_hydrogen_limit(cell):
    # A thin cell without carbon: ε → n_H σ_H d (1 − E_th / E), the probability of one n-p scatter above threshold.
    thin = N.NeutronDetector(depth="0.1 mm")
    e = 3.0
    n_sigma_d = 4.82e22 * float(N.hydrogen_cross_section(e)) * 1e-24 * 0.01
    expected = n_sigma_d * (1 - thin.proton_threshold() / e)
    assert thin.intrinsic_efficiency(e, carbon=False) == pytest.approx(expected, rel=1e-3)


def test_threshold_zeroes_the_efficiency_below_it(cell):
    e_th = cell.proton_threshold()
    assert cell.intrinsic_efficiency(0.9 * e_th) == 0.0
    assert cell.intrinsic_efficiency(1.1 * e_th) > 0.0
    higher = N.NeutronDetector(threshold="0.3 MeV")
    assert higher.intrinsic_efficiency(2.0) < cell.intrinsic_efficiency(2.0)


def test_efficiency_includes_the_solid_angle(cell):
    d, r = 510.0, 63.5
    omega = 2 * math.pi * (1 - d / math.hypot(d, r))
    assert cell.solid_angle() == pytest.approx(omega * 1e3, rel=1e-12)
    assert cell.efficiency(2.0) == pytest.approx(cell.intrinsic_efficiency(2.0) * omega / (4 * math.pi))
    arr = cell.intrinsic_efficiency(np.array([1.0, 2.0, 5.0]))
    assert arr.shape == (3,) and arr[1] == pytest.approx(cell.intrinsic_efficiency(2.0))


# -- time of flight -----------------------------------------------------------------------------------------------


def test_time_of_flight(cell):
    path = 1000.0
    v = math.sqrt(2 * 0.001 / N.NEUTRON_MEV) * N.C_MM_PER_NS  # 1 keV: non-relativistic
    assert N.time_of_flight(0.001, path) == pytest.approx(path / v, rel=1e-6)
    for e in (0.5, 3.0, 20.0):
        assert N.energy_from_tof(N.time_of_flight(e, path), path) == pytest.approx(e, rel=1e-10)
    assert cell.gamma_flash() == pytest.approx(cell.flight_path_mm() / N.C_MM_PER_NS)
    assert cell.time_of_flight(2.0) > cell.gamma_flash()
    with pytest.raises(ValueError):
        N.energy_from_tof(1.0, 1000.0)


def test_energy_resolution_non_relativistic_limit():
    # A thin cell: σ_E / E → 2 σ_t / t.
    thin = N.NeutronDetector(depth="0.001 mm", time_resolution="1 ns")
    e = 0.01
    t = thin.time_of_flight(e)
    expected = 2 * (1.0 / N.FWHM_PER_SIGMA) / t * e * N.FWHM_PER_SIGMA
    assert thin.energy_resolution(e) == pytest.approx(expected, rel=1e-4)


def test_misidentification_from_figure_of_merit():
    for fom in (0.5, 1.0, 1.5):
        tail = 1 - NormalDist().cdf(fom * N.FWHM_PER_SIGMA)
        assert N.misidentification_from_fom(fom) == pytest.approx(tail, rel=1e-9)
    assert N.misidentification_from_fom(0.0) == pytest.approx(0.5)


# -- cross-talk ---------------------------------------------------------------------------------------------------


def test_crosstalk_needs_twice_the_threshold(cell):
    e_th = cell.proton_threshold()
    assert cell.crosstalk([130, 0, 0], 1.9 * e_th) == 0.0
    assert cell.crosstalk([130, 0, 0], 3.0) > 0.0


def test_crosstalk_falls_with_distance_and_is_symmetric(cell):
    near = cell.crosstalk([130, 0, 0], 5.0, samples=100_000)
    others = [cell.crosstalk(off, 5.0, samples=100_000, seed=s)
              for s, off in enumerate(([-130, 0, 0], [0, 130, 0], [0, -130, 0]), start=2)]
    for o in others:
        assert o == pytest.approx(near, rel=0.05)
    far = cell.crosstalk([400, 0, 0], 5.0, samples=100_000)
    assert 0 < far < near / 3
    assert cell.crosstalk([3000, 0, 0], 5.0, samples=20_000) < 1e-3


def test_wall_crosstalk_sums_its_neighbours(cell):
    wall = N.NeutronWall(cell, cells=50, pitch="130 mm")
    one = cell.crosstalk([130, 0, 0], 5.0)
    assert wall.crosstalk(5.0) == pytest.approx(6 * one, rel=0.05)
    edge = N.NeutronWall(cell, cells=50, pitch="130 mm", edge_share=0.25)
    assert edge.crosstalk(5.0) == pytest.approx(0.75 * wall.crosstalk(5.0))


def test_wall_geometry(cell):
    wall = N.NeutronWall(cell, cells=50)
    assert wall.geometric_efficiency() == pytest.approx(50 * cell.solid_angle() / (4e3 * math.pi))
    c1 = math.cos(math.radians(wall.theta_max()))
    assert 2 * math.pi * (1 - c1) * 1e3 == pytest.approx(wall.solid_angle(), rel=1e-12)
    # an isotropic neutron is seen with the wall's efficiency
    rng = np.random.default_rng(5)
    c = 2 * rng.random(200_000) - 1
    dirs = np.stack([np.sqrt(1 - c * c), np.zeros_like(c), c], axis=1)
    mc = wall.detection_probability(np.full(len(c), 3.0), dirs).mean()
    assert mc == pytest.approx(wall.efficiency(3.0), rel=0.02)


# -- channels -----------------------------------------------------------------------------------------------------


def test_channel_names():
    assert N.channel_name(2, 1, 0) == "p2n"
    assert N.channel_name(1, 0, 1) == "αn"
    assert N.channel_name(1, 2, 0) == "2pn"
    for name in ("p2n", "αn", "2pn", "3n", "α2p3n"):
        assert N.channel_name(*N.parse_channel_name(name)) == name
    assert N.parse_channel_name("a3n") == (3, 0, 1)
    with pytest.raises(ValueError):
        N.parse_channel_name("2x")
    with pytest.raises(ValueError, match="same particles"):
        N.channels_from({"a": (1, 0, 0), "b": (1, 0, 0)})


CHANNELS = N.channels_from(["n", "2n", "3n", "pn", "p2n", "2pn", "αn"], gammas=8)

FULL = N.Detection(neutron_efficiency=0.25, proton_efficiency=0.6, alpha_efficiency=0.4, crosstalk=0.06,
                   gamma_efficiency=0.1, gamma_as_neutron=0.02, neutron_as_gamma=0.03,
                   proton_as_alpha=0.02, alpha_as_proton=0.05)


def test_analytic_rows_sum_to_one():
    m = N.selection_matrix(CHANNELS, FULL)
    assert m.assigned[-1] == N.UNASSIGNED
    assert np.allclose(m.matrix.sum(axis=1), 1.0, atol=1e-12)
    assert (m.matrix >= 0).all()


def test_analytic_binomial_closed_forms():
    eps = 0.3
    det = N.Detection(neutron_efficiency=eps)
    ch = N.channels_from(["n", "2n", "3n"])
    m = N.selection_matrix(ch, det)
    assert m.share("3n", "3n") == pytest.approx(eps ** 3)
    assert m.share("3n", "2n") == pytest.approx(3 * eps ** 2 * (1 - eps))
    assert m.share("3n", "n") == pytest.approx(3 * eps * (1 - eps) ** 2)
    assert m.share("3n", N.UNASSIGNED) == pytest.approx((1 - eps) ** 3)
    # p2n with protons and alphas mixed up
    det = N.Detection(neutron_efficiency=eps, proton_efficiency=0.5, proton_as_alpha=0.1)
    ch = N.channels_from(["2n", "p2n", "α2n"])
    m = N.selection_matrix(ch, det)
    assert m.share("p2n", "α2n") == pytest.approx(eps ** 2 * 0.5 * 0.1)
    assert m.share("p2n", "2n") == pytest.approx(eps ** 2 * 0.5)


def test_gamma_rays_counted_as_neutrons():
    g = 0.1 * 0.02
    det = N.Detection(gamma_efficiency=0.1, gamma_as_neutron=0.02)
    ch = N.channels_from(["n", "2n"], gammas=8)
    m = N.selection_matrix(ch, det)  # no neutron is ever seen, only γ rays posing as neutrons
    assert m.share("n", "n") == pytest.approx(8 * g * (1 - g) ** 7)
    assert m.share("n", "2n") == pytest.approx(28 * g ** 2 * (1 - g) ** 6)


def test_crosstalk_raises_the_apparent_2n_share_of_a_true_1n_channel():
    ch = N.channels_from(["n", "2n"])
    eps, c = 0.3, 0.08
    clean = N.selection_matrix(ch, N.Detection(neutron_efficiency=eps))
    talk = N.selection_matrix(ch, N.Detection(neutron_efficiency=eps, crosstalk=c))
    assert clean.share("n", "2n") == 0.0
    assert talk.share("n", "2n") == pytest.approx(eps * c)
    assert talk.share("n", "n") == pytest.approx(eps * (1 - c))
    # and the Monte Carlo sees the same
    events = N.isotropic_events(ch, weights=[1, 0], events=40_000, seed=11)
    mc = N.simulate_selection(events, ch, N.Detection(neutron_efficiency=eps, crosstalk=c), seed=12)
    assert mc.share("n", "2n") == pytest.approx(eps * c, abs=4 * mc.uncertainty()[0, 1])


def _within_statistics(mc, analytic, sigmas=4.5):
    # binomial errors from the expected shares, so cells with no expected events are judged too
    n = mc.events[:, None]
    err = np.sqrt(analytic.matrix * (1 - analytic.matrix) / n) + 1.0 / n
    z = np.abs(mc.matrix - analytic.matrix) / err
    return z.max()


def test_monte_carlo_agrees_with_the_analytic_matrix():
    """Item 75 'Detectors' done-when: the event generator's matrix agrees with the binomial expectation."""
    events = N.isotropic_events(CHANNELS, events=60_000, seed=21)
    mc = N.simulate_selection(events, CHANNELS, FULL, seed=22)
    m = N.selection_matrix(CHANNELS, FULL)
    assert mc.events.sum() == 60_000
    assert _within_statistics(mc, m) < 4.5


def test_monte_carlo_with_per_particle_efficiencies_agrees_after_averaging(cell):
    """Efficiencies that depend on each neutron's energy and direction (the wall's cone and the cell's curve)."""
    wall = N.NeutronWall(cell, cells=50)
    det = N.Detection(neutron_efficiency=wall.detection_probability, proton_efficiency=0.7, alpha_efficiency=0.5,
                      crosstalk=lambda e: 0.01 * np.clip(np.asarray(e) - 1.0, 0, None), proton_as_alpha=0.01,
                      gamma_efficiency=0.05, gamma_as_neutron=0.02)
    events = N.isotropic_events(CHANNELS, events=60_000, seed=31)
    with pytest.raises(ValueError, match="averaged"):
        N.selection_matrix(CHANNELS, det)
    m = N.selection_matrix(CHANNELS, det.averaged(events, CHANNELS))
    mc = N.simulate_selection(events, CHANNELS, det, seed=32)
    assert np.allclose(m.matrix.sum(axis=1), 1.0)
    assert _within_statistics(mc, m) < 4.5
    # the averaged neutron efficiency is the wall's efficiency over the events' spectrum
    av = det.averaged(events, CHANNELS)["2n"]
    assert 0.0 < av.neutron_efficiency < wall.geometric_efficiency()


def test_a_varying_number_of_gamma_rays_mixes_the_binomials():
    f = 0.1 * 0.02
    det = N.Detection(gamma_efficiency=0.1, gamma_as_neutron=0.02)
    # half the events have no γ ray, half have two: P(1 false neutron) = ½ · 2 f (1 − f)
    ch = N.channels_from({"n": (1, 0, 0), "2n": (2, 0, 0)}, gammas=(0.5, 0.0, 0.5))
    assert ch[0].gammas == (0.5, 0.0, 0.5)
    m = N.selection_matrix(ch, det)
    assert m.share("n", "n") == pytest.approx(0.5 * 2 * f * (1 - f))
    assert m.share("n", "2n") == pytest.approx(0.5 * f * f)
    assert np.allclose(m.matrix.sum(axis=1), 1.0)
    # a fixed number is the one-point distribution
    fixed = N.selection_matrix(N.channels_from(["n", "2n"], gammas=8), det)
    spread = N.selection_matrix(N.channels_from(["n", "2n"], gammas=(0,) * 8 + (1,)), det)
    assert np.allclose(fixed.matrix, spread.matrix)


def test_monte_carlo_agrees_with_the_analytic_matrix_with_thresholds_and_varying_gamma_rays():
    """Particle-detector thresholds (efficiencies that are zero below an energy) and a γ-ray number that varies from
    event to event: the analytic matrix from Detection.averaged and gamma_distributions over the same events is
    what the Monte Carlo converges to."""
    rng = np.random.default_rng(51)
    events = N.isotropic_events(CHANNELS, events=60_000, seed=51)
    for ev in events:
        ev["gammas"] = int(rng.poisson(10))
    det = replace(FULL, proton_efficiency=lambda e, d: 0.6 * (np.asarray(e) > 6.0),
                  alpha_efficiency=lambda e, d: 0.5 * (np.asarray(e) > 13.0), gamma_efficiency=0.2)
    m = N.selection_matrix(N.gamma_distributions(events, CHANNELS), det.averaged(events, CHANNELS))
    mc = N.simulate_selection(events, CHANNELS, det, seed=52)
    assert np.allclose(m.matrix.sum(axis=1), 1.0)
    assert _within_statistics(mc, m) < 4.5
    # the threshold removes the slow protons: the averaged efficiency is 0.6 times the share above 6 MeV
    p = [e for ev in events for k, e, _ in ev["particles"] if k == "p"]
    assert det.averaged(events, CHANNELS)["p2n"].proton_efficiency == pytest.approx(
        0.6 * np.mean(np.asarray(p) > 6.0), abs=0.01)
    # and using a fixed γ number instead of the events' own spread would be visibly wrong for γ-only fakes
    g = N.gamma_distributions(events, CHANNELS)
    assert g[0].gamma_counts.sum() == pytest.approx(1.0) and len(g[0].gamma_counts) > 15


def test_seeded_runs_are_identical():
    events = N.isotropic_events(CHANNELS, events=5_000, seed=41)
    a = N.simulate_selection(events, CHANNELS, FULL, seed=7)
    b = N.simulate_selection(events, CHANNELS, FULL, seed=7)
    assert np.array_equal(a.matrix, b.matrix)


def test_event_formats():
    ch = N.channels_from(["pn", "2n"])
    det = N.Detection(neutron_efficiency=1.0, proton_efficiency=1.0)
    z = (0.0, 0.0, 1.0)
    events = [
        [("n", 2.0, z), ("proton", 5.0, z)],                                                    # tuples
        {"particles": [{"kind": "n", "energy_mev": 2.0, "direction": z},
                       {"type": "neutron", "energy": 1.0, "direction": z}], "gammas": 3},      # dicts
        SimpleNamespace(particles=[SimpleNamespace(kind="n", energy_mev=1.0, direction=z),
                                   SimpleNamespace(kind="p", energy_mev=4.0, direction=z)]),   # objects
        [("n", 2.0, z), ("n", 2.0, z), ("n", 2.0, z)],                                         # 3n: no such channel
    ]
    mc = N.simulate_selection(events, ch, det)
    assert list(mc.events) == [2, 1]
    assert mc.share("pn", "pn") == 1.0 and mc.share("2n", "2n") == 1.0
    with pytest.raises(ValueError, match="unknown particle kind"):
        N.simulate_selection([[("deuteron", 1.0, z)]], ch, det)


def test_purity():
    ch = N.channels_from(["n", "2n"])
    m = N.selection_matrix(ch, N.Detection(neutron_efficiency=0.5))
    # 2n events land in "n" with 2·0.5·0.5 = 0.5, n events with 0.5: equal cross sections give purity 1/2
    assert m.purity({"n": 1.0, "2n": 1.0})["n"] == pytest.approx(0.5)
    assert m.purity({"n": 1.0, "2n": 1.0})["2n"] == pytest.approx(1.0)


# -- the detector as setup data -----------------------------------------------------------------------------------


def test_detector_fields_are_text_with_units():
    d = N.NeutronDetector(radius="6.35 cm", distance="1 m", threshold="50 keV", neutron_as_gamma=0.05)
    assert d.radius_mm == pytest.approx(63.5)
    assert d.distance_mm == pytest.approx(1000.0)
    assert d.threshold_mevee == pytest.approx(0.05)
    assert d.misidentification == (0.05, 0.01)
    assert d.check() == []
    assert set(N.TYPICAL_FIELDS) <= set(N.NeutronDetector.__dataclass_fields__)
    bad = N.NeutronDetector(radius="-1 mm", depth="1 MeV", neutron_as_gamma="150 %", material="water")
    problems = bad.check()
    assert len(problems) == 4
    det = N.Detection.from_detector(d, 0.2, gamma_efficiency=0.1)
    assert det.neutron_as_gamma == 0.05 and det.gamma_efficiency == 0.1
