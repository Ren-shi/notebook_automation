"""γ-detector response and calibration sources (backlog 53)."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment, response
from physim.nuclear.experiment import SetupError
from physim.nuclear.planner import Planner

CLOVER = {"name": "Clo", "model": "clover", "theta": "90 deg", "phi": "90 deg", "distance": "250 mm"}
LABR = {"name": "La", "model": "LaBr3_2x2", "theta": "90 deg", "phi": "-90 deg", "distance": "150 mm"}
DISC = {"name": "Disc", "theta": "45 deg", "phi": "90 deg", "distance": "120 mm", "radius": "35 mm"}


def setup(*gammas, **extra) -> Experiment:
    d = Experiment.example("coulex_ni58").to_dict()
    d["gamma_detectors"] = [dict(g) for g in gammas]
    d.update(extra)
    return Experiment.from_dict(d)


def resp(detector, **extra) -> response.Response:
    exp = setup(detector, **extra)
    return response.Response(exp, exp.gamma_detectors[0])


# -- attenuation ----------------------------------------------------------------------------------------------------


def test_attenuation_coefficients_are_the_nist_values():
    # NIST XCOM tables (Hubbell and Seltzer), cm²/g.
    assert response.mass_attenuation("Pb", 1.0) == pytest.approx(0.07102, rel=1e-4)
    assert response.mass_attenuation("Al", 1.0) == pytest.approx(0.06146, rel=1e-4)
    assert response.mass_attenuation("Ge", 1.25) == pytest.approx(0.05101, rel=1e-4)
    # Between tabulated points the interpolation is on log–log axes, and stays between its neighbours.
    lo, mid, hi = (float(response.mass_attenuation("Fe", e)) for e in (0.6, 0.662, 0.8))
    assert hi < mid < lo
    # Just above the K edge of lead (88 keV) the coefficient jumps up.
    assert response.mass_attenuation("Pb", 0.0881) > 3 * response.mass_attenuation("Pb", 0.0879)


def test_a_compound_is_the_mass_weighted_sum_of_its_elements():
    la, br = 138.905, 3 * 79.904
    expected = (la * response.mass_attenuation("La", 0.662) + br * response.mass_attenuation("Br", 0.662)) / (la + br)
    assert response.mass_attenuation("LaBr3", 0.662) == pytest.approx(float(expected), rel=1e-3)
    with pytest.raises(ValueError, match="no attenuation data for U"):
        response.mass_attenuation("U", 1.0)
    with pytest.raises(ValueError, match="5 keV to 20 MeV"):
        response.mass_attenuation("Pb", 30.0)


def test_doubling_an_absorber_squares_its_transmission():
    e = np.array([0.122, 0.662, 1.332])
    bare = resp(DISC).peak_efficiency(e)
    one = resp(dict(DISC, absorbers=[["Pb", "1 mm"]])).peak_efficiency(e)
    two = resp(dict(DISC, absorbers=[["Pb", "2 mm"]])).peak_efficiency(e)
    assert np.allclose(two / bare, (one / bare) ** 2, rtol=1e-12)
    mu = response.mass_attenuation("Pb", e) * 11.35  # per cm
    assert np.allclose(one / bare, np.exp(-mu * 0.1), rtol=1e-9)
    # The same thickness given as an areal density.
    areal = resp(dict(DISC, absorbers=[["Pb", "1135 mg/cm2"]])).peak_efficiency(e)
    assert np.allclose(areal, one, rtol=1e-9)


def test_the_chamber_wall_attenuates_detectors_outside_it():
    chamber = {"radius": "160 mm", "wall_thickness": "3 mm", "wall_material": "Al"}
    far = dict(DISC, distance="200 mm")
    inside, outside = resp(far), resp(far, chamber=chamber)
    ratio = float(outside.peak_efficiency(0.122) / inside.peak_efficiency(0.122))
    assert ratio == pytest.approx(float(response.transmission("Al", "3 mm", 0.122)), rel=1e-9)
    assert [a.origin for a in outside.absorbers] == ["chamber wall", "window"]


def test_an_absorber_needs_a_known_material():
    with pytest.raises(SetupError, match="absorber U: no attenuation data"):
        setup(dict(DISC, absorbers=[["U", "1 mm"]]))
    with pytest.raises(SetupError, match="is not a thickness"):
        setup(dict(DISC, absorbers=[["Pb", "1 keV"]]))
    with pytest.raises(SetupError, match="no γ-ray response"):
        setup(dict(DISC, material="NaI"))


# -- efficiency -----------------------------------------------------------------------------------------------------


def test_the_default_clover_has_the_published_efficiency_at_25_cm():
    """Four crystals of 21 to 22% relative efficiency each (Mirion's sheet): 4 × 0.215 × 1.2e-3 at 1332 keV and
    25 cm, less the few percent the end cap absorbs."""
    r = resp(CLOVER)
    eff = float(r.peak_efficiency(response.REFERENCE_MEV))
    assert eff == pytest.approx(4 * 0.215 * response.NAI_STANDARD * float(r.transmission(response.REFERENCE_MEV)),
                                rel=0.03)  # the four crystals sit a little off the axis
    assert 0.95 * 4 * 0.21 * response.NAI_STANDARD < eff < 4 * 0.22 * response.NAI_STANDARD
    one = float(r.peak_efficiency(response.REFERENCE_MEV, element=0))
    assert one == pytest.approx(eff / 4, rel=0.02)


def test_efficiency_falls_with_energy_and_low_energies_are_absorbed():
    r = resp(CLOVER)
    e = np.array([0.03, 0.122, 0.344, 0.779, 1.408, 2.614])
    eff = r.peak_efficiency(e)
    assert np.all(np.diff(eff[1:]) < 0)          # above 100 keV it only falls
    assert eff[0] < eff[1]                        # the end cap takes the X-rays
    # Between 344 and 1408 keV a coaxial germanium detector's efficiency falls roughly as E^-0.8.
    slope = math.log(eff[4] / eff[2]) / math.log(1.408 / 0.344)
    assert -1.0 < slope < -0.6
    assert np.all(r.total_efficiency(e) >= eff) and np.all(r.total_efficiency(e) <= sum(r.coverage))
    assert np.allclose(r.peak_to_total(e)[[0, 1]], [0.95, 0.18 * (0.122 / 1.332492) ** -0.68])


def test_a_bigger_or_longer_crystal_is_more_efficient():
    small, long_ = resp(dict(DISC, radius="25 mm")), resp(dict(DISC, radius="25 mm", crystal_diameter="50 mm",
                                                               crystal_length="90 mm", crystals=1))
    assert small.assumed_length and not long_.assumed_length
    assert float(long_.peak_efficiency(1.332)) > float(small.peak_efficiency(1.332))
    assert float(resp(DISC).peak_efficiency(1.332)) > float(small.peak_efficiency(1.332))
    assert resp(LABR).material == "LaBr3" and resp(DISC).material == "Ge"


def test_a_measured_curve_or_a_single_efficiency_replaces_the_model():
    curve = [["122 keV", "1.0 %"], ["344 keV", "0.5 %"], ["1408 keV", "0.15 %"]]
    r = resp(dict(DISC, efficiency_curve=curve))
    assert r.source == "curve"
    assert np.allclose(r.peak_efficiency([0.122, 0.344, 1.408]), [0.010, 0.005, 0.0015])
    between = float(r.peak_efficiency(0.2))
    assert between == pytest.approx(0.010 * (0.2 / 0.122) ** (math.log(0.5) / math.log(0.344 / 0.122)), rel=1e-9)
    # Beyond the last point it follows the model's shape from there.
    model = resp(DISC)
    assert float(r.peak_efficiency(2.0)) == pytest.approx(
        0.0015 * float(model.peak_efficiency(2.0) / model.peak_efficiency(1.408)), rel=1e-9)
    fixed = resp(dict(DISC, efficiency="2 %"))
    assert fixed.source == "fixed" and np.allclose(fixed.peak_efficiency([0.1, 1.0, 2.0]), 0.02)
    assert model.source == "model"
    with pytest.raises(SetupError, match="at least two points"):
        setup(dict(DISC, efficiency_curve=[["122 keV", "1 %"]]))
    # The curve and the absorbers survive a round trip through a setup file.
    exp = setup(dict(DISC, efficiency_curve=curve, absorbers=[["Cu", "0.5 mm"]]))
    again = Experiment.from_toml(exp.to_toml())
    assert again.gamma_detectors[0] == exp.gamma_detectors[0]


def test_resolution_depends_on_energy():
    ge = resp(CLOVER)
    # Mirion's sheet: 2.1 keV at 1332 keV and 1.05 keV at 122 keV.
    assert float(ge.fwhm(1.332492)) == pytest.approx(2.1e-3)
    assert float(ge.fwhm(0.122)) == pytest.approx(1.05e-3, rel=0.03)
    la = resp(LABR)
    # Saint-Gobain's note: 2.9% at 662 keV, 2.1% at 1332 keV, 1.6% at 2615 keV.
    for e, percent in ((0.662, 2.9), (1.332, 2.1), (2.615, 1.6)):
        assert 100 * float(la.fwhm(e)) / e == pytest.approx(percent, rel=0.1)


# -- the shape of the spectrum --------------------------------------------------------------------------------------


def test_the_response_to_one_gamma_ray_is_a_peak_a_continuum_and_escape_peaks():
    r = resp(CLOVER)
    edges = np.arange(0.0, 3.0, 0.0005)
    centres = (edges[:-1] + edges[1:]) / 2
    low = r.shape(0.662, edges)
    assert low.sum() == pytest.approx(1.0, abs=1e-3)
    peak = np.abs(centres - 0.662) < 0.005
    assert low[peak].sum() == pytest.approx(float(r.peak_to_total(0.662)), rel=0.02)  # with the continuum under it
    edge = response.compton_edge(0.662)
    assert edge == pytest.approx(0.4777, abs=1e-4)
    # Most of the rest lies below the Compton edge, a little between the edge and the peak.
    below = low[centres < edge + 0.003].sum()
    valley = low[(centres > edge + 0.005) & (centres < 0.655)].sum()
    assert below > 5 * valley > 0
    assert r.escape_fractions(0.662) == (0.0, 0.0)
    high = r.shape(2.614, edges)
    se, de = r.escape_fractions(2.614)
    assert 0 < se < de < 0.12
    for share, where in ((se, 2.614 - 0.511), (de, 2.614 - 1.022)):
        line = np.abs(centres - where) < 0.006
        flat = high[np.abs(centres - where - 0.02) < 0.006].sum()
        assert high[line].sum() - flat == pytest.approx(share, rel=0.05)


def test_the_compton_continuum_follows_klein_nishina():
    edges = np.linspace(0.0, 0.662, 400)
    h = response.compton_continuum(0.662, edges)
    assert h.sum() == pytest.approx(1.0)
    edge = response.compton_edge(0.662)
    centres = (edges[:-1] + edges[1:]) / 2
    assert h[centres > edge + 0.004].sum() == 0
    # The spectrum rises to the edge: its last bins below the edge hold more than the middle ones.
    inside = h[centres < edge - 0.004]
    assert inside[-1] > 1.5 * inside[len(inside) // 2]
    k = 0.662 / response.ELECTRON_MEV
    kn = lambda t: 2 + (t / 0.662) ** 2 / (k**2 * (1 - t / 0.662) ** 2) + (t / 0.662) / (1 - t / 0.662) * (  # noqa: E731
        t / 0.662 - 2 / k)
    assert inside[-1] / inside[10] == pytest.approx(kn(centres[len(inside) - 1]) / kn(centres[10]), rel=0.02)
    with_flat = response.compton_continuum(0.662, edges, multiple=0.2)
    assert with_flat[centres > edge + 0.004].sum() == pytest.approx(0.2, rel=0.03)


# -- calibration sources --------------------------------------------------------------------------------------------


def test_the_sources_carry_the_evaluated_lines():
    assert response.source_names() == ["22Na", "60Co", "88Y", "133Ba", "137Cs", "152Eu"]
    co = response.source("60Co")
    assert [round(1e3 * ln.energy_mev, 3) for ln in co.lines] == [1173.228, 1332.492]
    assert co.lines[1].intensity == pytest.approx(0.999826)
    assert co.half_life_s == pytest.approx(5.2711 * 365.2422 * 86400, rel=2e-3)
    eu = response.source("152Eu")
    strong = {round(1e3 * ln.energy_mev, 1): ln.intensity for ln in eu.strong()}
    assert strong[121.8] == pytest.approx(0.2841) and strong[1408.0] == pytest.approx(0.2085)
    assert strong[344.3] == pytest.approx(0.2659) and len(strong) == 11
    assert any(ln.kind == "x" for ln in eu.lines)
    with pytest.raises(ValueError, match="no decay data"):
        response.source("57Co")


def test_a_europium_run_gives_back_the_efficiency_that_was_put_in():
    exp = setup(CLOVER, DISC)
    run = response.source_run(exp, "152Eu", activity="37 kBq", time="1 h", seed=11)
    assert run.decays == pytest.approx(37e3 * 3600, rel=1e-4)
    for name in ("Clo", "Disc", "Clo B"):
        points = run.efficiency_points(name)
        assert len(points) >= 8, name
        for p in points:
            # Within statistics; the 964 keV line carries a weak neighbour 0.7 keV below it (1%).
            allowed = 4 * p["uncertainty"] + (0.012 if abs(p["energy_mev"] - 0.964) < 0.001 else 0.003) * p["true"]
            assert abs(p["efficiency"] - p["true"]) < allowed, (name, p)
    # The four crystals add up to the detector.
    crystals = sum(run.detector(f"Clo {c}") for c in "ABCD")
    assert np.array_equal(crystals, run.detector("Clo"))
    with pytest.raises(KeyError):
        run.detector("nothing")


def test_without_counting_noise_the_peak_areas_are_exact():
    run = response.source_run(setup(CLOVER), "152Eu", activity="37 kBq", time="1 h")
    run.spectra = run.expected
    for p in run.efficiency_points("Clo"):
        if abs(p["energy_mev"] - 0.964) > 0.001:
            assert p["efficiency"] == pytest.approx(p["true"], rel=3e-3), p


def test_the_cobalt_peak_to_total_is_the_value_put_in():
    exp = setup(CLOVER, LABR)
    run = response.source_run(exp, "60Co", activity="37 kBq", time="2 h", seed=5)
    ge = run.peak_to_total("Clo")
    assert abs(ge["value"] - ge["true"]) < 4 * ge["uncertainty"] + 0.002 * ge["true"]
    # For ⁶⁰Co the ratio is that of the response at its two lines.
    r = run.responses["Clo A"][0]
    assert ge["true"] == pytest.approx(float(np.mean(r.peak_to_total([1.173228, 1.332492]))), rel=0.01)
    # In LaBr₃ the 1173 keV peak stands on the Compton edge of the 1332 keV line, which the simple area misjudges.
    la = run.peak_to_total("La")
    assert la["value"] == pytest.approx(la["true"], rel=0.06)


def test_a_run_scales_with_activity_and_time_and_the_source_decays():
    exp = setup(CLOVER)
    short = response.source_run(exp, "60Co", "10 kBq", "10 min")
    long_ = response.source_run(exp, "60Co", "20 kBq", "30 min")
    assert long_.decays == pytest.approx(6 * short.decays, rel=1e-5)
    assert long_.expected["Clo A"].sum() == pytest.approx(6 * short.expected["Clo A"].sum(), rel=1e-5)
    # ⁸⁸Y (107 days) has lost half its activity after one half-life: the decays are A₀ T½ / ln 2 × ½.
    y = response.source_run(exp, "88Y", "1 kBq", f"{106.63 * 24} h")
    assert y.decays == pytest.approx(1e3 * response.source("88Y").half_life_s / math.log(2) / 2, rel=1e-3)
    assert not np.array_equal(response.source_run(exp, "60Co", seed=1).spectra["Clo A"],
                              response.source_run(exp, "60Co", seed=2).spectra["Clo A"])
    with pytest.raises(ValueError, match="no γ-ray detectors"):
        response.source_run(Experiment.example("alpha_on_gold"))


# -- in the planner -------------------------------------------------------------------------------------------------


def test_the_planner_counts_coincidences_with_the_response():
    p = Planner.example("coulex_ni58")
    r = p.rates()
    e0 = p.experiment.excitation.energy_mev
    eff = sum(float(response.Response(p.experiment, g).peak_efficiency(e0)) for g in p.experiment.gamma_detectors)
    assert r["gamma_efficiency"] == pytest.approx(eff) and r["gamma_efficiency_typical"]
    # Far below the geometric coverage that was used before.
    assert eff < 0.2 * sum(g.geometric_efficiency() for g in p.experiment.gamma_detectors)
    assert p.set("gamma detector 1", "absorbers", [["Pb", "5 mm"]])
    assert p.rates()["gamma_efficiency"] < eff
    curve = p.efficiency("gamma:0")
    assert curve["source"] == "model" and [a["origin"] for a in curve["absorbers"]] == ["absorber", "window"]
    assert len(curve["energy_mev"]) == len(curve["peak"]) == len(curve["total"])
    sel = p.selection("gamma:0")
    assert sel["gamma_energy_kev"] == pytest.approx(1e3 * e0) and sel["response"]["material"] == "Ge"
    assert sel["peak_efficiency"] == pytest.approx(float(np.interp(e0, curve["energy_mev"], curve["peak"])), rel=0.01)
    run = p.source_run("60Co", "10 kBq", "10 min")
    assert p.source_run("60Co", "10 kBq", "10 min") is run       # kept until the setup changes


def test_the_app_figures_and_fields():
    pytest.importorskip("plotly")
    from physim.nuclear import app

    p = Planner.example("coulex_ni58")
    run = p.source_run("152Eu", "37 kBq", "20 min")
    name = p.experiment.gamma_detectors[0].name
    fig = app.figure_efficiency(p, "gamma:0", run.efficiency_points(name))
    assert [t.name for t in fig.data] == ["full-energy peak", "any energy", "source run"]
    app.figure_source_spectrum(run, name).to_json()
    assert app._value("absorbers", "Pb 1 mm, Cu 0.5 mm") == [["Pb", "1 mm"], ["Cu", "0.5 mm"]]
    assert app._shown("absorbers", [["Pb", "1 mm"], ["Cu", "0.5 mm"]]) == "Pb 1 mm, Cu 0.5 mm"
    assert not p.set("gamma detector 1", "absorbers", app._value("absorbers", "lead"))
