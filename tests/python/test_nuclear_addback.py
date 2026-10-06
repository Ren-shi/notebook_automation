"""Add-back and Compton suppression for clovers (backlog 62)."""

import numpy as np
import pytest

from physim.nuclear import Experiment, guide, response
from physim.nuclear.detectors import Array
from physim.nuclear.events import simulate
from physim.nuclear.experiment import SetupError
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates
from physim.nuclear.scene import solids

EVENTS = 200_000
PER_EVENT = 60
E_REF = response.REFERENCE_MEV


def setup(**clover) -> Experiment:
    d = Experiment.example("coulex_ni58").to_dict()
    d["gamma_detectors"] = [{"name": "Clo", "model": "clover", "theta": "90 deg", "phi": "90 deg",
                             "distance": "100 mm", **clover}]
    return Experiment.from_dict(d)


def peak(g) -> tuple:
    m = g["counted"] & (g["deposited"] == g["energy_lab"])
    return float(g["weight"][m].sum()), int(m.sum())


def continuum(g) -> tuple:
    m = g["counted"] & (g["deposited"] < g["energy_lab"])
    return float(g["weight"][m].sum()), int(m.sum())


@pytest.fixture(scope="module")
def runs() -> dict:
    """The same particle events and γ-ray chain: bare, with add-back, and with a shield."""
    plain = setup()
    ev = simulate(plain, EVENTS, 1)
    out = {"plain": simulate_gammas(plain, EVENTS, seed=1, particle_events=ev, gammas_per_event=PER_EVENT)}
    out["addback"] = simulate_gammas(setup(addback=True), EVENTS, seed=1, particle_events=ev,
                                     gammas_per_event=PER_EVENT)
    out["shield"] = simulate_gammas(setup(shield="BGO"), EVENTS, seed=1, particle_events=ev,
                                    gammas_per_event=PER_EVENT)
    return out


def test_the_add_back_factor_at_1332_kev_is_the_value_put_in_and_within_the_published_range():
    plain, on = setup(), setup(addback=True)
    r0 = response.Response(plain, plain.gamma_detectors[0])
    for factor, exp in ((response.ADDBACK_FACTOR, on), (1.4, setup(addback=True, addback_factor=1.4))):
        r = response.Response(exp, exp.gamma_detectors[0])
        assert r.peak_efficiency(E_REF) / r0.peak_efficiency(E_REF) == pytest.approx(factor)
        assert r.peak_efficiency(E_REF, 2) / r0.peak_efficiency(E_REF, 2) == pytest.approx(factor)
        assert r.peak_to_total(E_REF) / r0.peak_to_total(E_REF) == pytest.approx(factor)
        # The bare response of the same detector is the crystals as they are.
        assert response.Response(exp, exp.gamma_detectors[0], bare=True).peak_efficiency(E_REF) == pytest.approx(
            r0.peak_efficiency(E_REF))
        assert r.describe()["addback"] and r.describe()["addback_factor"] == pytest.approx(factor)
    # Published for the EUROGAM clover (Duchêne et al. 1999): about 1.5 at 1.33 MeV.
    assert 1.3 <= response.ADDBACK_FACTOR <= 1.7
    r = response.Response(on, on.gamma_detectors[0])
    assert r.addback_factor(0.3) < r.addback_factor(E_REF) < r.addback_factor(3.0)
    assert np.all(r.addback_share(np.array([0.3, 1.0, 3.0])) <= 1.0)
    # The planned coincidences follow: the rates use the efficiency with add-back.
    assert Rates(on).gamma_efficiency()[0] / Rates(plain).gamma_efficiency()[0] == pytest.approx(
        float(r.addback_factor(on.excitation.energy_mev)))
    # A measured curve is the detector as it runs, with add-back: the bare response divides the factor out.
    curve = setup(addback=True, efficiency_curve=[["122 keV", "2 %"], ["1332 keV", "0.4 %"]])
    rc = response.Response(curve, curve.gamma_detectors[0])
    assert rc.peak_efficiency(E_REF) == pytest.approx(0.004, rel=1e-3)
    assert response.Response(curve, curve.gamma_detectors[0], bare=True).peak_efficiency(E_REF) == pytest.approx(
        0.004 / response.ADDBACK_FACTOR, rel=1e-3)


def test_add_back_returns_compton_events_to_the_peak_in_the_simulated_run(runs):
    g0, g = runs["plain"], runs["addback"]
    r = response.Response(setup(addback=True), setup(addback=True).gamma_detectors[0])
    # The same γ rays: add-back only changes what the crystals record of them.
    assert len(g) == len(g0) and np.array_equal(g["event"], g0["event"])
    assert np.array_equal(g["energy_lab"], g0["energy_lab"])
    p, n_peak = peak(g)
    p0, _ = peak(g0)
    expected = float(r.addback_factor(np.average(g0["energy_lab"], weights=g0["weight"])))
    added = n_peak - peak(g0)[1]
    assert added > 300
    assert p / p0 == pytest.approx(expected, rel=max(3 / np.sqrt(added), 0.03)), (p / p0, expected)
    # What was added back came out of the continuum, and nothing else changed.
    moved = g["deposited"] != g0["deposited"]
    assert moved.sum() == added
    assert np.all(g["deposited"][moved] == g["energy_lab"][moved])
    assert np.all(g0["deposited"][moved] < g0["energy_lab"][moved])
    assert np.array_equal(g["measured"][~moved], g0["measured"][~moved])
    assert continuum(g)[0] == pytest.approx(continuum(g0)[0] - float(g0["weight"][moved & g0["counted"]].sum()))
    # The γ ray goes to the crystal with the larger deposit: a neighbour of the first one when it changes.
    switched = g["crystal"] != g0["crystal"]
    assert switched.any() and np.all(moved[switched])
    assert np.all((g["crystal"][switched] - g0["crystal"][switched]) % 4 != 2)
    assert np.all(g0["deposited"][switched] < g0["energy_lab"][switched] / 2)
    # ... and the Doppler correction is made from that crystal.
    assert not np.allclose(g["corrected_recoil"][switched], g0["corrected_recoil"][switched])
    assert np.allclose(g["corrected_recoil"][~moved], g0["corrected_recoil"][~moved])
    assert "Add-back returned" in g.notes[-1]


def test_suppression_lowers_the_continuum_by_the_factor_and_leaves_the_peak(runs):
    g0, g = runs["plain"], runs["shield"]
    exp = setup(shield="BGO")
    r = response.Response(exp, exp.gamma_detectors[0])
    assert response.SUPPRESSION_FACTOR == 3.0 and r.suppression_factor(E_REF) == pytest.approx(3.0)
    # The peak is the same γ rays with the same weights; only the rest is thinned.
    assert peak(g) == peak(g0)
    assert len(g) < len(g0) and set(g["event"]) <= set(g0["event"])
    c, n = continuum(g)
    c0, n0 = continuum(g0)
    expected = float(r.suppression_factor(np.average(g0["energy_lab"], weights=g0["weight"])))
    assert c0 / c == pytest.approx(expected, rel=max(3 / np.sqrt(n), 0.03)), (c0 / c, expected)
    # In the response: the peak efficiency is unchanged, the continuum is divided by S, so P/T rises.
    r0 = response.Response(setup(), setup().gamma_detectors[0])
    assert r.peak_efficiency(E_REF) == pytest.approx(r0.peak_efficiency(E_REF))
    pt, pt0 = float(r.peak_to_total(E_REF)), float(r0.peak_to_total(E_REF))
    assert (1 - pt) / pt == pytest.approx((1 - pt0) / pt0 / 3.0)
    assert r.total_efficiency(E_REF) < r0.total_efficiency(E_REF)
    edges = np.linspace(0.0, 1.5, 1501)
    s, s0 = r.shape(E_REF, edges), r0.shape(E_REF, edges)
    low = edges[:-1] < 1.0
    # Per γ ray recorded, the shapes sum to 1; against the peak, the continuum is divided by S.
    assert s.sum() == pytest.approx(1.0, abs=2e-3) and s0.sum() == pytest.approx(1.0, abs=2e-3)
    assert s[low].sum() / pt == pytest.approx(s0[low].sum() / pt0 / 3.0, rel=0.02)
    assert r.describe()["shield"] == "BGO" and r.describe()["suppression_factor"] == pytest.approx(3.0)
    with_both = setup(addback=True, shield="BGO", suppression_factor=4.0)
    rb = response.Response(with_both, with_both.gamma_detectors[0])
    assert rb.peak_efficiency(E_REF) == pytest.approx(response.ADDBACK_FACTOR * r0.peak_efficiency(E_REF))
    assert rb.suppression_factor(E_REF) == pytest.approx(4.0)


def test_the_shield_is_a_blocking_volume_around_the_clover():
    plain, shielded = setup(), setup(shield="BGO")
    assert plain.gamma_detectors[0].housing()[1] == pytest.approx(101.0)
    assert shielded.gamma_detectors[0].housing()[1] == pytest.approx(101.0 + 2 * response.SHIELDS["BGO"])
    assert setup(shield="BGO", shield_thickness="30 mm").gamma_detectors[0].housing()[1] == pytest.approx(161.0)
    blocker = next(b for b in Array.from_experiment(shielded).blockers if "(housing)" in b.name)
    assert blocker.width == pytest.approx(151.0)
    parts = [p for s in solids(shielded) if s.key == "gamma:0" for p in s.parts]
    shield = next(p for p in parts if p.role == "shield")
    assert shield.hull and shield.w == pytest.approx(151.0) and shield.label == "BGO shield"
    assert not any(p.role == "shield" for s in solids(plain) for p in s.parts)


def test_with_both_switched_off_the_results_equal_those_without_the_fields(runs):
    g0 = runs["plain"]
    off = setup(addback=False)
    ev = simulate(off, EVENTS, 1)
    g = simulate_gammas(off, EVENTS, seed=1, particle_events=ev, gammas_per_event=PER_EVENT)
    for key in g.columns:
        assert np.array_equal(g[key], g0[key]), key
    # A run with the switches on, made plain, is the bare run too.
    both = setup(addback=True, shield="BGO")
    g2 = simulate_gammas(both, EVENTS, seed=1, particle_events=ev, gammas_per_event=PER_EVENT, plain=True)
    assert np.array_equal(g2["measured"], g0["measured"]) and np.array_equal(g2["crystal"], g0["crystal"])
    text = setup().to_toml()
    assert "addback" not in text and "shield" not in text


def test_the_setup_checks_and_the_round_trip():
    with pytest.raises(SetupError, match="add-back needs a clover"):
        d = Experiment.example("coulex_ni58").to_dict()
        d["gamma_detectors"][0]["addback"] = True
        Experiment.from_dict(d)
    with pytest.raises(SetupError, match="must be true or false"):
        setup(addback="yes")
    with pytest.raises(SetupError, match="shield 'NaI' is not known"):
        setup(shield="NaI")
    with pytest.raises(SetupError, match="must be at least 1"):
        setup(shield="BGO", suppression_factor=0.5)
    with pytest.raises(SetupError, match="must be a number without a unit"):
        setup(addback=True, addback_factor="1.5")
    with pytest.raises(SetupError, match="a shield needs 'housing_side'"):
        d = Experiment.example("coulex_ni58").to_dict()
        d["gamma_detectors"][0]["shield"] = "BGO"
        Experiment.from_dict(d)
    exp = setup(addback=True, addback_factor=1.45, shield="BGO", shield_thickness="20 mm", suppression_factor=3.5)
    again = Experiment.from_toml(exp.to_toml())
    gd = again.gamma_detectors[0]
    assert gd.addback is True and gd.addback_factor == 1.45 and gd.shield == "BGO"
    assert gd.shield_mm() == 20.0 and gd.suppression_factor == 3.5
    for f in ("addback", "addback_factor", "shield", "shield_thickness", "suppression_factor"):
        assert guide.help_for("gamma", f) is not None


def test_the_planner_and_the_figure_show_the_spectrum_with_and_without():
    pytest.importorskip("plotly")
    from physim.nuclear import app

    p = Planner(setup(addback=True, shield="BGO"))
    assert p.gamma_modes() == {"Clo": {"addback": True, "shield": "BGO"}}
    g = p.gamma_spectra(100_000, 1)
    s = g["spectra"]["Clo"]
    assert "plain" in s and len(s["plain"]["measured"]) == len(s["measured"])
    assert s["plain"]["measured"].sum() > s["measured"].sum()
    fig = app.figure_gamma_spectra(p, "Clo", 100_000, 1)
    names = [t.name for t in fig.data]
    assert any("without add-back and suppression" in n for n in names)
    fig.to_json()
    assert "plain" not in Planner.example("coulex_ni58").gamma_spectra(100_000, 1)["spectra"]["Ge90"]
    assert app._crystals_of({"model": "clover"}) == 4 and app._crystals_of({"radius": "35 mm"}) is None
    assert app._value("addback", "yes") is True and app._value("addback_factor", "1.4") == 1.4
