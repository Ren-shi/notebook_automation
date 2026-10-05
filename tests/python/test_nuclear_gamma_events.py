"""Particle–γ events, Doppler correction, coincidences and background (backlog 55)."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.events import simulate
from physim.nuclear.experiment import SetupError
from physim.nuclear.gamma import correlation_table, doppler_table
from physim.nuclear.gamma_events import ROOM_LINES, simulate_gammas
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates

EVENTS = 600_000


@pytest.fixture(scope="module")
def ni():
    exp = Experiment.example("coulex_ni58")
    return exp, simulate_gammas(exp, EVENTS, seed=3)


def peak(g, crystal, detector, key, e0, half=0.03):
    """Weighted mean and standard deviation of the energies of the full-energy events within ±half MeV of e0,
    and the error of the mean (the Compton continuum under the peak is left out, which a fit would do)."""
    m = g.select(crystal, detector) & (g["deposited"] == g["energy_lab"])
    x, w = g[key][m], g["weight"][m]
    inside = np.abs(x - e0) < half
    x, w = x[inside], w[inside]
    mean = np.average(x, weights=w)
    sd = math.sqrt(np.average((x - mean) ** 2, weights=w))
    err = sd * math.sqrt(np.sum(w**2)) / np.sum(w)
    return mean, sd, err, len(x)


def test_the_corrected_peak_sits_at_the_transition_energy_for_the_target(ni):
    exp, g = ni
    e0 = exp.excitation.energy_mev
    for crystal in ("Ge90", "Ge45", "Ge135"):
        raw = peak(g, crystal, None, "measured", e0)
        good = peak(g, crystal, None, "corrected_recoil", e0)
        wrong = peak(g, crystal, None, "corrected_projectile", e0, half=0.06)
        assert good[3] > 60
        assert abs(good[0] - e0) < 3 * good[2] + 0.5e-3, (crystal, good)
        # The correction narrows the peak; correcting for the wrong nucleus smears it instead.
        assert good[1] < raw[1] and wrong[1] > good[1]


def test_the_corrected_peak_sits_at_the_transition_energy_for_the_projectile():
    d = Experiment.example("coulex_ni58").to_dict()
    d["beam"]["energy"] = "70 MeV"
    d["reaction"] = {"type": "coulex", "excite": "projectile", "energy": "6.917 MeV", "multipolarity": "E2",
                     "b_up": "0.004 e2b2"}
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, EVENTS, seed=5)
    assert g.emitter == "16O"
    e0 = 6.917
    for crystal in ("Ge90", "Ge45"):
        raw = peak(g, crystal, None, "measured", e0, half=0.25)
        good = peak(g, crystal, None, "corrected_projectile", e0, half=0.25)
        assert good[3] > 30 and abs(good[0] - e0) < 3 * good[2] + 2e-3, (crystal, good)
        assert good[1] < 0.5 * raw[1]


def test_the_raw_width_agrees_with_the_analytic_doppler_broadening():
    """The analytic table (item 43) follows the scattered beam into each particle detector and spreads the γ
    rays evenly over each crystal; the events include the recoils that reach the detector too, and emit with
    the angular correlation. So the comparison is made with isotropic emission, on the events with the scattered
    beam detected, and on the true γ-ray energy before the crystal's resolution."""
    d = Experiment.example("coulex_ni58").to_dict()
    d["reaction"]["emission"] = "isotropic"
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, EVENTS, seed=3)
    table = {(r["particle_detector"], r["gamma_detector"]): r for r in doppler_table(exp)}
    ejectile = ~g.events.columns["recoil"][g["particle"]]
    checked = 0
    for det in g.events.detectors:
        for crystal in ("Ge90", "Ge45", "Ge135", "Ge0"):
            m = g.select(crystal, det) & ejectile
            if m.sum() < 150:
                continue
            x, w = g["energy_lab"][m], g["weight"][m]
            mean = np.average(x, weights=w)
            sd = math.sqrt(np.average((x - mean) ** 2, weights=w))
            row = table[(det, crystal)]
            assert mean * 1e3 == pytest.approx(row["mean_kev"], abs=1.5), (det, crystal)
            assert sd * 1e3 * 2.3548 == pytest.approx(row["doppler_fwhm_kev"], rel=0.25), (det, crystal)
            checked += 1
    assert checked >= 4


def test_event_counts_agree_with_the_analytic_rates(ni):
    exp, g = ni
    r = Rates(exp)
    e0 = exp.excitation.energy_mev
    factors = {(x["particle_detector"], x["gamma_detector"]): x["detector_factor"] for x in correlation_table(exp)}
    pulls = []
    for det in g.events.detectors:
        exc = r.rate(det, what="excitations")
        for i, gd in enumerate(exp.gamma_detectors):
            resp = g.crystals[i].response
            expected = exc * float(resp.total_efficiency(e0)) * factors[(det, gd.name)]
            rate, err = g.rate(gd.name, det)
            # 5%: the analytic factor comes from a coarse quadrature at the mid-target energy.
            pulls.append((rate - expected) / math.hypot(err, 0.05 * expected))
            assert abs(pulls[-1]) < 3, (det, gd.name, rate, expected)
    assert abs(np.mean(pulls)) < 1.0
    # The particle side is the generator's: its rates agree with the analytic ones as before.
    for det in g.events.detectors:
        rate, err = g.events.rate(det)
        assert abs(rate - r.rate(det)) < 3 * err + 0.01 * r.rate(det)


def test_random_coincidences_follow_the_singles_rates_and_the_window():
    d = Experiment.example("coulex_ni58").to_dict()
    d["run"]["coincidence_window"] = "200 ns"
    d["run"]["room_background"] = "5 /s"
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, 100_000, seed=2)
    assert g.window_s == pytest.approx(200e-9)
    for det in ("CD", "DSSD-L"):
        singles = sum(g.singles_rate["Ge90"].values())
        assert g.random_rate("Ge90", det) == pytest.approx(2 * 200e-9 * g.particle_rate[det] * singles)
        assert g.particle_rate[det] == pytest.approx(Rates(exp).rate(det))
    # The random spectrum carries the expected number of counts, and the Poisson draw agrees within statistics.
    t = exp.run.beam_time_s if hasattr(exp.run, "beam_time_s") else 24 * 3600
    h_true, edges = g.spectrum("Ge90", "DSSD-L", randoms=False)
    h_all, _ = g.spectrum("Ge90", "DSSD-L", randoms=True)
    expected = g.random_rate("Ge90", "DSSD-L") * t * g.live_fraction
    drawn = h_all.sum() - h_true.sum()
    assert expected > 100 and abs(drawn - expected) < 4 * math.sqrt(expected)
    assert g.random_spectrum("Ge90", "DSSD-L", edges).sum() == pytest.approx(g.random_rate("Ge90", "DSSD-L"), rel=0.02)
    assert g.gamma_random_rate("Ge90", "Ge45") == pytest.approx(
        2 * 200e-9 * sum(g.singles_rate["Ge90"].values()) * sum(g.singles_rate["Ge45"].values()))


def test_room_background_and_extra_lines_appear_in_the_singles():
    d = Experiment.example("coulex_ni58").to_dict()
    d["run"]["room_background"] = "10 /s"
    d["run"]["extra_lines"] = [["1274.5 keV", "2 /s"]]
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, 50_000, seed=1)
    s, edges = g.singles_spectrum["Ge90"], g.singles_edges
    centres = (edges[:-1] + edges[1:]) / 2
    assert g.singles_rate["Ge90"]["background"] > 12.0          # the peaks plus their continua
    for e_kev in (1460.820, 2614.511, 1274.5):
        window = np.abs(centres - e_kev * 1e-3) < 0.004
        beside = (np.abs(centres - e_kev * 1e-3) > 0.006) & (np.abs(centres - e_kev * 1e-3) < 0.012)
        assert s[window].sum() > 5 * s[beside].sum()
    assert any(line[1] == "40K" for line in ROOM_LINES)
    with pytest.raises(SetupError, match="is not a rate"):
        Experiment.from_dict({**d, "run": {**d["run"], "extra_lines": [["1274 keV", "2 keV"]]}})


def test_dead_time_and_thresholds():
    d = Experiment.example("coulex_ni58").to_dict()
    d["run"]["dead_time"] = "20 us"
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, 50_000, seed=1)
    total = sum(g.particle_rate.values()) + sum(sum(v.values()) for v in g.singles_rate.values())
    assert g.live_fraction == pytest.approx(1 / (1 + 20e-6 * total))
    assert g.counts("Ge90") == pytest.approx(g.rate("Ge90")[0] * g.events.beam_time_s * g.live_fraction)
    d["run"].pop("dead_time")
    d["gamma_detectors"][0]["threshold"] = "1.2 MeV"
    g2 = simulate_gammas(Experiment.from_dict(d), 50_000, seed=1)
    m = g2.select("Ge90", counted=False)
    assert np.all(g2["measured"][m][g2["counted"][m]] >= 1.2)
    assert g2["counted"][m].sum() < m.sum()


def test_the_same_setup_and_seed_give_the_same_gamma_rays():
    exp = Experiment.example("coulex_ni58")
    a, b = simulate_gammas(exp, 60_000, seed=4), simulate_gammas(exp, 60_000, seed=4)
    assert len(a) == len(b) > 0
    for name, col in a.columns.items():
        assert np.array_equal(col, b.columns[name]), name
    c = simulate_gammas(exp, 60_000, seed=5)
    assert not np.array_equal(a["measured"], c["measured"][:len(a)]) or len(c) != len(a)
    # Built on given particle events, the γ rays belong to those events.
    ev = simulate(exp, 60_000, seed=4)
    d = simulate_gammas(exp, 60_000, seed=4, particle_events=ev)
    assert d.events is ev and np.array_equal(d["event"], a["event"])
    assert np.all(ev.columns["event"][d["particle"]] == d["event"])


def test_gamma_rays_leave_the_moving_nucleus_with_the_correlation(ni):
    exp, g = ni
    # The emitter is the ⁵⁸Ni recoil: slow, so the shift is small; forward γ rays are blue-shifted.
    forward = g["theta"] < 60
    assert g["energy_lab"][forward].mean() > g.energy_mev > g["energy_lab"][g["theta"] > 120].mean()
    assert 0.005 < g["beta"].mean() < 0.03
    # With the backward ring, the crystal at 90° to the beam sees far fewer γ rays than those at 45° and 135°.
    counts = {name: g.counts(name, "CD") for name in ("Ge90", "Ge45", "Ge135")}
    assert counts["Ge90"] < 0.5 * counts["Ge45"] and counts["Ge90"] < 0.5 * counts["Ge135"]


def test_the_planner_and_the_app_figures(ni):
    pytest.importorskip("plotly")
    from physim.nuclear import app

    p = Planner.example("coulex_ni58")
    g = p.gamma_spectra(150_000, 1)
    assert g["available"] and set(g["gamma_detectors"]) == {"Ge90", "Ge45", "Ge135", "Ge0"}
    s = g["spectra"]["Ge90"]
    assert len(s["measured"]) == len(s["recoil"]) == len(s["edges"]) - 1 and s["counts"] > 0
    assert set(g["coincidences"]) == {"true", "random", "gamma_gamma_random"}
    assert p.gamma_events(150_000, 1) is p.gamma_events(150_000, 1)
    app.figure_gamma_spectra(p, "Ge45", 150_000, 1).to_json()
    assert not Planner.example("alpha_on_gold").gamma_spectra()["available"]


def test_root_export_carries_the_gamma_rays(tmp_path):
    uproot = pytest.importorskip("uproot")
    from physim.nuclear.rootio import read_root, write_root

    exp = Experiment.example("coulex_ni58")
    ev = simulate(exp, 60_000, seed=1)
    g = simulate_gammas(exp, 60_000, seed=1, particle_events=ev)
    path = write_root(ev, tmp_path / "e.root", exp, gammas=g)
    back = read_root(path)
    assert set(back["gammas"]) == set(g.columns) and len(back["gammas"]["crystal"]) == len(g)
    assert back["info"]["crystals"] == g.crystal_names
    with uproot.open(path) as f:
        assert "gamma_Ge90" in f
