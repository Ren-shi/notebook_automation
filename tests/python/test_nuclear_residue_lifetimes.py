"""Lifetimes of fusion-evaporation residues (lifetime wave 2B): level lifetimes and side-feeding times in the
residue cascade, residue transport through the plunger or the backing, and the predicted plunger components and
DSAM line shapes. Reference cases: the 50Cr three-foil plunger (Beckers et al., NIM A 1042, 167418 (2022)) and the
143Sm DSAM measurement (Rajbanshi et al., PLB 782, 143 (2018)), as the examples ``fusion_cr50_plunger`` and
``fusion_sm143_dsam`` set them up."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment, lifetimes
from physim.nuclear.experiment import SetupError
from physim.nuclear.levels import LN2
from physim.nuclear.residues import (SIDE_FEEDING_LEVELS, SIDE_FEEDING_Q_B, ChannelTable, ResidueCascade,
                                     dsam_lineshapes, lifetime_cascade, plunger_components, plunger_spectra,
                                     residue_scheme, simulate_residue_gammas)


@pytest.fixture(scope="module")
def cr50():
    return Experiment.example("fusion_cr50_plunger")


def _with(exp, **reaction):
    d = exp.to_dict()
    d["reaction"].update(reaction)
    return Experiment.from_dict(d)


# ---------------------------------------------------------------------------------------------------------------
# The two [reaction] fields


def test_side_feeding_fields_default_validate_and_round_trip(cr50):
    assert cr50.fusion.side_feeding_q is None and cr50.fusion.side_feeding_levels is None
    table = ChannelTable.from_experiment(cr50)
    assert (table.side_feeding_q_b, table.side_feeding_levels) == (SIDE_FEEDING_Q_B, SIDE_FEEDING_LEVELS) == (6.0, 5)
    exp = _with(cr50, side_feeding_q="3.5 b", side_feeding_levels=8)
    back = Experiment.from_dict(exp.to_dict())
    assert back.to_dict()["reaction"]["side_feeding_q"] == "3.5 b"
    assert back.to_dict()["reaction"]["side_feeding_levels"] == 8
    table = ChannelTable.from_experiment(back)
    assert (table.side_feeding_q_b, table.side_feeding_levels) == (3.5, 8)
    for bad, word in ((dict(side_feeding_q="6 MeV"), "side_feeding_q"), (dict(side_feeding_q="0 b"), "side_feeding_q"),
                      (dict(side_feeding_levels=2.5), "side_feeding_levels"),
                      (dict(side_feeding_levels=-1), "side_feeding_levels")):
        with pytest.raises(SetupError, match=word):
            _with(cr50, **bad)
    d = Experiment.example("coulex_ru98_plunger").to_dict()
    d["reaction"]["side_feeding_q"] = "6 b"
    with pytest.raises(SetupError, match="fusion-evaporation"):
        Experiment.from_dict(d)


# ---------------------------------------------------------------------------------------------------------------
# τ per level and the side-feeding band


def test_level_taus_come_from_the_scheme_and_unknown_ones_are_prompt(cr50):
    scheme = residue_scheme(cr50, "50Cr")
    cas = ResidueCascade(scheme)
    taus = cas.taus_ps
    assert taus[1] == pytest.approx(scheme.levels[1].half_life.value / LN2 * 1e12)   # 2+, ENSDF
    assert cas.tau_sources[1] == "half-life (ensdf)"
    assert math.isinf(taus[cas.ground])
    for i in cas.prompt_levels:            # γ rays but no τ: nan, never a made-up number
        assert np.isnan(taus[i]) and cas.tau_sources[i] is None and i in cas.branches
    assert set(cas.prompt_levels) == {i for i in cas.branches if cas.level_names[i]
                                      not in lifetimes.level_taus(scheme)}


def test_side_feeding_band_is_the_rotor_of_lineshape(cr50):
    scheme = residue_scheme(cr50, "50Cr")
    cas = ResidueCascade(scheme, side_feeding_q_b=6.0, side_feeding_levels=5)
    i = 1                                                       # 2+ at 783 keV, decays to the 0+ ground state
    band = cas.side_bands[i]
    assert [s for s, _e, _t in band] == [4, 6, 8, 10, 12]
    e_2 = next(t.energy.value for t in scheme.transitions if (t.initial, t.final) == (1, 0))
    a = e_2 / 6.0                                               # ħ²/2𝒥 from E(2+ → 0+) = 6A
    spin, e, tau = band[0]
    assert e == pytest.approx(a * 14)                            # E(4 → 2) = A (4·4 − 2)
    b_e2 = 5 / (16 * math.pi) * (6.0 * 100) ** 2 * 3 * 4 * 3 / (2 * 7 * 9)   # e² fm⁴
    assert tau == pytest.approx(lifetimes.tau_from_b(e, "E2", b_e2))
    assert cas.side_feeding_ps[i] == pytest.approx(sum(t for *_x, t in band))
    # A longer band takes longer, a smaller quadrupole moment longer still (as 1/Q²), and no band is prompt.
    longer = ResidueCascade(scheme, side_feeding_levels=10).side_feeding_ps[i]
    weaker = ResidueCascade(scheme, side_feeding_q_b=3.0).side_feeding_ps[i]
    assert longer > cas.side_feeding_ps[i]
    assert weaker == pytest.approx(4 * cas.side_feeding_ps[i])
    assert ResidueCascade(scheme, side_feeding_levels=0).side_feeding_ps[i] == 0.0
    assert np.isnan(cas.side_feeding_ps[cas.ground])


def test_timed_cascade_keeps_the_gamma_rays_and_adds_the_decay_times(cr50):
    scheme = residue_scheme(cr50, "50Cr")
    cas = ResidueCascade(scheme, side_feeding_levels=0)
    spins = np.full(20000, 4.0)                                   # enter at the 4+
    plain = cas.run(spins, np.random.default_rng(3))
    timed = cas.run(spins, np.random.default_rng(3), np.random.default_rng(9))
    for k in plain:
        assert np.array_equal(plain[k], timed[k])
    t = timed["time_ps"]
    four = (timed["initial"] == 2) & (timed["final"] == 1)
    two = (timed["initial"] == 1) & (timed["final"] == 0)
    tau4, tau2 = cas.taus_ps[2], cas.taus_ps[1]
    entered_at_4 = np.isin(timed["source"], timed["source"][four])
    assert t[four].mean() == pytest.approx(tau4, rel=0.05)
    assert t[two & entered_at_4].mean() == pytest.approx(tau4 + tau2, rel=0.05)


# ---------------------------------------------------------------------------------------------------------------
# Without a lifetime section nothing changes


def test_without_a_lifetime_section_the_gamma_rays_are_prompt(cr50):
    d = cr50.to_dict()
    d.pop("lifetime")
    exp = Experiment.from_dict(d)
    g = simulate_residue_gammas(exp, events=2000, seed=2)
    assert np.all(g["distance_um"] == 0) and not g["stopped"].any()
    assert np.array_equal(g["beta"], g.events["beta"][g["event"]])
    assert "_segment" not in g.columns
    with pytest.raises(ValueError, match="plunger"):
        simulate_residue_gammas(exp, events=100, distance_um=50.0)


# ---------------------------------------------------------------------------------------------------------------
# fusion_cr50_plunger: the three-foil plunger


DISTANCES = (17.0, 62.0, 200.0, 550.0)


@pytest.fixture(scope="module")
def cr50_runs(cr50):
    return {d: simulate_residue_gammas(cr50, events=60_000, seed=5, distance_um=d) for d in DISTANCES}


@pytest.fixture(scope="module")
def cr50_components(cr50):
    return plunger_components(cr50, DISTANCES, events=40_000, seed=11)


def _shares(g, initial, final, channel):
    """MC shares of one line's γ rays emitted before the degrader, between the foils and in the stopper (by the
    segment of rates.recoil_stack the residue was in), and their number."""
    m = (g["initial"] == initial) & (g["final"] == final) & (g["channel"] == channel)
    seg = g["_segment"][m]
    n = int(m.sum())
    return np.array([np.sum(seg <= 1), np.sum((seg == 2) | (seg == 3)), np.sum(seg >= 4)]) / n, n


def test_cr50_unshifted_share_of_the_2plus_follows_the_cascade(cr50, cr50_runs, cr50_components):
    comp = cr50_components
    assert comp["residue"] == "50Cr" and comp["kind"] == "differential"
    assert 0.035 < comp["beta_shifted"] < 0.045 and 0.01 < comp["beta_degraded"] < comp["beta_shifted"]
    channel = ChannelTable.from_experiment(cr50).index("4p2n")
    tau2 = comp["taus_ps"][1]
    beyond_single = 0.0
    for k, d in enumerate(DISTANCES):
        mc, n = _shares(cr50_runs[d], 1, 0, channel)
        predicted = np.array([comp["shifted"][1, k], comp["degraded"][1, k], comp["unshifted"][1, k]])
        sigma = np.sqrt(predicted * (1 - predicted) / n) + 0.005
        assert np.all(np.abs(mc - predicted) < 4 * sigma), (d, mc, predicted, n)
        # The 2+ is fed from the 4+ (τ 3.2 ps) and above: more of it survives to the degrader than a lone level.
        alone = math.exp(-comp["flight_ps"][k] / tau2)
        beyond_single = max(beyond_single, (mc[1] + mc[2] - alone) / math.sqrt(alone * (1 - alone) / n))
    assert beyond_single > 4
    # The cascade's own numbers: the unshifted share is lifetimes.unshifted at the residues' flight times, and it
    # falls with distance while the shifted share rises.
    assert np.all(np.diff(comp["unshifted"][1]) < 0) and np.all(np.diff(comp["shifted"][1]) > 0)
    cas = lifetime_cascade(cr50)
    one = lifetimes.unshifted(cas["taus_ps"], cas["feeding"], cas["populations"], comp["flight_ps"])[1]
    assert np.allclose(comp["shifted"][1], 1 - one, atol=0.03)


def test_cr50_three_components_appear_for_the_differential_plunger(cr50, cr50_runs):
    g = cr50_runs[62.0]
    channel = ChannelTable.from_experiment(cr50).index("4p2n")
    m = (g["initial"] == 1) & (g["final"] == 0) & (g["channel"] == channel) & g["counted"]
    seg = g["_segment"][m]
    theta = g["theta"][m]
    e = 1e3 * g["energy_lab"][m]
    groups = (seg <= 1, (seg == 2) | (seg == 3), seg >= 4)
    for backward in (True, False):
        ring = theta > 90 if backward else theta < 90
        means = [e[ring & grp].mean() for grp in groups]
        assert all(np.sum(ring & grp) > 30 for grp in groups)
        assert means[2] == pytest.approx(783.3, abs=1.5)           # unshifted: at rest in the stopper
        if backward:
            assert means[0] < means[1] < means[2] - 5               # fastest first: shifted, degraded, unshifted
        else:
            assert means[0] > means[1] > means[2] + 5
    assert g["stopped"][m].mean() > 0.1 and set(np.unique(g["distance_um"])) == {62.0}

    sp = plunger_spectra(cr50, [17.0, 200.0], events=20_000, seed=3)
    assert sp["theta_deg"] == [45.0, 142.3]
    line = next(x for x in sp["lines"] if (x["initial"], x["final"]) == (1, 0))
    centres = 0.5 * (line["edges_kev"][1:] + line["edges_kev"][:-1])
    for name in ("shifted", "degraded", "unshifted"):
        assert line[name].sum() > 0
    peak = {name: centres[np.argmax(line[name][0, 1])] for name in ("shifted", "degraded", "unshifted")}
    assert peak["shifted"] < peak["degraded"] < peak["unshifted"]   # at 142.3 deg
    assert peak["unshifted"] == pytest.approx(783.3, abs=1.5)
    # With distance, intensity moves from the degraded and unshifted peaks into the shifted one.
    assert line["shifted"][1].sum() > line["shifted"][0].sum()
    assert line["unshifted"][1].sum() < line["unshifted"][0].sum()
    assert np.allclose(line["total"].sum(axis=-1), line["per_residue"], rtol=0.02)


# ---------------------------------------------------------------------------------------------------------------
# fusion_sm143_dsam: line shapes in the gold backing


def _sm143(taus_ps, **reaction):
    """The 143Sm example with a small dipole band (illustrative, as the published band's 0.4–1 ps lifetimes): the
    143Sm high-spin scheme is not in the setup, so the test builds one. Each cascade enters at the top."""
    d = Experiment.example("fusion_sm143_dsam").to_dict()
    d["reaction"].update(reaction)
    for ch in d["reaction"]["channels"]:
        if ch["residue"] == "143Sm":
            ch["entry_spin"] = 21.5
    levels = [{"energy": "0 keV", "jpi": "13/2-"}]
    transitions = []
    for k, spin in enumerate((15, 17, 19, 21, 23), start=1):
        levels.append({"energy": f"{450.0 * k + 25.0 * k * k:g} keV", "jpi": f"{spin}/2-",
                       "half_life": f"{taus_ps[k - 1] * LN2:.6g} ps"})
        transitions.append({"from": k, "to": k - 1, "energy": f"{450.0 + 25.0 * (2 * k - 1):g} keV",
                            "intensity": 100, "multipolarity": "M1"})
    d.setdefault("levels", {})["143Sm"] = {"nuclide": "143Sm", "level": levels, "transition": transitions}
    return Experiment.from_dict(d)


def _width(line, ring):
    """Standard deviation (keV) of a line shape without resolution."""
    centres = 0.5 * (line["edges_kev"][1:] + line["edges_kev"][:-1])
    w = line["shapes"][ring] / line["shapes"][ring].sum()
    mean = w @ centres
    return math.sqrt(w @ (centres - mean) ** 2)


def test_sm143_lineshapes_broaden_and_shift_as_tau_shortens():
    out = {}
    for tau in (3.0, 0.5, 0.15):
        exp = _sm143([tau, 0.05, 0.05, 0.05, 0.05], side_feeding_levels=0)
        out[tau] = dsam_lineshapes(exp, theta_deg=[148.0, 90.0], n=3000, events=6000, resolution=False)
    first = {tau: next(x for x in r["lines"] if (x["initial"], x["final"]) == (1, 0)) for tau, r in out.items()}
    assert out[3.0]["residue"] == "143Sm" and 0.012 < out[3.0]["beta0"] < 0.02
    f = [first[tau]["attenuation"][0] for tau in (3.0, 0.5, 0.15)]
    assert 0 < f[0] < f[1] < f[2] < 1                               # F(τ) grows as τ shortens
    shift = [first[tau]["centroid_kev"][0] - first[tau]["energy_kev"] for tau in (3.0, 0.5, 0.15)]
    assert shift[0] > shift[1] > shift[2]                           # backward: shifted down, more for short τ
    assert _width(first[0.5], 0) > 1.5 * _width(first[3.0], 0)       # the line broadens out of the stopped peak
    # At 90 deg the line hardly moves (the second-order shift and the scattering only).
    assert abs(first[0.15]["centroid_kev"][1] - first[0.15]["energy_kev"]) < 0.25 * abs(shift[2])


def test_sm143_side_feeding_delays_the_line_shapes():
    taus = [0.4, 0.3, 0.3, 0.3, 0.3]
    prompt = dsam_lineshapes(_sm143(taus, side_feeding_levels=0), theta_deg=[148.0], n=3000, events=6000,
                             resolution=False)
    slow = dsam_lineshapes(_sm143(taus, side_feeding_q="0.5 b", side_feeding_levels=10), theta_deg=[148.0],
                           n=3000, events=6000, resolution=False)
    top = {r: next(x for x in r_["lines"] if x["initial"] == 5) for r, r_ in (("prompt", prompt), ("slow", slow))}
    assert slow["side_feeding_ps"][5] > 0.3 and prompt["side_feeding_ps"][5] == 0.0
    assert top["slow"]["attenuation"][0] < top["prompt"]["attenuation"][0] - 0.05
    # The delay carries down the cascade.
    low = {r: next(x for x in r_["lines"] if x["initial"] == 1) for r, r_ in (("prompt", prompt), ("slow", slow))}
    assert low["slow"]["attenuation"][0] < low["prompt"]["attenuation"][0]


def test_sm143_events_slow_and_stop_in_the_backing():
    exp = _sm143([0.5, 0.05, 0.05, 0.05, 0.05], side_feeding_levels=0)
    g = simulate_residue_gammas(exp, events=8000, seed=1, gammas_per_event=1)
    channel = ChannelTable.from_experiment(exp).index("5n")
    m = (g["channel"] == channel) & (g["initial"] == 1)
    assert m.sum() > 200
    stopped = g["stopped"][m]
    assert 0.1 < stopped.mean() < 0.9
    assert np.all(g["beta"][m][stopped] == 0)
    assert g["beta"][m].mean() < 0.8 * g.events["beta"][g["event"][m]].mean()
    assert np.all(g["distance_um"] == 0)
    with pytest.raises(ValueError, match="DSAM"):
        simulate_residue_gammas(exp, events=100, distance_um=10.0)
    with pytest.raises(ValueError, match="plunger"):
        plunger_components(exp)
