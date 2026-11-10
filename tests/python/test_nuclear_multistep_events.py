"""The simulated events take the level populations from the coupled equations: levels reached in two steps (the 4⁺
and the second 0⁺ state of ⁶²Ni) are populated and their γ rays are in the spectra with the strength the multi-step
solution gives; first order stays a switch; the rates and the events share one excitation probability."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import dataviews as dv
from physim.nuclear.coupled import CoupledChannels
from physim.nuclear.events import generator_config, simulate
from physim.nuclear.gamma import cascade_excitation, start_levels, unreached_levels
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.multistep import Multistep
from physim.nuclear.orientation import Excitation
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates, channels, excitation_probability, stack

#: Levels of ⁶²Ni in the example's scheme: the first 2⁺, the second 0⁺ (2048.55 keV) and the 4⁺ (2336.46 keV).
TWO, ZERO, FOUR = 1, 2, 4
LINE_4_2, LINE_0_2, LINE_2_0 = 1.16355, 0.87564, 1.17291  # MeV: differences of the level energies


def ni62(order=None, emission=None) -> Planner:
    p = Planner.example("coulex_ni62_pt194")
    d = p.draft
    if order is not None:
        d["reaction"]["order"] = order
    if emission is not None:
        d["reaction"]["emission"] = emission
    assert p._apply(d), p.problems
    return p


def test_many_angles_solved_together_are_the_single_solutions():
    exp = ni62().experiment
    first = cascade_excitation(exp, order="first")
    cc = CoupledChannels(first.beam, first.target, first.beam_energy, first.scheme, excite="projectile",
                         tolerance=1e-3)
    angles = [4.0, 30.0, 91.0, 150.0, 180.0]
    for theta, a in zip(angles, cc.solve_many(angles)):
        assert np.abs(a - cc.solve(theta)).max() < 1e-12, theta
    grid, states = cc.table(30.0)
    for theta, state in zip(grid, states):
        alone = cc.at(float(theta))
        assert state.direct == pytest.approx(alone.direct, rel=1e-9, abs=1e-15)
        assert state.population(TWO) == pytest.approx(alone.population(TWO), rel=1e-9)
    # On the table the events use, the forward angles where less than a millionth is excited stay in first order
    # (they take the longest to solve): nothing there differs by more than a part in a hundred million.
    grid, states = cc.table(2.0)
    weak = [s for s in states if sum(s.direct.values()) < 1e-6]
    assert 4 <= len(weak) <= 8 and all(s.direct[FOUR] == 0.0 for s in weak)
    assert weak[-1].direct[TWO] == pytest.approx(first.probability(TWO, weak[-1].theta_cm), rel=1e-9)
    assert weak[-1].direct[TWO] == pytest.approx(cc.probability(TWO, weak[-1].theta_cm), abs=1e-8)
    assert states[len(weak)].direct[FOUR] > 0.0


def test_all_orders_reach_the_levels_first_order_cannot_and_first_order_is_a_switch():
    exp = ni62().experiment
    cex = cascade_excitation(exp)
    assert isinstance(cex, CoupledChannels)
    assert cascade_excitation(exp) is cex  # solved once, shared by the events, the rates and the γ rays
    assert start_levels(cex) == [1, 2, 3, 4, 5] and unreached_levels(exp) == []
    back = cex.near(150.0)
    assert back.direct[FOUR] == pytest.approx(5.8e-3, rel=0.05) and back.direct[ZERO] == pytest.approx(8.1e-3, rel=0.05)
    # The two-step excitation leaves the 4⁺ state oriented: its γ rays are not isotropic.
    assert abs(back.coefficients(FOUR, TWO)[(2, 0)]) > 0.1

    first = ni62(order="first").experiment
    one = cascade_excitation(first)
    assert type(one) is Excitation and start_levels(one) == [1, 3]
    assert unreached_levels(first) == [2, 4, 5]
    g = simulate_gammas(first, 60_000, seed=4)
    assert set(np.unique(g["initial"]).tolist()) == {1, 3}
    p = ni62()
    assert not p.set("reaction", "order", "second") and p.set("reaction", "order", "first")


def test_the_rates_and_the_generator_share_the_excitation_probability():
    for order in ("all", "first"):
        exp = ni62(order=order).experiment
        layers = stack(exp)
        ch = next(c for c in channels(layers, exp) if c.excitation is not None)
        p = excitation_probability(exp, ch, layers)
        grid, states = cascade_excitation(exp).table()
        assert float(p.probability(grid[60])) == pytest.approx(sum(states[60].direct.values()))
        config, used = generator_config(exp)
        entry = next(c for c, k in zip(config["channels"], config["channel_of"]) if used[k].excitation is not None)
        assert np.allclose(entry["p_table"], p.probability(np.linspace(0.0, 180.0, 721)))
    # To all orders more is excited than in first order: the 2⁺ state is fed from the levels above it.
    all_orders, first = Rates(ni62().experiment), Rates(ni62(order="first").experiment)
    ratio = all_orders.rate("CD", what="excitations") / first.rate("CD", what="excitations")
    assert 1.02 < ratio < 1.05
    ev = simulate(ni62().experiment, 300_000, seed=2)
    label = next(c for c in ev.channels if "excited" in c)
    for det in ev.detectors:
        rate, err = ev.rate(det, counted=False, channel=label)
        expected = all_orders.rate(det, counted=False, what="excitations")
        assert abs(rate - expected) < 3 * err + 0.005 * expected, (det, rate, err, expected)


#: Half-width of the window a line is counted in, MeV: the Doppler-corrected lines are 12 to 20 keV wide (FWHM).
HALF = 0.025


def _line(g, initial: int, final: int, energy: float) -> np.ndarray:
    """The Doppler-corrected energies of the counted γ rays of one transition (the run knows which it was) that
    lie in the line's window."""
    x = g["corrected_projectile"]
    own = (g["initial"] == initial) & (g["final"] == final) & g["counted"]
    return x[own & (np.abs(x - energy) < HALF)]


def _detected(g, energy: float, half: float = HALF) -> float:
    """The share of the γ rays of ``energy`` that a crystal records within ``half`` of the full energy, summed over
    the crystals: the full-energy peak and the part of the continuum just below it."""
    fine = np.arange(0.0, 3.0, 0.0005)
    centres = (fine[:-1] + fine[1:]) / 2
    return sum(float(c.response.total_efficiency(energy, c.element))
               * float(c.response.shape(energy, fine)[np.abs(centres - energy) < half].sum()) for c in g.crystals)


def test_a_run_of_62ni_on_194pt_shows_the_4plus_and_0plus_lines_of_the_multistep_solution(tmp_path):
    """A beam run, every event simulated, with isotropic emission so that the efficiencies alone set the counts:
    the γ rays of 4⁺ → 2⁺ (1163.5 keV) and 0⁺ → 2⁺ (875.6 keV) against the multi-step yields of
    :class:`~physim.nuclear.multistep.Multistep` (its own quadrature over the detectors, two beam energies through
    the target) times what the crystals record in the line's window."""
    p = ni62(emission="isotropic")
    exp = p.experiment
    p.create_experiment(root=tmp_path)
    run = p.start_run("beam", duration="40 min", budget="40 min", seed=11)
    assert run.scale == 1.0
    g = run.gammas()
    seconds = run.real_s * g.live_fraction
    yields = Multistep(exp, exp.levels["beam"], "beam").yields()

    def predicted(transition, energy):
        return yields.total(transition) * _detected(g, energy) * seconds

    lines = {(TWO, 0): LINE_2_0, (FOUR, TWO): LINE_4_2, (ZERO, TWO): LINE_0_2}
    found = {t: _line(g, *t, e) for t, e in lines.items()}
    n = {t: len(x) for t, x in found.items()}
    assert n[(FOUR, TWO)] > 300 and n[(ZERO, TWO)] > 400, n
    for t, e in lines.items():
        # The Doppler correction for the projectile puts each line at its own energy: 4⁺ → 2⁺ at 1163.5 keV, not
        # at the 2⁺ line's 1172.9 keV.
        assert np.median(found[t]) == pytest.approx(e, abs=2e-3), t
        # In absolute counts; 3% for what the prediction leaves out (two γ rays of a cascade in one crystal, the
        # γ rays thrown along the emitter's velocity, the coarser quadrature of the multi-step yields).
        expected = predicted(t, e)
        assert abs(n[t] - expected) < 3 * math.sqrt(expected) + 0.03 * expected, (t, n[t], expected)
    # Relative to the 2⁺ line of the same run, where those cancel.
    for t in ((FOUR, TWO), (ZERO, TWO)):
        expected = predicted(t, lines[t]) / predicted((TWO, 0), LINE_2_0)
        measured = n[t] / n[(TWO, 0)]
        assert abs(measured - expected) < 3 * expected / math.sqrt(n[t]) + 0.01 * expected, (t, measured, expected)

    # The spectrum as the Data tab draws it. The 0⁺ → 2⁺ line stands clear of the others on the 2⁺ line's Compton
    # continuum (curved here, so a parabola through the bins on both sides): its counts above that are the
    # multi-step solution's.
    total = None
    for name in g.detector_names():
        view = dv.gamma_spectrum(run, name, correction="projectile", randoms="none", bins=1500, range=(0.0, 3.0))
        total = view["counts"] if total is None else total + view["counts"]
    centre = (view["edges"][:-1] + view["edges"][1:]) / 2
    off = centre - LINE_0_2
    peak = np.abs(off) < 0.020
    side = ((off > -0.100) & (off < -0.024)) | ((off > 0.024) & (off < 0.060))
    sigma = np.sqrt(np.maximum(total[side], 1.0))
    design = np.vander(off[side], 3) / sigma[:, None]
    coefficients, *_ = np.linalg.lstsq(design, total[side] / sigma, rcond=None)
    under = np.vander(off[peak], 3).sum(axis=0)
    net = total[peak].sum() - under @ coefficients
    err = math.sqrt(total[peak].sum() + under @ np.linalg.inv(design.T @ design) @ under)
    assert net > 5 * err, (net, err)
    expected = yields.total((ZERO, TWO)) * _detected(g, LINE_0_2, 0.020) * seconds
    assert abs(net - expected) < 3 * err, (net, err, expected)
    # The 4⁺ → 2⁺ line is 9.4 keV below the 2⁺ line, whose corrected width here is some 15 keV: its counts are
    # in the spectrum, inside the 2⁺ peak, and they are all there between 1.13 and 1.20 MeV.
    both = (centre > 1.130) & (centre < 1.200)
    x = g["corrected_projectile"]
    in_window = g["counted"] & (x > 1.130) & (x < 1.200)
    assert total[both].sum() == pytest.approx(in_window.sum(), abs=3)
    assert (in_window & (g["initial"] == FOUR)).sum() > 300

    # Predicted against measured: the peak the plan predicts is the 2⁺ state's own.
    for row in p.compare()["rows"]:
        if row["quantity"] in ("Rate", "Coincidences in the peak"):
            assert row["agree"], row


def test_a_setup_without_a_scheme_is_as_before():
    exp = Experiment.example("coulex_ni58")
    assert cascade_excitation(exp) is None
    layers = stack(exp)
    ch = next(c for c in channels(layers, exp) if c.excitation is not None)
    from physim.nuclear.coulex import Coulex

    assert isinstance(excitation_probability(exp, ch, layers), Coulex)
