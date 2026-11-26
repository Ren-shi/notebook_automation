"""Deorientation in vacuum in a setup's γ rays: the setting, the plan's angular correlation, the events (each level's own
time in vacuum, down a cascade), stored runs, and the GOSIA input."""

from __future__ import annotations

import copy
import json
import math
import re

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import gamma_events as ge
from physim.nuclear.deorientation import Parameters, TwoState
from physim.nuclear.experiment import SetupError
from physim.nuclear.gamma import (DepthOrientation, Faded, cascade_excitation, correlation_table, deorientation_of,
                                  deorientation_parameters, other_excitation)
from physim.nuclear.gamma_events import _InVacuum, simulate_gammas
from physim.nuclear.lifetime_plan import deorientation_warning
from physim.nuclear.multistep import Multistep
from physim.nuclear.planner import Planner
from physim.nuclear.orientation import attenuated, tensor_basis, tensors_of
from physim.nuclear.runs import ExperimentFolder, RunData, physics_changes, taken_setup


def setup(name: str, how=None, g=None, **levels) -> Experiment:
    raw = Experiment.example(name).to_dict()
    for key, value in (("deorientation", how), ("deorientation_g", g)):
        if value is not None:
            raw["reaction"][key] = value
    for role, scheme in levels.items():
        raw.setdefault("levels", {})[role] = scheme
    return Experiment.from_dict(json.loads(json.dumps(raw)))


def ni58(half_life="20e-12 s", how=None, g=None) -> Experiment:
    """The template of one state, with a level scheme that gives it a half-life: its γ rays follow no cascade."""
    raw = Experiment.example("coulex_ni58").to_dict()
    level = {"energy": raw["reaction"]["energy"], "jpi": "2+"}
    if half_life:
        level["half_life"] = half_life
    scheme = {"nuclide": "58Ni", "level": [{"energy": "0 keV", "jpi": "0+"}, level],
              "transition": [{"from": 1, "to": 0, "energy": raw["reaction"]["energy"], "multipolarity": "E2"}],
              "matrix_element": [{"from": 0, "to": 1, "multipolarity": "E2", "value": "26.36 efm2"}]}
    return setup("coulex_ni58", how, g, target=scheme)


def factors(exp, **more) -> dict:
    return {(r["particle_detector"], r["gamma_detector"], r["crystal"]): r["factor"] for r in correlation_table(exp, **more)}


# -- the setting ----------------------------------------------------------------------------------------------------


def test_it_is_on_unless_the_setup_says_otherwise():
    exp = Experiment.example("coulex_ni62_pt194")
    assert exp.excitation.deorientation is None and deorientation_of(exp) == "vacuum"
    assert deorientation_parameters(exp) == deorientation_parameters(exp, "other") == Parameters()
    p = Planner(exp)
    assert p.set("reaction", "deorientation", "none") and deorientation_of(p.experiment) == "none"
    assert deorientation_parameters(p.experiment) is None
    assert Experiment.from_toml(p.experiment.to_toml()).excitation.deorientation == "none"
    assert not p.set("reaction", "deorientation", "on")
    assert p.problems == ["reaction: deorientation must be 'vacuum' (GOSIA's two-state model) or 'none', got 'on'"]
    assert p.set("reaction", "deorientation", "vacuum") and p.set("reaction", "deorientation_g", 0.17)
    # The g factor is that of the nucleus the reaction names, as the rest of the section.
    assert deorientation_parameters(p.experiment).g_factor == 0.17
    assert deorientation_parameters(p.experiment, "other").g_factor is None
    assert Experiment.from_toml(p.experiment.to_toml()).excitation.deorientation_g == 0.17
    for wrong in (7, -6.5):
        assert not p.set("reaction", "deorientation_g", wrong) and "g factor" in p.problems[0]
    raw = Experiment.example("alpha_on_gold").to_dict()
    raw["reaction"] = {"type": "elastic", "deorientation": "none"}
    with pytest.raises(SetupError, match="deorientation only apply to type"):
        Experiment.from_dict(raw)


def test_a_lifetime_measurement_has_only_its_own():
    exp = Experiment.example("coulex_ru98_plunger")
    assert exp.lifetime.deorientation is not None and deorientation_of(exp) == "none"
    raw = exp.to_dict()
    del raw["lifetime"]["deorientation"]
    raw["reaction"]["deorientation"] = "vacuum"
    bare = Experiment.from_dict(copy.deepcopy(raw))
    assert deorientation_parameters(bare) is None and "No deorientation set" in deorientation_warning(bare)
    del raw["lifetime"]
    assert deorientation_parameters(Experiment.from_dict(raw)) == Parameters()


# -- density matrices -----------------------------------------------------------------------------------------------


def states(spin: float, n: int = 40) -> np.ndarray:
    d = int(2 * spin + 1)
    rng = np.random.default_rng(3)
    x = rng.normal(size=(n, d, 2)) + 1j * rng.normal(size=(n, d, 2))
    rho = x @ np.swapaxes(x.conj(), 1, 2)
    return rho / np.einsum("naa->n", rho)[:, None, None]


@pytest.mark.parametrize("spin", [0.0, 1.5, 2.0, 4.0])
def test_an_attenuated_state_keeps_its_trace_and_is_hermitian(spin):
    rho = states(spin)
    m = TwoState("194Pt", spin, ranks="all")
    g = m.g(np.linspace(0.0, 80.0, len(rho)), 0.05)
    out = attenuated(rho, spin, g)
    assert np.abs(np.einsum("naa->n", out) - 1).max() < 1e-12
    assert np.abs(out - np.swapaxes(out.conj(), 1, 2)).max() < 1e-12
    assert np.abs(out[0] - rho[0]).max() < 1e-12                              # no time: nothing happened
    for i in (7, 39):                                                          # each tensor times its G_k
        before, after = tensors_of(rho[i], spin), tensors_of(out[i], spin)
        assert all(after[k, q] == pytest.approx(before[k, q] * (g[k][i] if k else 1.0), abs=1e-12) for k, q in before)
    for still in (m.g(50.0, np.zeros(len(rho))), {}):                          # no speed; nothing asked for
        assert np.abs(attenuated(rho, spin, still) - rho).max() < 1e-12
    if spin > 1:
        from physim.nuclear.gamma_events import _attenuate

        two = {2: g[2], **({4: g[4]} if 4 in g else {})}                       # as the plunger's, rank by rank
        assert np.allclose(attenuated(rho, spin, two), _attenuate(rho, two), atol=1e-13)
        assert np.allclose(sum(tensor_basis(spin, k).reshape(2 * k + 1, -1).T @ tensor_basis(spin, k).reshape(2 * k + 1, -1)
                               for k in range(int(2 * spin) + 1)), np.eye(int(2 * spin + 1) ** 2), atol=1e-12)


# -- the plan -------------------------------------------------------------------------------------------------------


def test_what_never_flies_is_not_deoriented():
    """Without time or speed in vacuum, or a nucleus that stays in the target, the tensors are those of the feeding alone,
    which the solution gives: main's table, to the last digits."""
    exp = Experiment.example("coulex_ni62_pt194")
    for ex, level, final in ((cascade_excitation(exp), 1, 0), (cascade_excitation(exp), 4, 1), (other_excitation(exp), 1, 0)):
        orient = DepthOrientation(None, [ex])
        theta = np.linspace(5.0, 175.0, 30)
        want = orient.coefficients(level, final, theta, 237.0)
        lives = [0.0] + [3.0 * 2**n for n in range(len(ex.scheme.levels) - 1)]
        left = np.linspace(0.0, 3.0, 30)
        f = Faded(ex, level, lives, Parameters())
        for fade in (f.at(np.full(30, np.inf), np.full(30, 0.05)), f.at(left, np.zeros(30)),
                     Faded(ex, level, lives, Parameters(g_factor=0.0)).at(left, np.full(30, 0.05))):
            got = orient.coefficients(level, final, theta, 237.0, fade)
            assert set(got) == set(want)
            assert all(np.allclose(got[k], want[k], rtol=1e-11, atol=1e-12) for k in want)
        # With it what is left of the orientation is smaller, and the population stays.
        got = orient.coefficients(level, final, theta, 237.0, f.at(left, np.full(30, 0.05)))
        assert sum(np.abs(got[k]).sum() for k in want if k[0]) < 0.9 * sum(np.abs(want[k]).sum() for k in want if k[0])
        assert np.array_equal(got[0, 0], want[0, 0])
    assert Faded(ex, 1, [0.0] * len(ex.scheme.levels), Parameters()).off          # no level lives a while
    assert Faded(ex, 1, lives, Parameters(rate_per_ps=0.0)).off and not f.off


def test_one_level_is_deoriented_by_the_share_of_its_decays_outside():
    ex = other_excitation(Experiment.example("coulex_ni62_pt194"))
    tau, left, beta = 7.0, np.array([0.0, 0.5, 4.0, 200.0]), np.array([0.02, 0.03, 0.04, 0.05])
    lives = np.zeros(len(ex.scheme.levels))
    top = max(range(len(lives)), key=lambda n: ex.scheme.levels[n].energy.value)   # nothing feeds it
    lives[top] = tau
    fade = Faded(ex, top, lives, Parameters()).at(left, beta)
    g = TwoState("194Pt", ex.scheme.levels[top].spin).average(tau, beta)
    share = np.exp(-left / tau)
    for k in fade:
        assert fade[k][top] == pytest.approx(1 - share + share * g[k], rel=1e-10)


def test_the_feeding_is_what_drawn_decays_average_to():
    """The factors against decays drawn one by one: 4⁺ → 2⁺ → 0⁺ of ¹⁹⁴Pt, the level fed from above deoriented by its own
    time in vacuum after the other's, some nuclei still in the target when one or both decay; also with one that decays
    at once."""
    from physim.nuclear.angular import u_coefficient

    ex = other_excitation(Experiment.example("coulex_ni62_pt194"))
    levels = ex.scheme.levels
    four, two = next(n for n, lv in enumerate(levels) if lv.spin == 4), 1
    assert levels[two].spin == 2 and ex.decay(four, two).branch == pytest.approx(1.0)
    rng, n, beta, left = np.random.default_rng(11), 400_000, 0.05, 1.3
    for t4, t2 in ((2.0, 5.0), (3.0, 0.0), (0.0, 4.0), (3.0, 3.0)):
        lives = np.zeros(len(levels))
        lives[four], lives[two] = t4, t2
        fade = Faded(ex, two, lives, Parameters()).at(np.array([left]), np.array([beta]))
        a = rng.exponential(1.0, n) * t4
        b = a + rng.exponential(1.0, n) * t2
        g4 = TwoState("194Pt", 4).g(np.maximum(a - left, 0.0), beta)
        g2 = TwoState("194Pt", 2).g(np.where(b > left, b - np.maximum(a, left), 0.0), beta)
        for k in (2, 4):
            one = g4[k] * g2[k] * u_coefficient(k, 4, 2, 2)
            assert fade[k][four][0] == pytest.approx(one.mean(), abs=4 * one.std() / math.sqrt(n) + 1e-9), (t4, t2, k)


def test_the_table_takes_it_with_the_setups_half_lives_and_g_factor():
    """⁶²Ni (2 ps, v/c 0.07) keeps most of its pattern and ¹⁹⁴Pt (60 ps) hardly any; what is kept in memory for one setting
    is not given for another."""
    on, off = setup("coulex_ni62_pt194"), setup("coulex_ni62_pt194", "none")
    a, b = factors(on), factors(off)
    assert 0.90 < a["CD", "Ge45", ""] / b["CD", "Ge45", ""] < 0.96 and 0.97 < a["CD", "Ge0", ""] / b["CD", "Ge0", ""] < 1.05
    pt, flat = factors(on, nucleus="other"), factors(off, nucleus="other")
    assert 0.75 < pt["CD", "Ge45", ""] / flat["CD", "Ge45", ""] < 0.85
    # Towards isotropy: in every pair, between 1 (with the forward throw of some percent) and what it was.
    assert all(abs(pt[k] - 1) < abs(flat[k] - 1) for k in pt)
    # A larger moment fades more; that of the reaction's nucleus is not the other's.
    strong = setup("coulex_ni62_pt194", g=1.2)
    assert all(abs(x - 1) < abs(a[k] - 1) + 0.02 for k, x in factors(strong).items())
    assert factors(strong)["CD", "Ge45", ""] < a["CD", "Ge45", ""] - 0.02
    assert factors(strong, nucleus="other") == pt and factors(setup("coulex_ni62_pt194", g=0.0)) == b
    assert factors(on) == a and factors(off) == b
    # The rates and the coincidences expected follow it.
    from physim.nuclear.rates import Rates

    assert Rates(on).gamma_efficiency("CD")[0] < 0.99 * Rates(off).gamma_efficiency("CD")[0]


def test_one_state_is_deoriented_only_with_a_half_life():
    none = factors(ni58(how="none"))
    assert factors(ni58(None)) == factors(ni58(None, "none"))
    assert factors(setup("coulex_ni58")) == factors(setup("coulex_ni58", "none"))     # the template has no scheme
    long, short = factors(ni58("20e-12 s")), factors(ni58("0.65e-12 s"))
    key = ("CD", "Ge90", "") if ("CD", "Ge90", "") in none else next(iter(none))
    assert abs(long[key] - none[key]) > 4 * abs(short[key] - none[key]) > 0


# -- the events -----------------------------------------------------------------------------------------------------


def test_the_events_of_one_state_follow_the_plan():
    """Counts of a setup's events over those of isotropic ones, against the plan's, crystal by crystal, aimed or drawn over
    the whole sphere; the plan without deorientation is far from them."""
    exp = ni58()
    raw = exp.to_dict()
    raw["reaction"]["emission"] = "isotropic"
    flat = Experiment.from_dict(raw)
    got, even = (simulate_gammas(e, 60_000, seed=s, gammas_per_event=40) for e, s in ((exp, 21), (flat, 22)))
    whole = simulate_gammas(exp, 60_000, seed=23, aimed=False, gammas_per_event=10)
    plan, plain = factors(exp), factors(flat)
    off = factors(ni58(how="none"))

    def counted(g, crystal, detector):
        m = g.select(crystal, detector)
        per = np.bincount(g["particle"][m], weights=g["weight"][m])
        return per.sum(), math.sqrt((per**2).sum())

    for g in (got, whole):
        pulls = []
        for key, want in plan.items():
            (a, sa), (b, sb) = counted(g, key[1], key[0]), counted(even, key[1], key[0])
            err = (a / b) * math.hypot(sa / a, sb / b)
            pulls.append([(a / b - x / plain[key]) / err for x in (want, off[key])])
        on, without = np.array(pulls).T
        assert np.abs(on).max() < 4 and (without**2).sum() > (on**2).sum() + 40
    assert any("Deorientation in vacuum" in n for n in got.notes) and got["vacuum_ps"].max() > 20
    assert not any("Deorientation" in n for n in even.notes)


def test_a_cascade_hands_on_what_is_left_and_each_gamma_ray_has_its_own_time():
    """4⁺ → 2⁺ → 0⁺ of ¹⁹⁴Pt from one collision angle, in vacuum from the start: the second γ ray's mean pattern is that
    of the 2⁺ level fed by an attenuated 4⁺ and attenuated again. (Feeding with the unattenuated state, as without
    deorientation, is far outside.)"""
    ex = other_excitation(Experiment.example("coulex_ni62_pt194"))
    levels = ex.scheme.levels
    four = next(n for n, lv in enumerate(levels) if lv.spin == 4)
    n, t4, t2, beta = 60_000, 6.0, 1.0, 0.05
    theta = np.full(n, 130.0)
    rng = np.random.default_rng(5)
    src = np.repeat(np.arange(n), 2)
    initial, final = np.tile([four, 1], n), np.tile([1, 0], n)
    taken = (src, initial, final, np.ones(2 * n, dtype=bool))
    time = rng.exponential(1.0, 2 * n) * np.tile([t4, t2], n)
    lives = np.zeros(len(levels))
    lives[four], lives[1] = t4, t2
    vacuum = _InVacuum(ex.scheme, "194Pt", Parameters(), initial, time, np.full(2 * n, beta), lives, np.cumsum(time), None)
    known, deferred = ge._correlated(ex, theta, taken, np.full(n, four), rng, False, vacuum=vacuum)
    assert not deferred
    second = known[1::2]
    z = second[:, 2] ** 2
    p2 = (1.5 * z - 0.5).mean()
    err = (1.5 * z - 0.5).std() / math.sqrt(n)
    # ⟨P₂(cos θ)⟩ = a_20 / 5, from the collision's 4⁺ alone.
    tensors = ex.tensor_table(four)
    grid = ex.table()[0]
    pop = np.array([s.direct[four] for s in ex.table()[1]])
    t20 = np.interp(130.0, grid, tensors[2, 0].real) / np.interp(130.0, grid, pop)
    a2 = ex.distribution_coefficient(1, 0, 2)
    fade = Faded(ex, 1, lives, Parameters()).at(np.zeros(1), np.full(1, beta))[2][four][0]
    plain = Faded(ex, 1, lives, Parameters()).at(np.full(1, np.inf), np.full(1, beta))[2][four][0]
    assert p2 == pytest.approx(t20 * fade * a2 / 5, abs=4 * err)
    own = TwoState("194Pt", 2).average(t2, beta)[2] * plain                    # only the emitting level's own
    assert abs(t20 * own * a2 / 5 - p2) > 8 * err
    assert vacuum.looked > n and vacuum.below <= 1e-3 * vacuum.looked


def test_a_converted_step_takes_the_average_over_its_levels_life():
    ex = other_excitation(Experiment.example("coulex_ni62_pt194"))
    four = next(n for n, lv in enumerate(ex.scheme.levels) if lv.spin == 4)
    lives = np.zeros(len(ex.scheme.levels))
    lives[four] = 6.0
    way = (np.array([0.0, 0.4, np.inf, 0.4]), np.full(4, 0.05))
    vacuum = _InVacuum(ex.scheme, "194Pt", Parameters(), np.zeros(1, int), np.zeros(1), np.zeros(1), lives,
                       np.array([1.5]), lambda i: (way[0][i], way[1][i]))
    rho = states(4.0, 4)
    out = vacuum.converted(rho, four, np.arange(4), np.array([-1, -1, -1, 0]))
    g = TwoState("194Pt", 4, ranks="all").average(6.0, 0.05)
    p = np.array([1.0, math.exp(-0.4 / 6.0), 0.0, 1.0])     # never out; fed after it had left
    assert np.allclose(out, attenuated(rho, 4.0, {k: 1 - p + p * v for k, v in g.items()}), atol=1e-13)
    assert vacuum.converted(rho, 1, np.arange(4), np.full(4, -1)) is rho      # a level without a half-life


def test_a_nucleus_that_stays_in_matter_is_not_deoriented():
    """A thick backing stops every recoil: the events are those without, to the last digit."""
    raw = ni58(how="vacuum").to_dict()
    raw["target"]["backing"] = {"material": "Au", "thickness": "20 mg/cm2"}
    on = Experiment.from_dict(copy.deepcopy(raw))
    raw["reaction"]["deorientation"] = "none"
    off = Experiment.from_dict(raw)
    a, b = simulate_gammas(on, 20_000, seed=4), simulate_gammas(off, 20_000, seed=4)
    assert len(a) and not a["vacuum_ps"].any() and a["stopped"].any()
    assert all(np.array_equal(a.columns[k], b.columns[k]) for k in a.columns)


# -- stored runs ----------------------------------------------------------------------------------------------------


def _as_taken_before(run: RunData) -> RunData:
    """A run as physim stored it before it had deorientation: a setup that says nothing about it."""
    path = run.folder / "setup.toml"
    text = path.read_text(encoding="utf-8")
    assert len(re.findall(r'(?m)^deorientation = "none"$', text)) == 1
    path.write_text(re.sub(r"(?m)^deorientation = .*\n", "", text), encoding="utf-8")
    return RunData.load(run.folder)


def test_a_run_says_which_deorientation_its_events_have(tmp_path):
    exp = Experiment.example("coulex_ni58")
    folder = ExperimentFolder.create(exp, root=tmp_path)
    run = folder.take_run("beam", duration="2 s")
    assert exp.excitation.deorientation is None
    assert run.setup["reaction"]["deorientation"] == "vacuum"
    assert 'deorientation = "vacuum"' in (run.folder / "setup.toml").read_text(encoding="utf-8")
    assert run.changes(exp.to_dict()) == [] and physics_changes(exp.to_dict(), run.setup) == []
    assert run.changes(setup("coulex_ni58", "none").to_dict()) == ["reaction.deorientation: vacuum → none"]
    assert run.changes(setup("coulex_ni58", g=0.3).to_dict()) == ["reaction.deorientation_g: (none) → 0.3"]
    assert physics_changes(setup("coulex_ni58", "none").to_dict(), setup("coulex_ni58", "none", 0.3).to_dict()) == []
    source = folder.take_run("source", duration="10 s")
    assert "deorientation" not in source.setup["reaction"] and source.changes(exp.to_dict()) == []
    # A plunger has none of it in either case: nothing is written, and an old run's setup is the current one's physics.
    ru98 = Experiment.example("coulex_ru98_plunger")
    plunger = ExperimentFolder.create(ru98, root=tmp_path).take_run("beam", duration="0.1 s", distance_um=25.0)
    assert "deorientation" not in plunger.setup["reaction"] and plunger.changes(ru98.to_dict()) == []
    assert deorientation_of(plunger.experiment) == "none"


def test_a_run_stored_before_is_read_without_it_and_goes_on_without_it(tmp_path, monkeypatch):
    p = Planner(ni58())
    p.create_experiment(root=tmp_path)
    p.set("reaction", "deorientation", "none")
    old = _as_taken_before(p.start_run("beam", duration="2 s"))
    p.set("reaction", "deorientation", None)
    assert deorientation_of(p.experiment) == "vacuum"
    assert old.setup["reaction"]["deorientation"] == "none" == deorientation_of(old.experiment)
    assert taken_setup(old.folder).excitation.deorientation == "none"
    want = ["reaction.deorientation: none → vacuum"]
    assert old.changes(p.experiment.to_dict()) == want
    p.load_run(1)
    assert p.run_status()["stale"] and p.run_status()["changes"] == want and p.runs()[0]["stale"]
    from physim.nuclear.logbook import setup_history

    assert [h["changes"] for h in setup_history(p)] == [[], want]
    # Its analysis is that of its own setup: its correlation and the events of its simulated peak.
    a = p.analysis()
    own, now = ({k[:2]: v for k, v in factors(e).items()} for e in (old.experiment, p.experiment))
    assert a._correlation == pytest.approx(own, rel=1e-12) and a._correlation != pytest.approx(now, rel=1e-3)
    from physim.nuclear import analysis
    from physim.nuclear.analysis import Settings

    monkeypatch.setattr(analysis, "SHAPE_EVENTS", 20_000)
    monkeypatch.setattr(analysis, "SHAPE_EMISSIONS", 40_000)
    shape = a._simulated_gammas(Settings(), 0)
    assert shape is not None and not shape["vacuum_ps"].any()
    # Continued, it goes on without, and says so from then on; so does the same run taken again; a new one has it.
    longer = p.extend_run("1 s")
    assert longer.setup["reaction"]["deorientation"] == "none" and longer.summary["extended"] == 1
    assert 'deorientation = "none"' in (longer.folder / "setup.toml").read_text(encoding="utf-8")
    from physim.nuclear.runs import retake

    taker = retake(p.folder, p.folder.load_run(1))
    while not taker.done:
        taker.step()
    assert taker.write().setup["reaction"]["deorientation"] == "none"
    assert p.start_run("beam", duration="1 s").setup["reaction"]["deorientation"] == "vacuum"


# -- GOSIA ----------------------------------------------------------------------------------------------------------


def test_the_gosia_input_has_the_setups_deorientation():
    def vac(exp, role):
        lines = Multistep(exp, exp.levels[role], role).gosia_input("4pi").splitlines()
        at = lines.index("VAC,")
        end = lines.index("0 0", at)
        assert lines[end + 1:end + 3] == ["END,", ""]
        assert end - at <= 7, "GOSIA reads seven pairs at most, the closing 0 0 among them"
        return {int(x.split()[0]): float(x.split()[1]) for x in lines[at + 1:end]}

    exp = setup("coulex_ni62_pt194", g=-0.17)
    assert vac(exp, "beam") == pytest.approx({2: 0.02, 3: 0.0345, 4: 3.5, 5: 0.17, 6: 6e-6, 7: 0.6})
    assert vac(exp, "target")[5] == pytest.approx(78 / 194)
    for none in (setup("coulex_ni62_pt194", "none"), setup("coulex_ni62_pt194", g=0.0)):
        assert vac(none, "beam") == {3: 0.0}


# -- the analysis ---------------------------------------------------------------------------------------------------


def test_the_analysis_gives_b_e2_back_from_deoriented_events():
    """⁶²Ni on ¹⁹⁴Pt: B(E2) back within its errors; analysed without deorientation the same events read about 2.7 % less
    (12 samples of 400 000 reactions: −2.65 ± 0.02 %)."""
    from physim.nuclear.analysis import Analysis, Settings

    on, off = Experiment.example("coulex_ni62_pt194"), setup("coulex_ni62_pt194", "none")
    g = simulate_gammas(on, 150_000, seed=3)
    s = Settings(particle_gate="all")
    a, b = Analysis(on, g).run(s), Analysis(off, g).run(s)
    truth = on.excitation.b_up_e2fm
    assert abs(a.b_e2fm4 / truth - 1) < 4 * np.hypot(a.monte_carlo, a.statistical)
    assert b.b_e2fm4 / a.b_e2fm4 == pytest.approx(0.9735, abs=0.004)
