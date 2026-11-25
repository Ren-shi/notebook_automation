"""The E1 polarisation correction (GOSIA's form, manual eq. 3.27–3.28, after Alder and Winther, appendix J): every E2
coupling is multiplied by 1 − z a/r along the orbit, z = d E_cm A / Z² of the excited nucleus. It vanishes with d = 0,
its size matches a case worked by hand and GOSIA's; the setup's value (GOSIA's d = 0.0038 when it sets none) reaches the
events, the plan and the analysis; the analysis gives back the B(E2) a run was made with; and a run stored before there
was a correction is analysed without it."""

import copy
import json
import re

import numpy as np
import pytest

from physim.nuclear import Experiment, data
from physim.nuclear.analysis import Analysis, Settings
from physim.nuclear.coulex import (GOSIA_DIPOLE_POLARISATION, Coulex, orbit_integrals, polarisation_strength,
                                   polarised_integrals)
from physim.nuclear.coupled import CoupledChannels
from physim.nuclear.gamma import (cascade_excitation, e1_polarisation, excitation_of, other_excitation,
                                  state_transition)
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.multistep import Multistep
from physim.nuclear.orientation import Excitation, simple_scheme
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates, channels, coulex_for, depth_probability, stack
from physim.nuclear.runs import ExperimentFolder, RunData, physics_changes, taken_setup

D = GOSIA_DIPOLE_POLARISATION
NI58 = ("16O", "58Ni", 30.0)
NI62 = ("62Ni", "194Pt", 237.0)


def setup(name: str, d=None, **reaction) -> Experiment:
    raw = Experiment.example(name).to_dict()
    raw["reaction"].update(reaction)
    if d is not None:
        raw["reaction"]["e1_polarisation"] = d
    return Experiment.from_dict(raw)


# -- the solvers ----------------------------------------------------------------------------------------------------


def test_z_is_gosias_for_target_and_projectile_excitation():
    ni62, pt194 = data.nuclide("62Ni"), data.nuclide("194Pt")
    m1, m2 = ni62.nuclear_mass_mev, pt194.nuclear_mass_mev
    # GOSIA's subroutine LOAD, its two branches as it writes them (EP the beam energy, a2 and z2 the excited nucleus).
    target = D * 237.0 * 194 / (78**2 * (1 + m1 / m2))
    projectile = D * 237.0 * 62 / (28**2 * (1 + m1 / m2))
    c = {e: Coulex(*NI62, excite=e, energy=1.0, b_up=1.0, e1_polarisation=D) for e in ("target", "projectile")}
    # The collision's own E_cm is relativistic: 0.1 % below E A₂/(A₁ + A₂) here.
    assert c["target"].polarisation == pytest.approx(target, rel=2e-3)
    assert c["projectile"].polarisation == pytest.approx(projectile, rel=2e-3)
    assert c["projectile"].polarisation == pytest.approx(D * c["projectile"].e_cm * 62 / 28**2, rel=1e-12)
    assert c["projectile"].polarisation == pytest.approx(0.0539, abs=1e-4)
    assert polarisation_strength(0.0, 180.0, ni62) == polarisation_strength(None, 180.0, ni62) == 0.0
    # The same in every class, one z whatever the level.
    scheme = Experiment.example("coulex_ni62_pt194").levels["beam"]
    for cls in (Excitation, CoupledChannels):
        ex = cls(*NI62, scheme, excite="projectile", e1_polarisation=D)
        assert {p[2].polarisation for paths in ex.paths.values() for p in paths} == {c["projectile"].polarisation}


def test_without_a_d_nothing_changes_and_e1_and_e3_have_none():
    scheme = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2")
    for d in (0.0, None):
        assert Coulex(*NI58, energy=1.454, b_up=695.0, e1_polarisation=d).probability(150.0, exact=True) == \
            Coulex(*NI58, energy=1.454, b_up=695.0).probability(150.0, exact=True)
        assert Excitation(*NI58, scheme, e1_polarisation=d).probability(1, 150.0) == \
            Excitation(*NI58, scheme).probability(1, 150.0)
        a, b = CoupledChannels(*NI58, scheme, e1_polarisation=d), CoupledChannels(*NI58, scheme)
        assert a.polarisation == 0.0 and np.array_equal(a.solve(150.0), b.solve(150.0))
        assert all(np.array_equal(x, y) for x, y in zip(a.solve_many([60.0, 150.0]), b.solve_many([60.0, 150.0])))
    # GOSIA's QE multiplies the collision functions of λ = 2 alone.
    for mult, b in (("E1", 1e-3), ("E3", 1e4)):
        with_d = Coulex(*NI58, energy=3.0, multipolarity=mult, b_up=b, e1_polarisation=0.05)
        assert with_d.polarisation > 0
        assert with_d.probability(150.0, exact=True) == Coulex(*NI58, energy=3.0, multipolarity=mult,
                                                               b_up=b).probability(150.0, exact=True)
        s = simple_scheme("58Ni", "3.0 MeV", mult, b)
        assert np.array_equal(CoupledChannels(*NI58, s, e1_polarisation=0.05).solve(150.0),
                              CoupledChannels(*NI58, s).solve(150.0))


def test_the_polarised_orbit_integral_is_the_integral_with_the_factor():
    z, eps, xi = 0.3, 1.2, 0.7
    w = np.linspace(-12.0, 12.0, 400_001)
    r = eps * np.cosh(w) + 1
    e_iphi = ((np.cosh(w) + eps) + 1j * np.sqrt(eps**2 - 1) * np.sinh(w)) / r
    got = polarised_integrals(2, eps, xi, z)
    for mu in range(-2, 3):
        f = np.exp(1j * xi * (eps * np.sinh(w) + w)) * e_iphi**mu * (1 - z / r) / r**2
        assert got[mu + 2] == pytest.approx(np.sum(0.5 * (f[1:] + f[:-1])) * (w[1] - w[0]), abs=1e-7)
    assert polarised_integrals(2, eps, xi, 0.0) is orbit_integrals(2, eps, xi)
    assert polarised_integrals(3, eps, xi, z) is orbit_integrals(3, eps, xi)


def test_its_size_in_sudden_backscattering_is_the_hand_worked_one():
    """At 180° (ε = 1) only μ = 0 couples, and with ξ → 0 the first-order amplitude is ∫ dw / (cosh w + 1)² = 2/3; the
    correction takes z ∫ dw / (cosh w + 1)³ = 4z/15 from it, so the probability falls by (1 − 0.4 z)². A level at 1 keV
    has ξ = 5 × 10⁻⁴."""
    scheme = simple_scheme("58Ni", "0.001 MeV", "E2", "1 e2fm4")
    for d in (D, 0.05):
        c0, c1 = (Coulex(*NI58, energy=0.001, b_up=1.0, e1_polarisation=x) for x in (0.0, d))
        expected = (1 - 0.4 * c1.polarisation) ** 2
        assert c0.xi < 1e-3 and expected < 0.995
        assert c1.probability(180.0, exact=True) / c0.probability(180.0, exact=True) == pytest.approx(expected, abs=1e-6)
        e0, e1 = (Excitation(*NI58, scheme, e1_polarisation=x) for x in (0.0, d))
        assert e1.probability(1, 180.0) / e0.probability(1, 180.0) == pytest.approx(expected, abs=1e-6)
        k0, k1 = (CoupledChannels(*NI58, scheme, tolerance=1e-6, e1_polarisation=x) for x in (0.0, d))
        assert k1.probability(1, 180.0) / k0.probability(1, 180.0) == pytest.approx(expected, abs=1e-5)
        # Away from 180° the coupled equations still are first order for this weak state, the correction included.
        assert k1.probability(1, 90.0) / k0.probability(1, 90.0) == pytest.approx(
            c1.probability(90.0, exact=True) / c0.probability(90.0, exact=True), abs=1e-5)


def test_the_diagonal_coupling_has_it_as_second_order_theory_does():
    """The integration applies the factor as the nested quadrature of second-order perturbation theory does (another
    integration of the same couplings), with a diagonal element."""
    s = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0064 e2b2")
    s.set_matrix_element(1, 1, "E2", 26.4, source="user")
    p = {}
    for d in (0.0, 0.05):
        cc = CoupledChannels(*NI58, s, tolerance=1e-6, e1_polarisation=d)
        a = cc.solve(170.0)[cc._level_slices[1], 0]
        first, second = cc.second_order_check(1, 170.0)
        p[d] = (np.sum(np.abs(a) ** 2), np.sum(np.abs(first + second) ** 2))
    change, theory = p[0.05][0] / p[0.0][0], p[0.05][1] / p[0.0][1]
    assert change < 0.95
    assert change == pytest.approx(theory, abs=2e-3)


def test_the_probabilities_still_sum_to_one_and_the_angles_together_are_the_single_ones():
    scheme = Experiment.example("coulex_ni62_pt194").levels["beam"]
    cc = CoupledChannels(*NI62, scheme, excite="projectile", e1_polarisation=0.05)
    angles = [20.0, 90.0, 180.0]
    for theta, a in zip(angles, cc.solve_many(angles)):
        assert np.abs(a - cc.solve(theta)).max() < 1e-12
        assert np.sum(np.abs(a) ** 2) / a.shape[1] == pytest.approx(1.0, abs=1e-8)
        assert sum(cc.probabilities(theta).values()) == pytest.approx(1.0, abs=1e-8)
    # The forward angles that are left in first order have it too.
    weak = CoupledChannels(*NI62, scheme, excite="projectile", tolerance=1e-4, e1_polarisation=D)
    grid, states = weak.table(10.0)
    first = Excitation(*NI62, scheme, excite="projectile", e1_polarisation=D)
    assert states[0].direct[1] == pytest.approx(first.probability(1, 10.0), rel=1e-12)


def test_the_reorientation_is_weakened_less_than_the_excitation():
    """⟨2‖E2‖2⟩ = +6.6 e fm² of the template's ⁶²Ni 2⁺ state raises its population by 8 % at 170°. The correction
    weakens the diagonal coupling with the transitional one, so that effect shrinks by a percent or two of itself while
    the population falls by several; and it falls off one power of r faster than the excitation, so it weighs the
    backward angles more, as a smaller diagonal element does."""
    scheme = Experiment.example("coulex_ni62_pt194").levels["beam"]
    flat = copy.deepcopy(scheme)
    flat.set_matrix_element(1, 1, "E2", 0.0, source="user")
    p = {(d, s is flat): CoupledChannels(*NI62, s, excite="projectile", e1_polarisation=d).probabilities(170.0)[1]
         for d in (0.0, D) for s in (scheme, flat)}
    effect = {d: p[(d, False)] / p[(d, True)] - 1 for d in (0.0, D)}
    assert 0.07 < effect[0.0] < 0.09 and 0.97 < effect[D] / effect[0.0] < 1.0
    assert 0.93 < p[(D, False)] / p[(0.0, False)] < 0.97
    cc = {d: CoupledChannels(*NI62, scheme, excite="projectile", e1_polarisation=d) for d in (0.0, D)}
    ratio = {d: c.probabilities(160.0)[1] / c.probabilities(60.0)[1] for d, c in cc.items()}
    assert -0.02 < ratio[D] / ratio[0.0] - 1 < -0.01


# -- the setup ------------------------------------------------------------------------------------------------------


def test_the_setups_value_reaches_the_solutions_the_rates_and_the_gosia_input():
    p = Planner.example("coulex_ni62_pt194")
    assert p.experiment.excitation.e1_polarisation is None
    assert e1_polarisation(p.experiment) == GOSIA_DIPOLE_POLARISATION == 0.0038
    assert "DIP,3.8" in Multistep(p.experiment, p.experiment.levels["beam"], "beam").gosia_input("4pi").splitlines()
    for wrong in (5, -0.001, 0.06):
        assert not p.set("reaction", "e1_polarisation", wrong) and "DIP," in p.problems[0]
    assert p.set("reaction", "e1_polarisation", 0.005)
    on, off = p.experiment, setup("coulex_ni62_pt194", 0.0)
    assert Experiment.from_toml(on.to_toml()).excitation.e1_polarisation == 0.005
    assert Experiment.from_toml(off.to_toml()).excitation.e1_polarisation == 0.0 and e1_polarisation(off) == 0.0
    for exp, d in ((on, 0.005), (off, 0.0), (Experiment.example("coulex_ni62_pt194"), 0.0038)):
        # Each nucleus has its own z; the caches keep the solutions of each d apart.
        for ex in (cascade_excitation(exp), other_excitation(exp), cascade_excitation(exp, order="first")):
            assert ex.e1_polarisation == d
        z = cascade_excitation(exp)
        assert z.polarisation == pytest.approx(d * z.cx.e_cm * 62 / 28**2)
        pt = other_excitation(exp)
        assert pt.polarisation == pytest.approx(d * pt.cx.e_cm * 194 / 78**2)
        assert Multistep(exp, exp.levels["beam"], "beam").solver(230.0).e1_polarisation == d
        assert Multistep(exp, exp.levels["beam"], "beam", first_order=True).solver(230.0).e1_polarisation == d
        ch = next(c for c in channels(stack(exp), exp) if c.excitation is not None)
        assert depth_probability(exp, ch).nodes[0].solution.e1_polarisation == d
    # One state, no scheme: the first-order probability of the rates, and the orientation of its γ rays.
    for d, want in ((0.0, 0.0), (0.005, 0.005), (None, 0.0038)):
        one = setup("coulex_ni58", d)
        ch = next(c for c in channels(stack(one), one) if c.excitation is not None)
        c = coulex_for(one, ch)
        assert c.e1_polarisation == excitation_of(one).e1_polarisation == want
        assert c.polarisation == pytest.approx(want * c.e_cm * 58 / 28**2)
        assert all(n.e1_polarisation == want for n in depth_probability(one, ch).nodes)
    # GOSIA's input: the setup's d; none is a millionth of it, as GOSIA takes 0 as its default; or what is asked for.
    assert "DIP,5" in Multistep(on, on.levels["beam"], "beam").gosia_input("4pi").splitlines()
    assert "DIP,0.001" in Multistep(off, off.levels["beam"], "beam").gosia_input("4pi").splitlines()
    for dipole, line in ((20.0, "DIP,20"), (0.0, "DIP,0.001")):
        assert line in Multistep(on, on.levels["beam"], "beam").gosia_input("4pi", dipole=dipole).splitlines()


def test_only_a_coulomb_excitation_has_the_setting():
    raw = Experiment.example("alpha_on_gold").to_dict()
    raw["reaction"] = {"type": "elastic", "e1_polarisation": 0.0038}
    with pytest.raises(ValueError, match="e1_polarisation only apply"):
        Experiment.from_dict(raw)
    # A setup that cannot set it has GOSIA's default where it solves a collision (its foils, a level scheme).
    assert e1_polarisation(Experiment.example("alpha_on_gold")) == 0.0038


def test_the_rates_of_the_templates_fall_by_the_correction():
    for name, detector, low, high in (("coulex_ni58", "CD", 0.990, 0.996), ("coulex_ni62_pt194", "CD", 0.94, 0.96)):
        on, off = Rates(Experiment.example(name)), Rates(setup(name, 0.0))
        assert low < on.rate(detector, what="excitations") / off.rate(detector, what="excitations") < high


def test_the_events_are_made_with_it_and_the_analysis_gives_the_b_e2_back():
    on, off = Experiment.example("coulex_ni62_pt194"), setup("coulex_ni62_pt194", 0.0)
    # The events against the plan: the same generator, the same probability.
    g = simulate_gammas(on, 150_000, seed=3)
    r = Rates(on)
    ev = g.events
    for det in ("CD", "DSSD-L"):
        mine = ev["counted"] & (ev["detector"] == ev.detectors.index(det))
        excited = np.isin(ev["channel"], [i for i, lab in enumerate(ev.channels) if "excited" in lab])
        w = ev["weight"][mine & excited]
        want = r.rate(det, what="excitations")
        assert abs(w.sum() - want) < 4 * np.sqrt(np.sum(w**2)) + 2e-3 * want
        assert want / Rates(off).rate(det, what="excitations") == pytest.approx(0.95, abs=0.005)
    with_d, without = Analysis(on, g), Analysis(off, g)
    # At the B(E2) the run was made with the all-orders fit solves the equations of the events.
    level, ground = state_transition(on)
    truth = on.excitation.b_up_e2fm
    grid, y = with_d._scheme_yield(truth, level, 120.0, 170.0)
    events = cascade_excitation(on)
    expected = [events._populated(float(t), events._densities(events.solve(float(t)))).gamma_yield(level, ground)
                for t in grid]
    assert np.allclose(y[(level, ground)], expected, rtol=1e-5)
    s = Settings(particle_gate="all")
    a, b = with_d.run(s), without.run(s)
    # B(E2) back within its errors; analysed without the correction the same run reads about 7 % less.
    assert abs(a.b_e2fm4 / truth - 1) < 4 * np.hypot(a.monte_carlo, a.statistical)
    assert b.b_e2fm4 / a.b_e2fm4 == pytest.approx(0.93, abs=0.01)
    # The budget's term of the scheme's other matrix elements is solved with it too.
    assert with_d._scheme_responses()["yields"][0][(level, ground)].mean() < \
        0.97 * without._scheme_responses()["yields"][0][(level, ground)].mean()


# -- stored runs ----------------------------------------------------------------------------------------------------


def _as_taken_before(run: RunData) -> RunData:
    """A run as physim stored it before the correction: events without it, and a setup that says nothing about it."""
    path = run.folder / "setup.toml"
    text = path.read_text(encoding="utf-8")
    assert len(re.findall(r"(?m)^e1_polarisation = 0(\.0)?$", text)) == 1
    path.write_text(re.sub(r"(?m)^e1_polarisation = .*\n", "", text), encoding="utf-8")
    assert "e1_polarisation" not in path.read_text(encoding="utf-8")
    return RunData.load(run.folder)


def test_a_run_says_which_polarisation_its_events_have(tmp_path):
    exp = Experiment.example("coulex_ni58")
    folder = ExperimentFolder.create(exp, root=tmp_path)
    run = folder.take_run("beam", duration="2 s")
    # The setup sets none; the run's says the default it was taken with, and is the current setup's physics.
    assert exp.excitation.e1_polarisation is None and "e1_polarisation" not in (folder.path / "experiment.toml").read_text(
        encoding="utf-8")
    assert run.setup["reaction"]["e1_polarisation"] == run.experiment.excitation.e1_polarisation == 0.0038
    assert "e1_polarisation = 0.0038" in (run.folder / "setup.toml").read_text(encoding="utf-8")
    assert run.changes(exp.to_dict()) == [] and physics_changes(exp.to_dict(), run.setup) == []
    assert run.changes(setup("coulex_ni58", 0.0).to_dict()) == ["reaction.e1_polarisation: 0.0038 → 0.0"]
    # One that sets it says its own; a source run has no such events and an elastic setup no such setting.
    other = folder.take_run("beam", duration="2 s", experiment=setup("coulex_ni58", 0.005))
    assert other.setup["reaction"]["e1_polarisation"] == 0.005
    assert other.changes(exp.to_dict()) == ["reaction.e1_polarisation: 0.005 → 0.0038"]
    source = folder.take_run("source", duration="10 s")
    assert "e1_polarisation" not in source.setup["reaction"] and source.changes(exp.to_dict()) == []
    gold = ExperimentFolder.create(Experiment.example("alpha_on_gold"), root=tmp_path).take_run("beam", duration="1 s")
    assert "e1_polarisation" not in gold.setup["reaction"]


def test_a_run_stored_before_the_correction_is_read_without_it(tmp_path):
    p = Planner.example("coulex_ni62_pt194")
    p.create_experiment(root=tmp_path)
    p.set("reaction", "e1_polarisation", 0.0)
    old = _as_taken_before(p.start_run("beam", duration="10 min"))
    p.set("reaction", "e1_polarisation", None)  # the current setup: nothing set, so the default
    assert old.setup["reaction"]["e1_polarisation"] == 0.0 == e1_polarisation(old.experiment)
    assert taken_setup(old.folder).excitation.e1_polarisation == 0.0
    want = ["reaction.e1_polarisation: 0.0 → 0.0038"]
    assert old.changes(p.experiment.to_dict()) == want
    p.load_run(1)
    assert p.run_status()["stale"] and p.run_status()["changes"] == want and p.runs()[0]["stale"]
    from physim.nuclear.logbook import setup_history

    assert [h["changes"] for h in setup_history(p)] == [[], want]
    # Its analysis is that of its own setup, the simulated peak's too.
    a = Analysis(p._data_experiment(), old.gammas(), run=old)
    s = Settings(particle_gate="all")
    assert cascade_excitation(a.experiment).e1_polarisation == 0.0
    truth = p.experiment.excitation.b_up_e2fm
    r = a.run(s)
    assert abs(r.b_e2fm4 / truth - 1) < 4 * np.hypot(r.monte_carlo, r.statistical)
    assert p.analyse({"particle_gate": "all"})["b_e2fm4"] == pytest.approx(r.b_e2fm4, rel=1e-9)
    # With today's default the same run would read some 7 % too much.
    wrong = Analysis(p.experiment, old.gammas(), run=old).run(s)
    assert wrong.b_e2fm4 / r.b_e2fm4 == pytest.approx(1.075, abs=0.015)
    # Extended, it goes on without the correction, and says so from then on; so does the same run taken again.
    stream = old.summary["generator"]["probabilities"]
    longer = p.extend_run("1 min")
    assert longer.setup["reaction"]["e1_polarisation"] == 0.0 and longer.summary["extended"] == 1
    assert longer.summary["generator"]["probabilities"] == pytest.approx(stream, rel=1e-12)
    assert "e1_polarisation = 0" in (longer.folder / "setup.toml").read_text(encoding="utf-8")
    assert p.start_run("beam", duration="1 s").setup["reaction"]["e1_polarisation"] == 0.0038


def test_a_kept_result_from_before_the_correction_is_not_shown_as_the_setups(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    p._keep("multistep.json", {"x": 1}, p.experiment_folder)
    assert p._kept("multistep.json", p.experiment_folder) == {"x": 1}
    path = p.experiment_folder / "multistep.json"
    kept = json.loads(path.read_text(encoding="utf-8"))
    assert kept.pop("e1_polarisation") == 0.0038
    path.write_text(json.dumps(kept), encoding="utf-8")
    assert p._kept("multistep.json", p.experiment_folder) is None
    p.set("reaction", "e1_polarisation", 0.0)
    path.write_text(json.dumps(dict(kept, setup=p.experiment.to_dict())), encoding="utf-8")
    assert p._kept("multistep.json", p.experiment_folder) == {"x": 1}
