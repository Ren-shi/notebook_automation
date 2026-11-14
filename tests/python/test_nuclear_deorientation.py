"""Deorientation in vacuum (the two-state model with GOSIA's parameters): the attenuation itself, the setup's
switch, and the plan and the events with it."""

from __future__ import annotations

import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.angular import wigner_6j
from physim.nuclear.deorientation import (DEFAULTS, TwoState, atomic_spin, attenuate, charge_states, hard_core,
                                          rank_projectors)
from physim.nuclear.events import simulate
from physim.nuclear.experiment import SetupError
from physim.nuclear.gamma import correlation_table, deorientation_factors, deorientation_of, other_excitation
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.levels import Value
from physim.nuclear.orientation import simple_scheme, tensors_of

#: The mean life of the first 2⁺ level of ¹⁹⁴Pt (T½ = 41.7 ps in the template) and of ⁶²Ni (1.442 ps), ps.
PT_TAU, NI_TAU = 41.7 / math.log(2), 1.442 / math.log(2)


def switched(exp: Experiment, value: str = "vacuum", **target) -> Experiment:
    d = exp.to_dict()
    d["reaction"]["deorientation"] = value
    d["target"].update(target)
    return Experiment.from_dict(d)


def gosia_integral(z: int, a: int, spin: float, k: int, tau: float, beta: float) -> float:
    """GOSIA's time-integral coefficient, written out once more from its subroutines GKK and XSTATIC, number by
    number, for one rank and one speed (default parameters)."""
    h = 1 / (1 + (z**0.45 * 0.012008 / beta) ** 1.666667)
    qcen = z * h**0.6
    dq = math.sqrt(qcen * (1 - h)) / 2
    up_q, down_q = min(int(qcen + 3 * dq + 0.5), z), max(int(qcen - 3 * dq - 0.5), 1)
    norm = sum(math.exp(-((qcen - q) / dq) ** 2 / 2) for q in range(down_q, up_q + 1))

    def alpha(j):
        total, f = 0.0, abs(spin - j)
        while f <= spin + j + 1e-9:
            total += ((2 * f + 1) * wigner_6j(f, f, k, spin, spin, j)) ** 2 / (2 * j + 1)
            f += 1
        return total

    a_static = sum(alpha(atomic_spin(z - q)) * math.exp(-((qcen - q) / dq) ** 2 / 2) / norm
                   for q in range(down_q, up_q + 1))
    a_fluct = alpha(3.0)
    wsp = 4789.0 * (z / a) * (6e-6 * z * beta**0.6) / 3.0 * 3.5
    w2 = wsp * wsp * 3.0 * 4.0 / 3.0 * k * (k + 1)
    xlam = (1 - a_fluct) * (1 - math.exp(-w2 / (1 - a_fluct))) / 3.5
    up = (0.02 * tau * a_static + 1) / (tau * 0.02 + 1) * 0.0345 * tau + 1
    down = tau * (xlam + 0.0345) + 1
    alp = (math.sqrt(9 * xlam**2 + 8 * xlam * 3.5 * (w2 - xlam**2)) - 3 * xlam) / 4 / xlam / 3.5
    ccf = 1 + xlam * tau * (down - 2 * alp**2 * tau * 3.5) / ((down + alp * tau) * (down + 2 * alp * tau))
    return up / down * ccf


# -- the model ------------------------------------------------------------------------------------------------------


def test_the_defaults_are_gosias():
    assert (DEFAULTS.atomic_spin, DEFAULTS.width_per_ps, DEFAULTS.rate_per_ps, DEFAULTS.correlation_ps,
            DEFAULTS.field_tesla, DEFAULTS.field_power, DEFAULTS.g_factor) == (3.0, 0.02, 0.0345, 3.5, 600.0, 0.6, None)
    assert TwoState("194Pt", 2.0).g_factor == pytest.approx(78 / 194)
    # GOSIA's atomic ground-state spins by number of electrons (its function ATS), with its quirks: 20 and 72
    # electrons give 1, a bare ion and more than 96 electrons 0.
    assert [atomic_spin(n) for n in (0, 1, 2, 7, 8, 20, 26, 46, 47, 65, 72, 96, 97)] == [
        0.0, 0.5, 0.0, 1.5, 2.0, 1.0, 4.0, 0.0, 0.5, 7.5, 1.0, 2.0, 0.0]


def test_a_static_field_leaves_the_hard_core():
    # No electronic spin, no attenuation; and the 6j sums that make G_k^st(0) = 1 add up to 1.
    assert hard_core(2.0, 0.0, 2) == pytest.approx(1.0) and hard_core(2.0, 3.0, 0) == pytest.approx(1.0)
    for spin, j, k in ((2.0, 0.5, 2), (2.0, 3.0, 4), (4.0, 2.5, 6), (1.5, 1.0, 2)):
        total = 0.0
        fs = np.arange(abs(spin - j), spin + j + 0.5)
        for f in fs:
            for fp in fs:
                total += (2 * f + 1) * (2 * fp + 1) / (2 * j + 1) * wigner_6j(f, fp, k, spin, spin, j) ** 2
        assert total == pytest.approx(1.0)
        assert 0 < hard_core(spin, j, k) < 1


def test_the_charge_states_are_nikolaev_dmitriev():
    q, w = charge_states(78, 0.04)
    h = 1 / (1 + (78**0.45 * 0.012008 / 0.04) ** (5 / 3))
    assert w.sum() == pytest.approx(1.0)
    assert np.sum(q * w) == pytest.approx(78 * h**0.6, abs=0.1)        # Q₀ = 31.6 for platinum at β = 0.04
    assert (q[0], q[-1]) == (23, 39)


@pytest.mark.parametrize("nuclide,spin,tau,beta", [("194Pt", 2.0, PT_TAU, 0.04), ("62Ni", 2.0, NI_TAU, 0.05),
                                                    ("62Ni", 4.0, 1.27, 0.05), ("194Pt", 4.0, 5.3, 0.03)])
def test_the_attenuation_starts_at_one_and_averages_to_gosias(nuclide, spin, tau, beta):
    model = TwoState(nuclide, spin)
    ranks = model.ranks
    # At t = 0 nothing has happened yet, and it starts flat.
    assert model.g(0.0, beta)[0] == pytest.approx(np.ones(len(ranks)))
    assert model.g(1e-4, beta)[0] == pytest.approx(np.ones(len(ranks)), abs=1e-6)
    # Its average over the decay is GOSIA's closed form (manual 4.30–4.33) ...
    t = np.linspace(0.0, 40 * tau, 200_001)
    numeric = np.trapezoid(model.g(t, np.full(len(t), beta)) * np.exp(-t / tau)[:, None], t, axis=0) / tau
    assert model.integral(tau, beta)[0] == pytest.approx(numeric, abs=2e-5)
    # ... which is GOSIA's own arithmetic, rank by rank.
    nuc = {"194Pt": (78, 194), "62Ni": (28, 62)}[nuclide]
    for k in (2, 4):
        assert model.integral(tau, beta)[0][k] == pytest.approx(gosia_integral(*nuc, spin, k, tau, beta), rel=1e-6)
    # At long times what is left is the hard core of the static field, reached through the fluctuating state.
    late = model.g(1e5, beta)[0]
    assert late == pytest.approx(model.hard_limit(beta)[0], abs=1e-9)
    assert model.integral(1e7, beta)[0] == pytest.approx(model.hard_limit(beta)[0], abs=1e-4)
    assert np.all(late[1:] < model.alpha_static(beta)[0][1:])
    # It only ever attenuates, and more for a higher rank.
    g = model.g(t[::1000], np.full(len(t[::1000]), beta))
    assert np.all((g >= 0) & (g <= 1 + 1e-12))
    assert np.all(np.diff(model.integral(tau, beta)[0][::2]) < 0)


def test_the_template_levels():
    """How much the first 2⁺ levels of the 62Ni template are deoriented, as GOSIA would have them: the platinum
    level lives 60 ps and keeps a fifth of its alignment; nickel's lives 2 ps and keeps nine tenths."""
    pt = TwoState("194Pt", 2.0).integral(PT_TAU, 0.040)[0]
    ni = TwoState("62Ni", 2.0).integral(NI_TAU, 0.050)[0]
    assert pt[2] == pytest.approx(0.199, abs=0.002) and pt[4] == pytest.approx(0.078, abs=0.002)
    assert ni[2] == pytest.approx(0.916, abs=0.002) and ni[4] == pytest.approx(0.783, abs=0.002)


def test_attenuating_a_density_matrix_scales_its_tensors():
    rng = np.random.default_rng(1)
    for spin in (1.0, 2.0, 2.5, 4.0):
        d = int(2 * spin + 1)
        assert rank_projectors(spin).sum(axis=0) == pytest.approx(np.eye(d * d))
        x = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        rho = x @ x.conj().T
        rho /= np.trace(rho).real
        g = np.r_[1.0, rng.uniform(0.0, 1.0, d - 1)]
        out = attenuate(rho[None], spin, g[None])[0]
        before, after = tensors_of(rho, spin, k_max=d), tensors_of(out, spin, k_max=d)
        for (k, q), v in before.items():
            assert after[(k, q)] == pytest.approx(g[k] * v, abs=1e-12)
        assert np.allclose(out, out.conj().T)


# -- the setup -----------------------------------------------------------------------------------------------------


def test_the_switch_is_on_unless_set_to_none():
    exp = Experiment.example("coulex_ni62_pt194")
    assert exp.excitation.deorientation is None and deorientation_of(exp) == "vacuum"
    assert "deorientation" not in exp.to_dict()["reaction"]
    on = switched(exp)
    assert deorientation_of(on) == "vacuum" and on.to_dict()["reaction"]["deorientation"] == "vacuum"
    off = switched(exp, "none")
    assert deorientation_of(off) == "none" and off.to_dict()["reaction"]["deorientation"] == "none"
    assert deorientation_of(Experiment.from_dict(off.to_dict())) == "none"
    # Unset and "vacuum" are the same plan; "none" is another.
    assert correlation_table(exp, nucleus="other") == correlation_table(on, nucleus="other")
    assert correlation_table(exp, nucleus="other") != correlation_table(off, nucleus="other")
    with pytest.raises(SetupError, match="deorientation"):
        switched(exp, "gas")


# -- the plan ------------------------------------------------------------------------------------------------------


def test_the_plan_evens_out_the_long_lived_platinum_level():
    exp = switched(Experiment.example("coulex_ni62_pt194"), "none")
    on = switched(exp)

    def factors(e, nucleus):
        return {(r["particle_detector"], r["gamma_detector"]): r["detector_factor"]
                for r in correlation_table(e, nucleus=nucleus)}

    for nucleus, shrink in (("reaction", (0.75, 0.95)), ("other", (0.1, 0.4))):
        off_, on_ = factors(exp, nucleus), factors(on, nucleus)
        # The γ rays come out more evenly: the factors spread less about 1. (Not each one moves towards 1: ranks
        # 2 and 4 are attenuated differently, and a crystal where they partly cancel can move away from it.)
        spread = np.std(list(on_.values())) / np.std(list(off_.values()))
        assert shrink[0] < spread < shrink[1], (nucleus, spread)
    # Platinum behind the strip detectors: Ge90 sees 0.44 of even emission without deorientation, 0.93 with it.
    pt_off, pt_on = factors(exp, "other"), factors(on, "other")
    assert pt_off[("DSSD-L", "Ge90")] == pytest.approx(0.44, abs=0.02)
    assert pt_on[("DSSD-L", "Ge90")] == pytest.approx(0.93, abs=0.02)


def ni58(half_life=None, deorientation="vacuum", backing=None) -> Experiment:
    """The ⁵⁸Ni example (its 2⁺ state alone, no cascade) with the state's half-life ``half_life`` (s)."""
    d = Experiment.example("coulex_ni58").to_dict()
    d["reaction"]["deorientation"] = deorientation
    if backing:
        d["target"]["backing"] = {"material": backing[0], "thickness": backing[1]}
    if half_life is not None:
        s = simple_scheme("58Ni", "1.454 MeV", "E2", d["reaction"]["b_up"])
        s.levels[1].half_life = Value(half_life, source="user")
        d["levels"] = {"target": s.to_dict()}
    return Experiment.from_dict(d)


def detector_ratios(on: Experiment, off: Experiment, events: int, seed: int, nucleus: int = 0,
                    transition=(1, 0), gammas=None) -> list:
    """For each particle detector × γ-ray detector: (γ rays of ``transition`` counted with deorientation over
    without, on the same particle events and seeds; its statistical error; the plan's factor with over without).
    ``gammas`` = (with, without) reuses γ-ray events already made."""
    if gammas is None:
        ev = simulate(off, events, seed)
        gammas = simulate_gammas(on, particle_events=ev, seed=seed), simulate_gammas(off, particle_events=ev, seed=seed)
    g_on, g_off = gammas
    label = "reaction" if nucleus == 0 else "other"
    plan_on = {(r["particle_detector"], r["gamma_detector"]): r["detector_factor"]
               for r in correlation_table(on, nucleus=label)}
    plan_off = {(r["particle_detector"], r["gamma_detector"]): r["detector_factor"]
                for r in correlation_table(off, nucleus=label)}
    out = []
    for key, p0 in plan_off.items():
        sums = []
        for g in (g_on, g_off):
            m = (g.select(key[1], key[0]) & (g["nucleus"] == nucleus) & (g["initial"] == transition[0])
                 & (g["final"] == transition[1]))
            sums.append((g["weight"][m].sum(), math.sqrt(np.sum(g["weight"][m] ** 2))))
        (a, ea), (b, eb) = sums
        # The two samples share their particles; their γ rays are drawn with the same numbers only as long as the
        # draws do not depend on the weights (a cascade's acceptance loops part them). The error of the ratio is
        # that of two independent samples, an upper bound.
        out.append((key, a / b, a / b * math.hypot(ea / a, eb / b), plan_on[key] / p0))
    return out


def test_without_half_lives_nothing_is_deoriented():
    """``coulex_ni58`` has no level scheme, so no half-life: every decay happens as the nucleus is formed, and the
    γ rays with the switch on are those without it, one for one."""
    off = ni58(deorientation="none")
    on = ni58()
    ev = simulate(off, 20_000, 5)
    a, b = simulate_gammas(on, particle_events=ev, seed=5), simulate_gammas(off, particle_events=ev, seed=5)
    for key in ("particle", "crystal", "measured", "weight"):
        assert np.array_equal(a[key], b[key]), key
    assert any("no level of 58Ni has a half-life" in n for n in a.notes)
    assert correlation_table(on) == correlation_table(off)


def test_one_state_the_plan_and_the_events_agree_detector_by_detector():
    """The ⁵⁸Ni 2⁺ state given a half-life of 30 ps (its own is under 1 ps): the recoils leave the target at
    β ≈ 0.02 and decay in vacuum; the plan's factors change as the events' counts do, in every pair of detectors."""
    rows = detector_ratios(ni58(30e-12), ni58(30e-12, "none"), 300_000, seed=7)
    changed = 0
    for key, ratio, err, want in rows:
        assert abs(ratio - want) < 4 * err + 0.01, (key, ratio, err, want)
        changed += abs(want - 1) > 0.05
    assert changed >= 4


def test_cascades_the_plan_and_the_events_agree_crystal_by_crystal():
    """The ⁶²Ni on ¹⁹⁴Pt template, both nuclei followed down their cascades: the γ rays of each nucleus's first
    2⁺ → 0⁺ transition in each γ-ray detector with each particle detector change with deorientation as the plan's
    factors do. Platinum's change by up to a factor of two (Ge90 behind the strip detectors)."""
    exp = switched(Experiment.example("coulex_ni62_pt194"), "none")
    on = switched(exp)
    ev = simulate(exp, 150_000, 9)
    gammas = simulate_gammas(on, particle_events=ev, seed=9), simulate_gammas(exp, particle_events=ev, seed=9)
    for nucleus, transition in ((1, (1, 0)), (0, (1, 0))):
        rows = detector_ratios(on, exp, 0, 9, nucleus=nucleus, transition=transition, gammas=gammas)
        pulls = [(ratio - want) / err for _, ratio, err, want in rows]
        assert max(abs(p) for p in pulls) < 4, (nucleus, rows)
        assert abs(np.mean(pulls)) < 1.2, (nucleus, pulls)
        if nucleus == 1:
            assert max(want for *_, want in rows) > 1.6


def test_a_deoriented_run_counts_what_the_plan_predicts(tmp_path):
    """A run of the clover template with deorientation: the platinum peak, whose predicted coincidences drop by a
    fifth with it, is counted as the plan predicts in every particle detector (some 14 700 counts in two minutes,
    against 18 000 without deorientation)."""
    from physim.nuclear.planner import Planner

    p = Planner.example("coulex_ni62_pt194_clovers")
    assert p.set("reaction", "deorientation", "vacuum"), p.problems
    assert p.gamma()["deorientation"] == "vacuum"
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="1 h", budget="2 min")
    rows = [r for r in p.compare()["rows"] if r["quantity"] in ("Coincidences in the peak", "Rate")
            or r["quantity"].startswith("Coincidences in the 194Pt")]
    assert len(rows) == 9 and all(r["agree"] for r in rows), rows
    pt = [r for r in rows if r["quantity"].startswith("Coincidences in the 194Pt")]
    measured, predicted = sum(r["measured"] for r in pt), sum(r["predicted"] for r in pt)
    assert measured * p.run.real_s > 12_000
    assert measured / predicted == pytest.approx(1.0, abs=0.02)


def test_a_nucleus_that_stops_in_the_backing_is_not_deoriented():
    """The platinum recoils go forward into a thick gold backing and stop there before they decay: their table is
    as without deorientation, while the nickel, scattered backwards out of the front of the target, is still
    deoriented."""
    exp = Experiment.example("coulex_ni62_pt194")
    thick = switched(exp, backing={"material": "Au", "thickness": "40 mg/cm2"})
    plain = switched(thick, "none")
    pt = other_excitation(thick)
    from physim.nuclear import data
    from physim.nuclear.kinematics import TwoBody
    tb = TwoBody(pt.beam, pt.target, pt.beam_energy, excitation_mev=0.328, excite="recoil")
    assert deorientation_factors(thick, pt, "recoil", tb, data.nuclide("194Pt").nuclear_mass_mev,
                                 [0.0, PT_TAU * 1e-12]) == {1: pytest.approx(np.ones((90, 5)))}
    a, b = correlation_table(thick, nucleus="other"), correlation_table(plain, nucleus="other")
    assert [r["factor"] for r in a] == pytest.approx([r["factor"] for r in b], rel=1e-9)
    a, b = correlation_table(thick), correlation_table(plain)
    assert any(abs(x["factor"] - y["factor"]) > 0.01 for x, y in zip(a, b))
