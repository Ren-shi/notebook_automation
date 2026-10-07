"""Cascades (backlog item 71): each excitation decays through the level scheme's branches, so two populated levels
give true γ–γ coincidences in the numbers the populations, branches and efficiencies predict; one level gives
item 55's result unchanged."""

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import dataviews as dv
from physim.nuclear.events import simulate
from physim.nuclear.gamma import cascade_excitation
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.levels import Level, LevelScheme, Transition, Value
from physim.nuclear.orientation import simple_scheme
from physim.nuclear.planner import Planner

E1, E2 = 1454.2, 2775.4  # keV: the 2⁺₁ and 2⁺₂ states of ⁵⁸Ni


def v(x):
    return Value(x, source="user")


def two_states() -> LevelScheme:
    """⁵⁸Ni with two 2⁺ states excited from the ground state; the upper one decays to the lower (60%) or to the
    ground state (40%)."""
    s = LevelScheme("58Ni")
    s.levels = [Level(v(0.0), 0.0, 1), Level(v(E1), 2.0, 1), Level(v(E2), 2.0, 1)]
    s.transitions = [Transition(1, 0, v(E1), v(100.0), "E2"), Transition(2, 1, v(E2 - E1), v(60.0), "E2"),
                     Transition(2, 0, v(E2), v(40.0), "E2")]
    s.set_matrix_element(0, 1, "E2", 26.4)
    s.set_matrix_element(0, 2, "E2", 80.0)  # far stronger than in ⁵⁸Ni: about 9% of the excitations, for statistics
    return s


def setup(scheme) -> Experiment:
    p = Planner.example("coulex_ni58")
    d = p.draft
    d["detectors"] = [x for x in d["detectors"] if x["name"] == "CD"]
    d["gamma_detectors"] = [gd for gd in d["gamma_detectors"] if gd["name"] != "Ge135"]
    for gd in d["gamma_detectors"]:
        gd["distance"] = "60 mm"
    d["reaction"]["emission"] = "isotropic"  # the efficiencies alone set the coincidences
    if scheme is not None:
        d["levels"] = {"target": scheme.to_dict()}
    return Experiment.from_dict(d)


def test_one_level_gives_item_55s_result():
    plain = setup(None)
    one = setup(simple_scheme("58Ni", "1.4542 MeV", "E2", plain.excitation.b_up))
    assert cascade_excitation(one) is None
    a = simulate_gammas(plain, 50_000, 5)
    b = simulate_gammas(one, 50_000, 5)
    for k in ("measured", "crystal", "particle", "weight"):
        assert np.array_equal(a[k], b[k]), k


def test_two_levels_give_the_gamma_gamma_coincidences_the_branches_predict():
    exp = setup(two_states())
    cex = cascade_excitation(exp)
    assert cex is not None and set(cex.transitions()) == {(1, 0), (2, 1), (2, 0)}
    ev = simulate(exp, 400_000, 2)
    g = simulate_gammas(exp, particle_events=ev, gammas_per_event=20, seed=2)
    lines = {(2, 1): (E2 - E1) * 1e-3, (1, 0): E1 * 1e-3, (2, 0): E2 * 1e-3}
    assert {(int(i), int(f)) for i, f in zip(g["initial"], g["final"])} == set(lines)

    # The excitations with a detected particle that start in the upper level: by its share of the direct
    # excitation probability at each event's CM angle.
    c = ev.columns
    excited = np.flatnonzero(np.isin(c["channel"], [i for i, lab in enumerate(ev.channels) if "excited" in lab]))
    grid, states = cex.table()
    p1 = np.interp(c["theta_cm"][excited], grid, [s.direct.get(1, 0.0) for s in states])
    p2 = np.interp(c["theta_cm"][excited], grid, [s.direct.get(2, 0.0) for s in states])
    upper = float(np.sum(c["weight"][excited] * p2 / (p1 + p2)))   # per second
    branch = cex.decay(2, 1).gamma
    # The two γ rays in two different crystals (single-crystal detectors here): in one crystal they would sum
    # into one signal (backlog item 72) and count in neither peak.
    eff = {e: [gd.peak_efficiency(e, exp) for gd in exp.gamma_detectors] for e in (lines[(2, 1)], lines[(1, 0)])}
    a, b = eff[lines[(2, 1)]], eff[lines[(1, 0)]]
    pairs = sum(a[i] * b[j] for i in range(len(a)) for j in range(len(b)) if i != j)
    expected = upper * branch * pairs * 2

    def window(e):  # the Doppler-corrected full-energy peak
        return (e - 0.012, e + 0.012)

    rate, err = g.gamma_gamma(window(lines[(2, 1)]), window(lines[(1, 0)]), corrected="recoil")
    expected /= 2  # each cascade counts once
    assert rate > 0
    assert err < 0.25 * rate, "enough pairs to test"
    assert abs(rate - expected) < 3 * err + 0.05 * expected, (rate, err, expected)

    # Gating on the 1321 keV line selects the 1454 keV line, and not the 2775 keV one (which feeds the ground).
    corrected = g["corrected_recoil"]
    cas = g["cascade"]
    gated = set(cas[(np.abs(corrected - lines[(2, 1)]) < 0.012) & g["counted"]].tolist())
    with_gate = np.isin(cas, list(gated)) & g["counted"]
    near = lambda e: np.abs(corrected - e) < 0.012  # noqa: E731
    assert (with_gate & near(lines[(1, 0)])).sum() > 0
    assert (with_gate & near(lines[(2, 0)])).sum() == 0
    co = g.coincidences()
    assert sum(co["gamma_gamma_true"].values()) > 0


def test_a_run_keeps_the_cascades_and_the_data_tab_gates_on_a_line(tmp_path):
    p = Planner(setup(two_states()))
    p.create_experiment(root=tmp_path)
    run = p.start_run("beam", duration="10 min", budget="3 min")
    g = run.gammas()
    assert set(np.unique(g["initial"]).tolist()) >= {1, 2}
    lo, hi = (E2 - E1 - 12) * 1e-3, (E2 - E1 + 12) * 1e-3
    name = g.detector_names()[0]
    gated = dv.gamma_spectrum(run, name, correction="recoil", randoms="none", gamma_gate=(lo, hi), bins=300,
                              range=(0.5, 3.0))
    free = dv.gamma_spectrum(run, name, correction="recoil", randoms="none", bins=300, range=(0.5, 3.0))
    assert gated["counts"].sum() < free["counts"].sum()
    assert "γ gate" in gated["label"]
