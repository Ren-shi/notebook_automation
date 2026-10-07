"""Lifetimes (backlog item 73): a state that lives a few picoseconds decays partly after its nucleus has stopped in
a thick target and backing; the line shape's stopped share equals the share of decays after stopping, and a short
lifetime gives item 55's result."""

import numpy as np
import pytest

from physim.nuclear import Experiment, data
from physim.nuclear.events import simulate
from physim.nuclear.gamma_events import _clock, simulate_gammas
from physim.nuclear.levels import Value
from physim.nuclear.orientation import simple_scheme
from physim.nuclear.planner import Planner
from physim.nuclear.rates import _exit_paths, stack

E0 = 1.4542  # MeV


def setup(half_life=None, thickness="1 mg/cm2", backing=("Au", "20 mg/cm2")) -> Experiment:
    """The ⁵⁸Ni example with the CD only, close germanium detectors, isotropic emission, and the 2⁺ state's
    half-life ``half_life`` (s) in the setup's level scheme."""
    p = Planner.example("coulex_ni58")
    d = p.draft
    d["detectors"] = [x for x in d["detectors"] if x["name"] == "CD"]
    d["gamma_detectors"] = [gd for gd in d["gamma_detectors"] if gd["name"] != "Ge135"]
    for gd in d["gamma_detectors"]:
        gd["distance"] = "60 mm"
    d["reaction"]["emission"] = "isotropic"
    d["target"]["thickness"] = thickness
    if backing:
        d["target"]["backing"] = {"material": backing[0], "thickness": backing[1]}
    if half_life is not None:
        s = simple_scheme("58Ni", "1.4542 MeV", "E2", d["reaction"]["b_up"])
        s.levels[1].half_life = Value(half_life, source="user")
        d["levels"] = {"target": s.to_dict()}
    return Experiment.from_dict(d)


def test_the_stopped_share_is_the_share_of_decays_after_stopping():
    half = 1.0e-12  # s: comparable with the ~1 ps the recoils take to stop in the gold
    exp = setup(half)
    ev = simulate(exp, 600_000, 4)
    g = simulate_gammas(exp, particle_events=ev, gammas_per_event=150, seed=4)
    stopped = g["stopped"]
    assert 0.1 < stopped.mean() < 0.9

    # The share of decays after stopping, from the stopping time of each recoil and the exponential decay law.
    c = ev.columns
    excited = np.flatnonzero(np.isin(c["channel"], [i for i, lab in enumerate(ev.channels) if "excited" in lab]))
    from physim.nuclear.kinematics import TwoBody

    layers = stack(exp)
    ex = TwoBody("16O", "58Ni", 30.0, excitation_mev=E0, excite="recoil")
    rec = ex.recoil_for(c["theta_cm"][excited])
    e_rec = np.nan_to_num(np.asarray(rec.energy, dtype=float))
    th, ph = np.radians(np.asarray(rec.theta_lab)), np.radians(np.where(c["recoil"][excited], c["phi"][excited],
                                                                         c["phi"][excited] + 180.0))
    dirs = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    fwd, steps = _exit_paths(exp, layers, 0, c["depth"][excited], dirs)
    mass = data.nuclide("58Ni").nuclear_mass_mev + E0
    t_stop = np.full(len(excited), np.inf)
    e = e_rec.copy()
    elapsed = np.zeros(len(e))
    from physim.nuclear.rates import _after, stopping

    for j, path, forward in steps:
        m = fwd == forward
        mat = layers[j].material
        ge, gt = _clock("58Ni", mat, float(mat.density_g_cm3), mass)
        e_out = _after(stopping("58Ni", mat), e, np.asarray(path))
        stop_here = m & (e_out <= 0) & np.isinf(t_stop)
        t_stop[stop_here] = elapsed[stop_here] + np.interp(e[stop_here], ge, gt)
        go = m & np.isinf(t_stop)
        elapsed[go] += np.interp(e[go], ge, gt) - np.interp(e_out[go], ge, gt)
        e[go] = e_out[go]
    tau = half / np.log(2)
    after_stopping = np.where(np.isfinite(t_stop), np.exp(-t_stop / tau), 0.0)
    expected = float(np.average(after_stopping, weights=c["weight"][excited]))
    # Over the γ rays emitted (each excitation once), weighted as the events are.
    assert stopped.mean() == pytest.approx(expected, abs=0.08)

    # The line shape: in the forward detectors the stopped component sits at E0, the moving one 30 keV and more
    # above it.
    fwd_det = g.select("Ge0") | g.select("Ge45")
    full = fwd_det & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-9)
    raw = g["measured"][full]
    # In a forward detector every moving decay is shifted up, so what lies below E0 is the lower half of the
    # stopped peak (a Gaussian of the detector's resolution): twice its share is the stopped share.
    below = raw < E0
    share_in_line_shape = 2 * float(np.sum(g["weight"][full][below]) / np.sum(g["weight"][full]))
    # The line shape cannot tell a stopped nucleus from one that decays in the last, slow stretch of its slowing
    # down with a Doppler shift below the detector's resolution: those count as stopped in it.
    sigma = 2.5e-3 / 2.3548  # the germanium's resolution at 1.33 MeV, as a standard deviation (MeV)
    at_rest = g["stopped"][full] | (g["energy_lab"][full] - E0 < sigma)
    truth = float(np.sum(g["weight"][full][at_rest]) / np.sum(g["weight"][full]))
    n = full.sum()
    assert n > 150
    assert share_in_line_shape == pytest.approx(truth, abs=3 * np.sqrt(truth * (1 - truth) / n) + 0.03)


def test_a_short_lifetime_gives_item_55s_result():
    # No backing and a thin target: decaying at once or after the target is the same.
    flat = setup(None, thickness="0.2 mg/cm2", backing=None)
    short = setup(1e-16, thickness="0.2 mg/cm2", backing=None)
    a = simulate_gammas(flat, 300_000, 6, gammas_per_event=10)
    b = simulate_gammas(short, 300_000, 6, gammas_per_event=10)
    assert not b["stopped"].any()
    for name in ("Ge0", "Ge45"):
        ma, mb = a.select(name), b.select(name)
        ca = np.average(a["measured"][ma], weights=a["weight"][ma])
        cb = np.average(b["measured"][mb], weights=b["weight"][mb])
        assert cb == pytest.approx(ca, abs=0.002), name  # 2 keV on 1454
    assert len(b) == pytest.approx(len(a), rel=0.05)
