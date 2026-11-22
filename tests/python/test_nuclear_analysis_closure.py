"""Simulate, then analyse: the analysis gets B(E2) back on a thick target in inverse kinematics (Allmond et al.,
PRC 90, 034309 (2014), fixes 1, 2 and 7 of its verification).

62Ni at 111.6 MeV on 1.04 mg/cm² of carbon loses 38 % of its energy in the target, and the excitation probability
changes 150-fold through it; E_cm is only ~15 MeV against E* = 1.17 MeV, so the excited carbon recoil's CM angle
and lab-to-CM Jacobian differ from the elastic one's (by 4 % in the Jacobian). Before the analysis took ⟨P⟩/B
from the forward model it read B(E2) +9.7/+8.5/+8.8 % high (P at mid-target, elastic kinematics, before the forward
model integrated P over the depth), and +35 % once it did."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.analysis import GATE_OVERLAP_WARN, Analysis, Settings
from physim.nuclear.gamma_events import simulate_gammas

B_UP = 884.105  # e²fm⁴, ENSDF
#: Events of the weighted sample. With 2 000 000 the Monte Carlo uncertainty of B(E2) (Result.monte_carlo, the
#: scatter between seeds) is about 0.6 % for the inner ring, 0.35 % for the outer one and 0.3 % for both; the
#: excited share of the particles in the rings is known to 0.35 %, 0.22 % and 0.18 %.
EVENTS = 2_000_000


def allmond_setup() -> Experiment:
    """The 62Ni on natural-carbon setup of the Allmond 2014 verification (13C left out): carbon recoils in two
    full-azimuth rings at 14–28° and 28–44°, nine clovers at 21.75 cm sized to 2.24 % at 1173 keV."""
    ring = {"shape": "annular", "theta": "0 deg", "distance": "50 mm", "rings": 7, "thickness": "1000 um",
            "dead_layer": "0.5 um", "resolution": "500 keV", "threshold": "1 MeV"}
    clover = {"model": "clover", "distance": "217.5 mm", "crystals": 4, "crystal_diameter": "54.5131 mm",
              "crystal_length": "70 mm", "crystal_pitch": "49.0618 mm", "housing_side": "109.1237 mm",
              "window_gap": "5 mm", "resolution": "2.1 keV", "material": "Ge", "addback": True}
    places = [(90, 0), (90, 72), (90, 144), (90, 216), (90, 288), (132, 60), (132, 180), (132, 300), (154, 0)]
    return Experiment.from_dict({
        "schema": "physim.experiment/1",
        "title": "62Ni on natC (Allmond 2014)",
        "reaction": {"type": "coulex", "excite": "projectile", "energy": "1172.91 keV", "multipolarity": "E2",
                     "b_up": f"{B_UP} e2fm4"},
        "beam": {"nuclide": "62Ni", "energy": "111.6 MeV", "current": "5 ppA", "charge_state": 10},
        "target": {"material": "12C", "thickness": "1.03741 mg/cm2"},
        "run": {"beam_time": "24 h", "counts_wanted": 1000},
        "detectors": [dict(ring, name="Ring2", inner_radius="12.4664 mm", outer_radius="26.5855 mm", sectors=10),
                      dict(ring, name="Ring3", inner_radius="26.5855 mm", outer_radius="48.2844 mm", sectors=12)],
        "gamma_detectors": [dict(clover, name=f"G{t}_{p}", theta=f"{t} deg", phi=f"{p} deg") for t, p in places],
    })


@pytest.fixture(scope="module")
def allmond():
    exp = allmond_setup()
    return exp, Analysis(exp, simulate_gammas(exp, EVENTS, seed=11))


def test_the_expected_excitations_per_particle_are_those_of_the_events(allmond):
    """⟨P⟩/B comes from the forward cross sections (through the depth, with the excited recoil's kinematics): it
    is the share of the particles in the rings that the simulation excited, per unit B(E2), within the sample's
    own uncertainty."""
    exp, a = allmond
    ev = a.events
    c = ev.columns
    excited = np.isin(c["channel"], [i for i, lab in enumerate(ev.channels) if "excited" in lab])
    for det in ("Ring2", "Ring3", None):
        s = Settings(detectors=None if det is None else [det], particle_gate="all")
        rows = a._elastic_rows(s)
        w = c["weight"]
        share = w[rows & excited].sum() / w[rows].sum() / B_UP
        unc = math.sqrt(float((w[rows & excited] ** 2).sum())) / float(w[rows & excited].sum())
        p = a._mean_probability(s)
        assert p == pytest.approx(share, rel=max(3 * unc, 0.004)), (det, p, share, unc)


def test_b_e2_comes_back_per_ring_and_for_both(allmond):
    exp, a = allmond
    for det in ("Ring2", "Ring3", None):
        r = a.run(Settings(detectors=None if det is None else [det], particle_gate="all"))
        dev = r.b_e2fm4 / B_UP - 1
        assert r.monte_carlo < 0.007, (det, r.monte_carlo)
        # A residual of about +0.5 to +1 % remains, not from the Monte Carlo: seed 11 with 2 000 000 events gives
        # +0.9/+0.0/+0.2 % (MC 0.6/0.35/0.3 %), seed 7 with 4 000 000 gives +1.9/+0.6/+0.9 % (MC 0.4/0.25/0.2 %).
        # Of it, about 0.3 % is ⟨P⟩/B against the events' own excited share (the test above), and 0.3 to 0.6 % the
        # efficiency × angular-correlation factor against the γ rays the simulation detects per excitation (the
        # correlation table takes the kinematics of the middle of the target: gamma.correlation_table). So: 1 %
        # plus three Monte Carlo standard deviations.
        assert abs(dev) < 0.01 + 3 * r.monte_carlo, (det, 100 * dev, 100 * r.monte_carlo)
        assert r.gate == "all" and r.gate_warning is None and r.gate_overlap is None


def test_the_inelastic_gate_warns_when_the_groups_overlap(allmond):
    """Through 38 % energy loss the elastic and inelastic carbon groups overlap: the inelastic gate keeps a
    fraction of the excitations from the depths it happens to pick out, and the result says so and points to
    the gate "all". It does not switch by itself."""
    exp, a = allmond
    r = a.run(Settings(particle_gate="inelastic"))
    assert r.gate == "inelastic" and r.settings["particle_gate"] == "inelastic"
    assert r.gate_overlap > GATE_OVERLAP_WARN
    assert r.gate_warning and 'gate "all"' in r.gate_warning and "inelastic particle gate" in r.gate_warning
    assert r.notes[0] == r.gate_warning
    assert "overlap" in r.steps[0]["text"]
    # The gate still did what was asked: fewer γ rays than with every particle.
    assert r.lines[0]["gate_efficiency"] < 0.5


def test_no_warning_when_the_groups_are_apart():
    """A thin target and a fine detector: the groups stand apart and the inelastic gate is trusted."""
    d = Experiment.example("coulex_ni58").to_dict()
    d["target"]["thickness"] = "0.02 mg/cm2"
    d["detectors"] = [x for x in d["detectors"] if x["name"] == "DSSD-L"]
    for x in d["detectors"]:
        x["resolution"] = "20 keV"
    exp = Experiment.from_dict(d)
    a = Analysis(exp, simulate_gammas(exp, 60_000, seed=3))
    r = a.run(Settings(particle_gate="inelastic"))
    assert r.gate == "inelastic" and r.gate_overlap < GATE_OVERLAP_WARN and r.gate_warning is None, r.gate_overlap
