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
#: Events and seeds of the weighted samples. They are averaged: five samples of 400 000 know what one of 2 000 000
#: does (B(E2) to about 0.6 % in the inner ring, 0.35 % in the outer one and 0.3 % in both; the excited share of the
#: particles to 0.35 %, 0.22 % and 0.18 %) and need 3 GB at the most, where that one needs about 10 GB.
EVENTS = 400_000
SEEDS = (21, 22, 23, 24, 25)


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
    return exp, [Analysis(exp, simulate_gammas(exp, EVENTS, seed=seed)) for seed in SEEDS]


def mean_and_error(values, errors) -> tuple:
    """The mean of the samples' values, and its error from theirs."""
    return float(np.mean(values)), math.sqrt(float(np.sum(np.square(errors)))) / len(errors)


def excited_rows(a) -> np.ndarray:
    return np.isin(a.events.columns["channel"], [i for i, lab in enumerate(a.events.channels) if "excited" in lab])


def every_particle(det) -> Settings:
    return Settings(detectors=None if det is None else [det], particle_gate="all")


@pytest.mark.parametrize("det", ["Ring2", "Ring3", None])
def test_the_expected_excitations_per_particle_are_those_of_the_events(allmond, det):
    """⟨P⟩/B comes from the forward cross sections (through the depth, with the excited recoil's kinematics): it
    is the share of the particles in the rings that the simulation excited, per unit B(E2), within the samples'
    own uncertainty."""
    s = every_particle(det)
    dev, err = [], []
    for a in allmond[1]:
        rows, w = a._elastic_rows(s), a.events.columns["weight"]
        hit = rows & excited_rows(a)
        dev.append(a._mean_probability(s) / (w[hit].sum() / w[rows].sum() / B_UP) - 1)
        err.append(math.sqrt(float((w[hit] ** 2).sum())) / float(w[hit].sum()))
    mean, unc = mean_and_error(dev, err)
    assert abs(mean) < max(3 * unc, 0.004), (det, mean, unc)


@pytest.mark.parametrize("det", ["Ring2", "Ring3", None])
def test_the_efficiency_and_correlation_are_those_of_the_events(allmond, det):
    """The efficiency × angular-correlation factor is the full-energy γ rays the simulation detects per excitation.
    Over twelve samples it is 0.0 ± 0.4 %, 0.0 ± 0.2 % and +0.1 ± 0.2 % off (with the correlation table and the
    events' emitters at the middle of the target it was −0.9, −1.1 and −1.0 % off)."""
    s = every_particle(det)
    dev, err = [], []
    for a in allmond[1]:
        g, w = a.gammas, a.events.columns["weight"]
        keep, _ = a._gate(s)
        eff, corr = a._efficiency_and_correlation(keep)
        full = keep & g["counted"] & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-6)
        seen = g["weight"][full]
        dev.append(eff * corr / (seen.sum() / w[a._elastic_rows(s) & excited_rows(a)].sum()) - 1)
        err.append(math.sqrt(float((seen**2).sum())) / float(seen.sum()))
    mean, unc = mean_and_error(dev, err)
    assert unc < 0.006 and abs(mean) < 2.5 * unc, (det, mean, unc)


@pytest.mark.parametrize("det", ["Ring2", "Ring3", None])
def test_b_e2_comes_back_per_ring_and_for_both(allmond, det):
    """Within the Monte Carlo error alone. Over twelve samples −0.1 ± 0.4 %, +0.1 ± 0.2 % and −0.1 ± 0.2 %. The peak
    has the tails that the depth of the target gives it: with the area of a Gaussian alone it was −1.6, −3.2 and
    −2.7 %."""
    results = [a.run(every_particle(det)) for a in allmond[1]]
    dev, mc = mean_and_error([r.b_e2fm4 / B_UP - 1 for r in results], [r.monte_carlo for r in results])
    assert mc < 0.007, (det, mc)
    assert abs(dev) < 2.5 * mc, (det, 100 * dev, 100 * mc)
    assert all(r.gate == "all" and r.gate_warning is None and r.gate_overlap is None for r in results)
    assert all(r.fit.tails is not None and r.fit.tails["inside"] < 0.99 for r in results)


def test_the_inelastic_gate_warns_when_the_groups_overlap(allmond):
    """Through 38 % energy loss the elastic and inelastic carbon groups overlap: the inelastic gate keeps a
    fraction of the excitations from the depths it happens to pick out, and the result says so and points to
    the gate "all". It does not switch by itself."""
    a = allmond[1][0]
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
