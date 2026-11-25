"""Simulate, then analyse: the analysis gets B(E2) back on a thick target in inverse kinematics (Allmond et al.,
PRC 90, 034309 (2014), fixes 1, 2 and 7 of its verification).

62Ni at 111.6 MeV on 1.04 mg/cm² of carbon loses 38 % of its energy in the target, and the excitation probability
changes 150-fold through it; E_cm is only ~15 MeV against E* = 1.17 MeV, so the excited carbon recoil's CM angle
and lab-to-CM Jacobian differ from the elastic one's (by 4 % in the Jacobian). Before the analysis took ⟨P⟩/B
from the forward model it read B(E2) +9.7/+8.5/+8.8 % high (P at mid-target, elastic kinematics, before the forward
model integrated P over the depth), and +35 % once it did."""

import json
import math
from dataclasses import asdict, replace
from types import SimpleNamespace

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import analysis as analysis_module
from physim.nuclear.analysis import GATE_OVERLAP_WARN, Analysis, PeakFit, Settings
from physim.nuclear.coulex_reproduction import setup as reproduction_setup
from physim.nuclear.gamma_events import doppler_spread, reaction_depths, simulate_gammas

B_UP = 884.105  # e²fm⁴, ENSDF
#: Events and seeds of the weighted samples. They are averaged: five samples of 400 000 know what one of 2 000 000
#: does (B(E2) to about 0.6 % in the inner ring, 0.35 % in the outer one and 0.3 % in both; the excited share of the
#: particles to 0.35 %, 0.22 % and 0.18 %) and need 3 GB at the most, where that one needs about 10 GB.
EVENTS = 400_000
SEEDS = (21, 22, 23, 24, 25)


def allmond_setup() -> Experiment:
    """The 62Ni on natural-carbon setup of the Allmond 2014 verification (13C left out): carbon recoils in two
    full-azimuth rings at 14–28° and 28–44°, nine clovers at 21.75 cm sized to 2.24 % at 1173 keV.

    Without the E1 polarisation: the samples are those that the numbers below are of. (With it, 0.5 % fewer excitations,
    the samples of the same seeds are others. Of all that is tested here one number then leaves its band, the fitted
    area in ring 3, +0.28 ± 0.03 % in these five; over 21 samples, 21–25 and 301–316, it is +0.13 ± 0.04 % with and +0.07 ±
    0.04 % without. The closure with it is test_nuclear_allmond2014 and test_nuclear_polarisation.)"""
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
                     "b_up": f"{B_UP} e2fm4", "e1_polarisation": 0},
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


@pytest.mark.parametrize("det, beyond", [("Ring2", 0.008), ("Ring3", 0.020), (None, 0.017)])
def test_b_e2_comes_back_per_ring_and_for_both(allmond, det, beyond):
    """Within the Monte Carlo error alone. Over twelve samples −0.1 ± 0.4 %, +0.1 ± 0.2 % and −0.1 ± 0.2 %. The peak
    has the tails that the depth of the target gives it: with the area of a Gaussian alone it was −1.6, −3.2 and
    −2.7 %.

    The tails are used: every sample shows them, and over the samples the peak has more than ``beyond`` of itself
    outside ±3σ (a Gaussian has 0.27 %). One sample knows that share to 0.24, 0.16 and 0.13 % (the scatter of twelve:
    the inner ring has few entries there), so the bound is on the mean, three of its errors under the 1.18, 2.30 and
    1.87 % of these five with the segment centre of the Doppler correction as it is or as it is being corrected (1.24,
    2.35, 1.88 % and 1.18, 2.30, 1.87 %)."""
    results = [a.run(every_particle(det)) for a in allmond[1]]
    dev, mc = mean_and_error([r.b_e2fm4 / B_UP - 1 for r in results], [r.monte_carlo for r in results])
    assert mc < 0.007, (det, mc)
    assert abs(dev) < 2.5 * mc, (det, 100 * dev, 100 * mc)
    assert all(r.gate == "all" and r.gate_warning is None and r.gate_overlap is None for r in results)
    assert all(r.fit.tails is not None for r in results)
    assert np.mean([1 - r.fit.tails["inside"] for r in results]) > beyond


def full_energy_counts(a, keep) -> tuple:
    """The simulation's own counts in the full-energy peak, and their relative Monte Carlo error."""
    g = a.gammas
    w = g["weight"][keep & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-6)] * g.events.beam_time_s * g.live_fraction
    return float(w.sum()), math.sqrt(float((w**2).sum())) / float(w.sum())


@pytest.mark.parametrize("det", ["Ring2", "Ring3", None])
def test_the_fitted_area_is_the_full_energy_counts(allmond, det):
    """Round 3, fix 4. The peak has a flat top and tails that end: the share that two Gaussians have within ±3σ lost
    of it; the counts over the background of the whole window depend on it less. A wider Gaussian still took +0.26,
    +0.18 and +0.22 % from the background below the peak once the Doppler correction had the right speed; with a
    wider bell that ends it is +0.07 ± 0.10, +0.04 ± 0.05 and +0.04 ± 0.04 % (16 samples, seeds 301–316)."""
    s = replace(every_particle(det), subtract_randoms=False)
    dev = []
    for a in allmond[1]:
        true, _ = full_energy_counts(a, a._gate(s)[0])
        dev.append(a.run(s).area / true - 1)
    mean, err = float(np.mean(dev)), float(np.std(dev, ddof=1)) / math.sqrt(len(dev))
    assert abs(mean) < 0.001 + 2.5 * err and err < 0.003, (det, mean, err)


# -- decays inside the target -----------------------------------------------------------------------------------------


def with_half_life(text) -> Experiment:
    """The reproduction's setup (the two-level scheme, all orders) with that half-life of the 2⁺ state, or without."""
    d = reproduction_setup(levels="two").to_dict()
    level = d["levels"]["beam"]["level"][1]
    level.pop("half_life")
    if text is not None:
        level["half_life"] = text
    return Experiment.from_dict(d)


def spread_of(text) -> tuple:
    """(shifts in keV, share of each) of the peak in the inner ring, all crystals alike."""
    found = doppler_spread(with_half_life(text), {"Ring2": np.ones(36)}, 1e-4)
    return None if found is None else (1172.91 * found[0], found[1])


def test_the_spread_of_the_peak_follows_the_half_life():
    """Nuclei that decay inside the target are faster than the Doppler correction takes them, and the crystals are at
    90° and behind: a tail below the peak. None without a half-life; the peak of a level that lives long enough to
    leave is symmetric about the transition; one that decays at once is far below it."""
    assert spread_of(None) is None
    x, p = spread_of("1 us")
    tail = float(p[(x > -20) & (x < -7)].sum())
    assert p.sum() == pytest.approx(1) and abs(float((x * p).sum())) < 0.2 and tail < 0.02
    assert tail == pytest.approx(float(p[(x > 7) & (x < 20)].sum()), abs=0.005)
    y, q = spread_of("1.442 ps")
    assert float(q[(y > -20) & (y < -7)].sum()) > tail + 0.03
    assert -0.6 < float((y * q).sum()) < -0.3
    z, r = spread_of("1 fs")
    assert float((z * r).sum()) < -4
    # Through a thin target hardly any decays inside: what is left is as far above as below, and it ends.
    thin = with_half_life("1.442 ps").to_dict()
    thin["target"]["thickness"] = "0.01 mg/cm2"
    thin = Experiment.from_dict(thin)
    assert reaction_depths(thin) is None
    x, p = doppler_spread(thin, {"Ring2": np.ones(36)}, 1e-4)
    x = 1172.91 * x
    assert p.sum() == pytest.approx(1) and abs(float((x * p).sum())) < 0.2
    assert float(p[(x > -20) & (x < -7)].sum()) == pytest.approx(float(p[(x > 7) & (x < 20)].sum()), abs=0.005)


@pytest.fixture(scope="module")
def decaying():
    """Samples with the 2⁺ half-life (1.442 ps) in the level scheme."""
    exp = with_half_life("1.442 ps")
    return exp, [Analysis(exp, simulate_gammas(exp, EVENTS, seed=seed)) for seed in (11301, 11302, 11303, 11304)]


@pytest.mark.parametrize("det, scatter, of_tail", [("Ring2", 0.0031, 0.027), ("Ring3", 0.0017, 0.019),
                                                   (None, 0.0013, 0.016)])
def test_the_area_and_b_e2_come_back_with_the_half_life(decaying, det, scatter, of_tail):
    """The tail of the decays inside the target is in the area, and B(E2) comes back. ``scatter``: of the area of one
    sample; ``of_tail``: of the tail of four, relative to the fitted one (55 groups of four seeds, whose mean is right to
    0.6 %).

    Over 160 (seeds 17001–17160) the area is −0.08 ± 0.03, −0.03 ± 0.01 and −0.05 ± 0.01 % off (with a wider bell in place of
    the spread it was −0.40 ± 0.05 % in the inner ring), and B(E2) +0.15 ± 0.11, −0.06 ± 0.07 and −0.09 ± 0.06 % (with the
    orientation of the middle of the target in every event +0.36 ± 0.19, −0.29 ± 0.10 and −0.18 ± 0.09 % over sixty)."""
    exp, samples = decaying
    s = every_particle(det)
    area, b, mc, low, seen = [], [], [], 0.0, 0.0
    for a in samples:
        keep = a._gate(s)[0]
        true, _ = full_energy_counts(a, keep)
        fit = a.run(replace(s, subtract_randoms=False)).fit
        assert "spread" in fit.tails and not fit.tails["simulated"]
        area.append(fit.area / true - 1)
        # The tail of the fitted peak, and that of the events' own full-energy counts.
        g = a.gammas
        x = 1172.91 + np.arange(-19.975, -8, 0.05)
        curve = fit.curve(1e-3 * x) - fit.background_curve(1e-3 * x)
        low += 0.05e-3 / fit.tails["width"] * float(curve.sum()) / len(samples)
        tail = (keep & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-6)
                & (np.abs(1e3 * g["corrected_projectile"] - 1172.91 + 14) < 6))
        seen += float(g["weight"][tail].sum()) * g.events.beam_time_s * g.live_fraction / len(samples)
        r = a.run(s)
        b.append(r.b_e2fm4 / exp.excitation.b_up_e2fm - 1)
        mc.append(r.monte_carlo)
        assert 0 < r.budget["peak tails (simulated)"] < 0.02 and any("decay inside the target" in n for n in r.notes)
    assert low == pytest.approx(seen, rel=0.03 + 3 * of_tail), det
    assert abs(np.mean(area)) < 0.001 + 3 * scatter / math.sqrt(len(area)), (det, area)
    mean, unc = mean_and_error(b, mc)
    assert abs(mean) < 0.0025 + 3 * unc, (det, mean, unc)


# -- a run of few counts ----------------------------------------------------------------------------------------------


def as_run(a, hours: float, seed: int) -> Analysis:
    """A run of ``hours`` drawn from a weighted sample's γ rays: each as often as its rate says, all of weight one
    over the run's time. (The particles keep their weights: they are only counted.)"""
    g = a.gammas
    real = 3600.0 * hours
    picked = np.random.default_rng(seed).poisson(g["weight"] * real)
    rows = np.repeat(np.arange(len(g)), picked)
    cols = {k: np.asarray(v)[rows] for k, v in g.columns.items()}
    cols["weight"] = np.full(len(rows), 1.0 / real)
    run = SimpleNamespace(summary={"one_particle": True}, scale=1.0, real_s=real, duration_s=real)
    return Analysis(a.experiment, replace(g, columns=cols, events=replace(g.events, beam_time_s=real), seed=seed),
                    run=run)


@pytest.fixture(scope="module")
def hours(allmond):
    """Seven hours of beam: about 400 counts in the inner ring, too few to fix the tails."""
    return as_run(allmond[1][0], 7.0, 5)


def test_a_run_of_few_counts_takes_its_tails_from_a_simulation(hours, monkeypatch):
    """Round 3, fix 3: the tail shape comes from a simulation of the run's own setup, once; the step, a note and a line
    of the budget say so."""
    s = replace(every_particle("Ring2"), subtract_randoms=False)
    r = hours.run(s)
    assert 250 < r.area < 600
    t = r.fit.tails
    assert t is not None and t["simulated"] and r.fit.failed is None and r.fit.shape.startswith("Gaussian with simulated")
    step = next(x["text"] for x in r.steps if x["step"] == "peak fit")
    assert "simulation of this setup" in step and any("simulation of this setup" in n for n in r.notes)
    assert 0 < r.budget["peak tails (simulated)"] < 0.02
    assert r.systematic == pytest.approx(math.sqrt(sum(v**2 for v in r.budget.values())))
    assert PeakFit.of(json.loads(json.dumps(asdict(r.fit)))).tails["simulated"]
    # The simulated peak: many more counts, its own tails, and a width like the run's.
    sim = hours.simulated_peak(s)
    assert sim.tails is not None and not sim.tails["simulated"]
    assert t["share"] == sim.tails["share"] and r.fit.sigma == pytest.approx(sim.sigma, rel=0.2)
    assert r.budget["peak tails (simulated)"] >= analysis_module.SHAPE_KNOWN * t["simulated_loss"]
    # A second fit, and a change that leaves the peak alone, do not simulate again.
    monkeypatch.setattr(Analysis, "_simulated_gammas", lambda *a: pytest.fail("simulated again"))
    assert hours.run(s).area == r.area
    hours.run(replace(s, efficiency_unc=0.1))


def test_a_weighted_sample_is_not_simulated_again(allmond, monkeypatch):
    monkeypatch.setattr(Analysis, "_simulated_gammas", lambda *a: pytest.fail("simulated again"))
    r = allmond[1][0].run(every_particle("Ring3"))
    assert "peak tails (simulated)" not in r.budget and not r.fit.tails["simulated"]


def test_a_failed_fit_has_no_area(hours, monkeypatch):
    """A fit that ends at a bound is reported as failed: no area, no B(E2), and the result says why."""
    real = analysis_module.fit_tails

    def failing(*args, **kw):
        return replace(real(*args, **kw), failed="the peak's width ran to 2 times more than the setup gives", area=1e6)

    monkeypatch.setattr(analysis_module, "fit_tails", failing)
    r = hours.run(replace(every_particle("Ring2"), subtract_randoms=False))
    assert r.fit.failed and r.area == 0 and r.fit.area == 0 and r.b_e2fm4 != r.b_e2fm4 and not r.budget
    assert "failed" in next(x["text"] for x in r.steps if x["step"] == "peak fit")
    assert any(n.startswith("The peak fit failed") for n in r.notes)


# -- a run knows segments, not angles --------------------------------------------------------------------------------


def stored(a, setup=None) -> Analysis:
    """The particles of a weighted sample as a run keeps them: the counted ones, at their segment's mean angles."""
    from physim.nuclear.runs import _segment_directions

    c = a.events.columns
    rows = np.flatnonzero(c["counted"])
    cols = {k: np.asarray(v)[rows] for k, v in c.items()}
    cols["theta"], cols["phi"] = _segment_directions(a.experiment, cols)
    g = a.gammas
    run = SimpleNamespace(summary={"one_particle": True}, scale=1.0, real_s=1.0, duration_s=1.0)
    return Analysis(setup or a.experiment, replace(g, columns={k: np.asarray(v)[:0] for k, v in g.columns.items()},
                                                   events=replace(a.events, columns=cols)), run=run)


def exact(exp, det: str) -> float:
    """⟨P⟩/B of the recoils in one ring around the beam: the cross sections of each reaction, at 24 depths, integrated
    over the ring's polar angles."""
    from physim.nuclear.rates import beam_energy_at, channels, reaction_at, stack

    mm = lambda q: float(str(q).split()[0])  # noqa: E731
    d = next(x for x in exp.detectors if x.name == det)
    lo, hi = (math.degrees(math.atan(mm(r) / mm(d.distance))) for r in (d.inner_radius, d.outer_radius))
    x, w = np.polynomial.legendre.leggauss(24)
    xa, wa = np.polynomial.legendre.leggauss(100)
    th = lo + (xa + 1) / 2 * (hi - lo)
    omega = wa * np.sin(np.radians(th))
    layers = stack(exp)
    excited = every = 0.0
    for xi, wi in zip(x, w):
        e = float(beam_energy_at(exp, 0, [(xi + 1) / 2 * layers[0].thickness], layers)[0])
        for ch in channels(layers, exp):
            r = reaction_at(exp, ch, e, layers)
            p = r.kinematics.at_lab(th, "recoil")[0]
            part = wi * float(np.sum(omega * np.asarray(r.cross_section_cm(180.0 - p.theta_cm)) * p.jacobian))
            every += part
            excited += part if ch.excitation is not None else 0.0
    return excited / every / exp.excitation.b_up_e2fm


def with_rings(**changes) -> Experiment:
    d = allmond_setup().to_dict()
    for x in d["detectors"]:
        x.update(changes)
    return Experiment.from_dict(d)


@pytest.fixture(scope="module")
def one_ring():
    exp = with_rings(rings=1)
    return Analysis(exp, simulate_gammas(exp, 40_000, seed=3))


@pytest.mark.parametrize("det, at_the_mean", [("Ring2", 0.006), ("Ring3", 0.012)])
def test_a_run_averages_the_excitation_over_each_segment(allmond, one_ring, det, at_the_mean):
    """Round 3, fix 2: within 0.1 % of the cross sections integrated over the ring, with seven rings and with one. At the
    mean angle of one ring ⟨P⟩/B was +0.6 and +1.2 % off."""
    s = every_particle(det)
    for a in (allmond[1][0], one_ring):
        run, truth = stored(a), exact(a.experiment, det)
        assert run._mean_probability(s) == pytest.approx(truth, rel=1e-3)
        assert run._mean_probability(s, mode="forward") == pytest.approx(truth * B_UP, rel=1e-3)
    grid, ratio = run._forward_ratio("recoil")
    rows = run._elastic_rows(s)
    old = np.average(np.interp(run.events.columns["theta"][rows], grid, ratio),
                     weights=run.events.columns["weight"][rows]) / B_UP
    assert old / truth - 1 == pytest.approx(at_the_mean, abs=1e-3)


def test_the_segment_averages_are_kept_and_the_systematics_use_them(one_ring, monkeypatch):
    run = stored(one_ring)
    s = every_particle("Ring3")
    p = run._mean_probability(s)
    assert one_ring.__dict__.get("_segment_kept") is None      # a weighted sample keeps its own angles
    monkeypatch.setattr(Analysis, "_lab_sums", lambda *a, **k: pytest.fail("worked out again"))
    assert run._mean_probability(s) == p
    assert run._mean_probability(every_particle("Ring2")) > p  # other rings of the same run
    monkeypatch.undo()
    up = run._mean_probability(s, energy_scale=1.005, exact=True) / run._mean_probability(s, exact=True) - 1
    assert up == pytest.approx(0.0514, rel=0.02)


def deviation(a, det: str, recoil: bool) -> tuple:
    """The excited share of the particles of one kind in a detector: the run's expectation over what the events did,
    less one, and the error of the latter."""
    run = stored(a)
    c = run.events.columns
    rows = (c["detector"] == run.events.detectors.index(det)) & (c["recoil"] == recoil)
    hit, w = rows & excited_rows(run), c["weight"]
    err = math.sqrt(float((w[hit] ** 2).sum())) / float(w[hit].sum())
    return float(np.average(run._excited_share(rows), weights=w[rows]) / (w[hit].sum() / w[rows].sum())) - 1, err, run


def test_a_segment_partly_below_the_threshold():
    """Carbon recoils at 50° and beyond, a threshold of 12 MeV: those from deep in the target, which are hardly ever
    excited, are not counted. (At the mean angle, all depths counted, the share was 75 % too small.)"""
    d = allmond_setup().to_dict()
    d["detectors"][1].update(inner_radius="60 mm", outer_radius="180 mm", threshold="12 MeV", rings=1)
    exp = Experiment.from_dict(d)
    dev, err, run = deviation(Analysis(exp, simulate_gammas(exp, 300_000, seed=8)), "Ring3", True)
    assert err < 0.06 and abs(dev) < 3 * err, (dev, err)
    nodes = run._segment_nodes()
    cut = run._cut_by_threshold("recoil", np.flatnonzero(nodes["recoil"]))
    names = np.array(nodes["detector"])[nodes["recoil"]]
    assert cut[names == "Ring3"].all() and not cut[names == "Ring2"].any()


def test_a_segment_in_which_the_scattered_beam_turns_back():
    """62Ni scatters off carbon by at most 11.2°, where its cross section per solid angle diverges: a ring across that
    angle, both kinds of particle in it."""
    d = allmond_setup().to_dict()
    d["detectors"] = [dict(d["detectors"][0], inner_radius="6 mm", outer_radius="14 mm", distance="60 mm", rings=1)]
    exp = Experiment.from_dict(d)
    a = Analysis(exp, simulate_gammas(exp, 300_000, seed=9))
    for recoil in (False, True):
        dev, err, run = deviation(a, "Ring2", recoil)
        assert err < (0.05 if recoil else 0.01) and abs(dev) < 3 * err, (recoil, dev, err)


def test_segments_of_tilted_and_off_axis_detectors():
    """A tilted target, a strip, a disc and an annulus away from the beam axis, one segment each or a few."""
    d = Experiment.example("coulex_ni58").to_dict()
    d["target"]["tilt"] = "30 deg"
    face = {"thickness": "300 um", "resolution": "50 keV", "threshold": "1 MeV"}
    d["detectors"] = [dict(face, name="Strip", shape="rectangle", theta="100 deg", phi="20 deg", distance="60 mm",
                           width="60 mm", height="40 mm", strips_x=1, strips_y=1),
                      dict(face, name="Side", shape="rectangle", theta="60 deg", phi="200 deg", distance="50 mm",
                           width="50 mm", height="50 mm", strips_x=4, strips_y=1),
                      dict(face, name="Disc", shape="circle", theta="140 deg", phi="100 deg", distance="50 mm",
                           radius="25 mm"),
                      dict(face, name="Ring", shape="annular", theta="150 deg", phi="270 deg", distance="70 mm",
                           inner_radius="8 mm", outer_radius="30 mm", rings=1, sectors=1)]
    exp = Experiment.from_dict(d)
    a = Analysis(exp, simulate_gammas(exp, 400_000, seed=10))
    for det, recoil in (("Strip", False), ("Side", False), ("Side", True), ("Disc", False), ("Ring", False)):
        dev, err, _ = deviation(a, det, recoil)
        assert err < 0.06 and abs(dev) < 3 * err + 0.01, (det, recoil, dev, err)


# -- the assumed inputs in the budget ----------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def reorientation():
    """The Allmond setup with a level scheme of the one state and Q_s = 0.0 ± 0.1 eb."""
    from physim.nuclear.levels import quadrupole_factor
    from physim.nuclear.orientation import simple_scheme

    def setup(q_eb: float, unc=None) -> Experiment:
        scheme = simple_scheme("62Ni", "1172.91 keV", "E2", f"{B_UP} e2fm4")
        f = quadrupole_factor(2.0)
        scheme.set_matrix_element(1, 1, "E2", 100.0 * q_eb / f, unc=None if unc is None else 100.0 * unc / f)
        d = allmond_setup().to_dict()
        d["levels"] = {"beam": scheme.to_dict()}
        return Experiment.from_dict(d)

    exp = setup(0.0, 0.1)
    return Analysis(exp, simulate_gammas(exp, 100_000, seed=6)), setup


def test_the_target_thickness_term(reorientation):
    """Round 3, fix 12: ±0.0138 mg/cm² changes ⟨P⟩/B by −2.05 and +2.10 %, as setups of those thicknesses."""
    a, _ = reorientation
    s = replace(every_particle(None), thickness_unc=0.0138 / 1.03741)
    more, less = a._thickness_changes(s)
    assert more == pytest.approx(-0.0205, rel=0.1) and less == pytest.approx(0.0210, rel=0.1)
    for change, t in ((more, 1.03741 + 0.0138), (less, 1.03741 - 0.0138)):
        d = allmond_setup().to_dict()
        d["target"]["thickness"] = f"{t!r} mg/cm2"
        own = Analysis(Experiment.from_dict(d), a.gammas)._mean_probability(s)
        assert own / Analysis(allmond_setup(), a.gammas)._mean_probability(s) - 1 == pytest.approx(change, rel=0.02)


@pytest.mark.parametrize("det, up, down", [("Ring2", 0.0271, -0.0266), ("Ring3", 0.0191, -0.0188)])
def test_the_term_of_the_assumed_matrix_elements(reorientation, det, up, down):
    """Q_s ± 0.1 eb, as the coupled equations of setups with that moment give it."""
    a, setup = reorientation
    s = every_particle(det)
    r = a.run(s)
    (name, plus, minus), = a._matrix_element_changes(s, r.fit)["changes"]
    assert "Q_s of the 2+ state" in name and "0.1 eb" in name
    assert plus == pytest.approx(up, abs=1e-3) and minus == pytest.approx(down, abs=1e-3)
    gain = r.budget["efficiency"] / 0.05
    assert r.budget["matrix elements assumed"] == pytest.approx(max(plus, -minus) * gain)
    assert any("Q_s of the 2+ state" in n for n in r.budget_notes)
    assert "Q_s of the 2+ state" in next(x["text"] for x in r.steps if x["step"] == "normalisation")
    base, higher = (Analysis(setup(q), a.gammas)._mean_probability(s, mode="forward") for q in (0.0, 0.1))
    assert higher / base - 1 == pytest.approx(plus, abs=3e-4)
    # For a run it is the same response, averaged over the segments.
    run = stored(a)
    (_, p, m), = run._matrix_element_changes(s, r.fit)["changes"]
    assert (p, m) == (pytest.approx(plus, abs=2e-4), pytest.approx(minus, abs=2e-4))


def test_matrix_elements_without_uncertainties_have_no_term(reorientation):
    a, setup = reorientation
    b = Analysis(setup(0.0), a.gammas)
    r = b.run(every_particle("Ring3"))
    assert "matrix elements assumed" not in r.budget
    assert any("no uncertainty given for any of the other matrix elements" in n and "Q_s of the 2+ state" in n
               for n in r.budget_notes)


# -- the beam-energy systematic -------------------------------------------------------------------------------------


@pytest.mark.parametrize("det, expected", [("Ring2", 0.0525), ("Ring3", 0.0514), (None, 0.0517)])
def test_the_beam_energy_systematic_is_a_change_of_the_entry_energy(allmond, det, expected):
    """Round 3, fix 5: the same as a setup whose beam has 0.5 % more energy. Scaled at every depth instead, the
    systematic was 17 % too small (4.3 %), since the energy lost changes little and what is left is less."""
    exp, samples = allmond
    a, s = samples[0], every_particle(det)
    d = exp.to_dict()
    d["beam"]["energy"] = f"{exp.beam.energy_mev * 1.005!r} MeV"
    up = Analysis(Experiment.from_dict(d), a.gammas)
    real = up._mean_probability(s, exact=True) / a._mean_probability(s, exact=True) - 1
    app = a._mean_probability(s, energy_scale=1.005, exact=True) / a._mean_probability(s, exact=True) - 1
    assert app == pytest.approx(real, rel=0.02) and real == pytest.approx(expected, rel=0.02)


def test_feeding_is_named_only_when_there_is_a_level_to_feed_from():
    """Round 3, fix 8: a scheme of one excited state (with its static moment, so to all orders) has no feeding."""
    from physim.nuclear.orientation import simple_scheme

    scheme = simple_scheme("62Ni", "1172.91 keV", "E2", f"{B_UP} e2fm4")
    scheme.set_matrix_element(1, 1, "E2", 0.0, source="user")
    d = allmond_setup().to_dict()
    d["levels"] = {"beam": scheme.to_dict()}
    exp = Experiment.from_dict(d)
    r = Analysis(exp, simulate_gammas(exp, 60_000, seed=4)).run(every_particle(None))
    text = next(x["text"] for x in r.steps if x["step"] == "all orders")
    assert "feeding" not in text and "levels above" not in text and "coupled equations" in text


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


# -- another line of the nucleus in the window -------------------------------------------------------------------------


#: The seeds of the clovers' closure test (fixed before they were run).
BLENDED_SEEDS = (31, 32, 33, 34)


@pytest.fixture(scope="module")
def blended():
    """The shipped 62Ni on 194Pt example with clovers (the 4⁺ → 2⁺ line of 62Ni, 2.5 % of the 2⁺ → 0⁺ one, is 9.4 keV below
    it, at 2σ of the peak the clovers make; 237 MeV, 0.5 mg/cm² of 194Pt, one slice of the target)."""
    exp = Experiment.example("coulex_ni62_pt194_clovers")
    return exp, [Analysis(exp, simulate_gammas(exp, 200_000, seed=seed)) for seed in BLENDED_SEEDS]


@pytest.mark.parametrize("det, scatter", [("CD", 0.0010), (None, 0.0024)])
def test_a_peak_with_another_line_in_its_window_keeps_its_area(blended, det, scatter):
    """The area is the full-energy counts of the 2⁺ → 0⁺ and the 4⁺ → 2⁺ lines, not 4 % less (a Gaussian on a line,
    as the peak was fitted whenever the window held another line, lost the tails that the add-back gives it: a γ ray put down
    to the crystal beside is corrected with that one's angle);
    ``scatter`` is of one sample. Lines more than 3σ away are in the fit, not in the area. Over sixteen samples (seeds
    1–16): +0.07 ± 0.02 % (CD) and +0.09 ± 0.06 % (all); before, −3.58 ± 0.05 and −4.00 ± 0.07 %."""
    exp, samples = blended
    s = every_particle(det)
    dev = []
    for a in samples:
        keep = a._gate(s)[0]
        g = a.gammas
        full = keep & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-6)
        r = a.run(replace(s, subtract_randoms=False))
        assert r.fit.tails is not None and r.fit.tails["lines"]
        true = full & a._own & (((g["initial"] == 1) & (g["final"] == 0)) | ((g["initial"] == 4) & (g["final"] == 1)))
        counted = float(g["weight"][true].sum()) * g.events.beam_time_s * g.live_fraction
        dev.append(r.area / counted - 1)
    err = scatter / math.sqrt(len(dev))
    assert abs(np.mean(dev)) < 0.001 + 3 * err, (det, dev)


def test_the_strip_detectors_keep_their_area_with_few_entries(blended):
    """About 430 full-energy entries a sample, too few to show the tails that the add-back gives the peak: the fit has
    the setup's shape (the target is one slice) however few they are. Of 96 samples the area is +0.11 ± 0.07 % off and
    scatters by 0.70 %; with the tails only where the counts showed them it was −1.87 ± 0.23 and 2.3 %."""
    dev = []
    for det in ("DSSD-L", "DSSD-R"):
        s = replace(every_particle(det), subtract_randoms=False)
        for a in blended[1]:
            g = a.gammas
            r = a.run(s)
            assert "spread" in r.fit.tails and not r.fit.tails["simulated"]
            true = (a._gate(s)[0] & a._own & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-6)
                    & (((g["initial"] == 1) & (g["final"] == 0)) | ((g["initial"] == 4) & (g["final"] == 1))))
            dev.append(r.area / (float(g["weight"][true].sum()) * g.events.beam_time_s * g.live_fraction) - 1)
    assert abs(np.mean(dev)) < 0.001 + 3 * 0.0064 / math.sqrt(len(dev)) and np.std(dev, ddof=1) < 0.014, dev


def test_the_monte_carlo_error_counts_particles_and_the_fit(blended):
    """γ rays of one particle are one entry: each given twice, with half the weight, they tell no more. (By γ ray the
    error fell by √2.) The error is the fit's for that many entries, so it is above that of the peak's counts."""
    exp, samples = blended
    a, s = samples[0], replace(every_particle("DSSD-L"), subtract_randoms=False)
    fit = a.run(s).fit
    once = a._monte_carlo_uncertainty(s, a._gate(s)[0], fit)
    g = a.gammas
    twice = {k: np.repeat(np.asarray(v), 2) for k, v in g.columns.items()}
    twice["weight"] = twice["weight"] / 2
    b = Analysis(exp, replace(g, columns=twice))
    assert b._monte_carlo_uncertainty(s, b._gate(s)[0], fit) == pytest.approx(once, rel=1e-9)
    inside = a._gate(s)[0] & (np.abs(g["corrected_projectile"] - fit.centroid) < 3 * fit.sigma)
    w, p = g["weight"][inside], g["particle"][inside]
    by_gamma = math.sqrt(float((w**2).sum())) / float(w.sum())
    by_particle = math.sqrt(sum(float(w[p == q].sum()) ** 2 for q in np.unique(p))) / float(w.sum())
    assert by_gamma < by_particle < once < 1.5 * by_particle


@pytest.mark.parametrize("det", ["CD", None])
def test_b_e2_comes_back_from_a_peak_with_another_line_in_its_window(blended, det):
    """B(E2) comes back. Over sixteen samples (seeds 1–16) it was −6.19 ± 0.33 % (CD) and −6.48 ± 0.40 % (all detectors)
    off, the area 4 % short, B(E2) moving 1.43 times as much as the yield, and the 4⁺ → 2⁺ line counted as if it were seen as
    often as the 2⁺ → 0⁺ one; now −0.14 ± 0.36 and +0.11 ± 0.50 %."""
    exp, samples = blended
    results = [a.run(every_particle(det)) for a in samples]
    b = np.array([r.b_e2fm4 / exp.excitation.b_up_e2fm - 1 for r in results])
    mean, unc = mean_and_error(b, [r.monte_carlo for r in results])
    assert unc < 0.009 and abs(mean) < 2.5 * unc, (det, mean, unc)
    assert all(any("4+ → 2+" in x for x in r.lines[0]["blends"]) for r in results)
