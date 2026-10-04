"""Count rates, beam time, expected peaks and the Monte Carlo event generator (backlog item 39)."""

import math

import numpy as np
import pytest

from physim import _core
from physim.nuclear import Experiment, data
from physim.nuclear.detectors import Array
from physim.nuclear.events import COLUMNS, generator_config, simulate
from physim.nuclear.kinematics import TwoBody
from physim.nuclear.rates import MB_CM2, Rates, stack
from physim.nuclear.rutherford import Rutherford
from physim.nuclear.stopping import Stopping


def tilted(name="oxygen_on_lead_array", tilt="30 deg"):
    d = Experiment.example(name).to_dict()
    d["target"]["tilt"] = tilt
    return Experiment.from_dict(d)


# -- the Rust pieces against the Python ones -----------------------------------------------------------------------


@pytest.mark.parametrize("beam, target, energy, ejectile, excitation", [
    ("16O", "208Pb", 64.0, None, 0.0),
    ("4He", "197Au", 5.5, None, 0.0),
    ("208Pb", "1H", 1000.0, None, 0.0),         # inverse kinematics, double-valued
    ("2H", "3H", 10.0, "n", 0.0),                # Q = +17.6 MeV
    ("1H", "12C", 30.0, None, 4.4389),            # inelastic
])
def test_rust_kinematics_matches_python(beam, target, energy, ejectile, excitation):
    r = TwoBody(beam, target, energy, ejectile=ejectile, excitation_mev=excitation)
    th = np.linspace(0.0, 180.0, 181)
    for recoil in (False, True):
        theta, e = _core.nuclear_two_body([r.m1, r.m2, r.m3, r.m4], energy, th, recoil)
        ref = r.at_cm(th, "recoil" if recoil else "ejectile")
        # Kinetic energy is E − m: both sides round at the level of ε m, and platforms differ there.
        m = max(r.m1, r.m2, r.m3, r.m4)
        np.testing.assert_allclose(e, ref.energy, rtol=1e-11, atol=10 * np.finfo(float).eps * m)
        # A particle left at rest (an elastic recoil at θ* = 180°) has no direction: compare angles elsewhere.
        moving = ref.energy > 1e-9 * m
        np.testing.assert_allclose(theta[moving], ref.theta_lab[moving], atol=1e-9)


@pytest.mark.parametrize("ion, material, energy, thickness", [
    ("4He", "Au", 5.5, 0.5), ("16O", "Pb", 64.0, 0.2), ("1H", "Si", 10.0, 50.0), ("16O", "Si", 40.0, 5.0),
    ("208Pb", "Si", 8.0, 0.2),
])
def test_transport_tables_match_the_stopping_module(ion, material, energy, thickness):
    """The generator's formula σ² = S(E_out)² (W(E_in) − W(E_out)) solves the same straggling equation that
    Stopping.straggling integrates step by step."""
    st = Stopping(ion, material)
    tab = st.transport_table()
    assert np.all(np.diff(tab["w"]) >= 0)
    out, sigma = _core.nuclear_cross(tab, energy, thickness)
    assert out == pytest.approx(st.energy_after(energy, thickness), rel=1e-9)
    assert sigma == pytest.approx(st.straggling(energy, thickness), rel=0.02)
    assert _core.nuclear_cross(tab, energy, 10 * st.range(energy)) == (0.0, 0.0)


# -- analytic rates ---------------------------------------------------------------------------------------------


def test_rate_of_a_small_detector_is_flux_times_cross_section_times_solid_angle():
    exp = Experiment.example("alpha_on_gold")
    r = Rates(exp)
    layer = stack(exp)[0]
    n = layer.nuclides[(79, 197)]
    pps = exp.beam.particles_per_second
    assert pps == pytest.approx(1e-9 / 1.602176634e-19)
    for name in ("A30", "A60", "A135"):
        g = r.array[name]
        e_mid = Stopping("4He", "Au").energy_after(5.5, layer.thickness / 2)
        ruth = Rutherford("4He", "197Au", e_mid)
        sigma, _ = ruth.cross_section_lab(g.mean_theta())
        expect = pps * n * sigma * g.solid_angle() * 1e-3 * MB_CM2
        # The 2.5 mm discs span ±1.4°: the curvature of 1/sin⁴ adds about 0.5% at 30°.
        assert r.by_channel(name)[("197Au (target)", "ejectile")] == pytest.approx(expect, rel=0.01)


def test_rates_add_up_and_set_the_beam_time():
    exp = Experiment.example("oxygen_on_lead_array")
    r = Rates(exp)
    seg = r.per_segment("CD", counted=False)
    assert len(seg) == 16 * 24
    assert sum(seg.values()) == pytest.approx(r.rate("CD", counted=False), rel=1e-12)
    # The annular detector is symmetric about the beam: every sector of a ring sees the same rate.
    ring = [seg[(5, j)] for j in range(24)]
    assert np.ptp(ring) < 1e-6 * np.mean(ring)
    assert sum(r.by_channel("DSSD1").values()) == pytest.approx(r.rate("DSSD1", counted=False))
    t = r.beam_time_s
    assert t == 12 * 3600
    assert r.counts_in_run("DSSD3") == pytest.approx(r.rate("DSSD3") * t)
    assert r.beam_time_for("DSSD3") == pytest.approx(5000 / r.rate("DSSD3"))
    assert r.beam_time_for("DSSD3", counts=100) == pytest.approx(100 / r.rate("DSSD3"))
    assert r.relative_error("DSSD3") == pytest.approx(1 / math.sqrt(r.counts_in_run("DSSD3")))
    with pytest.raises(KeyError):
        r.rate("nope")


def test_warnings():
    w = Rates(Experiment.example("alpha_on_gold")).warnings()
    assert any("A20 counts" in x and "pile-up" in x for x in w)
    assert not any("A45 counts" in x for x in w)
    w = Rates(Experiment.example("oxygen_on_lead_array")).warnings(max_rate=1e9)
    assert any(x.startswith("12C (backing)") and "Coulomb barrier" in x for x in w)  # 16O on carbon: above it
    assert not any(x.startswith("208Pb") for x in w)  # 64 MeV 16O on lead: below the barrier
    d = Experiment.example("alpha_on_gold").to_dict()
    d["run"]["counts_wanted"] = 10**9
    w = Rates(Experiment.from_dict(d)).warnings()
    assert any("A135 collects" in x and "need" in x for x in w)


def test_distant_collisions_are_cut_by_energy():
    """Rutherford's cross section diverges for distant collisions, which send almost-90° recoils out with almost no
    energy: they are left out below the lowest threshold, so rates stay finite."""
    r = Rates(Experiment.example("alpha_on_gold"))
    assert r.min_energy == pytest.approx(0.2)
    for name, rate in r.per_detector(counted=False).items():
        assert 0 < rate < 1e5, name


# -- the event generator against the analytic rates and peaks ---------------------------------------------------


@pytest.mark.parametrize("make", [
    lambda: Experiment.example("alpha_on_gold"),
    lambda: Experiment.example("oxygen_on_lead_array"),
    lambda: tilted(),
], ids=["alpha_on_gold", "oxygen_on_lead_array", "tilted_30deg"])
def test_monte_carlo_rates_agree_with_analytic_rates(make):
    exp = make()
    r = Rates(exp)
    ev = simulate(exp, 2_000_000, seed=3)
    for name in ev.detectors:
        for counted in (False, True):
            a = r.rate(name, counted=counted)
            m, err = ev.rate(name, counted=counted)
            # Analytic "counted" rates judge the threshold from mean energies, so allow 2% on top of statistics.
            slack = 0.0 if not counted else 0.02 * a
            assert abs(m - a) < 4 * err + slack, (name, counted, a, m, err)
        for (channel, particle), a in r.by_channel(name).items():
            m, err = ev.rate(name, counted=False, channel=channel, particle=particle)
            assert abs(m - a) < 4 * err + 1e-3 * a, (name, channel, particle, a, m, err)


def test_peak_positions_and_widths():
    """Mean measured energy against kinematics plus mean energy loss; width against resolution, straggling,
    kinematic broadening and target thickness added in quadrature."""
    exp = Experiment.example("alpha_on_gold")
    r = Rates(exp)
    ev = simulate(exp, 4_000_000, seed=9)
    for name in ("A20", "A30", "A45", "A60"):
        pk = r.peaks(name)[0]
        assert (pk.channel, pk.particle) == ("197Au (target)", "ejectile")
        m = ev.select(name, particle="ejectile")
        x, w = ev["measured"][m], ev["weight"][m]
        mean = np.average(x, weights=w)
        sd = math.sqrt(np.average((x - mean) ** 2, weights=w))
        n = m.sum()
        assert abs(mean - pk.mean) < max(4 * sd / math.sqrt(n), 1e-3), name
        assert sd == pytest.approx(pk.sigma, rel=max(0.05, 4 / math.sqrt(2 * n))), name
        assert pk.components["resolution"] == pytest.approx(0.020 / 2.35482, rel=1e-4)
    # By hand: slow the beam to the middle of the foil, scatter to 30°, leave through the other half at 1/cos 30°.
    pk = r.peaks("A30")[0]
    st = Stopping("4He", "Au")
    k = TwoBody("4He", "197Au", st.energy_after(5.5, 0.25)).at_lab(30.0)[0].energy
    assert pk.mean == pytest.approx(st.energy_after(k, 0.25 / math.cos(math.radians(30))), abs=1e-3)


def test_peaks_from_a_backing_and_recoils():
    r = Rates(Experiment.example("oxygen_on_lead_array"))
    labels = {(pk.channel, pk.particle) for pk in r.peaks("DSSD1")}
    assert {("208Pb (target)", "ejectile"), ("208Pb (target)", "recoil"), ("12C (backing)", "ejectile"),
            ("12C (backing)", "recoil")} <= labels
    lead = r.peaks("DSSD1")[0]
    assert lead.channel == "208Pb (target)" and lead.particle == "ejectile"
    # A 50 mm strip detector at 150 mm spans ±9°: kinematic broadening dominates the whole-detector peak, but not
    # a single strip's.
    strip = [pk for pk in r.peaks("DSSD1", (8, 8)) if pk.particle == "ejectile" and pk.channel == lead.channel][0]
    assert lead.components["kinematic"] > 10 * strip.components["kinematic"]
    assert strip.sigma < lead.sigma / 5


# -- reproducibility and output ---------------------------------------------------------------------------------


def test_same_seed_same_events_and_pieces_add_up():
    exp = Experiment.example("oxygen_on_lead_array")
    a = simulate(exp, 50_000, seed=4)
    b = simulate(exp, 50_000, seed=4)
    c = simulate(exp, 50_000, seed=5)
    assert set(a.columns) == set(COLUMNS)
    for k in a.columns:
        assert np.array_equal(a[k], b[k]), k
    assert not np.array_equal(a["measured"], c["measured"])
    config, _ = generator_config(exp)
    first = _core.nuclear_events(config, 20_000, 4, total=50_000)
    rest = _core.nuclear_events(config, 30_000, 4, total=50_000, first=20_000)
    for k in a.columns:
        assert np.array_equal(np.concatenate([first[k], rest[k]]), a[k]), k


def test_columns_and_spectra():
    exp = Experiment.example("oxygen_on_lead_array")
    ev = simulate(exp, 200_000, seed=2)
    c = ev.columns
    assert np.all(np.diff(c["event"]) >= 0)
    assert np.all((c["detector"] >= 0) & (c["detector"] < len(ev.detectors)))
    cd = c["detector"] == ev.detectors.index("CD")
    assert c["segment_i"][cd].max() < 16 and c["segment_j"][cd].max() < 24
    assert np.all(c["depth"] <= sum(lay.thickness for lay in stack(exp)) + 1e-12)
    assert np.all(c["energy_face"] <= c["energy"] + 1e-12)
    assert np.all(c["deposited"] <= c["energy_face"] + 1e-12)
    assert np.all(c["measured"][c["counted"]] >= 0.3)
    # The backward annular detector sees only ejectiles, which went backwards.
    assert np.all(c["theta"][cd] > 90) and not np.any(c["recoil"][cd])
    h, edges = ev.spectrum("DSSD1", bins=100)
    assert h.sum() == pytest.approx(ev.counts("DSSD1"), rel=1e-9)
    assert len(edges) == 101
    assert ev.counts("DSSD1") == pytest.approx(ev.rate("DSSD1")[0] * 12 * 3600)
    with pytest.raises(ValueError):
        ev.select(particle="beam")
    with pytest.raises(ValueError):
        simulate(exp, 0)


def test_nothing_to_simulate_without_a_visible_detector():
    d = Experiment.example("alpha_on_gold").to_dict()
    d["detectors"] = [dict(d["detectors"][0], theta="0.2 deg", radius="0.1 mm", threshold="10 MeV")]
    with pytest.raises(ValueError, match="nothing to simulate"):
        generator_config(Experiment.from_dict(d))


def test_pictures():
    import matplotlib

    matplotlib.use("Agg")
    from physim.nuclear import plot

    ev = simulate(Experiment.example("oxygen_on_lead_array"), 100_000, seed=1)
    fig = plot.spectra(ev, ["DSSD1", "CD"])
    assert len(fig.axes) == 2
    ax = plot.theta_energy(ev)
    assert ax.get_xlabel().startswith("lab angle")


def test_species_and_masses_in_the_config():
    exp = Experiment.example("oxygen_on_lead_array")
    config, used = generator_config(exp)
    assert config["species"] == ["16O", "208Pb", "12C", "13C"]
    assert config["masses"][0] == pytest.approx(data.nuclide("16O").nuclear_mass_mev)
    assert len(config["tables"]) == 4 * config["materials"]
    assert [ch.label for ch in used] == ["208Pb (target)", "12C (backing)", "13C (backing)"]
    assert len(config["faces"]) == len(Array.from_experiment(exp).geometries)
