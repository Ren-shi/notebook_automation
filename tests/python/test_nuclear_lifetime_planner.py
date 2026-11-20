"""The planner's lifetime views (lifetime wave 3A): Planner.lifetimes, Planner.lifetime_spectra and
Planner.analyse_lifetimes, to the frozen shapes of planning/lifetimes.md.

The setups are the tests' own: the examples' detectors and beams with level schemes typed in here (no local ENSDF
copy is needed, as on CI) and lifetimes chosen here, so nothing pins the three lifetime examples' numbers. The
analysis is checked on synthetic runs drawn from a known τ, so it runs in seconds."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from physim.nuclear import Experiment, lifetimes
from physim.nuclear import lifetime_plan as lp
from physim.nuclear import transport
from physim.nuclear.gamma import doppler_energy
from physim.nuclear.gamma_events import crystals_of
from physim.nuclear.levels import LN2
from physim.nuclear.planner import Planner
from physim.nuclear.runs import DistanceSeries

PLAN_KEYS = {"available", "reason", "reaction", "nucleus", "method", "distances_um", "precision", "levels",
             "time_to_stop_ps", "problems", "warnings"}
LEVEL_KEYS = {"level", "energy_kev", "gamma_kev", "tau_ps", "tau_source", "beta", "method",
              "recommended_distances_um", "sensitive_um", "separation", "attenuation", "hours_per_distance"}
SPECTRA_KEYS = {"available", "reason", "method", "edges_kev", "spectra"}
ANALYSIS_KEYS = {"available", "reason", "level", "method", "tau_ps", "stat_ps", "budget", "points", "curve",
                 "b_down", "b_unit", "compare", "notes"}
POINT_KEYS = {"distance_um", "unshifted", "shifted", "error", "tau_ps", "tau_err_ps"}

#: The plunger distances of the Coulex setup (µm): log-spaced over 3–500 µm.
DISTANCES = (3.0, 8.0, 20.0, 45.0, 100.0, 220.0, 500.0)


def _coulex(b_up="2500 e2fm4", **lifetime):
    """The 98Ru-on-24Mg detectors and beam with this test's own B(E2), plunger distances and stopper; no
    deorientation unless given."""
    d = Experiment.example("coulex_ru98_plunger").to_dict()
    d["reaction"]["b_up"] = b_up
    d["lifetime"] = {"method": "rdds", "distances": [f"{x:g} um" for x in DISTANCES],
                     "stopper": {"kind": "stopper", "material": "Cu", "thickness": "20 mg/cm2"}, **lifetime}
    return Experiment.from_dict(d)


#: A yrast band typed in for 50Cr: (J^π, level keV, T½ ps, γ-ray keV to the level below). The half-lives are this
#: test's (one per method range), not ENSDF's.
CR50_BAND = (("0+", 0.0, None, None), ("2+", 783.3, 7.0, 783.3), ("4+", 1881.4, 2.0, 1098.1),
             ("6+", 3164.1, 0.9, 1282.7), ("8+", 4744.9, 0.2, 1580.8))


def _band(nuclide, rows, multipolarity="E2", source="user"):
    levels, transitions = [], []
    for k, (jpi, e, t_half, e_gamma) in enumerate(rows):
        lv = {"energy": f"{e:g} keV", "jpi": jpi, "source": source}
        if t_half is not None:
            lv["half_life"] = f"{t_half:g} ps"
        levels.append(lv)
        if e_gamma is not None:
            transitions.append({"from": k, "to": k - 1, "energy": f"{e_gamma:g} keV", "intensity": 100,
                                "multipolarity": multipolarity, "source": source})
    return {"nuclide": nuclide, "level": levels, "transition": transitions}


def _cr50():
    d = Experiment.example("fusion_cr50_plunger").to_dict()
    d.setdefault("levels", {})["50Cr"] = _band("50Cr", CR50_BAND)
    return Experiment.from_dict(d)


def _sm143(taus_ps=(0.6, 0.3, 0.2)):
    """The 143Sm-on-gold DSAM detectors and beam with a small dipole band typed in (lifetimes this test's)."""
    d = Experiment.example("fusion_sm143_dsam").to_dict()
    d["reaction"]["side_feeding_levels"] = 0
    for ch in d["reaction"]["channels"]:
        if ch["residue"] == "143Sm":
            ch["entry_spin"] = 19.5
    rows = [("13/2-", 0.0, None, None)]
    for k, spin in enumerate((15, 17, 19), start=1):
        rows.append((f"{spin}/2-", 500.0 * k + 30.0 * k * k, taus_ps[k - 1] * LN2, 500.0 + 30.0 * (2 * k - 1)))
    d.setdefault("levels", {})["143Sm"] = _band("143Sm", rows, multipolarity="M1")
    return Experiment.from_dict(d)


@pytest.fixture(scope="module")
def coulex():
    return Planner(_coulex())


@pytest.fixture(scope="module")
def cr50():
    return Planner(_cr50())


@pytest.fixture(scope="module")
def sm143():
    return Planner(_sm143())


# ---------------------------------------------------------------------------------------------------------------
# Planner.lifetimes


def test_methods_follow_the_lifetime_ranges():
    assert lp.method_for(0.005) == "outside"
    assert lp.method_for(0.05) == "dsam"
    assert lp.method_for(0.99) == "dsam"
    assert lp.method_for(1.0) == "rdds"
    assert lp.method_for(500.0) == "rdds"
    assert lp.method_for(5000.0) == "outside"
    assert lp.method_for(None) == "outside"
    # The one-level precision is best near x = 1.59 and worse on either side.
    assert lp.relative_precision(1.59) < lp.relative_precision(0.3) and lp.relative_precision(1.59) < \
        lp.relative_precision(5.0)


def test_coulex_plan_has_the_frozen_shape_and_the_setups_numbers(coulex):
    out = coulex.lifetimes()
    assert PLAN_KEYS <= set(out) and out["available"] and out["reaction"] == "coulex"
    assert out["nucleus"] == "98Ru" and out["method"] == "rdds" and out["precision"] == 0.1
    assert out["distances_um"] == list(DISTANCES)
    row = out["levels"][0]
    assert LEVEL_KEYS <= set(row)
    exp = coulex.experiment
    # τ from the setup's B(E2): B(E2↓) = B(E2↑)/5 for 0+ → 2+.
    tau = lifetimes.tau_from_b(1e3 * exp.excitation.energy_mev, "E2", exp.excitation.b_up_e2fm / 5)
    assert row["tau_ps"] == pytest.approx(tau, rel=1e-6) and "B(E2)" in row["tau_source"]
    assert row["method"] == "rdds"
    assert row["recommended_distances_um"] == lifetimes.recommend_distances(row["tau_ps"], row["beta"])
    lo, hi = lifetimes.sensitive_range(row["tau_ps"], row["beta"])
    assert row["sensitive_um"] == pytest.approx([lo, hi])
    assert row["attenuation"] is None and out["time_to_stop_ps"] > 0
    # Forward detectors see the line shifted up, backward ones down; 90° hardly at all.
    by_theta = {s["theta_deg"]: s for s in row["separation"]}
    assert by_theta[41.5]["separation_kev"] > 0 > by_theta[138.5]["separation_kev"]
    assert by_theta[41.5]["resolved"] and not by_theta[90.0]["resolved"]
    assert all(s["fwhm_kev"] > 0 for s in row["separation"])
    assert row["hours_per_distance"] > 0 and row["rate_per_s"] > 0


def test_beam_time_per_distance_goes_as_the_precision_squared(coulex):
    a = coulex.lifetimes(precision=0.1)["levels"][0]["hours_per_distance"]
    b = coulex.lifetimes(precision=0.05)["levels"][0]["hours_per_distance"]
    assert b == pytest.approx(4 * a, rel=1e-9)
    with pytest.raises(ValueError):
        coulex.lifetimes(precision=1.5)


def test_a_coulex_plunger_without_deorientation_is_warned_about_above_one_ps():
    long_lived = Planner(_coulex()).lifetimes()
    assert long_lived["levels"][0]["tau_ps"] > 1
    assert any("deorientation" in w for w in long_lived["warnings"])
    assert any("No deorientation" in w.text for w in Planner(_coulex()).warnings())  # the banner too
    given = Planner(_coulex(deorientation={"alpha2": 0.2, "gamma2": "0.1 /ps", "alpha4": 0.2, "gamma4": "0.3 /ps"}))
    assert not any("No deorientation" in w for w in given.lifetimes()["warnings"])
    assert not any("No deorientation" in w.text for w in given.warnings())
    short = Planner(_coulex(b_up="60000 e2fm4")).lifetimes()  # τ below a picosecond
    assert short["levels"][0]["tau_ps"] < 1
    assert not any("No deorientation" in w for w in short["warnings"])
    # ... and that level suits DSAM, which the plan says.
    assert short["levels"][0]["method"] == "dsam" and any("DSAM" in w for w in short["warnings"])


def test_fusion_plan_lists_the_band_studied_level_first(cr50):
    out = cr50.lifetimes()
    assert out["available"] and out["nucleus"] == "50Cr" and out["reaction"] == "fusion-evaporation"
    rows = {r["level"]: r for r in out["levels"]}
    assert out["levels"][0]["level"].startswith("2+")
    taus = {jpi: t * 1e0 / LN2 for jpi, _e, t, _g in CR50_BAND if t is not None}
    for jpi, tau in taus.items():
        row = next(r for name, r in rows.items() if name.startswith(jpi))
        assert row["tau_ps"] == pytest.approx(tau, rel=1e-6)
        assert row["tau_source"] == "half-life (user)"
        assert row["method"] == lp.method_for(tau)
    eight = next(r for name, r in rows.items() if name.startswith("8+"))
    assert eight["method"] == "dsam" and eight["recommended_distances_um"] == [] and eight["sensitive_um"] is None
    # Stronger lines feed lower levels in a yrast cascade: the 2+ line is the strongest.
    assert out["levels"][0]["rate_per_s"] >= max(r["rate_per_s"] for r in out["levels"][1:])
    one = cr50.lifetimes(level="4+")
    assert [r["level"] for r in one["levels"]] == [next(n for n in rows if n.startswith("4+"))]
    assert cr50.lifetimes(level="1098 keV")["levels"][0]["gamma_kev"] == pytest.approx(1098.1)
    with pytest.raises(ValueError):
        cr50.lifetimes(level="12+")


def test_dsam_plan_gives_f_on_the_backing_falling_with_tau(sm143):
    out = sm143.lifetimes()
    assert out["available"] and out["method"] == "dsam" and out["time_to_stop_ps"] > 0
    rows = out["levels"]
    assert all(r["method"] == "dsam" for r in rows)
    pairs = sorted((r["tau_ps"], r["attenuation"]) for r in rows)
    assert all(0 < f < 1 for _t, f in pairs)
    assert [f for _t, f in pairs] == sorted((f for _t, f in pairs), reverse=True)
    assert all(r["hours_per_distance"] is None for r in rows)
    assert any("known failures" in w for w in out["warnings"])


def test_a_setup_without_lifetimes_says_why():
    out = Planner.example("alpha_on_gold").lifetimes()
    assert PLAN_KEYS <= set(out) and not out["available"] and "Coulomb excitation" in out["reason"]
    assert out["levels"] == []


def test_a_fusion_setup_without_a_level_scheme_is_not_available(monkeypatch):
    from physim.nuclear import ensdf

    monkeypatch.setattr(ensdf, "available", lambda: False)
    monkeypatch.setattr("physim.nuclear.residues.residue_scheme", lambda exp, residue, *a, **k: None)
    out = Planner.example("fusion_cr50_plunger").lifetimes()
    assert not out["available"] and "level scheme" in out["reason"] and out["levels"] == []


# ---------------------------------------------------------------------------------------------------------------
# Planner.lifetime_spectra


def test_coulex_plunger_spectra_per_distance_and_ring(coulex):
    out = coulex.lifetime_spectra()
    assert SPECTRA_KEYS <= set(out) and out["available"] and out["method"] == "rdds"
    rings = sorted({s["theta_deg"] for s in out["spectra"]})
    assert rings == [41.5, 90.0, 138.5]
    assert len(out["spectra"]) == len(DISTANCES) * len(rings)
    edges = np.asarray(out["edges_kev"])
    e0 = out["gamma_kev"]
    for s in out["spectra"]:
        assert len(s["counts"]) == len(edges) - 1 and min(s["counts"]) >= 0
        assert s["unshifted_kev"] == pytest.approx(e0)
    fwd = [s for s in out["spectra"] if s["theta_deg"] == 41.5]
    shares = [s["unshifted_share"] for s in fwd]
    assert shares == sorted(shares, reverse=True) and shares[0] > 0.5 > shares[-1]
    assert fwd[0]["shifted_kev"] > e0
    # The line's counts in each ring: its rate times the ring's share of the efficiency times the time per distance.
    plan = coulex.lifetimes()["levels"][0]
    per = out["seconds_per_distance"]
    total = sum(sum(s["counts"]) for s in out["spectra"] if s["distance_um"] == DISTANCES[0])
    assert total == pytest.approx(plan["rate_per_s"] * per, rel=0.02)
    one = coulex.lifetime_spectra(distance_um=20.0)
    assert {s["distance_um"] for s in one["spectra"]} == {20.0}


def test_fusion_spectra_come_from_the_residues(cr50, sm143):
    out = cr50.lifetime_spectra(distance_um=42.0)
    assert out["available"] and "plunger_spectra" in out["model"]
    assert sorted({s["theta_deg"] for s in out["spectra"]}) == [45.0, 142.3]
    assert all(s["distance_um"] == 42.0 for s in out["spectra"])
    assert all(sum(s["counts"]) > 0 for s in out["spectra"])
    d = sm143.lifetime_spectra()
    assert d["available"] and d["method"] == "dsam" and "dsam_lineshapes" in d["model"]
    assert all(s["distance_um"] is None for s in d["spectra"])
    fwd = next(s for s in d["spectra"] if s["theta_deg"] == 57.0)
    back = next(s for s in d["spectra"] if s["theta_deg"] == 148.0)
    assert fwd["shifted_kev"] > d["gamma_kev"] > back["shifted_kev"]
    with pytest.raises(ValueError):
        sm143.lifetime_spectra(distance_um=10.0)


def test_spectra_without_a_lifetime_section_say_so():
    out = Planner.example("coulex_ni58").lifetime_spectra()
    assert SPECTRA_KEYS <= set(out) and not out["available"] and "set_lifetime" in out["reason"]


# ---------------------------------------------------------------------------------------------------------------
# Planner.analyse_lifetimes on synthetic runs


class _Gammas:
    """A particle run's γ rays as the lifetime analysis reads them (GammaEvents' surface)."""

    live_fraction = 1.0

    def __init__(self, measured_mev, crystal, crystals):
        self.columns = {"measured": measured_mev, "crystal": crystal, "weight": np.ones(len(measured_mev)),
                        "corrected_projectile": measured_mev, "corrected_recoil": measured_mev}
        self.crystals = crystals
        self.events = SimpleNamespace(beam_time_s=1.0)

    def select(self):
        return np.ones(len(self.columns["measured"]), dtype=bool)

    def __getitem__(self, key):
        return self.columns[key]


class _Run:
    def __init__(self, experiment, gammas, distance, kind="beam"):
        self.experiment, self._g, self.distance_um, self.kind = experiment, gammas, distance, kind
        self.real_s, self.duration_s, self.label = 10.0, 100.0, f"run-{distance}"

    def gammas(self, plain=False):
        return self._g


def _plunger_series(exp, tau_ps, beta, n=6000, seed=3, three_foil=None):
    """One synthetic run per distance: each crystal's γ rays shifted (decays before the stopper, at β) or not, with
    the share of a single level of mean life ``tau_ps``; ``three_foil`` = (β′, gap µm) adds the reduced-shifted
    component of a differential plunger."""
    rng = np.random.default_rng(seed)
    crystals = crystals_of(exp, plain=True)[0]
    theta = np.radians([float(exp.gamma_detectors[c.detector].theta.to("deg")) for c in crystals])
    e0 = 1e3 * exp.excitation.energy_mev
    sigma = 1.0
    lt_ = exp.lifetime
    runs = {}
    for d in lt_.distances_um:
        t1 = (d + lt_.offset_um) / (beta * lifetimes.C_UM_PER_PS)
        crystal = rng.integers(0, len(crystals), n)
        u1 = math.exp(-t1 / tau_ps)
        x = rng.random(n)
        e = np.full(n, e0)
        fly = x >= u1
        e[fly] = doppler_energy(e0, beta, np.cos(theta[crystal[fly]]))
        if three_foil is not None:
            b2, gap = three_foil
            u2 = math.exp(-(t1 + gap / (b2 * lifetimes.C_UM_PER_PS)) / tau_ps)
            mid = (x < u1) & (x >= u2)
            e[mid] = doppler_energy(e0, b2, np.cos(theta[crystal[mid]]))
        e = e + rng.normal(0.0, sigma, n)
        runs[float(d)] = [_Run(exp, _Gammas(1e-3 * e, crystal, crystals), d)]
    return DistanceSeries(exp.to_dict(), runs)


@pytest.mark.parametrize("method", [None, "decay_curve"])
def test_tau_comes_back_from_a_synthetic_distance_series(method):
    p = Planner(_coulex())
    tau = p.lifetimes()["levels"][0]["tau_ps"]
    beta = 0.045
    series = _plunger_series(p.experiment, tau, beta)
    out = p.analyse_lifetimes(series=series, method=method)
    assert ANALYSIS_KEYS <= set(out) and out["available"], out["reason"]
    assert out["method"] == (method or "ddcm") and out["level"] == "2"
    assert out["beta"] == pytest.approx(beta, rel=0.03)
    assert abs(out["tau_ps"] - tau) < max(3 * out["stat_ps"], 0.08 * tau)
    assert out["budget"]["statistics"] == pytest.approx(out["stat_ps"])
    assert len(out["points"]) == len(DISTANCES) and POINT_KEYS <= set(out["points"][0])
    shares = [pt["unshifted"] for pt in out["points"]]
    assert shares == sorted(shares, reverse=True)
    assert set(out["curve"]) >= {"x", "y", "kind"} and len(out["curve"]["x"]) == len(out["curve"]["y"]) > 0
    b, unc = out["b_down"]
    assert out["b_unit"] == "e2fm4" and b == pytest.approx(lifetimes.b_from_tau(
        1e3 * p.experiment.excitation.energy_mev, "E2", out["tau_ps"]), rel=1e-6) and unc > 0
    sources = {c["source"] for c in out["compare"]}
    assert "setup's Coulex B(E2)" in sources
    coulex_b = next(c for c in out["compare"] if c["source"] == "setup's Coulex B(E2)")
    assert coulex_b["tau_ps"] == pytest.approx(tau, rel=1e-6)
    assert coulex_b["b_down"] == pytest.approx(p.experiment.excitation.b_up_e2fm / 5)


def test_a_differential_plunger_is_analysed_by_the_three_foil_ddcm():
    exp = _coulex(stopper={"kind": "differential", "material": "Mg", "thickness": "1 mg/cm2"},
                  degrader2={"material": "Au", "thickness": "30 mg/cm2", "gap": "150 um"})
    p = Planner(exp)
    tau = p.lifetimes()["levels"][0]["tau_ps"]
    beta = 0.045
    b2 = lp.degraded_beta(lp.Subject(p)) * beta / lp.Subject(p).beta
    series = _plunger_series(exp, tau, beta, n=8000, three_foil=(b2, 150.0))
    out = p.analyse_lifetimes(series=series)
    assert out["available"], out["reason"]
    assert out["method"] == "ddcm_three_foil"
    assert abs(out["tau_ps"] - tau) < max(3 * out["stat_ps"], 0.1 * tau)
    pt = out["points"][0]
    assert {"fully_shifted", "reduced_shifted"} <= set(pt)
    assert pt["unshifted"] + pt["shifted"] == pytest.approx(1.0)


def test_analysis_without_runs_or_with_the_wrong_method_says_why(coulex, sm143):
    out = coulex.analyse_lifetimes()
    assert ANALYSIS_KEYS <= set(out) and not out["available"] and "No plunger runs" in out["reason"]
    assert coulex.analyse_lifetimes(method="lineshape")["reason"].startswith("A line-shape fit is for DSAM")
    assert not sm143.analyse_lifetimes()["available"]
    with pytest.raises(ValueError):
        coulex.analyse_lifetimes(method="napatau")


def test_dsam_tau_comes_back_from_a_synthetic_line_shape(sm143, monkeypatch):
    p = Planner(_sm143())
    spectra = p.lifetime_spectra()
    edges = np.asarray(spectra["edges_kev"])
    exp = p.experiment
    crystals = crystals_of(exp, plain=True)[0]
    ring_crystals = {}
    for k, c in enumerate(crystals):
        ring_crystals.setdefault(round(float(exp.gamma_detectors[c.detector].theta.to("deg")), 1), []).append(k)
    rng = np.random.default_rng(5)
    energies, crystal = [], []
    for s in spectra["spectra"]:
        expected = np.asarray(s["counts"], dtype=float)
        expected *= 40_000 / expected.sum()  # a run's worth of counts in each ring
        counts = rng.poisson(expected)
        e = np.repeat(edges[:-1], counts) + rng.random(int(counts.sum())) * np.repeat(np.diff(edges), counts)
        energies.append(e)
        crystal.append(rng.choice(ring_crystals[s["theta_deg"]], len(e)))
    energies, crystal = np.concatenate(energies), np.concatenate(crystal)

    class FusionGammas:
        columns = {"measured": 1e-3 * energies, "crystal": crystal, "corrected": 1e-3 * energies}

        def __init__(self):
            self.crystals = [c.name for c in crystals]

        def __len__(self):
            return len(energies)

        def __getitem__(self, key):
            return self.columns[key]

    run = _Run(exp, FusionGammas(), None, kind="fusion")
    monkeypatch.setattr(lp, "_dsam_runs", lambda planner, experiment: [run])
    out = p.analyse_lifetimes()
    assert out["available"], out["reason"]
    assert out["method"] == "lineshape" and out["points"] == []
    tau = p.lifetimes()["levels"][0]["tau_ps"]
    assert out["level"] == p.lifetimes()["levels"][0]["level"]
    assert abs(out["tau_ps"] - tau) < max(3 * out["total_ps"], 0.25 * tau)
    assert out["budget"]["stopping"] > 0
    assert out["curve"]["kind"].startswith("line shape") and len(out["curve"]["x"]) == len(out["curve"]["y"])
    assert out["b_unit"] == "mu_N2"


def test_planning_kernel_stops_the_recoil(sm143):
    s = lp.Subject(sm143)
    k, t_stop = lp.planning_kernel(s, n=200, times=40)
    assert t_stop > 0 and k.times_ps[-1] >= t_stop
    assert np.all(k.beta[-1] == 0)
    assert isinstance(k, transport.Kernel)


# ---------------------------------------------------------------------------------------------------------------
# 98Ru closure through the planner (lifetime wave 4B)

#: Distances of the closure series (µm): the example's range, fewer of them.
RU98_DISTANCES = (3.0, 10.0, 25.0, 60.0, 160.0, 450.0)


class _SimGammas(_Gammas):
    """A simulated distance's counted γ rays, their weights times ``scale`` (the same at every distance, as the
    reactions are: the efficiency ratio is fitted from a constant total)."""

    def __init__(self, g, scale):
        m = np.asarray(g["counted"], dtype=bool)
        super().__init__(np.asarray(g["measured"])[m], np.asarray(g["crystal"])[m], g.crystals)
        self.columns["weight"] = np.asarray(g["weight"], dtype=float)[m] * scale


def _ru98_series(exp, particle_events, seed):
    """One simulated run per distance from the same particle events (the beam and the recoils that reach the
    silicon detector, as a planner run has them), γ rays with their own seed per distance, in effective counts
    (Σw/Σw² of the first distance's γ rays, so that the spectra's errors are about those of counts)."""
    from physim.nuclear import gamma_events

    runs, scale = {}, None
    for k, d in enumerate(RU98_DISTANCES):
        g = gamma_events.simulate_gammas(exp, particle_events=particle_events, seed=seed + 101 * (k + 1),
                                         distance_um=d, gammas_per_event=1)
        if scale is None:
            w = np.asarray(g["weight"], dtype=float)[np.asarray(g["counted"], dtype=bool)]
            scale = w.sum() / (w * w).sum()
        runs[d] = [_Run(exp, _SimGammas(g, scale), d)]
    return DistanceSeries(exp.to_dict(), runs)


@pytest.fixture(scope="module")
def ru98_closure():
    """The coulex_ru98_plunger example (τ = 8.36 ps) with and without its deorientation: seeded simulated series
    from 200 000 reactions."""
    from physim.nuclear import events

    out = {}
    base = Experiment.example("coulex_ru98_plunger").to_dict()
    base["lifetime"]["distances"] = [f"{d:g} um" for d in RU98_DISTANCES]
    ev = events.simulate(Experiment.from_dict(base), 200_000, 7)
    for name in ("with", "without"):
        d = {**base, "lifetime": dict(base["lifetime"])}
        if name == "without":
            d["lifetime"].pop("deorientation")
        exp = Experiment.from_dict(d)
        out[name] = (Planner(exp), _ru98_series(exp, ev, 7))
    return out


@pytest.mark.parametrize("deorientation", ["with", "without"])
def test_ru98_tau_comes_back_from_a_simulated_series(ru98_closure, deorientation):
    """Planner.analyse_lifetimes on simulated ⁹⁸Ru runs gives τ back within its stated uncertainty, with and without
    deorientation in vacuum. Before wave 4B the DDCM came out ~10% long: Gaussian or window areas with one recoil
    speed (the runs' speeds spread over β ≈ 4–8%: the silicon detector sees the scattered beam and the recoils), no
    efficiency ratio and no deorientation correction."""
    p, series = ru98_closure[deorientation]
    tau = 8.36
    assert p.lifetimes()["levels"][0]["tau_ps"] == pytest.approx(tau, rel=0.01)
    out = p.analyse_lifetimes(series=series)
    assert out["available"], out["reason"]
    assert out["method"] == "ddcm"
    assert abs(out["tau_ps"] - tau) <= 2 * out["total_ps"], (out["tau_ps"], out["total_ps"])
    assert out["total_ps"] < 0.08 * tau
    notes = " ".join(out["notes"])
    assert "weighted by their own speed" in notes and "Efficiency of the shifted" in notes
    assert ("Deorientation in vacuum corrected" in notes) == (deorientation == "with")
    dc = p.analyse_lifetimes(series=series, method="decay_curve")
    assert abs(dc["tau_ps"] - tau) <= 2 * dc["total_ps"] + 0.03 * tau, (dc["tau_ps"], dc["total_ps"])
