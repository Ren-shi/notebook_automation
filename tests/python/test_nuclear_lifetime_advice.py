"""DSAM backing advice and the beam in the plunger's foils (lifetime wave 4A): ``lifetimes.dsam_sensitivity`` and
``slowing_time``, ``lifetime_plan.backing_candidates``/``backing_advice``/``stopper_beam``/``stopper_warnings``, their
place in ``Planner.lifetimes`` and the banner, and the MCP view.

Each test builds its own setup (typed level schemes, no ENSDF: CI has no local copy, and the stopper's level data is
typed in or hidden), and pins nothing on the three examples' values."""

from __future__ import annotations

import math

import numpy as np
import pytest

from physim.nuclear import Experiment, ensdf, lifetimes
from physim.nuclear import lifetime_plan as lp
from physim.nuclear.levels import LN2, LevelScheme
from physim.nuclear.planner import Planner


@pytest.fixture(autouse=True)
def _no_ensdf(monkeypatch, tmp_path):
    """No local ENSDF copy, as in CI."""
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))


def coulex_plunger(b_up="2500 e2fm4", current="2.5 enA"):
    """The 98Ru-on-24Mg plunger detectors and beam with this test's own B(E2), beam current, distances and a
    20 mg/cm2 Cu stopper (ENSDF is hidden, so the stopper's levels are this test's or none)."""
    d = Experiment.example("coulex_ru98_plunger").to_dict()
    d["reaction"]["b_up"] = b_up
    d["beam"]["current"] = current
    d["lifetime"] = {"method": "rdds", "distances": ["5 um", "20 um", "80 um", "300 um"],
                     "stopper": {"kind": "stopper", "material": "Cu", "thickness": "20 mg/cm2"},
                     "deorientation": {"alpha2": 0.2, "gamma2": "0.1 /ps", "alpha4": 0.2, "gamma4": "0.3 /ps"}}
    return Experiment.from_dict(d)


# -- physics ------------------------------------------------------------------------------------------------------

def test_sensitivity_and_slowing_time_on_the_exponential_kernel():
    alpha = 0.4
    k = lifetimes.exponential_kernel(0.02, alpha)
    for tau in (0.05, 0.2, 0.4, 1.0, 3.0):
        assert lifetimes.dsam_sensitivity(tau, k) == pytest.approx(alpha * tau / (alpha + tau) ** 2, rel=2e-3)
    # Largest (1/4) where the slowing-down time matches τ.
    taus = np.geomspace(0.04, 4.0, 41)
    s = [lifetimes.dsam_sensitivity(t, k) for t in taus]
    assert taus[int(np.argmax(s))] == pytest.approx(alpha, rel=0.15) and max(s) == pytest.approx(0.25, rel=2e-3)
    assert lifetimes.slowing_time(k) == pytest.approx(alpha, rel=1e-3)
    moving = lifetimes.SampledKernel(np.linspace(0.0, 2.0, 201), np.full((201, 1), 0.02))
    assert math.isinf(lifetimes.slowing_time(moving))
    with pytest.raises(ValueError):
        lifetimes.dsam_sensitivity(0.0, k)


# -- DSAM backing advice ------------------------------------------------------------------------------------------

def _band(nuclide, rows, multipolarity="M1"):
    levels, transitions = [], []
    for k, (jpi, e, t_half, e_gamma) in enumerate(rows):
        lv = {"energy": f"{e:g} keV", "jpi": jpi, "source": "user"}
        if t_half is not None:
            lv["half_life"] = f"{t_half:g} ps"
        levels.append(lv)
        if e_gamma is not None:
            transitions.append({"from": k, "to": k - 1, "energy": f"{e_gamma:g} keV", "intensity": 100,
                                "multipolarity": multipolarity, "source": "user"})
    return {"nuclide": nuclide, "level": levels, "transition": transitions}


def _dsam(taus_ps=(0.03, 0.3, 0.6), backing=("Au", "10 mg/cm2")):
    """The 143Sm-on-a-backing DSAM detectors and beam with this test's own dipole band and backing."""
    d = Experiment.example("fusion_sm143_dsam").to_dict()
    d["reaction"]["side_feeding_levels"] = 0
    for ch in d["reaction"]["channels"]:
        if ch["residue"] == "143Sm":
            ch["entry_spin"] = 19.5
    rows = [("13/2-", 0.0, None, None)]
    for k, spin in enumerate((15, 17, 19), start=1):
        rows.append((f"{spin}/2-", 500.0 * k + 30.0 * k * k, taus_ps[k - 1] * LN2, 500.0 + 30.0 * (2 * k - 1)))
    d.setdefault("levels", {})["143Sm"] = _band("143Sm", rows)
    if backing is None:
        d["target"].pop("backing", None)
    else:
        d["target"]["backing"] = {"material": backing[0], "thickness": backing[1]}
    return Experiment.from_dict(d)


@pytest.fixture(scope="module")
def gold():
    return Planner(_dsam())


def test_backing_candidates_stop_the_recoil_and_heavy_ones_stop_it_sooner(gold):
    subject = lp.Subject(gold)
    cands = lp.backing_candidates(subject)
    by = {c["material"]: c for c in cands}
    assert set(lp.BACKING_CANDIDATES) <= set(by)
    for c in cands:
        assert c["thickness_to_stop_mg_cm2"] > 0 and c["time_to_stop_ps"] > 0
        assert c["thickness_mg_cm2"] >= lp.BACKING_MARGIN * c["thickness_to_stop_mg_cm2"] * (1 - 1e-9)
        assert c["thickness_um"] == pytest.approx(10 * c["thickness_mg_cm2"] / c["density_g_cm3"])
        assert 0 < c["slowing_time_ps"] < c["time_to_stop_ps"]
    assert by["Au"]["setup"] and not by["Cu"]["setup"]
    # Dense, heavy backings stop the residue sooner than carbon.
    assert by["Au"]["time_to_stop_ps"] < by["C"]["time_to_stop_ps"]
    assert by["Au"]["slowing_time_ps"] < by["C"]["slowing_time_ps"]


def test_backing_advice_is_ranked_by_sensitivity_and_follows_tau(gold):
    cands = lp.backing_candidates(lp.Subject(gold))
    for tau in (0.03, 0.3, 3.0):
        adv = lp.backing_advice(tau, cands)
        s = [r["sensitivity"] for r in adv]
        assert s == sorted(s, reverse=True) and adv[0]["relative_sensitivity"] == pytest.approx(1.0)
        assert all(0 < r["attenuation"] < 1 and "kernel" not in r for r in adv)
        assert all(r["good"] == (r["relative_sensitivity"] >= lp.BACKING_GOOD) for r in adv)
    # A short life wants the fastest stopping, a long one the slowest: carbon is worst at 30 fs, best at 3 ps
    # among these.
    short = lp.backing_advice(0.03, cands)
    long_ = lp.backing_advice(3.0, cands)
    rank = {r["material"]: i for i, r in enumerate(short)}
    assert rank["C"] == len(short) - 1
    assert long_[0]["material"] == "C"
    # F falls with τ on each backing.
    au = [next(r for r in lp.backing_advice(t, cands) if r["material"] == "Au")["attenuation"] for t in (0.1, 0.5)]
    assert au[0] > au[1]


def test_the_plan_gives_backing_advice_per_dsam_level(gold):
    out = gold.lifetimes()
    assert out["method"] == "dsam"
    rows = [r for r in out["levels"] if r["tau_ps"] is not None]
    assert rows and all(r["backing_advice"] for r in rows)
    for r in rows:
        mine = next(c for c in r["backing_advice"] if c["setup"])
        assert mine["material"] == "Au"
        # The setup's F (its own kernel, more recoils) and the candidate's agree within the sampling.
        assert mine["attenuation"] == pytest.approx(r["attenuation"], abs=0.05)
    # The old keys are all still there.
    assert {"attenuation", "separation", "hours_per_distance", "recommended_distances_um"} <= set(rows[0])


def test_a_light_backing_for_a_short_life_is_warned_about():
    p = Planner(_dsam(taus_ps=(0.12, 0.15, 0.2), backing=("C", "5 mg/cm2")))
    out = p.lifetimes()
    first = out["levels"][0]
    mine = next(c for c in first["backing_advice"] if c["setup"])
    assert mine["material"] == "C" and not mine["good"]
    assert any("C backing gives" in w and "of the best" in w for w in out["warnings"])
    # A good backing is not warned about.
    gold = Planner(_dsam(taus_ps=(0.12, 0.15, 0.2))).lifetimes()
    assert not any("backing gives" in w for w in gold["warnings"])


def test_an_rdds_setup_gets_advice_only_for_its_dsam_levels():
    p = Planner(coulex_plunger(b_up="60000 e2fm4"))  # τ below a picosecond: DSAM's
    row = p.lifetimes()["levels"][0]
    assert row["method"] == "dsam" and row["backing_advice"]
    long_lived = Planner(coulex_plunger()).lifetimes()["levels"][0]
    assert long_lived["method"] == "rdds" and long_lived["backing_advice"] is None


# -- the beam in the plunger's foils ------------------------------------------------------------------------------

def _stopper_scheme(nuclide, gamma_kev, me="0.4 eb"):
    """A two-level scheme typed in for a stopper nucleus: an E2 from the ground state to one level at ``gamma_kev``."""
    problems = []
    scheme = LevelScheme.from_dict({"nuclide": nuclide, "level": [
        {"energy": "0 keV", "jpi": "3/2-", "source": "user"},
        {"energy": f"{gamma_kev:g} keV", "jpi": "1/2-", "half_life": "1 ps", "source": "user"}],
        "transition": [{"from": 1, "to": 0, "energy": f"{gamma_kev:g} keV", "intensity": 100,
                        "multipolarity": "E2", "source": "user"}]}, "test", problems)
    assert not problems, problems
    scheme.set_matrix_element(0, 1, "E2", me)
    return scheme


def test_the_beam_power_in_the_stopper_and_its_warning(monkeypatch):
    monkeypatch.setattr(lp, "_stopper_scheme", lambda exp, nuc: (None, "no level data in this test"))
    low = Planner(coulex_plunger(current="2.5 enA"))
    foils = low.lifetimes()["stopper"]
    assert [f["foil"] for f in foils] == ["stopper"]
    f = foils[0]
    exp = low.experiment
    # The 98Ru beam stops in 20 mg/cm² of Cu: it leaves all its energy there.
    assert f["beam_out_mev"] == 0 and f["deposited_mev"] == pytest.approx(f["beam_in_mev"])
    assert 0 < f["beam_in_mev"] < exp.beam.energy_mev
    assert f["power_w"] == pytest.approx(exp.beam.particles_per_second * f["deposited_mev"] * 1.602176634e-13)
    assert f["power_limit_w"] == lp.STOPPER_POWER_W and f["power_w"] < f["power_limit_w"]
    assert not any("deposits" in w for w in low.lifetimes()["warnings"])
    high = Planner(coulex_plunger(current="400 enA"))
    hot = high.lifetimes()["stopper"][0]
    assert hot["power_w"] == pytest.approx(160 * f["power_w"], rel=1e-9) and hot["power_w"] > lp.STOPPER_POWER_W
    text = next(w for w in high.lifetimes()["warnings"] if "deposits" in w)
    assert f"{lp.STOPPER_POWER_W:g} W" in text and "Cu stopper" in text
    assert any("deposits" in w.text and w.level == "warning" for w in high.warnings())  # the banner too


def test_the_stoppers_own_coulex_lines_near_a_studied_line(monkeypatch):
    near = {"63Cu": _stopper_scheme("63Cu", 653.5), "65Cu": _stopper_scheme("65Cu", 1500.0)}
    monkeypatch.setattr(lp, "_stopper_scheme", lambda exp, nuc: (near[nuc], "typed in") if nuc in near
                        else (None, "none"))
    p = Planner(coulex_plunger())
    f = p.lifetimes()["stopper"][0]
    by = {x["nuclide"]: x for x in f["lines"]}
    assert set(by) == {"63Cu", "65Cu"} and not f["missing"]
    assert by["63Cu"]["gamma_kev"] == pytest.approx(653.5)
    assert by["63Cu"]["yield_per_beam"] > 0 and by["63Cu"]["rate_per_s"] > 0
    warn = [w for w in p.lifetimes()["warnings"] if "Coulomb-excites" in w]
    assert len(warn) == 1 and "63Cu 653.5 keV" in warn[0] and "65Cu" not in warn[0] and "unshifted" in warn[0]
    assert any("Coulomb-excites" in w.text for w in p.warnings())
    # The same line far from every studied peak is listed but not warned about.
    far = {"63Cu": _stopper_scheme("63Cu", 900.0), "65Cu": _stopper_scheme("65Cu", 1500.0)}
    monkeypatch.setattr(lp, "_stopper_scheme", lambda exp, nuc: (far[nuc], "typed in"))
    q = Planner(coulex_plunger())
    assert q.lifetimes()["stopper"][0]["lines"]
    assert not any("Coulomb-excites" in w for w in q.lifetimes()["warnings"])


def test_coulex_yield_grows_with_the_matrix_element(monkeypatch):
    yields = []
    for me in ("0.2 eb", "0.4 eb"):
        s = {"63Cu": _stopper_scheme("63Cu", 900.0, me)}
        monkeypatch.setattr(lp, "_stopper_scheme", lambda exp, nuc, s=s: (s.get(nuc), "typed in" if nuc in s else
                                                                          "none"))
        yields.append(Planner(coulex_plunger()).lifetimes()["stopper"][0]["lines"][0]["yield_per_beam"])
    # First-order Coulomb excitation goes as the matrix element squared.
    assert yields[1] / yields[0] == pytest.approx(4.0, rel=0.05)


def test_a_stopper_without_level_data_is_said_not_guessed():
    p = Planner(coulex_plunger())  # ENSDF hidden
    f = p.lifetimes()["stopper"][0]
    assert not f["lines"] and {m["nuclide"] for m in f["missing"]} == {"63Cu", "65Cu"}
    assert any(w.startswith("No level data for the plunger's Cu stopper") for w in p.lifetimes()["warnings"])
    note = next(w for w in p.warnings() if w.text.startswith("No level data"))
    assert note.level == "note"


def test_the_stoppers_scheme_can_come_from_the_setup():
    exp = Experiment.example("coulex_ru98_plunger")
    exp.levels = dict(exp.levels, target=_stopper_scheme("63Cu", 700.0))
    scheme, source = lp._stopper_scheme(exp, "63Cu")
    assert scheme is exp.levels["target"] and "setup" in source
    assert lp._stopper_scheme(exp, "65Cu")[0] is None  # ENSDF hidden


def test_a_differential_plunger_has_two_foils_in_the_beam(monkeypatch):
    monkeypatch.setattr(lp, "_stopper_scheme", lambda exp, nuc: (None, "none"))
    p = Planner(Experiment.example("fusion_cr50_plunger"))
    foils = p.lifetimes()["stopper"]
    assert [f["foil"] for f in foils] == ["degrader", "degrader2"]
    # The beam crosses the first foil and enters the second with what it kept.
    assert foils[0]["beam_out_mev"] == pytest.approx(foils[1]["beam_in_mev"])
    assert all(f["deposited_mev"] == pytest.approx(f["beam_in_mev"] - f["beam_out_mev"]) for f in foils)


def test_dsam_and_no_lifetime_setups_have_no_foils():
    assert Planner(_dsam()).lifetimes()["stopper"] == []
    assert lp.stopper_warnings(Planner.example("alpha_on_gold")) == []


# -- MCP ------------------------------------------------------------------------------------------------------------

def test_the_mcp_view_shows_the_advice_and_the_foils(gold, monkeypatch):
    from physim.mcp.tools_lifetimes import plan_view

    view = plan_view(gold.lifetimes())
    adv = view["levels"][0]["backing_advice"]
    assert adv[0]["rank"] == 1 and "mg/cm2" in adv[0]["thickness_advised"]
    assert any(r.get("setup_backing") for r in adv)
    assert "beam_in_plunger_foils" not in view
    monkeypatch.setattr(lp, "_stopper_scheme", lambda exp, nuc: (None, "none"))
    rv = plan_view(Planner(coulex_plunger()).lifetimes())
    foil = rv["beam_in_plunger_foils"][0]
    assert "Cu stopper" in foil["foil"] and "stops in it" in foil["beam"] and "limit 0.5 W" in foil["power"]
    assert foil["no_level_data"] == ["63Cu", "65Cu"]
