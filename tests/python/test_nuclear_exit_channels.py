"""Exit channels from physim's own fusion and evaporation (backlog item 75, the wiring of Stages 2 and 3 into
Stage 1): the channel table they make, its bookkeeping against σ_fus, the planner filling the setup from it, the
Plan's fusion numbers, and the MCP tool."""

import json

import numpy as np
import pytest

from physim.nuclear.exit_channels import SOURCE, fusion_summary, physim_channels, thick_target
from physim.nuclear.experiment import Experiment
from physim.nuclear.fusion import Fusion
from physim.nuclear.planner import Planner

EXAMPLE = "fusion_dy158"
EVENTS = 4000


@pytest.fixture(scope="module")
def exp():
    return Experiment.example(EXAMPLE)


@pytest.fixture(scope="module")
def made(exp):
    return physim_channels(exp, events=EVENTS, seed=3)


def test_the_table_accounts_for_the_whole_fusion_cross_section(made):
    """Σσ of the channels plus the shares left out (fission, no particle, small channels) times σ_fus is σ_fus."""
    left = sum(made.left_out.values())
    assert made.table.total_mb + left * made.sigma_fusion_mb == pytest.approx(made.sigma_fusion_mb, rel=1e-9)
    assert all(ch.source == SOURCE for ch in made.table.channels)
    sigmas = [ch.cross_section_mb for ch in made.table.channels]
    assert sigmas == sorted(sigmas, reverse=True)


def test_sigma_fusion_is_the_fusion_codes_thick_target_average(exp, made):
    direct = Fusion("12C", "150Nd").thick_target(exp.beam.energy_mev, 1.0, "150Nd")
    assert made.sigma_fusion_mb == pytest.approx(direct.sigma_mb, rel=1e-9)
    assert made.sigma_fusion_mb == pytest.approx(thick_target(exp).sigma_mb, rel=1e-12)
    lo, hi = made.excitation_range_mev
    assert lo <= made.excitation_mev <= hi


def test_12c_on_150nd_at_60_mev_is_mostly_4n(made):
    """E* ≈ 49 MeV in ¹⁶²Dy: 4n strongest, then 5n (as in the example's own table and in PACE4 at this energy)."""
    labels = [ch.label for ch in made.table.channels]
    assert labels[:2] == ["4n", "5n"]
    first = made.table.channels[0]
    assert first.residue == "158Dy"
    assert first.entry_spin == int(first.entry_spin)             # even A: a whole spin
    assert made.table.channels[1].entry_spin % 1 == 0.5           # 157Dy: half-whole


def test_seeded_tables_are_identical(exp, made):
    again = physim_channels(exp, events=EVENTS, seed=3)
    assert again.table.to_setup() == made.table.to_setup()


def test_below_the_barrier_nothing_fuses(exp):
    d = exp.to_dict()
    d["beam"]["energy"] = "20 MeV"
    with pytest.raises(ValueError, match="nothing fuses"):
        physim_channels(Experiment.from_dict(d), events=100)


def test_heavy_beams_carry_the_high_spin_warning():
    d = Experiment.example(EXAMPLE).to_dict()
    d["beam"]["nuclide"] = "32S"
    d["beam"]["energy"] = "150 MeV"
    d["target"]["material"] = "92Mo"
    d["reaction"]["channels"] = [{"particles": "2n", "cross_section": "10 mb"}]
    d.pop("levels", None)
    made = physim_channels(Experiment.from_dict(d), events=500, seed=1)
    assert any("High spin" in w for w in made.warnings)


def test_fusion_summary_gives_the_barrier_and_an_excitation_function(exp):
    f = fusion_summary(exp)
    assert 40 < f["barrier_lab_mev"] < 55
    assert f["sigma_thick_mb"] < f["sigma_mb"]                    # the beam slows down through the target
    ef = f["excitation_function"]
    assert np.all(np.diff(ef["sigma_mb"]) >= 0)                   # rises with energy over this range
    assert f["e_out_mev"] < f["e_in_mev"] == exp.beam.energy_mev


def test_the_planner_fills_the_setup_from_physim(exp):
    p = Planner(Experiment.example(EXAMPLE))
    assert p.fill_channels(events=EVENTS, seed=3), p.problems
    chans = p.experiment.fusion.channels
    assert [c.particles for c in chans][:2] == ["4n", "5n"]
    assert all(c.source == SOURCE for c in chans)
    rows = p.residues()["rows"]
    assert rows[0]["residue"] == "158Dy" and rows[0]["source"] == SOURCE
    # The setup keeps its shape: it writes and reads back unchanged.
    again = Experiment.from_dict(p.experiment.to_dict())
    assert [(c.particles, c.cross_section_mb) for c in again.fusion.channels] == \
        [(c.particles, c.cross_section_mb) for c in chans]


def test_the_plan_carries_plain_fusion_numbers(exp):
    p = Planner(Experiment.example(EXAMPLE))
    pl = p.plan()
    f = pl["fusion"]
    assert f["compound"] == "162Dy" and "excitation_function" not in f
    json.dumps(f)
    assert Planner(Experiment.example("alpha_on_gold")).plan()["fusion"] is None


def test_the_mcp_tool_fills_the_channels(tmp_path):
    pytest.importorskip("mcp.server.mcpserver", reason="the MCP extra is not installed")
    from physim.mcp.session import Session
    from physim.mcp.tools_setup import SetupTools

    tools = SetupTools(Session(tmp_path))
    tools.new_experiment("Dy", from_example=EXAMPLE)
    out = tools.fill_channels(events=EVENTS, seed=3)
    assert out["ok"]
    assert out["physim"]["sigma_fusion_mb"] > 0
    assert out["residues"][0].startswith("4n -> 158Dy")


# -- the neutron wall and the channel selection ---------------------------------------------------------------


def _with_wall(**wall):
    d = Experiment.example(EXAMPLE).to_dict()
    d["neutron_wall"] = wall or {"cells": 50, "distance": "510 mm"}
    return d


def test_the_neutron_wall_reads_writes_and_takes_typical_values():
    exp = Experiment.from_dict(_with_wall())
    assert exp.to_dict()["neutron_wall"] == {"cells": 50, "distance": "510 mm"}
    wall = exp.neutron_wall.model()
    assert wall.cell.material == "EJ-301" and wall.pitch_mm == pytest.approx(127.0)
    assert 0 < wall.geometric_efficiency() < 1


def test_the_neutron_wall_reports_its_problems():
    with pytest.raises(Exception, match="scintillator"):
        Experiment.from_dict(_with_wall(cells=50, distance="510 mm", material="water"))
    with pytest.raises(Exception, match="whole sphere"):
        Experiment.from_dict(_with_wall(cells=5000, distance="510 mm"))
    d = Experiment.example("alpha_on_gold").to_dict()
    d["neutron_wall"] = {"cells": 10, "distance": "500 mm"}
    with pytest.raises(Exception, match="fusion-evaporation"):
        Experiment.from_dict(d)


def test_channel_selection_sorts_events_by_what_is_seen():
    from physim.nuclear.exit_channels import channel_selection

    r = channel_selection(Experiment.from_dict(_with_wall()), events=20_000, seed=2)
    m = np.array(r["matrix"])
    assert np.allclose(m.sum(axis=1), 1.0)
    assert r["true"] == ["3n", "4n", "5n", "p3n", "a2n"]
    i4, zero, one = r["true"].index("4n"), r["assigned"].index("γ"), r["assigned"].index("n")
    # 4n: each neutron is seen with about the wall's efficiency at its energy, so most events show none or one
    eps = r["neutron"]["efficiency_at_mean"]
    assert m[i4, zero] == pytest.approx((1 - eps) ** 4, abs=0.05)
    assert m[i4, one] > m[i4, r["assigned"].index("2n")]
    # without particle detectors at the right angles, protons and α particles are not seen
    assert r["charged_reached"] == {"p": 0.0, "a": 0.0}
    for c in r["contents"].values():
        assert sum(c["from"].values()) == pytest.approx(1.0)


def _with_silicon(**wall):
    """The example with its wall and two 300 µm silicon annular detectors, forward and backward (θ 17–50°, 130–163°)."""
    d = _with_wall(cells=50, distance="510 mm", **wall)
    for name, theta in (("CDf", "0 deg"), ("CDb", "180 deg")):
        d["detectors"].append({"name": name, "shape": "annular", "theta": theta, "distance": "40 mm",
                               "inner_radius": "12 mm", "outer_radius": "48 mm", "thickness": "300 um",
                               "dead_layer": "0.5 um"})
    return Experiment.from_dict(d)


def _z_max(r):
    """Largest |Monte Carlo − analytic| in units of the binomial error of the analytic share."""
    a, m = np.array(r["analytic"]), np.array(r["matrix"])
    n = np.array(r["events"])[:, None]
    return float((np.abs(m - a) / (np.sqrt(a * (1 - a) / n) + 1 / n)).max())


def test_the_new_wall_fields_read_write_and_check():
    exp = _with_silicon(gamma_efficiency="30 %", gamma_multiplicity=12, particle_threshold="0.5 MeV")
    w = exp.to_dict()["neutron_wall"]
    assert w["gamma_efficiency"] == "30 %" and w["gamma_multiplicity"] == 12 and w["particle_threshold"] == "0.5 MeV"
    assert exp.neutron_wall.model().cell.gamma_efficiency == "30 %"
    with pytest.raises(Exception, match="gamma_multiplicity"):
        Experiment.from_dict(_with_wall(cells=50, distance="510 mm", gamma_multiplicity=-1))
    with pytest.raises(Exception, match="particle_threshold"):
        Experiment.from_dict(_with_wall(cells=50, distance="510 mm", particle_threshold="-1 MeV"))


def test_particle_detector_thresholds_cut_protons_that_punch_through():
    """Evaporated protons from ¹⁶²Dy carry 10–17 MeV and punch through 300 µm of silicon, leaving about 2–3 MeV: a
    1 MeV threshold keeps them all, a 3 MeV one loses most; the α particles leave 17–25 MeV, above either."""
    from physim.nuclear.exit_channels import PARTICLE_THRESHOLD_MEV, channel_selection

    low = channel_selection(_with_silicon(), events=20_000, seed=4)
    high = channel_selection(_with_silicon(particle_threshold="3 MeV"), events=20_000, seed=4)
    assert low["thresholds_mev"]["CDf"] == PARTICLE_THRESHOLD_MEV and low["thresholds_mev"]["Monitor"] == 1.0
    assert high["thresholds_mev"]["CDf"] == 3.0 and high["thresholds_mev"]["Monitor"] == 1.0  # its own threshold
    for r in (low, high):
        assert 0.2 < r["charged_reached"]["p"] < 0.5
        assert r["charged_counted"]["a"] == pytest.approx(r["charged_reached"]["a"])
        assert np.allclose(np.array(r["matrix"]).sum(axis=1), 1.0)
        assert _z_max(r) < 4.5  # the Monte Carlo agrees with the analytic matrix of the same events
    assert low["charged_counted"]["p"] == pytest.approx(low["charged_reached"]["p"], abs=0.01)
    assert high["charged_counted"]["p"] < 0.6 * high["charged_reached"]["p"]
    def p_share(r):  # p3n events that show their proton and no neutron
        return r["matrix"][r["true"].index("p3n")][r["assigned"].index("p")]

    ratio = high["charged_counted"]["p"] / low["charged_counted"]["p"]
    assert p_share(high) == pytest.approx(p_share(low) * ratio, abs=0.04)


def test_gamma_rays_add_false_neutrons():
    """Each γ ray fires a cell with the wall's geometric efficiency × gamma_efficiency and becomes a neutron with
    gamma_as_neutron: a true 4n event that shows no neutron without γ rays shows false ones with 1 − (1 − f)^M."""
    from physim.nuclear.exit_channels import STATISTICAL_GAMMAS, channel_selection

    none = channel_selection(_with_silicon(gamma_multiplicity=0), events=20_000, seed=5)
    many = channel_selection(_with_silicon(gamma_multiplicity=20, gamma_as_neutron="5 %"), events=20_000, seed=5)
    auto = channel_selection(_with_silicon(), events=20_000, seed=5)
    assert none["gamma"]["multiplicity"] == 0 and many["gamma"]["multiplicity"] == 20
    assert auto["gamma"]["from_entry_spin"] and auto["gamma"]["multiplicity"] > STATISTICAL_GAMMAS
    f = many["gamma"]["false_neutron_per_gamma"]
    assert f == pytest.approx(many["neutron"]["geometric_efficiency"] * 0.40 * 0.05)
    def share(r, cls):
        return r["matrix"][r["true"].index("4n")][r["assigned"].index(cls)]

    # P(nothing seen) falls by the chance that no γ ray poses as a neutron, (1 − f)^20
    assert share(many, "γ") == pytest.approx(share(none, "γ") * (1 - f) ** 20, abs=0.012)
    assert share(many, "n") > share(none, "n") + 0.5 * 20 * f * share(none, "γ")
    for r in (none, many, auto):
        assert _z_max(r) < 4.5


def test_the_planner_adds_and_removes_the_wall_and_keeps_it_across_reaction_types():
    p = Planner(Experiment.example(EXAMPLE))
    assert p.add_neutron_wall(), p.problems
    assert p.experiment.neutron_wall.cells == 50
    assert p.set("reaction", "type", "elastic"), p.problems
    assert p.experiment.neutron_wall is None
    assert p.set("reaction", "type", "fusion-evaporation"), p.problems
    assert p.experiment.neutron_wall is not None
    assert p.remove_neutron_wall() and p.experiment.neutron_wall is None


def test_the_mcp_tools_add_a_wall_and_work_out_the_selection(tmp_path):
    pytest.importorskip("mcp.server.mcpserver", reason="the MCP extra is not installed")
    from physim.mcp.session import Session
    from physim.mcp.tools_setup import SetupTools

    tools = SetupTools(Session(tmp_path))
    tools.new_experiment("Dy wall", from_example=EXAMPLE)
    tools.set_value("neutron_wall", "cells", 50)
    out = tools.set_value("neutron_wall", "distance", "510 mm")
    assert out["ok"] and "cells 50" in out["neutron_wall"]
    r = tools.channel_selection(events=5000, seed=1)
    assert r["true"][0] == "3n" and "γ" in r["assigned"]
    json.dumps(r)
