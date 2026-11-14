"""Fusion-evaporation, Stage 1 (backlog item 75): residues from a channel table.

The "Done when" checks of Stage 1: the compound nucleus's excitation energy and speed equal the hand calculation
within 1 keV and 1e-6 c; the residues counted per channel equal the table × the luminosity within statistics; the
Doppler-corrected width of a residue line equals what the residues' opening angle and the crystal's size give;
seeded runs are bit-identical. Also: the setup section and its problems, the channel table other codes fill,
momentum and energy conservation event by event, the cascade and its alignment, and the planner, the app's model
and the MCP tools on the new reaction."""

import copy
import math

import numpy as np
import pytest

from physim.nuclear import data
from physim.nuclear.experiment import Experiment, ExitChannel, SetupError
from physim.nuclear.kinematics import TwoBody
from physim.nuclear.planner import Planner
from physim.nuclear.residues import (LIGHT, ChannelTable, ResidueCascade, _gaussian_tensors, compound,
                                     compound_profile, parse_particles, particles_label, residue_table,
                                     simulate_residue_gammas, simulate_residues)
from physim.nuclear.stopping import Stopping

EXAMPLE = "fusion_dy158"


@pytest.fixture(scope="module")
def exp():
    return Experiment.example(EXAMPLE)


@pytest.fixture(scope="module")
def events(exp):
    return simulate_residues(exp, events=200_000, seed=11)


def _setup(exp, **reaction):
    d = exp.to_dict()
    d["reaction"].update(reaction)
    return d


# -- the compound nucleus: exact kinematics --------------------------------------------------------------------


def test_excitation_and_speed_equal_the_hand_calculation():
    """E* = √s − M_CN and β = p_b / (E_b + m_t) from the AME2020 masses, worked out here step by step; and the same E*
    as E_cm + Q, with Q from the atomic mass excesses and E_cm from the two-body kinematics."""
    for beam, target, e in (("12C", "150Nd", 60.0), ("16O", "148Nd", 85.0), ("48Ca", "124Sn", 200.0)):
        mb, mt = data.nuclide(beam).nuclear_mass_mev, data.nuclide(target).nuclear_mass_mev
        b, t = data.nuclide(beam), data.nuclide(target)
        mc = data.nuclide((b.Z + t.Z, b.A + t.A)).nuclear_mass_mev
        p = math.sqrt(e**2 + 2 * e * mb)
        total = e + mb + mt
        roots = math.sqrt(total**2 - p**2)
        cn = compound(beam, target, e)
        assert abs(cn.excitation - (roots - mc)) < 1e-3          # 1 keV
        assert abs(cn.beta - p / total) < 1e-6
        # Q of the nuclei: the atomic Q of the mass excesses, less the binding energy the electrons gain when the
        # compound atom forms (an atomic energy, not the nucleus's: 68 keV for 162Dy).
        q_atomic = data.q_value([beam, target], [data.nuclide((b.Z + t.Z, b.A + t.A)).name])
        electrons = 1e-3 * (data.electron_binding_kev(b.Z) + data.electron_binding_kev(t.Z)
                            - data.electron_binding_kev(b.Z + t.Z))
        e_cm = TwoBody(beam, target, e).cm_energy_mev
        assert abs(cn.excitation - (e_cm + q_atomic + electrons)) < 1e-3
        assert abs(cn.q_value - (q_atomic + electrons)) < 1e-3
        assert abs(cn.e_cm - e_cm) < 1e-3
    # The example: 12C at 60 MeV on 150Nd makes 162Dy at 49.98 MeV, Q = −5.567 MeV.
    cn = compound("12C", "150Nd", "60 MeV")
    assert cn.nuclide == "162Dy"
    assert cn.excitation == pytest.approx(49.977, abs=2e-3)
    assert cn.q_value == pytest.approx(-5.567, abs=2e-3)
    # The non-relativistic speed v_b m_b / (m_b + m_t) is 1e-4 off: it is not the hand calculation within 1e-6 c.
    mb, mt = data.nuclide("12C").nuclear_mass_mev, data.nuclide("150Nd").nuclear_mass_mev
    v_beam = math.sqrt(1 - (mb / (60.0 + mb)) ** 2)
    assert 1e-6 < abs(v_beam * mb / (mb + mt) - cn.beta) < 1e-4 * 2


def test_excitation_through_the_target_follows_the_stopping(exp):
    """At the front face E* is that of the beam energy; at the back face that of the beam after the target's
    thickness, from the stopping powers; the spread grows with straggling."""
    prof = compound_profile(exp)
    e0 = exp.beam.energy_mev
    assert prof["beam_energy"][0] == pytest.approx(e0)
    assert abs(prof["excitation"][0] - compound("12C", "150Nd", e0).excitation) < 1e-3
    e_out = float(Stopping("12C", exp.target.material_data()).energy_after(e0, "1 mg/cm2"))
    assert abs(prof["beam_energy"][-1] - e_out) < 2e-3
    assert abs(prof["excitation"][-1] - compound("12C", "150Nd", e_out).excitation) < 2e-3
    assert prof["excitation"][-1] < prof["excitation"][0]
    assert prof["beam_sigma"][-1] > prof["beam_sigma"][0] > 0
    assert np.all(np.diff(prof["beta"]) < 0)


# -- the channel table and its residues ------------------------------------------------------------------------


def test_counts_per_channel_equal_the_table_times_the_luminosity(events):
    """Each channel's residues per second agree with σ × luminosity within statistics, and the total exactly."""
    ev = events
    assert ev.rate()[0] == pytest.approx(ev.expected(), rel=1e-9)
    for label in ev.channels:
        r, err = ev.rate(label)
        assert abs(r - ev.expected(label)) < 4 * err, label
    # The luminosity is beam particles per second × target nuclei per cm²: 5 pnA of 12C on 1 mg/cm2 of 150Nd.
    pps = 5e-9 / 1.602176634e-19
    atoms = 1e-3 * 6.02214076e23 / data.nuclide("150Nd").atomic_mass_u
    assert ev.luminosity == pytest.approx(pps * atoms, rel=2e-3)
    assert ev.counts("4n") == pytest.approx(ev.rate("4n")[0] * 48 * 3600)


def test_momentum_and_energy_are_conserved_event_by_event(events):
    """The light particles and the residue carry the compound nucleus's momentum and total energy."""
    ev, p = events, events.particles
    n = len(ev)
    masses = np.array([data.nuclide(LIGHT[k]).nuclear_mass_mev for k in LIGHT])
    m = masses[p["kind"]]
    pp = np.sqrt(p["energy"] ** 2 + 2 * p["energy"] * m)
    th, ph = np.radians(p["theta"]), np.radians(p["phi"])
    vec = pp[:, None] * np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=1)
    total = np.zeros((n, 3))
    np.add.at(total, p["event"], vec)
    energy = np.zeros(n)
    np.add.at(energy, p["event"], p["energy"] + m)
    table = ev.table
    residue_mass = np.array([data.nuclide(table.channels[c].residue).nuclear_mass_mev for c in ev["channel"]])
    m_res = residue_mass + ev["entry_excitation"]
    g = 1 / np.sqrt(1 - ev["beta"] ** 2)
    v = ev["_velocity"]
    total += (m_res * g)[:, None] * v
    energy += m_res * g
    m_cn = data.nuclide("162Dy").nuclear_mass_mev + ev["excitation"]
    g_cn = 1 / np.sqrt(1 - ev["beta_cn"] ** 2)
    p_cn = m_cn * g_cn * ev["beta_cn"]
    assert np.allclose(total[:, 0], 0.0, atol=1e-6 * p_cn.mean())
    assert np.allclose(total[:, 1], 0.0, atol=1e-6 * p_cn.mean())
    assert np.allclose(total[:, 2], p_cn, rtol=1e-8)
    assert np.allclose(energy, m_cn * g_cn, rtol=1e-12)
    # Every channel emits its particles: 4n events have four neutrons, a2n events an α and two neutrons.
    per_event = np.bincount(p["event"], minlength=n)
    for i, ch in enumerate(table.channels):
        assert np.all(per_event[ev["channel"] == i] == len(ch.particles))
    assert np.all(ev["entry_excitation"] >= 0)


def test_the_residues_opening_angle_follows_from_the_particles_momenta(events):
    """The residue's direction is the sum of the recoils: its mean squared angle to the beam is Σ ⟨p⊥²⟩ / P², with
    ⟨p⊥²⟩ = (2/3) ⟨p²⟩ = (2/3) 2μ⟨ε⟩ for each particle (isotropic emission)."""
    ev, p = events, events.particles
    i = ev.table.index("4n")
    rows = ev["channel"] == i
    sel = np.isin(p["event"], np.flatnonzero(rows))
    m_n = data.nuclide("n").nuclear_mass_mev
    mean_eps = p["energy_cm"][sel].mean()
    m_cn = data.nuclide("162Dy").nuclear_mass_mev
    p_cn = m_cn * ev["beta_cn"][rows].mean()
    psi2 = 4 * (2.0 / 3.0) * 2 * m_n * mean_eps / p_cn**2          # four neutrons
    measured = np.mean(np.radians(ev["theta"][rows]) ** 2)
    assert measured == pytest.approx(psi2, rel=0.03)
    # The evaporation spectrum ε e^(−ε/T): the neutrons' mean energy in the emitters' frames is about 2T.
    assert 1.5 < mean_eps < 3.5


def test_charged_particles_leave_above_their_barrier(events):
    p = events.particles
    kinds = list(LIGHT)
    for k in ("p", "a"):
        m = p["kind"] == kinds.index(k)
        assert m.any()
        assert p["energy_cm"][m].min() > 0.5 * p["energy_cm"][m].mean()
    assert p["energy_cm"][p["kind"] == kinds.index("a")].mean() > 12.0


def test_seeded_runs_are_bit_identical(exp):
    a = simulate_residues(exp, events=20_000, seed=5)
    b = simulate_residues(exp, events=20_000, seed=5)
    c = simulate_residues(exp, events=20_000, seed=6)
    for k in a.columns:
        assert np.array_equal(a.columns[k], b.columns[k]), k
    for k in a.particles:
        assert np.array_equal(a.particles[k], b.particles[k]), k
    assert not np.array_equal(a["depth"], c["depth"])
    ga = simulate_residue_gammas(exp, seed=5, residue_events=a)
    gb = simulate_residue_gammas(exp, seed=5, residue_events=b)
    assert len(ga) > 100
    for k in ga.columns:
        assert np.array_equal(ga.columns[k], gb.columns[k]), k


# -- the γ rays: Doppler correction, cascade, alignment -----------------------------------------------------------


def _one_crystal(exp, entry_spin=6):
    """The example with only the 4n channel entering at 6⁺ without side feeding, and one bare germanium crystal at
    90° with no resolution of its own: what is left of a line's width is the Doppler broadening."""
    d = exp.to_dict()
    d["reaction"]["channels"] = [{"residue": "158Dy", "particles": "4n", "cross_section": "500 mb",
                                  "entry_spin": entry_spin}]
    d["reaction"]["side_feeding"] = "0 %"
    d["levels"] = {"158Dy": d["levels"]["158Dy"]}
    d["beam"].pop("spot_size")
    d["gamma_detectors"] = [{"name": "Ge90", "theta": "90 deg", "phi": "0 deg", "distance": "200 mm",
                             "radius": "35 mm", "resolution": "0 keV"}]
    return Experiment.from_dict(d)


def test_doppler_corrected_width_equals_the_geometry():
    """At 90° the corrected line's width is E₀ β √(⟨ψ_x²⟩ + α²/4): ψ_x the residue's angle towards the crystal (the
    opening angle of the residues), α the crystal's half-angle (the 1D spread of a disc is α/2)."""
    e = _one_crystal(Experiment.example(EXAMPLE))
    ev = simulate_residues(e, events=200_000, seed=3)
    g = simulate_residue_gammas(e, seed=3, residue_events=ev, gammas_per_event=10)
    scheme = e.levels["158Dy"]
    six = next(i for i, lv in enumerate(scheme.levels) if lv.spin == 6)
    m = g.select("Ge90", full_energy=True) & (g["initial"] == six)
    x = g["corrected"][m]
    assert len(x) > 1500
    e0 = g["energy0"][m][0]
    v = ev["_velocity"]
    psi_x = v[:, 0] / np.linalg.norm(v, axis=1)
    alpha = math.atan2(35.0, 200.0)
    beta = ev["beta"].mean()
    expected = e0 * beta * math.sqrt(np.mean(psi_x**2) + alpha**2 / 4)
    assert x.std() == pytest.approx(expected, rel=0.05)
    assert x.mean() == pytest.approx(e0, abs=0.1 * expected)
    # Uncorrected, a crystal at 37° sees the line shifted up by about E₀ β cos 37°.
    assert g.beta_nominal == pytest.approx(compound("12C", "150Nd", ev["beam_energy"].mean()).beta, rel=0.01)


def test_cascade_yields_and_feeding(exp):
    """Without side feeding every yrast transition below the entry is crossed once (γ share 1/(1+α)); with side
    feeding the top one gets 1 − f and the lower levels share f."""
    scheme = exp.levels["158Dy"]
    cas = ResidueCascade(scheme, 0.3, 0.0)
    two = next(i for i, lv in enumerate(scheme.levels) if lv.spin == 2)
    t20 = next(t for t in scheme.transitions if t.initial == two)
    alpha = t20.conversion.value if t20.conversion is not None else 0.0
    ys = {(a, b): y for _, y, a, b in cas.line_yields(10)}
    assert ys[(two, t20.final)] == pytest.approx(1 / (1 + alpha), rel=1e-9)
    feed = ResidueCascade(scheme, 0.3, 0.25).feeding(10.5)
    ten = next(i for i, lv in enumerate(scheme.levels) if lv.spin == 10)
    assert feed[ten] == pytest.approx(0.75)
    assert sum(feed.values()) == pytest.approx(1.0)
    # The Monte Carlo cascade gives the same yields.
    rng = np.random.default_rng(1)
    out = cas.run(np.full(20_000, 10.0), rng)
    n20 = np.sum((out["initial"] == two) & (out["final"] == t20.final))
    assert n20 / 20_000 == pytest.approx(1 / (1 + alpha), abs=4 * math.sqrt(0.25 / 20_000))


def test_alignment_perpendicular_to_the_beam():
    """A spin perpendicular to the beam (small m about it) has B₂ < 0; a stretched E2 then has more γ rays along
    the beam than across it (a₂ = B₂ A₂ > 0). Down a cascade of stretched E2 transitions a₂ stays the same, as
    U₂ A₂ of each γ ray is A₂ of the one before: the yrast lines of a fusion-evaporation residue all show about the
    same anisotropy."""
    b = _gaussian_tensors(10.0, 0.3)
    assert b[2] < -0.5
    assert _gaussian_tensors(0.0, 0.3) == {2: 0.0, 4: 0.0}
    exp = Experiment.example(EXAMPLE)
    cas = ResidueCascade(exp.levels["158Dy"], 0.3, 0.0)
    out = cas.run(np.full(10, 10.0), np.random.default_rng(2))
    first = out["source"] == 0
    a2 = out["a2"][first]
    assert a2[0] > 0.2
    assert np.allclose(a2, a2[0], rtol=1e-9)
    assert a2[0] == pytest.approx(b[2] * -0.36861280877459485, rel=1e-6)  # B₂(10) A₂(10 → 8)
    # Isotropic without alignment information (a level of unknown spin) is the fall-back.
    assert cas.run(np.zeros(5), np.random.default_rng(1))["energy"].size == 0


# -- the setup section ----------------------------------------------------------------------------------------


def test_the_example_reads_and_writes_back_unchanged(exp):
    assert exp.reaction == "fusion-evaporation"
    assert [ch.residue_name for ch in exp.fusion.channels] == ["159Dy", "158Dy", "157Dy", "158Tb", "156Gd"]
    text = exp.to_toml()
    assert Experiment.from_toml(text).to_toml() == text
    assert "[[reaction.channels]]" in text and "[levels.158Dy]" in text


def test_particles_are_read_in_every_usual_spelling():
    assert parse_particles("4n") == ("n",) * 4
    assert parse_particles("p2n") == ("p", "n", "n")
    assert parse_particles("a n") == parse_particles("αn") == ("a", "n")
    assert particles_label(("n", "p", "n")) == "p2n"
    assert particles_label(("a", "n", "n")) == "α2n"
    for bad in ("", "4x", "n0", "2"):
        with pytest.raises(ValueError):
            parse_particles(bad)


@pytest.mark.parametrize("change, expected", [
    ({"channels": None}, "needs its exit channels"),
    ({"channels": []}, "channel table is empty"),
    ({"channels": [{"particles": "4x", "cross_section": "1 mb"}]}, "reaction channel 1 (4x): particles"),
    ({"channels": [{"particles": "4n", "cross_section": "300"}]}, "cross_section '300' has no unit"),
    ({"channels": [{"particles": "4n", "cross_section": "3 MeV"}]}, "cross_section '3 MeV' has the wrong unit"),
    ({"channels": [{"particles": "4n", "residue": "158Er", "cross_section": "3 mb"}]},
     "residue 158Er is not what 162Dy leaves after emitting 4n (158Dy)"),
    ({"channels": [{"particles": "4n", "cross_section": "3 mb"}, {"particles": "nnnn", "cross_section": "3 mb"}]},
     "the same particles as channel 1"),
    ({"channels": [{"particles": "4n", "cross_sectoin": "3 mb"}]}, "did you mean 'cross_section'"),
    ({"channels": [{"particles": "4n", "cross_section": "3 mb", "entry_spin": 7.3}]}, "entry_spin must be"),
    ({"alignment": 5.0}, "alignment must be above 0 and at most 2"),
])
def test_setup_problems_name_their_field(exp, change, expected):
    d = exp.to_dict()
    for k, v in change.items():
        if v is None:
            d["reaction"].pop(k)
        else:
            d["reaction"][k] = v
    with pytest.raises(SetupError) as err:
        Experiment.from_dict(d)
    assert any(expected in p for p in err.value.problems), err.value.problems


def test_fusion_fields_on_another_reaction_are_problems(exp):
    d = exp.to_dict()
    d["reaction"]["type"] = "elastic"
    d.pop("levels")
    with pytest.raises(SetupError) as err:
        Experiment.from_dict(d)
    assert any("only apply to type = \"fusion-evaporation\"" in p for p in err.value.problems)
    d = exp.to_dict()
    d["reaction"]["energy"] = "1 MeV"
    with pytest.raises(SetupError) as err:
        Experiment.from_dict(d)
    assert any("energy only apply to type = \"coulex\"" in p for p in err.value.problems)


def test_level_schemes_belong_to_a_residue(exp):
    d = exp.to_dict()
    d["levels"]["158Er"] = copy.deepcopy(d["levels"]["158Dy"])
    with pytest.raises(SetupError) as err:
        Experiment.from_dict(d)
    assert any("levels: unknown section '158Er'" in p for p in err.value.problems)
    d = exp.to_dict()
    d["levels"]["156Gd"] = copy.deepcopy(d["levels"]["158Dy"])
    with pytest.raises(SetupError) as err:
        Experiment.from_dict(d)
    assert any("is not the residue 156Gd" in p for p in err.value.problems)


# -- the seam for fusion.py and evaporation.py ------------------------------------------------------------------


def test_a_table_from_another_code_fills_the_setup_without_changing_its_shape(exp):
    t = ChannelTable.from_cross_sections("12C", "150Nd", {"4n": 480.0, "3n": 70.0, ("n",) * 5: 110.0, "p3n": 0.0},
                                         "physim fusion + evaporation")
    assert t.compound == "162Dy"
    assert t.labels == ["4n", "3n", "5n"]
    assert [ch.residue for ch in t.channels] == ["158Dy", "159Dy", "157Dy"]
    s = t.scaled_to(1000.0)
    assert s.total_mb == pytest.approx(1000.0)
    assert np.allclose(s.shares(), t.shares())
    d = exp.to_dict()
    d["reaction"]["channels"] = s.to_setup()
    d["levels"].pop("156Gd")  # a scheme belongs to a residue of the table
    e = Experiment.from_dict(d)
    back = ChannelTable.from_experiment(e)
    assert [ch.cross_section_mb for ch in back.channels] == pytest.approx([ch.cross_section_mb for ch in s.channels])
    assert {ch.source for ch in back.channels} == {"physim fusion + evaporation"}
    assert isinstance(e.fusion.channels[0], ExitChannel)


def test_residue_table_numbers(exp):
    t = residue_table(exp)
    assert t["compound"] == "162Dy"
    rows = {r["channel"]: r for r in t["rows"]}
    four = rows["4n"]
    assert four["rate_per_s"] == pytest.approx(520e-27 * t["luminosity"])
    assert four["lines"] and four["lines"][0]["detected_per_s"] > 0
    assert four["scheme"] == "setup"
    assert all(r["open"] for r in t["rows"])
    assert four["beam_time_s"] == pytest.approx(10_000 / four["lines"][0]["detected_per_s"])
    # Below the threshold a channel is closed and the table says so.
    d = exp.to_dict()
    d["reaction"]["channels"].append({"particles": "7n", "cross_section": "1 mb", "source": "placeholder"})
    t2 = residue_table(Experiment.from_dict(d))
    seven = [r for r in t2["rows"] if r["channel"] == "7n"][0]
    assert not seven["open"]
    assert any("7n" in n and "needs E*" in n for n in t2["notes"])
    assert any("placeholder" in n for n in t2["notes"])


# -- the planner, the app's model and the MCP tools ---------------------------------------------------------------


def test_planner_tabs_work_for_fusion_evaporation():
    p = Planner.example(EXAMPLE)
    assert not p.problems
    plan = p.plan()
    assert plan["residues"]["rows"][1]["residue"] == "158Dy"
    for tab in ("geometry", "kinematics", "rates", "energy_loss", "trajectories", "gamma"):
        p.tab(tab)
    assert p.spectra(events=20_000)["spectra"]
    assert p.gamma_spectra()["available"] is False
    g = p.residue_gamma_spectra(events=10_000)
    assert g["available"] and set(g["spectra"]) == {"Ge90", "Ge37", "Ge79", "Ge101", "Ge143"}
    assert g["spectra"]["Ge90"]["corrected"].sum() > 0


def test_planner_switches_to_fusion_evaporation_and_edits_channels():
    p = Planner.example("oxygen_on_lead_array")
    p.set("beam", "energy", "100 MeV")
    assert not p.problems
    p.set("reaction", "type", "fusion-evaporation")
    chans = p.draft["reaction"]["channels"]
    assert len(chans) == 1 and chans[0]["source"] == "placeholder" and chans[0]["particles"].endswith("n")
    assert not p.problems, p.problems
    assert any("placeholder" in w.text for w in p.warnings())
    assert p.set("reaction channel 1", "cross_section", "250 mb")
    assert p.experiment.fusion.channels[0].cross_section_mb == 250.0
    assert p.add_channel(particles="p2n", cross_section="20 mb")
    assert len(p.experiment.fusion.channels) == 2
    assert p.remove_channel("p2n")
    q = Planner.example(EXAMPLE)
    assert q.remove_channel("a2n") and "156Gd" not in q.draft["levels"] and not q.problems
    q.set("reaction", "type", "elastic")
    assert set(q.draft.get("levels", {})) == set() and not q.problems
    q.set("reaction", "type", "fusion-evaporation")
    assert "158Dy" in q.draft["levels"] and not q.problems
    p.set("reaction", "type", "elastic")
    assert p.draft["reaction"] == {"type": "elastic"}
    p.set("reaction", "type", "fusion-evaporation")
    assert p.draft["reaction"]["channels"][0]["cross_section"] == "250 mb"


def test_a_short_run_of_a_fusion_evaporation_setup(tmp_path):
    p = Planner.example(EXAMPLE)
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="5 s")
    assert p.run_status()["kind"] == "beam"
    assert p.spectra()["spectra"]["Monitor"]["counts"].sum() > 0


def test_the_mcp_tools_handle_fusion_evaporation(tmp_path):
    pytest.importorskip("mcp.server.mcpserver", reason="the MCP extra is not installed")
    from physim.mcp.session import Session
    from physim.mcp.tools_setup import SetupTools

    tools = SetupTools(Session(tmp_path))
    out = tools.new_experiment("fe", from_example=EXAMPLE)
    assert out["reaction"]["type"] == "fusion-evaporation"
    assert any(line.startswith("4n -> 158Dy") for line in out["residues"])
    out = tools.set_value("reaction channel 2", "cross_section", "400 mb")
    assert out["ok"]
    assert 'cross_section = "400 mb"' in tools.get_setup()["toml"]


def test_the_app_page_shows_the_residues_and_edits_the_channels():
    """The planner page on the example: the Plan tab's table of residues, and the channel table in the Target and
    reaction card, where a cross section is typed and applied."""
    pytest.importorskip("nicegui")
    import importlib

    app_tests = importlib.import_module("test_nuclear_app")

    async def scenario(user, p, on_loop):
        await user.open("/")
        await user.should_see("Beam time needed", retries=600)
        await user.should_see("Strongest lines", retries=600)
        assert "162Dy* at E* = " in app_tests._texts(user)
        await app_tests._open_card(user, "Target and reaction", "Cross section")
        boxes = app_tests._fields(user, "Cross section")
        assert sorted(str(b.value) for b in boxes) == ["10 mb", "110 mb", "25 mb", "35 mb", "520 mb"]
        four = next(b for b in boxes if str(b.value) == "520 mb")
        app_tests._type(user, four, "400 mb")
        await app_tests.asyncio_sleep(0.1)
        assert "1 change not applied" in app_tests._texts(user)
        await app_tests.asyncio_sleep(2.0)
        app_tests._press(user, "Apply")
        for _ in range(100):
            await app_tests.asyncio_sleep(0.05)
            if p.experiment.fusion.channels[1].cross_section_mb == 400.0:
                break
        assert p.experiment.fusion.channels[1].cross_section_mb == 400.0

    app_tests._page(EXAMPLE, scenario)
