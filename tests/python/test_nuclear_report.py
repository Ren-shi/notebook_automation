"""The beam-time report and data exports (backlog item 42): reproducible numbers, checked against hand calculations."""

import csv
import math

import pytest

from physim.nuclear import Experiment, data, report
from physim.nuclear.rates import stack
from physim.nuclear.rutherford import Rutherford
from physim.nuclear.stopping import Stopping

N_A = 6.02214076e23
E_CHARGE = 1.602176634e-19


@pytest.fixture(scope="module")
def alpha():
    return report.build(Experiment.example("alpha_on_gold"), seed=3, events=100_000, validation=False)


def test_same_setup_file_and_seed_give_the_same_numbers(tmp_path):
    exp = Experiment.example("oxygen_on_lead_array")
    path = tmp_path / "setup.toml"
    exp.save(path)
    a = report.build(Experiment.load(path), seed=7, events=50_000, validation=False)
    b = report.build(Experiment.load(path), seed=7, events=50_000, validation=False)
    assert a.numbers() == b.numbers()
    c = report.build(Experiment.load(path), seed=8, events=50_000, validation=False)
    assert [r["mc_rate_per_s"] for r in c.detectors] != [r["mc_rate_per_s"] for r in a.detectors]
    # Everything but the Monte Carlo columns is independent of the seed.
    strip = lambda rows: [{k: v for k, v in r.items() if not k.startswith("mc_")} for r in rows]  # noqa: E731
    assert strip(c.detectors) == strip(a.detectors)


def test_alpha_on_gold_against_a_hand_calculation(alpha):
    """1 pnA of 5.5 MeV α on 0.5 mg/cm² Au; a 2.5 mm disc at 100 mm and 45°."""
    pps = 1e-9 / E_CHARGE
    n = 0.5e-3 * N_A / data.nuclide("197Au").atomic_mass_u
    omega = 2 * math.pi * (1 - 100 / math.hypot(100, 2.5))  # sr
    loss = Stopping("4He", "Au").stopping_power(5.5) * 0.25  # to mid-target, MeV
    e_mid = 5.5 - loss
    # Rutherford in the lab at 45° by hand: CM energy, CM angle and Jacobian from non-relativistic kinematics.
    m1, m2 = 4.0026, 196.9666
    e_cm = e_mid * m2 / (m1 + m2)
    x = m1 / m2
    th = math.radians(45.0)
    th_cm = th + math.asin(x * math.sin(th))
    jac = (1 + x * x + 2 * x * math.cos(th_cm)) ** 1.5 / abs(1 + x * math.cos(th_cm))
    sigma_cm = (2 * 79 * 1.43996 / (4 * e_cm)) ** 2 / math.sin(th_cm / 2) ** 4 * 10  # mb/sr
    hand_rate = pps * n * sigma_cm * jac * omega * 1e-27
    row = {r["detector"]: r for r in alpha.detectors}["A45"]
    assert row["solid_angle_msr"] == pytest.approx(omega * 1e3, rel=1e-9)
    assert row["dsigma_domega_lab_mb_sr"] == pytest.approx(sigma_cm * jac, rel=2e-3)
    # Counted rate = the scattered α (the gold recoils, ~0.2 MeV at 45°, stop in the foil or fall below threshold).
    assert row["rate_per_s"] == pytest.approx(hand_rate, rel=0.01)
    assert row["rate_all_per_s"] > row["rate_per_s"]
    assert row["counts_in_run"] == pytest.approx(row["rate_per_s"] * 7200)
    assert row["beam_time_s"] == pytest.approx(10_000 / row["rate_per_s"])
    assert row["mc_rate_per_s"] == pytest.approx(row["rate_per_s"], abs=4 * row["mc_rate_error_per_s"])
    # Energy loss through the whole foil: S × t (S changes by under 1% across it).
    layer = alpha.energy_loss["layers"][0]
    assert layer["loss_mev"] == pytest.approx(Stopping("4He", "Au").stopping_power(5.45) * 0.5, rel=0.01)
    # Peak at 45°: kinematic factor at mid-target, then out through the other half at 1/cos 45°.
    k = ((m1 * math.cos(th) + math.sqrt(m2**2 - (m1 * math.sin(th)) ** 2)) / (m1 + m2)) ** 2
    e_out = k * e_mid - Stopping("4He", "Au").stopping_power(k * e_mid) * 0.25 / math.cos(th)
    peak = [p for p in alpha.peaks if p["detector"] == "A45" and p["particle"] == "ejectile"][0]
    assert peak["mean_MeV"] == pytest.approx(e_out, abs=2e-3)


def test_oxygen_on_lead_against_a_hand_calculation():
    rep = report.build(Experiment.example("oxygen_on_lead_array"), seed=1, events=50_000, validation=False)
    # 10 enA of 16O(6+) = 10e-9 / (6 e) particles per second on 200 µg/cm² of 208Pb.
    pps = 10e-9 / (6 * E_CHARGE)
    n = 0.2e-3 * N_A / data.nuclide("208Pb").atomic_mass_u
    row = {r["detector"]: r for r in rep.detectors}["DSSD3"]
    # A 50 × 50 mm detector at 150 mm: σ at its mean angle × Ω, within 3% (σ curves across ±9.5°).
    hand = pps * n * row["dsigma_domega_lab_mb_sr"] * row["solid_angle_msr"] * 1e-3 * 1e-27
    assert row["rate_per_s"] == pytest.approx(hand, rel=0.03)
    target, backing = rep.energy_loss["layers"]
    assert target["loss_mev"] == pytest.approx(Stopping("16O", "208Pb").stopping_power(63.8) * 0.2, rel=0.01)
    assert backing["energy_in_mev"] == pytest.approx(target["energy_out_mev"])
    # The kinematics table: elastic 16O on 208Pb at 105°.
    i = list(rep.kinematics["theta_lab_deg"]).index(105.0)
    from physim.nuclear.kinematics import TwoBody

    assert rep.kinematics["E_16O_on_208Pb_(target)_MeV"][i] == pytest.approx(
        TwoBody("16O", "208Pb", 64.0).at_lab(105.0)[0].energy, rel=1e-6)


def test_exports(tmp_path, alpha):
    paths = alpha.write(tmp_path, figures=False)
    names = {p.name for p in paths}
    assert {"report.html", "detectors.csv", "peaks.csv", "kinematics.csv", "energy_loss.csv", "strips.csv",
            "setup.toml"} <= names
    with open(tmp_path / "detectors.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 6 and "solid_angle_msr" in rows[0] and "rate_per_s" in rows[0]
    assert Experiment.load(tmp_path / "setup.toml").to_dict() == alpha.experiment.to_dict()
    text = (tmp_path / "report.html").read_text(encoding="utf-8")
    for heading in ("Warnings", "Detectors", "Kinematics", "Expected peaks", "Energy loss",
                    "Simulated spectra", "Models and their validation", "Data and references"):
        assert f"<h3>{heading}" in text, heading
    assert '<h2 id="setup">A. The setup</h2>' in text and "No run yet" in text  # the one report, without runs
    assert text.index('<h2 id="result">1. What this setup gives</h2>') < text.index('<h2 id="appendices"')
    for heading in ("Warnings", "Setup", "Detectors", "Kinematics", "Expected peaks", "Energy loss",
                    "Simulated spectra", "Models and their validation", "Data and references"):
        assert f"<h2>{heading}" in alpha.to_html(), heading  # the plan alone, as a page
    assert "pile-up" in text and "AME2020" in text and "@media print" in text
    assert text.count("data:image/png;base64,") == 4
    assert ">nan<" not in text.lower()


def test_figures_and_command_line(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    rc = report.main(["alpha_on_gold", "-o", str(tmp_path / "out"), "--events", "20000"])
    assert rc == 0
    for name in ("geometry", "coverage", "kinematics", "spectra"):
        assert (tmp_path / "out" / "figures" / f"{name}.png").stat().st_size > 10_000
        assert (tmp_path / "out" / "figures" / f"{name}.pdf").exists()


def _png_width(path):
    return int.from_bytes(path.read_bytes()[16:20], "big")


def test_figures_in_a_journal_style(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    # The default is Physical Review, one column: 86 mm at 600 dpi.
    rep = report.build(Experiment.example("alpha_on_gold"), events=20_000, validation=False)
    rep.write_figures(tmp_path / "prc")
    assert _png_width(tmp_path / "prc" / "kinematics.png") == pytest.approx(86 / 25.4 * 600, abs=2)
    figs = rep.figures()
    assert figs["kinematics"].axes[0].xaxis.label.get_fontfamily() == ["serif"]
    assert len(figs["spectra"].axes) == len(rep.experiment.detectors)  # one panel per detector

    rc = report.main(["alpha_on_gold", "-o", str(tmp_path / "nature"), "--events", "20000", "--journal", "nature",
                      "--width", "double"])
    assert rc == 0
    assert _png_width(tmp_path / "nature" / "figures" / "coverage.png") == pytest.approx(183 / 25.4 * 600, abs=2)
    with pytest.raises(ValueError, match="has no 'middle' width"):
        report.build(Experiment.example("alpha_on_gold"), journal="nature", width="middle")


def test_coulomb_excitation_in_the_report(tmp_path):
    rep = report.build(Experiment.example("coulex_ni58"), seed=1, events=20_000, validation=False)
    x = rep.excitation
    assert x["state"]["energy_kev"] == pytest.approx(1454.0)
    assert len(x["doppler"]) == 3 * 4
    paths = {p.name for p in rep.write(tmp_path, figures=False, root=False)}
    assert "gamma.csv" in paths
    text = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "<h3>Coulomb excitation" in text and "B(E2↑)" in text and "Doppler shift" in text
    assert "<h2>Coulomb excitation" in rep.to_html()
    # The inelastic peak sits below the elastic one by about the kinematic difference at the CD's angles: less
    # than E* for backscattered 16O, as the recoil takes part of Q. But less than that difference at mid-target:
    # P doubles from the back of the target to the front (30 → 27.9 MeV), so the excitations happen nearer the
    # front, where the beam is faster and the scattered 16O crosses less of the target on its way out (its
    # measured energy falls 4.1 MeV from front to back). That lifts the inelastic peak by about 0.24 MeV: the
    # difference is 0.63 MeV against 0.88 at mid-target (the event generator, cut into 16 slices, gives 0.62).
    from physim.nuclear.kinematics import TwoBody

    peaks = {p["channel"]: p for p in rep.peaks if p["detector"] == "CD" and p["particle"] == "ejectile"}
    el, inel = peaks["58Ni (target)"], peaks["58Ni (target), target excited to 1454 keV"]
    theta = rep.detectors[0]["theta_mean_deg"]
    e = 29.0  # MeV, mid-target
    kin = (TwoBody("16O", "58Ni", e).at_lab(theta)[0].energy
           - TwoBody("16O", "58Ni", e, excitation_mev=1.454).at_lab(theta)[0].energy)
    assert 0.6 * kin < el["mean_MeV"] - inel["mean_MeV"] < 0.85 * kin
    assert kin < 1.454
    assert report.build(Experiment.example("alpha_on_gold"), events=10_000, validation=False).excitation is None


def test_validation_status_in_the_report():
    rep = report.build(Experiment.example("alpha_on_gold"), events=10_000)
    rows = {r["capability"]: r for r in rep.register}
    assert rows["Rutherford cross section"]["literature"] == "pass"
    assert rows["Detector solid angles and response"]["tool"] == "none"
    assert rows["Count rates and beam time"]["tool"] in ("pending", "unchecked")


def test_root_file_reads_back_unchanged(tmp_path):
    pytest.importorskip("uproot")
    import numpy as np

    from physim.nuclear.events import simulate
    from physim.nuclear.rootio import read_root, write_root

    exp = Experiment.example("oxygen_on_lead_array")
    ev = simulate(exp, 100_000, seed=5)
    path = write_root(ev, tmp_path / "events.root", exp)
    back = read_root(path)
    assert set(back["events"]) == set(ev.columns)
    for k, v in ev.columns.items():
        assert np.array_equal(back["events"][k], np.asarray(v)), k
    assert Experiment.from_toml(back["setup"]).to_dict() == exp.to_dict()
    assert back["info"]["seed"] == 5 and back["info"]["n_events"] == 100_000
    for name in ev.detectors:
        counts, errors, edges = back["spectra"][name]
        assert counts.sum() == pytest.approx(ev.counts(name), rel=1e-9)
        w = ev["weight"][ev.select(name)] * ev.beam_time_s
        # Errors are the Monte Carlo statistical errors, √Σw² per bin.
        assert np.sqrt((errors**2).sum()) == pytest.approx(np.sqrt((w**2).sum()), rel=1e-9)
        assert edges[0] == 0.0 and len(edges) == 401
    # A Gaussian fit to the strongest peak recovers the analytic mean (as ROOT's Fit("gaus") would).
    from physim.nuclear.rates import Rates

    pk = Rates(exp).peaks("DSSD3")[0]
    counts, errors, edges = back["spectra"]["DSSD3"]
    centres = 0.5 * (edges[1:] + edges[:-1])
    sel = np.abs(centres - pk.mean) < 3 * pk.sigma
    mean = np.average(centres[sel], weights=counts[sel])
    assert mean == pytest.approx(pk.mean, abs=0.1)


def test_report_writes_the_root_file(tmp_path, alpha):
    pytest.importorskip("uproot")
    paths = alpha.write(tmp_path, figures=False)
    assert tmp_path / "events.root" in paths
    assert alpha.write(tmp_path / "no", figures=False, root=False)[-1].name == "setup.toml"


def test_the_app_downloads_the_report():
    pytest.importorskip("plotly")
    import io
    import zipfile

    from physim.nuclear.app import report_zip
    from physim.nuclear.planner import Planner

    data = report_zip(Planner.example("alpha_on_gold"), seed=1, events=20_000)
    names = zipfile.ZipFile(io.BytesIO(data)).namelist()
    assert {"report.html", "tables/detectors.csv", "setup.toml", "figures/spectra.pdf"} <= set(names)


# -- fusion-evaporation: the checks that apply and the limits (planning items F14, F15) ------------------------


def _fusion_variant(beam=None, energy=None, target=None):
    """The fusion_dy158 example with another beam, energy or target (one 2n channel when the system changes)."""
    d = Experiment.example("fusion_dy158").to_dict()
    if beam:
        d["beam"]["nuclide"] = beam
    if energy:
        d["beam"]["energy"] = energy
    if target:
        d["target"]["material"] = target
    if beam or target:
        b, t = data.nuclide(d["beam"]["nuclide"]), data.nuclide(d["target"]["material"])
        d["reaction"]["channels"] = [{"residue": data.nuclide((b.Z + t.Z, b.A + t.A - 2)).name, "particles": "2n",
                                      "cross_section": "100 mb"}]
        d.pop("levels", None)
    return Experiment.from_dict(d)


@pytest.fixture(scope="module")
def fusion_report():
    return report.build(Experiment.example("fusion_dy158"), seed=1, events=20_000)


def test_fusion_checks_are_picked_by_capability(fusion_report):
    from physim.nuclear import validation

    caps = report.fusion_capabilities(validation.CAPABILITIES)
    assert set(report.FUSION_MODELS) <= set(caps)
    checks = fusion_report.fusion["checks"]
    # Every check of those capabilities, and only those: one added to the suite appears without a change here.
    assert {c["name"] for c in checks} == {c.name for c in validation.CHECKS if c.capability in caps}
    rows = {r["capability"]: r for r in fusion_report.register}
    assert set(caps) <= set(rows)
    known = report.known_failures()
    for c in checks:
        assert c["status"] in ("pass", "pending", "known failure", "fail", "not run here")
        if c["status"] in ("known failure", "fail"):
            assert (c["status"] == "known failure") == (c["name"] in known)
        if c["kind"] == "tool" and c["checker"] == "PACE4":
            assert c["version"] == "4.36.0"
    pace = {c["name"]: c for c in checks}["fusion_vs_pace4"]
    # The PACE4 system nearest the run is the run's own; the suite's numbers for it are quoted, not the others'.
    assert pace["nearest"]["system"] == "12C + 150Nd" and pace["nearest"]["same"]
    assert "12C + 150Nd" in pace["nearest"]["quote"] and "48Ca" not in pace["nearest"]["quote"]


def test_fusion_section_in_the_report(fusion_report, tmp_path):
    import html

    page = fusion_report.to_html()
    section = page[page.index("<h2>Fusion and evaporation"):page.index("<h2>Energy loss")]
    fus = section[section.index("Fusion cross section"):section.index("<h3>Residues")]
    assert "σ_fus at the beam energy" in fus and "Checked against" in fus and "PACE4 4.36.0" in fus
    known = {c["name"] for c in fusion_report.fusion["checks"] if c["status"] == "known failure"}
    if "fusion_vs_pace4" in known:  # said plainly beside σ_fus, not only in the table at the end
        assert "Known failure" in fus
    assert "No coupled-channels fusion below the barrier" in fus  # attached, though not reached at 1.24 V_B
    res = section[section.index("<h3>Residues"):]
    assert "illustrative" in res and "The evaporation at high spin" in res
    assert "⚠" not in section  # the example is inside every limit
    models = page[page.index("<h2>Models and their validation"):page.index("<h2>Data and references")]
    assert "Fusion-evaporation checks" in models and "Limits of the fusion-evaporation models" in models
    for lim in report.LIMITS:
        assert html.escape(lim["title"]) in models
    # The app's report holds the same section in its plan.
    fusion_report.write(tmp_path, figures=False, root=False)
    text = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "<h3>Fusion and evaporation</h3>" in text and "Fusion-evaporation checks" in text
    assert ">nan<" not in text.lower()
    assert fusion_report.numbers()["fusion"]["summary"]["sigma_mb"] == pytest.approx(692.5, rel=0.01)


def test_fusion_limits_flag_the_runs_they_cover():
    from physim.nuclear.planner import Planner

    def flags(exp):
        p = Planner(exp)
        return {lim["key"]: lim for lim in report._fusion_data(p, exp)["limits"] if lim["flag"]}

    assert flags(_fusion_variant()) == {}
    below = flags(_fusion_variant(energy="45 MeV"))  # E_cm 41.7 MeV, V_B 44.7 MeV
    assert set(below) == {"coupled_channels"} and "lower bound" in below["coupled_channels"]["flag"]
    assert "near the fusion barrier" in flags(_fusion_variant(energy="52 MeV"))["coupled_channels"]["flag"]
    assert "inside the target" in flags(_fusion_variant(energy="49 MeV"))["coupled_channels"]["flag"]
    fast = flags(_fusion_variant(beam="4He", energy="50 MeV", target="165Ho"))  # 12.5 MeV/u
    assert set(fast) == {"incomplete_fusion"}
    heavy = flags(_fusion_variant(beam="48Ca", energy="200 MeV", target="124Sn"))
    assert set(heavy) == {"high_spin"} and heavy["high_spin"]["level"] == "note"  # the table is the setup's own
    qf = flags(_fusion_variant(beam="32S", energy="180 MeV", target="208Pb"))  # Z1Z2 = 1312
    assert qf["quasi_fission"]["level"] == "note"
    # The warning stands beside the number it bears on.
    rep = report.build(_fusion_variant(energy="45 MeV"), events=5_000, validation=False)
    page = rep.to_html()
    fus = page[page.index("<h3>Fusion cross section"):page.index("<h3>Residues")]
    assert "⚠" in fus and "lower bound" in fus
    # With a neutron wall: the selection matrix, and its simplifications beside it.
    p = Planner.example("fusion_dy158")
    p.add_neutron_wall()
    f = report._fusion_data(p, p.experiment)
    sel = {lim["key"]: lim for lim in f["limits"]}["selection"]
    assert sel["flag"] and "α" in sel["flag"]  # p3n and a2n: the particle identification matters
    rep = report.build(p.experiment, events=5_000, validation=False)
    assert "<h3>Channel selection" in rep.to_html() and "cross-talk" in rep.to_html()


def test_new_checks_appear_with_the_nearest_system():
    """A check the suite adds (a GEMINI++ comparison, here made up) shows with its version and the suite's numbers
    for the system nearest the run; a known failure is marked as such; other capabilities stay out."""
    from physim.nuclear.validation import Result

    gem = Result("Evaporation Monte Carlo (Hauser–Feshbach)", "tool", "evaporation_vs_gemini",
                 "GEMINI++ 3.3.0, ⁴He + ¹⁶⁵Ho, ¹²C + ¹⁵⁰Nd, ¹⁶O + ¹⁴⁴Sm", "30%", "fail", 0.4, 9,
                 "shares GEMINI++/physim: 4He + 165Ho: 165Tm 0.84/0.85; 12C + 150Nd: 157Dy 0.76/0.65, "
                 "158Dy 0.11/0.14; 16O + 144Sm: 157Yb 0.37/0.22")
    other = Result("Rutherford cross section", "tool", "x", "LISE++ 17", "1%", "pass", 0.1)
    measured = Result("Fusion cross section (Bass barrier, Wong)", "literature", "fusion_vs_data",
                      "Measured σ_ER: 16O + 142Nd, 16O + 150Nd (Visakh 2021)", "20%", "pass", 0.5, 30,
                      "worst per system: 16O + 142Nd 10%, 16O + 150Nd 17%")
    rows = report.fusion_checks([gem, other, measured], Experiment.example("fusion_dy158"),
                                known={"evaporation_vs_gemini"})
    assert [r["name"] for r in rows] == ["fusion_vs_data", "evaporation_vs_gemini"]
    g = rows[1]
    assert (g["checker"], g["version"], g["status"]) == ("GEMINI++", "3.3.0", "known failure")
    assert g["nearest"]["system"] == "12C + 150Nd" and g["nearest"]["same"]
    assert g["nearest"]["quote"] == "shares GEMINI++/physim: 12C + 150Nd: 157Dy 0.76/0.65, 158Dy 0.11/0.14"
    d = rows[0]  # no 12C data: the nearest measured system in beam and compound nucleus, its number only
    assert d["nearest"]["system"] == "16O + 150Nd" and not d["nearest"]["same"]
    assert d["nearest"]["quote"] == "worst per system: 16O + 150Nd 17%"


def test_known_failures_follow_the_suite():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location("_validation_tests",
                                                  Path(__file__).with_name("test_nuclear_validation.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert report.known_failures() == frozenset(mod.KNOWN_FAILURES)


def test_other_reactions_have_no_fusion_section():
    rep = report.build(Experiment.example("coulex_ni58"), events=5_000, validation=False)
    assert rep.fusion is None and "fusion" not in rep.numbers()
    parts = rep.parts(dpi=50)
    assert parts["fusion"] == "" and parts["fusion_plain"] == ""
    assert "Fusion-evaporation checks" not in parts["models"]
