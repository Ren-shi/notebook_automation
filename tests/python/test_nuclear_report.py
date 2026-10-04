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
    for heading in ("Warnings", "Setup", "Detectors", "Kinematics", "Expected peaks", "Energy loss",
                    "Simulated spectra", "Models and their validation", "Data and references"):
        assert f"<h2>{heading}" in text, heading
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
    assert {"report.html", "detectors.csv", "setup.toml", "figures/spectra.pdf"} <= set(names)
