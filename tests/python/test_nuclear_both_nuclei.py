"""Both nuclei of the collision are in the events: with a level scheme of the nucleus the ``[reaction]`` does not
name, its γ rays are simulated too; each nucleus's lines are sharp with its own Doppler correction; a search finds
the peaks, each is fitted and integrated, and B(E2) follows from the areas, to all orders. One particle stands for
one collision: elastic scattering is what the excitation leaves."""

import copy
import json

import numpy as np
import pytest

from physim.nuclear import Experiment, logbook
from physim.nuclear.analysis import Analysis, Result, Settings, find_peaks, fit_peak, identify
from physim.nuclear.events import simulate
from physim.nuclear.gamma import correlation_table, emitters, level_names, other_excitation, transition_lines
from physim.nuclear.gamma_events import _own_shares, crystals_of, other_lines, simulate_gammas, spectrum_top
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates
from physim.nuclear.response import FWHM_PER_SIGMA

#: The first transition of ¹⁹⁴Pt and of ⁶²Ni, keV, and their names in the template's schemes.
PT_LINE, NI_LINE = 328.473, 1172.91
PT_FIRST, NI_FIRST = "194Pt 2+₁ → 0+₁", "62Ni 2+₁ → 0+₁"
ELASTIC, EXCITED = "194Pt (target)", "194Pt (target), projectile excited to 1172.91 keV"


@pytest.fixture(scope="module")
def ni62():
    exp = Experiment.example("coulex_ni62_pt194")
    return exp, simulate_gammas(exp, 150_000, seed=3)


def without_target_scheme(exp) -> Experiment:
    d = exp.to_dict()
    del d["levels"]["target"]
    return Experiment.from_dict(d)


def test_the_setup_names_both_nuclei_and_their_lines():
    exp = Experiment.example("coulex_ni62_pt194")
    both = emitters(exp)
    assert [(em["nucleus"], em["nuclide"], em["role"], em["correction"]) for em in both] == [
        (0, "62Ni", "projectile", "projectile"), (1, "194Pt", "target", "recoil")]
    sol = other_excitation(exp)
    assert sol.excite == "target" and sol.scheme.nuclide == "194Pt"
    # The scheme holds the positive-parity levels up to the 6⁺ state: what E2 excitation reaches.
    assert level_names(sol.scheme.levels) == ["0+₁", "2+₁", "2+₂", "4+₁", "4+₂", "0+₂", "6+"]
    # Near the barrier the platinum is excited in almost every backward collision: no first-order matter.
    back = sol.near(160.0)
    assert 0.9 < sum(back.direct.values()) < 1.0 and back.gamma_yield(1, 0) > 0.8
    assert back.direct[6] > 0.02  # the 6⁺ state, three steps up
    lines = transition_lines(exp)
    pt = [ln for ln in lines if ln["nucleus"] == 1]
    assert {328.5, 293.6, 622.0, 482.8, 600.5} <= {round(ln["energy_kev"], 1) for ln in pt}
    assert all(ln["correction"] == "recoil" for ln in pt)
    first = next(ln for ln in pt if abs(ln["energy_kev"] - PT_LINE) < 0.1)
    assert first["label"] == PT_FIRST and (first["initial"], first["final"]) == (1, 0)
    assert any(ln["label"] == NI_FIRST and ln["correction"] == "projectile" for ln in lines)
    # Its γ rays count in the singles (and so in the random coincidences), line by line.
    rates = dict(other_lines(exp))
    assert max(rates, key=rates.get) == pytest.approx(PT_LINE * 1e-3) and all(v > 0 for v in rates.values())
    assert spectrum_top(exp) >= 1.1 * 2.3  # the top of the nickel scheme still sets the range
    # Without the scheme there is one nucleus, as before.
    alone = without_target_scheme(exp)
    assert other_excitation(alone) is None and len(emitters(alone)) == 1 and other_lines(alone) == []
    assert correlation_table(alone, nucleus="other") == []


def test_one_particle_stands_for_one_collision():
    """Elastic scattering is what the excitation leaves: elastic and excited together are Rutherford scattering,
    in the rates and in the events."""
    exp = Experiment.example("coulex_ni62_pt194")
    d = exp.to_dict()
    d["reaction"] = {"type": "elastic"}
    d.pop("levels", None)
    d["gamma_detectors"] = []
    rutherford = Rates(Experiment.from_dict(d))
    r = Rates(exp)
    ev = simulate(exp, 150_000, seed=4)
    assert ev.channels == [ELASTIC, EXCITED]
    for det in ev.detectors:
        by = {}
        for (channel, _), v in r.by_channel(det).items():
            by[channel] = by.get(channel, 0.0) + v
        # The share excited is large at these angles, and what is left is elastic.
        assert 0.1 < by[EXCITED] / (by[ELASTIC] + by[EXCITED]) < 0.3
        assert by[ELASTIC] + by[EXCITED] == pytest.approx(rutherford.rate(det, counted=False), rel=0.01)
        for channel in (ELASTIC, EXCITED):
            got, err = ev.rate(det, channel=channel)
            want = sum(v for (c, _), v in r.by_channel(det, counted=True).items() if c == channel)
            assert abs(got - want) < 4 * err + 0.01 * want, (det, channel, got, want, err)


def test_the_events_hold_the_other_nucleus_and_leave_the_first_as_it_was(ni62):
    exp, g = ni62
    own, other = g["nucleus"] == 0, g["nucleus"] == 1
    assert own.any() and other.sum() > own.sum()
    # The γ rays of the [reaction] nucleus are those of the setup without the other scheme, one for one.
    alone = simulate_gammas(without_target_scheme(exp), particle_events=g.events, seed=3)
    assert len(alone) == own.sum() and (alone["nucleus"] == 0).all()
    for key in ("particle", "crystal", "measured", "weight", "corrected_projectile", "initial", "final", "cascade"):
        assert np.array_equal(alone[key], g[key][own]), key
    # The other nucleus is excited on any collision, whether or not the first one was: some particles carry the
    # γ rays of both.
    ev = g.events
    elastic, excited = ev.channels.index(ELASTIC), ev.channels.index(EXCITED)
    assert set(np.unique(ev["channel"][g["particle"][other]])) == {elastic, excited}
    assert (ev["channel"][g["particle"][own]] == excited).all()
    assert np.intersect1d(g["particle"][own], g["particle"][other]).size > 0
    assert not np.intersect1d(g["cascade"][own], g["cascade"][other]).size
    # Their energies are the differences of the platinum levels, and the levels are counted in its scheme.
    levels = np.array([lv.energy.value for lv in exp.levels["target"].levels]) * 1e-3
    assert np.allclose(g["energy0"][other], levels[g["initial"][other]] - levels[g["final"][other]])
    assert set(np.unique(g["initial"][other])) <= set(range(1, 7))
    assert "194Pt, the other nucleus" in g.notes[-1] and "can excite both" in g.notes[-1]
    # Left out on request: the extension of a run that was taken without it.
    assert (simulate_gammas(exp, particle_events=ev, seed=3, other=False)["nucleus"] == 0).all()


def _in_line(g, mask, key, line_kev, half=2.5):
    e = 1e3 * g[key][mask & g["counted"]]
    w = g["weight"][mask & g["counted"]]
    return float(w[np.abs(e - line_kev) < half].sum())


def test_each_nucleus_is_sharp_with_its_own_correction(ni62):
    _, g = ni62
    own, other = g["nucleus"] == 0, g["nucleus"] == 1
    # Platinum: sharp when corrected for the recoil, smeared for the projectile and without a correction.
    recoil = _in_line(g, other, "corrected_recoil", PT_LINE)
    assert recoil > 5 * _in_line(g, other, "corrected_projectile", PT_LINE)
    assert recoil > 5 * _in_line(g, other, "measured", PT_LINE)
    # Nickel: the other way round.
    assert _in_line(g, own, "corrected_projectile", NI_LINE, 8.0) > 5 * _in_line(g, own, "corrected_recoil", NI_LINE, 8.0)


def test_only_the_gamma_rays_that_reach_a_crystal_are_drawn():
    """The aimed emission gives the γ rays of the whole-sphere one, in distribution: the same rates in all, per
    nucleus, forward and backward, and the same mean energy, on the same particles."""
    exp = Experiment.example("coulex_ni62_pt194")
    ev = simulate(exp, 12_000, seed=2)
    rows = {}
    for aimed in (True, False):
        out = []
        for seed in (1, 2, 3):
            g = simulate_gammas(exp, particle_events=ev, seed=seed, aimed=aimed, gammas_per_event=40, plain=True)
            row = []
            for nucleus in (0, 1):
                m = g["nucleus"] == nucleus
                w = g["weight"]
                row += [w[m].sum(), w[m & (g["theta"] < 80)].sum(), w[m & (g["theta"] > 100)].sum(),
                        float(np.average(g["energy_lab"][m], weights=w[m]))]
            out.append(row)
        rows[aimed] = np.array(out)
    a, b = rows[True], rows[False]
    err = np.hypot(a.std(axis=0, ddof=1), b.std(axis=0, ddof=1)) / np.sqrt(3)
    ratio = a.mean(axis=0) / b.mean(axis=0)
    assert np.all(np.abs(a.mean(axis=0) - b.mean(axis=0)) < 4 * err + 0.004 * b.mean(axis=0)), (ratio, err)
    assert np.all(np.abs(ratio - 1) < 0.05), ratio


def test_the_peak_search_finds_peaks_it_was_not_told_about():
    rng = np.random.default_rng(5)
    edges = np.arange(0.0, 1.5 + 1e-9, 0.5e-3)
    centres = (edges[:-1] + edges[1:]) / 2

    def fwhm(e):
        return np.full(np.shape(e), 2.5e-3)

    def gauss(centre, sigma, area):
        return area * np.exp(-0.5 * ((centres - centre) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi)) * 0.5e-3

    # A sloping level, a narrow peak at the resolution, a Doppler-broadened one, and a weak one below 5 σ.
    expected = (60 - 30 * centres + gauss(0.3285, 2.5e-3 / FWHM_PER_SIGMA, 4000)
                + gauss(1.1729, 15e-3 / FWHM_PER_SIGMA, 2500) + gauss(0.800, 2.5e-3 / FWHM_PER_SIGMA, 40))
    counts = rng.poisson(expected).astype(float)
    peaks = find_peaks(counts, edges, fwhm)
    found = sorted(round(1e3 * pk["energy"], 0) for pk in peaks)
    assert found == [328.0, 1173.0], peaks
    strong = max(peaks, key=lambda pk: pk["significance"])
    assert abs(1e3 * strong["energy"] - 328.5) < 0.3 and strong["counts"] == pytest.approx(4000, rel=0.1)
    assert all(pk["significance"] >= 5 for pk in peaks)
    assert find_peaks(rng.poisson(np.full(len(centres), 60.0)).astype(float), edges, fwhm) == []  # nothing in noise
    # Naming: a peak is a transition of the nucleus the correction is right for, else at most "near" another's.
    lines = [{"label": "194Pt 2+ → 0+", "nuclide": "194Pt", "energy_kev": 328.5, "correction": "recoil"},
             {"label": "62Ni 2+ → 0+", "nuclide": "62Ni", "energy_kev": 1172.9, "correction": "projectile"}]
    named = identify([dict(pk) for pk in peaks], lines, "recoil")
    by_energy = {round(1e3 * pk["energy"]): pk for pk in named}
    assert by_energy[328]["label"] == "194Pt 2+ → 0+" and by_energy[328]["line_kev"] == 328.5
    assert by_energy[1173]["label"] is None and by_energy[1173]["near"] == "62Ni 2+ → 0+"
    named = identify([dict(pk) for pk in peaks], lines, "projectile")
    assert {pk["label"] for pk in named} == {None, "62Ni 2+ → 0+"}


def test_a_scaled_sample_is_fitted_with_the_variance_of_its_entries():
    """A spectrum of few weighted entries scaled to many counts: taken for Poisson counts, the fit follows the
    scatter of the entries; given their variance, it finds the peak."""
    rng = np.random.default_rng(8)
    edges = np.arange(1.10, 1.25 + 1e-9, 0.4e-3)
    sigma = 6.5e-3
    energies = np.r_[rng.normal(1.1729, sigma, 1500), rng.uniform(1.10, 1.25, 500)]
    weight = rng.exponential(1.0, len(energies)) * 400.0
    counts, _ = np.histogram(energies, bins=edges, weights=weight)
    sumw2, _ = np.histogram(energies, bins=edges, weights=weight**2)
    good = fit_peak(counts, edges, 1.1729, sigma, 4.0, max_sigma=0.035, variance=sumw2)
    assert good.sigma == pytest.approx(sigma, rel=0.15)
    assert good.area == pytest.approx(weight[:1500].sum(), rel=0.08)
    # The uncertainty of the area is still that of the counts: what the beam time would give.
    assert good.area_unc == pytest.approx(np.sqrt(counts[np.abs((edges[:-1] + edges[1:]) / 2 - good.centroid)
                                                          < 3 * good.sigma].sum()), rel=0.3)


def test_both_peaks_are_found_integrated_and_give_b_e2(ni62):
    exp, g = ni62
    a = Analysis(exp, g)
    r = a.run(Settings(particle_gate="all"))
    assert [s["step"] for s in r.steps] == ["particle gate", "Doppler correction", "peak search", "peak fit",
                                            "random coincidences", "yield", "normalisation", "all orders"]
    # The search: each nucleus's transitions are named in its own spectrum and in no other.
    projectile, recoil = r.search["projectile"], r.search["recoil"]
    assert (projectile["nuclide"], recoil["nuclide"]) == ("62Ni", "194Pt")
    assert NI_FIRST in [pk["label"] for pk in projectile["peaks"]]
    assert {PT_FIRST, "194Pt 4+₁ → 2+₁", "194Pt 2+₂ → 2+₁"} <= {pk["label"] for pk in recoil["peaks"]}
    assert not any(pk["nuclide"] == "194Pt" for pk in projectile["peaks"])
    assert not any(pk["nuclide"] == "62Ni" for pk in recoil["peaks"])
    for sp in (projectile, recoil):  # what is left of the other nucleus's lines is no peak, or is said to be smeared
        assert all(pk["label"] or pk["fwhm_kev"] <= 0.03 * pk["energy_kev"] for pk in sp["peaks"])
    assert f"{NI_FIRST} at" in r.steps[2]["text"] and f"{PT_FIRST} at" in r.steps[2]["text"]
    # The two integrated peaks. The weighted sample is fitted with its own variance: the peak keeps its width.
    ni, pt = r.lines
    assert (ni["nuclide"], ni["correction"], pt["nuclide"], pt["correction"]) == ("62Ni", "projectile", "194Pt", "recoil")
    assert 12 < r.fit.fwhm_kev < 22
    assert ni["area"] == pytest.approx(r.area) and ni["yield"] == pytest.approx(r.area / (r.efficiency * r.correlation))
    assert ni["found"] is not None and pt["found"] is not None
    assert abs(pt["found"]["energy_kev"] - PT_LINE) < 1.0 and abs(1e3 * pt["fit"]["centroid"] - PT_LINE) < 1.0
    assert pt["area"] > 0 and pt["area_unc"] > 0 and pt["efficiency"] > ni["efficiency"]
    assert any("62Ni 4+ → 2+₁" in b for b in ni["blends"]) and pt["blends"] == []
    # The platinum peak gives the yield its level scheme does, to all orders.
    assert pt["ratio"] == pytest.approx(1.0, abs=0.06)
    assert 0.5 < pt["expected_per_particle"] < 1.0
    assert any("measured over expected" in n for n in r.notes)
    # B(E2) to all orders: below what first order reads from the same yield (the feeding, and the 4⁺ → 2⁺ line
    # in the window), and the value that was put in. The yield saturates, so the uncertainties are larger.
    first = a.run(Settings(particle_gate="all", order="first"))
    assert first.b_e2fm4 == pytest.approx(r.b_first_order) and "all orders" not in [s["step"] for s in first.steps]
    assert 0.90 < r.b_e2fm4 / r.b_first_order < 0.97
    text = r.steps[-1]["text"]
    assert "62Ni 4+ → 2+₁" in text and "saturate" in text
    assert r.b_e2fm4 == pytest.approx(exp.excitation.b_up_e2fm, rel=0.08)
    assert abs(first.b_e2fm4 / exp.excitation.b_up_e2fm - 1) > abs(r.b_e2fm4 / exp.excitation.b_up_e2fm - 1)
    assert 1.3 < r.statistical / first.statistical < 1.7 and 1.3 < r.budget["efficiency"] / 0.05 < 1.7
    # Normalised to the platinum peak, the particle count drops out: first order reads the Rutherford value over
    # the ratio of measured to expected.
    t = a.run(Settings(particle_gate="all", normalisation="target"))
    assert t.normalisation == "target" and f"measured {PT_FIRST} peak" in t.steps[-2]["text"]
    assert t.b_first_order == pytest.approx(r.b_first_order / pt["ratio"], rel=1e-6)
    assert t.b_e2fm4 == pytest.approx(exp.excitation.b_up_e2fm, rel=0.08)
    assert t.budget["reference B(E2)"] > 0.03 and t.statistical > 0
    # The result keeps all of it as plain data.
    again = Result.from_dict(copy.deepcopy(t.to_dict()))
    assert again.lines[1]["area"] == pytest.approx(pt["area"]) and again.search["recoil"]["peaks"]
    assert again.b_first_order == pytest.approx(t.b_first_order)
    page = logbook.analysis_html(t)
    assert 'alt="The peaks found"' in page and 'alt="The fitted peak of 194Pt"' in page
    assert "The peaks integrated" in page and "measured / expected" in page and "all orders" in page


def test_one_nucleus_is_searched_and_reported_the_same_way():
    exp = Experiment.example("coulex_ni58")
    r = Analysis(exp, simulate_gammas(exp, 150_000, seed=1)).run()
    assert len(r.lines) == 1 and r.lines[0]["correction"] == "recoil" and r.lines[0]["found"] is not None
    assert r.search["recoil"]["nuclide"] == "58Ni" and r.search["projectile"]["nuclide"] is None
    assert any(pk["label"] for pk in r.search["recoil"]["peaks"])
    # One state, first-order events: first order is all there is.
    assert r.b_first_order == r.b_e2fm4 and "all orders" not in [s["step"] for s in r.steps]
    with pytest.raises(ValueError, match="level scheme"):
        Analysis(exp, simulate_gammas(exp, 50_000, seed=1)).run(Settings(normalisation="target"))


def test_the_plan_predicts_the_other_peak():
    p = Planner.example("coulex_ni62_pt194")
    r = p._rates()
    other = r.other_peak()
    assert other["label"] == PT_FIRST and other["energy_kev"] == pytest.approx(PT_LINE)
    assert set(other["coincidence_per_s"]) == {"CD", "DSSD-L", "DSSD-R"}
    # Far more of it than of the nickel peak: the platinum is excited in nearly every backward collision.
    assert other["coincidence_per_s"]["CD"] > 8 * r.rate("CD", what="coincidences")
    plan = p.plan()
    assert plan["other_peak"]["label"] == PT_FIRST
    assert all(d["other_coincidence_per_s"] > 0 for d in plan["detectors"])
    assert Rates(without_target_scheme(p.experiment)).other_peak() is None
    assert Planner.example("coulex_ni58").plan()["other_peak"] is None


def test_a_run_holds_both_nuclei_and_its_analysis_shows_both_peaks(tmp_path):
    p = Planner.example("coulex_ni62_pt194")
    p.create_experiment(root=tmp_path)
    run = p.start_run("beam", duration="1 h", budget="20 s")
    g = run.gammas()
    assert set(np.unique(g["nucleus"])) == {0, 1}
    assert run.summary["nuclei"] == 2 and run.summary["one_particle"] is True and not run.one_nucleus
    with np.load(run.folder / "events.npz") as z:
        assert z["g_nucleus"].dtype == np.int8 and len(z["g_nucleus"]) == len(z["g_particle"])
    # Each particle of a run counts once: so does each γ ray, whichever nucleus sent it.
    assert np.ptp(g["weight"]) == 0
    # The coincidences in the peak that the plan predicts are the [reaction] nucleus's own; the other nucleus's
    # peak is counted beside them, and both are set against the plan.
    own_peak = sum(run.summary["peak_coincidences"].values())
    assert 0 < own_peak <= int((g["nucleus"] == 0).sum())
    assert sum(run.summary["other_peak_coincidences"].values()) > 5 * own_peak
    rows = {(row["quantity"], row["what"]): row for row in p.compare()["rows"]}
    theirs = rows[(f"Coincidences in the {PT_FIRST} peak", "CD")]
    assert theirs["agree"] and theirs["measured"] == pytest.approx(theirs["predicted"], rel=0.1)
    r = p.analyse({"normalisation": "target", "particle_gate": "all"})
    assert r["available"], r
    assert [ln["nuclide"] for ln in r["lines"]] == ["62Ni", "194Pt"]
    assert r["lines"][1]["ratio"] == pytest.approx(1.0, abs=0.1)
    assert r["b_e2fm4"] == pytest.approx(p.experiment.excitation.b_up_e2fm, rel=0.4)
    assert "particles in the same rings" in r["steps"][0]["text"]
    kept = p.run_analysis()
    assert kept.lines[1]["label"] == PT_FIRST and kept.search["recoil"]["peaks"]
    keys = [x["key"] for x in p.explanations(result=kept)]
    assert "b_e2" in keys and "yield" in keys

    # A run taken before the events held both nuclei: no column, no word in the summary, one nucleus's γ rays.
    with np.load(run.folder / "events.npz") as z:
        keep = z["g_nucleus"] == 0
        arrays = {k: (z[k][keep] if k.startswith("g_") else z[k]) for k in z.files if not k.endswith("_nucleus")}
    np.savez_compressed(run.folder / "events.npz", **arrays)
    summary = json.loads((run.folder / "summary.json").read_text(encoding="utf-8"))
    for key in ("nuclei", "other_peak_coincidences", "one_particle"):
        summary.pop(key, None)
    summary["gammas"] = int(keep.sum())
    (run.folder / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    p.load_run(1)
    old = p.run
    assert old.one_nucleus and (old.gammas()["nucleus"] == 0).all()
    assert [(s["number"], s["one_nucleus"]) for s in p.runs()] == [(1, True)]
    alone = Analysis(old.experiment, old.gammas(), run=old).run(Settings(particle_gate="all"))
    assert len(alone.lines) == 1 and alone.search["recoil"]["nuclide"] == "194Pt"
    assert "elastic" not in alone.steps[0]["text"]
    # Extended, it stays a run of one nucleus (its two parts would not go together); taken again, it has both,
    # with the same setup, seed and beam time, and the old run stays.
    n_gammas = len(old.gammas())
    extended = p.extend_run("1 h")
    assert extended.summary["nuclei"] == 1 and extended.one_nucleus and len(extended.gammas()) == n_gammas
    new = p.retake_run()
    assert new.number == 2 and new.summary["retake_of"] == 1 and new.summary["seed"] == old.summary["seed"]
    assert new.real_s == pytest.approx(old.real_s) and new.summary["nuclei"] == 2 and not new.one_nucleus
    assert set(np.unique(new.gammas()["nucleus"])) == {0, 1}
    assert [(s["number"], s["one_nucleus"]) for s in p.runs()] == [(1, True), (2, False)]


def test_the_clover_template_is_the_same_experiment_seen_by_four_clovers():
    exp = Experiment.example("coulex_ni62_pt194_clovers")
    single = Experiment.example("coulex_ni62_pt194")
    # Four clovers of four crystals, with add-back and a BGO shield each; the rest of the setup is the template's.
    assert [g.crystals for g in exp.gamma_detectors] == [4, 4, 4, 4]
    assert all(g.addback and g.shield == "BGO" for g in exp.gamma_detectors)
    assert (exp.beam, exp.target, exp.detectors, exp.levels) == (single.beam, single.target, single.detectors,
                                                                  single.levels)
    p, q = Planner(exp), Planner(single)
    # Nothing stands in the beam or overlaps: the checks hold the same note as the template's, and no warning.
    assert not [w for w in p.warnings() if w.level in ("error", "warning")]
    # The particle rates are the template's; the γ-ray efficiency is about twice (1.71% against 0.83%), so the
    # beam time for the counts wanted is well under half.
    assert p.plan()["detectors"][0]["rate_per_s"] == pytest.approx(q.plan()["detectors"][0]["rate_per_s"])
    assert p.selection()["gamma_efficiency"] == pytest.approx(0.0171, rel=0.03)
    assert p.selection()["gamma_efficiency"] > 1.9 * q.selection()["gamma_efficiency"]
    assert p.plan()["beam_time_needed_s"] < 0.45 * q.plan()["beam_time_needed_s"]
    # Both nuclei are in its events, in sixteen crystals.
    g = simulate_gammas(exp, events=60_000, seed=3)
    assert len(g.crystal_names) == 16 and set(np.unique(g["nucleus"])) == {0, 1}


def test_a_clover_run_counts_what_its_data_sheet_efficiency_says(tmp_path):
    # The faces of a clover's crystals overlap (discs 50 mm across, 45 mm apart): 3.1% of each face is covered
    # first by a neighbour, and a γ ray through it goes to that neighbour. A crystal that stands alone has all of
    # its face.
    exp = Experiment.example("coulex_ni62_pt194_clovers")
    shares = _own_shares(crystals_of(exp)[0])
    assert len(shares) == 16 and shares == pytest.approx(0.9688, abs=0.002)
    assert _own_shares(crystals_of(Experiment.example("coulex_ni62_pt194"))[0]) == pytest.approx(1.0)
    # The efficiency is the data sheet's all the same, in the events as in the plan: the run counts the platinum
    # peak the plan predicts. (With each crystal's efficiency spread over its whole face the run lost 3% of the
    # γ rays; platinum is the line to look at, with 15 000 counts in two minutes.)
    p = Planner(exp)
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="1 h", budget="2 min")
    rows = [r for r in p.compare()["rows"] if r["quantity"] == f"Coincidences in the {PT_FIRST} peak"]
    assert len(rows) == 3
    measured, predicted = sum(r["measured"] for r in rows), sum(r["predicted"] for r in rows)
    assert measured * p.run.real_s > 14_000   # deorientation in vacuum evens the platinum γ rays out
    assert measured / predicted == pytest.approx(1.0, abs=0.015)
