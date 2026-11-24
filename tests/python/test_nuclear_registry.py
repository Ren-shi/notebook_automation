"""The known physics of the nucleus view: whether a route applies to a nucleus, the checks between observables, what
a value is used for, and B(M1) and μ in the registry. Synthetic rows and no local ENSDF copy: seconds."""

import pytest

from physim.nuclear import compare, ensdf, nucleus, observables, results
from physim.nuclear import shellmodel_spaces as spaces
from physim.nuclear.analysis import Analysis
from physim.nuclear.experiment import Experiment
from physim.nuclear.observables import Route

#: ⁶²Ni's lowest levels as results.ensdf_levels gives them: (keV, 2J, parity, n).
NI62 = ((0.0, 0, 1, 1), (1172.91, 4, 1, 1), (2048.55, 0, 1, 2), (2336.5, 8, 1, 1))
FIRST, SECOND = "62Ni:2+1>0+1:E2", "62Ni:4+1>2+1:E2"
ENSDF_ORIGIN = results.reference("ENSDF, adopted levels and γ rays", 2025, where="NDS 204, 1 (2025)",
                                 evaluation=nucleus.ENSDF, dataset_date="202507")


@pytest.fixture(autouse=True)
def lab(tmp_path, monkeypatch):
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path / "experiments"))
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setattr(results, "ensdf_levels", lambda nuclide: NI62 if nuclide == "62Ni" else ())
    monkeypatch.setattr(nucleus, "decay", lambda nuclide, level: (None, None, "no γ ray known"))


def row(observable, key, value, unit, method="evaluation", evidence="evaluated", origin=ENSDF_ORIGIN, **more):
    return results.Row("62Ni", observable, key, value, unit, method, evidence, origin=origin, **more)


def ni62():
    """⁶²Ni as the evaluation has it, a shell-model B(E2), and simulated runs of everything."""
    jun45 = results.shell_model_job("JUN45", "f5pg9")
    run = results.from_run("62Ni-on-194Pt", "03-beam-24h")
    return [
        row("E_level", "62Ni:2+1", 1172.91, "keV", total=0.03), row("E_level", "62Ni:4+1", 2336.5, "keV", total=0.04),
        row("tau", "62Ni:2+1", 2.08, "ps", total=0.03), row("Qs", "62Ni:2+1", 5.0, "efm2", total=12.0),
        row("B(E2)", FIRST, 12.13, "W.u.", direction="down", total=0.18),
        row("B(E2)", SECOND, 21.0, "W.u.", direction="down", total=6.0),
        row("B(E2)", FIRST, 408.0, "e2fm4", "shell-model", "theory", jun45, direction="up"),
        row("B(E2)", FIRST, 918.6, "e2fm4", "coulex-yield", "projected", run, direction="up", total=65.0, truth=884.1),
        row("Qs", "62Ni:2+1", -27.0, "efm2", "coulex-fit", "projected", run, total=3.0, truth=-27.0),
        row("tau", "62Ni:2+1", 0.3, "ps", "dsam", "projected", run, total=0.1, truth=0.3)]


# -- whether a route applies ----------------------------------------------------------------------------------------


def test_a_target_or_a_beam():
    coulex = Route("coulex-yield", "coulex")
    fits, why = nucleus.applies("62Ni", "B(E2)", FIRST, coulex)
    assert fits is True and "3.63 % of natural Ni" in why and "target or the beam" in why
    fits, why = nucleus.applies("68Ni", "B(E2)", "68Ni:2+1>0+1:E2", coulex)  # no test for it: not a False
    assert fits is None and "not in the natural element" in why and "not checked" in why
    fits, why = nucleus.applies("68Ni", "B(E2)", "68Ni:2+1>0+1:E2", Route("nrf"))
    assert fits is None and "needs a target" in why
    # The yield analysis is of the excitation from the ground state; the fit takes several steps.
    fits, why = nucleus.applies("62Ni", "B(E2)", SECOND, coulex)
    assert fits is False and "fit to all orders" in why
    assert nucleus.applies("62Ni", "B(E2)", SECOND, Route("coulex-fit", "coulex"))[0] is True
    assert nucleus.applies("62Ni", "B(E2)", FIRST, Route("inelastic-scattering", "inelastic"))[0] is None
    assert nucleus.applies("62Ni", "mu", "62Ni:2+1", Route("moment")) == (
        None, "physim has no test for this route: nothing is checked")
    assert nucleus.applies("62Ni", "Qs", "62Ni:0+1", Route("coulex-fit", "coulex"))[0] is False  # no Q_s for J = 0
    assert nucleus.applies("62Ni", "B(E2)", FIRST, Route("evaluation"))[0] is None  # no local ENSDF copy


def test_the_lifetime_against_the_range_of_the_method():
    plunger, dsam = Route("rdds", "lifetimes"), Route("dsam", "lifetimes")
    rows = ni62()
    fits, why = nucleus.applies("62Ni", "B(E2)", FIRST, plunger, rows)
    assert fits is True and why.startswith("τ = 2.08 ps (ENSDF adopted") and "1 ps to 1 ns" in why
    fits, why = nucleus.applies("62Ni", "tau", "62Ni:2+1", dsam, rows)
    assert fits is False and "above the range" in why and "a plunger" in why and "10 fs to 1 ps" in why
    # A simulated run's τ (0.3 ps here) and a calculated one decide nothing.
    alone = [r for r in rows if r.observable == "tau" and r.evidence == "projected"]
    alone.append(row("tau", "62Ni:2+1", 4.5, "ps", "shell-model", "theory", results.shell_model_job("JUN45")))
    fits, why = nucleus.applies("62Ni", "tau", "62Ni:2+1", dsam, alone)
    assert fits is None and "not tested" in why
    # A published value by the method is not said to be out of its reach.
    kenn = row("tau", "62Ni:2+1", 2.01, "ps", "dsam", "measured", results.reference("Kenn et al.", 2001), total=0.12)
    route = next(r for r in nucleus.routes("62Ni", "tau", "62Ni:2+1", [kenn], every=rows + [kenn])
                 if r["method"] == "dsam")
    assert (route["status"], route["applies"]) == ("done", None)
    assert route["reason"].startswith("a value by this method is here") and route["ids"] == [kenn.id]


def test_whether_a_shell_model_calculation_fits(monkeypatch):
    route = Route("shell-model")
    monkeypatch.setattr(spaces, "available_memory_bytes", lambda: 16 * 1024**3)
    fits, why = nucleus.applies("24Mg", "B(E2)", "24Mg:2+1>0+1:E2", route)
    assert fits is True and "sd space" in why and "half of the 16 GB" in why
    monkeypatch.setattr(spaces, "available_memory_bytes", lambda: 0)
    fits, why = nucleus.applies("24Mg", "B(E2)", "24Mg:2+1>0+1:E2", route)
    assert fits is False and "Too big" in why
    fits, why = nucleus.applies("158Dy", "B(E2)", "158Dy:2+1>0+1:E2", route)
    assert fits is False and "none of the interactions" in why
    assert nucleus.applies("24Mg", "E_level", "24Mg:?@1368", route)[0] is None  # no spin: nothing to size


# -- the checks between observables ---------------------------------------------------------------------------------


def test_the_checks_are_plain_ratios_of_independent_rows():
    rows = ni62()
    results.add_points([r for r in rows if r.evidence != "projected"])
    v = nucleus.observable_view("62Ni", "B(E2)", FIRST, ledger=rows)
    e42, b42, rotor, theory = v["checks"]
    assert all(set(c) == {"label", "text", "how", "reference", "rows", "value"} and c["how"] for c in v["checks"])
    projected = {r.id for r in rows if r.evidence == "projected"}
    assert not any(set(c["rows"]) & projected for c in v["checks"])
    assert e42["value"] == pytest.approx(2336.5 / 1172.91) and e42["rows"] == [rows[1].id, rows[0].id]
    assert b42["value"] == pytest.approx(21.0 / 12.13)  # ENSDF's own two numbers, both ↓ in W.u.
    assert theory["value"] == pytest.approx(408 / 884.1, abs=1e-3) and theory["label"].startswith("Shell model, JUN45")
    # The rotor reading is the Coulex analysis's own, to the digit, and the evaluated Q_s stands beside it.
    b = compare.canonical(rows[4]).value
    analysis = Analysis.__new__(Analysis)
    analysis.experiment = Experiment.example("coulex_ni62_pt194")
    shape = analysis._shape(b, 62)
    assert rotor["value"] == shape["qs_2plus_efm2"] < 0
    assert observables.rotor_reading(b, "62Ni")["beta2"] == shape["beta2"]
    assert "Evaluated: Q_s(2⁺₁) = 5 ± 12 e fm² (ENSDF adopted" in rotor["text"] and rows[3].id in rotor["rows"]
    # The same from the side of Q_s and of the energies; the next transition has its own; the report's rows are used.
    assert [c["label"] for c in nucleus.observable_view("62Ni", "Qs", "2+1", ledger=rows)["checks"]] == [
        e42["label"], rotor["label"]]
    assert nucleus.observable_view("62Ni", "E_level", "4+1")["checks"] == [e42]
    assert [c["value"] for c in nucleus.observable_view("62Ni", "B(E2)", SECOND)["checks"]] == [
        e42["value"], b42["value"]]
    assert nucleus.observable_view("62Ni", "B(E2)", FIRST, ledger=[])["checks"] == []
    with pytest.raises(ValueError, match="closed against its own truth"):
        compare.ratio(rows[7], rows[4])


def test_a_moment_without_a_sign_says_so():
    rows = ni62()
    rows[3] = row("Qs", "62Ni:2+1", 5.0, "efm2", total=12.0, detail="ENSDF gives no sign")
    rotor = nucleus.observable_view("62Ni", "Qs", "2+1", ledger=rows)["checks"][-1]
    assert "ENSDF gives no sign for it" in rotor["text"]


# -- the registry ---------------------------------------------------------------------------------------------------


def test_what_a_value_is_used_for_and_the_new_observables():
    for name, o in observables.OBSERVABLES.items():
        used = nucleus.observable_view("62Ni", name, FIRST.replace("E2", o.multipolarity or "")
                                       if o.belongs_to == "transition" else "2+1")["used_for"]
        assert used == [{"text": u.text, "reference": u.reference} for u in o.used_for] and used
    assert (observables.get("B(M1)").unit, observables.get("mu").unit) == ("mu_N2", "mu_N")
    # 1 W.u. of B(M1)↓ is 1.790 μ_N² whatever the nucleus, and B(M1)↑ = B(M1)↓ between two levels of one spin.
    for nuclide in ("62Ni", "158Dy"):
        assert observables.to_canonical("B(M1)", 1.0, "W.u.", "down", nuclide=nuclide, j_lower=2,
                                        j_upper=2) == pytest.approx(1.790, abs=1e-3)
    point = nucleus.add_data_point("62Ni", "mu", "2+1", 0.33, "mu_N", uncertainty=0.05, method="moment",
                                   reference={"authors": "A. Author", "year": 2001})
    v = nucleus.observable_view("62Ni", "mu", "2+1")
    assert v["label"] == "μ(2+₁)" and [r["status"] for r in v["routes"] if r["method"] == "moment"] == ["done"]
    assert point.id in v["cells"]
    # A published inelastic-scattering value counts for its route, which "other" could not.
    nucleus.add_data_point("62Ni", "B(E2)", FIRST, 850.0, "e2fm4", direction="up", uncertainty=60.0,
                           method="inelastic-scattering", reference={"authors": "B. Author", "year": 1970})
    routes = {r["method"]: r["status"] for r in nucleus.observable_view("62Ni", "B(E2)", FIRST)["routes"]}
    assert routes["inelastic-scattering"] == "done" and "other" not in routes
