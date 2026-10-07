"""The Physics tab (backlog item 70): the register as a reference, and every explained number pointing at the
section that shows its physics."""

from pathlib import Path

from physim.nuclear import workbench as wb
from physim.nuclear.planner import REGISTER, Planner

DOCS = Path(__file__).resolve().parents[2] / "docs"


def test_sections_follow_the_register_and_their_pages_and_tests_exist():
    keys = [k for k, _, _, _ in wb.PHYSICS]
    assert keys == ["data", "kinematics", "stopping", "rutherford", "detectors", "rates", "coulex"]
    for key, title, page, tests in wb.PHYSICS:
        assert (DOCS / f"{page}.md").is_file(), page
        for t in tests:
            assert (Path(__file__).parent / t).is_file(), t
    # Every page of the planner's register is a section's page.
    pages = {page for _, _, page, _ in wb.PHYSICS}
    assert {p for t, p in REGISTER.items() if t not in ("report", "spectra")} <= pages | {"physics-register/README"}


def test_every_explained_number_opens_a_physics_section(tmp_path):
    from physim.nuclear import logbook
    from physim.nuclear.record import analysis_explanations, explanations

    sections = {k for k, _, _, _ in wb.PHYSICS}
    p = Planner.example("coulex_ni58")
    keys = [x.key for x in explanations(p)]
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="5 min", budget="1 min")
    out = p.analyse(settings={"particle_gate": "all"})
    if out["available"]:
        keys += [x.key for x in analysis_explanations(p, p.analysis().history[-1])]
    keys += [x.key for x in logbook.run_explanations(p)]
    for k in keys:
        assert wb.physics_for(k) in sections, k
    # The specific ones go where their physics is, not to the default.
    assert wb.physics_for("solid_angle") == "detectors" and wb.physics_for("rate") == "rates"
    assert wb.physics_for("weisskopf") == "data" and wb.physics_for("doppler") == "coulex"
    # And every setup card's consequence line.
    for d in wb.decisions(p):
        assert wb.physics_for(d.key) in sections
    assert wb.physics_for("target") == "stopping" and wb.physics_for("detector:1") == "detectors"
