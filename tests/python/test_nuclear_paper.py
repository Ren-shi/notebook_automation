"""Paper-ready figures in a journal's style (backlog item 49)."""

import pytest

from physim.nuclear import paper
from physim.nuclear.planner import Planner

pytest.importorskip("matplotlib")


def _options(p, name):
    if name == "sweep":
        return {"sweep": p.sweep("beam energy", ["4 MeV", "5 MeV"], detector=p.experiment.detectors[0].name)}
    return {"events": 20_000} if name == "spectra" else {}


@pytest.mark.parametrize("example", Planner.examples())
def test_every_figure_draws_and_has_its_data(example):
    p = Planner.example(example)
    for name in paper.FIGURES:
        options = _options(p, name)
        fig = paper.figure(p, name, **options)
        assert fig.axes, name
        lines = paper.data_csv(p, name, **options).splitlines()
        assert lines[0].startswith("series,")
        # The excitation figure has data only for a Coulomb-excitation setup.
        assert len(lines) > 1 or (name == "gamma" and p.experiment.excitation is None), name


@pytest.mark.parametrize("key", paper.JOURNALS)
def test_figures_have_the_journal_width_and_type(key):
    p = Planner.example("alpha_on_gold")
    style = paper.JOURNALS[key]
    assert "single" in style.widths and "double" in style.widths
    for width, mm in style.widths.items():
        fig = paper.figure(p, "kinematics", key, width)
        assert fig.get_figwidth() * 25.4 == pytest.approx(mm)
        ax = fig.axes[0]
        assert ax.xaxis.label.get_fontsize() == pytest.approx(style.size)
        assert ax.xaxis.label.get_fontfamily() == [style.family]
        assert fig.get_facecolor()[:3] == (1.0, 1.0, 1.0)


def test_curves_differ_without_colour():
    fig = paper.figure(Planner.example("coulex_ni58"), "kinematics")
    styles = [str(line.get_linestyle()) for line in fig.axes[0].get_lines()]
    assert len(styles) >= 4 and len(set(styles)) == len(styles)


def test_export_formats():
    p = Planner.example("alpha_on_gold")
    pdf = paper.export(p, "kinematics", "pdf", "physical_review", "single")
    assert pdf.startswith(b"%PDF")
    svg = paper.export(p, "kinematics", "svg", "nature", "double").decode()
    assert "<svg" in svg and "lab angle" in svg  # the text stays text
    png = paper.export(p, "energy_loss", "png", "elsevier", "middle")
    assert png.startswith(b"\x89PNG")
    # 140 mm at 600 dpi
    assert int.from_bytes(png[16:20], "big") == pytest.approx(140 / 25.4 * 600, abs=2)
    assert paper.preview(p, "kinematics").startswith(b"\x89PNG")


def test_wrong_choices_are_named():
    p = Planner.example("alpha_on_gold")
    with pytest.raises(ValueError, match="unknown journal"):
        paper.figure(p, "kinematics", "physics_today")
    with pytest.raises(ValueError, match="has no 'middle' width"):
        paper.figure(p, "kinematics", "nature", "middle")
    with pytest.raises(ValueError, match="unknown figure"):
        paper.figure(p, "histogram")
    with pytest.raises(ValueError, match="unknown format"):
        paper.export(p, "kinematics", "eps")
    with pytest.raises(ValueError, match="needs sweep"):
        paper.figure(p, "sweep")
    assert paper.width_mm("nature", 120) == 120.0
