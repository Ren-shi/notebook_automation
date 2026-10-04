import copy
import math
import pytest

from physim import nuclear as nu
from physim.nuclear import Experiment, Quantity, SetupError

# Shipped inside the package, so these tests also run against an installed wheel.
EXAMPLES = sorted(nu.EXAMPLES.glob("*.toml"))


def minimal():
    """The smallest valid setup, as the structure of a setup file."""
    return {
        "schema": nu.SCHEMA,
        "title": "test",
        "beam": {"nuclide": "4He", "energy": "5 MeV", "current": "1 pnA"},
        "target": {"material": "Au", "thickness": "1 mg/cm2"},
        "detectors": [{"shape": "circle", "radius": "5 mm", "theta": "30 deg", "distance": "100 mm",
                       "thickness": "300 um"}],
        "run": {"beam_time": "1 h"},
    }


def problems_of(data):
    with pytest.raises(SetupError) as info:
        Experiment.from_dict(data)
    return info.value.problems


# -- quantities and names ---------------------------------------------------------------------------------------


def test_quantity_parse_and_convert():
    q = Quantity.parse("5.5 MeV")
    assert q == Quantity(5.5, "MeV") and str(q) == "5.5 MeV" and q.kind == "energy"
    assert Quantity.parse("200 ug/cm2").to("mg/cm2") == pytest.approx(0.2)
    assert Quantity.parse("500 µm").to("mm") == pytest.approx(0.5)
    assert Quantity.parse("90 deg").to("rad") == pytest.approx(math.pi / 2)
    assert Quantity.parse("1.5e-3 m").canonical == pytest.approx(1.5)
    assert Quantity.parse("30 min").to("h") == pytest.approx(0.5)
    assert str(Quantity(40.0, "mm")) == "40 mm" and str(Quantity(0.1, "%")) == "0.1 %"
    with pytest.raises(ValueError, match="has no unit"):
        Quantity.parse("40")
    with pytest.raises(ValueError, match="unknown unit 'furlong'"):
        Quantity.parse("3 furlong")
    with pytest.raises(ValueError, match="cannot convert"):
        Quantity(1, "MeV").to("mm")


def test_nuclide_names():
    assert nu.parse_nuclide("4He") == (2, 4)
    assert nu.parse_nuclide("He-4") == (2, 4) == nu.parse_nuclide("alpha")
    assert nu.parse_nuclide("208Pb") == (82, 208)
    assert nu.parse_nuclide("p") == (1, 1) and nu.parse_nuclide("d") == (1, 2)
    with pytest.raises(ValueError, match="not an element"):
        nu.parse_nuclide("4Xx")
    with pytest.raises(ValueError, match="not possible"):
        nu.parse_nuclide("1Au")
    with pytest.raises(ValueError, match="cannot read"):
        nu.parse_nuclide("gold")


def test_material_names():
    assert nu.parse_material("Au") == [(79, None, 1.0)]
    assert nu.parse_material("197Au") == [(79, 197, 1.0)]
    assert nu.parse_material("CD2") == [(6, None, 1.0), (1, 2, 2.0)]
    assert nu.parse_material("Mylar") == [(6, None, 10.0), (1, None, 8.0), (8, None, 4.0)]
    assert nu.parse_material("13CH4") == [(6, 13, 1.0), (1, None, 4.0)]
    # Two-letter symbols starting with D or T are elements, not deuterium or tritium.
    assert nu.parse_material("TiO2") == [(22, None, 1.0), (8, None, 2.0)]
    assert nu.parse_material("Te") == [(52, None, 1.0)] and nu.parse_material("Ds") == [(110, None, 1.0)]
    for bad in ("", "gold", "Xx", "C-12", "CH2x"):
        with pytest.raises(ValueError):
            nu.parse_material(bad)


# -- reading, writing and round trips ---------------------------------------------------------------------------


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.name)
def test_examples_round_trip_exactly(path, tmp_path):
    text = path.read_text(encoding="utf-8")
    exp = Experiment.load(path)
    # The examples are kept in the canonical layout, so writing them back gives the same text.
    assert exp.to_toml() == text
    out = tmp_path / "again.toml"
    exp.save(out)
    assert Experiment.load(out) == exp
    assert out.read_text(encoding="utf-8") == text


def test_examples_exist():
    assert nu.example_names() == ["alpha_on_gold", "coulex_ni58", "oxygen_on_lead_array"]
    assert Experiment.example("alpha_on_gold") == Experiment.load(nu.EXAMPLES / "alpha_on_gold.toml")
    with pytest.raises(ValueError, match="no example setup 'nope'; available: alpha_on_gold"):
        Experiment.example("nope")


def test_derived_values():
    exp = Experiment.example("oxygen_on_lead_array")
    b = exp.beam
    assert (b.Z, b.A, b.charge) == (8, 16, 6)
    assert b.energy_mev == pytest.approx(64.0)  # 4 MeV/u x 16
    assert b.particle_current_pna == pytest.approx(10 / 6)  # 10 electrical nA at charge 6+
    assert b.particles_per_second == pytest.approx(10e-9 / 6 / 1.602176634e-19)
    cd = exp.detectors[-1]
    assert cd.position_mm() == pytest.approx((0.0, 0.0, -40.0), abs=1e-12)
    assert cd.facing_direction() == pytest.approx((0.0, 0.0, 1.0), abs=1e-12)
    d2 = exp.detectors[1]  # theta 60, phi 180, 150 mm
    assert d2.position_mm() == pytest.approx((-150 * math.sin(math.pi / 3), 0.0, 75.0), abs=1e-9)
    assert d2.detector_material == "Si"
    gold = Experiment.from_dict(minimal())
    assert gold.beam.charge == 2 and gold.beam.energy_mev == 5.0


def test_position_and_facing_placement():
    data = minimal()
    data["detectors"][0] = {"shape": "rectangle", "width": "40 mm", "height": "40 mm", "thickness": "1 mm",
                            "position": ["0 mm", "50 mm", "100 mm"], "facing": [0, 0, -1]}
    exp = Experiment.from_dict(data)
    det = exp.detectors[0]
    assert det.position_mm() == (0.0, 50.0, 100.0)
    assert det.facing_direction() == (0.0, 0.0, -1.0)
    assert Experiment.from_toml(exp.to_toml()) == exp


def test_python_edits_are_validated(tmp_path):
    exp = Experiment.from_dict(minimal())
    exp.detectors[0].theta = "60 deg"
    exp.validate()
    exp.save(tmp_path / "ok.toml")
    assert Experiment.load(tmp_path / "ok.toml").detectors[0].theta == Quantity(60, "deg")
    exp.detectors[0].distance = "-40 mm"
    with pytest.raises(SetupError, match=r"detector 1: distance must be positive, got '-40 mm'"):
        exp.save(tmp_path / "bad.toml")
    assert not (tmp_path / "bad.toml").exists()


# -- invalid input: every problem names its field ---------------------------------------------------------------


def _set(path, value):
    def edit(d):
        *parents, last = path
        for p in parents:
            d = d[p]
        if value is _DELETE:
            del d[last]
        else:
            d[last] = value
    return edit


_DELETE = object()

INVALID = [
    # (edit, expected text in the problem)
    (_set(("schema",), "physim.experiment/9"), "schema 'physim.experiment/9' is not supported"),
    (_set(("schema",), _DELETE), "'schema' is missing"),
    (_set(("title",), ""), "'title' is missing or empty"),
    (_set(("beam",), _DELETE), "the [beam] section is missing"),
    (_set(("beam", "nuclide"), "4Xx"), "beam: nuclide 'Xx' in '4Xx' is not an element symbol"),
    (_set(("beam", "energy"), 5.5), 'beam: energy 5.5 has no unit; write it as text with a unit, e.g. "5.5 MeV"'),
    (_set(("beam", "energy"), "-5 MeV"), "beam: energy must be positive, got '-5 MeV'"),
    (_set(("beam", "energy"), "5 mm"), "beam: energy '5 mm' has the wrong unit; use one of MeV"),
    (_set(("beam", "current"), "1 nA"), "beam: current unknown unit 'nA'"),
    # Text without a unit, or unreadable text: the example unit suits the field.
    (_set(("beam", "energy"), "5.5"), "beam: energy '5.5' has no unit; write it with one, e.g. '5.5 MeV'"),
    (_set(("beam", "energy"), "abc"), "beam: energy cannot read 'abc' as a number with a unit, e.g. '1 MeV'"),
    (_set(("run", "beam_time"), "12"), "run: beam_time '12' has no unit; write it with one, e.g. '12 h'"),
    (_set(("beam", "energy"), "5 MeVV"), "beam: energy unknown unit 'MeVV' in '5 MeVV'; use one of MeV,"),
    (_set(("beam", "charge_state"), 3), "beam: charge_state must be between 1 and 2 for 4He"),
    (_set(("beam", "enrgy"), "5 MeV"), "beam: unknown field 'enrgy'; did you mean 'energy'?"),
    (_set(("target", "material"), "gold"), "target: material cannot read 'gold' as a material"),
    (_set(("target", "thickness"), "0 mg/cm2"), "target: thickness must be positive"),
    (_set(("target", "tilt"), "90 deg"), "target: tilt must be between -90 and 90 deg"),
    (_set(("target", "backing"), {"material": "C"}), "target.backing: 'thickness' is missing"),
    (_set(("detectors",), []), "detectors must be written as [[detectors]] sections"),
    (_set(("detectors",), _DELETE), "no detectors: add at least one [[detectors]] section"),
    (_set(("detectors", 0, "distance"), "-40 mm"), "detector 1: distance must be positive, got '-40 mm'"),
    (_set(("detectors", 0, "theta"), "200 deg"), "detector 1: theta must be between 0 and 180 deg"),
    (_set(("detectors", 0, "theta"), _DELETE), "detector 1: 'theta' is missing (or give the detector a position)"),
    (_set(("detectors", 0, "shape"), "hexagon"), "detector 1: shape 'hexagon' is not known"),
    (_set(("detectors", 0, "radius"), _DELETE), "detector 1: a circle detector needs 'radius'"),
    (_set(("detectors", 0, "width"), "5 mm"), "detector 1: 'width' is not used by a circle detector"),
    (_set(("detectors", 0, "position"), ["0 mm", "0 mm", "10 mm"]),
     "detector 1: give either position or theta/phi/distance, not both"),
    (_set(("detectors", 0, "facing"), [0, 0, 0]), "detector 1: facing must not be the zero vector"),
    (_set(("detectors", 0, "thickness"), _DELETE), "detector 1: 'thickness' is missing"),
    (_set(("run", "beam_time"), "0 h"), "run: beam_time must be positive"),
    (_set(("run", "counts_wanted"), 0.5), "run: counts_wanted must be a whole number"),
    (_set(("reaction",), {"type": "transfer"}), "reaction: type 'transfer' is not available yet"),
    (_set(("extra",), {}), "unknown section 'extra'"),
]


@pytest.mark.parametrize("edit, expected", INVALID, ids=[e for _, e in INVALID])
def test_invalid_input_names_the_field(edit, expected):
    data = minimal()
    edit(data)
    problems = problems_of(data)
    assert any(expected in p for p in problems), problems


def test_annular_radii_and_named_detector_labels():
    data = minimal()
    data["detectors"].append({"name": "CD", "shape": "annular", "inner_radius": "40 mm", "outer_radius": "10 mm",
                              "theta": "180 deg", "distance": "40 mm", "thickness": "300 um"})
    data["detectors"].append({"name": "CD", "shape": "circle", "radius": "5 mm", "theta": "90 deg",
                              "distance": "50 mm", "thickness": "300 um"})
    problems = problems_of(data)
    assert "detector 2 (CD): inner_radius must be smaller than outer_radius" in problems
    assert "detector 3 (CD): name 'CD' is already used by detector 2" in problems


def test_all_problems_are_reported_together():
    data = minimal()
    data["beam"]["energy"] = "-1 MeV"
    data["target"]["thickness"] = "1 parsec"
    data["detectors"][0]["distance"] = "-1 mm"
    with pytest.raises(SetupError) as info:
        Experiment.from_dict(data)
    assert len(info.value.problems) == 3
    assert str(info.value).startswith("the setup has 3 problems:")


def test_not_toml():
    with pytest.raises(SetupError, match="not valid TOML"):
        Experiment.from_toml("schema = \n")


def test_minimal_setup_is_untouched_by_reading():
    data = minimal()
    before = copy.deepcopy(data)
    Experiment.from_dict(data)
    assert data == before
