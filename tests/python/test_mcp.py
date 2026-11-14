"""The planner's MCP server (physim mcp): the experiment and setup tools on saved experiments, what goes back to
the model, refusals as messages, and the server spoken to by a client, in process and over stdio."""

import asyncio
import json
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("mcp.server.mcpserver", reason="the MCP extra is not installed")

from mcp.client import Client  # noqa: E402
from mcp.client.stdio import StdioServerParameters  # noqa: E402

from physim.mcp.server import build_server, message  # noqa: E402
from physim.mcp.session import Session, ToolProblem  # noqa: E402
from physim.mcp.shape import LIMIT, plain  # noqa: E402
from physim.mcp.tools_setup import SetupTools  # noqa: E402
from physim.nuclear.planner import Planner  # noqa: E402


@pytest.fixture
def tools(tmp_path):
    return SetupTools(Session(tmp_path))


def _data(result):
    """A tool's result as the client sees it: the structured content, or the JSON text."""
    assert not result.is_error, result.content[0].text
    if result.structured_content is not None:
        return result.structured_content
    return json.loads(result.content[0].text)


# -- what goes back to the model ---------------------------------------------------------------------------------


def test_plain_summarises_long_arrays_and_keeps_short_ones():
    out = plain({"long": np.arange(LIMIT + 1, dtype=float), "short": np.arange(3), "list": list(range(LIMIT + 5)),
                 "nan": float("nan"), "n": np.int64(4), "x": np.float32(1.234567891), "flag": np.bool_(True)})
    assert out["long"] == {"summary": "array", "length": LIMIT + 1, "min": 0.0, "max": float(LIMIT),
                           "mean": LIMIT / 2}
    assert out["short"] == [0, 1, 2]
    assert out["list"]["summary"] == "array"
    assert out["nan"] is None and out["n"] == 4 and out["flag"] is True
    assert out["x"] == pytest.approx(1.23457)
    json.dumps(out)


# -- experiments -------------------------------------------------------------------------------------------------


def test_tools_need_an_open_experiment(tools):
    with pytest.raises(ToolProblem, match="No experiment is open"):
        tools.get_setup()


def test_new_experiment_from_an_example_is_saved_and_listed(tools, tmp_path):
    out = tools.new_experiment(name="Ni test", from_example="coulex_ni58")
    assert out["experiment"]["name"] == "Ni test"
    folder = Path(out["experiment"]["folder"])
    assert folder.parent == tmp_path and (folder / "experiment.toml").is_file()
    assert out["beam"]["energy"] == "30 MeV" and out["reaction"]["type"] == "coulex"
    assert out["detectors"][0].startswith("1. CD: shape annular")
    assert out["problems"] == [] and out["current_run"] is None
    json.dumps(out)
    listed = tools.list_experiments()
    assert [e["name"] for e in listed["experiments"]] == ["Ni test"]
    assert listed["experiments"][0]["runs"] == 0


def test_new_experiment_needs_exactly_one_start(tools, tmp_path):
    with pytest.raises(ToolProblem, match="exactly one"):
        tools.new_experiment()
    with pytest.raises(ToolProblem, match="No example"):
        tools.new_experiment(from_example="nope")
    setup = tmp_path / "mine.toml"
    Planner.example("alpha_on_gold").save(setup)
    out = tools.new_experiment(from_setup_file=str(setup))
    assert out["beam"]["nuclide"] == "4He"


def test_list_examples_gives_titles(tools):
    names = [e["name"] for e in tools.list_examples()["examples"]]
    assert names == Planner.examples()
    assert all(e["title"] for e in tools.list_examples()["examples"])


def test_open_experiment_by_name_folder_or_path(tools, tmp_path):
    folder = Path(tools.new_experiment(name="Gold foil", from_example="alpha_on_gold")["experiment"]["folder"])
    other = SetupTools(Session(tmp_path))
    for key in ("Gold foil", "gold foil", folder.name.lower(), str(folder)):
        assert other.open_experiment(key)["experiment"]["name"] == "Gold foil"
    with pytest.raises(ToolProblem, match="'Gold foil'"):
        other.open_experiment("silver")


# -- editing the setup -------------------------------------------------------------------------------------------


def test_an_edit_is_saved_and_matches_the_planner(tools, tmp_path):
    out = tools.new_experiment(name="Gold", from_example="alpha_on_gold")
    folder = out["experiment"]["folder"]
    out = tools.set_value("beam", "energy", "6 MeV")
    assert out["ok"] and out["beam"]["energy"] == "6 MeV"
    reopened = Planner.open_experiment(folder)
    assert reopened.experiment.beam.energy_mev == pytest.approx(6.0)
    direct = Planner.example("alpha_on_gold")
    direct.set("beam", "energy", "6 MeV")
    assert [(w.level, w.text) for w in direct.warnings()] == [(w["level"], w["text"]) for w in out["warnings"]]


def test_an_invalid_edit_reports_problems_and_keeps_the_saved_setup(tools, tmp_path):
    folder = tools.new_experiment(name="Gold", from_example="alpha_on_gold")["experiment"]["folder"]
    out = tools.set_value("beam", "energy", "-5 MeV")
    assert not out["ok"] and out["problems"] and "note" in out
    saved = Planner.open_experiment(folder).experiment
    assert saved.beam.energy_mev == pytest.approx(5.5)
    assert tools.get_setup()["problems"] == out["problems"]
    assert tools.set_value("beam", "energy", "5 MeV")["ok"]


def test_detectors_by_name_or_number(tools):
    tools.new_experiment(from_example="coulex_ni58")
    out = tools.set_value("detector CD", "distance", "35 mm")
    assert out["ok"] and "distance 35 mm" in out["detectors"][0]
    assert "distance 35 mm" in tools.set_value("detector 1", "distance", "35 mm")["detectors"][0]
    n = len(out["detectors"])
    out = tools.duplicate_detector("DSSD-L", {"theta": "60 deg"})
    assert len(out["detectors"]) == n + 1 and "DSSD-L-2" in out["detectors"][-1]
    out = tools.remove_detector(n + 1)
    assert len(out["detectors"]) == n
    out = tools.add_detector({"name": "A75", "shape": "circle", "theta": "75 deg", "distance": "100 mm",
                              "radius": "2.5 mm", "thickness": "300 um"})
    assert out["ok"] and out["detectors"][-1].startswith(f"{n + 1}. A75")
    with pytest.raises(ToolProblem, match="No detector named"):
        tools.set_value("detector nope", "distance", "1 mm")
    with pytest.raises(ToolProblem, match="Unknown shape"):
        tools.add_detector({"name": "X", "shape": "hexagon"})
    with pytest.raises(ToolProblem, match="count from 1"):
        tools.remove_detector(0)


def test_gamma_detectors(tools):
    tools.new_experiment(from_example="coulex_ni58")
    n = len(tools.get_setup()["toml"].split("[[gamma_detectors]]")) - 1
    out = tools.duplicate_gamma_detector(1, {"phi": "45 deg"})
    assert len(out["gamma_detectors"]) == n + 1
    out = tools.remove_gamma_detector(n + 1)
    assert len(out["gamma_detectors"]) == n


def test_replace_setup_checks_first(tools):
    tools.new_experiment(from_example="alpha_on_gold")
    before = tools.get_setup()["toml"]
    bad = tools.replace_setup(before.replace('energy = "5.5 MeV"', 'energy = "lots"'))
    assert not bad["ok"] and bad["problems"]
    assert tools.get_setup()["toml"] == before
    assert not tools.replace_setup("this is [not toml")["ok"]
    good = tools.replace_setup(before.replace('energy = "5.5 MeV"', 'energy = "7 MeV"'))
    assert good["ok"] and good["beam"]["energy"] == "7 MeV"


def test_level_schemes_and_the_excited_state(tools):
    out = tools.new_experiment(from_example="coulex_ni62_pt194")
    assert set(out["level_schemes"]) >= {"target"}
    out = tools.use_state("target", 1)
    assert out["ok"] and out["reaction"]["excite"] == "target"
    assert out["reaction"]["energy"].startswith("328.473")
    assert out["levels"]["target"]["levels"][1]["jpi"] == "2+"
    out = tools.set_matrix_element("target", 0, 1, "E2", "1.3 eb")
    assert out["ok"]


# -- the server --------------------------------------------------------------------------------------------------


def test_messages_drop_python_quoting():
    assert message(KeyError("no detector named 'X'")) == "no detector named 'X'"
    assert "count" in message(IndexError("there are 3 detectors")).lower()


def test_the_server_lists_and_calls_the_tools(tmp_path):
    server = build_server(tmp_path)

    async def talk():
        async with Client(server) as c:
            listed = {t.name: t for t in (await c.list_tools()).tools}
            assert set(listed) == {name for name, _ in SetupTools.TOOLS}
            assert listed["list_experiments"].annotations.read_only_hint
            assert not listed["set_value"].annotations.read_only_hint
            assert listed["replace_setup"].annotations.destructive_hint
            assert "section" in listed["set_value"].input_schema["properties"]
            assert "text with a unit" in listed["set_value"].description

            refused = await c.call_tool("get_setup", {})
            assert refused.is_error and "No experiment is open" in refused.content[0].text
            out = _data(await c.call_tool("new_experiment", {"from_example": "coulex_ni58"}))
            assert out["beam"]["energy"] == "30 MeV"
            out = _data(await c.call_tool("set_value", {"section": "detector 2", "field": "rings", "value": 8}))
            assert "problems" in out
            refused = await c.call_tool("remove_detector", {"which": "nope"})
            assert refused.is_error and refused.content[0].text.endswith(": no detector named 'nope'")

    asyncio.run(talk())
    assert server.session.planner is not None


def test_physim_mcp_over_stdio(tmp_path):
    params = StdioServerParameters(command=sys.executable,
                                   args=["-m", "physim", "mcp", "--experiments", str(tmp_path)])

    folders = []

    async def talk():
        async with Client(params) as c:
            names = {t.name for t in (await c.list_tools()).tools}
            assert "new_experiment" in names
            out = _data(await c.call_tool("new_experiment", {"name": "Over stdio", "from_example": "alpha_on_gold"}))
            assert out["experiment"]["name"] == "Over stdio"
            folders.append(out["experiment"]["folder"])
            out = _data(await c.call_tool("set_value", {"section": "beam", "field": "energy", "value": "6 MeV"}))
            assert out["ok"]

    asyncio.run(talk())
    assert Planner.open_experiment(folders[0]).experiment.beam.energy_mev == pytest.approx(6.0)


def test_physim_mcp_without_the_extra(monkeypatch, capsys):
    import physim.mcp as cli

    monkeypatch.delitem(sys.modules, "physim.mcp.server", raising=False)
    for name in [m for m in sys.modules if m == "mcp" or m.startswith(("mcp.", "mcp_types"))]:
        monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, "mcp", None)
    assert cli.main([]) == 2
    assert "physim-engine[mcp]" in capsys.readouterr().err
