"""The MCP server's lifetime tools (lifetime wave 3A): planning a lifetime measurement, setting up the plunger
through the planner's draft methods, the predicted spectra without their bins, the analysis's refusal without runs,
and the tools as the server lists them. Quantities come back as text with units, like the other tools."""

import asyncio

import pytest

pytest.importorskip("mcp.server.mcpserver", reason="the MCP extra is not installed")

from mcp.client import Client  # noqa: E402

from physim.mcp.server import build_server  # noqa: E402
from physim.mcp.session import Session, ToolProblem  # noqa: E402
from physim.mcp.tools_lifetimes import LifetimeTools  # noqa: E402
from physim.mcp.tools_setup import SetupTools  # noqa: E402


@pytest.fixture
def tools(tmp_path):
    session = Session(tmp_path)
    setup = SetupTools(session)
    setup.new_experiment(name="plunger", from_example="coulex_ru98_plunger")
    return LifetimeTools(session)


def test_the_plan_comes_back_as_text_with_units(tools):
    out = tools.plan_lifetime_measurement()
    assert out["available"] and out["nucleus"] == "98Ru" and out["method"] == "rdds"
    level = out["levels"][0]
    assert level["tau"].endswith(" ps") and "B(E2)" in level["tau_source"]
    assert level["method"] == "rdds" and level["recoil_beta"].endswith("% of c")
    assert all(d.endswith(" um") for d in level["recommended_distances"])
    assert " to " in level["sensitive_range"] and level["sensitive_range"].endswith(" um")
    assert set(level["peak_separation"]) == {"41.5 deg", "90 deg", "138.5 deg"}
    assert "resolved" in level["peak_separation"]["41.5 deg"]
    assert level["beam_time_per_distance"].split()[-1] in ("s", "min", "h")
    assert out["distances"][0].endswith(" um") and out["offset"] == "8 um"
    assert isinstance(out["warnings"], list) and isinstance(out["problems"], list)


def test_planning_another_method_changes_nothing(tools):
    before = tools.session.require().draft
    out = tools.plan_lifetime_measurement(method="dsam")
    assert "as if the method were dsam" in out["note"]
    assert tools.session.require().draft == before
    with pytest.raises(ToolProblem):
        tools.plan_lifetime_measurement(method="fast timing")


def test_set_plunger_goes_through_the_draft_methods(tools):
    out = tools.set_plunger(distances=["5 um", 20, "100 um"], offset="2 um", kind="stopper", material="Au",
                            thickness="25 mg/cm2", deorientation={})
    assert out["ok"] and out["problems"] == []
    lt = out["lifetime"]
    assert lt["distances"] == ["5 um", "20 um", "100 um"] and lt["offset"] == "2 um"
    assert lt["stopper"] == {"kind": "stopper", "material": "Au", "thickness": "25 mg/cm2"}
    assert "deorientation" not in lt
    assert any("No deorientation" in w["text"] for w in out["warnings"])  # in the summary's banner
    exp = tools.session.require().experiment
    assert exp.lifetime.distances_um == [5.0, 20.0, 100.0]
    # Only what is given changes; a foil half given takes the rest from the section.
    out = tools.set_plunger(thickness="30 mg/cm2")
    assert out["lifetime"]["stopper"]["thickness"] == "30 mg/cm2" and out["lifetime"]["stopper"]["material"] == "Au"
    out = tools.set_plunger(kind="differential", material="Mg", thickness="1 mg/cm2",
                            second_foil={"material": "Au", "thickness": "30 mg/cm2", "gap": "150 um"})
    assert out["ok"], out["problems"]
    assert out["lifetime"]["degrader2"]["gap"] == "150 um"
    with pytest.raises(ToolProblem):
        tools.set_plunger(second_foil={"material": "Au"})
    out = tools.set_plunger(remove=True)
    assert out["ok"] and out["lifetime"] is None and tools.session.require().experiment.lifetime is None


def test_an_invalid_plunger_is_kept_with_its_problems(tools):
    out = tools.set_plunger(kind="stopper", material="Cu", thickness="0.5 mg/cm2")
    assert not out["ok"] and any("does not stop" in p for p in out["problems"])
    assert "not valid" in out["note"]


def test_spectra_leave_out_the_bins(tools):
    out = tools.lifetime_spectra(distance="25 um")
    assert out["available"] and out["method"] == "rdds" and "analytic" in out["model"]
    assert {r["distance"] for r in out["spectra"]} == {"25 um"}
    row = out["spectra"][0]
    assert row["shifted_at"].endswith(" keV") and row["counts_in_line"] >= 0 and 0 <= row["unshifted_share"] <= 1
    assert "counts" not in row and out["energy_range"].endswith(" bins")
    assert tools.lifetime_spectra(distance=25.0)["spectra"][0]["distance"] == "25 um"


def test_analysis_without_runs_says_what_to_do(tools):
    out = tools.analyse_lifetimes()
    assert not out["available"] and "No plunger runs" in out["reason"]
    with pytest.raises(ToolProblem):
        tools.analyse_lifetimes(method="napatau")


def test_the_server_lists_the_lifetime_tools(tmp_path):
    server = build_server(tmp_path)

    async def talk():
        async with Client(server) as c:
            listed = {t.name: t for t in (await c.list_tools()).tools}
            assert {name for name, _ in LifetimeTools.TOOLS} <= set(listed)
            assert listed["plan_lifetime_measurement"].annotations.read_only_hint
            assert listed["lifetime_spectra"].annotations.read_only_hint
            assert listed["analyse_lifetimes"].annotations.read_only_hint
            assert not listed["set_plunger"].annotations.read_only_hint
            assert "distances" in listed["set_plunger"].input_schema["properties"]
            refused = await c.call_tool("plan_lifetime_measurement", {})
            assert refused.is_error and "No experiment is open" in refused.content[0].text

    asyncio.run(talk())
