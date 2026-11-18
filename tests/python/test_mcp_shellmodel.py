"""The MCP server's shell-model tools: options and size estimates, the memory refusal, jobs from the cache, results
compared with a level scheme, applying values without overriding ENSDF or derived ones, importing a table, and the
tools as the server lists them. One real KSHELL run at the end, skipped without KSHELL."""

import asyncio
import json
import os
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp.server.mcpserver", reason="the MCP extra is not installed")

from mcp.client import Client  # noqa: E402
from mcp.client.stdio import StdioServerParameters  # noqa: E402

from physim.external import kshell  # noqa: E402
from physim.mcp.server import build_server  # noqa: E402
from physim.mcp.session import Session, ToolProblem  # noqa: E402
from physim.mcp.tools_setup import SetupTools  # noqa: E402
from physim.mcp.tools_shellmodel import ShellModelTools, parse_states, states_text  # noqa: E402
from physim.nuclear import ensdf, shellmodel  # noqa: E402
from physim.nuclear import shellmodel_jobs as jobs  # noqa: E402
from physim.nuclear import shellmodel_spaces as spaces  # noqa: E402
from physim.nuclear.planner import Planner  # noqa: E402

REF = Path(__file__).parents[1] / "reference" / "nuclear" / "shellmodel"
LOGS = REF / "kshell" / "si28_usdb" / "logs"
TABLE = REF / "import" / "si28_usdb_table.csv"
GB = 1024**3
SI28_STATES = "0+3, 2+3, 4+2"


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    """A shell-model cache of the test's own, and plenty of free memory unless a test says otherwise."""
    root = tmp_path / "shellmodel"
    monkeypatch.setenv("PHYSIM_SHELLMODEL", str(root))
    monkeypatch.setattr(spaces, "available_memory_bytes", lambda: 16 * GB)
    return root


@pytest.fixture
def fake_install(tmp_path, monkeypatch):
    """A KSHELL 'install' that is never run: the tools only need its build name for the job ID."""
    folder = tmp_path / "kshell"
    folder.mkdir()
    install = kshell.KshellInstall(folder, folder / "kshell.exe", folder / "transit.exe", None, None, "test build\n")
    monkeypatch.setattr(kshell, "find_install", lambda: install)
    return install


@pytest.fixture
def tools(tmp_path):
    return ShellModelTools(Session(tmp_path / "experiments"))


def card(nucid, kind, cols=(), cont=" "):
    """One 80-column ENSDF record: ``cols`` are (first column counting from 1, text)."""
    line = list(f"{nucid:<5}{cont} {kind} ".ljust(80))
    for start, text in cols:
        line[start - 1:start - 1 + len(text)] = text
    return "".join(line)


@pytest.fixture
def ensdf_si28(tmp_path, monkeypatch):
    """A local 'ENSDF' folder with a hand-typed ²⁸Si scheme (values close to the evaluated ones)."""
    def c(kind, cols, cont=" "):
        return card(" 28SI", kind, cols, cont)

    lines = [c(" ", [(10, "ADOPTED LEVELS, GAMMAS")]),
             c("L", [(10, "0.0"), (22, "0+"), (40, "STABLE")]),
             c("L", [(10, "1779.03"), (20, "11"), (22, "2+"), (40, "0.475 PS"), (50, "17")]),
             c("G", [(10, "1779.0"), (20, "11"), (22, "100"), (32, "E2")]),
             c("G", [(10, "BE2W=13.2 5")], cont="2"),
             c("L", [(10, "4617.86"), (20, "4"), (22, "4+"), (40, "34 FS")]),
             c("G", [(10, "2838.3"), (22, "100"), (32, "E2")])]
    folder = tmp_path / "ensdf"
    folder.mkdir()
    (folder / "ensdf.028").write_text("\n".join(lines))
    monkeypatch.setenv(ensdf.ENV, str(folder))
    ensdf._adopted.cache_clear()
    yield
    ensdf._adopted.cache_clear()


@pytest.fixture
def si28_experiment(tools, ensdf_si28):
    """An open experiment with a ²⁸Si target whose level scheme came from (the fake) ENSDF."""
    setup = SetupTools(tools.session)
    setup.new_experiment(name="Si28 coulex", from_example="coulex_ni58")
    assert setup.set_value("target", "material", "28Si")["ok"]
    out = setup.lookup_levels("target", "5 MeV")
    assert out["ok"]
    return out["experiment"]["folder"]


def _cached_job(install, eff_charge=None):
    """Put the ²⁸Si/USDB fixture logs into the cache as a finished job; returns the job ID."""
    req = shellmodel.ShellModelRequest("28Si", "usdb", parse_states(SI28_STATES), ("E2", "M1"), eff_charge)
    key = jobs.job_key(req, install.version)
    folder = jobs.job_folder(key)
    folder.mkdir(parents=True)
    result = shellmodel.read_folder(LOGS, jobs.resolved(req))
    (folder / "results.json").write_text(json.dumps({**result.to_dict(), "job": key}), encoding="utf-8")
    return key


# -- reading the arguments ---------------------------------------------------------------------------------------


def test_states_are_read_as_j_parity_and_how_many():
    assert parse_states(SI28_STATES) == ((0, 1, 3), (4, 1, 3), (8, 1, 2))
    assert parse_states("0+3 2+3") == ((0, 1, 3), (4, 1, 3))
    assert parse_states("2+, 3/2-2; 1-x2") == ((4, 1, 1), (3, -1, 2), (2, -1, 2))
    assert parse_states(["0+3", "2+ 1"]) == ((0, 1, 3), (4, 1, 1))
    assert states_text(parse_states("0+3, 3/2-2")) == "0+3, 3/2-2"


@pytest.mark.parametrize("text, words", [
    ("2", "cannot read '2'"), ("0 3", "cannot read '0'"), ("two plus", "cannot read 'two'"),
    ("2/2+", "odd numerator"), ("3/4+", "odd numerator"), ("0+0", "at least one state"),
    ("0+3, 0+1", "asked for twice"), ("", "no states"), (" , ", "no states"), (5, "not 5"),
])
def test_unreadable_states_are_refused_with_what_to_write(text, words):
    with pytest.raises(ToolProblem, match=words):
        parse_states(text)


def test_requests_that_cannot_be_made_are_refused(tools):
    with pytest.raises(ToolProblem, match="2J must be even"):
        tools.estimate_shell_model("28Si", "usdb", "3/2+2")
    with pytest.raises(ToolProblem, match="unknown interaction 'nope'.*shell_model_options"):
        tools.estimate_shell_model("28Si", "nope", SI28_STATES)
    with pytest.raises(ToolProblem, match="does not fit in the sd space"):
        tools.estimate_shell_model("208Pb", "usdb", "0+1")
    with pytest.raises(ToolProblem, match="cannot read 'Si' as a nuclide"):
        tools.estimate_shell_model("Si", "usdb", "0+1")
    with pytest.raises(ToolProblem, match="open an experiment"):
        tools.estimate_shell_model("target", "usdb", "0+1")


# -- options and estimates ---------------------------------------------------------------------------------------


def test_options_list_the_interactions_with_charges_citations_and_a_typical_size(tools):
    out = tools.shell_model_options("28Si")
    json.dumps(out)
    names = [o["interaction"] for o in out["interactions"]]
    assert {"usda", "usdb"} <= set(names) and "kb3g" not in names
    usdb = out["interactions"][names.index("usdb")]
    assert (usdb["space"], usdb["core"], usdb["valence_protons"], usdb["valence_neutrons"]) == ("sd", "16O", 6, 6)
    assert usdb["parities"] == ["+"] and usdb["m_scheme_dimension"] == {"+": 93710}
    assert usdb["effective_charges"] == {"e_p": 1.5, "e_n": 0.5, "unit": "e"}
    assert "Phys. Rev. C" in usdb["citation"] and usdb["defaults_from"]
    typical = usdb["typical_request"]
    assert typical["states"] == SI28_STATES and typical["fits"] is True
    assert typical["dimensions"] == {"0+": 93710, "2+": 81122, "4+": 52175}
    assert typical["limit"] == "8 GB" and out["memory_limit"].startswith("8 GB (half of the 16 GB")
    odd = tools.shell_model_options("27Al")
    assert odd["interactions"][0]["typical_request"]["states"] == "1/2+2, 3/2+2, 5/2+2"
    assert tools.shell_model_options("208Pb")["interactions"] == []


def test_options_default_to_the_open_experiments_target_or_beam(tools, si28_experiment):
    assert tools.shell_model_options()["nuclide"] == "28Si"
    beam = tools.shell_model_options(role="beam")
    assert beam["nuclide"] == "16O" and "usdb" not in [o["interaction"] for o in beam["interactions"]]
    with pytest.raises(ToolProblem, match="role"):
        tools.shell_model_options(role="ejectile")


def test_the_estimate_gives_dimensions_memory_and_the_limit(tools, monkeypatch):
    out = tools.estimate_shell_model("28Si", "usdb", SI28_STATES, threads=6)
    assert out["fits"] and out["lanczos_vectors"] == "in memory"
    assert out["dimensions"]["0+"] == 93710 and out["limit"] == "8 GB" and out["free_ram_now"] == "16 GB"
    assert "6 threads" in out["rough_time"] and "rough" in out["rough_time"]
    monkeypatch.setattr(spaces, "available_memory_bytes", lambda: 20 * 1024**2)  # 20 MB free: limit 10 MB
    small = tools.estimate_shell_model("28Si", "usdb", SI28_STATES)
    assert not small["fits"] and "Too big" in small["reason"] and small["limit"] == "10 MB"


# -- running: the memory rule and the cache ----------------------------------------------------------------------


def test_a_run_over_half_the_free_memory_is_refused_and_starts_nothing(tools, fake_install, monkeypatch, cache):
    monkeypatch.setattr(spaces, "available_memory_bytes", lambda: 20 * 1024**2)
    submitted = []
    monkeypatch.setattr(jobs, "submit", lambda req, **kw: submitted.append(req) or "0123456789abcdef")
    with pytest.raises(ToolProblem, match="Refused: Too big.*Nothing was started.*explicitly"):
        tools.run_shell_model("28Si", "usdb", SI28_STATES)
    assert submitted == [] and jobs.list_jobs() == []

    # Only on the user's explicit instruction: allow_over_limit goes into the request.
    monkeypatch.setattr(jobs, "status", lambda job: {"job": job, "state": "queued", "message": "Queued."})
    out = tools.run_shell_model("28Si", "usdb", SI28_STATES, allow_over_limit=True)
    assert out["job"] == "0123456789abcdef" and out["state"] == "queued"
    assert submitted[0].allow_over_limit and "explicit" in out["warning"]
    assert out["estimate"]["fits"] is False and "shell_model_status" in out["next"]


def test_a_run_needs_kshell(tools, monkeypatch):
    monkeypatch.setattr(kshell, "find_install", lambda: None)
    with pytest.raises(ToolProblem, match="PHYSIM_KSHELL.*import_shell_model"):
        tools.run_shell_model("28Si", "usdb", SI28_STATES)


def test_a_run_in_the_cache_returns_its_results_at_once(tools, fake_install, monkeypatch):
    job = _cached_job(fake_install)
    monkeypatch.setattr(spaces, "available_memory_bytes", lambda: 20 * 1024**2)  # a cache hit is never refused
    out = tools.run_shell_model("28Si", "USDB", "4+2, 0+3, 2+3", transitions="M1, E2", eff_charge="1.5, 0.5")
    assert out["job"] == job and out["cached"] and out["state"] == "done"
    assert [s["state"] for s in out["states"]][:3] == ["0+1", "2+1", "4+1"]
    json.dumps(out)


def test_status_list_and_cancel(tools, fake_install, cache):
    job = _cached_job(fake_install)
    st = tools.shell_model_status(job, wait_seconds=5)
    assert st["state"] == "done" and job in st["next"]
    listed = tools.list_shell_model_jobs()
    assert [j["job"] for j in listed["jobs"]] == [job]
    assert tools.cancel_shell_model(job)["state"] == "done"  # a finished job is left as it is
    for bad in ("0123456789abcdef", "not-a-job"):
        with pytest.raises(ToolProblem, match="list_shell_model_jobs"):
            tools.shell_model_status(bad)
        with pytest.raises(ToolProblem, match="list_shell_model_jobs"):
            tools.cancel_shell_model(bad)


def test_a_failed_job_says_why(tools, cache):
    job = "fedcba9876543210"
    folder = cache / job
    folder.mkdir(parents=True)
    (folder / "status.json").write_text(json.dumps({
        "job": job, "state": "failed", "message": "kshell.exe failed", "log_tail": "error: out of memory",
        "request": {"nuclide": "28Si", "interaction": "usdb", "states": [[0, 1, 3]], "transitions": ["E2"]}}))
    st = tools.shell_model_status(job)
    assert st["state"] == "failed" and st["log_tail"] == "error: out of memory"
    assert st["request"]["states"] == "0+3"
    with pytest.raises(ToolProblem, match="(?s)failed.*out of memory"):
        tools.shell_model_results(job)


# -- results, the comparison with a level scheme, and apply ------------------------------------------------------


def test_results_are_numbers_with_units(tools, fake_install):
    job = _cached_job(fake_install)
    out = tools.shell_model_results(job)
    ref = shellmodel.read_folder(LOGS)
    assert out["nuclide"] == "28Si" and out["interaction"] == "usdb"
    assert out["effective_charges"] == {"e_p": 1.5, "e_n": 0.5, "unit": "e"}
    two = next(s for s in out["states"] if s["state"] == "2+1")
    assert two["ex_kev"] == pytest.approx(1000 * ref.state(4, 1, 1).ex_mev, abs=0.1)
    assert two["q_s_efm2"] == pytest.approx(ref.state(4, 1, 1).q_s, rel=1e-5)
    e2 = next(t for t in out["transitions"] if t["kind"] == "E2" and (t["from"], t["to"]) == ("2+1", "0+1"))
    t = ref.transition("E2", (4, 1, 1), (0, 1, 1))
    assert e2["b_down_e2fm4"] == pytest.approx(t.b_down, rel=1e-5)
    assert e2["matrix_element_efm2"] == pytest.approx(t.mred, rel=1e-5) and e2["sign_known"]
    assert any(t["kind"] == "M1" and "b_down_mu_N2" in t for t in out["transitions"])
    assert not any(t["from"] == t["to"] for t in out["transitions"])  # moments are with the states
    json.dumps(out)


def test_results_compared_with_the_experiments_level_scheme(tools, fake_install, si28_experiment):
    job = _cached_job(fake_install)
    out = tools.shell_model_results(job, role="target")
    cmp = out["comparison"]
    pairs = {p["state"]: p for p in cmp["pairs"]}
    assert pairs["2+1"]["level"] == 1 and pairs["2+1"]["known_kev"] == pytest.approx(1779.03)
    assert pairs["2+1"]["known_source"] == "ensdf"
    assert pairs["2+1"]["delta_kev"] == pytest.approx(pairs["2+1"]["shell_kev"] - 1779.03, abs=0.1)
    assert pairs["0+2"]["level"] is None and "new" in pairs["0+2"]["note"]
    assert cmp["rms_over_states"] == 2 and cmp["rms_deviation_kev"] > 0
    assert any("kept" in line and "replace=True" in line for line in cmp["apply_would"])
    # The preview changed nothing.
    assert tools.session.planner.levels()["target"]["matrix_elements"][0]["source"] == "derived"
    with pytest.raises(ToolProblem, match="is for 28Si, but the beam nuclide is 16O"):
        tools.shell_model_results(job, role="beam")


def test_apply_fills_gaps_and_keeps_ensdf_and_derived_values(tools, fake_install, si28_experiment):
    job = _cached_job(fake_install)
    out = tools.apply_shell_model(job, "target")
    assert out["ok"] and out["changes"] and out["source"] == job
    t = out["levels"]["target"]
    assert [lv["energy_kev"] for lv in t["levels"]] == [0.0, pytest.approx(1779.03), pytest.approx(4617.86)]
    assert {lv["source"] for lv in t["levels"]} == {"ensdf"}
    first = next(m for m in t["matrix_elements"] if (m["from"], m["to"]) == (0, 1))
    assert first["source"] == "derived"
    q = next(m for m in t["matrix_elements"] if m["from"] == m["to"] == 1)
    assert q["source"] == "shell model" and q["q_efm2"] > 0
    assert any("kept" in c for c in out["changes"])
    saved = Planner.open_experiment(si28_experiment).levels()["target"]
    assert {(m["from"], m["to"]): m["source"] for m in saved["matrix_elements"]}[(1, 1)] == "shell model"

    # replace=True: the user's explicit choice overwrites the derived value.
    out = tools.apply_shell_model(job, "target", which="E2", replace=True)
    first = next(m for m in out["levels"]["target"]["matrix_elements"] if (m["from"], m["to"]) == (0, 1))
    ref = shellmodel.read_folder(LOGS).transition("E2", (4, 1, 1), (0, 1, 1))
    assert first["source"] == "shell model" and abs(first["value_efm"]) == pytest.approx(abs(ref.mred))
    assert {lv["source"] for lv in out["levels"]["target"]["levels"]} == {"ensdf"}  # energies were not chosen

    with pytest.raises(ToolProblem, match="role"):
        tools.apply_shell_model(job, "projectile")
    with pytest.raises(ToolProblem, match="cannot apply 'bogus'"):
        tools.apply_shell_model(job, "target", which="bogus")


def test_apply_needs_an_open_experiment_and_finished_results(tools, fake_install):
    job = _cached_job(fake_install)
    with pytest.raises(ToolProblem, match="No experiment is open"):
        tools.apply_shell_model(job, "target")
    with pytest.raises(ToolProblem, match="neither a shell-model job ID"):
        tools.shell_model_results("no such thing")


# -- import ------------------------------------------------------------------------------------------------------


def test_import_a_table_then_apply_it_by_its_handle(tools, si28_experiment):
    out = tools.import_shell_model(str(TABLE), role="target")
    assert out["handle"] == str(TABLE.resolve())
    assert out["imported_from"]["format"] == "physim shell-model table"
    assert out["effective_charges"] == {"e_p": 1.31, "e_n": 0.46, "unit": "e"}
    two = next(s for s in out["states"] if s["state"] == "2+1")
    assert two["ex_kev"] == 1932.0 and two["q_s_efm2"] == pytest.approx(18.5)
    e2 = next(t for t in out["transitions"] if (t["from"], t["to"]) == ("2+1", "0+1"))
    assert e2["b_down_e2fm4"] == pytest.approx(78.2) and e2["matrix_element_efm2"] == pytest.approx(-19.7738)
    assert out["comparison"]["pairs"][1]["known_kev"] == pytest.approx(1779.03)
    assert not any(n.startswith("imported from") for n in out["notes"])

    applied = tools.apply_shell_model(out["handle"], "target", which="Q, 2+1")
    q = next(m for m in applied["levels"]["target"]["matrix_elements"] if m["from"] == m["to"] == 1)
    assert q["source"] == "shell model" and q["q_efm2"] == pytest.approx(18.5, rel=1e-3)
    assert "1.31" in q["note"]


def test_import_refuses_what_it_cannot_read(tools, tmp_path):
    with pytest.raises(ToolProblem, match="No such file"):
        tools.import_shell_model(str(tmp_path / "missing.csv"))
    bad = tmp_path / "bad.csv"
    bad.write_text("# format: physim shell-model table 1\n# nuclide: 28Si\n# interaction: usdb\n"
                   "kind,state,final,value,unit\nEx,2+1,,1.9,parsecs\n")
    with pytest.raises(ToolProblem, match="Cannot import"):
        tools.import_shell_model(str(bad))
    lpt = tmp_path / "si28.lpt"
    lpt.write_text("NuShellX output\n")
    with pytest.raises(ToolProblem, match="NuShellX"):
        tools.import_shell_model(str(lpt))


# -- the server --------------------------------------------------------------------------------------------------


def test_the_server_lists_the_shell_model_tools(tmp_path):
    server = build_server(tmp_path)

    async def talk():
        async with Client(server) as c:
            listed = {t.name: t for t in (await c.list_tools()).tools}
            assert {name for name, _ in ShellModelTools.TOOLS} <= set(listed)
            assert {name for name, _ in SetupTools.TOOLS} <= set(listed)
            assert listed["shell_model_options"].annotations.read_only_hint
            assert listed["shell_model_results"].annotations.read_only_hint
            assert not listed["run_shell_model"].annotations.read_only_hint
            assert listed["apply_shell_model"].annotations.destructive_hint
            assert "explicitly" in listed["apply_shell_model"].description
            assert "explicitly" in listed["run_shell_model"].description
            assert "half of the RAM free" in listed["run_shell_model"].description
            assert "states" in listed["estimate_shell_model"].input_schema["properties"]

            out = (await c.call_tool("import_shell_model", {"path": str(TABLE)}))
            assert not out.is_error, out.content[0].text
            refused = await c.call_tool("estimate_shell_model", {"nuclide": "28Si", "interaction": "usdb",
                                                                 "states": "2"})
            assert refused.is_error and "cannot read '2'" in refused.content[0].text

    asyncio.run(talk())
    assert "Shell-model calculations" in server.instructions


def test_shell_model_tools_over_stdio(tmp_path):
    # The whole environment, so the server imports this physim (PYTHONPATH) and uses the test's cache.
    params = StdioServerParameters(command=sys.executable, env=dict(os.environ),
                                   args=["-m", "physim", "mcp", "--experiments", str(tmp_path)])

    async def talk():
        async with Client(params) as c:
            names = {t.name for t in (await c.list_tools()).tools}
            assert {name for name, _ in ShellModelTools.TOOLS} <= names
            out = await c.call_tool("shell_model_options", {"nuclide": "28Si"})
            assert not out.is_error and "usdb" in out.content[0].text

    asyncio.run(talk())


# -- a real run --------------------------------------------------------------------------------------------------


def test_si28_usdb_for_real_through_the_tools(tools, monkeypatch):
    monkeypatch.undo()  # the real free memory (the cache fixture faked it); a cache of this test's own below
    monkeypatch.setenv("PHYSIM_SHELLMODEL", str(Path(tools.session.root).parent / "shellmodel-real"))
    if kshell.find_install() is None:
        pytest.skip("KSHELL is not installed (set PHYSIM_KSHELL)")
    out = tools.run_shell_model("28Si", "usdb", SI28_STATES)
    job = out["job"]
    if not out.get("cached"):
        for _ in range(10):  # up to 10 minutes; the run takes about half a minute on a desktop
            st = tools.shell_model_status(job, wait_seconds=60)
            if st["state"] not in ("queued", "running"):
                break
        assert st["state"] == "done", st
    res = tools.shell_model_results(job)
    ref = shellmodel.read_folder(LOGS)
    for want in ref.states:
        got = next(s for s in res["states"] if s["state"] == want.label)
        assert got["ex_kev"] == pytest.approx(1000 * want.ex_mev, abs=1.0)
    e2 = next(t for t in res["transitions"] if t["kind"] == "E2" and (t["from"], t["to"]) == ("2+1", "0+1"))
    assert e2["b_down_e2fm4"] == pytest.approx(ref.transition("E2", (4, 1, 1), (0, 1, 1)).b_down, rel=1e-3)
