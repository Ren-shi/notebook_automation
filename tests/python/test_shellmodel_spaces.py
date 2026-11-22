"""Shell-model valence spaces, interactions, partition files, dimensions and memory estimates.

The oracles are KSHELL's own tools, run in wave 0 (``tests/reference/nuclear/shellmodel/``): 14 partition files
from ``gen_partition.py`` with no truncation and ``count_dim.exe``'s dimensions for each, the ``.snt`` files, and
the ²⁸Si and ⁴⁸Cr logs (vector counts and run times).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import pytest

from physim.nuclear import shellmodel_spaces as sms

REF = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "shellmodel"
PTN_DIR = REF / "partitions"
SNT_DIR = REF / "snt"
PTN_FILES = sorted(PTN_DIR.glob("*.ptn"))
GB = 1024**3


def _case(ptn: Path):
    """(interaction, z_val, n_val, parity) from a fixture partition file's name and header."""
    _, inter, _ = ptn.stem.split("_")
    m = re.search(r"Z=(\d+)\s+N=(\d+)\s+parity=([+-]1)", ptn.read_text())
    return inter, int(m.group(1)), int(m.group(2)), int(m.group(3))


def _count_dim(ptn: Path) -> tuple[dict[int, int], dict[int, int]]:
    """2M → M-scheme dimension and 2M → J-scheme dimension from count_dim.exe's output."""
    m_dims, j_dims = {}, {}
    for line in ptn.with_suffix(".count_dim.txt").read_text().splitlines():
        words = line.split()
        if len(words) >= 4 and words[0] == "dim.":
            m_dims[int(words[1])] = int(words[2])
            j_dims[int(words[1])] = int(words[3])
    return m_dims, j_dims


def test_fixtures_present():
    assert len(PTN_FILES) == 14


# --- spaces and interactions ------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["ckpot", "usdb", "kb3g", "gxpf1a", "jun45", "sdpf-mu"])
def test_built_in_space_matches_snt(name):
    spec = sms.INTERACTIONS[name]
    space = sms.read_space(SNT_DIR / spec.snt)
    assert space.name == name
    assert (space.core_z, space.core_n) == (spec.space.core_z, spec.space.core_n)
    assert space.orbits == spec.space.orbits


def test_read_space_details():
    space = sms.read_space(SNT_DIR / "sdpf-mu.snt")
    assert (space.core_z, space.core_n) == (8, 8)
    assert len(space.proton_orbits) == len(space.neutron_orbits) == 7
    assert space.orbits[0] == sms.Orbit(0, 2, 5, -1)
    assert space.orbits[-1] == sms.Orbit(1, 1, 1, 1)
    assert space.orbits[3].label == "p 0f7/2"
    assert space.capacity(-1) == space.capacity(1) == 32
    assert space.core_label == "16O"
    # kb3.snt is not in the table but uses the same fp space as kb3g
    assert sms.read_space(SNT_DIR / "kb3.snt").orbits == sms.FP_SPACE.orbits


def test_read_space_rejects_garbage(tmp_path):
    bad = tmp_path / "bad.snt"
    bad.write_text("! nothing\n 1 1 8 8\n 1 0 2 5 1\n")
    with pytest.raises(ValueError, match="bad.snt"):
        sms.read_space(bad)


def test_snt_namelist_defaults():
    assert sms.read_snt_defaults(SNT_DIR / "gxpf1a.snt")["eff_charge"] == (1.5, 0.5)
    jun = sms.read_snt_defaults(SNT_DIR / "jun45.snt")
    assert jun["eff_charge"] == (1.5, 1.1)
    assert jun["gl"] == (1.0, 0.0)
    assert jun["gs"] == (3.910, -2.678)
    mu = sms.read_snt_defaults(SNT_DIR / "sdpf-mu.snt")
    assert mu["eff_charge"] == (1.35, 0.35) and mu["hw_type"] == 2
    assert sms.read_snt_defaults(SNT_DIR / "usdb.snt") == {}


@pytest.mark.parametrize("name", ["ckpot", "usdb", "kb3g", "gxpf1a", "jun45", "sdpf-mu"])
def test_interaction_defaults_follow_the_snt(name):
    """Where the .snt file recommends inputs, the table uses them; elsewhere 1.5/0.5 and free g-factors."""
    spec = sms.INTERACTIONS[name]
    snt = sms.read_snt_defaults(SNT_DIR / spec.snt)
    assert spec.eff_charge == snt.get("eff_charge", (1.5, 0.5))
    assert spec.gl == snt.get("gl", (1.0, 0.0))
    assert spec.gs == snt.get("gs", (5.585, -3.826))
    assert spec.hw_type == snt.get("hw_type", spec.hw_type)
    assert spec.citation and spec.defaults_from


def test_table_covers_the_plan():
    assert set(sms.INTERACTIONS) == {"ckpot", "usda", "usdb", "kb3g", "gxpf1a", "jun45", "sdpf-mu"}
    assert "USDA" in sms.INTERACTIONS["usdb"].notes          # the header's mislabel is recorded
    assert sms.interaction("USDB.snt") is sms.INTERACTIONS["usdb"]
    with pytest.raises(ValueError, match="unknown interaction"):
        sms.interaction("nosuch")


def test_hw_type_and_beta_cm_as_kshells_helper_chooses():
    hw = {n: i.hw_type for n, i in sms.INTERACTIONS.items()}
    assert hw == {"ckpot": 1, "usda": 2, "usdb": 2, "kb3g": 1, "gxpf1a": 1, "jun45": 1, "sdpf-mu": 2}
    beta = {n for n, i in sms.INTERACTIONS.items() if i.beta_cm}
    assert beta == {"sdpf-mu"}
    assert sms.INTERACTIONS["sdpf-mu"].beta_cm == 10.0
    # the ²⁸Si input KSHELL's helper wrote in wave 0
    text = (REF / "kshell" / "si28_usdb" / "inputs" / "Si28_usdb_0.input").read_text()
    assert re.search(r"hw_type\s*=\s*2\b", text) and re.search(r"beta_cm\s*=\s*0\.0\b", text)


# --- options ----------------------------------------------------------------------------------------------------


def test_options_si28():
    opts = {o.interaction: o for o in sms.options("28Si")}
    assert set(opts) == {"usda", "usdb", "sdpf-mu"}
    assert (opts["usdb"].z_val, opts["usdb"].n_val, opts["usdb"].parities) == (6, 6, (1,))
    assert opts["usdb"].dimension_by_parity == {1: 93710}
    assert opts["sdpf-mu"].parities == (1, -1)
    assert opts["usdb"].spec is sms.INTERACTIONS["usdb"]


def test_options_other_nuclei():
    # the fixtures named c12/c13 hold 2+2 and 2+3 valence nucleons above 4He: 8Be and 9Be
    assert [(o.interaction, o.z_val, o.n_val, o.parities) for o in sms.options("9Be")] == [("ckpot", 2, 3, (-1,))]
    assert [(o.interaction, o.z_val, o.n_val, o.parities) for o in sms.options("13C")] == [("ckpot", 4, 5, (-1,))]
    cr48 = {o.interaction: o for o in sms.options("Cr-48")}
    # 48Cr also fits (16 + 20 nucleons) in the sd-pf space; its dimension says how hopeless that is
    assert set(cr48) == {"kb3g", "gxpf1a", "sdpf-mu"}
    assert cr48["kb3g"].parities == (1,) and cr48["kb3g"].dimension_by_parity == {1: 1_963_461}
    assert cr48["sdpf-mu"].dimension_by_parity[1] > 1e12
    # 64Ge: 4 + 4 above 56Ni in f5pg9, 12 + 12 in the fp shell, 24 + 24 in sd-pf
    assert {o.interaction for o in sms.options("64Ge")} == {"kb3g", "gxpf1a", "jun45", "sdpf-mu"}
    assert sms.options("208Pb") == []
    assert all(o.interaction != "usdb" for o in sms.options("16O"))   # no valence nucleons
    with pytest.raises(ValueError):
        sms.options("not a nuclide")


def test_parities():
    assert sms.parities(sms.SDPF_SPACE, 0, 6) == (1, -1)
    assert sms.parities(sms.P_SPACE, 2, 2) == (1,)
    assert sms.parities(sms.P_SPACE, 2, 3) == (-1,)
    assert sms.parities(sms.F5PG9_SPACE, 4, 4) == (1, -1)
    assert sms.parities(sms.SD_SPACE, 13, 0) == ()


# --- dimensions and partitions against KSHELL -------------------------------------------------------------------


@pytest.mark.parametrize("ptn", PTN_FILES, ids=lambda p: p.stem)
def test_dimensions_match_count_dim(ptn):
    inter, z_val, n_val, parity = _case(ptn)
    space = sms.INTERACTIONS[inter].space
    m_dims, j_dims = _count_dim(ptn)
    assert m_dims
    dims = sms.dimensions(space, z_val, n_val, parity)
    # count_dim also prints a 0 row at the largest 2M of the other parity (ge64 negative parity, 2M = 52)
    assert dims == {m2: d for m2, d in m_dims.items() if d}
    assert all(dims.get(m2, 0) == d for m2, d in m_dims.items())
    for m2, j in j_dims.items():
        assert sms.j_dimension(space, z_val, n_val, parity, m2) == j


def test_ni56_full_fp_dimension():
    dims = sms.dimensions(sms.FP_SPACE, 8, 8, 1)
    assert dims[0] == 1_087_455_228
    assert sms.dimensions(sms.FP_SPACE, 8, 8, -1) == {}


@pytest.mark.parametrize("ptn", PTN_FILES, ids=lambda p: p.stem)
def test_partitions_match_gen_partition(ptn):
    inter, z_val, n_val, parity = _case(ptn)
    spec = sms.INTERACTIONS[inter]
    text = sms.partitions(spec.space, z_val, n_val, parity, snt_name=spec.snt)
    expected = ptn.read_text()   # the fixtures have Windows line endings
    assert text.splitlines() == expected.splitlines()


def _ptn_sets(text: str):
    """Proton partitions, neutron partitions and pairs (as occupation tuples) of a .ptn text."""
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    n_p, n_n = (int(w) for w in lines[1].split())
    prot = [tuple(int(w) for w in ln.split()[1:]) for ln in lines[2:2 + n_p]]
    neut = [tuple(int(w) for w in ln.split()[1:]) for ln in lines[2 + n_p:2 + n_p + n_n]]
    n_pairs = int(lines[2 + n_p + n_n])
    pairs = {(prot[int(a) - 1], neut[int(b) - 1])
             for a, b in (ln.split() for ln in lines[3 + n_p + n_n:3 + n_p + n_n + n_pairs])}
    assert len(pairs) == n_pairs
    return set(prot), set(neut), pairs


@pytest.mark.parametrize("ptn", PTN_FILES, ids=lambda p: p.stem)
def test_partitions_match_as_sets(ptn):
    inter, z_val, n_val, parity = _case(ptn)
    spec = sms.INTERACTIONS[inter]
    assert _ptn_sets(sms.partitions(spec.space, z_val, n_val, parity)) == _ptn_sets(ptn.read_text())


def test_partitions_default_header_and_errors():
    assert sms.partitions(sms.P_SPACE, 2, 2, 1).startswith("# partition file of p.snt  Z=2  N=2  parity=+1\n")
    with pytest.raises(ValueError, match="no - parity"):
        sms.partitions(sms.SD_SPACE, 2, 2, -1)
    with pytest.raises(ValueError, match="cannot hold"):
        sms.partitions(sms.SD_SPACE, 13, 2, 1)


# --- memory -----------------------------------------------------------------------------------------------------


def test_available_memory_is_read_from_the_os():
    free = sms.available_memory_bytes()
    assert isinstance(free, int) and free > 64 * 1024**2


def test_lanczos_parameters_as_kshells_helper():
    assert sms.lanczos_parameters(3) == (10, 200)
    assert sms.lanczos_parameters(32) == (48, 200)
    assert sms.lanczos_parameters(200) == (300, 350)


@pytest.mark.parametrize("log, n_eigen", [("log_Si28_usdb_j0p.txt", 3), ("log_Si28_usdb_j8p.txt", 2)])
def test_vectors_on_disk_mode_match_kshell_logs(log, n_eigen):
    text = (REF / "kshell" / "si28_usdb" / "logs" / log).read_text(errors="replace")
    count = int(re.search(r"Memory / process is:.*?x\s+(\d+)\s+=", text).group(1))
    assert sms.vectors_in_memory(n_eigen, 1) == count


def test_vectors_in_memory_mode_cover_kshell_log():
    """⁴⁸Cr, 32 states, vectors in memory: KSHELL held 216 vectors; the estimate must not count fewer."""
    text = (REF / "kshell" / "cr48_kb3" / "log_cr48_m0.txt").read_text(errors="replace")
    count = int(re.search(r"Memory / process is:.*?x\s+(\d+)\s+=", text).group(1))
    assert count == 216
    assert sms.vectors_in_memory(32, 0) >= count


@dataclass(frozen=True)
class FakeRequest:
    """Stands in for session 1C's ShellModelRequest (only the attributes the estimate reads)."""

    nuclide: str
    interaction: str
    states: tuple
    threads: int | None = 12


SI28 = FakeRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3), (8, 1, 2)))
CR48 = FakeRequest("48Cr", "kb3g", ((0, 1, 1),))
NI56 = FakeRequest("56Ni", "gxpf1a", ((0, 1, 1),))


def test_limit_is_half_of_free_memory_users_example(monkeypatch):
    """16 GB installed, 10 GB in use: 6 GB free, so a run may use 3 GB."""
    monkeypatch.setattr(sms, "available_memory_bytes", lambda: 16 * GB - 10 * GB)
    est = sms.estimate(SI28)
    assert est.free_bytes == 6 * GB
    assert est.limit_bytes == 3 * GB
    assert est.fits and est.mode_lv_hdd == 0
    assert "3 GB" in est.reason and "6 GB" in est.reason


def test_si28_estimate():
    est = sms.estimate(SI28, free_bytes=6 * GB)
    assert est.dimension_by_state == {(0, 1): 93710, (4, 1): 81122, (8, 1): 52175}
    assert est.interaction == "usdb" and (est.z_val, est.n_val) == (6, 6)
    # 204 + 3 vectors of 93 710 doubles, times 1.2, plus the fixed overhead
    assert est.memory_bytes == int(93710 * 8 * 206 * 1.2) + sms.FIXED_OVERHEAD_BYTES
    assert est.disk_bytes == (3 * 93710 + 3 * 81122 + 2 * 52175) * 8
    assert est.memory_bytes_on_disk_mode < est.memory_bytes_in_memory_mode
    # measured 27 s for the whole set (12 threads): the rough guess is within a factor of two
    assert 13 < est.seconds_guess < 54
    json.dumps(est.to_dict())
    assert est.to_dict()["states"][1] == {"j2": 4, "parity": 1, "label": "2+", "dimension": 81122}


def test_cr48_goes_to_disk_under_3_gb_and_memory_with_more():
    tight = sms.estimate(CR48, free_bytes=6 * GB)
    assert tight.dimension_by_state == {(0, 1): 1_963_461}
    assert tight.memory_bytes_in_memory_mode > 3 * GB          # 200 Lanczos vectors of 15.7 MB do not fit
    assert tight.fits and tight.mode_lv_hdd == 1
    assert tight.disk_bytes >= 200 * 1_963_461 * 8
    assert "on disk" in tight.reason
    roomy = sms.estimate(CR48, free_bytes=16 * GB)
    assert roomy.fits and roomy.mode_lv_hdd == 0 and roomy.limit_bytes == 8 * GB
    assert 36 < roomy.seconds_guess < 144                         # measured 72 s


def test_ni56_does_not_fit():
    est = sms.estimate(NI56, free_bytes=6 * GB)
    assert est.dimension_by_state == {(0, 1): 1_087_455_228}
    assert not est.fits and est.mode_lv_hdd == 1
    assert est.memory_bytes > 30 * GB
    assert "Too big" in est.reason and "cluster" in est.reason


def test_threads_scale_the_time():
    one = sms.estimate(FakeRequest("48Cr", "kb3g", ((0, 1, 1),), threads=1), free_bytes=16 * GB)
    twelve = sms.estimate(CR48, free_bytes=16 * GB)
    assert 4 < one.seconds_guess / twelve.seconds_guess < 6.5     # measured 401 s / 72 s = 5.6


def test_estimate_notes_missing_states_and_rejects_bad_requests():
    est = sms.estimate(FakeRequest("28Si", "usdb", ((28, 1, 2), (30, 1, 1))), free_bytes=6 * GB)
    assert est.dimension_by_state == {(28, 1): 1, (30, 1): 0}
    assert "14+ (1 exists, 2 asked)" in est.reason and "15+ (0 exist, 1 asked)" in est.reason
    with pytest.raises(ValueError, match="2J must be even"):
        sms.estimate(FakeRequest("28Si", "usdb", ((1, 1, 1),)), free_bytes=GB)
    with pytest.raises(ValueError, match="does not fit"):
        sms.estimate(FakeRequest("208Pb", "usdb", ((0, 1, 1),)), free_bytes=GB)
    with pytest.raises(ValueError, match="unknown interaction"):
        sms.estimate(FakeRequest("28Si", "nosuch", ((0, 1, 1),)), free_bytes=GB)
    with pytest.raises(ValueError, match="no states"):
        sms.estimate(FakeRequest("28Si", "usdb", ()), free_bytes=GB)


def test_estimate_accepts_an_interaction_object_and_negative_parity():
    req = FakeRequest("22O", sms.INTERACTIONS["sdpf-mu"], ((2, -1, 2), (0, 1, 1)))
    est = sms.estimate(req, free_bytes=6 * GB)
    m_neg, _ = _count_dim(PTN_DIR / "o22_sdpf-mu_n.ptn")
    m_pos, _ = _count_dim(PTN_DIR / "o22_sdpf-mu_p.ptn")
    assert est.dimension_by_state == {(2, -1): m_neg[2], (0, 1): m_pos[0]}
    assert est.interaction == "sdpf-mu" and (est.z_val, est.n_val) == (0, 6)


def test_an_effective_charge_with_no_valence_particles_is_warned_about():
    """62Ni in JUN45 (56Ni core): 0 valence protons, 6 neutrons, so e_p has no effect."""
    (warn,) = sms.charge_warnings("62Ni", "jun45")
    assert warn.startswith("e_p = 1.5 e (the interaction's default) has no effect") and "56Ni core" in warn
    (given,) = sms.charge_warnings("62Ni", "jun45", (1.3, 0.6))
    assert given.startswith("e_p = 1.3 e (given)")
    assert sms.charge_warnings("62Ni", "jun45", (0.0, 0.6)) == []        # a zero charge is harmless
    assert sms.charge_warnings("28Si", "usdb") == []                     # both species have valence particles
    (o22,) = sms.charge_warnings("22O", "usdb")                           # 0 protons above 16O
    assert o22.startswith("e_p") and "no valence protons" in o22
    est = sms.estimate(FakeRequest("62Ni", "jun45", ((0, 1, 1),)), free_bytes=6 * GB)
    assert est.warnings == sms.charge_warnings("62Ni", "jun45") and est.to_dict()["warnings"] == est.warnings
    assert sms.estimate(SI28, free_bytes=6 * GB).warnings == []
