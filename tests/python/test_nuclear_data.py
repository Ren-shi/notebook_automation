import re
from pathlib import Path

import pytest

from physim.nuclear import Experiment, Quantity, SetupError, data

N_A = 6.02214076e23


# -- atomic masses ----------------------------------------------------------------------------------------------


def test_ame2020_table():
    table = data.nuclides()
    assert len(table) == 3558
    assert data.nuclide("12C").mass_excess_kev == 0.0  # the definition of u
    assert data.nuclide("4He").mass_excess_kev == pytest.approx(2424.91587)
    assert data.nuclide((82, 208)).mass_excess_kev == pytest.approx(-21748.519)  # ± 1.148 keV
    n = data.nuclide("n")
    assert (n.Z, n.A, n.name) == (0, 1, "n") and n.mass_excess_kev == pytest.approx(8071.31806)
    assert not data.nuclide("208Pb").estimated
    assert sum(x.estimated for x in table) > 1000  # far from stability, masses come from systematics
    with pytest.raises(ValueError, match="not in the AME2020 mass table"):
        data.nuclide("150Au")


def test_masses_in_other_units():
    he4 = data.nuclide("4He")
    assert he4.atomic_mass_u == pytest.approx(4.002603254, abs=1e-9)
    assert he4.atomic_mass_mev == pytest.approx(3728.4013, abs=1e-3)
    # The alpha particle mass (CODATA 2018: 3727.3794066 MeV): two electron masses less, plus 79.0 eV binding.
    assert he4.nuclear_mass_mev == pytest.approx(3727.3794, abs=2e-4)
    # The proton (CODATA 2018: 938.27208816 MeV), from the hydrogen atom (13.6 eV binding).
    assert data.nuclide("p").nuclear_mass_mev == pytest.approx(938.272088, abs=1e-6)


# Q-values as tabulated in nuclear data compilations, each quoted to the precision given (MeV).
PUBLISHED_Q = [
    (["2H", "3H"], ["4He", "n"], 17.589),  # D-T fusion
    (["2H", "2H"], ["3He", "n"], 3.269),
    (["2H", "2H"], ["3H", "p"], 4.033),
    (["2H", "3He"], ["4He", "p"], 18.353),
    (["n", "6Li"], ["3H", "4He"], 4.783),
    (["n", "p"], ["2H"], 2.224566),  # deuteron binding energy
    (["238U"], ["234Th", "4He"], 4.270),  # alpha decay
    (["12C", "4He"], ["16O"], 7.162),  # 12C(alpha, gamma)16O
]


@pytest.mark.parametrize("entrance, exit, q", PUBLISHED_Q, ids=lambda x: "+".join(x) if isinstance(x, list) else x)
def test_q_values_match_published_to_1_kev(entrance, exit, q):
    assert data.q_value(entrance, exit) == pytest.approx(q, abs=1e-3)


def test_q_value_checks_conservation_and_excitation():
    assert data.q_value(["4He", "197Au"], ["4He", "197Au"]) == 0.0
    assert data.q_value(["16O", "208Pb"], ["16O", "208Pb"], excitation_mev=6.13) == pytest.approx(-6.13)
    with pytest.raises(ValueError, match="does not conserve charge"):
        data.q_value(["2H", "3H"], ["4He", "p"])
    with pytest.raises(ValueError, match="does not conserve mass number"):
        data.q_value(["2H", "3H"], ["4He"])


# -- elements ---------------------------------------------------------------------------------------------------


def _within_standard_weight(molar_mass, weight):
    interval = re.match(r"^\[([\d.]+),([\d.]+)\]$", weight)
    if interval:
        return float(interval[1]) <= molar_mass <= float(interval[2]), 0.0
    m = re.match(r"^([\d.]+)\((\d+)\)$", weight)
    decimals = len(m[1].split(".")[1]) if "." in m[1] else 0
    sigma = int(m[2]) * 10.0**-decimals
    deviation = abs(molar_mass - float(m[1])) / sigma
    return deviation <= 1.0, deviation


def test_atomic_weights_from_isotopes_agree_with_nist():
    """Abundances (NIST, 2009) x AME2020 masses reproduce the standard atomic weights.

    All elements with a natural composition land inside the published interval or within 1 sigma, except Se and Hg
    (1.5 and 2.4 sigma): their standard atomic weights were revised in 2013, after the isotopic compositions used.
    """
    checked = 0
    for z in range(1, 93):
        el = data.element(z)
        if not el.isotopes:
            continue
        ok, deviation = _within_standard_weight(el.molar_mass, el.standard_atomic_weight)
        if el.symbol in ("Se", "Hg"):
            assert deviation < 3.0, el
        else:
            assert ok, (el.symbol, el.molar_mass, el.standard_atomic_weight)
        checked += 1
    assert checked == 84  # Z = 1-92 except Tc, Pm and Po to Ac


def test_element_properties():
    au = data.element("Au")
    assert (au.Z, au.name, au.density_g_cm3, au.I_eV) == (79, "Gold", 19.32, 790.0)
    assert au.isotopes == ((197, 1.0),)
    assert data.element(14).I_eV == 173.0
    tc = data.element("Tc")
    assert tc.isotopes == ()
    with pytest.raises(ValueError, match="no natural isotopic composition; name an isotope, e.g. '98Tc'"):
        tc.molar_mass
    assert data.element("Og").density_g_cm3 is None  # NIST's table stops at uranium


# -- materials and areal densities (hand calculations) ----------------------------------------------------------


def test_gold():
    au = data.material("Au")
    assert au.areal_density_mg_cm2("1 um") == pytest.approx(1.932)  # 19.32 g/cm3 x 1e-4 cm
    assert au.thickness_um("1.932 mg/cm2") == pytest.approx(1.0)
    atoms = au.atoms_per_cm2("1 mg/cm2")
    assert atoms == {(79, 197): pytest.approx(1e-3 * N_A / 196.96657, rel=1e-6)}  # 3.0574e18


def test_carbon_and_silicon():
    c = data.material("C")
    assert c.density_g_cm3 == 1.7  # NIST's graphite value; real foils vary, so set density when it matters
    atoms = c.atoms_per_cm2("1 mg/cm2")
    assert sum(atoms.values()) == pytest.approx(1e-3 * N_A / 12.0107, rel=1e-4)  # 5.014e19
    assert atoms[(6, 13)] / sum(atoms.values()) == pytest.approx(0.0107)
    si = data.material("Si")
    assert si.areal_density_mg_cm2("300 um") == pytest.approx(2.33 * 0.03 * 1e3)  # 69.9 mg/cm2


def test_polyethylene_and_deuterated_polyethylene():
    ch2 = data.material("CH2")
    assert ch2.density_g_cm3 is None
    assert ch2.molar_mass == pytest.approx(12.0107 + 2 * 1.00794, rel=1e-4)
    by_z = {}
    for (z, _), n in ch2.atoms_per_cm2("1 mg/cm2").items():
        by_z[z] = by_z.get(z, 0) + n
    assert by_z[6] == pytest.approx(1e-3 * N_A / 14.0266, rel=1e-4)  # 4.293e19
    assert by_z[1] == pytest.approx(2 * by_z[6])
    cd2 = data.material("CD2", density="1.06 g/cm3")
    assert cd2.molar_mass == pytest.approx(12.0107 + 2 * 2.014102, rel=1e-4)
    d_atoms = cd2.atoms_per_cm2("1 mg/cm2")[(1, 2)]
    assert d_atoms == pytest.approx(2e-3 * N_A / 16.0389, rel=1e-4)  # 7.509e19
    assert cd2.areal_density_mg_cm2("10 um") == pytest.approx(1.06)
    with pytest.raises(ValueError, match="density of CH2 is not known"):
        ch2.areal_density_mg_cm2("10 um")


def test_named_compounds_and_isotopes():
    mylar = data.material("Mylar")
    assert mylar.density_g_cm3 == 1.38 and mylar.I_eV == 78.7
    by_z = mylar.elements
    # Mylar is C10H8O4: atom ratios from NIST's mass fractions.
    assert by_z[1] / by_z[6] == pytest.approx(0.8, rel=1e-3) and by_z[8] / by_z[6] == pytest.approx(0.4, rel=1e-3)
    assert data.material("CsI").density_g_cm3 == 4.51
    assert "Water, Liquid" in data.compound_names()
    pb = data.material("208Pb")
    assert pb.atoms == ((82, 208, 1.0),)
    # Same atoms per volume as natural lead, so the density scales with the atomic mass.
    assert pb.density_g_cm3 == pytest.approx(11.35 * 207.9767 / 207.2169, rel=1e-4)
    c13 = data.Material.enriched("C", {13: 99, 12: 1})
    assert c13.atoms == ((6, 12, 0.01), (6, 13, 0.99))
    assert c13.density_g_cm3 == pytest.approx(1.7 * (0.01 * 12 + 0.99 * 13.003355) / 12.0107, rel=1e-4)
    with pytest.raises(ValueError, match="no natural isotopic composition"):
        data.material("Tc")


# -- the setup file uses the data -------------------------------------------------------------------------------


def _setup(**target):
    return {
        "schema": "physim.experiment/1", "title": "t",
        "beam": {"nuclide": "4He", "energy": "5 MeV", "current": "1 pnA"},
        "target": {"material": "Au", "thickness": "1 mg/cm2", **target},
        "detectors": [{"shape": "circle", "radius": "5 mm", "theta": "30 deg", "distance": "100 mm",
                       "thickness": "300 um"}],
        "run": {"beam_time": "1 h"},
    }


def test_setup_material_checks():
    exp = Experiment.from_dict(_setup(material="CD2", thickness="10 um", density="1.06 g/cm3"))
    assert exp.target.material_data().areal_density_mg_cm2(exp.target.thickness) == pytest.approx(1.06)
    assert exp.target.density == Quantity(1.06, "g/cm3")
    assert Experiment.from_toml(exp.to_toml()) == exp
    with pytest.raises(SetupError, match="target: thickness '10 um' is a length, but the density of CD2 is not known"):
        Experiment.from_dict(_setup(material="CD2", thickness="10 um"))
    with pytest.raises(SetupError, match="target: material Tc has no natural isotopic composition"):
        Experiment.from_dict(_setup(material="Tc"))
    d = _setup()
    d["beam"]["nuclide"] = "250Au"
    with pytest.raises(SetupError, match="beam: nuclide 250Au is not in the AME2020 mass table"):
        Experiment.from_dict(d)
    d["beam"]["nuclide"] = "n"
    with pytest.raises(SetupError, match="beam: nuclide must be charged"):
        Experiment.from_dict(d)
    si = Experiment.from_dict(_setup()).detectors[0].material_data()
    assert si.areal_density_mg_cm2("300 um") == pytest.approx(69.9)


# -- provenance -------------------------------------------------------------------------------------------------


def test_every_data_file_has_its_source_recorded():
    folder = Path(data.DATA)
    sources = (folder / "SOURCES.md").read_text(encoding="utf-8")
    files = [p.name for p in folder.iterdir() if p.suffix in (".txt", ".csv")]
    assert sorted(files) == ["mass_1.mas20.txt", "nist_astar.csv", "nist_compounds.csv", "nist_elements.csv",
                             "nist_isotopes.csv", "nist_pstar.csv"]
    for name in files:
        assert f"`{name}`" in sources, name
    assert "Creative Commons Attribution 3.0" in sources and "U.S. Secretary of Commerce" in sources
