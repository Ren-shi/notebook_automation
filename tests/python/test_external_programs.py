"""physim.external.programs: where an outside program lives (environment variable, ~/.physim/programs.toml, the
copy bundled with the installed planner, the repository's tools/ folder), checked on temporary folders."""

from pathlib import Path

import pytest

from physim.external import programs


@pytest.fixture
def places(tmp_path, monkeypatch):
    """A fake home, install folder and repository, each with its own tools/kshell/, and no PHYSIM_KSHELL."""
    home, install, repo = tmp_path / "home", tmp_path / "install", tmp_path / "repo"
    for root in (install, repo):
        (root / "tools" / "kshell").mkdir(parents=True)
    home.mkdir()
    monkeypatch.delenv("PHYSIM_KSHELL", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(programs, "install_folder", lambda: install)
    monkeypatch.setattr(programs, "repo_folder", lambda: repo)
    return {"home": home, "install": install, "repo": repo, "tmp": tmp_path}


def test_lookup_order(places, monkeypatch):
    install, repo, tmp = places["install"], places["repo"], places["tmp"]
    assert programs.find("kshell") == (install / "tools" / "kshell").resolve()  # bundled before the repository

    (install / "tools" / "kshell").rmdir()
    assert programs.find("kshell") == (repo / "tools" / "kshell").resolve()

    own = tmp / "my kshell"
    own.mkdir()
    programs.set_path("kshell", own)
    assert programs.find("kshell") == own.resolve()  # the settings file before the bundled copy

    env = tmp / "cluster build"
    env.mkdir()
    monkeypatch.setenv("PHYSIM_KSHELL", str(env))
    assert programs.find("kshell") == env.resolve()  # the environment variable before everything

    assert [where for where, _ in programs.candidates("kshell")] == [
        "PHYSIM_KSHELL", str(programs.settings_file()), "bundled", "repository"]


def test_missing_places_are_passed_over(places, monkeypatch):
    tmp = places["tmp"]
    monkeypatch.setenv("PHYSIM_KSHELL", str(tmp / "nowhere"))
    programs.set_path("kshell", tmp / "gone")
    assert programs.find("kshell") == (places["install"] / "tools" / "kshell").resolve()
    for root in (places["install"], places["repo"]):
        (root / "tools" / "kshell").rmdir()
    assert programs.find("kshell") is None
    assert programs.find("gosia") is None


def test_a_setting_naming_the_exe_stands_for_its_folder(places, monkeypatch):
    exe = places["tmp"] / "build" / "kshell.exe"
    exe.parent.mkdir()
    exe.write_bytes(b"")
    monkeypatch.setenv("PHYSIM_KSHELL", str(exe))
    assert programs.find("kshell") == exe.parent.resolve()


def test_set_path_keeps_other_programs_and_can_forget(places):
    tmp = places["tmp"]
    a, b = tmp / "a", tmp / "b"
    a.mkdir()
    b.mkdir()
    programs.set_path("kshell", a)
    programs.set_path("gosia", b)
    text = programs.settings_file().read_text(encoding="utf-8")
    assert programs.settings_file() == places["home"] / ".physim" / "programs.toml"
    assert f'kshell = "{a.resolve().as_posix()}"' in text and f'gosia = "{b.resolve().as_posix()}"' in text
    assert programs.find("gosia") == b.resolve()

    programs.set_path("kshell", None)
    assert programs.find("kshell") == (places["install"] / "tools" / "kshell").resolve()
    assert programs.find("gosia") == b.resolve()


def test_a_broken_settings_file_is_passed_over(places):
    path = programs.settings_file()
    path.parent.mkdir(parents=True)
    path.write_text("kshell = = nonsense", encoding="utf-8")
    assert programs.find("kshell") == (places["install"] / "tools" / "kshell").resolve()


def test_names_are_checked():
    for bad in ("", "../kshell", "k shell", "1kshell"):
        with pytest.raises(ValueError):
            programs.find(bad)


def test_install_and_repo_folders_are_recognised(tmp_path, monkeypatch):
    installed = tmp_path / "physim planner" / "python" / "Lib" / "site-packages" / "physim"
    installed.mkdir(parents=True)
    monkeypatch.setattr(programs, "_package_folder", lambda: installed)
    assert programs.install_folder() == tmp_path / "physim planner"
    assert programs.repo_folder() is None

    checkout = tmp_path / "checkout" / "python" / "physim"
    checkout.mkdir(parents=True)
    monkeypatch.setattr(programs, "_package_folder", lambda: checkout)
    assert programs.install_folder() is None
    assert programs.repo_folder() is None  # no pyproject.toml: not a checkout
    (tmp_path / "checkout" / "pyproject.toml").write_text("")
    assert programs.repo_folder() == tmp_path / "checkout"


def test_this_checkout_is_found_as_the_repository():
    checkout = Path(__file__).resolve().parents[2]
    if not Path(programs.__file__).resolve().is_relative_to(checkout / "python"):
        pytest.skip("physim is installed, not imported from this checkout")
    assert programs.repo_folder() == checkout and programs.install_folder() is None
