"""The pure-Python parts of scripts/build_kshell.py (checksums, the rename patch, DLL discovery, the log comparison)
and the KSHELL bundling in installer/windows/build.py. Nothing here compiles or runs KSHELL."""

from __future__ import annotations

import importlib.util
import io
import sys
import tarfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
LOGS = REPO / "tests" / "reference" / "nuclear" / "shellmodel" / "kshell" / "si28_usdb" / "logs"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # dataclasses look their module up there
    spec.loader.exec_module(module)
    return module


bk = _load("build_kshell", REPO / "scripts" / "build_kshell.py")
installer = _load("installer_build", REPO / "installer" / "windows" / "build.py")


# --------------------------------------------------------------------------------------------- checksums


def _archive_for(path: Path, md5: bool = True):
    return bk.Archive(path.name, "https://example.invalid/x", bk.file_digest(path, "sha256"),
                      bk.file_digest(path, "md5") if md5 else None)


def test_verify_accepts_matching_checksums(tmp_path):
    f = tmp_path / "a.tar.gz"
    f.write_bytes(b"kshell")
    bk.verify(f, _archive_for(f))


def test_verify_rejects_wrong_sha256_and_md5(tmp_path):
    f = tmp_path / "a.tar.gz"
    f.write_bytes(b"kshell")
    good = _archive_for(f)
    with pytest.raises(ValueError, match="SHA-256"):
        bk.verify(f, bk.Archive(f.name, good.url, "0" * 64))
    with pytest.raises(ValueError, match="MD5"):
        bk.verify(f, bk.Archive(f.name, good.url, good.sha256, "0" * 32))


def test_find_archive_searches_subfolders_and_skips_bad_copies(tmp_path):
    (tmp_path / "bad").mkdir()
    (tmp_path / "good" / "deeper").mkdir(parents=True)
    good = tmp_path / "good" / "deeper" / "kshell.tar.gz"
    good.write_bytes(b"right")
    (tmp_path / "bad" / "kshell.tar.gz").write_bytes(b"wrong")
    assert bk.find_archive(_archive_for(good), [tmp_path]) == good
    assert bk.find_archive(bk.Archive("kshell.tar.gz", "", "0" * 64), [tmp_path]) is None


def test_pinned_checksums_are_well_formed():
    assert len(bk.KSHELL.sha256) == 64 and bk.KSHELL.sha256.startswith("ce0659c3")
    assert len(bk.OPENBLAS.sha256) == 64 and bk.OPENBLAS.md5 == "5047115a08a7d39471375a25994108c5"


def test_the_vendored_release_is_the_pinned_one():
    assert bk.VENDORED_KSHELL.is_file()
    bk.verify(bk.VENDORED_KSHELL, bk.KSHELL)
    assert bk.get_archive(bk.KSHELL, None, Path("unused"), bk.VENDORED_KSHELL) == bk.VENDORED_KSHELL


def test_archives_folder_overrides_the_vendored_copy(tmp_path):
    local = tmp_path / "local" / "cpc"
    local.mkdir(parents=True)
    copy = local / "kshell.tar.gz"
    copy.write_bytes(bk.VENDORED_KSHELL.read_bytes())
    assert bk.get_archive(bk.KSHELL, tmp_path / "local", tmp_path / "dl", bk.VENDORED_KSHELL) == copy
    # a folder without it falls back to the vendored copy; without that either, it is an error, not a download
    assert bk.get_archive(bk.KSHELL, tmp_path / "dl", tmp_path / "dl", bk.VENDORED_KSHELL) == bk.VENDORED_KSHELL
    with pytest.raises(SystemExit, match="not found under"):
        bk.get_archive(bk.KSHELL, tmp_path / "dl", tmp_path / "dl", tmp_path / "missing.tar.gz")


def test_a_changed_vendored_copy_is_an_error(tmp_path):
    bad = tmp_path / "kshell.tar.gz"
    bad.write_bytes(b"\x1f\x8b not the release")
    with pytest.raises(ValueError, match="SHA-256"):
        bk.get_archive(bk.KSHELL, None, tmp_path / "dl", bad)


def test_a_web_page_instead_of_the_archive_fails_clearly(tmp_path, monkeypatch):
    class Page(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(bk.urllib.request, "urlopen", lambda req, timeout: Page(b"<!DOCTYPE html><html>"))
    with pytest.raises(SystemExit, match="not a .gz archive .* probably a web page"):
        bk.get_archive(bk.KSHELL, None, tmp_path / "dl", tmp_path / "missing.tar.gz")
    assert not (tmp_path / "dl" / "kshell.tar.gz").exists()  # nothing half-downloaded is left behind

    def refuse(req, timeout):
        raise bk.urllib.error.HTTPError(req.full_url, 403, "Forbidden", {}, None)

    monkeypatch.setattr(bk.urllib.request, "urlopen", refuse)
    with pytest.raises(SystemExit, match="403"):
        bk.download(bk.KSHELL, tmp_path / "dl")


def test_unpack_tar_refuses_paths_outside_the_folder(tmp_path):
    archive = tmp_path / "evil.tar"
    with tarfile.open(archive, "w") as tar:
        info = tarfile.TarInfo("../evil.txt")
        info.size = 1
        tar.addfile(info, io.BytesIO(b"x"))
    with pytest.raises(ValueError, match="unexpected member"):
        bk.unpack_tar(archive, tmp_path / "out")
    assert not (tmp_path / "evil.txt").exists()


# --------------------------------------------------------------------------------------------- the patch

SAMPLE = """module lanczos
  implicit none
contains
  subroutine restart(fn)
    character(len=*) :: fn
    call rename( trim(fn)//'1', &
         trim(fn)//'2' )
    CALL RENAME(trim(fn)//'3', trim(fn)//'4')
    call renamed_thing(fn)
  end subroutine restart
end module lanczos
"""


def test_patch_replaces_every_rename_call_and_adds_the_helper():
    text, n = bk.patch_rename_text(SAMPLE, "lanczos")
    assert n == 2
    assert text.count("call rename_replace(") == 2
    assert "call renamed_thing(fn)" in text  # other names starting with "rename" are left alone
    helper = text.index("subroutine rename_replace(old, new)")
    assert text.index("end subroutine restart") < helper < text.index("end module lanczos")
    assert "call rename(old, new)" in text  # the helper itself calls the intrinsic once
    # POSIX semantics: renaming a file to itself, or a missing file, does nothing
    assert "if (trim(old) == trim(new)) return" in text
    assert "if (.not. ex) return" in text


def test_patch_is_idempotent():
    once, _ = bk.patch_rename_text(SAMPLE, "lanczos")
    twice, n = bk.patch_rename_text(once, "lanczos")
    assert (twice, n) == (once, 0)


def test_patch_refuses_unexpected_sources():
    with pytest.raises(ValueError, match="end module bp_io"):
        bk.patch_rename_text(SAMPLE, "bp_io")
    with pytest.raises(ValueError, match="no 'call rename"):
        bk.patch_rename_text("module bp_io\ncontains\nend module bp_io\n", "bp_io")


def test_patch_sources_writes_files_and_a_diff(tmp_path):
    (tmp_path / "lanczos.f90").write_text(SAMPLE)
    (tmp_path / "bp_io.F90").write_text(SAMPLE.replace("lanczos", "bp_io"))
    counts, diff = bk.patch_sources(tmp_path)
    assert counts == {"lanczos.f90": 2, "bp_io.F90": 2}
    assert "--- a/kshell/src/lanczos.f90" in diff and "+++ b/kshell/src/bp_io.F90" in diff
    assert "+    call rename_replace( trim(fn)//'1', &" in diff
    assert "rename_replace" in (tmp_path / "bp_io.F90").read_text()


# --------------------------------------------------------------------------------------------- DLLs

OBJDUMP = """
kshell.exe:     file format pei-x86-64

The Import Tables (interpreted .idata section contents)
 vma:            Hint    Time      Forward  DLL       First
\tDLL Name: libgcc_s_seh-1.dll
\tvma:  Hint/Ord Member-Name Bound-To
\t3c8a8     115  __register_frame_info
\tDLL Name: libgfortran-5.dll
\tDLL Name: KERNEL32.dll
\tDLL Name: api-ms-win-crt-heap-l1-1-0.dll
\tDLL Name: libopenblas.dll
\tDLL Name: libgfortran-5.dll
"""


def test_parse_objdump_imports():
    assert bk.parse_objdump_imports(OBJDUMP) == [
        "libgcc_s_seh-1.dll", "libgfortran-5.dll", "KERNEL32.dll", "api-ms-win-crt-heap-l1-1-0.dll",
        "libopenblas.dll"]
    assert bk.parse_objdump_imports("no imports here") == []


def _fake_tree(tmp_path):
    mingw, blas = tmp_path / "mingw" / "bin", tmp_path / "openblas" / "bin"
    mingw.mkdir(parents=True)
    blas.mkdir(parents=True)
    imports = {
        "kshell.exe": ["libgfortran-5.dll", "libgomp-1.dll", "KERNEL32.dll", "api-ms-win-crt-heap-l1-1-0.dll",
                       "libopenblas.dll"],
        "transit.exe": ["libgfortran-5.dll", "libopenblas.dll"],
        "libgfortran-5.dll": ["libquadmath-0.dll", "libgcc_s_seh-1.dll", "ADVAPI32.dll"],
        "libgomp-1.dll": ["libdl.dll", "libwinpthread-1.dll"],
        "libquadmath-0.dll": ["libgcc_s_seh-1.dll"],
        "libgcc_s_seh-1.dll": ["libwinpthread-1.dll"],
        "libopenblas.dll": ["msvcrt.dll"],
    }
    for name in ["libgfortran-5.dll", "libgomp-1.dll", "libquadmath-0.dll", "libgcc_s_seh-1.dll", "libdl.dll",
                 "libwinpthread-1.dll"]:
        (mingw / name).write_bytes(b"")
    (blas / "libopenblas.dll").write_bytes(b"")
    system = {"kernel32.dll", "advapi32.dll", "msvcrt.dll"}
    return mingw, blas, (lambda p: imports.get(Path(p).name, [])), (lambda n: n.lower() in system)


def test_collect_dlls_follows_imports_recursively(tmp_path):
    mingw, blas, imports, system = _fake_tree(tmp_path)
    found = bk.collect_dlls([Path("kshell.exe"), Path("transit.exe")], imports, [blas, mingw], system)
    assert sorted(found) == ["libdl.dll", "libgcc_s_seh-1.dll", "libgfortran-5.dll", "libgomp-1.dll",
                             "libopenblas.dll", "libquadmath-0.dll", "libwinpthread-1.dll"]
    assert found["libopenblas.dll"] == blas / "libopenblas.dll"
    assert found["libgomp-1.dll"] == mingw / "libgomp-1.dll"


def test_collect_dlls_reports_a_missing_dll(tmp_path):
    mingw, blas, imports, system = _fake_tree(tmp_path)
    (mingw / "libdl.dll").unlink()
    with pytest.raises(FileNotFoundError, match="libdl.dll .needed by libgomp-1.dll"):
        bk.collect_dlls([Path("kshell.exe")], imports, [blas, mingw], system)


def test_api_sets_are_never_shipped():
    assert bk.is_api_set("api-ms-win-crt-stdio-l1-1-0.dll")
    assert bk.is_api_set("EXT-MS-WIN-foo.dll")
    assert not bk.is_api_set("libgomp-1.dll")


# --------------------------------------------------------------------------------------------- licences


def test_leading_comment_extracts_the_licence_header():
    header = ("/***************\n * Copyright (c) 2025, The OpenBLAS Project\n * All rights reserved.\n"
              " * Redistribution and use in source and binary forms\n * ****/\n#ifndef CBLAS_H\n/* later */\n")
    text = bk.leading_comment(header)
    assert text.startswith("Copyright (c) 2025, The OpenBLAS Project\nAll rights reserved.\n")
    assert "Redistribution and use" in text and "later" not in text and "*" not in text


def test_gpl_text_is_recognised(tmp_path):
    assert bk.is_gpl3("  GNU GENERAL PUBLIC LICENSE\n  Version 3, 29 June 2007\n...")
    assert not bk.is_gpl3("GNU GENERAL PUBLIC LICENSE\nVersion 2, June 1991")
    bad = tmp_path / "gpl2.txt"
    bad.write_text("GNU GENERAL PUBLIC LICENSE\nVersion 2, June 1991")
    with pytest.raises(SystemExit):
        bk.find_gpl_text(bad, None)


def test_the_committed_gpl_text_is_the_default_and_complete():
    assert bk.find_gpl_text(None, None) == bk.GPL_TEXT
    text = bk.GPL_TEXT.read_text(encoding="utf-8")
    assert bk.is_gpl3(text)
    for part in ("END OF TERMS AND CONDITIONS", "How to Apply These Terms to Your New Programs",
                 "Public License instead of this License."):
        assert part in text


# --------------------------------------------------------------------------------------------- log comparison


def _logs():
    return bk.read_logs(LOGS)


def test_parse_fixture_logs():
    logs = _logs()
    e = bk.parse_energies(logs["log_Si28_usdb_j0p.txt"])
    assert sorted(e) == [1, 2, 3]
    assert e[1] == (pytest.approx(-135.86073, abs=1e-5), 0, 1)
    t = bk.parse_transitions(logs["log_Si28_usdb_tr_j0p_j4p.txt"])
    assert t[("E2", 0, 1, 4, 1)] == (499.2211, 99.8442, 0.0)
    m1 = bk.parse_transitions(logs["log_Si28_usdb_tr_j4p_j4p.txt"])
    assert m1[("M1", 4, 1, 4, 1)][2] == 1.0312  # the diagonal line carries the moment
    assert bk.transition_sections(logs["log_Si28_usdb_tr_j0p_j0p.txt"]) == ["E2", "M1"]
    assert bk.parse_transitions(logs["log_Si28_usdb_tr_j0p_j0p.txt"]) == {}


def test_compare_identical_logs_passes():
    logs = _logs()
    assert bk.compare_logs(logs, logs) == []


def test_compare_catches_a_shifted_energy_and_a_changed_b():
    logs = _logs()
    new = dict(logs)
    new["log_Si28_usdb_j0p.txt"] = logs["log_Si28_usdb_j0p.txt"].replace("-135.86073", "-135.85873")  # 2 keV
    new["log_Si28_usdb_tr_j0p_j4p.txt"] = logs["log_Si28_usdb_tr_j0p_j4p.txt"].replace("499.2211", "499.9000")
    problems = bk.compare_logs(new, logs)
    assert any("log_Si28_usdb_j0p.txt: state 1" in p for p in problems)
    assert any("('E2', 0, 1, 4, 1)" in p for p in problems)
    assert len(problems) == 2


def test_compare_tolerates_rounding_and_sign_of_mred():
    logs = _logs()
    new = dict(logs)
    # 0.5 keV, 0.05 % in B, and an opposite phase of the wave function (Mred's sign) are all fine
    new["log_Si28_usdb_j0p.txt"] = logs["log_Si28_usdb_j0p.txt"].replace("-135.86073", "-135.86123")
    new["log_Si28_usdb_tr_j0p_j4p.txt"] = (logs["log_Si28_usdb_tr_j0p_j4p.txt"]
                                           .replace("499.2211", "499.4700").replace("-22.3433", " 22.3433"))
    assert bk.compare_logs(new, logs) == []


def test_compare_reports_missing_logs_and_collapse():
    logs = _logs()
    new = {k: v for k, v in logs.items() if k != "log_Si28_usdb_j8p.txt"}
    new["log_Si28_usdb_j4p.txt"] = logs["log_Si28_usdb_j4p.txt"] + "\n possible JJ collapse\n"
    problems = bk.compare_logs(new, logs)
    assert "log_Si28_usdb_j8p.txt: missing" in problems
    assert any("collapse" in p for p in problems)


# --------------------------------------------------------------------------------------------- installer


def _kshell_folder(tmp_path, *, licence=True):
    folder = tmp_path / "kshell"
    (folder / "source").mkdir(parents=True)
    for name in ["kshell.exe", "transit.exe", "README.txt", "source/kshell.tar.gz", "libopenblas.dll",
                 "LICENSE-OpenBLAS.txt"]:
        (folder / name).write_bytes(b"")
    (folder / "VERSION").write_text("KSHELL 2019 CPC release + physim rename patch 1\nbuilt: today\n")
    if licence:
        (folder / "COPYING").write_text("GNU GENERAL PUBLIC LICENSE\nVersion 3, 29 June 2007\n")
    return folder


def test_installer_bundles_a_complete_kshell(tmp_path):
    source, out = _kshell_folder(tmp_path), tmp_path / "bundle"
    out.mkdir()
    assert installer.kshell_problems(source) == []
    assert installer.bundle_kshell(source, out)
    assert (out / "tools" / "kshell" / "kshell.exe").is_file()
    installer.write_notices(out, True)
    notices = (out / "THIRD-PARTY-NOTICES.txt").read_text(encoding="utf-8")
    assert "GNU General Public License, version 3" in notices
    assert "KSHELL 2019 CPC release + physim rename patch 1" in notices


def test_installer_warns_and_skips_without_kshell_or_its_licence(tmp_path, capsys):
    out = tmp_path / "bundle"
    out.mkdir()
    assert not installer.bundle_kshell(tmp_path / "absent", out)
    assert "KSHELL not bundled" in capsys.readouterr().out
    assert not installer.bundle_kshell(_kshell_folder(tmp_path, licence=False), out)
    assert "COPYING missing" in capsys.readouterr().out
    assert not (out / "tools").exists()
    installer.write_notices(out, False)
    assert "KSHELL" not in (out / "THIRD-PARTY-NOTICES.txt").read_text(encoding="utf-8")
