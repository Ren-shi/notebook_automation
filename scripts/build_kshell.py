"""Build KSHELL (the 2019 release published with Comput. Phys. Commun.) for physim's shell-model calculations.

    python scripts/build_kshell.py                         # verify, patch, build, check
    python scripts/build_kshell.py --archives <dir>        # take the archives from a local folder instead
    python scripts/build_kshell.py --check-only            # run the 28Si check on an existing tools/kshell/
    python scripts/build_kshell.py --compare <dir> <dir>   # compare two folders of KSHELL logs (Linux vs Windows)

Steps: take ``kshell.tar.gz`` from ``third_party/kshell/`` (vendored: Mendeley Data refuses automated downloads) and
on Windows download the OpenBLAS 0.3.34 binaries from GitHub, and verify both by checksum; unpack into
``build/kshell/``; apply physim's rename patch (below); compile ``kshell.exe``, ``transit.exe`` and
``count_dim.exe`` with gfortran (OpenMP, no MPI); and assemble ``tools/kshell/`` with the executables, every DLL
they need (Windows), the interaction files (``snt/``), the source archive with the patch and this script, the
licence texts and a ``VERSION`` file. Then the 28Si/USDB case from the test fixtures is run with the new
executables in a clean folder, with nothing but Windows itself on the PATH (so a missing DLL shows up), and its
energies and B(E2)/B(M1) are compared with the fixture logs: energies to 1 keV, reduced probabilities to 0.1 %.

KSHELL is GPLv3 and stays a separate program next to physim (MIT): it is run as an executable, never linked.

The rename patch. KSHELL's Lanczos restart with vectors on disk (``mode_lv_hdd = 1``, its default) renames temporary
files onto names that exist. POSIX ``rename`` replaces the target; the MinGW one fails silently on Windows, KSHELL
reads stale vectors back and stops with "possible JJ collapse". The patch adds ``rename_replace(old, new)`` to
``lanczos.f90`` and ``bp_io.F90``, which behaves as POSIX does (replaces ``new``; does nothing when ``old`` is
missing or equals ``new``, both of which happen in KSHELL), and calls it instead of ``rename``. It is applied on
every platform: on Linux it changes nothing.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import difflib
import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional

REPO = Path(__file__).resolve().parents[1]
WINDOWS = sys.platform == "win32"


@dataclass(frozen=True)
class Archive:
    name: str
    url: str
    sha256: str
    md5: Optional[str] = None


#: The KSHELL release physim bundles (GPLv3): N. Shimizu, T. Mizusaki, Y. Utsuno and Y. Tsunoda,
#: Comput. Phys. Commun. 244, 372 (2019), Mendeley Data kgzdz6ryyk.
KSHELL = Archive(
    "kshell.tar.gz",
    "https://data.mendeley.com/public-files/datasets/kgzdz6ryyk/files/7ad64d9c-5df2-4027-a7f2-768befe2894a/"
    "file_downloaded",
    "ce0659c3b39dc8c8ed623f6ebe1e10766445fc72372b6e33e952a578c7227526",
)
#: The copy of that release kept in the repository, used unless --archives names another folder. The URL above is
#: only a last resort: Mendeley Data refuses non-browser clients (403, or a web page instead of the archive).
VENDORED_KSHELL = REPO / "third_party" / "kshell" / KSHELL.name
#: OpenBLAS (BSD-3-Clause) for LAPACK/BLAS on Windows. The MD5 is the one published with the release.
OPENBLAS_VERSION = "0.3.34"
OPENBLAS = Archive(
    f"OpenBLAS-{OPENBLAS_VERSION}-x64.zip",
    f"https://github.com/OpenMathLib/OpenBLAS/releases/download/v{OPENBLAS_VERSION}/"
    f"OpenBLAS-{OPENBLAS_VERSION}-x64.zip",
    "e9cb6134541f36c27346d5fc5995652f060fba227cebbbabcbda5a5a44d7c76b",
    md5="5047115a08a7d39471375a25994108c5",
)

#: Raise when the patch changes, so builds (and caches keyed on VERSION) can tell them apart.
PATCH_LEVEL = 1

#: KSHELL's modules in dependency order (from its Makefile), then the three programs.
MODULES = [
    "constant.f90", "model_space.f90", "class_stopwatch.F90", "rotation_group.f90", "harmonic_oscillator.f90",
    "sp_matrix_element.f90", "operator_jscheme.f90", "operator_mscheme.f90", "interaction.f90", "lib_matrix.F90",
    "partition.F90", "wavefunction.F90", "bridge_partitions.F90", "bp_io.F90", "bp_expc_val.F90", "lanczos.f90",
    "bp_block.F90", "block_lanczos.F90",
]
PROGRAMS = {"kshell.exe": "kshell.F90", "transit.exe": "transit.F90", "count_dim.exe": "count_dim.f90"}

FFLAGS = ["-O3", "-fopenmp"]
WINDOWS_LDFLAGS = ["-Wl,--stack,40000000"]  # KSHELL's own advice for Windows (Cygwin) with OpenMP
LINUX_LIBS = ["-llapack", "-lblas"]

FIXTURES = REPO / "tests" / "reference" / "nuclear" / "shellmodel" / "kshell" / "si28_usdb"
E_TOL_MEV = 0.001  # 1 keV
B_REL_TOL = 0.001  # 0.1 %
B_ABS_TOL = 2e-4  # the logs print four decimals


def say(*args) -> None:
    print(*args, flush=True)


# ------------------------------------------------------------------------------------------------ archives


def file_digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify(path: Path, archive: Archive) -> None:
    """Raise ValueError unless ``path`` has the archive's SHA-256 (and MD5, when one is pinned)."""
    got = file_digest(path, "sha256")
    if got != archive.sha256:
        raise ValueError(f"{path}: SHA-256 {got}, expected {archive.sha256}")
    if archive.md5:
        got = file_digest(path, "md5")
        if got != archive.md5:
            raise ValueError(f"{path}: MD5 {got}, expected {archive.md5}")


def find_archive(archive: Archive, folders: Iterable[Path]) -> Optional[Path]:
    """The first file named like the archive under the folders (searched recursively) whose checksums match."""
    for folder in folders:
        for path in sorted(Path(folder).rglob(archive.name)):
            try:
                verify(path, archive)
                return path
            except ValueError as e:
                say("  skipping", e)
    return None


#: The first bytes of the archive types fetched here, to tell an archive from an error page.
MAGIC = {".gz": b"\x1f\x8b", ".zip": b"PK\x03\x04"}


def download(archive: Archive, folder: Path) -> Path:
    """Fetch the archive and verify it. Fails with a clear message, and leaves nothing behind, when the server
    refuses or sends something else (a web page instead of the archive)."""
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / archive.name
    say("downloading", archive.url)
    req = urllib.request.Request(archive.url, headers={"User-Agent": "physim-build"})
    try:
        with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as f:
            shutil.copyfileobj(r, f)
    except OSError as e:  # urllib.error.URLError and HTTPError are OSErrors
        dest.unlink(missing_ok=True)
        raise SystemExit(f"downloading {archive.name} failed: {e}. Pass --archives <folder with a copy>.")
    head = dest.read_bytes()[:4]
    magic = MAGIC.get(Path(archive.name).suffix)
    try:
        if magic and not head.startswith(magic):
            raise ValueError(f"{archive.url} sent something that is not a {Path(archive.name).suffix} archive "
                             f"(starts with {head!r}; probably a web page)")
        verify(dest, archive)
    except ValueError as e:
        dest.unlink(missing_ok=True)
        raise SystemExit(f"{e}. Pass --archives <folder with a copy>.")
    return dest


def get_archive(archive: Archive, local: Optional[Path], download_dir: Path, vendored: Optional[Path] = None) -> Path:
    """The archive, checked: from --archives when it has a copy, else the vendored copy, else earlier downloads,
    else the network (the archives' folder, when given, must then have it)."""
    path = find_archive(archive, [local]) if local and local.is_dir() else None
    if path:
        pass
    elif vendored and vendored.is_file():
        verify(vendored, archive)  # a changed vendored copy is an error, never silently replaced by a download
        path = vendored
    elif local:
        raise SystemExit(f"{archive.name} with SHA-256 {archive.sha256} not found under {local}")
    else:
        if vendored:
            say(f"WARNING: {vendored} is missing; trying {archive.url}")
        path = find_archive(archive, [download_dir] if download_dir.is_dir() else []) or download(archive,
                                                                                                 download_dir)
    say(f"{archive.name}: {path} (checksum ok)")
    return path


def unpack_tar(path: Path, dest: Path) -> None:
    with tarfile.open(path) as tar:
        for m in tar.getmembers():
            if m.name.startswith(("/", "\\")) or ".." in Path(m.name).parts or not (m.isfile() or m.isdir()):
                raise ValueError(f"{path}: unexpected member {m.name!r}")
        try:
            tar.extractall(dest, filter="data")
        except TypeError:  # Python without extraction filters; the members were checked above
            tar.extractall(dest)


def unpack_zip(path: Path, dest: Path) -> None:
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            if name.startswith(("/", "\\")) or ".." in Path(name).parts or ":" in name:
                raise ValueError(f"{path}: unexpected member {name!r}")
        z.extractall(dest)


# ------------------------------------------------------------------------------------------------ the patch

RENAME_HELPER = """
  subroutine rename_replace(old, new)
    ! physim patch: rename old to new, replacing new if it exists (POSIX behaviour; the Windows rename fails
    ! otherwise). Like POSIX, a missing old or old == new changes nothing (KSHELL relies on both).
    character(len=*), intent(in) :: old, new
    logical :: ex
    integer :: lun
    if (trim(old) == trim(new)) return
    inquire(file=old, exist=ex)
    if (.not. ex) return
    inquire(file=new, exist=ex)
    if (ex) then
       open(newunit=lun, file=new, status='old')
       close(lun, status='delete')
    end if
    call rename(old, new)
  end subroutine rename_replace

"""

_RENAME_CALL = re.compile(r"\bcall\s+rename\s*\(", re.IGNORECASE)

#: The files patched and the module each helper goes into.
PATCHED = {"lanczos.f90": "lanczos", "bp_io.F90": "bp_io"}


def patch_rename_text(text: str, module: str) -> tuple[str, int]:
    """Replace every ``call rename(`` with ``call rename_replace(`` and add the helper before ``end module``.

    Returns the new text and the number of calls replaced. Text already patched comes back unchanged with 0. Raises
    ValueError when the module end is missing or there is no rename call to patch (the source is not the one this
    patch was written for)."""
    if "rename_replace" in text:
        return text, 0
    ends = list(re.finditer(rf"^[ \t]*end[ \t]+module[ \t]+{re.escape(module)}\b", text, re.IGNORECASE | re.MULTILINE))
    if not ends:
        raise ValueError(f"'end module {module}' not found")
    new, n = _RENAME_CALL.subn("call rename_replace(", text)
    if n == 0:
        raise ValueError(f"no 'call rename(' in module {module}")
    ends = list(re.finditer(rf"^[ \t]*end[ \t]+module[ \t]+{re.escape(module)}\b", new, re.IGNORECASE | re.MULTILINE))
    i = ends[-1].start()
    return new[:i] + RENAME_HELPER + new[i:], n


def patch_sources(src: Path) -> tuple[dict[str, int], str]:
    """Patch the files in place; return the calls replaced per file and a unified diff of the change."""
    counts, diff = {}, []
    for name, module in PATCHED.items():
        path = src / name
        old = path.read_text(encoding="latin-1")
        new, n = patch_rename_text(old, module)
        with open(path, "w", encoding="latin-1", newline="") as f:
            f.write(new)
        counts[name] = n
        diff += difflib.unified_diff(old.splitlines(keepends=True), new.splitlines(keepends=True),
                                     f"a/kshell/src/{name}", f"b/kshell/src/{name}")
    return counts, "".join(diff)


# ------------------------------------------------------------------------------------------------ building


def compiler_info(fc: str) -> tuple[str, list[int]]:
    banner = subprocess.run([fc, "--version"], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    version = subprocess.run([fc, "-dumpfullversion", "-dumpversion"], capture_output=True, text=True,
                             check=True).stdout.strip()
    return banner.strip(), [int(x) for x in re.findall(r"\d+", version)[:3]]


def compile_kshell(src: Path, obj: Path, fc: str, fflags: list[str], link: list[str], jobs: int) -> dict[str, Path]:
    """Compile the modules in order, then link the three programs. Returns the executables."""
    obj.mkdir(parents=True, exist_ok=True)
    base = [fc, *fflags, "-I", src, "-J", obj, "-I", obj]
    t0 = time.time()
    for name in MODULES:
        say("  compiling", name)
        r = subprocess.run([str(c) for c in [*base, "-c", src / name, "-o", obj / (Path(name).stem + ".o")]],
                           capture_output=True, text=True)
        if r.returncode:  # gfortran warns a lot about KSHELL's BLAS calls; show its output only on failure
            say(r.stdout + r.stderr)
            raise SystemExit(f"compiling {name} failed")
    objects = [obj / (Path(n).stem + ".o") for n in MODULES]

    def build(exe: str, main: str) -> subprocess.Popen:
        return subprocess.Popen([str(c) for c in [*base, "-o", obj / exe, src / main, *objects, *link]])

    say("  compiling and linking", ", ".join(PROGRAMS))
    procs = [build(exe, main) for exe, main in PROGRAMS.items()] if jobs > 1 else []
    if not procs:
        for exe, main in PROGRAMS.items():
            if build(exe, main).wait():
                raise SystemExit(f"building {exe} failed")
    for (exe, _), p in zip(PROGRAMS.items(), procs):
        if p.wait():
            raise SystemExit(f"building {exe} failed")
    say(f"  built in {time.time() - t0:.0f} s")
    return {exe: obj / exe for exe in PROGRAMS}


# ------------------------------------------------------------------------------------------------ DLLs

_DLL_NAME = re.compile(r"^\s*DLL Name:\s*(\S+)\s*$", re.MULTILINE)


def parse_objdump_imports(text: str) -> list[str]:
    """The DLL names in ``objdump -p`` output, in order, without repeats."""
    seen: dict[str, None] = {}
    for name in _DLL_NAME.findall(text):
        seen.setdefault(name, None)
    return list(seen)


def is_api_set(name: str) -> bool:
    """Windows API sets (api-ms-win-*, ext-ms-*): virtual names the loader resolves, never files to ship."""
    return name.lower().startswith(("api-ms-win-", "ext-ms-"))


def collect_dlls(binaries: Iterable[Path], imports_of: Callable[[Path], list[str]], search_dirs: Iterable[Path],
                 is_system: Callable[[str], bool]) -> dict[str, Path]:
    """Every DLL the binaries need, recursively, that is not part of Windows: name (lower case) -> file.

    ``imports_of`` gives a binary's imported DLL names; ``is_system`` says whether a name belongs to Windows. A DLL
    found neither in ``search_dirs`` nor in Windows raises FileNotFoundError naming who wanted it."""
    dirs = [Path(d) for d in search_dirs]
    found: dict[str, Path] = {}
    queue = [(Path(b), None) for b in binaries]
    while queue:
        binary, _ = queue.pop(0)
        for name in imports_of(binary):
            key = name.lower()
            if key in found or is_api_set(name) or is_system(name):
                continue
            path = next((d / name for d in dirs if (d / name).is_file()), None)
            if path is None:
                raise FileNotFoundError(f"{name} (needed by {binary.name}) is not in {', '.join(map(str, dirs))}")
            found[key] = path
            queue.append((path, binary))
    return found


def windows_system_dll(name: str) -> bool:
    system32 = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32"
    return (system32 / name).is_file()


def objdump_imports(objdump: str) -> Callable[[Path], list[str]]:
    def imports(path: Path) -> list[str]:
        out = subprocess.run([objdump, "-p", str(path)], capture_output=True, text=True, check=True).stdout
        return parse_objdump_imports(out)
    return imports


# ------------------------------------------------------------------------------------------------ licences


def is_gpl3(text: str) -> bool:
    return "GNU GENERAL PUBLIC LICENSE" in text and "Version 3, 29 June 2007" in text


#: The verbatim GPL v3 text kept in the repository (KSHELL's tarball has none), copied to tools/kshell/COPYING.
GPL_TEXT = REPO / "installer" / "windows" / "licenses" / "GPL-3.0.txt"


def find_gpl_text(explicit: Optional[Path], fc_dir: Optional[Path]) -> Optional[Path]:
    """The GPL v3 text: the given file, else the repository's copy; only if that has gone missing, the system's
    common licences (Debian, Ubuntu) or GCC's COPYING3 in the compiler's tree."""
    if explicit:
        if explicit.is_file() and is_gpl3(explicit.read_text(encoding="utf-8", errors="replace")):
            return explicit
        raise SystemExit(f"{explicit} is not the GPL v3 text")
    candidates = [GPL_TEXT, Path("/usr/share/common-licenses/GPL-3")]
    if fc_dir:
        root = fc_dir.parent
        for sub in ("share/licenses", "share/doc", "licenses"):
            if (root / sub).is_dir():
                candidates += sorted((root / sub).rglob("COPYING3"))
    for path in candidates:
        try:
            if path.is_file() and is_gpl3(path.read_text(encoding="utf-8", errors="replace")):
                return path
        except OSError:
            continue
    return None


def leading_comment(text: str) -> str:
    """The first C block comment of a file, without the comment markers and the leading ' * '."""
    m = re.search(r"/\*(.*?)\*/", text, re.DOTALL)
    if not m:
        return ""
    lines = [re.sub(r"^\s*\*+ ?", "", ln).rstrip() for ln in m.group(1).splitlines()]
    lines = [ln for ln in lines if ln.strip("*").strip() or ln == ""]
    return "\n".join(lines).strip() + "\n"


# ------------------------------------------------------------------------------------------------ assembling


def readme_text(windows: bool, dlls: Iterable[str]) -> str:
    dll_list = ", ".join(sorted(dlls)) or "none"
    return """KSHELL for physim
=================

These programs are KSHELL, a shell-model code by N. Shimizu, T. Mizusaki, Y. Utsuno and Y. Tsunoda,
Comput. Phys. Commun. 244, 372 (2019), from the release published with that paper (Mendeley Data kgzdz6ryyk).
Please cite the paper when you use results from KSHELL.

KSHELL is free software under the GNU General Public License, version 3 or later (COPYING). It is a separate
program: physim (MIT licence) only writes its input files, runs the executables and reads their output.
The complete source is in source/: the original release (kshell.tar.gz), physim's patch to it
(physim-rename.patch, explained in PATCH.txt) and the script that built these files (build_kshell.py).

Contents
  kshell.exe, transit.exe, count_dim.exe   built from the source above (see VERSION)
  snt/                                     the interaction files of the release (GPLv3, as above)
  source/                                  the corresponding source
  VERSION                                  release, patch, compiler, libraries and build date
""" + (f"""  *.dll                                    run-time libraries: {dll_list}

Run-time libraries (Windows)
  libopenblas.dll    OpenBLAS {OPENBLAS_VERSION}, BSD 3-clause licence (LICENSE-OpenBLAS.txt), https://www.openblas.net
  libgfortran, libgomp, libquadmath, libgcc_s
                     GCC run-time libraries, GPL v3 with the GCC Runtime Library Exception
  libwinpthread      MinGW-w64 winpthreads, MIT and BSD licences
  libdl              dlfcn-win32, MIT licence (present only with some MinGW builds)
""" if windows else """
On Linux the programs use the system's LAPACK and BLAS libraries.
""")


PATCH_TEXT = f"""physim's patch to KSHELL (patch level {PATCH_LEVEL})
==========================================

File: source/physim-rename.patch (a unified diff against kshell.tar.gz; apply with "patch -p1" in the folder
that contains kshell/).

Why. KSHELL's Lanczos restart with vectors on disk (mode_lv_hdd = 1, its default) renames temporary files onto
names that exist. POSIX rename replaces the target; on Windows the MinGW rename fails silently, so KSHELL read
stale vectors back and stopped with "possible JJ collapse".

What. A subroutine rename_replace(old, new) in lanczos.f90 and bp_io.F90 that behaves as POSIX rename does:
it replaces new when new exists, and does nothing when old is missing or equals new (KSHELL relies on both).
Every "call rename(" in those two files calls it instead. Nothing else is changed. With the patch, disk and
memory modes give identical eigenvalues; on Linux the results are unchanged.
"""


def assemble(out: Path, exes: dict[str, Path], dlls: dict[str, Path], src_root: Path, kshell_archive: Path,
             diff: str, version: str, gpl: Optional[Path], openblas_licence: Optional[str]) -> None:
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    for name, path in exes.items():
        shutil.copy2(path, out / name)
    for path in dlls.values():
        shutil.copy2(path, out / path.name)
    (out / "snt").mkdir()
    for snt in sorted((src_root / "snt").glob("*.snt")):
        shutil.copy2(snt, out / "snt" / snt.name)
    (out / "source").mkdir()
    shutil.copy2(kshell_archive, out / "source" / KSHELL.name)
    (out / "source" / "physim-rename.patch").write_bytes(diff.encode("latin-1"))
    shutil.copy2(Path(__file__), out / "source" / "build_kshell.py")
    (out / "PATCH.txt").write_text(PATCH_TEXT, encoding="utf-8")
    (out / "VERSION").write_text(version, encoding="utf-8")
    (out / "README.txt").write_text(readme_text(bool(dlls), [p.name for p in dlls.values()]), encoding="utf-8")
    if gpl:
        shutil.copy2(gpl, out / "COPYING")
    else:
        (out / "COPYING.TODO").write_text(
            "TODO: the GNU GPL v3 text belongs here as COPYING. KSHELL's release has none and the repository's copy\n"
            "(installer/windows/licenses/GPL-3.0.txt) was missing. Restore it or rebuild with --gpl-text <file>.\n"
            "installer/windows/build.py does not bundle KSHELL without it.\n",
            encoding="utf-8")
    if openblas_licence:
        (out / "LICENSE-OpenBLAS.txt").write_text(openblas_licence, encoding="utf-8")


def version_text(banner: str, flags: list[str], link: list[str], blas: str, counts: dict[str, int]) -> str:
    now = _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0)
    patched = ", ".join(f"{name}: {n} call(s)" for name, n in counts.items())
    return (f"KSHELL 2019 CPC release + physim rename patch {PATCH_LEVEL}\n"
            f"release: {KSHELL.name} (Comput. Phys. Commun. 244, 372 (2019); Mendeley Data kgzdz6ryyk), "
            f"sha256 {KSHELL.sha256}\n"
            f"patch: physim rename patch {PATCH_LEVEL} ({patched}); source/physim-rename.patch\n"
            f"compiler: {banner}\n"
            f"flags: {' '.join(flags)}; link: {' '.join(Path(x).name if os.sep in x else x for x in link)}\n"
            f"blas: {blas}\n"
            f"platform: {platform.system()} {platform.release()} {platform.machine()}\n"
            f"built: {now.isoformat().replace('+00:00', 'Z')}\n")


def build(args) -> Path:
    fc = shutil.which(args.fc) or args.fc
    if not Path(fc).exists() and not shutil.which(fc):
        raise SystemExit(f"gfortran not found ({args.fc}); pass --fc")
    fc_dir = Path(fc).resolve().parent
    banner, ver = compiler_info(fc)
    say("compiler:", banner)
    flags = list(FFLAGS) + (["-fallow-argument-mismatch"] if ver and ver[0] >= 10 else [])

    work = args.build_dir.resolve()
    downloads = work / "downloads"
    kshell_tar = get_archive(KSHELL, args.archives, downloads, VENDORED_KSHELL)
    for sub in ("src", "openblas", "obj"):
        shutil.rmtree(work / sub, ignore_errors=True)
    unpack_tar(kshell_tar, work / "src")
    root = work / "src" / "kshell"
    counts, diff = patch_sources(root / "src")
    say("patched:", ", ".join(f"{k} ({v} rename calls)" for k, v in counts.items()))

    openblas_licence = None
    if WINDOWS:
        ob_zip = get_archive(OPENBLAS, args.archives, downloads)
        unpack_zip(ob_zip, work / "openblas")
        ob = work / "openblas"
        # The DLL's import library: the static libopenblas.a does not link with UCRT gfortran (__imp__cprintf).
        link = [str(ob / "lib" / "libopenblas.dll.a"), *WINDOWS_LDFLAGS]
        blas = f"OpenBLAS {OPENBLAS_VERSION} ({OPENBLAS.name}, sha256 {OPENBLAS.sha256})"
        openblas_licence = leading_comment((ob / "include" / "cblas.h").read_text(encoding="utf-8", errors="replace"))
        if "Redistribution and use" not in openblas_licence:
            openblas_licence = None
    else:
        link = args.libs.split() if args.libs else list(LINUX_LIBS)
        blas = "system: " + " ".join(link)

    exes = compile_kshell(root / "src", work / "obj", fc, flags, link, args.jobs)

    dlls: dict[str, Path] = {}
    if WINDOWS:
        objdump = shutil.which("objdump", path=str(fc_dir)) or shutil.which("objdump")
        if not objdump:
            raise SystemExit("objdump not found next to gfortran or on the PATH")
        dlls = collect_dlls(exes.values(), objdump_imports(objdump), [work / "openblas" / "bin", fc_dir],
                            windows_system_dll)
        say("DLLs:", ", ".join(p.name for p in dlls.values()))

    gpl = find_gpl_text(args.gpl_text, fc_dir)
    version = version_text(banner, flags, link, blas, counts)
    out = args.out.resolve()
    assemble(out, exes, dlls, root, kshell_tar, diff, version, gpl, openblas_licence)
    say(f"wrote {out}")
    if not gpl:
        say(f"WARNING: {GPL_TEXT} missing: tools/kshell/COPYING.TODO written; rebuild with --gpl-text <file>")
    if WINDOWS and not openblas_licence:
        say("WARNING: no OpenBLAS licence text found in the archive")
    return out


# ------------------------------------------------------------------------------------------------ the check

_STATE = re.compile(r"^\s*(\d+)\s+<H>:\s+(-?\d+\.\d+)\s+<JJ>:\s+(-?\d+\.\d+)\s+J:\s+(\d+)/2\s+prty\s+(-?\d+)",
                    re.MULTILINE)
_SECTION = re.compile(r"^\s*(E\d|M\d)\s+transition", re.IGNORECASE)
_ROW = re.compile(r"^\s*(\d+)\(\s*(\d+)\)\s+(-?\d+\.\d+)\s+(\d+)\(\s*(\d+)\)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
                  r"\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")


def parse_energies(text: str) -> dict[int, tuple[float, int, int]]:
    """kshell.exe log: state index -> (energy MeV, 2J, parity); the last listing wins."""
    return {int(m[1]): (float(m[2]), int(m[4]), int(m[5])) for m in _STATE.finditer(text)}


def parse_transitions(text: str) -> dict[tuple, tuple[float, float, float]]:
    """transit.exe log: (kind, 2Jf, nf, 2Ji, ni) -> (B->, B<-, Mom.). Mred is left out: its sign follows the
    arbitrary phases of the wave functions, which may differ between builds."""
    rows, kind = {}, None
    for line in text.splitlines():
        m = _SECTION.match(line)
        if m:
            kind = m[1].upper()
            continue
        m = _ROW.match(line)
        if m and kind:
            rows[(kind, int(m[1]), int(m[2]), int(m[4]), int(m[5]))] = (float(m[9]), float(m[10]), float(m[11]))
    return rows


def transition_sections(text: str) -> list[str]:
    """The multipolarities a transit.exe log lists (a 0+ to 0+ log has E2 and M1 sections with no rows)."""
    return [m[1].upper() for m in map(_SECTION.match, text.splitlines()) if m]


def b_close(a: float, b: float) -> bool:
    return abs(a - b) <= max(B_REL_TOL * abs(b), B_ABS_TOL)


def compare_logs(new: dict[str, str], ref: dict[str, str]) -> list[str]:
    """Problems found comparing logs of the same name (kshell logs by energy, transit logs by B and moments)."""
    problems = []
    for name, ref_text in sorted(ref.items()):
        if name not in new:
            problems.append(f"{name}: missing")
            continue
        text = new[name]
        if "collapse" in text.lower():
            problems.append(f"{name}: KSHELL reports a collapse")
        if "<H>:" in ref_text:
            e_new, e_ref = parse_energies(text), parse_energies(ref_text)
            if not e_ref:
                problems.append(f"{name}: no states in the reference")
            for i, (e, j2, p) in e_ref.items():
                if i not in e_new:
                    problems.append(f"{name}: state {i} missing")
                elif e_new[i][1:] != (j2, p) or abs(e_new[i][0] - e) > E_TOL_MEV:
                    problems.append(f"{name}: state {i} {e_new[i]} vs {(e, j2, p)}")
        else:
            t_new, t_ref = parse_transitions(text), parse_transitions(ref_text)
            if transition_sections(text) != transition_sections(ref_text):
                problems.append(f"{name}: sections {transition_sections(text)} vs {transition_sections(ref_text)}")
            for key, vals in t_ref.items():
                if key not in t_new:
                    problems.append(f"{name}: {key} missing")
                elif not all(b_close(a, b) for a, b in zip(t_new[key], vals)):
                    problems.append(f"{name}: {key} {t_new[key]} vs {vals}")
    return problems


def read_logs(folder: Path) -> dict[str, str]:
    return {p.name: p.read_text(encoding="latin-1") for p in sorted(folder.glob("log_*.txt"))}


def clean_env(threads: Optional[int]) -> dict[str, str]:
    """The environment for the check: on Windows nothing but Windows on the PATH, so the DLLs next to the
    executables must be complete."""
    env = dict(os.environ)
    if WINDOWS:
        root = os.environ.get("SystemRoot", r"C:\Windows")
        keep = {"SYSTEMROOT", "SYSTEMDRIVE", "WINDIR", "TEMP", "TMP", "USERPROFILE", "COMSPEC",
                "NUMBER_OF_PROCESSORS", "PROCESSOR_ARCHITECTURE", "PATHEXT", "LOCALAPPDATA", "APPDATA"}
        env = {k: v for k, v in os.environ.items() if k.upper() in keep}
        env["PATH"] = os.pathsep.join([rf"{root}\System32", root, rf"{root}\System32\Wbem"])
    env.update({"OMP_STACKSIZE": "64m", "GFORTRAN_UNBUFFERED_PRECONNECTED": "y"})
    if threads:
        env["OMP_NUM_THREADS"] = str(threads)
    return env


def run_program(exe: Path, input_name: str, log: Path, cwd: Path, env: dict) -> int:
    with open(log, "wb") as f:
        return subprocess.run([str(exe), input_name], cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT).returncode


def check(out: Path, fixtures: Path, work: Optional[Path], threads: Optional[int]) -> bool:
    """Run the 28Si/USDB case (vectors on disk, as in the fixtures, then the 0+ again in memory) and compare."""
    out = out.resolve()
    for exe in ("kshell.exe", "transit.exe"):
        if not (out / exe).is_file():
            raise SystemExit(f"{out / exe} not found: build first")
    work = Path(tempfile.mkdtemp(prefix="kshell-check-")) if work is None else work.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    inputs = fixtures / "inputs"
    for p in inputs.iterdir():
        shutil.copy2(p, work / p.name)
    shutil.copy2(out / "snt" / "usdb.snt", work / "usdb.snt")
    env = clean_env(threads)
    say(f"28Si/USDB check in {work}" + (" (PATH: Windows only)" if WINDOWS else ""))

    ok = True
    t0 = time.time()
    runs = [("kshell.exe", f"Si28_usdb_{j}.input", f"log_Si28_usdb_j{j}p.txt") for j in (0, 4, 8)]
    runs += [("transit.exe", f"Si28_usdb_{a}_{b}.input", f"log_Si28_usdb_tr_j{a}p_j{b}p.txt")
             for a, b in ((0, 0), (0, 4), (4, 4), (4, 8), (8, 8))]
    # The 0+ again with the Lanczos vectors in memory: both storage modes must agree.
    mem = (inputs / "Si28_usdb_0.input").read_text()
    mem = re.sub(r"mode_lv_hdd\s*=\s*\d+", "mode_lv_hdd = 0", mem).replace("Si28_usdb_j0p.wav", "mem_j0p.wav")
    (work / "Si28_usdb_0_mem.input").write_text(mem)
    runs.append(("kshell.exe", "Si28_usdb_0_mem.input", "mem_log_Si28_usdb_j0p.txt"))
    for exe, inp, log in runs:
        code = run_program(out / exe, inp, work / log, work, env)
        say(f"  {exe} {inp}: exit {code}" + (" (a DLL is missing)" if code == 0xC0000135 else ""))
        if code:
            ok = False
            tail = (work / log).read_text(encoding="latin-1", errors="replace").splitlines()[-15:]
            say("    " + "\n    ".join(tail))
    say(f"  ran in {time.time() - t0:.0f} s")

    ref = read_logs(fixtures / "logs")
    new = read_logs(work)
    problems = compare_logs(new, ref)
    mem_log = {"log_Si28_usdb_j0p.txt": (work / "mem_log_Si28_usdb_j0p.txt").read_text(encoding="latin-1")}
    problems += ["memory mode: " + p for p in compare_logs(mem_log, {"log_Si28_usdb_j0p.txt":
                                                                         ref["log_Si28_usdb_j0p.txt"]})]
    n_states = sum(len(parse_energies(t)) for n, t in ref.items() if "<H>:" in t)
    n_trans = sum(len(parse_transitions(t)) for n, t in ref.items() if "<H>:" not in t)
    for p in problems:
        say("  FAIL", p)
    ok = ok and not problems
    say(f"  {'PASS' if ok else 'FAIL'}: {n_states} energies (1 keV) and {n_trans} transition rows "
        f"(B to 0.1 %) against {fixtures / 'logs'}; disk and memory modes")
    return ok


# ------------------------------------------------------------------------------------------------ main


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archives", type=Path,
                    help="folder with local copies of the archives (searched recursively); a KSHELL archive found "
                         "there is used instead of the vendored one, and OpenBLAS is not downloaded")
    ap.add_argument("--out", type=Path, default=REPO / "tools" / "kshell")
    ap.add_argument("--build-dir", type=Path, default=REPO / "build" / "kshell")
    ap.add_argument("--fc", default="gfortran", help="the gfortran to use")
    ap.add_argument("--libs", help=f"LAPACK/BLAS link flags on Linux (default {' '.join(LINUX_LIBS)})")
    ap.add_argument("--gpl-text", type=Path,
                    help=f"a file with the GPL v3 text, copied as COPYING (default {GPL_TEXT.relative_to(REPO)})")
    ap.add_argument("--jobs", type=int, default=3, help="programs linked at once")
    ap.add_argument("--threads", type=int, help="OMP_NUM_THREADS for the check (default: all)")
    ap.add_argument("--work", type=Path, help="folder for the check's runs (default: a new temporary folder)")
    ap.add_argument("--fixtures", type=Path, default=FIXTURES)
    ap.add_argument("--no-check", action="store_true", help="build only")
    ap.add_argument("--check-only", action="store_true", help="run the 28Si check on --out without building")
    ap.add_argument("--compare", nargs=2, type=Path, metavar=("NEW", "REF"),
                    help="compare the log_*.txt files of two folders and exit")
    args = ap.parse_args(argv)

    if args.compare:
        problems = compare_logs(read_logs(args.compare[0]), read_logs(args.compare[1]))
        for p in problems:
            say("FAIL", p)
        say("PASS" if not problems else f"{len(problems)} difference(s)")
        return 1 if problems else 0
    out = args.out if args.check_only else build(args)
    if args.no_check:
        return 0
    return 0 if check(out, args.fixtures, args.work, args.threads) else 1


if __name__ == "__main__":
    raise SystemExit(main())
