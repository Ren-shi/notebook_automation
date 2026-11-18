"""Assemble the Windows planner bundle: the official embeddable Python, the physim wheel with the app's packages,
the icon, and KSHELL for shell-model calculations when tools/kshell/ has been built (scripts/build_kshell.py; a
separate GPLv3 program, copied to tools/kshell/ in the bundle; without it the bundle is made with a warning).
installer/windows/planner.iss then packs the folder into an installer.

    python installer/windows/build.py --wheel dist/physim_engine-0.1.0-cp39-abi3-win_amd64.whl --out build/planner

Run it with a CPython of the same minor version as PYTHON_VERSION (it installs the packages with the host's pip and
precompiles them). See backlog/44-one-click-installer.md for why the bundle is built this way.
"""

from __future__ import annotations

import argparse
import compileall
import hashlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

#: The embeddable distribution from python.org. Its python.exe and pythonw.exe are signed by the Python Software
#: Foundation, so the bundle contains no unsigned executables of ours.
PYTHON_VERSION = "3.12.10"
EMBED_URL = "https://www.python.org/ftp/python/{v}/python-{v}-embed-amd64.zip"

#: Large parts of the bundled packages that the planner never loads.
PRUNE = ["numpy/_core/tests", "numpy/tests", "numpy/*/tests", "matplotlib/tests", "matplotlib/mpl-data/sample_data",
         "pandas/tests", "plotly/tests", "pip", "setuptools"]

#: Pure-Python dependencies published only as source. The bundle is installed with --only-binary (it targets another
#: platform), so these are built into wheels first; without proxy_tools, pip falls back to pywebview 3.4, which
#: NiceGUI's native mode cannot import.
SOURCE_ONLY = ["proxy_tools"]


def fetch_python(dest: Path) -> None:
    url = EMBED_URL.format(v=PYTHON_VERSION)
    print("downloading", url)
    with urllib.request.urlopen(url, timeout=120) as r:
        data = r.read()
    print(f"  {len(data) / 1e6:.1f} MB, sha256 {hashlib.sha256(data).hexdigest()}")
    zipfile.ZipFile(io.BytesIO(data)).extractall(dest)


def enable_site_packages(python: Path) -> None:
    """The embeddable Python ignores site-packages until its ._pth file says otherwise."""
    pth = next(python.glob("python3*._pth"))
    lines = [ln for ln in pth.read_text().splitlines() if ln.strip() and ln.strip() != "#import site"]
    pth.write_text("\n".join(lines + ["Lib\\site-packages", "import site", ""]))
    (python / "Lib" / "site-packages").mkdir(parents=True, exist_ok=True)


def install(wheel: Path, site: Path) -> None:
    major, minor = PYTHON_VERSION.split(".")[:2]
    with tempfile.TemporaryDirectory() as local:  # outside the bundle, which planner.iss packs whole
        subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--disable-pip-version-check",
                        "--wheel-dir", local, *SOURCE_ONLY], check=True)
        cmd = [sys.executable, "-m", "pip", "install", "--no-cache-dir", "--disable-pip-version-check",
               "--target", str(site), "--platform", "win_amd64", "--python-version", f"{major}.{minor}",
               "--implementation", "cp", "--only-binary=:all:", "--find-links", local,
               f"{wheel.resolve()}[app,root,window]"]
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)


def prune(site: Path) -> None:
    for pattern in PRUNE:
        for p in site.glob(pattern):
            shutil.rmtree(p, ignore_errors=True)
    for p in list(site.rglob("__pycache__")):
        shutil.rmtree(p, ignore_errors=True)
    for p in site.glob("bin"):  # console scripts written for the host Python
        shutil.rmtree(p, ignore_errors=True)


def make_icon(path: Path) -> bool:
    """A hyperbolic orbit around a nucleus, at the sizes Windows asks for. Needs Pillow; without it the shortcut
    shows Python's icon."""
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        print("Pillow not installed: no icon")
        return False
    import math

    n, k = 256, 4  # drawn at 4x and scaled down, for smooth edges
    img = Image.new("RGBA", (n * k, n * k), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    d.rounded_rectangle([8 * k, 8 * k, (n - 8) * k, (n - 8) * k], radius=48 * k, fill=(30, 41, 59, 255))
    # A repulsive Coulomb orbit: the hyperbola branch x = a(cosh w + ε), y = a√(ε² − 1) sinh w from the nucleus.
    eps, a, cx, cy = 1.6, 22.0, 168.0, 128.0
    for i in range(-800, 801):  # a round brush along the curve (wide polylines show seams)
        w = i / 400
        x = cx - a * (math.cosh(w) + eps)
        y = cy + a * math.sqrt(eps**2 - 1) * math.sinh(w)
        d.ellipse([(x - 7) * k, (y - 7) * k, (x + 7) * k, (y + 7) * k], fill=(96, 165, 250, 255))
    d.ellipse([(cx - 24) * k, (cy - 24) * k, (cx + 24) * k, (cy + 24) * k], fill=(248, 113, 113, 255))
    img = img.resize((n, n), Image.LANCZOS)
    img.save(path, sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    return True


#: Files a KSHELL folder from scripts/build_kshell.py must have before it is bundled: the programs physim runs,
#: the build record, and the licence texts (KSHELL is GPLv3, OpenBLAS BSD-3: shipping them without is not allowed).
KSHELL_REQUIRED = ["kshell.exe", "transit.exe", "VERSION", "COPYING", "README.txt", "source/kshell.tar.gz"]


def kshell_problems(folder: Path) -> list[str]:
    """Why a KSHELL folder cannot be bundled (empty when it can)."""
    if not folder.is_dir():
        return [f"{folder} not found (build it with: python scripts/build_kshell.py)"]
    problems = [f"{name} missing" for name in KSHELL_REQUIRED if not (folder / name).is_file()]
    if (folder / "libopenblas.dll").is_file() and not (folder / "LICENSE-OpenBLAS.txt").is_file():
        problems.append("LICENSE-OpenBLAS.txt missing")
    return problems


def bundle_kshell(source: Path, out: Path) -> bool:
    """Copy KSHELL (a separate GPLv3 program, run as an executable) to <bundle>/tools/kshell/ when it is there and
    complete; otherwise warn and leave it out (the planner then finds no bundled KSHELL, as before)."""
    problems = kshell_problems(source)
    if problems:
        print("WARNING: KSHELL not bundled: " + "; ".join(problems))
        return False
    dest = out / "tools" / "kshell"
    shutil.copytree(source, dest)
    print("KSHELL bundled:", (source / "VERSION").read_text(encoding="utf-8").splitlines()[0])
    return True


def write_notices(out: Path, kshell: bool) -> None:
    """THIRD-PARTY-NOTICES.txt, shown by the installer before it installs (planner.iss, InfoBeforeFile)."""
    text = ["Third-party software in the physim planner",
            "===========================================",
            "",
            "physim itself is under the MIT licence (the licence page before this one). The installer also contains",
            "separately licensed software:",
            "",
            "Python and the packages the planner uses",
            "  python\\ is the official embeddable Python from python.org (Python Software Foundation licence,",
            "  python\\LICENSE.txt) with packages under their own licences (python\\Lib\\site-packages\\*.dist-info).",
            ""]
    if kshell:
        version = (out / "tools" / "kshell" / "VERSION").read_text(encoding="utf-8").splitlines()[0]
        text += ["KSHELL (tools\\kshell\\), for shell-model calculations",
                 f"  {version}.",
                 "  KSHELL is by N. Shimizu, T. Mizusaki, Y. Utsuno and Y. Tsunoda, Comput. Phys. Commun. 244, 372",
                 "  (2019). It is free software under the GNU General Public License, version 3 or later",
                 "  (tools\\kshell\\COPYING), and comes WITHOUT ANY WARRANTY. It is a separate program: physim runs",
                 "  its executables and reads their output, and is not linked with it. The complete source, physim's",
                 "  small patch to it and the build script are in tools\\kshell\\source\\.",
                 "  Its run-time libraries are OpenBLAS (BSD licence, tools\\kshell\\LICENSE-OpenBLAS.txt) and the GCC",
                 "  and MinGW-w64 run-time libraries (GPL v3 with the GCC Runtime Library Exception; MIT/BSD):",
                 "  see tools\\kshell\\README.txt.",
                 ""]
    (out / "THIRD-PARTY-NOTICES.txt").write_text("\n".join(text), encoding="utf-8", newline="\r\n")


def check(out: Path) -> None:
    """The bundled Python finds physim and the app's packages in the bundle, and nothing outside it."""
    code = ("import sys, physim, nicegui, plotly, matplotlib, uproot, physim.nuclear.app as a; "
            "assert all(p.lower().startswith(sys.prefix.lower()) for p in sys.path if p), sys.path; "
            "print('bundle ok:', physim.__version__, physim.__file__)")
    env = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON")}
    subprocess.run([str(out / "python" / "python.exe"), "-c", code], check=True, env=env)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--wheel", type=Path, required=True, help="the physim wheel for win_amd64")
    ap.add_argument("--out", type=Path, default=Path("build/planner"))
    ap.add_argument("--kshell", type=Path, default=REPO / "tools" / "kshell",
                    help="the KSHELL folder from scripts/build_kshell.py, bundled when present")
    args = ap.parse_args(argv)
    if sys.version_info[:2] != tuple(int(x) for x in PYTHON_VERSION.split(".")[:2]):
        ap.error(f"run with Python {PYTHON_VERSION.rsplit('.', 1)[0]} (it precompiles the bundled packages)")

    out = args.out
    shutil.rmtree(out, ignore_errors=True)
    python = out / "python"
    fetch_python(python)
    enable_site_packages(python)
    site = python / "Lib" / "site-packages"
    install(args.wheel, site)
    prune(site)
    compileall.compile_dir(str(site), quiet=1, workers=0)
    make_icon(out / "planner.ico")
    shutil.copy(REPO / "LICENSE", out / "LICENSE.txt")
    write_notices(out, bundle_kshell(args.kshell, out))
    check(out)
    size = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print(f"bundle: {size / 1e6:.0f} MB in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
