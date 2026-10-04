"""Assemble the Windows planner bundle: the official embeddable Python, the physim wheel with the app's packages,
and the icon. installer/windows/planner.iss then packs the folder into an installer.

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
import urllib.request
import zipfile
from pathlib import Path

#: The embeddable distribution from python.org. Its python.exe and pythonw.exe are signed by the Python Software
#: Foundation, so the bundle contains no unsigned executables of ours.
PYTHON_VERSION = "3.12.10"
EMBED_URL = "https://www.python.org/ftp/python/{v}/python-{v}-embed-amd64.zip"

#: Large parts of the bundled packages that the planner never loads.
PRUNE = ["numpy/_core/tests", "numpy/tests", "numpy/*/tests", "matplotlib/tests", "matplotlib/mpl-data/sample_data",
         "pandas/tests", "plotly/tests", "pip", "setuptools"]


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
    cmd = [sys.executable, "-m", "pip", "install", "--no-cache-dir", "--disable-pip-version-check",
           "--target", str(site), "--platform", "win_amd64", "--python-version", f"{major}.{minor}",
           "--implementation", "cp", "--only-binary=:all:", f"{wheel.resolve()}[app,root]"]
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
    shutil.copy(Path(__file__).resolve().parents[2] / "LICENSE", out / "LICENSE.txt")
    check(out)
    size = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print(f"bundle: {size / 1e6:.0f} MB in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
