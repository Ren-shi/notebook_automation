"""GEMINI++ reference results for physim's evaporation Monte Carlo (validation step B5).

GEMINI++ (R. J. Charity) is a second, independent statistical-decay code beside PACE4. LISE++ ships a GUI build of
version 3.3.0 as ``Gemini.exe``. It has no silent mode: given a ``.gemini`` file on the command line it opens the file
in its window and waits for Execute. ``scripts/gemini_drive.ps1`` presses Execute through Windows UI Automation (no
mouse, no user), reads the results window's tabs and closes the program. GEMINI++ itself writes ``<name>.html`` (the
summary: residue cross sections, mean multiplicities) and ``<name>.cs4`` next to the input. The mean CM energies of
n, p and α are only in the results window's second tab, which the driver saves as ``<name>.tab1.txt``.

The five systems are the PACE4 comparison of ``docs/theory/evaporation.md``. GEMINI++ runs in fusion mode at the same
lab energy, so it computes the compound nucleus's E* itself (AME2016 masses); the file records it beside PACE4's. Its
spin distribution is a Fermi function, σ_ℓ ∝ (2ℓ+1)/(1 + exp((ℓ − L0)/Δ)), with Δ = 2 ħ as PACE4's L-diffuseness.
PACE4's σ_J is the same form (Bass model), so L0 is fitted to PACE4's σ_J. GEMINI++ takes the critical spin as an
integer and derives its own L0 from it (about 0.15 ħ lower); short probe runs pick the integer whose L0 is nearest
the fit, within 0.5 ħ. GEMINI++'s random seed is fixed: identical inputs give identical results.

Work happens in a scratch copy of the LISE++ install; nothing is written into the install. Run on Windows:

    python scripts/gemini_reference.py ["D:/Program Files/LISEcute/LISEcute"] [--work DIR] [--events 20000]

One CSV per system goes to ``tests/reference/nuclear/gemini_<system>.csv`` in the format of that folder's README.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import math
import re
import shutil
import subprocess
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "gemini_drive.ps1"
OUT = ROOT / "tests" / "reference" / "nuclear"
DEFAULT_LISE = Path(r"D:\Program Files\LISEcute\LISEcute")
DIFFUSENESS = 2.0  # ħ, PACE4's L-diffuseness
# GEMINI++'s evaporation mode: 0 = widths (and kinetic energies) from Weisskopf, spin and ℓ from Hauser-Feshbach;
# 1 = Hauser-Feshbach; 2 = switches between 0 and 1 by the ratio of rotational to thermal energy. In this GUI build
# mode 1 hangs (and later crashes, 0xC0000005) on some event of 16O + 144Sm (between the 1500th and 2000th), and of
# 32S + 92Mo and 48Ca + 124Sn within the first events at PACE4's spins (also with GEMINI's own Bass spins, and for
# 48Ca still at critical spin 50); mode 0 hangs on 16O + 144Sm too. Mode 2 finishes but silently loses events
# (4He + 165Ho: all 20; 16O + 144Sm: 18863 of 20000; 32S: 40%; 48Ca: 48%), so it is not used. Mode 1 it is; systems
# it cannot run get no file.
EVAP_MODES = {0: "0 (Weisskopf widths, Hauser-Feshbach spins)", 1: "1 (Hauser-Feshbach)",
              2: "2 (switches between Weisskopf and Hauser-Feshbach by the ratio of rotational to thermal energy)"}
EVAP_MODE = 1

# PACE4 4.36.0 runs of the same systems (scripts/pace4_reference.py, tests/reference/nuclear/pace4_<system>.csv):
# E* in MeV and σ_J in mb for J = 0, 1, 2, ... (beam and target spins 0, so J = ℓ).
SYSTEMS = [
    dict(name="4He_165Ho_50MeV", beam=(2, 4), target=(67, 165), e_lab=50.0, pace4_ex=47.618,
         pace4_sigma_j=[3.4, 10, 17, 24, 31, 38, 45, 52, 58, 65, 72, 79, 85, 92, 98, 100, 110, 110, 110, 100, 95, 83,
                        67, 51, 37, 26, 17, 11, 7.2, 4.6, 2.9, 1.8, 1.2, 0.72, 0.45, 0.28]),
    dict(name="12C_150Nd_70MeV", beam=(6, 12), target=(60, 150), e_lab=70.0, pace4_ex=59.317,
         pace4_sigma_j=[0.91, 2.7, 4.6, 6.4, 8.2, 10, 12, 14, 15, 17, 19, 21, 23, 25, 26, 28, 30, 32, 34, 36, 37, 39,
                        41, 43, 44, 46, 47, 48, 48, 47, 45, 41, 35, 29, 22, 16, 11, 7.4, 4.8, 3.1, 2, 1.2, 0.77, 0.48,
                        0.3, 0.19, 0.12]),
    dict(name="16O_144Sm_90MeV", beam=(8, 16), target=(62, 144), e_lab=90.0, pace4_ex=52.46,
         pace4_sigma_j=[0.56, 1.7, 2.8, 3.9, 5.1, 6.2, 7.3, 8.4, 9.6, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23,
                        24, 25, 26, 28, 29, 30, 31, 32, 33, 33, 34, 34, 34, 32, 30, 27, 22, 18, 13, 9.2, 6.2, 4.1,
                        2.6, 1.7, 1, 0.65, 0.41, 0.25, 0.16, 0.097, 0.06]),
    dict(name="32S_92Mo_150MeV", beam=(16, 32), target=(42, 92), e_lab=150.0, pace4_ex=63.382,
         pace4_sigma_j=[0.25, 0.75, 1.2, 1.7, 2.2, 2.7, 3.2, 3.7, 4.2, 4.7, 5.2, 5.7, 6.2, 6.7, 7.2, 7.7, 8.2, 8.7,
                        9.2, 9.7, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20,
                        21, 21, 22, 22, 22, 23, 23, 23, 23, 23, 23, 21, 19, 17, 14, 10, 7.6, 5.2, 3.5, 2.2, 1.4, 0.9,
                        0.56, 0.35, 0.21, 0.13, 0.082, 0.05]),
    dict(name="48Ca_124Sn_200MeV", beam=(20, 48), target=(50, 124), e_lab=200.0, pace4_ex=70.983,
         pace4_sigma_j=[0.13, 0.39, 0.66, 0.92, 1.2, 1.4, 1.7, 2, 2.2, 2.5, 2.8, 3, 3.3, 3.6, 3.8, 4.1, 4.3, 4.6, 4.9,
                        5.1, 5.4, 5.7, 5.9, 6.2, 6.4, 6.7, 7, 7.2, 7.5, 7.8, 8, 8.3, 8.6, 8.8, 9.1, 9.3, 9.6, 9.9, 10,
                        10, 11, 11, 11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16,
                        16, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 17, 16, 15, 13, 10, 7.8, 5.7, 3.9, 2.6, 1.6,
                        1, 0.65, 0.4, 0.25, 0.15, 0.095, 0.058, 0.036]),
]


def fermi_l0(sigma_j: list[float], dif: float = DIFFUSENESS) -> float:
    """Half-point L0 of the Fermi distribution N(2J+1)/(1+exp((J-L0)/dif)) that fits σ_J best (least squares)."""

    def residual(l0: float) -> float:
        shape = [(2 * j + 1) / (1.0 + math.exp((j - l0) / dif)) for j in range(len(sigma_j))]
        norm = sum(s * f for s, f in zip(sigma_j, shape)) / sum(f * f for f in shape)
        return sum((s - norm * f) ** 2 for s, f in zip(sigma_j, shape))

    lo, hi = 0.0, float(len(sigma_j))
    g = (math.sqrt(5) - 1) / 2
    for _ in range(100):  # golden-section search
        a, b = hi - g * (hi - lo), lo + g * (hi - lo)
        if residual(a) < residual(b):
            hi = b
        else:
            lo = a
    return 0.5 * (lo + hi)


def gemini_input(system: dict, events: int, critical_spin: float, evap_mode: int = EVAP_MODE) -> str:
    (zp, ap), (zt, at) = system["beam"], system["target"]
    # The [Common], [Spin] and [Compound] keys follow LISE++'s example files\examples\40Ca_112Sn.gemini. Which
    # calculation runs follows the selected tab of the window (the driver keeps "Fusion reaction"); the
    # [Compound] section is for the other tab and is not used here.
    return "\n".join([
        "[Common]", "mode=1", "AME=true", "IMF=true", "IMFenh=false", f"optEvap={evap_mode}", "",
        "[Projectile]", f"A={ap}", f"Z={zp}", f"E={system['e_lab']:g}", f"num_events={events}", "",
        "[Target]", f"A={at}", f"Z={zt}", "",
        "[Spin]", "S=0", f"l0={int(critical_spin)}", f"dif={DIFFUSENESS:g}", "",
        "[Compound]", f"A={ap + at}", f"Z={zp + zt}", f"Ex={system['pace4_ex']:g}", "J=30", "num_casc=300", "",
    ])


class _Tables(HTMLParser):
    """All table rows of an HTML page, as lists of cell texts."""

    def __init__(self):
        super().__init__()
        self.rows: list[list[str]] = []
        self._row: list[str] | None = None
        self._cell: list[str] | None = None

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self._row = []
        elif tag in ("td", "th") and self._row is not None:
            if self._cell is not None:  # GEMINI++ leaves some <td> unclosed
                self.handle_endtag(tag)
            self._cell = []

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self._cell is not None:
            self._row.append(" ".join("".join(self._cell).split()))
            self._cell = None
        elif tag == "tr" and self._row is not None:
            self.rows.append(self._row)
            self._row = None

    def handle_data(self, data):
        if self._cell is not None:
            self._cell.append(data)


def _value(rows: list[list[str]], label: str) -> str:
    for row in rows:
        for i, cell in enumerate(row[:-1]):
            if cell.strip().lower() == label.lower():
                return row[i + 1]
    raise ValueError(f"'{label}' not in the GEMINI++ output")


def _number(text: str) -> float:
    return float(re.match(r"\s*([-+0-9.eE]+)", text).group(1))


def parse_html(path: Path) -> dict:
    t = _Tables()
    t.feed(path.read_text(encoding="utf-8", errors="replace"))
    rows = t.rows
    out = {
        "excitation_mev": _number(_value(rows, "Compound nucleus Excitation energy (MeV)")),
        "sigma_fus_mb": _number(_value(rows, "Fusion cross section")),
        "l0": _number(_value(rows, "L0")),
        "critical_spin": _number(_value(rows, "Critical spin")),
        "bass_l": _number(_value(rows, "Bass L")),
        "imf_events": int(_value(rows, "Intermediate Mass Fragments")),
        "fission_events": int(_value(rows, "Symmetric Fission")),
        "residue_events": int(_value(rows, "Residual Nuclei")),
        "total_events": int(_value(rows, "TOTAL")),
    }
    try:
        out["gamma_energy_mev"] = _number(_value(rows, "average energy in gamma rays"))
    except ValueError:  # absent when no event ended as a residue
        out["gamma_energy_mev"] = float("nan")
    m = re.search(r"num_events=(\d+)", path.with_suffix(".gemini").read_text()) if path.with_suffix(
        ".gemini").exists() else None
    if m and out["total_events"] != int(m.group(1)):
        raise EventsLost(f"{path.name}: GEMINI++ accounts for {out['total_events']} of {m.group(1)} events")
    for row in rows:  # "Version 3.3.0" | build date
        for i, cell in enumerate(row[:-1]):
            if cell.startswith("Version "):
                out["version"] = f"{cell.split()[1]} of {row[i + 1]}"
    residues, z = [], None
    for row in rows:
        if len(row) >= 6 and re.fullmatch(r"\d+ [A-Z][a-z]?", row[1]) and re.fullmatch(r"\d+", row[2]):
            if row[0].strip():
                z = int(row[0])
            residues.append((z, int(row[1].split()[0]), int(row[2]), float(row[4])))
    out["residues"] = sorted(residues, key=lambda r: -r[2])
    # the .cs4 file beside it has E* to 1 keV: "!Gemini Zcomp Ncomp Energy CSfus ..."
    cs4 = path.with_suffix(".cs4")
    if cs4.exists():
        line = [x for x in cs4.read_text().splitlines() if x.startswith("!Gemini")][0]
        out["excitation_mev"] = float(line.split()[3])
    return out


def parse_cm_spectra(path: Path, events: int) -> dict:
    """Multiplicities and mean CM energies of n, p and α from the 'C.M. spectra of emitted particles' table.

    Its columns are Neut, Prot, Alpha and Gamma; the Gamma column histograms the summed γ energy per event (its
    total is the number of events with γ rays), so only the particle columns are used."""
    text = path.read_text(encoding="utf-8-sig")
    block = text[text.index("C.M. spectra of emitted particles"):]
    cells = [" ".join(c.split()) for c in re.split("[\ufdd0\ufdd1]", block)]  # also folds non-breaking spaces

    def row(label: str) -> list[float]:
        i = cells.index(label)
        return [float(c) if c.strip() else 0.0 for c in cells[i + 1:i + 5]]

    header = cells[cells.index("Ex(MeV)") + 1:cells.index("Ex(MeV)") + 5]
    if header != ["Neut", "Prot", "Alpha", "Gamma"]:
        raise ValueError(f"unexpected C.M. spectra columns {header}")
    total, mean = row("Total"), row("Average Energy")
    keys = ("n", "p", "a", "gamma")
    return {"multiplicity": {k: v / events for k, v in zip(keys, total)},
            "mean_cm_energy_mev": dict(zip(keys, mean))}


HIDDEN_DESKTOP = "physim_gemini"


def run_on_hidden_desktop(cmdline: str, timeout_s: float) -> int:
    """Run a command on a separate Windows desktop of this window station, which is never shown.

    Its windows, and those of the programs it starts (they inherit the desktop), never appear on the user's screen
    or take the focus, and UI Automation called from that desktop sees them. Returns the exit code."""
    import ctypes
    from ctypes import wintypes

    class STARTUPINFOW(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("lpReserved", wintypes.LPWSTR), ("lpDesktop", wintypes.LPWSTR),
                    ("lpTitle", wintypes.LPWSTR), ("dwX", wintypes.DWORD), ("dwY", wintypes.DWORD),
                    ("dwXSize", wintypes.DWORD), ("dwYSize", wintypes.DWORD), ("dwXCountChars", wintypes.DWORD),
                    ("dwYCountChars", wintypes.DWORD), ("dwFillAttribute", wintypes.DWORD),
                    ("dwFlags", wintypes.DWORD), ("wShowWindow", wintypes.WORD), ("cbReserved2", wintypes.WORD),
                    ("lpReserved2", ctypes.c_void_p), ("hStdInput", wintypes.HANDLE),
                    ("hStdOutput", wintypes.HANDLE), ("hStdError", wintypes.HANDLE)]

    class PROCESS_INFORMATION(ctypes.Structure):
        _fields_ = [("hProcess", wintypes.HANDLE), ("hThread", wintypes.HANDLE),
                    ("dwProcessId", wintypes.DWORD), ("dwThreadId", wintypes.DWORD)]

    user32 = ctypes.WinDLL("user32", use_last_error=True)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    user32.CreateDesktopW.restype = wintypes.HANDLE
    user32.CreateDesktopW.argtypes = [wintypes.LPCWSTR, wintypes.LPCWSTR, ctypes.c_void_p, wintypes.DWORD,
                                      wintypes.DWORD, ctypes.c_void_p]
    kernel32.CreateProcessW.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, ctypes.c_void_p, ctypes.c_void_p,
                                        wintypes.BOOL, wintypes.DWORD, ctypes.c_void_p, wintypes.LPCWSTR,
                                        ctypes.POINTER(STARTUPINFOW), ctypes.POINTER(PROCESS_INFORMATION)]
    desktop = user32.CreateDesktopW(HIDDEN_DESKTOP, None, None, 0, 0x01FF, None)  # DESKTOP_ALL_ACCESS
    if not desktop:
        raise OSError(ctypes.get_last_error(), "CreateDesktop failed")
    try:
        si = STARTUPINFOW(cb=ctypes.sizeof(STARTUPINFOW), lpDesktop=HIDDEN_DESKTOP)
        pi = PROCESS_INFORMATION()
        buf = ctypes.create_unicode_buffer(cmdline)
        if not kernel32.CreateProcessW(None, buf, None, None, False, 0x08000000, None, None,  # CREATE_NO_WINDOW
                                       ctypes.byref(si), ctypes.byref(pi)):
            raise OSError(ctypes.get_last_error(), "CreateProcess failed")
        try:
            if kernel32.WaitForSingleObject(pi.hProcess, int(timeout_s * 1000)) != 0:
                kernel32.TerminateProcess(pi.hProcess, 1)
                return -1
            code = wintypes.DWORD()
            kernel32.GetExitCodeProcess(pi.hProcess, ctypes.byref(code))
            return code.value
        finally:
            kernel32.CloseHandle(pi.hThread)
            kernel32.CloseHandle(pi.hProcess)
    finally:
        user32.CloseDesktop(desktop)


def _kill_gemini(exe: Path) -> None:
    """Stop any Gemini.exe left running from the scratch copy (after a hang or a killed driver)."""
    ps = (f"Get-Process Gemini -ErrorAction SilentlyContinue | Where-Object {{ $_.Path -eq '{exe}' }} "
          "| Stop-Process -Force")
    subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True)


def run_gemini(exe: Path, folder: Path, name: str, text: str, visible: bool = False) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.glob(f"{name}.*"):
        old.unlink()
    (folder / f"{name}.gemini").write_text(text, encoding="ascii")
    log = folder / f"{name}.driver.log"
    args = ["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(DRIVER), "-Exe", str(exe),
            "-InputFile", str(folder / f"{name}.gemini"), "-Tab", "Fusion reaction", "-TimeoutSeconds", "1800",
            "-LogFile", str(log)]
    for attempt in range(3):  # a run that failed (window did not come up, program crashed) is started afresh
        if log.exists():
            log.unlink()
        if visible:
            code = subprocess.run(args, capture_output=True, text=True).returncode
        else:
            code = run_on_hidden_desktop(subprocess.list2cmdline(args), timeout_s=2000)
        _kill_gemini(exe)
        message = log.read_text(encoding="utf-8-sig").strip() if log.exists() else "(no driver log)"
        if code == 0 and (folder / f"{name}.html").exists():
            print("  " + message)
            return
        print(f"  attempt {attempt + 1} failed: {message}")
        if "hung" in message or "exited during the run" in message:
            break  # GEMINI++'s seed is fixed: the same input fails the same way again
    raise GeminiFailed(f"GEMINI++ run {name} failed: {message}")


class GeminiFailed(RuntimeError):
    pass


class EventsLost(GeminiFailed):
    """The run finished, but its outcomes (residues + fission + IMF) do not add up to the events asked for."""


def scratch_install(lise: Path, work: Path) -> Path:
    copy = work / "LISEcute"
    if not (copy / "Gemini.exe").exists():
        print(f"copying {lise} -> {copy}")
        shutil.copytree(lise, copy, dirs_exist_ok=True)
    return copy / "Gemini.exe"


def write_csv(system: dict, res: dict, cm: dict, events: int, l0_fit: float, crit: float, text: str,
              note: str = "") -> Path:
    (zp, ap), (zt, at) = system["beam"], system["target"]
    sym = system["name"].split("_")
    n, m = cm["mean_cm_energy_mev"], cm["multiplicity"]
    path = OUT / f"gemini_{system['name']}.csv"
    lines = [
        f"# tool: GEMINI++ {res['version']} (GUI build in LISE++ (LISEcute), Gemini.exe), in a scratch copy of the "
        "install; Execute pressed through Windows UI Automation by scripts/gemini_drive.ps1 on a hidden desktop (no "
        "silent mode exists)",
        f"# produced by: scripts/gemini_reference.py (Claude Code), {_dt.date.today().isoformat()}",
        f"# system: {sym[0]} + {sym[1]} at {system['e_lab']:g} MeV (lab) -> compound nucleus Z = {zp + zt}, "
        f"A = {ap + at}",
        "# settings: fusion mode at PACE4's lab energy (E* from GEMINI's own AME2016 masses); spin option 0 (critical "
        f"spin input) with Fermi diffuseness {DIFFUSENESS:g} hbar; the integer critical spin (GEMINI takes no",
        f"#   decimals) whose L0 is nearest the half-point of a Fermi fit to PACE4's sigma_J (fit L0 = {l0_fit:.2f}); "
        f"evaporation mode {EVAP_MODES[EVAP_MODE]};",
        "#   masses AME2016; IMF emission on (GEMINI default), not enhanced;",
        f"#   {events} events; GEMINI's random seed is fixed (identical inputs give identical results). "
        "Multiplicities and mean CM energies from GEMINI's 'C.M. spectra of emitted particles' table",
        "# input file: " + " ".join(line for line in text.splitlines() if line),
        f"# excitation_mev: {res['excitation_mev']:g}",
        f"# pace4_excitation_mev: {system['pace4_ex']:g}",
        f"# sigma_fus_mb: {res['sigma_fus_mb']:g}",
        f"# spin: critical_spin_hbar={crit:.2f},l0_hbar={res['l0']:g},l0_fit_to_pace4_hbar={l0_fit:.3f},"
        f"diffuseness_hbar={DIFFUSENESS:g},bass_l_hbar={res['bass_l']:g}",
        f"# mean_cm_energy_mev: n={n['n']:g},p={n['p']:g},a={n['a']:g}",
        f"# multiplicity: n={m['n']:g},p={m['p']:g},a={m['a']:g}",
        f"# gamma_energy_per_event_mev: {res['gamma_energy_mev']:g}",
        f"# outcomes: residues={res['residue_events']},imf={res['imf_events']},fission={res['fission_events']}",
        f"# cascades: {events}",
    ]
    if note:
        lines.append(f"# note: {note}")
    lines.append("Z,A,events,sigma_mb")
    lines += [f"{z},{a},{ev},{sig:g}" for z, a, ev, sig in res["residues"]]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run_system(system: dict, exe: Path, runs: Path, args) -> Path:
    """Probe the critical spin, run GEMINI++, write the CSV."""
    l0_fit = fermi_l0(system["pace4_sigma_j"])
    # GEMINI reads the critical spin as an integer (a decimal, in the file or the window, becomes 0) and puts
    # its L0 about 0.15 hbar below it. Probe the nearest integers and keep the one whose L0 is closest.
    guess = round(l0_fit + 0.15)
    probes = {}
    for crit in (guess - 1, guess, guess + 1):
        try:
            run_gemini(exe, runs, f"{system['name']}_probe", gemini_input(system, 20, crit), visible=args.visible)
            probes[crit] = parse_html(runs / f"{system['name']}_probe.html")["l0"]
        except GeminiFailed:
            continue
    if not probes:
        raise GeminiFailed(f"{system['name']}: every probe run failed")
    crit = min(probes, key=lambda c: abs(probes[c] - l0_fit))
    # GEMINI++ hangs, then crashes, on some event of some inputs (see EVAP_MODE). Its seed is fixed and the events
    # come in a fixed order, so a shorter run is the same events up to that point: if a run hangs, fall back to
    # fewer events and say so in the file.
    hung = []
    for events in [args.events] + [n for n in (10000, 5000, 2000, 1500, 1000) if n < args.events]:
        text = gemini_input(system, events, crit)
        try:
            run_gemini(exe, runs, system["name"], text, visible=args.visible)
            break
        except GeminiFailed:
            hung.append(events)
    else:
        raise GeminiFailed(f"{system['name']}: runs of {', '.join(map(str, hung))} events all hung or crashed")
    note = (f"runs of {', '.join(map(str, hung))} events hung GEMINI++ (Gemini.exe stops using CPU, later crashes "
            "with 0xC0000005); the largest event count that finished is used" if hung else "")
    res = parse_html(runs / f"{system['name']}.html")
    cm = parse_cm_spectra(runs / f"{system['name']}.tab1.txt", events)
    path = write_csv(system, res, cm, events, l0_fit, crit, text, note)
    top = ", ".join(f"{a}({z}) {ev / events:.3f}" for z, a, ev, _ in res["residues"][:3])
    print(f"  E*={res['excitation_mev']} L0={res['l0']} (fit {l0_fit:.2f}) sigma_fus={res['sigma_fus_mb']} mb; "
          f"{top}; <e> {cm['mean_cm_energy_mev']}; M {cm['multiplicity']} -> {path.name}")
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("lise", nargs="?", type=Path, default=DEFAULT_LISE, help="LISE++ install folder (read only)")
    ap.add_argument("--work", type=Path, default=Path(tempfile.gettempdir()) / "physim_gemini")
    ap.add_argument("--events", type=int, default=20000)
    ap.add_argument("--only", nargs="*", help="system names to run (default: all five)")
    ap.add_argument("--visible", action="store_true",
                    help="run GEMINI++ on the normal desktop (its windows show and take the focus) instead of a "
                         "hidden one")
    args = ap.parse_args()
    if sys.platform != "win32":
        sys.exit("GEMINI++ from LISE++ is a Windows program")
    exe = scratch_install(args.lise, args.work)
    runs = args.work / "runs"
    failed = {}
    for system in SYSTEMS:
        if args.only and system["name"] not in args.only:
            continue
        print(system["name"])
        try:
            run_system(system, exe, runs, args)
        except GeminiFailed as failure:
            failed[system["name"]] = str(failure)
            print(f"  NO FILE: {failure}")
    if failed:
        print("GEMINI++ could not run: " + "; ".join(failed))
        sys.exit(1)

if __name__ == "__main__":
    main()
