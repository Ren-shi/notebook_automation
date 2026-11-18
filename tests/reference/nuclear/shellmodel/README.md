# KSHELL fixtures (shell-model wave 0)

Real KSHELL inputs and outputs, so the parser, partition and size-estimate code can be tested without running
KSHELL. Produced on 18 November 2026 on Windows 10 (Ryzen 5 3600, 12 threads) with the KSHELL 2019 release
published with Comput. Phys. Commun. 244, 372 (2019) (Mendeley Data `kgzdz6ryyk`, `kshell.tar.gz`, SHA-256
`ce0659c3…7526`, GPLv3), built with MinGW-w64 gfortran 16.2 (`-O3 -fopenmp -fallow-argument-mismatch`) against
OpenBLAS 0.3.34, with physim's rename patch (below). KSHELL ver.4 (2023-07-14), built the same way, gave
byte-identical results in the summary for every ²⁸Si case.

Files are kept byte for byte (`.gitattributes`: `-text`): KSHELL's logs on Windows mix CRLF and LF line endings,
and the parsers must read them as written.

## The Windows rename patch

KSHELL's Lanczos restart with vectors on disk (`mode_lv_hdd = 1`, its default) renames temporary files onto names
that exist. POSIX `rename` replaces the target; on Windows it fails, KSHELL reads stale vectors back and stops with
"possible JJ collapse". The patch adds `rename_replace(old, new)` to `lanczos.f90` and `bp_io.F90` that behaves
as POSIX does (replaces `new`; does nothing when `old` is missing or `old == new`). With it, disk and memory
modes give identical eigenvalues. Without it, only `mode_lv_hdd = 0` works.

## Contents

- `kshell/si28_usdb/` — ²⁸Si in the sd shell with USDB: lowest three 0⁺, three 2⁺, two 4⁺ (J-projected,
  `is_double_j = .true.`, so `mtot` is 2J), and E2/M1 transitions between them.
  - `inputs/`: the namelist inputs exactly as KSHELL ver.4's `kshell_ui.py` wrote them (`Si28_usdb_<2J>.input`
    for kshell.exe, `Si28_usdb_<2Ji>_<2Jf>.input` for transit.exe) and the partition file `Si28_usdb_p.ptn`.
  - `logs/`: `log_Si28_usdb_j<2J>p.txt` (kshell.exe), `log_Si28_usdb_tr_j<2J>p_j<2J>p.txt` (transit.exe), with
    the default effective charges e_p = 1.5, e_n = 0.5 and free g-factors; `summary_Si28_usdb.txt` is KSHELL's
    `collect_logs.py` output for the same logs.
  - `logs_eff_1.31_0.46/`: the transit runs again with e_p = 1.31, e_n = 0.46 (the input is
    `transit_template.input` with the two file names filled in), for comparison with the paper below.
- `kshell/cr48_kb3/` — ⁴⁸Cr in the fp shell with KB3: KSHELL's own test input (`test/cr48.sh`: block Lanczos,
  32 states at M = 0, no J projection, vectors in memory) and its log. 6.7 min on 12 threads.
- `partitions/` — partition files from KSHELL ver.4's `gen_partition.py` with no truncation (answer `0`), named
  `<nucleus>_<interaction>_<p|n>.ptn` (positive or negative parity), and `count_dim.exe`'s output for each
  (`*.count_dim.txt`: M-scheme and J-scheme dimension for every 2M). 14 cases over ckpot, usdb, kb3g, gxpf1a,
  jun45 and sdpf-mu, odd and even, both parities where the space allows.
- `benchmarks.csv` — the published numbers below as a reference file (case, logs folder, quantity, state,
  value as printed, unit, source with paper and table), read by the validation check
  `physim.nuclear.validation.shellmodel_vs_published`. It allows one unit in the last digit the source prints
  (dimensions exact). Added in wave 2 with the ⁴⁸Cr excited states from the same paper's Table 2 (exact:
  2⁺₁ −32.148, 4⁺₁ −31.130, 6⁺₁ −29.555, 0⁺₂ −28.564 MeV; KSHELL −32.14801, −31.13038, −29.55518, −28.56458).
  B(E2; 2⁺₁→0⁺₁), 78.20 against 78.3, is the closest to the limit (99% of it); the other ²⁸Si numbers agree within
  half a unit. An oscillator length 0.03% longer (B(E2) × 1.0013, Q_s × 1.0006) would put all of them within
  half a unit; the paper does not state the constants it used.
- `si28_ensdf.csv` — ENSDF's adopted ²⁸Si levels below 10 MeV and the B(E2) physim derives from them
  (`LevelScheme.from_ensdf("28Si", max_energy_kev=10000)`, ENSDF file dated 2013-10), so the shell model can be
  compared with experiment (the check `shellmodel_dimensions` reports it as context) without the local copy.
- `snt/` — the interaction files those cases use, from the same KSHELL 2019 distribution (GPLv3, see below).
  Note: the comment header of `usdb.snt` says "USDA" (copied from `usda.snt`), but its numbers are USDB's
  (single-particle energies 2.1117, −3.9257, −3.2079 MeV).

## Benchmarks these outputs pass

| Case | Quantity | KSHELL here | Published | Source |
|---|---|---|---|---|
| ²⁸Si USDB | E(0⁺₁) relative to ¹⁶O | −135.861 MeV | −135.9 MeV | Frycz et al., arXiv:2404.14506 (2024), Table 3, ANTOINE |
| ²⁸Si USDB, e_p 1.31, e_n 0.46 | B(E2; 2⁺₁→0⁺₁) | 78.20 e²fm⁴ | 78.3 e²fm⁴ | same, Table 5 |
| | B(E2; 4⁺₁→2⁺₁) | 110.02 e²fm⁴ | 110 e²fm⁴ | same, Table 5 |
| | Q_s(2⁺₁) | +18.50 e fm² | +18.5 e fm² | same, Table 2 |
| ²⁸Si sd shell | M-scheme dimension, M = 0 | 93 710 | 93 710 | standard |
| ⁴⁸Cr KB3 | E(0⁺₁) | −32.95367 MeV | −32.953 MeV | "The Density Matrix Renormalization Group and the Nuclear Shell Model", arXiv:0808.1303 (exact full-fp result quoted there) |
| ⁴⁸Cr fp shell | M-scheme dimension, M = 0 | 1 963 461 | 1 963 461 | same |
| ⁵⁶Ni fp shell | M-scheme dimension, M = 0 | 1 087 455 228 | 1 087 455 228 | standard |

Conventions seen in the logs (for the parser): transit.exe's `Mred.` is the signed reduced matrix element
⟨f‖M(Eλ)‖i⟩ and `B(EM )->` = Mred² / (2J_i + 1) with i the right-hand (second) state, which is physim's
convention in `levels.py`; `Mom.` on a diagonal line is the spectroscopic moment (Q_s in e fm², μ in μ_N).
kshell.exe prints `<eQ>` (Q_s with the effective charges) and the g-factor inputs for every state.

## Licence of these files

`snt/*.snt` and the partition files' format come from KSHELL (Copyright 2019 Noritaka Shimizu), distributed under
the GNU General Public License v3 or later (header of `kshell.F90`; Mendeley Data lists the release as GPLv3).
They are kept here as test data under that licence, separately from physim's own MIT-licensed code. The logs are
KSHELL's output for the inputs above. Cite N. Shimizu, T. Mizusaki, Y. Utsuno and Y. Tsunoda, Comput. Phys.
Commun. 244, 372 (2019) when using results from KSHELL.

Remade by `tools/kshell-build/assemble_fixtures.py` (local, git-ignored, wave 0); wave 1's
`scripts/build_kshell.py` takes over building and running.
