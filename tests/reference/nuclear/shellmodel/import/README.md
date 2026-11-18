# Shell-model import fixtures (wave 2B)

Files for `physim.nuclear.shellmodel_import`, which reads shell-model results run elsewhere. The KSHELL logs and
KSHELL's `collect_logs.py` summary it also reads are the wave 0 fixtures one folder up (`../kshell/`, see
`../README.md`); nothing is copied here.

- `si28_usdb_table.csv` — a filled example of physim's shell-model table format: ²⁸Si with USDB and effective
  charges e_p = 1.31, e_n = 0.46, typed by hand from this repository's own KSHELL run (the benchmark table of
  `../README.md` for E(0⁺₁), B(E2) and Q_s(2⁺₁); `../kshell/si28_usdb/logs/summary_Si28_usdb.txt` for the
  excitation energies; the transit logs in `../kshell/si28_usdb/logs_eff_1.31_0.46/` for the signed matrix
  elements and μ(2⁺₁)). Its notes give the published ANTOINE numbers (Frycz et al., arXiv:2404.14506) these agree
  with. The tests check it against the logs it was typed from.
- `table_template.csv` — the template a user copies: every header entry, every kind and unit explained in `##`
  comments, with the same ²⁸Si numbers as `si28_usdb_table.csv` written in other units (keV, e b, e²b², W.u.,
  a transition written upwards). The tests check that it reads to the same values, which tests the unit
  conversions and the turning round of transitions.

No NuShellX output is here: none was available to test a reader against (searched this machine and the
repository on 18 November 2026), so `import_path` raises `NotImplementedError` for `.lpt` files. A NuShellX
reader needs real output of one calculation (the `.lpt` level list and the transition output, with the inputs that
made them) before it is written.
