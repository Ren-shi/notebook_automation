# Shell-model calculations

physim can run a nuclear shell-model calculation for a nucleus in your experiment and put the results into its
level scheme: level energies, B(E2), quadrupole moments Q_s, and (as numbers to read) B(M1) and magnetic moments
μ. The values then feed the Coulomb-excitation yields and the γ-ray lists exactly as ENSDF values or values you
typed do, each marked with the source **shell model** and a note naming the interaction and effective charges.

## What it does and does not do

- **physim drives KSHELL; it does not solve the shell model itself.** KSHELL (N. Shimizu et al., below) is an
  outside program. physim writes its inputs (partition files and namelists), runs it out of sight in the
  background, reads its logs, and records where every number came from.
- **Laptop-sized valence spaces.** The p, sd and fp shells, f5pg9 and sd-pf, with the interactions listed
  [below](#interactions). A big calculation belongs on a cluster; run it there and [import](#importing-results)
  the results.
- **At most half of the free RAM.** Before a run starts, physim estimates its memory need and compares it with
  half of the RAM that is free at that moment (not the installed RAM: 16 GB installed with 10 GB in use leaves
  6 GB free, so the limit is 3 GB). A run over the limit is refused unless you explicitly allow it, and the free
  memory is checked again before each J^π step.
- **No truncation.** Every partition of the valence nucleons is kept. The size estimate tells you whether that
  fits.
- **Shell-model values never override other values unless you ask.** In a level scheme they fill gaps only: an
  ENSDF, derived or user value is kept unless you explicitly ask for it to be replaced. (Values marked
  "assumed", physim's own placeholders, are replaced.)
- **Signs are kept.** Shell-model matrix elements carry real relative signs, unlike ENSDF's, and are stored with
  them.
- Results are kept in a cache, so the same calculation is never run twice. One KSHELL calculation runs at a time;
  others wait in a queue.

States are written J, parity and an index: **2+1** is the lowest 2⁺ state, **3/2-2** the second 3/2⁻. When you
*ask for* states, the number is how many: **0+3, 2+3, 4+2** means the lowest three 0⁺, three 2⁺ and two 4⁺ states.

## Where KSHELL comes from

The Windows installer of the planner bundles a KSHELL build (in `tools\kshell\` of the install folder). physim
looks for KSHELL in this order and takes the first folder that exists:

1. the environment variable `PHYSIM_KSHELL`;
2. `~/.physim/programs.toml`, a line `kshell = "C:/path/to/kshell"`, written by
   `physim.external.programs.set_path("kshell", folder)`;
3. the copy bundled with the installed planner;
4. `tools/kshell/` in a source checkout (build it with `python scripts/build_kshell.py`).

So you can point physim at your own KSHELL build (an MPI build on a cluster, say) without changing anything else.

## From Python

```python
from physim.nuclear import shellmodel, shellmodel_jobs as jobs, shellmodel_spaces as spaces
from physim.nuclear.shellmodel import ShellModelRequest

spaces.options("28Si")                  # the interactions that apply: valence nucleons, parities, dimensions

# States are (2J, parity, how many): the lowest three 0+, three 2+ and two 4+.
request = ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3), (8, 1, 2)), transitions=("E2", "M1"))
spaces.estimate(request)                # dimensions, memory, disk, the limit, fits, a rough time

job = jobs.submit(request)              # returns at once; the calculation runs in a background process
jobs.wait(job)                          # or poll jobs.status(job)
result = jobs.result(job)               # a ShellModelResult

result.state(4, +1, 1).ex_mev                              # the 2+1: excitation energy, Q_s, mu
result.transition("E2", (4, 1, 1), (0, 1, 1)).b_down       # B(E2; 2+1 -> 0+1) in e2fm4
```

To put values into an experiment's level scheme, read the scheme from ENSDF first and then apply the result
through the planner, which saves the setup:

```python
from physim.nuclear.planner import Planner

p = Planner.open_experiment("Si28 coulex")
p.lookup_levels("target", "5 MeV")
p.apply_shell_model("target", result)                       # fills gaps; keeps ENSDF, derived and user values
p.apply_shell_model("target", result, which=["Q", "2+1"])   # only the quadrupole moment of the 2+1
p.apply_shell_model("target", result, which=["E2"], replace=True)  # only on your explicit choice
```

`which` takes quantities (`"energy"`, `"E2"`, `"E1"`, `"Q"`), states (`"2+1"`) and `"new"` (add the states the
scheme has no level for); `shellmodel.match(result, scheme)` shows which shell-model state is paired with which
known level (by J^π, then energy order). Without a level scheme for the role, the calculated states become the
levels, the lowest as the ground state.

Jobs live in `~/.physim/shellmodel/<job>/` (the environment variable `PHYSIM_SHELLMODEL` moves the folder): the
request, the status, `results.json` and KSHELL's inputs and logs. Wave functions are deleted when a job ends. The
cache is kept under 256 MB (`PHYSIM_SHELLMODEL_CACHE_MB`), removing the jobs used longest ago first.

## From an LLM, through the MCP server

physim's MCP server (`physim mcp`, with the `mcp` extra: `pip install "physim-engine[mcp]"`; add it to Claude Code
with `claude mcp add physim -- physim mcp`) has tools for the whole cycle, so you can ask "run a shell-model
calculation for the target with USDB and put the quadrupole moments into the level scheme":

| Tool | Does |
|---|---|
| `shell_model_options(nuclide?, role?)` | the interactions that apply (by default to the open experiment's target), with valence nucleons, parities, default effective charges and where they come from, citations, and a size estimate for a typical request against the memory limit now |
| `estimate_shell_model(nuclide, interaction, states)` | dimensions, memory, disk, the limit (half of the free RAM now), whether it fits, a rough time |
| `run_shell_model(nuclide, interaction, states, transitions?, eff_charge?, allow_over_limit?)` | starts a background job and returns its ID (or the cached result) |
| `shell_model_status(job, wait_seconds?)` | queued, running (which step), done, failed (with the end of KSHELL's log) or cancelled |
| `list_shell_model_jobs()`, `cancel_shell_model(job)` | housekeeping |
| `shell_model_results(job, role?)` | the states and transitions as numbers with units; with a role, compared with the experiment's level scheme and what applying would change or keep |
| `apply_shell_model(job_or_path, role, which?, replace?)` | puts chosen values into the level scheme and saves the experiment |
| `import_shell_model(path, role?)` | reads results run elsewhere (below) |

The two rules are written into the tools' descriptions, which is what the model reads: `allow_over_limit` and
`replace` are set only when you explicitly say so.

## From the app

The panel **Shell model: calculate or import the levels of …**, beside the level scheme of the target and of the
beam, runs, imports and applies calculations in the planner app (see {doc}`planner-app`). Every result read or
applied is also copied into the nucleus's own results as theory rows (level energies, B(E2), B(M1), Q_s and μ),
so it stays on the nucleus's page after the job cache is cleaned (see {doc}`results`).

## Interactions

The interactions physim knows, with the papers to cite (as given in the headers of KSHELL's interaction files) and
the default effective charges (e_p, e_n). g-factors are the free-nucleon values unless stated.
`shell_model_options` gives the same list for a nucleus.

| Name | Space | Core | Default e_p, e_n | Cite |
|---|---|---|---|---|
| `ckpot` | p (0p3/2, 0p1/2) | ⁴He | 1.5, 0.5 (physim's generic default) | S. Cohen and D. Kurath, Nucl. Phys. 73, 1 (1965) |
| `usda` | sd (0d5/2, 1s1/2, 0d3/2) | ¹⁶O | 1.5, 0.5 (generic) | W. A. Richter, S. Mkhize and B. A. Brown, Phys. Rev. C 78, 064302 (2008); the fit itself: B. A. Brown and W. A. Richter, Phys. Rev. C 74, 034315 (2006) |
| `usdb` | sd | ¹⁶O | 1.5, 0.5 (generic) | as USDA (the `usdb.snt` header says "USDA", copied from `usda.snt`; its numbers are USDB's) |
| `kb3g` | fp (0f7/2, 1p3/2, 0f5/2, 1p1/2) | ⁴⁰Ca | 1.5, 0.5 (generic) | A. Poves, J. Sánchez-Solano, E. Caurier and F. Nowacki, Nucl. Phys. A 694, 157 (2001) |
| `gxpf1a` | fp | ⁴⁰Ca | 1.5, 0.5 (from the `.snt` file) | M. Honma, T. Otsuka, B. A. Brown and T. Mizusaki, Eur. Phys. J. A 25, Suppl. 1, 499 (2005) |
| `jun45` | f5pg9 (1p3/2, 0f5/2, 1p1/2, 0g9/2) | ⁵⁶Ni | 1.5, 1.1, and g_s 0.7 × free (from the `.snt` file) | M. Honma, T. Otsuka, T. Mizusaki and M. Hjorth-Jensen, Phys. Rev. C 80, 064323 (2009) |
| `sdpf-mu` | sd-pf | ¹⁶O | 1.35, 0.35 (from the `.snt` file) | Y. Utsuno, T. Otsuka, B. A. Brown, M. Honma, T. Mizusaki and N. Shimizu, Phys. Rev. C 86, 051301(R) (2012) |

For sd-pf, KSHELL's Lawson term (β_cm = 10) removes centre-of-mass excitations. Give other effective charges with
`eff_charge` (they change B(E2) and Q_s, not the energies).

## Importing results

Results run elsewhere are as usable as results run here (`shellmodel_import.import_path(path)` in Python,
`import_shell_model` through the MCP server). physim reads:

- **a folder of KSHELL logs** (or one log, or a list of them): the `kshell.exe` and `transit.exe` logs of a
  calculation, from any KSHELL build;
- **KSHELL's summary** (`summary_*.txt` from `collect_logs.py`) when the logs are not there: energies to 1 keV
  and B values to KSHELL's rounding, without signs, Q_s or μ;
- **a table** in physim's format, for numbers from a paper or a collaborator's email.

NuShellX output is not read yet (no real output files were available to test a reader against): type the numbers
into a table meanwhile.

The table is a CSV file with one value per row, each with its unit. The header lines say what produced the
numbers:

```text
# format: physim shell-model table 1
# nuclide: 28Si
# interaction: USDB
# effective charges: 1.31, 0.46
# reference: A. Author et al., Journal 1, 2 (2024), Table 3
kind,state,final,value,unit
Ex,0+1,,0.0,MeV
Ex,2+1,,1.779,MeV
Q,2+1,,+18.5,e fm2
B(E2),2+1,0+1,78.3,e2fm4
M(E2),4+1,2+1,-31.47,e fm2
```

Kinds: `Ex` and `E` (excitation and absolute energy, MeV or keV), `T`, `Q` (e fm2 or e b), `mu` (mu_N), `B(E2)`,
`B(M1)`, `B(E1)` (B(state → final), also in W.u.) and `M(E2)`, `M(M1)`, `M(E1)` (the signed reduced matrix
element ⟨final‖M‖state⟩, Bohr–Mottelson convention). A transition given as a B value only has no known sign: it is
kept as the positive root and listed as such. Every state a transition names needs an `Ex` or `E` row. A template
is in `tests/reference/nuclear/shellmodel/import/table_template.csv`; the full description is in the docstring of
`physim.nuclear.shellmodel_import`.

## Validation

How KSHELL, as physim builds, runs and reads it, compares with published calculations with the same Hamiltonian
(²⁸Si/USDB against ANTOINE, ⁴⁸Cr/KB3 against exact energies) and with physim's own dimension count is in the
[physics register](physics-register/README.md), row "Shell-model calculations (KSHELL)". The comparison with
experiment (²⁸Si/USDB against ENSDF) is reported there as information, not as a test: the shell model is a
model.

## Licence and citation

KSHELL is free software under the **GNU General Public License, version 3**. physim bundles the 2019 release
published with the paper below (Mendeley Data `kgzdz6ryyk`), built for Windows with one small patch (file renames
that replace their target, as on POSIX; without it, runs with Lanczos vectors on disk fail on Windows). It is a
separate program that physim runs, so physim's own MIT licence is unchanged. The source archive is in
`third_party/kshell/kshell.tar.gz` and in the installed `tools\kshell\source\`, and the licence text in
`installer/windows/licenses/GPL-3.0.txt`.

If you publish numbers from a shell-model calculation run with physim, cite KSHELL and the interaction:

- N. Shimizu, T. Mizusaki, Y. Utsuno and Y. Tsunoda, "Thick-restart block Lanczos method for large-scale
  shell-model calculations", Comput. Phys. Commun. 244, 372 (2019);
- the interaction's paper, from the table above.
