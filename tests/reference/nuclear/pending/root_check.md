# Request: open a physim ROOT file in ROOT

For backlog item 42 ("the ROOT file opens in ROOT: TTree::Draw and a histogram fit work"). physim writes ROOT files
with uproot, and the automated tests read them back with uproot. This check needs a real ROOT installation.

1. Make the file (any machine with physim and uproot):

   ```bash
   python -m physim.nuclear.report alpha_on_gold -o alpha-report --seed 1
   ```

2. In ROOT:

   ```cpp
   TFile f("alpha-report/events.root");
   f.ls();                                   // events, spectrum_A20 ... spectrum_A135, setup, info
   events->Draw("measured", "weight*7200*(detector==1 && counted)");   // A30, counts in the 2 h run
   spectrum_A30->Fit("gaus");
   ```

3. Record the result in `tests/reference/nuclear/root_check.md`:
   - ROOT version;
   - whether everything opened without errors;
   - the fitted mean and sigma.

   physim expects a mean of 5.3495 MeV and σ = 13.5 keV, from `Rates(...).peaks("A30")`.
