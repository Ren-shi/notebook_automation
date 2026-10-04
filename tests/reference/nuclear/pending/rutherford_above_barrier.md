# Request: elastic scattering above the Coulomb barrier (literature)

For backlog item 37. The aim is to show, with measured data, that the cross section leaves Rutherford's formula
where physim's grazing-angle warning says it should.

**Suggested source:** G. W. Farwell and H. E. Wegner, "Elastic Scattering of Intermediate-Energy Alpha Particles by
Heavy Nuclei", Phys. Rev. **95**, 1212 (1954), doi:10.1103/PhysRev.95.1212. They measured α scattering on Ag, Ta, Pb
and Th at 60° (lab) from 13 to 42 MeV. Their abstract says the cross section follows Rutherford up to a critical
energy E₀ that rises with Z, then falls as exp[−K(E − E₀)] with K ≈ 0.26 MeV⁻¹.

**What is needed:** for at least Pb (and ideally all four targets), the measured ratio σ/σ_Rutherford against beam
energy at 60°, or E₀ for each target, transcribed from the paper's tables or figures. Save it as
`tests/reference/nuclear/rutherford_farwell_wegner.csv` with columns `target,theta_lab_deg,e_lab_mev,ratio_to_rutherford`
(or `target,theta_lab_deg,e0_mev`), with a `# source:` line naming the table or figure.

**physim's prediction to compare:** with the default interaction radius R = 1.2 (A₁^⅓ + A₂^⅓) + 2 fm, the grazing
angle of α + ²⁰⁸Pb drops below 60° lab (61.0° CM) at **E ≈ 32.5 MeV**, so physim warns from there on. The test will check that E₀ lies near
that energy and adjust the documented default radius if the data say otherwise.
