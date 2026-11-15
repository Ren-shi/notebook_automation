# Request: a published Coulomb-excitation measurement

For backlog item 43. physim's first-order semiclassical Coulomb excitation (`physim.nuclear.coulex`) passes its
closed-form and engine-orbit checks, and agrees with GOSIA (`coulex_gosia.csv`, first order, and
`coulex_gosia_all_orders.csv`, all orders; both from `coulex_gosia.inp`). It still needs comparing with data.

Any first-order case with a quoted cross section or excitation probability against angle, with its beam, energy,
target, state and B(E2). Classic choices are ⁵⁸,⁶⁰Ni or ¹¹⁴Cd with α particles or ¹⁶O. Transcribe the table into
`tests/reference/nuclear/coulex_literature.csv` with the columns
`beam,target,beam_energy_mev,state_mev,b_e2_e2b2,theta_cm_deg,probability,dsigma_cm_mb_sr,uncertainty`, and record
the reference (authors, journal, table) in the header.
