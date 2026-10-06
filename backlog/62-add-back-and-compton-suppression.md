# 62 · Add-back and Compton suppression for clovers

**Priority:** P3 · **Size:** M · **Area:** Experiment workbench

**Status: Done.**

> **Done** (`addback`, `addback_factor`, `shield`, `shield_thickness` and `suppression_factor` on the γ-ray
> detector in `experiment.py`; `Response.addback_factor/addback_share/suppression_factor` and
> `Response(..., bare=True)` in `response.py`; `_addback_and_suppression` and `simulate_gammas(..., plain=True)`
> in `gamma_events.py`; the shield as a solid in `scene.py`; `Planner.gamma_modes` and the "plain" spectra in
> `gamma_spectra`; the two switches and `CLOVER_FIELDS` in the app; `tests/python/test_nuclear_addback.py`;
> guide section "Add-back and Compton suppression" in `docs/nuclear-setup.md`).
> - **Add-back:** of the γ rays that leave a Compton deposit, the share (F(E) − 1) P/T ÷ (continuum share)
>   returns to the peak, F(E) = 1 + (F₀ − 1)(E / 1332 keV)^0.8 with F₀ = `addback_factor`; the deposit is summed
>   with the neighbour's (one of the two crystals beside it round the square, since the scattered photon is
>   not transported) and the γ ray goes to the crystal with the larger deposit, which the Doppler correction
>   then uses. The resolution of a summed event is that of two crystals.
> - **Suppression:** the shield keeps every deposit that is not the full energy with probability 1/S(E),
>   S(E) = 1 + (S₀ − 1)(E / 1332 keV)^0.3 with S₀ = `suppression_factor`; the peak is unchanged. The shield's
>   wall (25 mm unless given) widens the housing, which the rates, the events and the scene take as a blocking
>   volume.
> - **Consistency:** the event chain starts from the bare crystals and applies both to the deposits with their
>   own random numbers, so the γ rays are the same with and without; the response model applies the same
>   factors to the efficiency, the peak-to-total ratio and the spectrum shapes (singles, randoms, the room
>   background), so the planned coincidences and the analysis see the add-back too. A measured efficiency
>   curve is taken as the detector runs, with add-back.
> - **Results** (the ⁵⁸Ni example with a clover at 100 mm, 200 000 reactions, 60 emissions each):
>   - The add-back factor of the response at 1332 keV is the value put in (1.5 by default, 1.4 when given)
>     exactly, and the simulated peak grows by F(1454 keV) = 1.54 within statistics; the default lies in the
>     published range for clovers (about 1.5 at 1.33 MeV, Duchêne et al. 1999).
>   - The shield leaves the peak's γ rays and weights identical and lowers the continuum by S(1454 keV) = 3.05
>     within statistics; the response's P/T rises from 0.17 to 0.38 at 1454 keV.
>   - With both off, every column of the simulated γ rays equals item 55's.
> - **Not done here, or to confirm:**
>   - The exponents of the two energy dependences and the shield's thickness are typical values, not from a
>     document; a measured add-back factor or suppressed P/T for the user's clover should replace the defaults.
>   - The scattered photon is not transported (the neighbour is drawn at random), the shield's own spectrum is
>     not made, and a shield does not attenuate the γ rays on their way in (its front is taken as open).
>   - Add-back across two clovers, and the add-back of escape peaks, are left out.

## Why
Clovers are normally used with add-back, which sums the energies of neighbouring crystals and raises the
full-energy efficiency at high energy. Many also sit inside BGO shields that reject Compton-scattered events. Both
change the peak areas and the background. The first version works without them.

## Scope
- **Add-back:** a parametrised probability, depending on energy, that a γ ray scatters from one crystal into a
  neighbour. With add-back on, those events return to the full-energy peak.
  - The Doppler correction then uses the crystal with the larger energy deposit.
- **Compton suppression:** a BGO shield as a blocking volume around the clover (item 51), with a suppression
  factor that depends on energy and lowers the continuum.
- **In the app:** both are switches for each clover, with the spectrum shown with and without.
- **Left out:** transport of the scattered photon; the shield's own spectrum.

## Depends on
53, 55.

## Done when
- The add-back factor (the ratio of peak areas with and without) at 1.33 MeV equals the value put in, and lies
  within the range published for clovers.
- Suppression lowers the continuum by the factor put in and leaves the peak area unchanged.
- With both switched off, the results equal those of item 55.
